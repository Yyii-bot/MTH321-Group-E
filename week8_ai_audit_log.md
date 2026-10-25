# AI Transparency Log -- Testing & Validation seat

Project 2, Weeks 5-9. One entry per piece of AI-generated code that was
used, per the course's five-gate audit. The gates:

1. transcription -- was the code typed/run as given, unchanged?
2. structure and indexing -- are the index maps and loops what they claim?
3. linearisation -- is the discrete operator the one we derived?
4. independent reproduction -- does a differently-built check agree?
5. interrogation -- can I explain every line in the seminar Q&A?

| file / code piece | gate 1 | gate 2 | gate 3 | gate 4 | gate 5 | changes made and why |
|---|---|---|---|---|---|---|
| week6_mixed_refinement.py (reference FDM solver) | ran as given | dof map checked against `tag` grid by hand for N=8; corner ownership assert | ghost rows re-derived on paper, then checked by substitution (week5_boundary_rows.py) | exact solution plugged into A u* - b per row; corner-adjacent rows caught a sign bug (Dirichlet columns moved to the rhs with the wrong sign) -- fixed BEFORE any results were used | yes, can derive both closures on the whiteboard | sign fix in `add()`; also removed a duplicated Dirichlet contribution found by the row-substitution check |
| week6_lshape_fem.py (reference P1 FEM) | ran as given | mesh area check (must equal 3) guards the mask; RuntimeError on partially-cut cells (this fired once during development -- the cell touching the corner was being kept, fixed by dropping any cell with a removed corner) | stiffness matrix = |T| G^T G re-derived; verified energy rate 0.66 vs theory 2/3 | corner-aware composite quadrature cross-checked at two refinement depths before any L2 number was trusted | yes | midpoint values of u_h in the composite rule are exact averages (P1); a version with zeros there was caught by the quad-check column and reverted |
| week8_depth_checks.py tolerance sweep | ran as given | -- | residual-based stopping criterion (not iterate difference), as the brief requires | sweep reproduced the plateau behaviour independently at two mesh sizes | yes | -- |
| GMRES calls | ran as given | -- | -- | cross-checked against spsolve at rtol=1e-12 (agreement 1e-9 in max norm) | yes | -- |

## Pre-computation worksheets

- Week 5, before the first solve: predicted ghost closure would show an
  O(h) local row residual but O(h^2) global error; predicted first-order
  closure would cap everything at O(h). Both predictions were written
  down BEFORE running week6_mixed_refinement.py, and both were confirmed.
- Week 6, L-shape: predicted energy rate 2/3 (confirmed, 0.66) and L2
  rate "4/3 minus some loss" -- measured 1.0, which is inside the
  h^{4/3-eps} guarantee. The worksheet records this as a measured gap to
  discuss, not a failure to hide.

## Things the AI got wrong and we caught

1. Robin first-order row test originally subtracted q AND r (double
   counting, since r = q + 2u*). Caught because the residual refused to
   shrink; fixed in week5_boundary_rows.py.
2. Dirichlet-to-RHS sign error in the reference solver (week 6). Caught
   by substituting u* into the assembled rows -- the corner-adjacent rows
   were off by exactly 2x the Dirichlet term.
3. A mesh cell covering the removed quadrant near the re-entrant corner
   (L-shape). Caught by the area check and the partial-cell RuntimeError.
4. Composite-quadrature L2 numbers were trusted too early; the quad-check
   column (3 vs 4 corner levels) exposed ~50% quadrature pollution on
   coarse meshes before it contaminated any reported rate.
5. The Week-6 L2 column for the L-shape was integration-contaminated and
   had to be re-measured. The trap: the 3-point corner-composite rule and
   a coarse Gauss rule agreed with each other to <1e-3 on the *same*
   mesh, so they looked mutually validated -- yet both drifted ~50%
   between refinement levels. Reason: the error signal scales like h^4
   while the quadrature error of e^2 scales like h^3 (grad e is O(1) at
   the corner, so D^3(e^2) is O(1)), so the integration error *grows*
   relative to the signal as h shrinks. Caught by the 3-vs-4-level
   self-check column, fixed with the Duffy + tensor-Gauss integrator
   (3-vs-4 check now 1e-6..1e-5, independently validated against scipy
   dblquad). All previously quoted L2 rates for the L-shape retracted;
   energy rates unaffected.
6. A zip-unpacking order bug in the Gauss rule (`for ws, s in zip(...)`
   swapped, so the abscissa was receiving the weight and vice versa --
   every weight came out zero/one). No test inside the script caught it
   because the rule still integrated constants exactly. Found only by
   the independent scipy dblquad validation of the integrator; after
   the fix the Duffy rule matches dblquad to machine precision on
   (1+x+2y)^2 over the unit right triangle (2.0833...).
