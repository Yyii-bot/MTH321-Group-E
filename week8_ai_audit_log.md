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
| week6_lshape_fem.py (reference P1 FEM) | ran as given | mesh area check (must equal 3) guards the mask; RuntimeError on partially-cut cells (this fired once during development -- the cell touching the corner was being kept, fixed by dropping any cell with a removed corner) | stiffness matrix = \|T\| G^T G re-derived; verified energy rate 0.64-0.66 vs theory 2/3 | three rounds, each caught by a differently-built check (items 5, 8, 9 below): composite 3-point vs composite 3-point; Duffy vs dblquad; final table vs the Algorithm branch's E5 CSV at identical N_dof | yes | corner node made Dirichlet; corner-exact analytic integral + Duffy n=8 replaced all low-order rules; both changes forced by cross-implementation disagreement, not by aesthetics |
| week8_depth_checks.py tolerance sweep | ran as given | -- | residual-based stopping criterion (not iterate difference), as the brief requires | sweep reproduced the plateau behaviour independently at two mesh sizes | yes | -- |
| GMRES calls | ran as given | -- | -- | cross-checked against spsolve at rtol=1e-12 (agreement 1e-9 in max norm) | yes | -- |

## Pre-computation worksheets

- Week 5, before the first solve: predicted ghost closure would show an
  O(h) local row residual but O(h^2) global error; predicted first-order
  closure would cap everything at O(h). Both predictions were written
  down BEFORE running week6_mixed_refinement.py, and both were confirmed.
- Week 6, L-shape: predicted energy rate 2/3 (confirmed, 0.64-0.66) and
  L2 rate "4/3 minus some loss". The worksheet records the measurement
  history honestly: first read ~1.0 (quadrature pollution), then 1.25-1.32
  (free corner node), final 1.31-1.34 -- the textbook 4/3 within a third
  of an order at every level once the measurement is clean.

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
   self-check column, fixed with the Duffy + tensor-Gauss integrator.
   (This integrator was itself superseded on 2026-10-27 by the
   corner-exact + Duffy-n8 protocol of items 8-9; the retraction it
   triggered stays in the record.)
6. A zip-unpacking order bug in the Gauss rule (`for ws, s in zip(...)`
   swapped, so the abscissa was receiving the weight and vice versa --
   every weight came out zero/one). No test inside the script caught it
   because the rule still integrated constants exactly. Found only by
   the independent scipy dblquad validation of the integrator; after
   the fix the Duffy rule matches dblquad to machine precision on
   (1+x+2y)^2 over the unit right triangle (2.0833...).
7. Two prediction misses while building `week8_iterative_ladder.py`,
   both caught by running the numbers instead of trusting the prediction.
   (a) Pre-computation worksheet predicted CG would visibly break down on
   the nonsymmetric mixed matrix; measured: it converged cleanly (748
   iters vs 436 on the SPD matrix at N=128, clean error). The report
   now says "unjustified in theory, survives mild nonsymmetry in
   practice, quantified" instead of the predicted failure. (b) The
   first version ran Jacobi/GS at N=256 (O(N^2) iterations, ~300k)
   and hit the 300 s time budget; the released version caps the
   classical rungs at N=128, states that openly in the output, and
   runs N=256 only for SOR/CG/GMRES+ILU. Also a real code bug in the
   first COO assembler (lists of arrays instead of concatenated flat
   arrays) -- caught the moment the theory table ran, fixed.
8. **The re-entrant corner was left as a free unknown** in the L-shape
   FEM (a modelling slip, not a coding bug). It surfaced only because
   the corner-exact integral ported from the Algorithm branch asserts
   u_h(corner) = 0 and the assertion fired with u_h(corner) ~ 0.29 on
   coarse meshes. A free corner node violates the zero trace of u* on
   the cut edges; it also inflated the smooth leg's apparent "corner
   pollution" (an artificial ~5/3 cap on the L2 rate that we had
   already published in the depth-comparison narrative). Fixing it
   moved the N_dof sequence to exactly the Algorithm branch's
   (161/705/2945/12033), the singular-leg L2 to textbook 4/3, and the
   smooth leg to a flat 2.00 with theta ~ 2.0. The published 5/3 claim
   is retracted in writing in the report and milestones.
9. **Low-order quadrature away from the corner biased the singular-leg
   L2 by ~7x.** The revision-2 Duffy rule with composite corner
   refinement self-checked at 1e-6..1e-5 and still read 2.05e-2 at n=32
   where the certified value is 3.02e-3. Lesson, sharper than item 5's:
   polynomial self-checks cannot see a systematic bias on a singular
   integrand -- every rule we compared against shared the same blind
   spot. What closed the question was not a better self-check but a
   head-to-head against a differently-built implementation (the
   Algorithm branch's E5 table): once the corner node was fixed and
   their corner-exact integral was ported, our table reproduces theirs
   to three+ significant figures at identical N_dof, confirmed
   independently by brute-force 4x-refined Duffy n=8 on the whole mesh
   (agreement 3e-11 on the corner elements). Final protocol: corner
   elements exact radial/adaptive-angular, everything else Duffy n=8
   (their intorder=8), energy by 3-point edge-midpoint; self-check
   column is Duffy n=8 vs n=6 at 1e-7.
