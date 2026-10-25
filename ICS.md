# Individual Contribution Statement

**Project 2 -- Numerical Solution of Elliptic BVPs, Weeks 5-9**
Team direction: Direction 2 (mixed Dirichlet/Neumann/Robin data, FDM,
two boundary closures) + Open Direction 1 (L-shaped domain with adaptive
refinement); all four depth comparisons.
My seat: **Testing & Validation**.

> Final version. Teammate names in the "What I did not do" section are
> left as role placeholders because this file is mine to sign, not
> theirs.

## My responsibility

I owned the question "can the numbers in the report be believed". That
meant the manufactured-solution pipeline, the four-mesh refinement
protocol with named norms and fitted orders, separating algebraic error
from discretisation error, and hunting the edge cases that everyone else
was too busy building things to look for.

## What I actually did

**Week 5.** Verified all manufactured data before any solver existed:
re-derived f = -2 exp(x+y) with a stencil instead of trusting my own
calculus, checked boundary data signs (q > 0 nails the outward normal),
unit-tested the L-shape branch convention including the quadrant-IV
points, verified f = 0, the two zero-trace rays and the nonzero outer
data. Verified all four right-edge discrete rows by substituting sampled
values of u*, and added a corner-ownership assertion so the Dirichlet
corners can never be overwritten by a boundary row. Caught one bug in
my own test code while doing this (the Robin row was double-counting q).

**Week 6.** Built a small *reference* FDM solver (separate from the main
implementation) so the tables could exist early and the main solver
would have something independent to disagree with. It caught a real
sign error: Dirichlet columns were being moved to the rhs with the wrong
sign, and only the row-substitution check found it -- the scaled
residuals were machine-zero the whole time, which is exactly why the
course says residuals alone prove nothing. Produced the four-mesh tables
for all four closure/BC combinations: first-order closure lands at
order ~1 in every norm, ghost closure at ~2 including at the boundary
nodes. Split boundary vs. interior error (the first-order closure's
defect dominates globally). Quantified the matrix nonsymmetry
(~0.12) and recorded why that kills CG as an option. Ran the
tolerance-sweep separation of algebraic vs. discretisation error.
For the L-shape I built a reference P1 FEM on exactly-represented
meshes: energy error 0.66 vs h (theory 2/3), and the L2 rate first
measured at ~1.0. That L2 number turned out to be quadrature pollution
(see Week 8/9), not physics; the re-measured rate is 1.25-1.32, rising
toward the h^{4/3} guarantee. The episode also showed that plain 3-point
quadrature mis-measures the L2 norm by up to 50% near the corner, so
every L2 number we report now comes from the Duffy-transform
tensor-Gauss integrator with a 3-vs-4-level cross-check.

**Week 7.** Assembled the clinic readiness sheet; froze the
discretisation/norm/tolerance triple and the matched-settings protocol
for the depth comparisons. Designed the effectivity-index study
(eta / true energy error on every adaptive level) that our adaptive
claim rests on.

**Week 8.** Ran the tolerance sweep used as the "tolerance
contamination" evidence (GMRES at rtol=1e-8 still corrupts the error by
~80%; rtol=1e-10 is the first clean level). Ran the direct-vs-iterative
comparison at matched accuracy: no crossover in the tested range, direct
wins by 19x-214x, reported as such. Completed the five-gate AI audit
for every AI-generated solver piece and wrote the audit log. Reproduced
every table from a clean shell.

**Week 9.** Ran the final audit end to end: clean-shell reproduction of
every script (all table numbers reproduced, timings excepted), the
Ctrl-F number-tracing pass over the report, the final residual audit
(1e-15..1e-14 scaled on every mesh), and the package freeze. Finished
the smooth-vs-nonsmooth correction from Week 8: re-measured the L-shape
L2 column with the Duffy + tensor-Gauss integrator and completed the
theta/rho decomposition, giving the final quantified ladder
2.0 -> 5/3 -> 4/3. Wrote the Testing & Validation report, finalised this
ICS. Wrote the diff-script scope for the main-solver cross-check so it
is a 30-minute job the day the interface lands.

## What I did not do

- The main FDM/FEM implementation belongs to [Algorithm Implementation
  teammate]; my reference solvers are validation scaffolding, clearly
  labelled as such, not the deliverable solver.
- Weak forms, truncation derivations and regularity theory belong to
  [Mathematical Theory teammate]; I checked their discrete rows against
  manufactured values but did not derive them.
- Report/slides typesetting belongs to [Visualization & Report
  teammate]; I supplied the verification tables, residual audits and
  figure-ready data, and proofread every accuracy claim.
- Repo, meetings and weekly records belong to [PM teammate].

## The part I would defend hardest in the Q&A

The ghost-closure story: O(h) local defect, O(h^2) global error, and the
evidence chain that shows both (row substitution, then solution-error
tables). Also the scaled-residual argument, and the honest L2 gap on the
L-shape -- measuring 1.0 where the optimistic quote is 4/3, and being
able to say why that's still consistent with theory.

## Reflection

The moment that rewired how I think about verification was finding the
Dirichlet sign error while every residual sat at 1e-15: a number can be
exactly self-consistent and still wrong, which is why "two low-order
methods agree with each other" turned out to mean nothing for the L2
column either. Catching my own wrong numbers twice in one project was
not embarrassing, it was the job -- the evidence chain is only worth
anything if it can retract a number as loudly as it quotes one.

---

Name: Rong Yin   Date: 2026-10-25
