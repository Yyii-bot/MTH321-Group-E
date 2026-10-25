# Week 7 -- clinic readiness (Testing & Validation)

## Clinic readiness sheet

**Discretisation / norm / tolerance triple (the thing the clinic exists to
vet):**

| item | our answer |
|---|---|
| problem | -Delta u = f, u* = exp(x+y); mixed BC on unit square (Direction 2); singular harmonic benchmark on L-domain (open direction 1) |
| discretisation | FDM five-point stencil, two right-edge closures (one-sided first-order, ghost-point second-order); P1 FEM on exact L-shape triangulation |
| named norms | discrete L2 (sqrt(h^2 sum e^2)), nodal Linf (all nodes), right-edge nodal Linf, FEM: L2 and energy seminorm via corner-aware composite quadrature |
| refinement protocol | h, h/2, h/4, h/8 minimum; we run 5-6 levels (N=32..256 / n=8..256); orders fitted pairwise and stated per mesh pair |
| algebraic tolerance | direct solve (spsolve) -- no tolerance; backward scaled residual reported per mesh (1e-15..1e-14); iterative runs must pass the tolerance sweep before their numbers count |
| residual reporting | scaled residual D^{-1}(b - Au), D_ii = max_j |A_ij|; raw residual reported too, to show why scaling matters (PDE rows ~h^-2, boundary rows ~h^-1) |

**Checks attached as evidence:**

1. Manufactured data verified independently (Week 5 script output).
2. Boundary rows verified by substitution; corner ownership asserted.
3. Four-mesh tables for all four closure/BC combinations, boundary and
   interior errors separated, observed orders as predicted.
4. Algebraic-vs-discretisation error separated on a representative mesh
   (direct backward residual + GMRES tolerance sweep; the error plateaus
   where it should).
5. Matrix structure: spy plots, ||A - A^T||/||A|| ~ 0.12, solver choice
   (sparse direct) justified by nonsymmetry.
6. L-shape: branch convention unit-tested; f = 0 verified; zero-trace
   rays verified; energy rate 0.66 (theory 2/3); L2 measured 1.0,
   discussed against the h^{4/3-eps} guarantee; quadrature adequacy
   demonstrated (3 vs 4 corner composite levels).

**Known weaknesses to raise at the clinic (better here than in the report):**

- L2 rate on uniform L-shape meshes sits at ~1.0, not 4/3. We believe
  this is the honest measured behaviour in the reachable range; we will
  NOT claim 4/3.
- The adaptive estimator still needs its effectivity evidence (below).

## Adaptive study: effectivity plan (innovation point)

The claim "adaptive beats uniform" is only as good as the estimator.
Since u* is known, every adaptive level gets:

    eff = eta / ||u* - u_h||_energy

Acceptance: eff stays in a bounded band (say 0.2..5) across levels.
`week7_effectivity.py` reads the adaptive driver's CSV (level, n_dof,
energy_error, eta) and produces the table + figure; energy_error MUST be
computed with the same corner-aware quadrature as the uniform study.

Marking rule (to freeze this week): Dörfler marking with theta = 0.5 on
a gradient-jump (residual-based, f = 0 so no volume term) indicator;
newest-vertex bisection or red-refinement -- team decision, record it.

Uniform-vs-adaptive comparison protocol (freeze this week):

- same solver tolerance, same error functional (energy), same
  quadrature depth;
- plot error vs N_dof for both sequences on one log-log figure;
- uniform reference slopes: energy ~ N_dof^{-1/3}; a working adaptive
  sequence should approach N_dof^{-1/2} -- measured, not assumed;
- cost columns: N_dof, nnz, solve wall time per level.

## Plan for Week 8

Run the two depth comparisons we own the evidence for (direct vs
iterative at matched accuracy -- done, see week8 notes; FDM vs FEM gate);
support the smooth-vs-nonsmooth comparison; assemble the AI audit log
(five gates) for every AI-written solver piece.

## Blocker

Adaptive driver CSV does not exist until the FEM seat wires the loop --
needs the mesh-refinement code by mid Week 8 for the effectivity numbers
to land in time.
