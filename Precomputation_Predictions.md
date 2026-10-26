# Pre-computation Predictions

**Status:** recorded before the full E1–E6 experiment sweeps. These predictions are hypotheses to reconcile with output, not acceptance criteria.

## E1 — Core mixed-boundary systems

1. The one-sided closure will approach first order in right-boundary nodal error; its global error will be limited by that closure.
2. The ghost closure will approach second order in global discrete \(L2_h\) and \(L_\infty\) errors for both Neumann and Robin data, although the full eliminated boundary PDE row has a lower-order local residual than the centred derivative formula alone.
3. All four matrices will be nonsymmetric and have positive diagonals and nonpositive off-diagonal entries. The right corners will have no unknown rows.
4. Direct-solver scaled residuals will be near round-off and much smaller than discretisation errors.

## E2 — Algebraic solvers

1. With the same scaled-residual target, Jacobi and Gauss–Seidel iteration counts will grow approximately like \(N^2\); suitably tuned red-black SOR should grow more slowly.
2. Sparse SuperLU is expected to be fastest on the moderate grids tested. No crossover will be claimed unless measured times actually show one.
3. Iterate-difference stopping may terminate with a larger true residual than residual stopping on fine grids.

## E3 — Matched two-dimensional FDM/FEM benchmark

1. Both methods will show second-order common-node \(L2_h\) convergence on the two smooth data sets.
2. FDM and P1 FEM will have different sparsity and assembly costs despite identical mesh scale and interior degrees of freedom.
3. Any error-constant difference will be reported without attributing it to unmatched data, norms, or tolerances.

## E4 — Smooth versus singular L-shaped solutions

1. The smooth harmonic control will approach second-order \(L2_h\) and nodal convergence.
2. The singular solution will have a reduced nodal rate near \(2/3\) and a reduced global \(L2_h\) rate; errors will concentrate near the re-entrant corner.
3. The discrete energy error will converge more slowly than the discrete \(L2_h\) error.

## E5 — Graded FEM meshes

1. Uniform P1 energy error will scale approximately like \(N_{dof}^{-1/3}\).
2. Corner-graded meshes should improve the energy-error-versus-DOF slope toward \(N_{dof}^{-1/2}\), provided the outer mesh is also refined and quadrature error is controlled.
3. A radial exponent sweep will be measured rather than assuming that the theoretical threshold \(\mu=3/2\) is automatically optimal in a finite mesh family.

## E6 — Optional pure Neumann audit

1. The ungauged stiffness matrix will have one numerical constant null vector.
2. Low-order data quadrature will produce a measurable discrete compatibility defect; the correction will reduce \(1^Tb\) to round-off.
3. The mean gauge will produce a unique finite solution with the prescribed mean.
