# Implementation-versus-Library Boundary

## Principle

The project implements the numerical choices whose correctness is being assessed—grid ownership, equations, boundary rows, stopping metrics, stationary iterations, error norms, and verification—and delegates general-purpose sparse storage, mature factorisations, and standard finite-element infrastructure. This is deliberate: the brief permits justified library routines and assesses correct application and the ability to explain and modify the code, not unnecessary reimplementation of every algorithm.

## Implemented directly in the project

| Component | Source | Why direct implementation matters |
|---|---|---|
| Exact solutions, forcing, Neumann/Robin data, singular-angle branch | `code/mixed_bvp/data.py` | Makes every analytic identity and sign auditable; prevents hidden benchmark-specific corrections in assembly. |
| Mixed-grid unknown ownership and `k=(j-1)N+(i-1)` | `core.py` | Protects Dirichlet corners structurally and makes transposition/index errors testable. |
| Five-point assembly, strong Dirichlet elimination, both right-edge closures | `core.py`, `dirichlet.py` | These are the core assessed discretisations; every coefficient and RHS sign is visible and modifiable. |
| Jacobi, red-black Gauss–Seidel, and red-black SOR | `solvers.py` | Demonstrates the actual update equations, one red+black sweep per iteration, vectorised colour operations, and residual stopping. |
| Row-scaled residual, componentwise backward error, named norms, pairwise orders, asymptotic fits | `metrics.py` | Prevents residual/error and area/line norm conflation; keeps every reported metric reproducible. |
| Independent exact local residuals and structural audits | `verification.py`, `core.py` | Checks the assembled equations without using the numerical solution as its own oracle. |
| L-shaped active grid and unique-edge discrete energy error | `dirichlet.py`, `metrics.py` | Makes geometry ownership and singular-corner treatment explicit. |
| Pure-Neumann load correction and mean-gauged saddle-point construction | `pure_neumann.py` | Demonstrates the constant nullspace, discrete compatibility, and gauge rather than hiding them in a solver call. |
| Singular-corner error integration | `fem.py` | Uses analytic radial powers and adaptive angular quadrature; avoids evaluating an infinite exact gradient at the origin. |
| Experiment protocols and reproducible evidence | `experiments/e1_...` through `e6_...` | Keeps meshes, tolerances, norms, quantities of interest, cost measures, and output files explicit. |

## Delegated to established libraries

| Library routine | Used for | Justification and limits |
|---|---|---|
| NumPy array operations | Vectorised updates, norms, small diagnostic spectra | Avoids Python loops over unknowns in solvers and metrics. Dense eigendecompositions are restricted to explicitly labelled small meshes. |
| `scipy.sparse` COO/CSR | Sparse matrix storage, conversion, matvec, absolute-value operations | Standard, tested sparse infrastructure; the project still constructs every core stencil coefficient itself. Production solves never densify. |
| `scipy.sparse.linalg.spsolve` / SuperLU | Reference sparse direct solution | A mature nonsymmetric factorisation is more reliable and efficient than a handwritten LU. Backward error, fill, conditioning, and an independent GMRES solve are still audited. |
| `scipy.sparse.linalg.gmres` | Tight-tolerance algebraic cross-check | Provides a differently implemented solver oracle for the same assembled system. It is not presented as a tuned preconditioned performance study. |
| `scipy.sparse.linalg.onenormest` and `splu` | Condition and LU-fill diagnostics | Separates algebraic trustworthiness and direct-solver storage cost from discretisation error. |
| `scipy.sparse.csgraph.connected_components` | Matrix connectivity audit | Independent graph infrastructure checks accidental isolated unknowns. |
| `scipy.integrate.quad` | Adaptive smooth angular integrals at singular-corner elements | Radial powers are integrated analytically by the project; only the remaining smooth 1D integral is delegated. |
| scikit-fem mesh/P1 basis/assembly/condensation/refinement | Independent P1 FEM comparison, graded meshes, pure-Neumann extension | Supplies standard element and mesh infrastructure for a genuinely different discretisation. The weak form, exact data, boundary conditions, gauges, errors, and mesh-quality checks remain explicit in project code. |
| pandas and Matplotlib | CSV evidence and figures | Output infrastructure only; neither determines numerical acceptance. |

Conjugate gradients are deliberately **not** used for the mixed-boundary core: the assembled matrices are measurably nonsymmetric. The pure-Dirichlet FDM/P1 matrices may be symmetric, but using the same SuperLU protocol there keeps the matched comparison controlled.

## Evidence that the code can be explained and modified

1. **Change the physics without changing assembly.** `CoreData` accepts different exact, forcing, Dirichlet, Neumann, Robin, and `beta` functions. The manufactured-data tests check the new identities independently.
2. **Change the boundary closure.** `assemble_core(..., closure="one_sided")` and `closure="ghost"` share the grid and data but construct different rows. Their coefficients and RHS signs have separate tests.
3. **Change the iterative algorithm.** The SOR update is two sparse colour operations; setting `omega=1` gives red-black Gauss–Seidel. Invalid colouring is rejected rather than silently producing a wrong sweep.
4. **Change the acceptance metric.** Residual scaling, backward error, error norms, forward algebraic error, and successive error/QoI changes are separate functions or fields. The E2 tolerance selector can be tightened or given a different explicit budget.
5. **Change the geometry.** Square and L-shaped grids use explicit active-node maps; tests reject an unknown next to an inactive L-shaped node and check radial polygon corners.
6. **Change FEM resolution or grading.** Uniform, red-graded, and radial-graded families call the same P1 solve/error routines. Mesh-quality metrics are reported alongside errors.

The viva code map in `Code_Review_and_Five_Gate_Audit.md` links these modifications to the equations and measured evidence. Unit tests are regression checks, not a substitute for that explanation.
