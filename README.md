# Mixed-Boundary BVP: Role 3 Algorithm Implementation

This is the code-only Algorithm Implementation deliverable for Direction 2. It contains the sparse finite-difference core, direct and iterative solvers, independent verification, controlled numerical experiments, tests, and one-command reproduction. No report or TeX installation is required to run or review the implementation.

The implementation covers mixed Dirichlet/Neumann and Dirichlet/Robin data, one-sided and ghost-point right-edge closures, matched two-dimensional FDM/P1-FEM benchmarks, smooth and singular L-shaped-domain problems, graded meshes, and a pure-Neumann compatibility/nullspace extension.

## 1. Model and swappable data

The highlighted core model is the unit-diffusion Poisson equation

\[
-\Delta u=f\quad\text{on }(0,1)^2.
\]

The left, bottom, and top edges, including all four corners, are Dirichlet-owned. The open right edge has either

\[
\partial_n u=q\quad\text{or}\quad \partial_n u+\beta u=r.
\]

The manufactured data are centralised in `code/mixed_bvp/data.py`:

- `u_exact=exp(x+y)` and `f=-2*exp(x+y)`;
- `q(y)=exp(1+y)`;
- `beta=2` and `r=q+beta*u_exact=3*exp(1+y)`.

`CoreData` keeps the forcing, Dirichlet data, Neumann/Robin data, exact solution, and Robin parameter separate from mesh construction and assembly. A subclass or replacement data object can change the physics without changing `assemble_core`. The highlighted model has unit diffusion, so no unnecessary variable-coefficient abstraction is introduced.

A minimal direct solve is:

```bash
PYTHONPATH=code python - <<'PY'
from mixed_bvp.core import assemble_core
from mixed_bvp.solvers import solve_direct

system = assemble_core(16, "neumann", "ghost")
result = solve_direct(system.A, system.b, system.row_kind)
print(result.converged, result.residuals["scaled_relative_inf"])
PY
```

## 2. Environment

The recorded full run used Python 3.12.12 and the exact package versions in `requirements.txt`.

```bash
python -m pip install -r requirements.txt
```

Required packages are NumPy, SciPy, pandas, Matplotlib, and scikit-fem. PyAMG is not required: Jacobi, red-black Gauss–Seidel, and red-black SOR are implemented directly in `code/mixed_bvp/solvers.py`.

## 3. Reproduce every result and figure

Run from the extracted project root:

```bash
python code/run_all.py
```

The full driver runs E1–E6 and regenerates:

- all numerical tables in `results/*.csv`;
- all 21 figures in `figures/*.png`;
- software versions, command scope, experiment list, and prediction-file hash in `results/run_manifest.json`.

The formal manifest must contain all six experiments and `quick=false`. Wall times depend on hardware, BLAS/LAPACK, and concurrent load; iteration counts, matrix structure, residuals, errors, and fitted orders are the portable numerical evidence.

A smaller smoke run and selective drivers are available for debugging:

```bash
python code/run_all.py --quick
python code/run_all.py --experiments E1 E3
python code/run_all.py --experiments E2
```

Quick and selective runs overwrite the same output folders. They are not formal evidence; rerun the full command before using the outputs for assessment.

## 4. Run the tests

```bash
PYTHONPATH=code python -m unittest discover -s code/tests -v
```

The 22 tests check:

- five-point and boundary-row coefficients, RHS signs, and corner ownership;
- independent exact local residuals and representative convergence;
- red-black colouring, stationary-solver agreement, and small-mesh spectral radii;
- tolerance-sweep forward error, error stagnation, and midpoint-QoI stagnation;
- P1 linear reproduction, common-node ordering, and scaled FDM/P1 stiffness equivalence;
- L-shaped geometry, singular traces, maximum-principle behaviour, and finite corner integration;
- radial polygon geometry and the pure-Neumann nullspace, compatibility correction, and mean gauge.

## 5. Explicit indexing and sparse assembly

The core unknowns are `(i,j)` for `1<=i<=N` and `1<=j<=N-1`, giving exactly `N(N-1)` unknowns. The zero-based matrix map is implemented explicitly as

\[
k=(j-1)N+(i-1).
\]

Right corners never receive a natural-boundary row. Interior rows use the positive `-Delta` five-point stencil. Dirichlet neighbours are strongly eliminated with `b_k <- b_k-a_kl*g_l`, which gives the `+g/h^2` contribution for a standard neighbour.

Assembly constructs every coefficient in COO format and converts to CSR before solution. Production solves do not densify the matrix. Dense iteration matrices and eigenspectra occur only in explicitly labelled small-mesh diagnostics.

## 6. Solver and verification conventions

### Algebraic residuals

For `s_i=sum_j |A_ij|`, the primary iterative stopping metric is

\[
\rho_{\mathrm{scal}}=
\frac{\|D_s^{-1}(b-Ax)\|_\infty}
{\|D_s^{-1}b\|_\infty}.
\]

Iterate differences are diagnostic only. Direct solves additionally record componentwise backward error, conditioning/fill diagnostics, and a tight row-scaled GMRES cross-check. Conjugate gradients are not applied to the nonsymmetric mixed-boundary core matrices.

### Forward accuracy and stagnation

A small residual alone is not accepted as proof of forward accuracy. E2 measures

\[
\eta_{\mathrm{alg}}=
\frac{\|x_{\mathrm{iter}}-x_{\mathrm{direct}}\|_\infty}
{E_{\mathrm{direct},\infty}}.
\]

At every tolerance it also records `L2h`, `Linf`, and the physical quantity of interest `u_h(1,1/2)`. Successive changes in these errors and the QoI are normalised by the corresponding direct-solve discretisation-error scales. A reporting tolerance is selected only when both conditions hold:

1. forward algebraic error is at most 1% of direct discretisation error;
2. every successive `L2h`, `Linf`, and QoI change is at most 1% of its direct error scale.

The loosest first tolerance cannot pass the stagnation gate because it has no predecessor. `E2_accuracy_target_summary.csv` records the selected tolerance, scaled residual, iterations, cost, and both gates. Baseline common-residual cost and matched dual-gate cost are deliberately separate outputs.

### Error norms and observed orders

`code/mixed_bvp/metrics.py` defines all norms in one place:

- `L2h=sqrt(h^2*sum(e_k^2))` on the declared unknown set;
- `Linf=max|e_k|`;
- separate interior and right-edge norms;
- `L2gamma_right=sqrt(h*sum_right(e_k^2))` as a line norm, not an area norm;
- discrete energy error from unique active east/north edge differences of numerical error.

Refinement tables contain at least four meshes. Observed orders include pairwise rates and last-three-level log-log fits with `R2`, fit levels, and maximum natural-log deviation. Refinement figures state the fitted order and fit deviation rather than relying only on nominal reference slopes.

## 7. Experiment map

| Experiment | Purpose | Main outputs |
|---|---|---|
| E1 | Neumann/Robin data with one-sided and ghost-point closures on four meshes | `E1_core_convergence.csv`, `E1_asymptotic_fits.csv`, matrix and independent local-row audits |
| E2 | Sparse SuperLU versus Jacobi, red-black Gauss–Seidel, and red-black SOR | `E2_solver_comparison.csv`, `E2_accuracy_target_sweep.csv`, `E2_accuracy_target_summary.csv`, matched cost and stagnation figures |
| E3 | Separate matched 2D five-point FDM/P1-FEM Dirichlet benchmark | `E3_matched_FDM_FEM.csv`, `E3_stencil_equivalence.csv`, common-node errors, costs, and sparsity |
| E4 | Smooth versus re-entrant-singular harmonic solutions on the L-shaped domain | `E4_Lshape_smooth_vs_singular.csv`, `E4_asymptotic_fits.csv`, error maps |
| E5 | Uniform, red-graded, and radially graded P1 FEM for the singular solution | `E5_graded_mesh_convergence.csv`, `E5_error_vs_dof_fits.csv`, mesh-quality evidence |
| E6 | Pure-Neumann compatibility correction, constant nullspace, and mean gauge | `E6_pure_Neumann_audit.csv`, `E6_small_mesh_nullspace.csv` |

The matched FDM/FEM study is intentionally a separate pure-Dirichlet benchmark so that different boundary formulations cannot be mistaken for an FDM-versus-FEM effect. The L-shaped smooth control distinguishes solver failure from regularity-induced rate loss. Grading results are assessed against actual DOFs and mesh quality, not only an outer mesh parameter.

## 8. Code architecture

- `code/mixed_bvp/data.py`: swappable exact solutions, forcing, and boundary data.
- `code/mixed_bvp/core.py`: explicit index map, COO/CSR assembly, strong Dirichlet elimination, and matrix audits.
- `code/mixed_bvp/solvers.py`: SuperLU, Jacobi, red-black Gauss–Seidel/SOR, fill/conditioning diagnostics, and GMRES cross-check.
- `code/mixed_bvp/metrics.py`: residuals, backward error, named norms, order fits, and discrete energy error.
- `code/mixed_bvp/verification.py`: independent exact local residuals and corner/index audits.
- `code/mixed_bvp/dirichlet.py`: aligned square and L-shaped FDM grids.
- `code/mixed_bvp/fem.py`: real 2D P1 meshes, solves, graded families, and smooth/singular error integration.
- `code/mixed_bvp/pure_neumann.py`: compatibility correction and mean-gauged saddle-point solve.
- `code/mixed_bvp/experiments/`: separated E1–E6 protocols and shared output utilities.
- `code/run_all.py`: complete reproduction entry point and manifest writer.

Assembly, solution, and verification are separate functions. Geometry, physical data, solver tolerances, and experiment selection are not hard-coded together.

## 9. Code-facing documentation and audits

- `docs/Precomputation_Predictions.md`: unchanged hypothesis record hashed by the full driver.
- `docs/Code_Review_and_Five_Gate_Audit.md`: checklist findings, resolutions, viva code map, and the exact transcription, structure-and-indexing, linearisation, independent-reproduction, and interrogation gates.
- `docs/Implementation_and_Library_Boundary.md`: direct implementation versus justified library delegation.
- `docs/AI_Transparency_Log.md`: substantive AI-assisted code changes, reasons, and independent evidence.
- `results/reproduction_audit.json`: clean-extraction full-run and deterministic-result comparison, generated after the final reproduction.

No personal commit history, meeting record, human external review, or teammate contribution is invented in these records.
