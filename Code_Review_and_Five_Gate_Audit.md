# Code Review and Five-Gate AI Audit

**Scope:** Direction 2, Role 3 — Algorithm Implementation release  
**Review basis:** supplied BVP checklist and project brief, analytic substitution, manual source inspection, 22 unit tests, formal E1–E6 evidence, and clean-archive reproduction  
**Review mode:** AI-output self-audit with automated checks; no human external review or certification is implied  
**Attachment:** this record accompanies `AI_Transparency_Log.md`, which records what was changed in AI output and why.

## Release record under the brief's exact five gates

| Required gate | Interrogation and acceptance evidence | Result |
|---|---|---|
| **transcription** | Every continuous sign and discrete row was checked against the stated model. Exact identities give `f=-2e^(x+y)`, `q=e^(1+y)`, and `r=q+beta*u=3e^(1+y)` with `beta=2`. All corners are Dirichlet-owned. One-sided rows replace the boundary PDE; ghost rows retain it and eliminate the auxiliary value. Independent substitution checks the complete assembled `A*u_exact-b`, including RHS signs and the doubled inward ghost coupling. | **Pass.** `test_manufactured_forcing_identities`, boundary coefficient/RHS tests, `E1_local_row_consistency.csv`. |
| **structure and indexing** | The unknown set has `N(N-1)` entries and the explicit zero-based map `k=(j-1)N+(i-1)`. Production assembly is COO followed by CSR. Audits check index permutation, connectivity, diagonal/off-diagonal signs, nonsymmetry, and zero right-corner unknowns. FDM/FEM coordinates are compared in matrix-unknown order. Red-black colours derive from the same coordinates. L-shaped active neighbours and radial polygon corners are tested. | **Pass.** `test_index_map_and_corner_protection`, common-node and stencil-identity tests, geometry tests, `E1_matrix_audits.csv`. |
| **linearisation** | The core Poisson problem and all submitted extensions are linear. There is no nonlinear iteration, Newton step, Jacobian approximation, or nonlinear branch to audit. The gate is therefore passed as a **linear-system integrity check**, not by inventing a Jacobian test: assembly constructs one explicit `A,b`; every solver targets that same system; the pure-Neumann saddle point is also linear. Any future nonlinear extension would require a new Fréchet/Jacobian and linearised-solve audit before use. | **Pass — linear problem; no linearisation approximation introduced.** Direct/GMRES and stationary-solver agreement tests. |
| **independent reproduction** | Verification does not rely on regenerating the same formula with the same code. Closed-form local residuals check coefficients independently; SuperLU is cross-checked by tight GMRES; matched FDM and P1 FEM provide a second discretisation on a separate controlled benchmark; smooth/singular controls test regularity claims. The packaged project is also extracted into a clean directory and rerun from its documented entry point. Deterministic CSV fields are compared independently of wall times. | **Pass.** `E1_algebraic_error_audit.csv`, `E3_stencil_equivalence.csv`, `results/reproduction_audit.json`, full-run manifest, 22 tests. |
| **interrogation** | Adversarial review asked: Can a small residual still contaminate the answer? Can a second-order centred interface condition coexist with a first-order complete boundary-row residual? Do FDM/FEM really have different sparsity? Does grading restore every norm? Is radial `mu=3/2` automatically optimal? Does pure Neumann need more than removing one row? The answers come from measurements, and falsified predictions are retained. | **Pass after corrections.** Tolerance/stagnation sweeps, prediction reconciliation, mesh-quality evidence, nullspace/gauge audits, and explicit limitations. |

No open blocker or major code finding remains in the release audit. Scientific limitations remain disclosed; they are not reclassified as passes or hidden as code defects.

## A. Can it run?

- [x] `python code/run_all.py` is the end-to-end entry point for E1–E6.
- [x] `python code/run_all.py --quick` provides a smaller smoke run.
- [x] `python code/run_all.py --experiments E1 E3` supports selective reproduction.
- [x] `README.md` states dependencies, installation, tests, model/data contracts, and run commands.
- [x] `requirements.txt` records the tested package versions.
- [x] Every generated figure has an experiment driver; no report file is needed for reproduction.
- [x] `results/run_manifest.json` records versions, command, complete experiment scope, and prediction-file hash.
- [x] Clean-extraction reproduction checks instructor-style execution rather than relying on the development working directory.

## B. Is the model right?

- [x] The discrete operator represents positive `-Delta`, with forcing and boundary signs matching the continuous specification.
- [x] The right-edge outward derivative is `u_x(1,y)`; positive `q` provides a sign check.
- [x] Robin data satisfy `r=q+beta*u_exact` with `beta=2`.
- [x] All corners, including both right corners, are Dirichlet-owned.
- [x] Exact solutions and data are centralised and swappable.
- [x] Pure-Neumann compatibility uses outward fluxes on every edge and the constant nullspace is measured before gauging.

## C. Is the discretisation right?

- [x] Interior rows implement the stated five-point stencil.
- [x] Dirichlet neighbours are eliminated by `b_k <- b_k-a_kl*g_l`, giving the correct `+g/h^2` contribution.
- [x] Independent manufactured-solution residuals audit the complete rows, not only the solution.
- [x] FDM/FEM common-node coordinates and the scaled stiffness identity are checked.
- [x] The aligned L-shaped grid excludes inactive nodes and needs no cut-cell correction.
- [x] Singular traces, branch, finite corner integration, and mesh incidence are tested.
- [x] Solvers use vectorised sparse operations rather than Python loops over individual unknowns.

## D. Is the solver honest?

- [x] Row-scaled residual stopping is primary for every iterative method; iterate differences are diagnostic only.
- [x] PDE and boundary residuals are separated, and componentwise backward error is recorded.
- [x] Nonconvergence is returned explicitly.
- [x] The tolerance sweep checks both forward algebraic error and successive changes in reported `L2h`, `Linf`, and `u_h(1,1/2)`.
- [x] The reporting tolerance is selected only when the 1% forward target and 1% stagnation target both pass.
- [x] Loose `1e-6` runs are identified as accuracy-contaminated and excluded from convergence claims.
- [x] Direct solves report backward errors, condition/fill diagnostics, and a tight GMRES cross-check.
- [x] CG is not applied to nonsymmetric core matrices.

**Evidence:** `E2_accuracy_target_sweep.csv`, `E2_accuracy_target_summary.csv`, `E2_tolerance_stagnation.png`, `E2_stopping_criterion_audit.csv`, and `E1_algebraic_error_audit.csv`.

## E. Are the numbers believable?

- [x] Neumann and Robin variants each have four-mesh tables for both closures.
- [x] Norm weights and unknown sets are defined in `README.md` and implemented in one `metrics.py` module.
- [x] Interior and right-boundary errors are separated.
- [x] Observed orders use the last three levels and include `R2` and maximum log deviation.
- [x] Refinement figures display fitted order and fit deviation, not only a nominal reference slope.
- [x] FDM/FEM comparisons use identical data, nodes, DOFs, common-node norms, and direct-solver protocol.
- [x] Smooth and singular L-shaped cases use the same aligned FDM solver and four refinements.
- [x] Graded refinement is assessed against DOFs and mesh quality, not only an outer `N`.
- [x] Reduced or unexpected rates and falsified predictions are reported rather than replaced by nominal theory.

## F. Is it maintainable and explainable?

- [x] Assembly, solve, metrics, verification, and experiment drivers are separate.
- [x] Hand-implemented methods and delegated library routines are explicitly justified in `Implementation_and_Library_Boundary.md`.
- [x] Geometry, physical data, and solver tolerance selection are not hard-coded together.
- [x] Dense iteration matrices and eigenspectra are restricted to small-mesh diagnostics.
- [x] SOR omega values are measured in a scan rather than treated as exact model constants.
- [x] Quick-mode tolerances, representative meshes, and kappa annotations follow the selected experiment parameters.
- [x] The viva code map below ties derivations, code decisions, and measured evidence to reviewable source locations.

## Findings and resolutions

| Severity | Finding | Resolution and regression evidence | Status |
|---|---|---|---|
| Blocker during development | FEM error integration broadcast `values[mesh.t]` with the wrong axis order | Changed to `values[mesh.t].T`; linear reproduction and singular integration tests pass | Closed |
| Blocker during development | `unknown_ij()` followed array order rather than matrix-unknown order | Sort coordinates by `k_of`; common-node and stencil-identity tests added | Closed |
| Major during development | Recursive singular-corner quadrature was slow and unreliable | Analytic radial integration plus adaptive angular quadrature; no corner gradient evaluation | Closed |
| Major during release review | Residual tolerance `1e-10` did not ensure the 1% forward target for Jacobi/GS at `N=128` | Added tolerance refinement and forward-error ratios | Closed |
| Major during rubric review | Forward error alone did not demonstrate error/QoI stagnation | Added successive L2h/Linf/midpoint-QoI changes and dual-gate reporting-tolerance selection | Closed |
| Major during evidence review | Weak red grading did not approach intended error-versus-DOF slopes | Retained kappa scan, selected measured `kappa=6`, disclosed remaining Linf limitation | Closed |
| Minor | Sparse off-diagonal diagnostic could include implicit zeros | Audit uses stored off-diagonal data | Closed |
| Minor | Plot drivers assumed formal tolerance/mesh values or omitted kappa labels | Use measured values and maximum selected `N`; actual kappa annotations | Closed |
| Minor | Refinement figures omitted fit deviation | Added fitted `p` and `delta_log` to E1/E3/E4/E5/E6 figures | Closed |
| Minor | Invalid/dead plotting code and unused imports reduced clarity | Removed; compile and unit tests pass | Closed |
| Prediction finding | FDM/FEM sparsity was expected to differ | Measured identical patterns and a scaled stencil identity; prediction rejected | Reconciled |
| Prediction finding | Iterate-difference stopping was expected to show false convergence | Measured later stopping and smaller residual in the configured comparison | Reconciled |

## Viva code map

| Question the implementation must answer | Source and answer anchor |
|---|---|
| Why are there `N(N-1)` unknowns, and why can no right corner receive a Neumann/Robin row? | `index_map` in `core.py`: unknowns have `1<=i<=N`, `1<=j<=N-1`; all corners retain `k_of=-1`. |
| Why does Dirichlet elimination add `g/h^2`? | In `assemble_core`, the unknown equation contains `-g/h^2` on the left for a Dirichlet neighbour; moving it to the right gives `+g/h^2`. Coefficient/RHS tests check the sign. |
| Where do the ghost-row doubled inward coefficient and Robin diagonal come from? | `core.py` and `verification.py`: ghost elimination gives west coefficient `-2/h^2`, diagonal `(4+2*h*beta)/h^2`, and RHS `f+2*r/h` on the right edge. |
| How can the centred derivative condition be second order while the complete ghost row is first order? | `exact_local_residuals` independently evaluates the entire eliminated row. Its measured residual order is about 0.949; global solution-error order is measured separately rather than inferred from the derivative condition. |
| Why is conjugate gradients inappropriate for the core matrices? | `E1_matrix_audits.csv` records nonzero nonsymmetry for all core assemblies. SuperLU is the direct reference, GMRES is an algebraic cross-check, and stationary solvers use the same nonsymmetric system. |
| Why is a small scaled residual not a forward-accuracy proof? | `metrics.py`, `solvers.py`, and E2: the sweep also measures `||x_iter-x_direct||inf/E_direct,inf` and successive `L2h`, `Linf`, and midpoint-QoI changes. Both 1% gates must pass. |
| Why did Jacobi require a tighter reporting tolerance than its forward-only choice? | `E2_accuracy_target_sweep.csv`: at `1e-12` the normalised successive `L2h` change is still about 1.130%, so the scan extends to `1e-13` without relaxing the criterion. |
| Why is FDM/FEM compared on a separate Dirichlet problem, and why can their sparsity match? | `e3_fdm_fem.py` removes boundary-treatment confounds. On the chosen right-isosceles triangulation the P1 hypotenuse coupling vanishes; `E3_stencil_equivalence.csv` measures zero pattern mismatches and scaled stiffness agreement. |
| What distinguishes singular regularity loss from a defective L-shaped solver? | E4 uses the same aligned solver for smooth and singular harmonic data. The smooth control remains about second order, while singular errors concentrate near the re-entrant corner at reduced measured rates. |
| Does grading restore every norm, and why is radial `mu=3/2` not assumed optimal? | E5 measures errors against actual DOFs and records angles/quality ratios. Red grading improves H1/L2 but not every norm; radial meshes become anisotropic and the best tested finite-mesh exponent is measured, not prescribed. |
| Why does pure Neumann require both compatibility correction and a gauge? | `pure_neumann.py` and E6: the constant nullspace requires `1^T b_c=0` for solvability, while the mean gauge selects one solution and imposes the exact continuous mean. |
| What can be modified without rewriting the project? | `CoreData` changes physics; `closure` changes boundary rows; `omega` changes SOR; `accuracy_target` changes the explicit E2 budget; grid/mesh functions change geometry. |

## Final verification commands

```bash
python -m pip install -r requirements.txt
python code/run_all.py
PYTHONPATH=code python -m unittest discover -s code/tests -v
```

Expected test-suite result for this release:

```text
Ran 22 tests
OK
```

The complete run writes a `quick=false` manifest containing E1, E2, E3, E4, E5, and E6. Quick/selective runs are smoke/debug tools, not the formal evidence set. The clean-archive comparison excludes wall-time fields from numerical equality checks because timings are environment-dependent; iteration counts, errors, residuals, orders, matrix structure, and other deterministic fields are compared.

## Actionable practice to retain

Keep the independent, coefficient-level verification chain as the project evolves: exact local residuals, common-node ordering assertions, named norms, prediction hashes, tolerance-to-forward-error checks, and successive error/QoI stagnation checks. These checks detected issues that a plausible plot or a small residual alone would have missed.
