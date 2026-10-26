# AI Transparency Log

**Scope.** AI assistance was used for implementation drafting, numerical-experiment design, debugging, and code-facing documentation. AI output was treated as a draft, not as evidence. Acceptance required analytic substitution, structural checks, independent numerical checks, or reproducible experiments. The attached five-gate audit is `Code_Review_and_Five_Gate_Audit.md`.

**Prediction control.** `Precomputation_Predictions.md` was retained as the hypothesis record. Its SHA-256 is written to `results/run_manifest.json` by the full driver. Outputs were not used to rewrite those predictions.

## Substantive changes to AI-assisted output

| Stage | AI-output issue or unsupported suggestion | What was changed, and why | Independent evidence |
|---|---|---|---|
| FEM error evaluation | The first element-value broadcasting pattern did not match scikit-fem's axis convention. | Replaced `values[mesh.t]` with `values[mesh.t].T` so vertex values align with element quadrature. | Linear-field reproduction and finite singular-corner integration tests. |
| Common-node FDM/FEM comparison | `np.nonzero` returned array order, not matrix-unknown order; this could compare values at different physical nodes. | Sorted coordinates through `k_of` and added a runtime coordinate-order assertion. | `test_fdm_fem_common_node_coordinates`; `E3_stencil_equivalence.csv`. |
| FDM/FEM prediction | The draft prediction expected different matrix sparsity. | Rejected the prediction after deriving and measuring vanishing P1 hypotenuse coupling on the chosen triangulation. | Zero pattern mismatches; relative scaled-stencil difference at round-off. |
| Singular energy integration | Recursive dyadic corner subdivision was slow and unreliable. | Integrated all radial powers analytically and delegated only the remaining smooth angular integral to adaptive quadrature. No exact gradient is evaluated at the corner. | `test_corner_polar_error_integration_is_finite_and_accurate_l2`. |
| Sparse matrix audit | An off-diagonal maximum could include implicit sparse zeros and obscure the intended sign diagnostic. | Restricted the diagnostic to stored off-diagonal `.data`. | Matrix sign/connectivity tests and `E1_matrix_audits.csv`. |
| Algebraic accuracy | The initial draft treated scaled-residual tolerance `1e-10` as sufficient for matched accuracy. At `N=128`, Jacobi/GS forward algebraic errors were 8.41%/4.21% of discretisation error. | Added successively tighter tolerances and the explicit target `||x_iter-x_direct||inf / E_direct,inf <= 0.01`. | `E2_accuracy_target_sweep.csv`; loose-tolerance contamination is disclosed. |
| Brief compliance | A forward-error ratio alone did not demonstrate that the reported discretisation error and physical quantity of interest stop changing. During the formal run, Jacobi at `1e-12` still had an `L2h` successive change of 1.130% of the direct error scale. | Added `L2h`, `Linf`, and `u_h(1,1/2)` at every tolerance; measured successive changes against direct discretisation-error scales; selected reporting tolerances only when both forward and 1% stagnation gates pass. Extended the formal scan to `1e-13` rather than relaxing the failed criterion. | `E2_tolerance_stagnation.png`, `E2_accuracy_target_summary.csv`, and the new solver regression test. |
| Cost comparison | A common residual baseline was initially presented too close to a matched-accuracy comparison. | Kept baseline iteration/cost data, but added a separate accuracy-and-stagnation matched cost figure and selected-tolerance table. | `E2_cost_comparison.png` versus `E2_matched_accuracy_cost.png`. |
| SOR parameter choice | The pure-Dirichlet model omega was initially a tempting “optimal” value for a different matrix. | Used it only to centre a measured scan; selected the value with the smallest converged iteration count. | `E2_omega_scan.csv`; small-mesh spectral audit. |
| Iterate stopping | The prediction suggested that iterate-difference stopping might stop too early. | Reported the measured negative result: in the configured audit it stopped later with a smaller residual. Residual stopping remains primary. | `E2_stopping_criterion_audit.csv`; unchanged predictions and reconciliation in the code-review audit. |
| Grading claims | Early wording implied that corner grading should restore all error norms and that the theoretical radial exponent should be optimal. | Limited claims to measured H1/L2 improvement; disclosed the lower graded Linf exponent and strong radial anisotropy. | `E5_error_vs_dof_fits.csv`, `E5_graded_mesh_convergence.csv`, mesh-quality fields. |
| Refinement graphics | Initial plots showed curves and reference slopes but did not state fit quality on the figure. | Added last-three-level fitted order `p` and maximum natural-log fit deviation `delta_log` to E1/E3/E4/E5/E6 refinement graphics. | Regenerated figures and matching `*_asymptotic_fits.csv` / `E5_error_vs_dof_fits.csv`. |
| Experiment generality | Draft drivers contained formal-run assumptions in representative mesh and kappa plotting logic; reference slopes were anchored to arbitrary coarse maxima far above the data. | Used the maximum selected mesh for error maps, the measured tolerance for plots, actual kappa annotations, and reference slopes anchored to the finest value of the intended target curve. | Quick and selective driver paths; formal full-run regeneration and figure inspection. |
| Code quality | Dead plotting fragments, an invalid groupby fragment, and unused imports reduced reviewability. | Removed invalid/dead code and cleaned imports rather than preserving non-working branches. | `compileall` and the complete unit-test suite. |
| AI-audit naming | The first internal audit used descriptive gates that did not match the course brief. | Rewrote the release record under the exact required gates: transcription; structure and indexing; linearisation; independent reproduction; interrogation. | `Code_Review_and_Five_Gate_Audit.md`. |

## Claims deliberately not made

- No claim that a visually plausible solution verifies the implementation.
- No claim that a small residual alone proves forward accuracy.
- No claim that FDM/FEM sparsity differs when the measured patterns are identical.
- No claim of a direct/iterative crossover outside the measured evidence.
- No claim of optimal convergence in every norm or on anisotropic radial meshes.
- No claim of nonlinear Newton linearisation: the submitted core and extensions are linear.
- No claim of human external review, personal meetings, commits, or team activity not present in the supplied record.

The submitting student must be able to explain and modify these changes. The code-review audit and its viva code map identify the corresponding source and evidence; no unrecorded personal or team activity is claimed.
