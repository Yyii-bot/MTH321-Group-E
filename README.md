# E_validation -- Testing & Validation seat, Project 2

Everything the Testing & Validation seat owns, split by week. All scripts
run on the team venv (numpy / scipy / matplotlib only):

    python week5_verify_analytic.py
    python week5_boundary_rows.py
    python week6_mixed_refinement.py
    python week6_lshape_fem.py
    python week8_depth_checks.py        # slow: GMRES sweeps + N=256 timings
    python week8_smooth_vs_nonsmooth.py  # depth comparison 4: both legs, one command
    python week7_effectivity.py <adaptive.csv>   # once the adaptive driver exists

## Files

| file | week | what it is |
|---|---|---|
| week5_verify_analytic.py | 5 | manufactured data checks (direction 2 + L-shape), no solver involved |
| week5_boundary_rows.py | 5 | the four right-edge rows verified by substituting u* samples; corner-ownership helper |
| week5_milestone.md | 5 | weekly milestone text (completed / plan / blocker) |
| week6_mixed_refinement.py | 6 | reference FDM solver + four-mesh tables, scaled residuals, asymmetry, GMRES tolerance sweep |
| week6_lshape_fem.py | 6 | reference P1 FEM on the exact L-shape mesh, uniform refinement, corner-aware quadrature, corner/field error split |
| week6_milestone.md | 6 | weekly milestone text with the actual measured numbers |
| week7_clinic_readiness.md | 7 | clinic sheet: discretisation/norm/tolerance triple, evidence list, weaknesses, adaptive protocol |
| week7_effectivity.py | 7 | effectivity index (eta / true error) table + figure for the adaptive sequence |
| week8_depth_checks.py | 8 | matched-settings gate, tolerance sweep, direct-vs-iterative at matched accuracy, FDM-vs-FEM gate |
| week8_smooth_vs_nonsmooth.py | 8 | depth comparison 4: smooth vs. singular on the same FEM/meshes/norms; Duffy L2 integrator + theta/rho split |
| week8_ai_audit_log.md | 8 | five-gate AI audit records for the report appendix |
| week8_milestone.md | 8 | weekly milestone text |
| week9_final_audit.md | 9 | the pre-submission checklist (run it top to bottom) |
| week9_milestone.md | 9 | weekly milestone text: final audit executed, package frozen |
| report_testing_validation.md | 9 | the verification & validation report behind every quoted number |
| results_verified.md | 6-9 | every real, reproducible number with its rerun command -- the only file to quote from |
| needs_from_team.md | 7 | table of missing inputs, who provides them, what they block, fallbacks |
| ICS.md | 9 | Individual Contribution Statement draft |

## Conventions used throughout

- We solve **-Delta u = f** (the test bank's sign convention); for
  u* = exp(x+y) that means f = -2 exp(x+y). The problem pack prints
  "Delta u = f" once but its data matches -Delta u = f.
- Residuals are reported **row-scaled** (D_ii = max_j |A_ij|) wherever
  PDE rows (h^-2) and boundary rows (h^-1) coexist. Raw residuals are
  printed once, to justify the scaling.
- All L-shape L2 errors use the Duffy-transform tensor-Gauss integrator
  from week8_smooth_vs_nonsmooth.py (validated against scipy dblquad);
  energy errors use the 3-point rule. The corner 3-vs-4-level check must
  stay below ~5% before a number is quoted. The L2 values printed by
  week6_lshape_fem.py are quadrature-inflated and superseded -- do not
  quote them.
- The reference solvers in week6 are validation scaffolding. They are
  deliberately independent of the main implementation; if the two
  disagree, that is a finding.

## Headline measured numbers (Week 6-8 runs)

- Direction 2: first-order closure -> order ~1 in every norm (both BCs);
  ghost closure -> order ~2 in every norm including the right-edge Linf.
- Scaled residuals 1e-15..1e-14 on all meshes; ||A - A^T||/||A|| ~ 0.12
  (both closures) -> sparse direct, no CG.
- Tolerance contamination: GMRES rtol 1e-8 still shifts Linf error
  ~80%; first clean level 1e-10 (N=128).
- Direct vs. iterative: no crossover for N <= 256 (direct faster by
  19x-214x, unpreconditioned GMRES).
- L-shape, uniform meshes: energy p = 0.66 vs h (theory 2/3), p = 0.33
  vs N_dof; L2 p = 1.25-1.32 vs h (Duffy quadrature, self-check 1e-5),
  rising toward the h^{4/3} guarantee with refinement. The earlier
  ~1.0 L2 rate came from quadrature pollution and is retracted
  (see results_verified.md, section 4).
- Smooth vs. nonsmooth (same solver, meshes, norms): smooth leg L2
  1.65-1.70 / energy 1.00; singular leg L2 1.25-1.32 / energy 0.66.
  Theta/rho split shows the smooth leg's own L2 is pollution-capped
  near 5/3 -- the quantified ladder 2.0 -> 5/3 -> 4/3.
