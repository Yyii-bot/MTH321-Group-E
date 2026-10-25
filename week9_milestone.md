# Week 9 -- milestone notes (Testing & Validation)

## What I did this week

Final audit week. Everything on this list was executed, not just planned.

**1. Clean-shell reproduction pass.** Reran every script from a plain
shell on the venv: `week5_verify_analytic.py`, `week5_boundary_rows.py`,
`week6_mixed_refinement.py`, `week6_lshape_fem.py`,
`week8_depth_checks.py`, `week8_smooth_vs_nonsmooth.py`. Every table
number reproduced; wall-clock timings moved a little between runs
(they are timings), orders and error values did not.

**2. Number-tracing pass.** Every number in the report traces:
report -> `results_verified.md` -> table -> script -> command. The
Ctrl-F pass found one stale sentence in the Week 6 milestone and the
audit checklist (the pre-correction "~1.0" L2 rate); both now carry the
correction note pointing at the Duffy re-measurement.

**3. The smooth-vs-nonsmooth correction landed.** The L2 column for the
L-shape was re-measured with the Duffy + tensor-Gauss integrator and the
theta/rho decomposition is in. Final quoted ladder: 2.0 (interpolation)
-> 5/3 (smooth FEM, measured 1.65-1.70) -> 4/3 (singular FEM, measured
1.25-1.32 and rising). Recorded in `results_verified.md` section 6 and
in the audit log as items 5-6.

**4. Final residual audit.** Scaled residuals on all Direction-2 tables:
1e-15..1e-14 on every mesh and closure/BC combination. Raw residuals
shown once to justify the scaling. The nonsymmetry number
(||A - A^T||/||A|| ~ 0.12) re-confirmed at N=64.

**5. Package freeze for submission.** `report_testing_validation.md`
written as the verification evidence behind the report numbers;
`week9_final_audit.md` checkboxes resolved (done items checked, pending
items explicitly marked as team dependencies, not silently dropped);
ICS finalised; README updated so a stranger can rerun everything.

## Plan for the remaining days

Hand the package to the report owner; stay on call for the cross-check
the day the main solver interface lands (the diff script is a 30-minute
job, already scoped in needs_from_team.md). Rehearse the clinic Q&A list
at the bottom of the audit checklist.

## Blocker

None on my side. The three open items (main-solver cross-check,
adaptive CSV, FDM-vs-FEM) all sit with teammates and are tracked in
needs_from_team.md with fallbacks stated.
