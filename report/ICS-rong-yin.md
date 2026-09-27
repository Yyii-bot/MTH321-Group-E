# Individual Contribution Statement (ICS) — Rong Yin (Sylvia-OvO11)

**Project:** Numerical Solution of ODE Initial Value Problems — Topic ③ Chemical
kinetics: Robertson problem
**Role:** ⑤ Testing & Validation
**Team:** MTH321 Group E (team members per PM's roster)

## 1. My ownership in one sentence

I own whether our results can be believed: I constructed the team's reference
solution, checked every headline number in the report against it, hunted the
edge cases, and made the whole pipeline reproducible from a single command.

## 2. What I actually did (by week)

**Week 1 — oracle and test plan.**
- Built the three-level reference chain (`code/reference.py`): Radau at
  rtol = 1e-11/atol = 1e-13 repeated at tolerances 100× tighter, retaining only
  digits confirmed by both runs (y(40) confirmed to 12/17/12 digits);
  cross-validated against LSODA and BDF (agreement ≤ 1.4e-11 over the grid);
  sanity-checked against the problem pack's orientation values.
- Wrote the validation plan (`report/00_Plan_and_Analysis.md`) mapping every
  problem-pack validation requirement to a concrete experiment.
- Wrote my independent implementations of the three solvers + Newton
  (`code/solvers.py`) so that validation does not depend on the code it checks
  (self-consistency check per Brief §2).

**Week 2 — unit tests vs oracle + first edge cases.**
- Implemented the 8-test PASS/FAIL suite (`code/validation_tests.py`): y(40)
  componentwise vs oracle, mass conservation, non-negativity, RK4/Euler order
  on the logistic sanity problem before touching Robertson, the first-step
  Newton residual-rise-then-quadratic-decay trace, eigenvalue check values,
  oracle self-agreement. Result: 8/8 PASS.
- Edge-case hunting (`code/edge_cases.py`): explicit Euler and RK4 step sweeps
  verifying the frozen-Jacobian stability limits 5.9e-4 / 8.2e-4 (including the
  honest discrepancy at h = 6e-4: stable but negative); Newton-tolerance sweep
  discovering that eᵀJ = 0 makes every Newton iterate conserve mass exactly —
  so a 28 %-wrong solution can conserve perfectly; adaptive-tolerance sweep
  showing local-error control ≠ global-error guarantee.

**Week 3 — convergence orders + code review.**
- Ran the order studies (`code/convergence_study.py`): implicit Euler
  0.99 ± 0.00, explicit Euler 1.00 ± 0.00, RK4 5.09 ± 0.15 (claimed as
  "≥ design order", with the narrow usable window and oracle measurement floor
  stated explicitly).
- Compiled and ran the code-review checklist (`report/code_review_checklist.md`)
  over all validation code (Brief §7 + visualisation guide) in preparation for
  the Tutorial-3 review; results presented at the tutorial.

**Week 4 — final check and submission support.**
- y₂(40) tolerance scan cross-checked against the reference at five accuracy
  levels; cost-vs-accuracy at a matched error target (7e-5), including the
  negative result that step-doubling costs ~1.6× fixed stepping here.
- Built `code/run_all.py` + README; verified end-to-end reproduction from a
  clean figures directory (5/5 steps, 8/8 tests); compiled the report-number
  traceability table so every quoted figure maps to a figure or printed output.
- Reviewed the team's final report numbers against my outputs; fed the
  validated numbers to the Visualization & Report owner.

## 3. What I can explain (viva readiness)

- Why the reference solution is credible (three-level chain) and why the
  problem pack's numbers are check values, not answers.
- Why the stiffness ratio must exclude the structural zero eigenvalue, and why
  the initial state is degenerate.
- Why hλ overlays are frozen-Jacobian diagnostics and how the step sweeps test
  them (including the h = 6e-4 discrepancy).
- Why Newton's tolerance barely matters for the invariant here (eᵀJ = 0) but
  matters enormously for accuracy.
- Every line of my solvers, including the Newton update sign and the
  step-doubling controller logic.

## 4. Honest limits

- The RK4 observed order (5.09 ± 0.15) is steeper than 4; "exactly 4" is not
  identifiable between the stability limit and the oracle measurement floor.
  We report "at least design order" rather than forcing a cleaner number.
- My solvers are a validation scaffold; the team's production solvers are owned
  by the Algorithm Implementation role. I checked their outputs against my
  independent implementations (they agree to the quoted digits), but the
  production code is not mine.
- Wall-time figures are single-machine, single-run measurements; nfev/Newton
  counts are the robust cost measure and wall time is indicative only.

## 5. AI transparency (for the report's AI log appendix)

I used an AI assistant as a programming aid throughout, following the course
guidance ("plan first, then let AI help step by step"): I wrote the validation
plan, experiment designs and acceptance thresholds myself first; AI assistance
was used for individual coding steps (boilerplate, plotting layouts, debugging).
Every AI-produced numerical claim was re-derived or re-measured by my own code
before entering the report — all quoted numbers come from `run_all.py` output
on this project, not from AI output. Where AI-generated plotting or solver code
diverged from the plan (e.g. a wrong axis sum in the conservation check, an
over-complex failure-point plotting expression), I corrected it and re-ran the
affected tests. I can explain and modify every line of the final code.

## 6. Signature

Name: Rong Yin  ·  Date: 2026-09-27  ·  Role: ⑤ Testing & Validation
