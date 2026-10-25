# Weekly Milestone — Week 3

**Author:** Rong Yin · **Role:** ⑤ Testing & Validation · **Topic:** ③ Robertson problem

Format: what we completed / what we plan next / current blocker.

---

## Week 3 — Convergence orders + code review

**Completed**
- Convergence-order studies (`code/validation/convergence_study.py`, fig06):
  implicit Euler **0.99 ± 0.00** (theory 1), explicit Euler **1.00 ± 0.00**
  (theory 1, measured inside its stability window h <= 5e-4), RK4
  **5.09 ± 0.15** (theory 4). Every panel follows visualisation-guide Rule 1:
  log-log axes, theoretical-slope reference line, polyfit slope +/- standard
  error, 95% confidence band, and a round-off-floor annotation. The team agreed
  to phrase the RK4 claim as "at least design order", with the narrow usable
  window (stability limit 8.2e-4 above, oracle measurement floor ~1e-13 below)
  stated explicitly as a limitation.
- Finished the edge-case experiments started in Week 2
  (`code/validation/edge_cases.py`):
  - E4 Newton-tolerance sweep (fig04): because the conserved direction
    e=(1,1,1)^T satisfies e^T J = 0, EVERY Newton iterate preserves the linear
    invariant exactly - the measured mass defect is 0 to 1e-15 for any Newton
    tolerance from 1e-1 down to 1e-14 - yet with newton_tol >= 1e-2 the solution
    itself is wrong by 28% at t=40. Quantitative proof that conservation alone
    is weak evidence, exactly as the problem pack warns. Production runs use
    newton_tol = 1e-12, far past the point (<= 1e-4) where the global error
    stops changing.
  - E5 adaptive-tolerance sweep (fig05): for tol >= 1e-4 the step-doubling
    controller never fires (0 rejections, h stays 0.1) yet the achieved global
    error is 1.7e-4 > tol; local-error control is not a global-error guarantee.
    As tol tightens, accepted/rejected steps (1023/8 at tol=1e-6, 10163/12 at
    1e-8) confirm the predicted transient behaviour: h shrinks by two orders of
    magnitude through the stiff layer and grows back afterwards.
- Code-review checklist (`code/validation/CODE_REVIEW_CHECKLIST.md`) compiled
  from Brief Sections 2/7 and the visualisation guide, run over the validation
  code (all domains pass) and applied to the team's implementation at the
  Tutorial-3 review; two minor issues found in team code (timestamped figure
  file names, missing README) were logged as tickets and fixed.

**Plan next (Week 4)**
- y2(40) multi-tolerance cross-check against the reference; cost-vs-accuracy
  at a matched error target; end-to-end `run_all.py` reproduction check from a
  clean state; report-number traceability table; hand validated numbers to the
  Visualization & Report owner.

**Blocker**
- None.
