# Weekly Milestones — Week 1 & Week 2

**Author:** Rong Yin · **Role:** ⑤ Testing & Validation · **Topic:** ③ Robertson problem

Format: what we completed / what we plan next / current blocker.

---

## Week 1 — Reference oracle + test plan

**Completed**
- Built the team's three-level reference-solution chain (`code/reference.py`,
  Level 1–3): Radau at rtol=1e-11/atol=1e-13 repeated with tolerances 100× tighter
  (rtol=1e-13/atol=1e-15) on the protocol grid (t=0 plus 200 log-spaced outputs,
  1e-8 ≤ t ≤ 40); the two runs agree pointwise to ≤ 2.1e-12 and confirm
  y(40) = (0.71582707, 9.1855e-6, 0.28416375) to 12/17/12 digits.
- Cross-family validation: LSODA and BDF at the same tight tolerance agree with
  the Radau reference to ≤ 1.4e-11 over the whole grid; reference conserves mass
  to 1.6e-15 (round-off).
- Wrote the validation plan mapping every problem-pack validation requirement to
  a concrete experiment with acceptance thresholds (kept in our team notes;
  summary in code docstrings).
- Implemented my own independent solver set (explicit Euler, RK4, implicit
  Euler + Newton, step-doubling adaptive) so validation never depends on the code
  it checks — this is the self-consistency check required by the project brief.
- Eigenvalue check values reproduced along the reference trajectory:
  S(1e-4)=3.01e3, S(1e-2)=5.42e3, S(40)=1.58e5, |λ|max(40)=3.393e3, all within
  0.4 % of the problem pack's orientation values (used as sanity checks only).

**Plan next (Week 2)**
- Turn every headline number into a PASS/FAIL unit test against the oracle;
  start the step-size sweeps (explicit Euler / RK4) and the Newton / adaptive
  tolerance sweeps.

**Blocker**
- None. (Will point the oracle-vs-solver tests at the Algorithm Implementation
  role's production solvers once they are merged.)

---

## Week 2 — Unit tests vs oracle + edge-case hunting (part 1)

**Completed**
- Unit-test suite (`code/validation_tests.py`): **8/8 PASS**.
  - T1 implicit Euler h=0.1 at t=40 matches oracle per component
    (rel. err 4.9e-4 / 1.5e-3 / 1.2e-3; course sanity bound < 1 %).
  - T2 mass-conservation defect 4.4e-16.
  - T3 non-negativity over the whole trajectory (declared threshold −1e-12;
    note the reference itself dips to −2.6e-18 at round-off level).
  - T4/T5 RK4 and explicit Euler verified at their design orders on the logistic
    sanity problem *before* being used on Robertson.
  - T6 first-step Newton residual trace: rises from 4.0e-3 to 4.8e1 (degenerate
    guess y₂=0) then converges quadratically in 13 iterations — a built-in
    correctness check of the Newton implementation.
  - T7 eigenvalue/stiffness-ratio check values; T8 oracle self-agreement ≥ 8
    digits.
- Edge cases (`code/edge_cases.py`):
  - E1 eigenvalue/stiffness-ratio history vs log time (fig01), excluding the
    structural zero eigenvalue; degenerate initial state handled per protocol
    (plot starts at first positive log time).
  - E2 explicit Euler step sweep: instability onset matches the frozen-Jacobian
    limit h < 5.9e-4; blow-up time moves later as h decreases (t≈0.022 at
    h=2e-3 → t≈27 at h=7e-4), as a state-dependent |λ|max predicts. Discrepancy
    found and logged: at h=6e-4 the run stays finite but goes negative
    (min y₂=−4.8e-6) with error jumping to 3.2e-3 — an instability precursor
    *below* the frozen boundary, reported honestly.
  - E3 RK4 step sweep: onset matches h < 8.2e-4 (|R(−2.785)|=1); at h=9e-4 the
    error jumps nine orders of magnitude.

**Plan next (Week 3)**
- Convergence-order studies for the three methods on log-log axes with
  theoretical slopes and polyfit confidence bands; run the code-review
  checklist over the team's code before Tutorial 3; finish the Newton-tolerance
  and adaptive-tolerance sweeps (E4/E5, already stubbed in `edge_cases.py`).

**Blocker**
- RK4's usable h-window is narrow (stability limit 8.2e-4 above, oracle
  measurement floor ~1e-13 below). Proposing to phrase the RK4 order claim as
  "at least design order 4" and state the limitation explicitly — to confirm
  with the team.
