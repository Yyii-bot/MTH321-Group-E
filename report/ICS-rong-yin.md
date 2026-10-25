# ICS & AI Transparency Log — Entries by Rong Yin (Role E: Testing & Validation)

> Compiled to the official format of ICS.pdf (Part A team log + Part B personal
> statement). Part A is a **team document**: the entries below are my
> contribution to it; the Project Manager merges all five members' entries into
> one log for the report appendix. Part B is mine alone.
> Git traceability: all my artefacts are committed on branch
> `testing-validation` of `Yyii-bot/MTH321-Group-E` (commits `aa01217`,
> `9a7cf6a`, `02f61b7`, and this update).

---

## Part A — Team AI Transparency Log (Rong Yin's entries to be merged)

### A1. AI tools used (my row)

| Tool | Version / Date | Primary use |
|---|---|---|
| Kimi AI assistant (coding assistant in Kimi Work) | Sep 2026 | Scaffold code for the validation suite (solvers, oracle chain, test harness, plotting boilerplate); debugging; drafting figure layouts. Every numerical claim it produced was re-measured by my own code before use. |

### A2. Significant AI-assisted claim audits

**Item 1 — "The mass-conservation defect of implicit Euler is controlled by the
Newton residual."**
- Prompt/claim: while building the Newton-tolerance sweep, the plan assumed the
  defect would scale with the Newton tolerance (the generic statement in the
  problem pack).
- AI output matched that assumption; risk: accepting a plausible-sounding
  general statement as a theorem for this system.
- Possible falsifier: a sweep showing the defect is independent of tolerance.
- **Independent verification and evidence:** I ran the sweep anyway
  (`code/validation/edge_cases.py`, experiment E4; `figures/validation/fig04`).
  The measurement falsified the generic claim: the defect is 0 to 1e-15 for
  every tolerance from 1e-1 to 1e-14. Deriving why (e^T J = 0 implies
  e^T(I−hJ)^{-1} = e^T, so every Newton iterate preserves the linear invariant)
  turned the wrong generic claim into a stronger system-specific theorem, now
  in the report (Section 4.4).
- Final decision: claim replaced; the report quotes the measured table and the
  derivation, not the AI's (and my own initial) assumption.

**Item 2 — "RK4 converges at order 4.00 on [0, 40]."**
- Risk: fitting all error points including the oracle measurement floor
  (~1e-13, below which one measures the reference's own tolerance jitter, not
  the method) gives a meaningless slope (my first fit gave 3.42 ± 0.70).
- Possible falsifier: error increasing as h decreases near the floor.
- **Independent verification:** error table reproduced twice
  (`code/validation/convergence_study.py`); the floor behaviour is visible in
  the raw numbers. Claim rephrased to "RK4 is at least its design order on this
  problem; exactly 4.00 is not identifiable between the stability limit 8.2e-4
  and the measurement floor".
- Final decision: limitation stated in the report and in my milestones.

**Item 3 — Reference-solution digits.**
- AI-assisted code produced a Radau oracle; risk: treating solver output as
  ground truth.
- **Independent verification:** dual-tolerance repeat (100×), cross-family
  LSODA/BDF runs, and the problem-pack orientation values
  (`code/validation/reference.py`, Level 1–3). Only digits confirmed by all
  checks are quoted (y(40) confirmed to 12/17/12 digits).
- Final decision: accepted with the stated digit limits; this chain is what the
  report's Section 4.1 describes.

### A3. Cross-reference

All A2 items correspond to Part B items 1–4 below; evidence locations are the
files/figures named there.

### A4. Intellectual ownership declaration (my row)

| Name | Signature | Team-owned artefact / contribution | Trusted dependency and interface |
|---|---|---|---|
| Rong Yin | Rong Yin | `code/validation/` package: independent solver set (Euler/RK4/implicit Euler+Newton/adaptive), three-level reference oracle, 8-test PASS/FAIL suite, edge-case experiments, convergence study, final check, run_all.py; milestones Week 1–4; code-review checklist | SciPy `solve_ivp` (Radau/BDF/LSODA) used only to construct the reference solution; interface tested by the dual-tolerance repeat and cross-family comparison (reference.py). NumPy linear algebra inside Newton. |

---

## Part B — Individual Contribution Statement

**Individual Contribution Statement — Project 1 (ODE IVP)**

**Name:** Rong Yin (2363466) · GitHub: Sylvia-OvO11
**Role this project:** E — Testing & Validation

**The specific things I did this project:**

1. Built the team's reference-solution credibility chain
   (`code/validation/reference.py`): Radau at rtol=1e-11/atol=1e-13 repeated
   with tolerances 100× tighter on the protocol grid (t=0 plus 200 log-spaced
   outputs), cross-validated against LSODA and BDF, sanity-checked against the
   problem pack; y(40) confirmed to 12/17/12 digits. All my commits on
   `testing-validation`.
2. Implemented an independent solver set (`code/validation/solvers.py`):
   explicit Euler, RK4, implicit Euler with Newton (`newton_solve_ie`), and a
   step-doubling adaptive controller, with a cost ledger (`WorkCounter`) —
   deliberately independent from the production solvers so that validation
   never depends on the code it checks (self-consistency check, Brief §2).
3. Wrote the unit-test suite (`code/validation/validation_tests.py`): 8 PASS/FAIL
   tests covering y(40) vs oracle, mass conservation, non-negativity with a
   declared threshold, RK4/Euler order checks on a non-stiff sanity problem
   before touching Robertson, the first-step Newton residual-rise-then-
   quadratic-decay trace, eigenvalue check values, and oracle self-agreement.
   Result: 8/8 PASS (log: `figures/validation/run_all_output.txt`).
4. Ran the edge-case experiments (`code/validation/edge_cases.py`): explicit
   Euler and RK4 step sweeps confirming the frozen-Jacobian stability limits
   5.9e-4 / 8.2e-4 (fig02, fig03; including the honest h=6e-4 "stable but
   negative" discrepancy); Newton-tolerance sweep discovering the e^T J = 0
   exact-invariant-preservation property and the "conserves mass yet is 28%
   wrong" failure mode (fig04); adaptive-tolerance sweep showing local-error
   control is not a global-error guarantee (fig05).
5. Ran the convergence studies (`code/validation/convergence_study.py`, fig06):
   observed orders 0.99 ± 0.00 (IE), 1.00 ± 0.00 (EE), 5.09 ± 0.15 (RK4),
   each with theoretical slope, polyfit standard error and 95% band; flagged
   and documented the RK4 measurement-floor limitation.
6. Ran the pre-submission final check (`code/validation/final_check.py`):
   y2(40) cross-checked at five accuracy levels against the oracle; cost-vs-
   accuracy at the matched target err≈7e-5 with the negative result that
   step-doubling costs ~1.6× the f-evaluations of fixed stepping (fig07);
   compiled the report-number traceability table used to cross-check the team
   report's Section 4 numbers against my independent runs (all headline numbers
   corroborated).
7. Verified reproducibility and reviewed code: `python run_all.py` end-to-end
   from a clean state (5/5 steps, 8/8 tests); code-review checklist per the
   official Tutorial-3 format; weekly milestones Week 1–4 committed under
   `milestones/`.
8. Cross-checked the final team report against my reference: the quoted
   reference y(40), stiffness data S(10^-4)=2996, S(10^-2)=5423, S(40)=1.5840e5,
   |λ|max(40)=3392.79, both explicit thresholds, the τN-sweep invariant finding
   and the y2(40) cross-check table all agree with my independent measurements.

**What I would do differently:**
I would push for the production code to be merged into the shared repo earlier.
As it stands, my review of the production solvers had to be conducted through
the report and my independent re-implementation rather than by reading and
running the team's `methods.py`/`adaptive.py` directly from git.

**Anything I should have asked for help with earlier:**
The intermittent GitHub connectivity from my machine cost several retry loops;
asking the PM for a preferred push window or mirror earlier would have saved
time. Also worth asking earlier: whether the team wanted the validation suite
merged into the production `code/` tree or kept namespaced — I chose
`code/validation/` on my branch and would confirm at merge.
