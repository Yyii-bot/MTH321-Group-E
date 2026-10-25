# Code-Review Checklist (official Tutorial-3/4 format)

**Reviewer:** Rong Yin (Role E, Testing & Validation)
**Code under review:** (i) team production code as described in the report
(Table 3: methods.py, model.py, adaptive.py, experiments.py, validation_tests.py,
run_all.py, make_figures.py); (ii) my validation package `code/validation/`.
**Scope note (honest status):** at the time of this review the production code
was not yet merged into the shared repo (only README.md on the other branches),
so items on the production code were verified **indirectly** — by reproducing
its published numbers with my independent implementation (`code/validation/`,
branch `testing-validation`) and by checking internal consistency of the
report's tables. Items marked ● verified directly on code in the repo;
items marked ◐ verified indirectly via independent reproduction; items marked
○ pending direct review at merge time.

## A. Can it run?

- ◐ `python code/run_all.py` runs end-to-end without errors — my namespaced
  `code/validation/run_all.py` verified 5/5 steps from a clean state; the
  production `run_all.py` is claimed equivalent in the report (Section 3.2) and
  must be re-verified once merged.
- ● README explains how to install deps and run (`code/validation/README.md`:
  dependencies, per-module commands, outputs). Production-tree README pending
  at merge.

## B. Is the model right?

- ● The RHS f(t, y) matches the report's equations — spot-checked all three
  components of `robertson_model.robertson_rhs` against report eq. (2),
  including the convention that the B+B channel carries exactly 3e7·y2² with no
  factor 2. Two–three terms as instructed: ∂f1 = −0.04y1 + 1e4y2y3 ✓;
  ∂f2 destruction term −1e4y2y3 − 3e7y2² ✓; ∂f3 = +3e7y2² ✓.
- ● Parameters defined in a findable constants block (`robertson_model.py`,
  K1/K2/K3 with comments).
- ● Scales noted: stiffness ratio computed and plotted (fig01);
  |λ|max(40)=3393 and S(40)=1.6e5 stated.

## C. Are the methods correctly implemented?

- ● Explicit Euler is y + h f(t, y), one line (`solvers.euler_step`).
- ● RK4 stages use the correct arguments (k2 at t+h/2 with y+h·k1/2, etc.,
  `solvers.rk4_step`); order 4 confirmed on the logistic sanity problem (T4)
  and on Robertson above the measurement floor (5.09 ± 0.15 ≥ 4).
- ● Implicit Euler Newton solves w − y_n − h f(t_{n+1}, w) = 0, not a
  simplified version (`solvers.newton_solve_ie`).
- ● Newton matrix J_F = I − h J_f with the correct sign; the strongest
  available check is behavioural: from the degenerate guess y2=0 the residual
  must rise first and then decay quadratically — measured exactly so
  (4.0e-3 → 4.8e1 → 5.6e-17 in 13 iterations, test T6). A wrong sign cannot
  reproduce this trace.
- ◐ No hard-coded step sizes where the code claims adaptivity — production
  controller (I-control law, report eq. (27)) differs from my step-doubling
  prototype; acceptance/rejection counts in report Table 6 are consistent with
  a genuine controller (rejections 3–4 across tolerances).
- ◐ Adaptive controller estimates an error and uses it to change h — report
  fig. 5 (right) shows the normalised estimate tracking the acceptance
  threshold from below; direct code check at merge.

## D. Are the numbers believable?

- ● Comparison against a reference: every quoted number traced to the oracle
  (three-level chain; dual-tolerance repeat; LSODA/BDF cross-check).
- ● Observed order consistent with claims: 0.99/1.00 within 0.01 of theory for
  both first-order methods; RK4 claim correctly limited to "≥ design order".
- ● Invariant conservation checked separately from accuracy (T2/T3; E4
  demonstrated that conservation alone is weak: exact conservation coexists
  with a 28% wrong solution at loose Newton tolerance).
- ● Every "stable" claim backed by an experiment, not an assertion: step sweeps
  with per-run stable/non-negative/accurate status (fig02/fig03, report Table 7
  uses the same three-way criteria).

## E. Is it maintainable?

- ● Functions are small and named by what they do (one page per solver;
  `WorkCounter` ledger isolates cost accounting).
- ● f(t, y) is swappable: all integrators take `f` (and `jac`) as callables;
  demonstrated by running the same solvers on the logistic problem (T4/T5).
- ● Magic numbers named (K1/K2/K3; RK4_NEG_LIMIT computed, not hard-coded).
- ● Comments explain why (e.g. the e^T J = 0 invariant derivation in E4).
- ● No dead code in `code/validation/` (unused-import sweep done at Week 4).

## F. Findings

| # | Severity | Location | What's wrong | Suggested fix |
|---|---|---|---|---|
| 1 | major (process) | repo, branches other than testing-validation | Production code is not in the shared repo yet; the submission checklist requires the instructor to run `python code/run_all.py` from the submitted zip | PM coordinates the merge of all role branches into the submission tree before packaging the code zip |
| 2 | minor | report §4.1 vs pack protocol | The report's reference grid uses 241 log points; the pack asks for "at least 200" — compliant, but the number should be stated where the grid is defined | one-line note in §4.1 (already lists "241 logarithmically spaced times" — verify wording at final read) |
| 3 | minor | my `code/validation/README.md` | Figure-output path differs between standalone layout and namespaced repo layout (explained in code, not in README) | add the two-layout note to README at merge |

## G. Top 3 things to fix before Seminar 2

1. Merge all role branches into the submission tree and re-run
   `python code/run_all.py` once on a fresh clone (finding 1).
2. Re-run my `validation_tests.py` against the **production** solvers once
   merged (currently run against my independent re-implementation).
3. Confirm the report's AI Transparency Log (Part A) merges all five members'
   entries — currently only mine is drafted.

## Actionable positive (keep doing)

The team's habit of stating acceptance criteria per run (the
stable / non-negative / accurate triple in Table 7) instead of a vague "it
works" is exactly what the problem pack asks; keep that discipline in any
future project.
