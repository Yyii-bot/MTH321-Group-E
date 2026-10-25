# Weekly Milestone — Week 4

**Author:** Rong Yin · **Role:** ⑤ Testing & Validation · **Topic:** ③ Robertson problem

Format: what we completed / what we plan next / current blocker.

---

## Week 4 — Final check, reproducibility, hand-off

**Completed**
- y2(40) tolerance scan (`code/validation/final_check.py` F1), five accuracy
  levels all compared against the team-built oracle (true y2(40) = 9.18553e-6,
  confirmed to 17 digits):
  fixed h=0.4 -> 9.23917406e-6 (err 5.4e-8); h=0.1 -> 9.19906765e-6 (1.4e-8);
  h=0.0125 -> 9.18723117e-6 (1.7e-9); adaptive tol=1e-6 -> 9.18824714e-6
  (2.7e-9); adaptive tol=1e-8 -> 9.18581299e-6 (2.8e-10). The report quotes
  only digits stable across this table.
- Cost-vs-accuracy at a matched error target (`final_check.py` F2, fig07):
  cost measure named first (nfev / Newton iterations / wall time; Jacobian
  assembly and factorisation, plotting and interpreter startup not counted and
  common to both). At err@40 ~ 7e-5: fixed h=0.02 gives err 6.98e-5 with
  nfev=6008; step-doubling adaptive (tol=1e-6) gives err 6.98e-5 with
  nfev=9396. **Negative result reported honestly:** naive step doubling costs
  ~1.6x the f-evaluations of a well-chosen fixed step on this problem; its
  value is not saving work but not needing to know h in advance.
- End-to-end reproduction: `python run_all.py` verified from a clean figures
  directory - 5/5 steps succeed, 8/8 unit tests PASS, 7 figures regenerated
  under fixed file names; the full log is stored at
  `figures/validation/run_all_output.txt` and every number quoted in the report
  is copied from it.
- Report-number traceability table (`final_check.py` F3 output): each headline
  claim mapped to its figure and output line; the Visualization & Report owner
  cross-checked the final report against this table.
- Code package tidied for submission: `code/validation/` (model, solvers,
  reference chain, tests, edge cases, convergence study, final check, run_all,
  README, code-review checklist), `figures/validation/` fig01-07,
  milestones Week 1-4.
- Individual Contribution Statement (ICS) and the Testing & Validation report
  section drafted (validated numbers only), ready for the report appendix and
  the team LaTeX document respectively.

**Plan next**
- Support the Seminar Q&A (viva prep on the oracle chain, the e^T J = 0
  invariant finding, the frozen-Jacobian limits and the negative cost result is
  done).

**Blocker**
- None.
