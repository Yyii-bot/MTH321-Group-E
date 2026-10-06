# Submission Checklist — Status from Role E (Testing & Validation)

Cross-checked against `submission_checklist.pdf` (end-of-project list). Overall
report status is owned by the Visualization & Report role / PM; this file maps
every validation-relevant item to its evidence on branch `testing-validation`.

## Three-part package

| Item | Status | Evidence |
|---|---|---|
| Team report PDF ≥ 8 pages, AI log + one ICS per member as appendices | ● report body 14 pages (Group E, `Robertson_ODE_Report_TeamE_正文.pdf`); AI log Part A entries drafted per member (`report/ICS-rong-yin.md` Part A — needs merging with other four members by PM); ICS Part B drafted (`report/ICS-rong-yin.md`) | appendices pending assembly |
| Code zip with README | ● my package has README + run_all; **blocker:** merge all role branches before zipping | `code/validation/` |
| Slides PDF | ○ owned by PM/Viz role | — |

## Report checks (validation-owned)

| Item | Status | Evidence |
|---|---|---|
| Numbers checked against own constructed oracle, with the method stated | ● §4.1 states the Radau dual-tolerance chain and digit-retention rule; my independent chain agrees | `code/validation/reference.py`, §4.1, Table 4 |
| Observed order stated for ≥ 2 methods (log-log slope) | ● report Table 5 (PR test) + Table 7/Fig. 6 (Robertson); my independent orders 0.99/1.00/5.09(≥4) corroborate | `convergence_study.py`, fig06 |
| Cost measure stated; two-method comparison at matched accuracy | ● §4.2/4.4 (matched error 0.7423, per-step cost ratio 5.16–5.25; work–precision fig. 8); my independent matched-target comparison (err≈7e-5) agrees qualitatively | `final_check.py`, fig07 |
| Adaptive controller demonstrated (step varies, tolerance respected) | ● report fig. 5 (six decades of h; estimate hugs threshold); my E5 sweep agrees | `edge_cases.py`, fig05 |
| Stiffness discussion (eigenvalue spread / ratio) | ● §1.2/§2.3/fig. 2; S(40)=1.5840e5 and |λ|max(40)=3392.79 match my fig01 to all quoted digits | `edge_cases.py`, fig01 |
| Conservation / non-negativity / full-state error as separate diagnostics | ● §4.2/4.4 three-way criteria; the e^T J = 0 exact-preservation theorem and "conserves but wrong" experiment are in §4.4 | `edge_cases.py` E4, fig04 |
| AI Transparency Log present | ◐ my Part A entries drafted; team merge pending | `report/ICS-rong-yin.md` |
| ICS per member, traceable to git | ● mine traces to commits on `testing-validation` | `report/ICS-rong-yin.md` |

## Code checks (my package)

| Item | Status |
|---|---|
| `python code/run_all.py` regenerates every figure | ● 5/5 steps, 8/8 tests, fig01–07 from clean state (log committed) |
| README: dependencies / how to run / how to reproduce | ● `code/validation/README.md` |
| No absolute paths | ● relative paths only; FIGDIR auto-adapts to namespaced layout |
| f(t, y) model-agnostic | ● all solvers take callables; demonstrated on logistic (T4/T5) |
| Seeds fixed | N/A — deterministic ODE project, no random elements |
| Figures reproducible, no stale png | ● all fig01–07 regenerated in the Week 3-4 commit |
| .gitignore covers .venv/__pycache__/build | ● committed at Week 1-2 |

## Report-number cross-check (my independent runs vs the report)

| Claim in report | My independent measurement | Verdict |
|---|---|---|
| y(40) = (0.7158270687, 9.185534765e-6, 0.2841637457) | (0.7158270687, 9.1855347646e-6, 0.2841637457) | match to all quoted digits |
| S(1e-4)=2996, S(1e-2)=5423, S(40)=1.5840e5; \|λ\|max(40)=3392.79 | 3010, 5422, 1.5840e5; 3392.8 | match (pack values 3.0e3/5.4e3) |
| EE threshold ≈5.9e-4; min y2 = −4.84e-6 at h=6e-4 | unstable ≥6e-4; min y2 = −4.842e-6 at h=6e-4 | match |
| RK4 threshold 8.2e-4 | 1e-11 error at 8.2e-4; 4e-2 at 9e-4 | match |
| τN sweep: conserves mass yet solution wrong | defect ≤1e-15 for all τN; 28 % error at τN≥1e-2 | same finding |
| y2(40) cross-check table | my five-level table agrees digit-for-digit at equal h | match |

## Last sanity question (the checklist's own)

Weakest link: the production code is not yet in the shared repo, so the
instructor-runs-the-code check cannot be demonstrated from git history alone.
Fix first: PM merges branches → one fresh-clone `run_all.py` run → zip.
