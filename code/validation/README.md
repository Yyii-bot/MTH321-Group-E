# Robertson Problem - Testing & Validation Code Package

Validation code for Topic 3 (Chemical kinetics: Robertson problem), owned by
Role 5 (Testing & Validation).

## How to run

Dependencies: Python 3.9+, `numpy`, `scipy`, `matplotlib` (already in the
course venv).

```bash
cd code
python run_all.py   # reproduce all 7 figures + all console output
```

Individual modules can also be run standalone:

| Command | Content | Week |
|---------|---------|------|
| `python reference.py` | reference-solution chain L1-L3 (Radau dual-tolerance, LSODA/BDF cross-check) | Week 1 |
| `python validation_tests.py` | 8 unit tests, all PASS/FAIL | Week 2 |
| `python edge_cases.py` | edge hunting: eigenvalue history, Euler/RK4 step sweeps, Newton/adaptive tolerance sweeps | Week 2-3 |
| `python convergence_study.py` | convergence orders for three methods (log-log + theory slopes + 95% bands) | Week 3 |
| `python final_check.py` | y2(40) tolerance scan, cost vs accuracy, report-number traceability | Week 4 |

Figures are written to `../figures/fig01..fig07_*.png` (dpi=150, fixed file
names).

## What is the "oracle"?

The course deliberately provides no benchmark numbers (Brief, Section 2). An
oracle is a trusted reference solution that we construct **ourselves** and
endorse with a credibility chain, which all other methods are then checked
against. Ours: tight-tolerance Radau, repeated with tolerances 100x apart
(only digits confirmed by both runs are quoted), cross-validated against LSODA
and BDF, and sanity-checked against the problem pack's orientation values.

## File overview

- `robertson_model.py` - Robertson RHS, analytic Jacobian, reduced Jacobian,
  eigenvalue/stiffness-ratio helpers, protocol grid
- `solvers.py` - explicit Euler / RK4 / implicit Euler + Newton / step-doubling
  adaptive implicit Euler (with WorkCounter cost ledger)
- `reference.py` - oracle construction (three-level chain)
- `validation_tests.py` - unit-test suite (PASS/FAIL per quoted number)
- `edge_cases.py` - step too large / bad Newton start / loose adaptive
  tolerance / negative concentrations / "conserved but wrong"
- `convergence_study.py` - order verification
- `final_check.py` - pre-submission checks and number traceability
- `run_all.py` - end-to-end reproduction

## Relation to the team's production code

The solvers in this package are a validation scaffold re-implemented
independently (same mathematical definitions as the production solvers owned by
the Algorithm Implementation role). Agreement of two independent
implementations in the fine-mesh limit is the self-consistency check required
by the project brief. The report quotes only numbers confirmed by both.
