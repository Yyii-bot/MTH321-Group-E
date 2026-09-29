# Robertson Stiff ODE Project — Code Package

Reproduces every experiment, table, and figure in the report
*Numerical Solution of the Robertson Chemical Kinetics Problem*.

## Layout
| file | purpose |
|---|---|
| `model.py` | Robertson right-hand side, full/reduced Jacobians, Prothero–Robinson test problem |
| `methods.py` | explicit Euler, classical RK4, implicit Euler with Newton iteration |
| `adaptive.py` | adaptive implicit Euler with step-doubling error control |
| `experiments.py` | drivers for each experiment (returns plain data) |
| `run_all.py` | runs all experiments, writes `results/` (JSON + NPZ) |
| `make_figures.py` | reads `results/`, produces `figures/fig_*.pdf/.png` |
| `validation_tests.py` | PASS/FAIL unit tests (T1–T8) |

## Quick start
```bash
pip install -r requirements.txt
python validation_tests.py   # ~1 min; exits non-zero on failure
python run_all.py            # ~5–8 min; fills results/
python make_figures.py       # produces the report figures
```

## Notes
- Python 3.10+; NumPy 2.x, SciPy 1.11+, Matplotlib 3.8+.
- The reference solution is built with `scipy.integrate.solve_ivp`
  (Radau, rtol=1e-12, atol=1e-14) and cross-checked by a run with
  tolerances tightened by a factor of 100; only digits stable under that
  tightening are quoted.
- Core integrators are implemented by the team; SciPy is used only as an
  independent oracle and as a production-solver benchmark (Radau, BDF).
