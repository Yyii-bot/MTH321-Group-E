# Testing & Validation Package — Project 2 (Elliptic BVPs)

Code-and-evidence deliverable for the Testing & Validation seat (Rong Yin).
This package answers one question for the team report: **can every quoted
number survive the instructor rerunning the code?** It contains the
manufactured-data checks, two deliberately independent reference solvers,
the four depth-comparison studies, the five-gate AI audit log, and the
ledger (`results_verified.md`) that the report is allowed to quote from.
One statement per number, one rerun command per statement.

## Problems and data

Both problems use the test-bank sign convention **−Δu = f**.

**Direction 2 — mixed boundary data** on the unit square. Dirichlet on
x = 0, y = 0, y = 1 (corners included); the right edge x = 1 carries
Neumann ∂ₙu = q or Robin ∂ₙu + βu = r (β = 2), closed either with the
one-sided first-order row or the ghost-point row. Manufactured data
(centralised in `week5_verify_analytic.py`):

- u* = exp(x+y),  f = −2 exp(x+y)
- q(y) = exp(1+y),  r(y) = q + βu* = 3 exp(1+y)

**Open direction 1 — L-shaped domain** (−1,1)² minus one unit quadrant,
area exactly 3. Singular solution u* = r^{2/3} sin(2φ/3) with
φ ∈ [0, 3π/2] (unit-tested in week 5, including the quadrant-IV branch);
smooth comparison leg u_s = exp(x+y). P1 FEM on exactly-represented
uniform meshes, plus the smooth-vs-singular depth study on identical
meshes and norms.

## Environment

Recorded runs: Python 3.12, NumPy, SciPy, Matplotlib (the course venv).
No other dependencies; nothing here needs scikit-fem, pandas, or PyAMG.

## How to reproduce

Run each script from this directory; outputs are console tables intended
to be pasted into `results_verified.md`.

**Fast (seconds):**

```bash
python week5_verify_analytic.py     # manufactured data, all green expected
python week5_boundary_rows.py       # right-edge rows: O(h) local residuals
python week6_mixed_refinement.py    # Direction-2 convergence tables + asymmetry
python week6_lshape_fem.py          # L-shape energy rates (its printed L2 column is superseded, see Conventions)
python week8_smooth_vs_nonsmooth.py # depth comparison 4, both legs + Duffy integrator validation
python week7_effectivity.py <adaptive.csv>   # once the adaptive CSV exists
```

**Slow (minutes):**

```bash
python week8_depth_checks.py         # GMRES sweeps + N=256 timings, ~2-3 min
python week8_iterative_ladder.py     # Jacobi/GS at N=128 need 85k/43k iterations, ~2 min
```

Wall times move between runs and machines; iteration counts, residuals,
errors and fitted orders are the portable evidence.

## Verification gates (all must stay green before a number is quoted)

- scaled residual ≤ 1e-10 at any iterative stop (rtol = 1e-8 was measured
  to leave the solution error 80% wrong — week8_depth_checks.py);
- L-shape L2 quadrature self-check (corner levels 3 vs 4) below ~5%;
- mesh area exactly 3 at every refinement level;
- Dirichlet corner-ownership assertion passes in the assembler;
- any number in the report traces to `results_verified.md` → table →
  script → command.

## File map

| Path | Week | Contents |
|---|---|---|
| `week5_verify_analytic.py` | 5 | manufactured data checks, no solver involved |
| `week5_boundary_rows.py` | 5 | four right-edge rows verified by substitution; corner-ownership helper |
| `week6_mixed_refinement.py` | 6 | reference FDM solver; four-mesh tables; scaled residuals; asymmetry; tolerance probe |
| `week6_lshape_fem.py` | 6 | reference P1 FEM on exact L-shape meshes; energy rates; mesh audit |
| `week7_clinic_readiness.md` | 7 | frozen discretisation/norm/tolerance triple; evidence list; adaptive protocol |
| `week7_effectivity.py` | 7 | effectivity index (η / true error) for the adaptive sequence, CSV in → table + figure |
| `week8_depth_checks.py` | 8 | matched-settings gate; tolerance sweep; direct-vs-iterative; FDM-vs-FEM gate |
| `week8_smooth_vs_nonsmooth.py` | 8 | depth comparison 4; Duffy L2 integrator (validated vs `scipy.dblquad`); θ/ρ split |
| `week8_iterative_ladder.py` | 8 | Jacobi/GS/SOR/CG/GMRES+ILU ladder on the Direction-2 matrix; SOR ω_opt theory check |
| `week8_ai_audit_log.md` | 8 | five-gate audit of every AI-generated solver piece, incl. caught bugs |
| `week9_final_audit.md` | 9 | pre-submission checklist with `[x]` / `[team]` status per item |
| `report_testing_validation.md` | 9 | the V&V narrative behind the report's verification section |
| `milestones.md` | 5–9 | weekly milestones, one section per week |
| `results_verified.md` | 6–9 | every quotable number with its rerun command — quote nothing else |
| `needs_from_team.md` | 7 | missing inputs, who provides them, what they block, fallbacks |
| `ICS.md` | 9 | Individual Contribution Statement (official Part B format) |

## Conventions and status notes

- Residuals are **row-scaled** (D_ii = max_j |A_ij|) wherever PDE rows
  (h⁻²) and boundary rows (h⁻¹) coexist; raw residuals are printed once
  to justify the scaling.
- Nonsymmetry is quoted in the **inf norm**: ‖A−Aᵀ‖∞/‖A‖∞ ≈ 0.12 at
  N = 64 for both closures. (The Algorithm Implementation seat reports a
  Frobenius-norm ratio, a different functional — name the norm wherever
  either is quoted.)
- `week6_lshape_fem.py`'s printed L2 column is quadrature-inflated and
  **superseded** by the Duffy re-measurement in
  `week8_smooth_vs_nonsmooth.py`; its energy rates and mesh code remain
  authoritative (the week-8/9 scripts import them).
- The week-6 reference solvers are validation scaffolding, deliberately
  independent of the main implementation. The cross-check is no longer
  pending: on 2026-10-26 the main solver's ghost-Neumann errors agreed
  with this package's reference solver to six significant figures at
  N = 32/64/128 (L2 and Linf).
- One open discrepancy is recorded in the audit trail (item in
  `results_verified.md` section 4 note and `needs_from_team.md`): the
  main implementation's corner-exact L2 quadrature reports singular-leg
  errors ~3× smaller than this package's Duffy values at identical
  degrees of freedom, same fitted orders. Absolute L2 values from the
  two packages must not be mixed in one table until that head-to-head
  is closed.
