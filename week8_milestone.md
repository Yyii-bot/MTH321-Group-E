# Week 8 -- milestone notes (Testing & Validation)

## What I did this week

Closed the verification evidence for the four depth comparisons and the
reproduction check.

**1. Matched-settings gate (`week8_depth_checks.py`).** Every comparison
in the report must share problem, norm, refinement protocol and solver
tolerance. The gate function refuses to compare runs that differ in any
of these -- it caught its first mismatch in the demo (L2 vs Linf), which
is the point.

**2. Tolerance sweep (contamination evidence).** N=128, ghost Neumann,
GMRES at rtol = 1e-2 .. 1e-12 vs. the spsolve reference
(Linf error 1.35e-5):

| rtol | scaled residual | Linf error | contaminated? |
|---|---|---|---|
| 1e-2 | 5.3e-3 | 5.35e+00 | YES |
| 1e-4 | 2.8e-5 | 1.27e-01 | YES |
| 1e-6 | 3.7e-7 | 1.21e-03 | YES |
| 1e-8 | 3.9e-9 | 2.45e-05 | YES |
| 1e-10 | 3.8e-11 | 1.36e-05 | no |
| 1e-12 | 3.4e-13 | 1.35e-05 | no |

The error is still visibly wrong at rtol=1e-8 even though the residual
looks tiny. Stopping criterion used throughout: scaled residual, never
iterate difference.

**3. Direct vs. iterative at matched accuracy (no crossover in range).**

| N | spsolve s | GMRES s (matched tol) | speedup |
|---|---|---|---|
| 64 | 0.009 | 0.17 | 19x |
| 128 | 0.046 | 9.1 | 198x |
| 256 | 0.274 | 58.5 | 214x |

GMRES without a preconditioner never gets close in the tested range --
reported as "no crossover observed", which the brief explicitly accepts
as a finding. Caveat recorded in the script: timings with a
preconditioner would need the preconditioner named.

**4. FDM vs. FEM gate.** The matched-settings gate functions and passed
its demo self-check, but no real FEM numbers have gone through it yet --
stated honestly, the gate is *ready*, the comparison is *pending* both
main solvers.

**5. AI audit log.** `week8_ai_audit_log.md` records the five-gate audit
for every AI-generated solver piece, including the bugs the audit caught
(Dirichlet RHS sign; Robin row double-counting; L-corner mesh cell) and
the pre-computation worksheets. Goes in the report appendix.

**6. Reproduction check.** Fresh shell, plain `python` from the venv:
`week5_*.py`, `week6_*.py`, `week8_depth_checks.py` rerun end-to-end;
every table number reproduced to the last printed digit. Consolidated
into `results_verified.md` (only real numbers, each with its rerun
command).

**7. Smooth-vs-nonsmooth, final leg (`week8_smooth_vs_nonsmooth.py`).**
Both legs now run from one script on the same reference FEM, the same
L-shaped meshes (n = 8..128) and the same error functionals. L2 is
measured with a Duffy-transform tensor-Gauss integrator that was
validated against `scipy.integrate.dblquad` on a known integral before
any error norm was trusted. Measured: smooth leg (u_s = exp(x+y))
L2 order 1.70 -> 1.65, energy exactly 1.00; singular leg L2 1.25 ->
1.32, energy 0.66. Splitting the smooth-leg error into theta (discrete
pollution, exact formula) and rho (interpolation, Duffy-integrated)
gives the story: rho decays at exactly 2.00 while theta decays at
~1.55-1.66 and dominates, so even the smooth problem's L2 rate on the
L-domain is capped near 5/3 by corner pollution. That retracts and
corrects the Week-6 L2 numbers (the old corner-composite rule was
quadrature-inflated; energy numbers unchanged and still valid). The
quantified ladder 2.0 -> 5/3 -> 4/3 is the content of this depth
comparison; it is in results_verified.md section 6.

## Plan for Week 9

Final audit: trace every report number back to a manufactured check or a
refinement table; final residual audit; finish ICS; verify the zip runs
on a clean machine.

## Blocker

None.
