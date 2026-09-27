# Testing & Validation — Report Section (draft for the team LaTeX report)

> For integration into the team report. All numbers come from the fixed output
> of `code/run_all.py` (stored in `figures/validation/run_all_output.txt`);
> figure numbers correspond to `figures/validation/fig01...fig07`.
> Author: Rong Yin (Sylvia-OvO11), Role 5 Testing & Validation.

## V.1 Reference solution: construction and credibility chain

The course deliberately provides no benchmark numbers, so every claim below is
checked against a reference solution we constructed ourselves. The oracle is built
in three levels (`reference.py`):

- **Level 1 (primary).** `scipy.integrate.solve_ivp` with `method='Radau'` on the
  protocol grid (t = 0 plus 200 logarithmically spaced outputs, 1e-8 ≤ t ≤ 40),
  computed twice with tolerances 100× apart (rtol = 1e-11/atol = 1e-13 vs
  rtol = 1e-13/atol = 1e-15). The maximum pointwise difference over the whole grid
  is 2.1e-12, and at t = 40 the two runs agree to 12, 17 and 12 decimal digits in
  the three components. Only digits confirmed by both runs are quoted:
  **y(40) = (0.71582707, 9.18553e-6, 0.28416375)**.
- **Level 2 (cross-family).** LSODA and BDF at the same tight tolerance differ from
  the Radau reference by at most 8.9e-12 (LSODA) and 1.4e-11 (BDF) over the whole
  grid — three independent algorithm families agree.
- **Level 3 (sanity).** The reference conserves mass to 1.6e-15 (round-off) and
  agrees with the problem pack's orientation value to relative errors
  4.4e-8 / 3.8e-6 / 1.6e-7. Its minimum component is −2.6e-18: even the reference
  is negative at round-off level, which is why the non-negativity diagnostic in
  this project uses an explicit declared threshold of −1e-12.

## V.2 Unit tests against the oracle

`validation_tests.py` encodes every headline number as a PASS/FAIL check
(8/8 PASS at submission). Highlights: implicit Euler with h = 0.1 matches the
oracle at t = 40 with relative errors 4.9e-4 / 1.5e-3 / 1.2e-3 (the 1% sanity
bound from the course scaffold); mass conservation defect 4.4e-16; the first
Newton solve of the integration (from the degenerate guess y₂ = 0 at t = 0)
raises the residual from 4.0e-3 to 4.8e1 before converging quadratically in
13 iterations — the residual-rise-then-quadratic-decay trace is a built-in check
that the Newton implementation is correct (a wrong Jacobian cannot reproduce it).

## V.3 Frozen-Jacobian predictions vs. step sweeps

Fig. 1 (fig01) shows the two nonzero eigenvalues and the stiffness ratio
S(t) = max|Re λ| / min|Re λ| over the nonzero spectrum along the reference
trajectory. The full Jacobian always has a structural zero eigenvalue (the
conserved direction), and the initial state is degenerate (both nonzero
eigenvalues vanish at y = (1,0,0)), so the plot starts at the first positive
logarithmic time as the protocol requires. We find S(1e-4) = 3.01e3,
S(1e-2) = 5.42e3, |λ|max(40) = 3.393e3 and S(40) = 1.58e5, matching the
problem pack's check values to ≤ 0.4 %.

The hλ overlay is a local frozen-Jacobian diagnostic, not a nonlinear stability
proof, so we tested it with step sweeps (Figs. 2–3):

- **Explicit Euler** (Fig. 2, fig02): frozen estimate h < 2/|λ|max(40) = 5.9e-4.
  Measured: h ≥ 7e-4 blows up, and the blow-up time moves *later* as h decreases
  (t ≈ 0.022 at h = 2e-3; t ≈ 27 at h = 7e-4) — exactly what a state-dependent
  eigenvalue predicts, since |λ|max grows with t. At h = 6e-4 the run stays
  finite but develops an oscillatory overshoot: min y₂ = −4.8e-6 and the error
  jumps to 3.2e-3, an instability precursor *below* the frozen boundary — a
  discrepancy we report rather than hide. From h = 5.9e-4 downward the runs are
  stable, non-negative, and first-order accurate.
- **RK4** (Fig. 3, fig03): frozen estimate h < 2.785/|λ|max = 8.2e-4. Measured:
  at h = 9e-4 the error jumps nine orders of magnitude (4.3e-2 vs 1.0e-11 at
  8.2e-4); at and below the boundary errors fall to 1e-11…1e-13.

## V.4 Conservation, non-negativity and the "conserved but wrong" trap

The problem pack insists these diagnostics be reported separately, and our
Newton-tolerance sweep (Fig. 4, fig04) shows why: for this system the conserved
direction e = (1,1,1)ᵀ satisfies eᵀJ = 0, hence eᵀ(I − hJ)⁻¹ = eᵀ, so **every
Newton iterate preserves the linear invariant algebraically** — the measured mass
defect is 0 to 1.1e-15 for *any* Newton tolerance from 1e-1 down to 1e-14.
Yet with Newton tol = 1e-1 or 1e-2 the solution itself is wrong by 28 % at t = 40.
A method can conserve perfectly and be completely inaccurate; conservation alone
is weak evidence, exactly as the problem pack warns. For the production runs we
use Newton tol = 1e-12, far past the point (≤ 1e-4) where the global error stops
changing — time-discretisation error dominates throughout.

## V.5 Observed convergence orders

Fig. 6 (fig06) shows error-vs-h on log-log axes with the theoretical slope, the
polyfit slope ± standard error and a 95 % confidence band (visualisation guide
Rule 1): implicit Euler 0.99 ± 0.00 (theory 1); explicit Euler 1.00 ± 0.00
(theory 1, measured inside its stability window h ≤ 5e-4); RK4 5.09 ± 0.15.
For RK4 the usable window between the stability limit (8.2e-4) and the oracle's
measurement floor (~1e-13, below which we measure the reference's own tolerance
jitter) is narrow, and within it the effective slope is slightly steeper than 4
because the error is dominated by the stiff transient. We therefore claim "RK4 is
at least its design order on this problem" and state explicitly that "exactly
4.00" is not identifiable here — an honest limitation.

## V.6 Adaptive controller under loose tolerances

Fig. 5 (fig05) sweeps the step-doubling controller's tolerance. For tol ≥ 1e-4
the controller never fires (0 rejections, h stays at 0.1) because the local
error estimate e ≈ 3.6e-5 at t = 0 is below the budget — yet the achieved global
error is 1.7e-4, *larger* than tol. Local-error control is not a global-error
guarantee; the report quotes global errors measured against the oracle, never
the controller's own estimate. As tol tightens, the accepted/rejected step counts
(1023/8 at tol = 1e-6; 10163/12 at 1e-8) confirm the predicted transient behaviour:
h shrinks by two orders of magnitude through the stiff layer and grows back after
it relaxes.

## V.7 y₂(40) at several accuracies (cross-checked)

Problem pack requirement: the final y₂ reported by the implicit method at several
tolerances, cross-checked against our own reference (true y₂(40) = 9.18553e-6
confirmed to 17 digits):

| run | y₂(40) | error vs oracle |
|---|---|---|
| implicit Euler h = 0.4 | 9.23917406e-6 | 5.4e-8 |
| implicit Euler h = 0.1 | 9.19906765e-6 | 1.4e-8 |
| implicit Euler h = 0.0125 | 9.18723117e-6 | 1.7e-9 |
| adaptive tol = 1e-6 | 9.18824714e-6 | 2.7e-9 |
| adaptive tol = 1e-8 | 9.18581299e-6 | 2.8e-10 |

Only the digits stable across this table are quoted in the conclusions.

## V.8 Cost vs accuracy at a matched error target

Following the brief's cost protocol we name the cost measure first:
right-hand-side evaluations (nfev), Newton iterations and wall time; Jacobian
assembly/factorisation, plotting and interpreter startup are not counted and are
common to both variants. At the matched target err(40) ≈ 7e-5 (Fig. 7, fig07):

| method | err@40 | nfev | Newton iters | wall |
|---|---|---|---|---|
| implicit Euler fixed h = 0.02 | 6.98e-5 | 6008 | 4008 | 0.098 s |
| implicit Euler step-doubling, tol = 1e-6 | 6.98e-5 | 9396 | 6303 | 0.150 s |

**Negative result, reported honestly:** naive step-doubling costs ~1.6× the
f-evaluations of a well-chosen fixed step on this problem (each accepted step
solves three implicit steps; the factor is 3 near the transient and →1 on the
relaxed tail). Its value is not saving work — it is not needing to know h in
advance and crossing the transient automatically.

## V.9 Reproducibility and code review

`python code/run_all.py` regenerates all seven figures and every number quoted
above from a machine with numpy/scipy/matplotlib installed (verified from a
clean figures directory; full log in `figures/run_all_output.txt`). The
code-review checklist (`report/code_review_checklist.md`) was run against the
validation code before Tutorial 3 and against the team's implementation before
submission; every figure passes the visualisation-guide pre-submission checklist.
