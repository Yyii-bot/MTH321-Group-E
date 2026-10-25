# Weekly milestones -- Testing & Validation (Rong Yin)

Project 2, Weeks 5-9. One file, one section per week.

## Week 5

Ran every analytic check *before* any solver existed, so that when numbers
go wrong later we know whether to blame the maths or the code.

- `week5_verify_analytic.py` -- all green:
  - Direction 2: f = -2 exp(x+y) confirmed against a fine 5-point stencil
    (not by hand-waving derivatives); boundary data sampled on all four
    edges; q = exp(1+y) > 0 confirms the outward normal on x=1 points
    the right way; Robin data r = 3 exp(1+y) satisfies r = q + 2u*.
  - L-shape: branch convention phi = theta - pi/2 in [0, 3pi/2] unit-tested
    at five points including the nasty quadrant IV; u* = 0 on both rays
    (machine zero); outer boundary data nonzero (so the problem isn't
    trivially u = 0); Delta u* = 0 verified by stencil at points with
    r >= 0.3; u* >= 0 over 1500 random domain samples.
- `week5_boundary_rows.py` -- the four right-edge rows (first-order and
  ghost, Neumann and Robin) verified by substituting sampled u* values:
  - first-order rows: residual shrinks like h (ratio 1.99 per halving),
    as a first-order closure must.
  - ghost rows: residual shrinks like h *locally* -- this is expected
    (LeVeque FDML), the O(h^2) payoff must show up in the solution error,
    which is exactly what the Week 6 tables will test.
  - corner ownership helper: asserts the four Dirichlet corners are never
    overwritten by a boundary row. Runs against the real tag grid.

Flag for the team: the problem pack prints the PDE as "Delta u = f" but
then gives f = -2 exp(x+y), which is consistent with -Delta u = f (the
test bank's convention). We adopt -Delta u = f. Worth a line in the
report so the clinic doesn't trip on it.

*Plan for Week 6:* four-mesh refinement tables for all four closure/BC
combinations, with boundary-vs-interior errors split out, scaled
residuals, and the matrix asymmetry numbers. Same for the L-shape on
uniform meshes. *Blocker: none.*

## Week 6

Produced the first full refinement tables, and (true to the role) found
two real bugs in the reference solver before any number was trusted:

1. Dirichlet columns were moved to the rhs with the wrong sign in the
   assembly -- the corner-adjacent rows were off by exactly 2x the
   Dirichlet term. Caught by plugging u* into the assembled rows
   (A u* - b per row), fixed, only then were tables produced.
2. Robin first-order test row double-counted q (r already contains q + 2u*).
   Caught because the residual refused to shrink under refinement.

**Direction 2 results (u* = exp(x+y), N = 32..256):**

| closure | BC | L2 order | Linf order | right-edge Linf order |
|---|---|---|---|---|
| first-order | Neumann | 1.07 / 1.03 / 1.02 | 1.03 / 1.01 / 1.01 | 1.03 / 1.01 / 1.01 |
| ghost | Neumann | 2.05 / 2.03 / 2.01 | 2.00 / 2.00 / 2.00 | 2.00 / 2.00 / 2.00 |
| first-order | Robin | 1.05 / 1.03 / 1.01 | 1.01 / 1.01 / 1.00 | 1.01 / 1.01 / 1.00 |
| ghost | Robin | 2.06 / 2.03 / 2.02 | 2.00 / 2.00 / 2.00 | 2.00 / 2.00 / 2.00 |

(pairs are per mesh doubling 32->64->128->256; fitted only on the levels shown)

Take-aways for the report:

- The closure, not the PDE stencil, sets the global rate: first-order
  closure caps everything at ~1, ghost recovers the full 2 even at the
  boundary nodes themselves. So the boundary closure DOES limit the
  global error for the first-order variant.
- Scaled residuals stay at 1e-15..1e-14 on every mesh (row scaling:
  D_ii = max_j |A_ij|; raw residuals grow like h^-2 and would mislead).
- ||A - A^T||/||A|| ~ 0.12-0.125 for both closures: nonsymmetric, so CG
  is unjustified and sparse direct is the honest choice.
- Algebraic vs. discretisation error separated on N=128 (ghost Neumann):
  spsolve scaled residual 3.4e-15; GMRES at rtol=1e-8 still shifts the
  Linf error by ~80%, at 1e-12 agrees to 1e-9. The tolerance has to be
  ~7 orders below the discretisation error before it stops mattering.

**L-shape results (uniform refinement, n = 8..128, reference P1 FEM):**

- Energy error: p = 0.64..0.66 vs h (theory 2/3), p = 0.28..0.33 vs
  N_dof (theory 1/3). Clean match.
- Mesh sanity: total area must equal exactly 3; any partially-cut cell is
  a hard error. The L-shape is represented exactly (cut follows grid lines).

> Correction (2026-10-25), recorded per the honesty rule: the L2 numbers
> first reported this week (~1.0 order) were later found to be
> quadrature-inflated. The L2 signal (~h^4) is smaller than the O(h^3)
> quadrature error of any low-order rule applied to e^2 (grad e is O(1)
> pointwise, so D^3(e^2) = O(1)). Re-measured with the Duffy-transform
> tensor Gauss integrator (validated against scipy.integrate.dblquad;
> self-check 1e-5), the L-shape L2 errors are ~3x smaller and the order
> is 1.25-1.32 rising toward 4/3, not ~1.0. See `results_verified.md`
> section 4 and `week8_smooth_vs_nonsmooth.py`. Energy numbers were
> never affected and are unchanged.

*Plan for Week 7:* clinic readiness sheet; effectivity-index design for
the adaptive study; freeze the matched-settings protocol for the four
depth comparisons. *Blocker: none.*

## Week 7

- Assembled the clinic readiness sheet (`week7_clinic_readiness.md`):
  the discretisation/norm/tolerance triple is frozen, the evidence list
  is attached, known weaknesses are stated out loud, and the adaptive
  protocol (marking rule, effectivity measurement, matched settings
  across uniform and adaptive sequences) is written down before any
  adaptive numbers exist.
- Designed the effectivity-index study (eta / true energy error on every
  adaptive level) that our adaptive claim rests on; the analysis script
  `week7_effectivity.py` is ready and waiting for the adaptive driver's
  CSV.
- Froze the matched-settings protocol for the depth comparisons:
  problem, norm, refinement protocol and solver tolerance must be
  identical between any two runs being compared.

*Plan for Week 8:* run the tolerance sweep as contamination evidence,
the direct-vs-iterative cost comparison, the five-gate AI audit, and the
reproduction pass. *Blocker: none.*

## Week 8

Closed the verification evidence for the four depth comparisons and the
reproduction check.

1. **Matched-settings gate (`week8_depth_checks.py`).** Every comparison
   in the report must share problem, norm, refinement protocol and solver
   tolerance. The gate function refuses to compare runs that differ in any
   of these -- it caught its first mismatch in the demo (L2 vs Linf),
   which is the point.
2. **Tolerance sweep (contamination evidence).** N=128, ghost Neumann,
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
3. **Direct vs. iterative at matched accuracy (no crossover in range).**
   Unpreconditioned GMRES at matched tol never gets close for N <= 256
   (direct faster by 19x-214x across runs) -- reported as "no crossover
   observed", which the brief explicitly accepts as a finding. Any timing
   claim with a preconditioner must name it; ours has none, a deliberate
   worst case.
4. **FDM vs. FEM gate.** The matched-settings gate functions passed their
   demo self-check, but no real FEM numbers have gone through it yet --
   stated honestly, the gate is *ready*, the comparison is *pending* both
   main solvers.
5. **AI audit log.** `week8_ai_audit_log.md` records the five-gate audit
   for every AI-generated solver piece, including the bugs the audit
   caught (Dirichlet RHS sign; Robin row double-counting; L-corner mesh
   cell) and the pre-computation worksheets. Goes in the report appendix.
6. **Reproduction check.** Fresh shell, plain `python` from the venv:
   `week5_*.py`, `week6_*.py`, `week8_depth_checks.py` rerun end-to-end;
   every table number reproduced to the last printed digit. Consolidated
   into `results_verified.md` (only real numbers, each with its rerun
   command).
7. **Smooth-vs-nonsmooth, final leg (`week8_smooth_vs_nonsmooth.py`).**
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

*Plan for Week 9:* final audit -- trace every report number back to a
manufactured check or a refinement table; final residual audit; finish
ICS; verify the zip runs on a clean machine. *Blocker: none.*

## Week 9

Final audit week. Everything on this list was executed, not just planned.

1. **Clean-shell reproduction pass.** Reran every script from a plain
   shell on the venv: `week5_verify_analytic.py`,
   `week5_boundary_rows.py`, `week6_mixed_refinement.py`,
   `week6_lshape_fem.py`, `week8_depth_checks.py`,
   `week8_smooth_vs_nonsmooth.py`. Every table number reproduced;
   wall-clock timings moved a little between runs (they are timings),
   orders and error values did not.
2. **Number-tracing pass.** Every number in the report traces:
   report -> `results_verified.md` -> table -> script -> command. The
   Ctrl-F pass found one stale sentence in the audit checklist (the
   pre-correction "~1.0" L2 rate); it now carries the correction note
   pointing at the Duffy re-measurement.
3. **The smooth-vs-nonsmooth correction landed.** The L2 column for the
   L-shape was re-measured with the Duffy + tensor-Gauss integrator and
   the theta/rho decomposition is in. Final quoted ladder: 2.0
   (interpolation) -> 5/3 (smooth FEM, measured 1.65-1.70) -> 4/3
   (singular FEM, measured 1.25-1.32 and rising). Recorded in
   `results_verified.md` section 6 and in the audit log as items 5-6.
4. **Final residual audit.** Scaled residuals on all Direction-2 tables:
   1e-15..1e-14 on every mesh and closure/BC combination. Raw residuals
   shown once to justify the scaling. The nonsymmetry number
   (||A - A^T||/||A|| ~ 0.12) re-confirmed at N=64.
5. **Package freeze for submission.** `report_testing_validation.md`
   written as the verification evidence behind the report numbers;
   `week9_final_audit.md` checkboxes resolved (done items checked, pending
   items explicitly marked as team dependencies, not silently dropped);
   ICS finalised; README updated so a stranger can rerun everything.
6. **Repo handover.** Package pushed to the group repo's
   testing-validation branch; previous project files preserved under
   their original paths, its root README under `archive/`.

*Remaining days:* hand the package to the report owner; stay on call for
the cross-check the day the main solver interface lands (the diff script
is a 30-minute job, already scoped in needs_from_team.md). Rehearse the
clinic Q&A list at the bottom of the audit checklist.

*Blocker:* none on my side. The three open items (main-solver
cross-check, adaptive CSV, FDM-vs-FEM) all sit with teammates and are
tracked in needs_from_team.md with fallbacks stated.
