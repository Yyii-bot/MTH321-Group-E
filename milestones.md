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
- ||A - A^T||_F/||A||_F = 3.95e-02 (ghost) / 3.93e-02 (first-order) at
  N=64, Frobenius norm, one shared definition with the Algorithm branch
  (their E1 CSV carries the same two values): nonsymmetric, so CG is
  unjustified on paper and sparse direct is the honest choice.
- Algebraic vs. discretisation error separated on N=128 (ghost Neumann):
  spsolve scaled residual 3.4e-15; GMRES at rtol=1e-8 still shifts the
  Linf error by ~80%, at 1e-12 agrees to 1e-9. The tolerance has to be
  ~7 orders below the discretisation error before it stops mattering.

**L-shape results (uniform refinement, n = 8..128, reference P1 FEM):**

- Energy error: p = 0.64..0.66 vs h (theory 2/3), p = 0.28..0.32 vs
  N_dof (theory 1/3). Clean match, stable across every correction below.
- Mesh sanity: total area must equal exactly 3; any partially-cut cell is
  a hard error. The L-shape is represented exactly (cut follows grid lines).
- N_dof sequence 33 / 161 / 705 / 2945 / 12033, identical to the
  Algorithm branch's E5 table (this became true only after the
  corner-node fix below).

> Corrections (honesty rule), final state 2026-10-27: the L2 numbers
> first reported this week (~1.0 order) were quadrature-inflated; the
> second measurement (Duffy tensor Gauss, still a free re-entrant corner
> node) read 1.25-1.32 with absolute values ~3x the certified ones. Two
> real defects were then found and fixed: (i) the re-entrant corner is a
> Dirichlet node -- a free corner node violates the zero trace of u* on
> the cut edges; (ii) low-order quadrature away from the corner is
> replaced by Duffy n=8, and the corner elements by the exact
> radial/adaptive-angular integral ported from the Algorithm branch's
> fem_impl. Final table: L2 1.31-1.34 vs h (textbook 4/3), absolute
> values matching the Algorithm branch's E5 to three+ significant figures
> at identical N_dof. See `results_verified.md` section 4 and audit log
> items 5, 8, 9.

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
   worst case. **Extended 2026-10-25 by `week8_iterative_ladder.py`**
   (results_verified.md section 7): the full lecture ladder (Jacobi, GS,
   SOR, CG, GMRES+ILU) on the same matrix at the same stop criterion
   still loses to direct on every rung -- GMRES+ILU is the closest at
   1.5-1.9x SLOWER, because the ILU setup costs about what the direct LU
   costs at these 2D sizes. SOR's omega_opt formula is verified on its
   own model problem (270 iters at the theory value vs 310 at the best
   grid point, N=64). The surprise: CG converged on the nonsymmetric
   matrix (||A-A^T||_F/||A||_F = 3.95e-02) at ~1.7x the SPD iteration
   count -- unjustified in theory, survives mild nonsymmetry in practice,
   and we say both halves of that sentence.
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
   Both legs run from one script on the same reference FEM, the same
   L-shaped meshes (n = 8..128) and the same error functionals, with the
   corner-exact + Duffy-n8 integrator shared with week6. Final measured
   numbers (2026-10-27): smooth leg (u_s = exp(x+y)) L2 order 2.00 flat,
   energy 1.00; singular leg L2 1.31 -> 1.34, energy 0.64 -> 0.66. The
   theta/rho split shows rho at exactly 2.00 and theta now also ~2.0
   (1.89 -> 2.00), an order of magnitude below ||e|| -- with the corner
   Dirichlet-fixed there is no pollution cap; the earlier "~5/3 cap" was
   the free-corner-node artifact (audit log item 8), retracted in
   writing. The content of this depth comparison is the textbook
   regularity ladder measured with one measuring stick: smooth 2.0/1.0
   vs singular 4/3/2/3 (L2/energy). In results_verified.md section 6.

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
   pointing at the re-measurement history.
3. **The L-shape corrections landed (final pass 2026-10-27).** After the
   corner-node fix and the corner-exact + Duffy-n8 integrator, the L2
   column reads 1.31-1.34 vs h (textbook 4/3) and reproduces the
   Algorithm branch's E5 table to three+ significant figures at identical
   N_dof (161/705/2945/12033). The smooth leg of the depth comparison
   moved to a flat 2.00 with theta ~ 2.0 -- the old "5/3 pollution cap"
   is retracted as a free-corner artifact. Final quoted contrast: smooth
   2.0 (L2) / 1.0 (energy) vs singular 4/3 / 2/3. Recorded in
   `results_verified.md` sections 4, 6, 8 and audit log items 5, 8, 9.
4. **Final residual audit.** Scaled residuals on all Direction-2 tables:
   1e-15..1e-14 on every mesh and closure/BC combination. Raw residuals
   shown once to justify the scaling. The nonsymmetry numbers
   (||A - A^T||_F/||A||_F = 3.95e-02 / 3.93e-02 at N=64) re-confirmed
   and now identical in definition and value to the Algorithm branch's
   E1 CSV.
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
