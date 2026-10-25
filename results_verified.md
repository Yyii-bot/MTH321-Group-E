# Verified results (all real, all reproducible)

Every number on this page was produced by a script in this folder that
runs as-is on the team venv. Reproduce any line with the command in the
last column. Nothing here depends on the team's main solver; where a
number comes from a *reference* implementation it says so.

Generated 2026-10-21; last full rerun and correction pass 2026-10-25.
If a script changes, rerun it and update this file.

## 1. Direction 2 -- mixed boundary data (reference FDM solver)

Source: `week6_mixed_refinement.py` (this is my reference solver, not the
team's main one -- cross-check against the main solver is still pending).
Problem: -Delta u = f on (0,1)^2, u* = exp(x+y), f = -2 exp(x+y).
Reproduce: `python week6_mixed_refinement.py`

### Observed convergence orders (pairwise, 32->64->128->256)

| BC | closure | L2 | Linf | right-edge Linf |
|---|---|---|---|---|
| Neumann | first-order | 1.07 / 1.03 / 1.02 | 1.03 / 1.01 / 1.01 | 1.03 / 1.01 / 1.01 |
| Neumann | ghost | 2.05 / 2.03 / 2.01 | 2.00 / 2.00 / 2.00 | 2.00 / 2.00 / 2.00 |
| Robin | first-order | 1.05 / 1.03 / 1.01 | 1.01 / 1.01 / 1.00 | 1.01 / 1.01 / 1.00 |
| Robin | ghost | 2.06 / 2.03 / 2.02 | 2.00 / 2.00 / 2.00 | 2.00 / 2.00 / 2.00 |

### Supporting quantities (representative values)

| quantity | value | notes |
|---|---|---|
| scaled residual, all meshes/all combos | 1e-15 .. 1e-14 | row scaling D_ii = max_j |A_ij| |
| raw residual growth | ~h^-2 | shown once, justifies the scaling |
| ||A - A^T||/||A|| (inf), N=64 ghost | 1.25e-01 | same for first-order: 1.23e-01 |
| nnz, N=64 ghost | 19906 | |

### Algebraic vs. discretisation error, N=128 ghost Neumann

Reproduce: same script, bottom section.

| solve | wall s | scaled residual | Linf error |
|---|---|---|---|
| spsolve | 0.048 | 3.3e-15 | 1.35e-05 |
| GMRES rtol=1e-4 | 2.5 | 3.6e-05 | 1.27e-01 |
| GMRES rtol=1e-8 | 5.1 | 3.0e-09 | 2.45e-05 |
| GMRES rtol=1e-12 | 8.0 | 3.4e-13 | 1.35e-05 |

## 2. Tolerance sweep (contamination threshold)

Source: `week8_depth_checks.py`. N=128, ghost Neumann.
Reproduce: `python week8_depth_checks.py`

| rtol | scaled residual | Linf error | contaminated |
|---|---|---|---|
| 1e-2 | 5.3e-03 | 5.35e+00 | YES |
| 1e-4 | 2.8e-05 | 1.27e-01 | YES |
| 1e-6 | 3.7e-07 | 1.21e-03 | YES |
| 1e-8 | 3.9e-09 | 2.45e-05 | YES |
| 1e-10 | 3.8e-11 | 1.36e-05 | no |
| 1e-12 | 3.4e-13 | 1.35e-05 | no |

Threshold for this mesh: rtol must be << 1e-8; 1e-10 is the first clean level.

## 3. Direct vs. iterative at matched accuracy (no crossover)

Source: `week8_depth_checks.py` (unpreconditioned GMRES -- deliberate
worst case; any timing claim in the report must name the preconditioner).

| N | spsolve s | GMRES s (matched tol=1e-12) | speedup |
|---|---|---|---|
| 64 | 0.009 | 0.17 | 19x |
| 128 | 0.046 | 9.1 | 198x |
| 256 | 0.274 | 58.5 | 214x |

Finding: no crossover in the tested range; reported as such.

## 4. L-shape, uniform meshes (reference P1 FEM)

Source: `week6_lshape_fem.py` for the energy error and mesh/solve;
**the L2 numbers in this section were re-measured with the Duffy-Gauss
integrator in `week8_smooth_vs_nonsmooth.py` (correction, 2026-10-25)**.
The original 3-point corner-composite L2 values printed by
week6_lshape_fem.py were quadrature-inflated (its own 3-vs-4-level check
was ~9%); the Duffy values carry a 1e-5-level check and an independent
scipy dblquad validation of the integrator. Energy errors are unaffected
(their integrand always converged cleanly) and unchanged below.
Reproduce: `python week6_lshape_fem.py` (energy), then
`python week8_smooth_vs_nonsmooth.py` (L2, both legs, one command).

| n | h | N_dof | L2 error (Duffy, corrected) | energy error | L2 order (vs h) | energy order (vs h) |
|---|---|---|---|---|---|---|
| 8 | 0.25 | 34 | 5.64e-02 | 3.60e-01 | -- | -- |
| 16 | 0.125 | 162 | 2.37e-02 | 2.31e-01 | 1.25 | 0.64 |
| 32 | 0.0625 | 706 | 9.67e-03 | 1.47e-01 | 1.29 | 0.65 |
| 64 | 0.03125 | 2946 | 3.90e-03 | 9.31e-02 | 1.31 | 0.66 |
| 128 | 0.015625 | 12034 | 1.57e-03 | 5.88e-02 | 1.32 | 0.66 |

| quantity | value |
|---|---|
| energy order vs N_dof (all levels) | 0.28 .. 0.33 (theory 1/3) |
| L2 quadrature self-check (corner levels 3 vs 4, Duffy) | 1e-6 .. 1e-5 |
| mesh area sanity | exactly 3 (= 4 - 1) at every level |

Interpretation (measured): energy matches the regularity-aware theory
(2/3 vs h, 1/3 vs N_dof). L2 sits at 1.25-1.32 and rises with
refinement -- consistent with the h^{4/3-eps} guarantee and visibly
approaching 4/3, unlike the quadrature-contaminated first measurement
(~1.0). The uniform-mesh energy reference slope vs N_dof is ~0.33 --
the number the adaptive sequence has to beat.

## 5. Week 5 analytic checks (all green)
Reproduce: `python week5_verify_analytic.py`, `python week5_boundary_rows.py`

- f = -2 exp(x+y) confirmed by stencil at 4 sample points (rtol 1e-5).
- q = exp(1+y) > 0 -- outward normal sign consistent.
- r = 3 exp(1+y) = q + 2u* -- Robin data consistent.
- L-shape: phi in [0, 3pi/2] at 5 test points incl. quadrant IV; u* = 0
  on both rays (machine zero); outer boundary data nonzero; Delta u* ~ 0
  (stencil, r >= 0.3); u* >= 0 over 1500 random samples.
- First-order rows: residual ratio 1.99 per halving (O(h)).
- Ghost rows: residual ratio ~1.98 per halving (O(h) locally -- the
  O(h^2) claim is about the global error, verified in section 1).
- Corner ownership assertion: passes on a correct tag grid.

## 6. Depth comparison 4: smooth vs. nonsmooth (same solver, meshes, norms)

Source: `week8_smooth_vs_nonsmooth.py` -- one command reproduces both
legs. L2 by Duffy-transform tensor Gauss (integrator validated against
`scipy.integrate.dblquad` to machine precision; singular leg
corner-refined with a 3-vs-4-level check of 1e-6..1e-5).
Reproduce: `python week8_smooth_vs_nonsmooth.py`

Same reference P1 FEM, same L-shaped meshes (n = 8..128), same error
functionals.  Smooth leg: u_s = exp(x+y), f = -2 exp(x+y).  Nonsmooth
leg: u* = r^{2/3} sin(2 phi/3), f = 0.

| leg | L2 order (pairwise) | energy order (pairwise) |
|---|---|---|
| smooth | 1.70 / 1.67 / 1.66 / 1.65 | 1.00 / 1.00 / 1.00 / 1.00 |
| singular | 1.25 / 1.29 / 1.31 / 1.32 | 0.63 / 0.65 / 0.66 / 0.66 |

Decomposition of the smooth-leg error, e = theta - rho with
theta = u_h - u_I (pollution) and rho = u_e - u_I (interpolation):

| n | ||e||_L2 | ||theta|| (exact) | ||rho|| (Duffy) |
|---|---|---|---|
| 8 | 3.15e-02 | 2.90e-02 | 2.92e-02 |
| 16 | 9.69e-03 | 9.91e-03 | 7.33e-03 |
| 32 | 3.04e-03 | 3.24e-03 | 1.84e-03 |
| 64 | 9.61e-04 | 1.04e-03 | 4.59e-04 |
| 128 | 3.05e-04 | 3.29e-04 | 1.15e-04 |

||rho|| decays at exactly 2.00 per halving (the textbook interpolation
rate, measured).  ||theta|| decays at ~1.55-1.66 and dominates -- so on
the L-domain even the SMOOTH problem's L2 rate is pollution-capped near
5/3 (= h^{1+s}, s -> 2/3: the dual-problem regularity limit), and the
singular leg degrades one more step toward 4/3.  The ladder
2.0 (interpolation only) -> 5/3 (smooth FEM) -> 4/3 (singular FEM) is
the quantified content of this depth comparison; the energy column
(1.0 -> 2/3) is its textbook counterpart.

## 7. Iterative method ladder (Jacobi / GS / SOR / CG / GMRES+ILU)

Source: `week8_iterative_ladder.py` -- the course lecture's method ladder
run on the SAME discrete system as section 1 (same assembler, ghost
Neumann, same scaled-residual stop 1e-10). Extends the Week-8
"direct vs. iterative" depth item beyond unpreconditioned GMRES.
Reproduce: `python week8_iterative_ladder.py` (about 2 min: Jacobi/GS
at N=128 need 85k/43k iterations). Wall times are single runs and
vary a few percent between runs; iteration counts are stable.

### 7a. The ladder, mixed BC (ghost Neumann), u* = exp(x+y)

| N | method | wall s | iters | Linf error at stop |
|---|---|---|---|---|
| 64 | spsolve (direct) | 0.010 | -- | 5.41e-05 |
| 64 | Jacobi | 1.475 | 23250 | 5.42e-05 |
| 64 | Gauss-Seidel | 3.431 | 11620 | 5.42e-05 |
| 64 | SOR w=1.900 | 0.186 | 590 | 5.41e-05 |
| 64 | CG | 0.017 | 425 | 5.41e-05 |
| 64 | GMRES+ILU | 0.015 | 5 | 5.41e-05 |
| 128 | spsolve (direct) | 0.047 | -- | 1.35e-05 |
| 128 | Jacobi | 38.387 | 85630 | 1.40e-05 |
| 128 | Gauss-Seidel | 57.081 | 42810 | 1.40e-05 |
| 128 | SOR w=1.900 | 3.456 | 2480 | 1.36e-05 |
| 128 | CG | 0.512 | 748 | 1.35e-05 |
| 128 | GMRES+ILU | 0.088 | 6 | 1.35e-05 |
| 256 | spsolve (direct) | 0.272 | -- | 3.38e-06 |
| 256 | SOR w=1.900 | 55.824 | 9360 | 3.60e-06 |
| 256 | CG | 1.322 | 1351 | 3.40e-06 |
| 256 | GMRES+ILU | 0.475 | 9 | 3.38e-06 |

Jacobi/GS at N=256 not run (O(N^2) iterations; the 64->128 counts show
the growth explicitly). SOR omega for the mixed system swept at N=64:
best 1.900 (590 iters).

Measured findings:
1. **No crossover on any rung, including GMRES+ILU** (closest: 1.5-1.9x
   SLOWER than direct at every N). The ILU setup costs about what the
   direct LU costs at these 2D sizes, and the preconditioned GMRES then
   needs only 5-9 iterations -- the preconditioner is excellent, the
   total still loses. Extrapolation beyond N=256 is stated, not tested.
2. **Classical methods scale as theory says**: Jacobi iterations grow
   23250 -> 85630 per halving (x3.7 ~ N^2 log); GS is consistently a
   factor ~2 better in iterations but slower per iteration here.
3. **CG on the nonsymmetric matrix converged** (425/748/1351 iters,
   scaled res ~1e-11, clean errors) even though ||A-A^T||/||A|| = 0.12
   makes it unjustified in theory. Theory-vs-practice nuance, measured
   not assumed: mild nonsymmetry slowed CG by ~1.7x vs the SPD
   Dirichlet matrix (748 vs 436 iters at N=128) but did not break it.

### 7b. SOR omega_opt on the lecture's model problem (pure-Dirichlet square)

| N | h | theory omega_opt = 2/(1+sin(pi h)) | measured best omega | iters at best | iters at theory value |
|---|---|---|---|---|---|
| 64 | 0.01562 | 1.906455 | 1.900 (0.025-spaced grid) | 310 | 270 |
| 128 | 0.00781 | 1.952093 | (sweep at 64 only) | -- | 530 |

The theory value is verified on its own model problem: 270 iterations at
omega = 1.906455, better than any grid point, and omega > omega_opt
diverges (aborted by a 100x residual-growth detector). On the
nonsymmetric mixed system the same formula is a starting guess only --
the sweep found 1.900 there too, but the lecture's own caveat (the
formula assumes the standard model) is stated in the report.

### 7c. Two theory spot-checks

| check | result |
|---|---|
| CG vs CG+Jacobi-precond, Dirichlet square N=64 | 221 vs 221 iterations -- identical, confirming D^-1A = (h^2/4)A is a constant scaling, so the Jacobi preconditioner buys nothing here (lecture point 13, verified) |
| CG on SPD Dirichlet matrix, N=128 | 436 iters, clean |
| CG on nonsymmetric mixed matrix, N=128 | 748 iters, scaled res 5.7e-11, Linf error 1.35e-05 -- converged, see 7a finding 3 |

## Not yet real (do not quote as results)

- Adaptive effectivity numbers, adaptive-vs-uniform slopes -- awaiting
  the adaptive driver's CSV (see needs_from_team.md).
- FDM-vs-FEM comparison numbers -- awaiting both main solvers.
- Cross-check of these reference numbers against the team's main solver.
- Everything in week9_final_audit.md is a checklist, not data.
