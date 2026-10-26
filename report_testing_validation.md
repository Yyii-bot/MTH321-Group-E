# Testing & Validation Report -- Project 2

**Numerical Solution of Elliptic BVPs, Weeks 5-9**
Team direction: Direction 2 (mixed Dirichlet/Neumann/Robin data, FDM, two
boundary closures) + Open Direction 1 (L-shaped domain, regularity-aware
interpretation, adaptive refinement); all four depth comparisons.
Seat: **Testing & Validation** (Rong Yin)

This report is the verification and validation evidence behind every number
the team report quotes. Each number below can be reproduced with the command
given in `results_verified.md`, on a clean machine with numpy / scipy /
matplotlib installed. Nothing here depends on the team's main solver; where a
number comes from a reference implementation, it says so, and every
cross-check against the main solver is recorded with its result rather than
assumed.

Sign convention throughout: we solve **-Delta u = f** (the test bank's
convention). For u* = exp(x+y) that gives f = -2 exp(x+y). The problem pack
prints "Delta u = f" once, but its own data matches -Delta u = f.

---

## 1. Verification of the manufactured data (Week 5, before any solver)

Every piece of problem data was re-derived independently of my own calculus:

- f = -2 exp(x+y) confirmed by applying the 5-point stencil to u* at four
  sample points (agreement to rtol 1e-5). Script: `week5_verify_analytic.py`.
- Boundary data signs checked by sampling: q = exp(1+y) > 0 on x=1 fixes the
  outward normal direction; r = 3 exp(1+y) = q + 2u* confirms the Robin
  data is consistent with the chosen u*.
- All four right-edge discrete rows (Neumann/Robin x first-order/ghost)
  were verified by substituting sampled values of u* into the row and
  watching the residual shrink like the derived local order. Script:
  `week5_boundary_rows.py`. First-order rows: residual ratio ~1.99 per
  halving (O(h) local). Ghost rows: ratio ~1.98 (O(h) local -- the O(h^2)
  claim for ghost is about the *global* error and is verified on the
  solution in Section 2, not on the row).
- The L-shape branch convention (phi in [0, 3pi/2]) was unit-tested
  including the quadrant-IV points where atan2 alone gets it wrong;
  u* = r^{2/3} sin(2phi/3) vanishes on both rays (machine zero), is
  nonnegative over 1505 random samples, is harmonic (stencil residual
  ~1e-7 for r >= 0.3), and the outer boundary data is nonzero.
- A corner-ownership assertion in the assembler guarantees a boundary row
  can never overwrite a Dirichlet corner node.

One bug was caught at this stage, in my own test code: the Robin row test
double-counted q (subtracting both q and r when r = q + 2u* already
contains it). The residual refused to shrink and exposed it.

## 2. Direction 2: mixed boundary data (reference FDM solver, Week 6)

A deliberately small reference FDM solver (`week6_mixed_refinement.py`),
independent of the main implementation, produced four-mesh refinement
tables (N = 32..256) for all four closure/BC combinations. The exact
solution u* = exp(x+y) is known, so the error is measured, not assumed.

**Observed orders (pairwise), from `results_verified.md` section 1:**

| BC | closure | L2 | Linf | right-edge Linf |
|---|---|---|---|---|
| Neumann | first-order | 1.07 / 1.03 / 1.02 | 1.03 / 1.01 / 1.01 | 1.03 / 1.01 / 1.01 |
| Neumann | ghost | 2.05 / 2.03 / 2.01 | 2.00 / 2.00 / 2.00 | 2.00 / 2.00 / 2.00 |
| Robin | first-order | 1.05 / 1.03 / 1.01 | 1.01 / 1.01 / 1.00 | 1.01 / 1.01 / 1.00 |
| Robin | ghost | 2.06 / 2.03 / 2.02 | 2.00 / 2.00 / 2.00 | 2.00 / 2.00 / 2.00 |

Two findings matter beyond the rates themselves:

1. **A small algebraic residual alone does not validate the boundary
   discretisation.** The row-substitution check caught a Dirichlet
   sign error (columns moved to the rhs with the wrong sign) while the
   scaled residuals were machine-zero (~1e-15) on every mesh. Residuals
   measure consistency of the *assembly*, not correctness of the *operator*.

2. **Row scaling is mandatory when PDE and boundary rows coexist.**
   Raw residuals grow like h^-2 while scaled residuals (D_ii = max_j
   |A_ij|) sit at 1e-15..1e-14 on every mesh and combination. The raw
   residual is shown once in the script output to justify the scaling.

The matrix is nonsymmetric: ||A - A^T||_F/||A||_F = 3.95e-02 (ghost) and
3.93e-02 (first-order) at N=64, in the Frobenius norm and with exactly the
definition the Algorithm branch uses (their E1 CSV carries the same two
values). That rules CG out on paper; the reference tables use sparse direct
(spsolve), with CG's actual behaviour studied honestly in Depth 2b below.

## 3. L-shaped domain: regularity-aware interpretation (Weeks 6, 8)

Reference P1 FEM on exactly-represented L-shape meshes
(`week6_lshape_fem.py`; the same mesh/solve code the smooth-vs-nonsmooth
study imports). The mesh sanity check (total area exactly 3 = 4 - 1) guards
against partially-cut corner cells. The re-entrant corner (0,0) is a
Dirichlet node: it lies on the two cut edges, where u* has zero trace.

**Measured on five uniform levels (n = 8..128), from
`results_verified.md` section 4.** Both error columns come from one
measuring stick: elements touching the corner use the exact radial /
adaptive-angular integral ported from the Algorithm branch's fem_impl;
all other elements use Duffy tensor Gauss at n=8 (their intorder=8) for
L2 and the 3-point edge-midpoint rule for the energy seminorm.

| n | h | N_dof | L2 error | energy error | L2 order | energy order |
|---|---|---|---|---|---|---|
| 8 | 0.25 | 33 | 1.8812e-02 | 1.9289e-01 | -- | -- |
| 16 | 0.125 | 161 | 7.5922e-03 | 1.2402e-01 | 1.31 | 0.64 |
| 32 | 0.0625 | 705 | 3.0235e-03 | 7.9193e-02 | 1.33 | 0.65 |
| 64 | 0.03125 | 2945 | 1.1972e-03 | 5.0324e-02 | 1.34 | 0.65 |
| 128 | 0.015625 | 12033 | 4.7313e-04 | 3.1878e-02 | 1.34 | 0.66 |

- The **energy rate matches regularity-aware theory**: 0.64-0.66 vs h
  (theory 2/3, driven by H^{1+2/3-eps} regularity at the re-entrant
  corner) and 0.28..0.32 vs N_dof (theory 1/3). This uniform-mesh slope
  vs N_dof is the reference the adaptive sequence must beat.
- The **L2 rate is 1.31-1.34 vs h** -- the textbook 4/3 rate within a
  third of an order at every level, as the Aubin-Nitsche h^{4/3-eps}
  guarantee suggests it should be once the measurement is clean.
- **Cross-validation against the team's main solver (their E5 table):**
  identical N_dof (161/705/2945/12033) and L2 errors agreeing to three
  or more significant figures on every shared level (their 7.59e-3 /
  3.02e-3 / 1.20e-3 / 4.73e-4 vs our 7.5922e-3 / 3.0235e-3 / 1.1972e-3 /
  4.7313e-4). Two independent implementations agreeing on a singular
  problem to this precision is the strongest single piece of evidence in
  this package.

**Two corrections, recorded here because the honesty rule applies to the
validator first.** An earlier version of this table reported L2 ~1.0 from
a 3-point corner-composite quadrature; a second version reported 1.25-1.32
with Duffy tensor Gauss but a *free* re-entrant corner node. Both defects
were real: (i) low-order quadrature of e^2 has error O(h^3), the same
order as the h^4 signal (grad e is O(1) near the corner, so D^3(e^2) is
O(1)), so the integration error grows relative to the signal exactly when
the mesh refines -- two such rules can agree with each other to <1e-3 and
both be wrong; (ii) a free corner node violates the zero trace of u* on
the cut edges, leaves u_h(corner) ~ 0.29 on coarse meshes, and quietly
degrades both legs' rates. The final protocol removes both: corner-exact
integral (ported, certified by its own adaptive error estimate at 1e-11)
plus Duffy n=8 elsewhere, with the corner node Dirichlet-fixed. The
self-check column (Duffy n=8 vs n=6) is now 1e-7, and the absolute values
reproduce the Algorithm branch's E5 table. All energy numbers were never in
doubt across these corrections (their integrand converged cleanly) and
only moved within the last digit. Full record in `week8_ai_audit_log.md`,
items 5, 8, 9.

The corner region (r < 0.25) carries 91-92% of the squared L2 error at
every level -- the quantitative motivation for the team's adaptive
refinement around the re-entrant corner.

## 4. Depth comparisons (Weeks 7-8)

The protocol (frozen in `week7_clinic_readiness.md`): every comparison
shares problem, norm, refinement protocol and solver tolerance, and a
gate function (`week8_depth_checks.py`) refuses to compare runs that
differ in any of these -- it caught its first mismatch in the demo.

**Depth 1 -- tolerance contamination.** GMRES sweep at N=128, ghost
Neumann, against the spsolve reference (Linf error 1.35e-5):

| rtol | scaled residual | Linf error | contaminated |
|---|---|---|---|
| 1e-2 | 5.3e-3 | 5.35e+00 | YES |
| 1e-4 | 2.8e-5 | 1.27e-01 | YES |
| 1e-6 | 3.7e-7 | 1.21e-03 | YES |
| 1e-8 | 3.9e-9 | 2.45e-05 | YES |
| 1e-10 | 3.8e-11 | 1.36e-05 | no |
| 1e-12 | 3.4e-13 | 1.35e-05 | no |

The residual at rtol=1e-8 looks excellent (3.9e-9) and the error is
still 80% too large. Stopping criterion used everywhere after this:
scaled residual at rtol <= 1e-10, never iterate difference.

**Depth 2 -- direct vs. iterative at matched accuracy.** Unpreconditioned
GMRES at rtol=1e-12 vs spsolve, N = 64..256: direct is faster by 19x-214x
depending on run and mesh (see `results_verified.md` section 3 for the
recorded run; wall times vary per run, the finding does not). **No
crossover observed in the tested range**, reported as such. Any timing
claim with a preconditioner must name the preconditioner; ours has none,
a deliberate worst case.

**Depth 2b -- the full iterative ladder (lecture methods) on the same
system.** `week8_iterative_ladder.py` (`results_verified.md` section 7)
runs Jacobi, Gauss-Seidel, SOR, CG and GMRES+ILU on the same ghost-Neumann
matrix, same stop (scaled residual <= 1e-10), N = 64..256. Every rung
loses to the direct solve. The classical rungs scale exactly as the
lecture predicts -- Jacobi needs 23250 -> 85630 iterations per halving
(O(N^2 log)), GS half of that, SOR with a swept omega (1.900 on this
system) 590 -> 9360 -- but even the best rung, GMRES+ILU, is 1.5-1.9x
SLOWER than direct at every N: the ILU setup costs about what the direct
LU costs at these 2D sizes, and the preconditioned GMRES then needs only
5-9 iterations. The no-crossover finding survives preconditioning.
Two theory checks from the lecture are verified with real numbers:
SOR's omega_opt = 2/(1+sin(pi h)) on its own model problem (pure
Dirichlet square) takes 270 iterations at N=64 -- fewer than any value
on the 0.025-spaced sweep grid (best grid point 1.900: 310), and any
omega above it diverges; and CG with the Jacobi preconditioner converges
in exactly the same 221 iterations as plain CG on the 5-point Dirichlet
matrix, confirming D^-1A is a constant scaling there so the Jacobi
preconditioner buys nothing. One honest surprise: CG *converged* on our
nonsymmetric mixed matrix (||A-A^T||_F/||A||_F = 3.95e-02), at 748 vs
436 iterations on the SPD Dirichlet matrix at N=128, with a clean error.
The theory says CG requires symmetry; mild nonsymmetry slowed it by
~1.7x but did not break it. We report both halves of that sentence.

**Depth 3 -- FDM vs. FEM.** The matched-settings gate is written, tested,
and passed its self-check; the real comparison awaits both main solvers
and is reported as pending, not faked.

**Depth 4 -- smooth vs. nonsmooth (same solver, same meshes, same
norms).** `week8_smooth_vs_nonsmooth.py` runs both legs on the reference
P1 FEM and the same L-shaped meshes (n = 8..128), with the section-3
measuring stick:

| leg | L2 order (pairwise) | energy order (pairwise) |
|---|---|---|
| smooth (u_s = exp(x+y)) | 2.00 / 2.00 / 2.00 / 2.00 | 0.99 / 1.00 / 1.00 / 1.00 |
| singular (r^{2/3} sin(2phi/3)) | 1.31 / 1.33 / 1.34 / 1.34 | 0.64 / 0.65 / 0.65 / 0.66 |

The decomposition e = theta - rho (theta = u_h - u_I discrete pollution,
computed exactly; rho = u_I - u exact interpolation error) turns this
from a table into an explanation: rho decays at exactly 2.00 per halving,
and theta now also decays like h^2 (1.89 -> 2.00) while staying more than
an order of magnitude below ||e|| -- with correct boundary data, even on
the L-domain the smooth problem keeps its full O(h^2) L2 rate. (An
earlier draft reported a "pollution cap near 5/3" for the smooth leg;
that was the free corner node artifact, retracted in writing -- audit log
item 8.) The singular leg degrades exactly one regularity step: L2 to
4/3, energy to 2/3. The quantified content of this depth comparison is
the textbook regularity ladder measured on one solver with one measuring
stick: **smooth 2.0 / 1.0 vs singular 4/3 / 2/3** (L2 / energy).

## 5. Reproducibility and audit

- Every number in this report traces to `results_verified.md`, which
  carries the rerun command for each table. All scripts were rerun from a
  clean shell during the Week 9 audit (2026-10-25) and reproduced; the
  2026-10-27 correction pass reran the two FEM scripts again and every
  quoted value above is from that final run.
- The five-gate AI audit (transcription / structure / linearisation /
  independent reproduction / interrogation) for every AI-generated solver
  piece, including the caught bugs and the two-stage L2 correction, is in
  `week8_ai_audit_log.md` and goes in the report appendix.
- Pre-computation worksheets: the ghost-closure O(h)-local/O(h^2)-global
  prediction and the L-shape rate predictions were written down before
  the runs; both predictions were confirmed.

## 6. What is still open (stated, not hidden)

- Adaptive effectivity numbers -- pending the adaptive driver's CSV; the
  effectivity script (`week7_effectivity.py`) is ready.
- FDM vs. FEM comparison -- pending both main solvers.
- These items are listed in `results_verified.md` under "Not yet real"
  and must not be quoted as results.

---

Rong Yin, Testing & Validation seat, 2026-10-27
