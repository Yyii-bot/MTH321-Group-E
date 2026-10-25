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
number comes from a reference implementation, it says so, and the
cross-check against the main solver is listed as pending rather than
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

The matrix is nonsymmetric: ||A - A^T||/||A|| ~ 0.12 for both closures
at N=64. That rules CG out; the reference tables use sparse direct
(spsolve), with GMRES studied separately below.

## 3. L-shaped domain: regularity-aware interpretation (Weeks 6, 8)

Reference P1 FEM on exactly-represented L-shape meshes
(`week6_lshape_fem.py`; reuses the same mesh/solve code as the
smooth-vs-nonsmooth study). The mesh sanity check (total area exactly 3 =
4 - 1) guards against partially-cut corner cells.

**Measured on five uniform levels (n = 8..128), from
`results_verified.md` section 4:**

| n | h | N_dof | L2 error (Duffy, corrected) | energy error | L2 order | energy order |
|---|---|---|---|---|---|---|
| 8 | 0.25 | 34 | 5.64e-02 | 3.60e-01 | -- | -- |
| 16 | 0.125 | 162 | 2.37e-02 | 2.31e-01 | 1.25 | 0.64 |
| 32 | 0.0625 | 706 | 9.67e-03 | 1.47e-01 | 1.29 | 0.65 |
| 64 | 0.03125 | 2946 | 3.90e-03 | 9.31e-02 | 1.31 | 0.66 |
| 128 | 0.015625 | 12034 | 1.57e-03 | 5.88e-02 | 1.32 | 0.66 |

- The **energy rate matches regularity-aware theory**: 0.66 vs h (theory
  2/3, driven by H^{1+2/3-eps} regularity at the re-entrant corner) and
  0.28..0.33 vs N_dof (theory 1/3). This uniform-mesh slope vs N_dof is
  the reference the adaptive sequence must beat.
- The **L2 rate is 1.25-1.32 and rises with refinement**, consistent
  with the h^{4/3-eps} guarantee and visibly approaching 4/3.

**Correction, recorded here because the honesty rule applies to the
validator first.** An earlier version of this table reported L2 ~1.0,
from a 3-point corner-composite quadrature. That integrator was
quadrature-inflated: two low-order rules agreed with each other to <1e-3
on the same mesh (so they *looked* mutually validated) yet both drifted
~50% between refinement levels, because the error signal scales like h^4
while the quadrature error of e^2 scales like h^3 (grad e is O(1) near
the corner, so D^3(e^2) is O(1)) -- the integration error grows relative
to the signal exactly when the mesh refines. The column was re-measured
with a Duffy-transform tensor-Gauss integrator carrying a 3-vs-4-level
self-check of 1e-6..1e-5 and independently validated against
`scipy.integrate.dblquad` on a known integral. All energy numbers were
never in doubt (their integrand converged cleanly under the same check)
and are unchanged. Full record in `week8_ai_audit_log.md`, items 5-6.

The corner region carries a disproportionate share of the error, which is
the quantitative motivation for the team's adaptive refinement around the
re-entrant corner.

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

**Depth 3 -- FDM vs. FEM.** The matched-settings gate is written, tested,
and passed its self-check; the real comparison awaits both main solvers
and is reported as pending, not faked.

**Depth 4 -- smooth vs. nonsmooth (same solver, same meshes, same
norms).** `week8_smooth_vs_nonsmooth.py` runs both legs on the reference
P1 FEM and the same L-shaped meshes (n = 8..128):

| leg | L2 order (pairwise) | energy order (pairwise) |
|---|---|---|
| smooth (u_s = exp(x+y)) | 1.70 / 1.67 / 1.66 / 1.65 | 1.00 / 1.00 / 1.00 / 1.00 |
| singular (r^{2/3} sin(2phi/3)) | 1.25 / 1.29 / 1.31 / 1.32 | 0.63 / 0.65 / 0.66 / 0.66 |

The decomposition e = theta - rho (theta = u_h - u_I discrete pollution,
computed exactly; rho = u_I - u exact interpolation error, Duffy-integrated)
turns this from a table into an explanation: rho decays at exactly 2.00
per halving (the textbook interpolation rate), while theta decays at
~1.55-1.66 and dominates. So even the *smooth* problem's L2 rate on the
L-domain is capped near 5/3 (= h^{1+s}, s -> 2/3, the dual-problem
regularity limit), and the singular leg degrades one more step toward 4/3.
The quantified ladder **2.0 -> 5/3 -> 4/3** is the content of this depth
comparison; the energy column (1.0 -> 2/3) is its textbook counterpart.

## 5. Reproducibility and audit

- Every number in this report traces to `results_verified.md`, which
  carries the rerun command for each table. All scripts were rerun from a
  clean shell during the Week 9 audit (2026-10-25) and reproduced.
- The five-gate AI audit (transcription / structure / linearisation /
  independent reproduction / interrogation) for every AI-generated solver
  piece, including five caught bugs, is in `week8_ai_audit_log.md` and
  goes in the report appendix.
- Pre-computation worksheets: the ghost-closure O(h)-local/O(h^2)-global
  prediction and the L-shape rate predictions were written down before
  the runs; both predictions were confirmed.

## 6. What is still open (stated, not hidden)

- Cross-check of the reference numbers against the team's main solver --
  pending the main solver interface (see `needs_from_team.md`).
- Adaptive effectivity numbers -- pending the adaptive driver's CSV; the
  effectivity script (`week7_effectivity.py`) is ready.
- FDM vs. FEM comparison -- pending both main solvers.
- These items are listed in `results_verified.md` under "Not yet real"
  and must not be quoted as results.

---

Rong Yin, Testing & Validation seat, 2026-10-25
