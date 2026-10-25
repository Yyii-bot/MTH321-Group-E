# Week 6 -- milestone notes (Testing & Validation)

## What I did this week

Produced the first full refinement tables, and (true to the role) found
two real bugs in the reference solver before any number was trusted:

1. Dirichlet columns were moved to the rhs with the wrong sign in the
   assembly -- the corner-adjacent rows were off by exactly 2x the
   Dirichlet term. Caught by plugging u* into the assembled rows
   (A u* - b per row), fixed, only then were tables produced.
2. Robin first-order test row double-counted q (r already contains q + 2u*).
   Caught because the residual refused to shrink under refinement.

## Direction 2 results (u* = exp(x+y), N = 32..256)

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
  is unjustified and sparse direct is the honest choice. Spy plots in the
  report should show the doubled inward coupling of the ghost rows.
- Algebraic vs. discretisation error separated on N=128 (ghost Neumann):
  spsolve scaled residual 3.4e-15; GMRES at rtol=1e-8 still shifts the
  Linf error by ~80%, at 1e-12 agrees to 1e-9. The tolerance has to be
  ~7 orders below the discretisation error before it stops mattering.

## L-shape results (uniform refinement, n = 8..256, reference P1 FEM)

- Energy error: p = 0.64..0.66 vs h (theory 2/3), p = 0.28..0.33 vs
  N_dof (theory 1/3). Clean match.
- L2 error: p ~ 1.0-1.05 vs h, stable from n=16 to n=128. Below the optimistic
  h^{4/3}, inside the guaranteed h^{4/3-eps}. The corner region (r < 0.25)
  carries roughly a quarter of the squared L2 error -- the rest decays at
  ~h^1 as well, so the global L2 is limited everywhere, not just at the
  corner. This is a measured result to discuss, not a bug to hide.
- Quadrature: plain 3-point rules mis-measure the L2 norm badly near the
  corner (the halo r in [0.25, 0.5] still has a singular-ish integrand);
  composite refinement out to r < 0.5 (4 levels) with a 3-vs-4-level
  cross-check (~9% on L2, ~0.2% on energy) is what all reported L-shape
  numbers are based on.

> Correction (2026-10-25), recorded per the honesty rule: the L2 numbers
> above were later found to be quadrature-inflated.  The L2 signal
> (~h^4) is smaller than the O(h^3) quadrature error of any low-order
> rule applied to e^2 (grad e is O(1) pointwise, so D^3(e^2) = O(1)).
> Re-measured with the Duffy-transform tensor Gauss integrator
> (validated against scipy.integrate.dblquad; self-check 1e-5), the
> L-shape L2 errors are ~3x smaller and the order is 1.25-1.32 rising
> toward 4/3, not ~1.0.  See `results_verified.md` section 4 and
> `week8_smooth_vs_nonsmooth.py`.  Energy numbers were never affected.
- Mesh sanity: total area must equal exactly 3; any partially-cut cell is
  a hard error. The L-shape is represented exactly (cut follows grid lines).

## Plan for Week 7

Clinic readiness sheet; effectivity-index design for the adaptive study;
freeze the matched-settings protocol for the four depth comparisons.

## Blocker

None.
