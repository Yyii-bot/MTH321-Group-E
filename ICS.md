# Individual Contribution Statement — Project 2 (Elliptic BVPs)

**Name:** Rong Yin
**Role this project:** Testing & Validation

## The specific things I did this project

1. Verified all manufactured data before any solver existed
   (`week5_verify_analytic.py`, `week5_boundary_rows.py`): re-derived
   f = -2 exp(x+y) with a five-point stencil, sign-checked q and r on
   all edges, unit-tested the L-shape branch convention including
   quadrant-IV points, and verified all four right-edge discrete rows by
   substituting sampled u* values into the assembled rows.

2. Built the reference FDM solver (`week6_mixed_refinement.py`),
   deliberately independent of the main implementation, and caught the
   Dirichlet-to-RHS sign error by plugging u* into A u* - b per row
   while every scaled residual sat at 1e-15 — the reason the report
   says "a small algebraic residual alone does not validate the boundary
   discretisation". Produced the four-mesh convergence tables for all
   four closure/BC combinations (first-order ~1 everywhere, ghost ~2
   everywhere).

3. Built the reference P1 FEM on exactly-represented L-shape meshes
   (`week6_lshape_fem.py`): energy order 0.66 vs h (theory 2/3), mesh
   area == 3 hard check, corner-aware quadrature with a 3-vs-4-level
   cross-check column.

4. Found and retracted my own wrong numbers: the Week-6 L2 column was
   quadrature-inflated (two low-order rules agreed with each other to
   <1e-3 but both drifted ~50% across refinement levels). Re-measured
   with a Duffy-transform tensor-Gauss integrator
   (`week8_smooth_vs_nonsmooth.py`), validated against scipy dblquad,
   giving the corrected L2 order 1.25-1.32 rising toward 4/3. The
   retraction is written down in results_verified.md section 4 and
   audit-log item 5, not silently replaced.

5. Ran the smooth-vs-nonsmooth depth comparison on the same solver,
   meshes and norms: smooth leg L2 1.65-1.70 / energy 1.00, singular leg
   L2 1.25-1.32 / energy 0.66, with the theta/rho decomposition showing
   the smooth leg is itself pollution-capped near 5/3 — the quantified
   ladder 2.0 -> 5/3 -> 4/3 (`week8_smooth_vs_nonsmooth.py`,
   results_verified.md section 6).

6. Quantified solver-tolerance contamination (`week8_depth_checks.py`):
   GMRES at rtol=1e-8 still corrupts the solution error by ~80% with an
   apparently excellent residual; 1e-10 is the first clean level. Froze
   the scaled-residual stopping criterion the whole package now uses.

7. Ran the iterative method ladder on the Direction-2 matrix
   (`week8_iterative_ladder.py`): Jacobi/GS/SOR/CG/GMRES+ILU at the same
   stop criterion, N = 64..256. No rung beats the sparse direct solve —
   even GMRES+ILU is 1.5-1.9x slower — so the no-crossover finding
   survives preconditioning. Verified the lecture's SOR omega_opt
   formula on its own model problem (270 iterations at the theory value,
   N=64) and recorded the surprise that CG converged on the nonsymmetric
   matrix (748 vs 436 SPD iterations) instead of the predicted breakdown.

8. Wrote and maintained the evidence ledger (`results_verified.md`, the
   only file the report quotes numbers from), the five-gate AI audit log
   (`week8_ai_audit_log.md`, including five caught bugs and two recorded
   prediction misses), the clinic readiness sheet, the effectivity-index
   script awaiting the adaptive CSV, and the final pre-submission audit
   (`week9_final_audit.md`). Every claim above traces to a commit and a
   rerun command.

## What I would do differently

Budget the quadrature validation before the first L2 table, not after:
the corner quadrature failure cost a written retraction and a re-run of
the L-shape study. A cheap manufactured-integral test of the integrator
on day one would have caught it in minutes.

## Anything I should have asked for help with earlier

The main solver interface for the cross-check. I scoped the diff script
in Week 7 but should have pushed harder for a callable interface by Week
8 — the "reference solver, cross-check pending" label on my tables could
have been removed before the deadline instead of after it.

---

Name: Rong Yin   Date: 2026-10-26
