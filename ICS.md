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
   (`week6_lshape_fem.py`): energy order 0.64-0.66 vs h (theory 2/3),
   L2 order 1.31-1.34 vs h (theory 4/3), mesh area == 3 hard check, and
   the final corner-exact + Duffy-n8 error integrator whose N_dof
   sequence (161/705/2945/12033) and absolute L2 values match the
   Algorithm branch's E5 table to three+ significant figures on every
   shared level.

4. Found and retracted my own wrong numbers, twice, in writing. Round
   one: the Week-6 L2 column was quadrature-inflated (two low-order
   rules agreed with each other to <1e-3 but both drifted ~50% across
   refinement levels); re-measured with a Duffy tensor-Gauss integrator.
   Round two (2026-10-27): the re-entrant corner had been left as a free
   unknown, violating the zero trace of u* on the cut edges, and the
   off-corner quadrature was still low-order — fixing both moved the
   singular leg to textbook 4/3 and the smooth leg to a flat 2.00,
   retracting the published "~5/3 pollution cap" as an artifact. Final
   table reproduces the main implementation's E5 numbers. The
   retractions are in results_verified.md section 4 and audit-log
   items 5, 8, 9, not silently replaced.

5. Ran the smooth-vs-nonsmooth depth comparison on the same solver,
   meshes and norms: smooth leg L2 2.00 flat / energy 1.00, singular leg
   L2 1.31-1.34 / energy 0.64-0.66, with the theta/rho decomposition
   showing both interpolation and pollution decay at ~2.0 for smooth
   data — the measured contrast is the textbook regularity ladder
   2.0 / 1.0 (smooth) vs 4/3 / 2/3 (singular) for L2 / energy
   (`week8_smooth_vs_nonsmooth.py`, results_verified.md section 6).

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
   matrix (748 vs 436 SPD iterations; ||A-A^T||_F/||A||_F = 3.95e-02)
   instead of the predicted breakdown.

8. Wrote and maintained the evidence ledger (`results_verified.md`, the
   only file the report quotes numbers from), the five-gate AI audit log
   (`week8_ai_audit_log.md`, including the caught bugs, two recorded
   prediction misses, and the two-round L2 correction), the clinic
   readiness sheet, the effectivity-index script awaiting the adaptive
   CSV, and the final pre-submission audit (`week9_final_audit.md`).
   Every claim above traces to a commit and a rerun command.

## What I would do differently

Budget the quadrature and boundary-node audits before the first L2
table, not after: the corner quadrature failure cost a written
retraction and two re-runs of the L-shape study, and the free corner
node hid inside a plausible-looking rate for a week. Two cheap habits
would have caught both on day one: a manufactured-integral test of the
integrator, and a head-to-head against any independently-built solver
before publishing absolute values.

## Anything I should have asked for help with earlier

The main solver interface for the cross-check. I scoped the diff script
in Week 7 but should have pushed harder for a callable interface by Week
8 — the head-to-head against the Algorithm branch's E5 table on
2026-10-27 is what finally settled my L2 numbers, and it could have
happened a week earlier.

---

Name: Rong Yin   Date: 2026-10-27
