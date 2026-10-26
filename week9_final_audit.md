# Week 9 -- final audit checklist (Testing & Validation)

Run this top to bottom before the common deadline (Monday Week 10, 23:59).
Every box is a claim in the report that must survive the instructor
running the code.

Status: executed 2026-10-25 from a clean shell; correction pass re-run
2026-10-27 (all FEM numbers below are from that final run). `[x]` =
verified by my scripts this week. `[team]` = depends on a teammate
deliverable, tracked in needs_from_team.md with a stated fallback; it is
*not* checked, and the report says "pending" for it.

## Direction 2 (mixed boundary data)

- [x] f = -2 exp(x+y) independently verified (week5 script, attached output)
- [x] q, r, g sampled and sign-checked (q > 0 on x=1)
- [x] four boundary rows derived AND substituted (week5_boundary_rows.py)
- [x] Dirichlet corner ownership asserted in the assembler
- [x] four-mesh tables, all four closure/BC combos, named norms,
      boundary/interior split, fitted orders stated per mesh pair
- [x] scaled residual reported on every mesh; raw residual shown once to
      justify the scaling
- [x] ||A - A^T||_F/||A||_F (3.95e-02 / 3.93e-02 at N=64) + spy plot;
      solver choice justified by nonsymmetry; definition and values
      identical to the Algorithm branch's E1 CSV
- [x] algebraic-vs-discretisation separation on a representative mesh
      (direct backward residual AND tolerance sweep)
- [x] the sentence "a small algebraic residual alone does not validate
      the boundary discretisation" appears in the report
      (report_testing_validation.md, section 2)

## Open direction 1 (L-shape + adaptive)

- [x] branch convention unit-tested incl. quadrant IV points
- [x] f = 0, zero-trace rays, nonzero outer data verified
- [x] mesh represents the L exactly (area == 3 check; no cut cells)
- [x] re-entrant corner is a Dirichlet node (zero trace on the cut edges);
      N_dof sequence identical to the Algorithm branch (161/705/2945/12033)
- [x] >= 4 uniform levels; errors vs h AND vs N_dof
- [x] energy rate ~2/3 vs h reported (measured 0.64-0.66); L2 measured
      rate 1.31-1.34 vs h = textbook 4/3 within a third of an order at
      every level (final integrator: corner-exact + Duffy n=8; the
      earlier ~1.0 and 1.25-1.32 rounds were our quadrature and free
      corner node, retracted in writing, audit log items 5, 8, 9)
- [x] cross-validation: L2 errors reproduce the Algorithm branch's E5
      table to three+ significant figures at identical N_dof
- [x] quadrature adequacy demonstrated: Duffy n=8 vs n=6 self-check at
      1e-7, corner angular quadrature's own estimate at 1e-11..1e-8,
      corner integral verified against 6x-refined brute-force Duffy at 3e-11
- [team] adaptive: marking rule stated; effectivity index table (eta/err)
      attached; improvement measured vs N_dof, uniform slope ~1/3 taken
      as the reference; mesh pictures show corner-concentrated refinement
- [team] same solver tolerance + same quadrature depth across uniform and
      adaptive sequences (matched settings) -- script ready, needs the CSV

## Depth comparisons (all four)

- [x] matched settings asserted per comparison (problem/norm/protocol/tol)
- [x] cost measure stated and measured (wall time, nnz, N_dof)
- [x] direct vs. iterative: no-crossover finding stated explicitly
- [team] FDM vs. FEM: same manufactured problem, same norms, common-node
      Linf, both N_dof stated -- gate written and self-tested, numbers
      need both main solvers
- [x] smooth vs. nonsmooth: same solver, degradation explained via
      regularity (theta/rho decomposition; measured ladder smooth
      2.0 / 1.0 vs singular 4/3 / 2/3 for L2 / energy), not waved away

## Hygiene

- [x] every report number traced to a table/script output (Ctrl-F pass:
      number -> table -> script -> command, via results_verified.md)
- [x] five-gate AI audit complete, log attached (appendix)
- [x] pre-computation worksheets filed; predictions written before runs
- [team] `python code/run_all.py` reproduces every figure on a clean env
      (team-level script owned by Algo+Viz; my package reproduces via
      the per-script commands in the README instead)
- [x] README in the zip: dependencies, how to run, what each script is
- [x] ICS written and honest (see ICS.md, signed 2026-10-27)

## If the instructor cross-questions me, I can explain:

- why the ghost closure is O(h^2) globally despite an O(h) local defect
- why scaled residuals are mandatory here (h^-2 vs h^-1 row scales)
- why CG is unjustified on paper (nonsymmetry, quantified in the Frobenius
  norm: 3.95e-02 at N=64) and what it actually did (converged, 748 vs 436
  SPD iterations, clean error)
- why the uniform-mesh L2 rate now measures 1.31-1.34 and what the two
  earlier rounds (~1.0, then 1.25-1.32) were: quadrature pollution and a
  free re-entrant corner node, both found by cross-implementation
  head-to-head and fixed, retracted in writing (audit log items 5, 8, 9)
- how the effectivity index validates (or invalidates) the adaptive loop
- where the tolerance contamination threshold is and how it was found
- why "two quadrature rules agree with each other" is not evidence that
  either is right, and why the check that finally worked was a
  differently-built implementation, not a better self-check (audit log
  items 5, 9)
