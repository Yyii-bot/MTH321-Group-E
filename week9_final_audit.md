# Week 9 -- final audit checklist (Testing & Validation)

Run this top to bottom before the common deadline (Monday Week 10, 23:59).
Every box is a claim in the report that must survive the instructor
running the code.

Status: executed 2026-10-25 from a clean shell. `[x]` = verified by my
scripts this week. `[team]` = depends on a teammate deliverable, tracked
in needs_from_team.md with a stated fallback; it is *not* checked, and
the report says "pending" for it.

## Direction 2 (mixed boundary data)

- [x] f = -2 exp(x+y) independently verified (week5 script, attached output)
- [x] q, r, g sampled and sign-checked (q > 0 on x=1)
- [x] four boundary rows derived AND substituted (week5_boundary_rows.py)
- [x] Dirichlet corner ownership asserted in the assembler
- [x] four-mesh tables, all four closure/BC combos, named norms,
      boundary/interior split, fitted orders stated per mesh pair
- [x] scaled residual reported on every mesh; raw residual shown once to
      justify the scaling
- [x] ||A - A^T||/||A|| + spy plot; solver choice justified by nonsymmetry
- [x] algebraic-vs-discretisation separation on a representative mesh
      (direct backward residual AND tolerance sweep)
- [x] the sentence "a small algebraic residual alone does not validate
      the boundary discretisation" appears in the report
      (report_testing_validation.md, section 2)

## Open direction 1 (L-shape + adaptive)

- [x] branch convention unit-tested incl. quadrant IV points
- [x] f = 0, zero-trace rays, nonzero outer data verified
- [x] mesh represents the L exactly (area == 3 check; no cut cells)
- [x] >= 4 uniform levels; errors vs h AND vs N_dof
- [x] energy rate ~2/3 vs h reported (measured 0.66); L2 measured rate
      reported honestly: 1.25-1.32 vs h rising toward the h^{4/3}
      guarantee (Duffy re-measurement, 2026-10-25; the earlier ~1.0 was
      quadrature pollution and is retracted in writing, see audit log 5)
- [x] quadrature adequacy demonstrated: Duffy 3-vs-4-level check
      1e-6..1e-5, integrator independently validated against scipy dblquad
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
      regularity (theta/rho decomposition; ladder 2.0 -> 5/3 -> 4/3),
      not waved away

## Hygiene

- [x] every report number traced to a table/script output (Ctrl-F pass:
      number -> table -> script -> command, via results_verified.md)
- [x] five-gate AI audit complete, log attached (appendix)
- [x] pre-computation worksheets filed; predictions written before runs
- [team] `python code/run_all.py` reproduces every figure on a clean env
      (team-level script owned by Algo+Viz; my package reproduces via
      the per-script commands in the README instead)
- [x] README in the zip: dependencies, how to run, what each script is
- [x] ICS written and honest (see ICS.md, signed 2026-10-25)

## If the instructor cross-questions me, I can explain:

- why the ghost closure is O(h^2) globally despite an O(h) local defect
- why scaled residuals are mandatory here (h^-2 vs h^-1 row scales)
- why CG is unjustified (nonsymmetry, quantify it)
- why the uniform-mesh L2 rate measures 1.25-1.32 and not the optimistic
  4/3, what h^{4/3-eps} means, and why the first measurement said 1.0
  (quadrature pollution -- audit log item 5)
- how the effectivity index validates (or invalidates) the adaptive loop
- where the tolerance contamination threshold is and how it was found
- why "two quadrature rules agree with each other" is not evidence that
  either is right (audit log item 5)
