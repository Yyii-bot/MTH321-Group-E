# What I need from the team (Testing & Validation seat)

Status of every external input my evidence chain is missing. "Blocks" =
which report claim stays unbacked without it. Fallback = what I do if it
never arrives.

| # | item | from | format / spec | needed by | blocks | fallback if missing |
|---|---|---|---|---|---|---|
| 1 | Sign convention agreement: -Delta u = f, f = -2 exp(x+y) | Algo (Maths confirm) | one-line reply + code uses it | Week 5 (now) | matched-settings between solvers; FDM vs FEM comparison | I keep the bank convention and document the discrepancy; comparisons stay internal to my reference code |
| 2 | Main FDM solver interface: callable as solve(N, bc, closure) -> nodal grid U | Algo | function + 2-min usage note | Week 6-7 | cross-check of every Direction-2 table number against an independent implementation | report numbers as "reference solver, cross-check pending" -- weaker but honest |
| 3 | Cross-check run of (2) vs my reference solver | me, once (2) lands | diff report: ||U_main - U_ref||_inf per (N, bc, closure), expect < 1e-8 | Week 7 (clinic) | the word "independent" in the verification chain | none -- do it the day (2) lands, ~10 min |
| 4 | Adaptive driver CSV: one row per level `level, n_dof, energy_error, eta` | Algo + FEM (Maths picks indicator) | CSV; energy_error MUST use corner-aware composite quadrature (same rule as week6_lshape_fem.py) | Week 8 mid | the entire adaptive-validation claim (effectivity index, uniform-vs-adaptive figure, "adaptive wins" statement) | drop the adaptive claim to a proposal; report uniform study only -- loses the innovation point |
| 5 | Estimator + marking decision: indicator type (gradient-jump?), Dörfler theta | Maths proposes, team agrees | one line in the Week-7 milestone | Week 7 clinic | (4) cannot start | I proposed theta = 0.5 in week7_clinic_readiness.md; team can override |
| 6 | FEM solver on the L-shape (scikit-fem or own) callable for uniform + adaptive meshes | FEM seat | function + mesh/uh output format | Week 8 | FDM-vs-FEM depth item; independent reproduction of my L-shape rates | my mini-FEM is labelled "reference"; FDM-vs-FEM item reported as "not completed in range" |
| 7 | run_all.py regenerating every report figure | Algo + Viz | runs end-to-end from repo root | Week 9 | Dimension C (reproducibility); my Week-9 audit item "clean env rerun" | zip my scripts with README as a partial reproduction package |
| 8 | Report drafts (verification section source) | Viz | shared file | Week 8-9 | number-tracing pass (every figure traces to a table) | I trace against results_verified.md only |
| 9 | L2-norm head-to-head: their corner-exact quadrature vs my Duffy values on the singular leg | me + Algo | same mesh, same solution, both error functionals printed side by side | before report freeze | mixing absolute L2 values from the two packages in one table (they differ ~3x at identical dof; fitted orders agree) | report the two norm definitions separately and say why; never average them |

Status updates 2026-10-26: (1) confirmed — their data.py matches the bank
convention. (2)+(3) closed — the main solver's ghost-Neumann errors agree
with my reference solver to six significant figures at N = 32/64/128 (L2
and Linf), so the Direction-2 cross-check is done. (9) is new from that
review: singular-leg L2 absolute values differ ~3.2x at identical dof
(theirs corner-exact adaptive angular quadrature ~1e-14, mine Duffy +
tensor Gauss n=3), orders agree (both rising toward 4/3). Open.

## What I do NOT need

- Their data to re-derive my Week 5-6 numbers: those stand on my scripts.
- Benchmark numbers from the instructor: the course deliberately hands
  none out; my manufactured checks replace them.
- A preconditioner study: the direct-vs-iterative finding is reported as
  "no crossover with unpreconditioned GMRES", which is complete as stated.

## Offers (things I can produce for them)

| offer | for whom | effort |
|---|---|---|
| Auto cross-check script: main vs reference solver diff report | Algo | 30 min once interface (2) exists |
| Corner-aware error-norm function they can import for the FEM driver | FEM seat | already exists: integrate_errors in week6_lshape_fem.py |
| clinic readiness Q&A rehearsal | whole team | 30 min before Seminar 3 |
