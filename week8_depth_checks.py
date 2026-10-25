"""
Week 8 -- matched-settings checks for the depth comparisons.

Every depth comparison in the report has to run under matched settings:
same problem, same norms, same refinement protocol, stated tolerances,
stated cost measure.  This module is the checklist made code: each
function fails loudly when two runs cannot be compared.

Also contains the tolerance-sweep driver used to separate algebraic from
discretisation error for iterative solvers (the "tolerance contamination"
hunt), and the cost measure (wall time + iteration count + nnz) recorded
for the direct-vs-iterative comparison.

Run:  python week8_depth_checks.py
"""
import time
import numpy as np
import scipy.sparse.linalg as spla

from week6_mixed_refinement import solve_mixed, errors, scaled_residual


# ----------------------------------------------------------------------
# matched-settings gate
# ----------------------------------------------------------------------
def assert_comparable(run_a, run_b):
    """Two run dicts are comparable only if these fields match."""
    keys = ("problem", "norm", "refinement", "solver_tol")
    for k in keys:
        if run_a[k] != run_b[k]:
            raise AssertionError(
                f"runs differ in '{k}': {run_a[k]!r} vs {run_b[k]!r} -- "
                f"this comparison would be invalid, fix the settings first")
    return True


# ----------------------------------------------------------------------
# tolerance sweep: where does the algebraic error stop mattering?
# ----------------------------------------------------------------------
def tolerance_sweep(N=128, bc="neumann", closure="ghost", tols=(1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12)):
    """Solve the same system at successively tighter GMRES tolerances.

    The reported error should drop with tol and then plateau; the plateau
    is the discretisation error and the tol where it starts is the
    tolerance contamination threshold for this mesh.  Stopping at a looser
    tol than that silently corrupts every error number in the report.
    """
    U, A, b, dof, tag, _ = solve_mixed(N, bc, closure)
    u_direct = U[tag != "dirichlet"]
    l2_ref, linf_ref, _ = errors(U, N)

    print(f"tolerance sweep, N={N} ({bc}, {closure}); "
          f"direct-solve reference: L2={l2_ref:.4e} Linf={linf_ref:.4e}")
    print(f"  {'rtol':>8} {'scaled res':>12} {'wall s':>8} {'L2 err':>12} {'Linf err':>12} {'contaminated?':>14}")
    for tol in tols:
        t0 = time.perf_counter()
        u, info = spla.gmres(A, b, rtol=tol, maxiter=50000)
        wall = time.perf_counter() - t0
        if info != 0:
            print(f"  {tol:8.0e} {'--':>12} {wall:8.2f}  GMRES did not converge (info={info})")
            continue
        U[tag != "dirichlet"] = u
        l2, linf, _ = errors(U, N)
        _, sc = scaled_residual(A, b, u)
        bad = "YES" if abs(linf - linf_ref) > 0.5 * linf_ref else "no"
        print(f"  {tol:8.0e} {sc:12.2e} {wall:8.2f} {l2:12.4e} {linf:12.4e} {bad:>14}")
    print("  -> anything flagged 'YES' is tolerance contamination:")
    print("     the solver error was comparable to the discretisation error.")


# ----------------------------------------------------------------------
# direct vs. iterative at matched accuracy (cost measure: wall time)
# ----------------------------------------------------------------------
def direct_vs_iterative(Ns=(64, 128, 256), target="discretisation"):
    """For each mesh, find the GMRES tolerance that matches the direct
    solve's accuracy (scaled residual ~ machine level is the strict
    match; usually 1e-10..1e-12), and record the cost of both.

    The deliverable number is the crossover: if GMRES-at-matched-accuracy
    is always slower, say so explicitly in the report -- 'no crossover in
    the tested range' is a legitimate finding, the brief asks for it.
    """
    print("| N | direct s | GMRES s (matched tol) | GMRES rtol used | iterations | speedup direct/GMRES |")
    print("|---|---|---|---|---|---|")
    for N in Ns:
        U, A, b, dof, tag, t_dir = solve_mixed(N, "neumann", "ghost")
        u_dir = U[tag != "dirichlet"]
        _, sc_dir = scaled_residual(A, b, u_dir)

        # tighten until the GMRES solution matches the direct one closely
        for tol in (1e-6, 1e-8, 1e-10, 1e-12):
            t0 = time.perf_counter()
            u, info = spla.gmres(A, b, rtol=tol, maxiter=100000)
            wall = time.perf_counter() - t0
            if info == 0 and np.max(np.abs(u - u_dir)) < 1e-8:
                break
        iters = "n/a"  # gmres in this scipy version does not return the count;
                       # wrap the solver to count if the report needs it
        print(f"| {N} | {t_dir:.3f} | {wall:.3f} | {tol:.0e} | {iters} | {wall/t_dir:.1f}x |")
    print("note: GMRES here has no preconditioner; that is a deliberate worst")
    print("case.  If the report claims solver timings, state the preconditioner")


# ----------------------------------------------------------------------
# FDM vs. FEM gate: same manufactured problem, same norms, common nodes
# ----------------------------------------------------------------------
def fdm_vs_fem_gate(fdm_run, fem_run):
    """fdm_run / fem_run: dicts with keys 'u*', 'norm', 'h', 'ndof', 'l2',
    'linf_common', 'solver_tol'.  Common-node Linf requires evaluating the
    FEM solution at FDM nodes -- the team FEM driver must provide that."""
    assert_comparable(fdm_run, fem_run)
    for k in ("l2", "linf_common", "ndof"):
        if k not in fdm_run or k not in fem_run:
            raise KeyError(f"both runs must report '{k}'")
    print("FDM vs FEM at matched settings:")
    print(f"  h = {fdm_run['h']},  N_dof FDM = {fdm_run['ndof']},  N_dof FEM = {fem_run['ndof']}")
    print(f"  L2:   FDM {fdm_run['l2']:.4e}   FEM {fem_run['l2']:.4e}")
    print(f"  Linf at common nodes: FDM {fdm_run['linf_common']:.4e}   FEM {fem_run['linf_common']:.4e}")
    ratio = fdm_run["l2"] / fem_run["l2"]
    print(f"  L2 ratio FDM/FEM = {ratio:.2f}  (state the mesh scale before attributing this to the method)")


if __name__ == "__main__":
    # demo: the matched-settings gate catching a real mismatch
    try:
        assert_comparable({"problem": "mixed2", "norm": "L2", "refinement": "h/2", "solver_tol": 1e-12},
                          {"problem": "mixed2", "norm": "Linf", "refinement": "h/2", "solver_tol": 1e-12})
    except AssertionError as e:
        print("matched-settings gate works, caught:", e)

    tolerance_sweep()

    print()
    direct_vs_iterative()
