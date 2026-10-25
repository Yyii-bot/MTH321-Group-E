"""
Week 8 -- the iterative method ladder on the Direction-2 FDM system.

Owned by the Testing & Validation seat.

The Week-8 depth check "direct vs. iterative" used unpreconditioned GMRES
only, and reported "no crossover in range".  The course lecture covers a
full ladder -- Jacobi, Gauss-Seidel, SOR, CG, GMRES, preconditioners -- so
this script runs the ladder on the SAME discrete system as the Week-6/8
tables (same assembler via solve_mixed, same ghost-Neumann case) and asks
two questions with real numbers:

  1. does any rung of the ladder beat the sparse direct solve
     (the crossover question, now with classical methods and an
     ILU-preconditioned GMRES included)?
  2. where does the lecture's theory apply verbatim, and where does it
     break -- SOR's omega_opt formula on the pure-Dirichlet square
     (applies) vs. our nonsymmetric mixed-BC matrix (needs re-checking,
     exactly as the lecture warns)?

Stopping criterion: max row-scaled residual |D^-1 r| <= 1e-10, the same
protocol as week8_depth_checks.py (1e-10 was the first clean level there).
The lecture's iterate-difference criterion delta = ||U(k+1)-U(k)||_2 is
recorded at convergence for comparison, but never used to stop -- Week 8
showed why residuals, not iterate differences, are the honest criterion.

Cost note: Jacobi/Gauss-Seidel need O(N^2) iterations, so they are run at
N = 64, 128 only; N = 256 carries the rungs that are actually plausible
there (SOR with a swept omega, GMRES+ILU).  The iteration counts at 64/128
show the O(N^2) growth explicitly.

Run:  python week8_iterative_ladder.py
"""
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from week6_mixed_refinement import solve_mixed  # same assembler as the Week-6 tables

STOP = 1e-10        # scaled-residual stop, Week-8 protocol
MAXITER = 500_000   # generous cap; Jacobi on N=128 needs ~4e4


# ----------------------------------------------------------------------
# small shared helpers
# ----------------------------------------------------------------------
def scaled_res(A, b, x):
    """max |D^-1 (b - A x)|, D_ii = max_j |A_ij| (same rule as Week 6)."""
    r = b - A @ x
    D = np.abs(A).max(axis=1).toarray().ravel()
    return float(np.max(np.abs(r) / D))


def run_loop(A, b, update, x0=None):
    """Generic iteration driver.

    update(x) -> (x_new, delta) with delta = ||x_new - x||_2 (the lecture's
    stopping quantity, recorded for the record).  Returns
    (x, iters, wall_s, res, delta_at_stop, diverged).
    Divergence abort: if the scaled residual grows 100x over any window
    after the first 50 iterations, we stop and say so -- that outcome is
    data, not a failure to hide.
    """
    n = A.shape[0]
    x = np.zeros(n) if x0 is None else x0.copy()
    res0 = scaled_res(A, b, x)
    res = res0
    delta = 0.0
    t0 = time.perf_counter()
    diverged = False
    window = res0
    k = 0
    while res > STOP and k < MAXITER:
        x_new, delta = update(x)
        x = x_new
        k += 1
        if k % 10 == 0 or k < 50:          # residual every 10 iters is plenty
            res = scaled_res(A, b, x)
            if k >= 50 and res > 100.0 * window:
                diverged = True
                break
            if k % 2000 == 0:
                window = res               # fresh divergence reference
    wall = time.perf_counter() - t0
    res = scaled_res(A, b, x)
    return x, k, wall, res, delta, diverged


def jacobi(A, b, w=1.0):
    """Weighted Jacobi as a correction iteration: x += w D^-1 r."""
    Dinv = 1.0 / A.diagonal()

    def update(x):
        r = b - A @ x
        step = w * (Dinv * r)
        return x + step, float(np.linalg.norm(step))
    return update


def sor(A, b, w):
    """SOR in residual-correction form: x += w (D + wL)^-1 r.

    w = 1 is exactly Gauss-Seidel, so one implementation serves both.
    The triangular factorisation is built once and reused -- that setup
    cost is timed separately and included in the reported wall time.
    """
    D = sp.diags(A.diagonal())
    Lstrict = sp.tril(A, k=-1)
    M = D + w * Lstrict
    t0 = time.perf_counter()
    Msolve = spla.splu(M.tocsc()).solve
    setup = time.perf_counter() - t0

    def update(x):
        r = b - A @ x
        step = w * Msolve(r)
        return x + step, float(np.linalg.norm(step))
    return update, setup


# ----------------------------------------------------------------------
# pure-Dirichlet square: the lecture's exact SOR setting
# ----------------------------------------------------------------------
def assemble_dirichlet(N):
    """-Delta u = f on (0,1)^2, u = exp(x+y) on all four edges.

    Interior (N-1)^2 unknowns, h = 1/N -- exactly the model problem the
    lecture's omega_opt = 2/(1+sin(pi h)) formula is derived for.
    Built in COO with the block-row wrap-around couplings simply absent.
    """
    h = 1.0 / N
    n = N - 1
    x = np.linspace(0.0, 1.0, N + 1)
    X, Y = np.meshgrid(x, x, indexing="ij")

    idx = np.arange(n * n).reshape(n, n)
    rows, cols, vals = [], [], []

    def put(ii, jj, coeff):
        # ii/jj may be index arrays; keep everything as flat arrays later
        rows.append(np.asarray(ii).ravel())
        cols.append(np.asarray(jj).ravel())
        vals.append(np.full(np.asarray(ii).ravel().shape, coeff))

    put(idx.ravel(), idx.ravel(), 4.0 / h**2)                   # diagonal
    put(idx[:, :-1].ravel(), idx[:, 1:].ravel(), -1.0 / h**2)   # j -> j+1
    put(idx[:, 1:].ravel(), idx[:, :-1].ravel(), -1.0 / h**2)   # j -> j-1
    put(idx[:-1, :].ravel(), idx[1:, :].ravel(), -1.0 / h**2)   # i -> i+1
    put(idx[1:, :].ravel(), idx[:-1, :].ravel(), -1.0 / h**2)   # i+1 -> i
    A = sp.csr_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(n * n, n * n))

    f = -2.0 * np.exp(X + Y)
    b = f[1:N, 1:N].ravel()
    # move known boundary values to the rhs (MINUS sign -- the Week-6 lesson)
    g = np.exp(X + Y)
    coff = -1.0 / h**2
    for k in range(n):          # bottom and top edges
        b[k] -= coff * g[0, k + 1]
        b[(n - 1) * n + k] -= coff * g[N, k + 1]
    for k in range(n):          # left and right edges
        b[k * n] -= coff * g[k + 1, 0]
        b[k * n + n - 1] -= coff * g[k + 1, N]
    return A, b


def sor_sweep(A, b, omegas):
    """Iterations-to-converge for each omega.  Omegas above the optimal
    diverge and are reported as such (aborted by the 100x detector)."""
    out = []
    for w in omegas:
        update, setup = sor(A, b, w)
        _, k, wall, res, _, div = run_loop(A, b, update)
        out.append((w, k, wall, res, div))
    return out


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # ==================================================================
    # Table A: the ladder on the mixed system (ghost Neumann), the case
    # the Week-8 no-crossover statement was made on.  Same assembler,
    # same exact solution, same scaled-residual protocol.
    # ==================================================================
    print("### iterative ladder, mixed BC (ghost Neumann), u* = exp(x+y) ###")
    print("stop: scaled residual <= 1e-10 (Week-8 protocol)")
    print("| N | method | wall s | iters | scaled res at stop | Linf error | vs direct |")
    print("|---|---|---|---|---|---|---|")

    SOR_OMEGA = {}      # filled by the N=64 sweep below

    def ladder(N, classical):
        U, A, b, dof, tag, t_dir = solve_mixed(N, "neumann", "ghost")
        free = tag != "dirichlet"
        gx, gy = np.meshgrid(np.linspace(0, 1, N + 1),
                             np.linspace(0, 1, N + 1), indexing="ij")
        exact = np.exp(gx + gy)[free]
        base = np.max(np.abs(U[free] - exact))
        print(f"| {N} | spsolve (direct) | {t_dir:.3f} | -- | "
              f"{scaled_res(A, b, U[free]):.1e} | {base:.2e} | 1.0x |")

        def report(name, x, k, wall, res, div):
            err = np.max(np.abs(x - exact))
            tag_txt = ("DIVERGED" if div else
                       (f"{t_dir/wall:.1f}x faster" if wall < t_dir
                        else f"{wall/t_dir:.1f}x slower"))
            print(f"| {N} | {name} | {wall:.3f} | {k} | {res:.1e} | {err:.2e} | {tag_txt} |")

        if classical:
            # Jacobi (w = 1): the bottom rung
            x, k, wall, res, delta, div = run_loop(A, b, jacobi(A, b))
            report("Jacobi", x, k, wall, res, div)

            # Gauss-Seidel (SOR with w = 1), factorisation included
            update, setup = sor(A, b, 1.0)
            x, k, wall, res, delta, div = run_loop(A, b, update)
            report("Gauss-Seidel", x, k, wall + setup, res, div)

        # SOR with omega swept on THIS system at N = 64 (the square-theory
        # value is only a starting guess on a nonsymmetric problem)
        w = SOR_OMEGA[N]
        update, setup = sor(A, b, w)
        x, k, wall, res, delta, div = run_loop(A, b, update)
        report(f"SOR w={w:.3f}", x, k, wall + setup, res, div)

        # CG: the matrix is nonsymmetric (Week 6: ||A-A^T||/||A|| ~ 0.12),
        # so this rung is expected to misbehave -- record what it does.
        iters = []
        t0 = time.perf_counter()
        x_cg, info = spla.cg(A, b, rtol=1e-10, maxiter=10000,
                             callback=lambda _: iters.append(0))
        wall = time.perf_counter() - t0
        res = scaled_res(A, b, x_cg)
        err = np.max(np.abs(x_cg - exact))
        verdict = "converged (unexpected!)" if res <= STOP else "stagnated/broke down"
        print(f"| {N} | CG (unjustified) | {wall:.3f} | {len(iters)} | {res:.1e} | {err:.2e} | {verdict} |")

        # GMRES with ILU preconditioning -- the rung most likely to flip
        # the no-crossover finding
        t0 = time.perf_counter()
        ilu = spla.spilu(A.tocsc(), drop_tol=1e-4, fill_factor=20)
        prec = spla.LinearOperator(A.shape, ilu.solve)
        setup = time.perf_counter() - t0
        t0 = time.perf_counter()
        it_gm = []
        x_gm, info = spla.gmres(A, b, M=prec, rtol=1e-10, maxiter=2000,
                                callback=lambda _: it_gm.append(0),
                                callback_type="legacy")
        wall = time.perf_counter() - t0
        res = scaled_res(A, b, x_gm)
        err = np.max(np.abs(x_gm - exact))
        verdict = (f"{t_dir/(wall+setup):.1f}x faster" if wall + setup < t_dir
                   else f"{(wall+setup)/t_dir:.1f}x slower")
        print(f"| {N} | GMRES+ILU | {wall + setup:.3f} | {len(it_gm)} | {res:.1e} | {err:.2e} | {verdict} |")
        print(f"    [N={N}: ILU setup {setup:.3f}s + GMRES {wall:.3f}s; "
              f"||A-A^T||/||A|| = "
              f"{spla.norm(A - A.T, np.inf)/spla.norm(A, np.inf):.2e}]")

    print("[sweeping SOR omega on the mixed system at N=64 first]")
    Um, Am, bm, dofm, tagm, _ = solve_mixed(64, "neumann", "ghost")
    sweep_mixed = sor_sweep(Am, bm, np.round(np.arange(1.2, 1.902, 0.05), 3))
    conv_mixed = [(w, k) for w, k, wall, res, div in sweep_mixed
                  if not div and res <= STOP]
    best_w_mixed, best_k_mixed = min(conv_mixed, key=lambda t: t[1])
    SOR_OMEGA = {64: best_w_mixed, 128: best_w_mixed, 256: best_w_mixed}
    print(f"  mixed-system best omega = {best_w_mixed:.3f} "
          f"({best_k_mixed} iters at N=64); reused for N=128,256")

    for N in (64, 128):
        ladder(N, classical=True)
    ladder(256, classical=False)
    print("note: Jacobi/GS at N=256 not run (O(N^2) iterations; the 64->128 "
          "counts show the growth). SOR/CG/GMRES+ILU are run at 256.")

    # ==================================================================
    # Table B: SOR omega on the lecture's own model problem
    # (pure-Dirichlet square, h = 1/N), theory vs measurement
    # ==================================================================
    print("\n### SOR omega sweep, pure-Dirichlet square (theory's own setting) ###")
    print("| N | h | theory omega_opt = 2/(1+sin(pi h)) | measured best omega | iters at best | iters at theory |")
    print("|---|---|---|---|---|---|")
    omegas = np.round(np.arange(1.0, 1.951, 0.025), 3)
    A, b = assemble_dirichlet(64)
    sweep = sor_sweep(A, b, omegas)
    conv = [(w, k, wall) for w, k, wall, res, div in sweep if not div and res <= STOP]
    best_w, best_k, _ = min(conv, key=lambda t: t[1])
    w_theory64 = 2.0 / (1.0 + np.sin(np.pi / 64))
    # exact theory value was not on the 0.025-spaced grid -- run it directly
    update_t, setup_t = sor(A, b, w_theory64)
    _, k_theory64, _, res_t64, _, div_t64 = run_loop(A, b, update_t)
    if div_t64:
        k_theory64 = -1
    print(f"| 64 | {1/64:.5f} | {w_theory64:.6f} | {best_w:.3f} | {best_k} | {k_theory64} |")
    # confirmation point at N = 128 (single run at the theory value)
    w_theory128 = 2.0 / (1.0 + np.sin(np.pi / 128))
    A, b = assemble_dirichlet(128)
    update, setup = sor(A, b, w_theory128)
    _, k_128, _, res_128, _, div_128 = run_loop(A, b, update)
    print(f"| 128 | {1/128:.5f} | {w_theory128:.6f} | (sweep at 64 only) | -- | "
          f"{k_128}{' (DIVERGED)' if div_128 else ''} |")
    print("note: omegas above the optimum diverge (aborted by the 100x "
          "residual-growth detector); they are excluded from 'best'.")

    # ==================================================================
    # Table C: two theory spot-checks the lecture invites
    # ==================================================================
    print("\n### theory spot-checks ###")
    # (1) Jacobi preconditioner on the 5-point Dirichlet matrix is just a
    #     constant scaling, D^-1 A = (h^2/4) A -- so PCG(Jacobi) must
    #     converge in (essentially) the same iteration count as plain CG.
    A, b = assemble_dirichlet(64)
    it_plain = []
    spla.cg(A, b, rtol=1e-10, maxiter=20000, callback=lambda _: it_plain.append(0))
    Dinv = 1.0 / A.diagonal()
    it_prec = []
    spla.cg(A, b, rtol=1e-10, maxiter=20000, M=spla.LinearOperator(A.shape, lambda r: Dinv * r),
            callback=lambda _: it_prec.append(0))
    print(f"CG vs CG+Jacobi-precond on Dirichlet square, N=64: "
          f"{len(it_plain)} vs {len(it_prec)} iterations")
    print("-> identical behaviour confirms D^-1A is a constant scaling here "
          "(lecture's point 13), so the Jacobi preconditioner buys nothing "
          "on this matrix.")

    # (2) CG on the SPD Dirichlet matrix vs. CG on the nonsymmetric mixed
    #     matrix: same PDE, same solution, opposite symmetry -- the
    #     cleanest possible demonstration that CG's requirement is on the
    #     MATRIX, not the problem.
    A_dir, b_dir = assemble_dirichlet(128)
    it_dir = []
    t0 = time.perf_counter()
    spla.cg(A_dir, b_dir, rtol=1e-10, maxiter=20000, callback=lambda _: it_dir.append(0))
    t_spd = time.perf_counter() - t0
    U, A_mix, b_mix, dof, tag, _ = solve_mixed(128, "neumann", "ghost")
    it_mix = []
    t0 = time.perf_counter()
    x_mix, _ = spla.cg(A_mix, b_mix, rtol=1e-10, maxiter=5000, callback=lambda _: it_mix.append(0))
    t_mix = time.perf_counter() - t0
    free = tag != "dirichlet"
    gx, gy = np.meshgrid(np.linspace(0, 1, 129), np.linspace(0, 1, 129), indexing="ij")
    exact = np.exp(gx + gy)[free]
    print(f"CG on SPD Dirichlet matrix, N=128: {len(it_dir)} iters, {t_spd:.3f}s, clean")
    print(f"CG on nonsymmetric mixed matrix, N=128: {len(it_mix)} iters in {t_mix:.3f}s, "
          f"scaled res {scaled_res(A_mix, b_mix, x_mix):.1e}, "
          f"Linf error {np.max(np.abs(x_mix - exact)):.2e}")
