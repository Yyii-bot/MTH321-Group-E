"""
Week 6 -- four-mesh refinement tables for the mixed-boundary problem.

Owned by the Testing & Validation seat.  This file deliberately contains a
small *reference* FDM solver (mixed Dirichlet/Neumann/Robin, two boundary
closures) so that

  * the refinement tables can be produced before the main implementation
    is finished, and
  * the main solver has something independent to disagree with -- if the
    two solvers give different numbers, that IS the finding.

Solver: -Delta u = f on (0,1)^2, u* = exp(x+y), f = -2*exp(x+y).
  Dirichlet on x=0, y=0, y=1 (corners included),
  right edge x=1 in a Neumann or a Robin variant (beta = 2),
  closed either with the one-sided first-order row or the ghost-point row.

The matrix is nonsymmetric for both closures (boundary rows are scaled
differently from PDE rows, and ghost elimination doubles the inward
coupling), so we use a sparse direct solve -- CG would need justification
we cannot give.

Everything printed below is meant to be pasted into the weekly milestone
and, later, into the report's verification section.

Run:  python week6_mixed_refinement.py
"""
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


# ----------------------------------------------------------------------
# reference solver
# ----------------------------------------------------------------------
def solve_mixed(N, bc, closure, beta=2.0):
    """Assemble and solve on an (N+1)x(N+1) nodal grid.

    bc:      'neumann' or 'robin'   (right edge)
    closure: 'first' (one-sided) or 'ghost'
    Returns U (full (N+1,N+1) grid, Dirichlet nodes filled exactly),
    plus the system (A, b, dof map) for residual checks.
    """
    h = 1.0 / N
    x = np.linspace(0.0, 1.0, N + 1)
    y = np.linspace(0.0, 1.0, N + 1)
    X, Y = np.meshgrid(x, y, indexing="ij")

    # label every node; the four corners stay 'dirichlet' no matter what
    tag = np.full((N + 1, N + 1), "dirichlet", dtype=object)
    tag[1:N + 1, 1:N] = "pde"
    tag[N, 1:N] = "bc1" if closure == "first" else "bc2"

    dof = -np.ones((N + 1, N + 1), dtype=int)
    free = tag != "dirichlet"
    dof[free] = np.arange(free.sum())

    def ustar(xx, yy):
        return np.exp(xx + yy)

    rows, cols, vals = [], [], []
    rhs = np.zeros(free.sum())

    def add(i, j, coeff):
        # unknown neighbour -> matrix entry; Dirichlet node -> move its
        # column contribution to the rhs with a MINUS sign:
        #   A_ff u_f = b - A_fd g   (this sign flip is a classic bug)
        if tag[i, j] == "dirichlet":
            rhs[k] -= coeff * ustar(x[i], y[j])
        else:
            rows.append(k); cols.append(dof[i, j]); vals.append(coeff)

    for i in range(1, N + 1):
        for j in range(1, N):
            k = dof[i, j]
            if tag[i, j] == "pde":
                # five-point stencil for -Delta (positive definite form)
                add(i, j, 4.0 / h**2)
                add(i + 1, j, -1.0 / h**2)
                add(i - 1, j, -1.0 / h**2)
                add(i, j + 1, -1.0 / h**2)
                add(i, j - 1, -1.0 / h**2)
                rhs[k] += -2.0 * np.exp(x[i] + y[j])          # f = -2 exp(x+y)
            elif tag[i, j] == "bc1":
                # one-sided first-order closure: PDE row is *replaced*
                add(i, j, 1.0 / h)
                add(i - 1, j, -1.0 / h)
                if bc == "robin":
                    rows.append(k); cols.append(k); vals.append(beta)
                    rhs[k] += 3.0 * np.exp(1.0 + y[j])        # r = 3 exp(1+y)
                else:
                    rhs[k] += np.exp(1.0 + y[j])              # q = exp(1+y)
            else:  # bc2: five-point equation at the boundary node with the
                # ghost value eliminated (see week5_boundary_rows.py)
                if bc == "neumann":
                    add(i, j, 4.0 / h**2)
                    add(i - 1, j, -2.0 / h**2)                # doubled inward coupling
                    add(i, j - 1, -1.0 / h**2)
                    add(i, j + 1, -1.0 / h**2)
                    rhs[k] += -2.0 * np.exp(1.0 + y[j]) + 2.0 * np.exp(1.0 + y[j]) / h
                else:
                    add(i, j, (4.0 + 2.0 * beta * h) / h**2)  # beta lands on the diagonal
                    add(i - 1, j, -2.0 / h**2)
                    add(i, j - 1, -1.0 / h**2)
                    add(i, j + 1, -1.0 / h**2)
                    rhs[k] += -2.0 * np.exp(1.0 + y[j]) + 2.0 * 3.0 * np.exp(1.0 + y[j]) / h

    A = sp.csr_matrix((vals, (rows, cols)), shape=(free.sum(),) * 2)
    t0 = time.perf_counter()
    u = spla.spsolve(A, rhs)
    solve_time = time.perf_counter() - t0

    U = ustar(X, Y)
    U[free] = u
    return U, A, rhs, dof, tag, solve_time


# ----------------------------------------------------------------------
# error norms -- names matter, these go in the report verbatim
# ----------------------------------------------------------------------
def errors(U, N):
    h = 1.0 / N
    x = np.linspace(0.0, 1.0, N + 1)
    y = np.linspace(0.0, 1.0, N + 1)
    X, Y = np.meshgrid(x, y, indexing="ij")
    E = U - np.exp(X + Y)                       # zero on Dirichlet nodes by elimination

    l2 = np.sqrt(h**2 * np.sum(E**2))           # discrete L2, trapezoidal-equal weights
    linf = np.max(np.abs(E))                    # nodal L-infinity over all nodes
    bnd = np.max(np.abs(E[N, 1:N]))             # nodal L-infinity on the right edge only
    return l2, linf, bnd


def scaled_residual(A, b, u):
    """Row-scaled residual D^{-1}(b - A u), D_ii = max_j |A_ij|.

    Raw residuals are meaningless across meshes here: PDE rows scale like
    h^-2, first-order boundary rows like h^-1.  Always report the scaled
    version (we print both so the report can show why).
    """
    r = b - A @ u
    D = np.abs(A).max(axis=1).toarray().ravel()
    return np.max(np.abs(r)), np.max(np.abs(r) / D)


def asymmetry(A):
    """Relative Frobenius asymmetry ||A - A^T||_F / ||A||_F.

    Definition shared with the Algorithm branch (their core.py), so the
    two reports quote one number.  On the N=64 ghost-Neumann matrix this
    gives ~0.039; the inf-norm ratio we used earlier gave ~0.125 purely
    because the norm choice inflates single-row entries.
    """
    d = A - A.T
    return np.sqrt(d.power(2).sum() / A.power(2).sum())


def observed_order(errs):
    """Pairwise orders from a list of errors on h, h/2, h/4, ..."""
    return [np.log2(errs[k] / errs[k + 1]) for k in range(len(errs) - 1)]


# ----------------------------------------------------------------------
# main: the four required tables + the algebra-vs-discretisation separation
# ----------------------------------------------------------------------
if __name__ == "__main__":
    Ns = [32, 64, 128, 256]

    for bc in ("neumann", "robin"):
        for closure in ("first", "ghost"):
            print(f"\n### {bc.capitalize()}, {closure}-order closure ###")
            print("| N | h | L2 error | Linf error | right-edge Linf | scaled residual | raw residual |")
            print("|---|---|---|---|---|---|---|")
            errs = []
            for N in Ns:
                U, A, b, dof, tag, _ = solve_mixed(N, bc, closure)
                l2, linf, bnd = errors(U, N)
                raw, sc = scaled_residual(A, b, U[tag != "dirichlet"])
                errs.append((l2, linf, bnd))
                print(f"| {N} | {1/N:.5f} | {l2:.4e} | {linf:.4e} | {bnd:.4e} | {sc:.2e} | {raw:.2e} |")

            for name, idx in (("L2", 0), ("Linf", 1), ("right-edge Linf", 2)):
                orders = observed_order([e[idx] for e in errs])
                print(f"  observed {name} orders: " +
                      ", ".join(f"p={p:.2f}" for p in orders) +
                      "  (fit on " + "/".join(f"{Ns[k]}->{Ns[k+1]}" for k in range(len(Ns) - 1)) + ")")

    # --- matrix structure evidence (one representative case) ---
    print("\n### matrix structure, N=64, ghost Neumann ###")
    U, A, b, dof, tag, _ = solve_mixed(64, "neumann", "ghost")
    print(f"  nnz = {A.nnz},  ||A - A^T||_F/||A||_F = {asymmetry(A):.3e}")
    U, A, b, dof, tag, _ = solve_mixed(64, "neumann", "first")
    print(f"  same quantity, first-order closure = {asymmetry(A):.3e}")

    # --- separate algebraic from discretisation error ---
    # Direct solve: report backward residual (spsolve has no tolerance knob).
    # Then re-solve with GMRES at loose and tight tolerances and check the
    # *solution error* stops changing once the tolerance is well below the
    # discretisation error.
    print("\n### algebra vs. discretisation, N=128, ghost Neumann ###")
    N = 128
    U, A, b, dof, tag, t_direct = solve_mixed(N, "neumann", "ghost")
    u_dir = U[tag != "dirichlet"]
    free = tag != "dirichlet"
    x = np.linspace(0, 1, N + 1); y = np.linspace(0, 1, N + 1)
    X, Y = np.meshgrid(x, y, indexing="ij")
    exact_free = np.exp(X + Y)[free]
    _, sc_dir = scaled_residual(A, b, u_dir)
    print(f"  spsolve:               wall {t_direct:.3f}s  scaled residual {sc_dir:.2e}  "
          f"solution Linf error {np.max(np.abs(u_dir - exact_free)):.4e}")
    for tol in (1e-4, 1e-8, 1e-12):
        t0 = time.perf_counter()
        u_gm, info = spla.gmres(A, b, rtol=tol, maxiter=20000)
        t_gm = time.perf_counter() - t0
        _, sc_gm = scaled_residual(A, b, u_gm)
        print(f"  GMRES rtol={tol:.0e}:     wall {t_gm:.3f}s  scaled residual {sc_gm:.2e}  "
              f"info={info}  solution Linf error {np.max(np.abs(u_gm - exact_free)):.4e}  "
              f"||u_gmres - u_direct||_inf {np.max(np.abs(u_gm - u_dir)):.2e}")
    print("  -> once rtol is small enough, the error line should stop moving:")
    print("     that plateau IS the discretisation error, and the residual")
    print("     column alone would never tell you that.")
