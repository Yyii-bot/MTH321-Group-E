"""
Week 5 -- check the right-edge boundary rows by hand (then by script).

The problem pack asks for the first-order and ghost-point rows to be
verified by substituting sampled values of u*.  Doing it by hand once is
worth it; doing it for a grid of h values is what this script is for.

What we expect (and what the clinic will ask about):
  * first-order row:  (u*(N,j) - u*(N-1,j))/h - q_j  is O(h), and this
    defect DOES limit the global error to first order.
  * ghost row: the eliminated boundary equation has an O(h) *local*
    residual -- that is normal (LeVeque, FDML, elliptic chapter).  The
    defect stays confined to the boundary layer and the *global* solution
    error is still O(h^2).  So don't read the second-order claim off this
    table; it has to show up in the Week 6 solution-error refinement table.

Run:  python week5_boundary_rows.py
"""
import numpy as np

U = lambda x, y: np.exp(x + y)          # manufactured solution, direction 2
Q = lambda y: np.exp(1.0 + y)           # Neumann data on x=1 (outward normal (1,0))
R = lambda y: 3.0 * np.exp(1.0 + y)     # Robin data, beta = 2: r = q + 2 u*
F = lambda x, y: -2.0 * np.exp(x + y)   # forcing, -Delta u = f convention
BETA = 2.0


def first_order_residual(N, yj, kind):
    """Substitute exact nodal values into the one-sided boundary row.

    kind = 'neumann': (U_N - U_{N-1})/h - q_j
    kind = 'robin':   (U_N - U_{N-1})/h + beta U_N - r_j
    Both should shrink like O(h) as the grid refines.
    """
    h = 1.0 / N
    uN, uW = U(1.0, yj), U(1.0 - h, yj)
    if kind == "neumann":
        res = (uN - uW) / h - Q(yj)      # (U_N - U_{N-1})/h = q_j
    else:
        # Robin row is (U_N - U_{N-1})/h + beta U_N = r_j.
        # Note: NO separate q term -- r already contains q + beta u*.
        res = (uN - uW) / h + BETA * uN - R(yj)
    return res


def ghost_row_residual(N, yj, kind):
    """Substitute exact values into the *ghost-eliminated* boundary equation.

    Start from the five-point equation at (N,j) with the ghost value
        U_{N+1} = U_{N-1} + 2 h q_j              (Neumann)
        U_{N+1} = U_{N-1} + 2 h (r_j - beta U_N) (Robin)
    substituted in.  The remaining equation is what actually lands in the
    matrix -- this residual should shrink like O(h^2) (scaled by f).

    Row (for neumann):  (4U_N - 2U_{N-1} - U_{N,j-1} - U_{N,j+1})/h^2 = f + 2q/h
    """
    h = 1.0 / N
    u = U(1.0, yj)
    uW = U(1.0 - h, yj)
    uS, uN = U(1.0, yj - h), U(1.0, yj + h)

    if kind == "neumann":
        lhs = (4 * u - 2 * uW - uS - uN) / h**2
        rhs = F(1.0, yj) + 2 * Q(yj) / h
    else:  # robin: ghost carries (r_j - beta*U_N), so beta moves onto the diagonal
        lhs = ((4 + 2 * BETA * h) * u - 2 * uW - uS - uN) / h**2
        rhs = F(1.0, yj) + 2 * R(yj) / h
    return lhs - rhs


def show(label, residual_fn, expect_power):
    print(f"-- {label} --")
    print(f"   {'N':>5} {'h':>8} {'residual':>14} {'res/h^p':>12}")
    prev = None
    for N in (16, 32, 64, 128):
        h = 1.0 / N
        r = abs(residual_fn(N, 0.37))          # fixed sample point y = 0.37
        ratio = "" if prev is None else f"{prev / r:6.2f}x"
        print(f"   {N:>5} {h:8.4f} {r:14.6e} {ratio:>12}")
        prev = r
    print(f"   (expect the residual to shrink like h^{expect_power}, i.e. ratio -> {2**expect_power})\n")


def corner_ownership(tag, N):
    """The four corners belong to the Dirichlet edges and must never be
    overwritten by a Neumann or Robin row.  Call this on whatever array the
    assembler uses to label nodes (values like 'dirichlet', 'pde', 'bc')."""
    corners = {"(0,0)": tag[0, 0], "(0,N)": tag[0, N], "(N,0)": tag[N, 0], "(N,N)": tag[N, N]}
    bad = {k: v for k, v in corners.items() if v != "dirichlet"}
    if bad:
        raise AssertionError(f"corners overwritten by boundary rows: {bad}")
    return True


if __name__ == "__main__":
    show("first-order Neumann row  (expect O(h), ratio -> 2)",
         lambda N, y: first_order_residual(N, y, "neumann"), 1)
    show("first-order Robin row    (expect O(h), ratio -> 2)",
         lambda N, y: first_order_residual(N, y, "robin"), 1)
    # ghost rows: O(h) *local* residual is the expected behaviour; the
    # second-order payoff shows up in the solution error, not here.
    show("ghost Neumann row        (expect O(h) locally, ratio -> 2; global error is still O(h^2))",
         lambda N, y: ghost_row_residual(N, y, "neumann"), 1)
    show("ghost Robin row          (expect O(h) locally, ratio -> 2; global error is still O(h^2))",
         lambda N, y: ghost_row_residual(N, y, "robin"), 1)

    # corner ownership demo on a fake tag array (replace with the real one)
    N = 8
    tag = np.full((N + 1, N + 1), "dirichlet", dtype=object)
    tag[1:N + 1, 1:N] = "pde"
    tag[N, 1:N] = "bc"
    print("corner ownership on a correct tag grid:", corner_ownership(tag, N))
