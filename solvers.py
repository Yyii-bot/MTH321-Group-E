"""Sparse direct and stationary solvers with honest residual stopping."""

from __future__ import annotations

from dataclasses import dataclass, field
import time
import warnings
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .metrics import residual_metrics, row_scaling


@dataclass
class SolveResult:
    x: np.ndarray
    method: str
    converged: bool
    iterations: int
    wall_time: float
    residuals: dict
    omega: float | None = None
    final_iterate_change: float = np.nan
    history: list = field(default_factory=list)
    message: str = ""


def solve_direct(A, b, row_kind=None):
    """Solve with SciPy's sparse SuperLU driver and report backward errors."""
    A = A.tocsr()
    b = np.asarray(b, dtype=float)
    start = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("error", spla.MatrixRankWarning)
        try:
            x = spla.spsolve(A, b)
        except spla.MatrixRankWarning as exc:
            raise ValueError("The assembled matrix is singular") from exc
    elapsed = time.perf_counter() - start
    return SolveResult(
        x=np.asarray(x),
        method="direct",
        converged=True,
        iterations=0,
        wall_time=elapsed,
        residuals=residual_metrics(A, b, x, row_kind),
    )


def _scaled_relative(A, b, x, scale, rhs_scale):
    r = b - A @ x
    rho = float(np.max(np.abs(r) / scale) / rhs_scale)
    return r, rho


def _validate_red_black(A, colours):
    colours = np.asarray(colours)
    if colours.shape != (A.shape[0],) or not np.all((colours == 0) | (colours == 1)):
        raise ValueError("colours must contain one red/black label per unknown")
    for colour in (0, 1):
        ix = np.flatnonzero(colours == colour)
        block = A[ix][:, ix].tocsr()
        off = block - sp.diags(block.diagonal(), format="csr")
        if off.nnz and np.max(np.abs(off.data)) > 1e-13:
            raise ValueError("The supplied colouring has same-colour couplings")
    return colours.astype(np.int8)


def solve_jacobi(
    A,
    b,
    tol=1e-10,
    maxiter=100000,
    x0=None,
    row_kind=None,
    stop_criterion="residual",
    iterate_tol=1e-8,
    store_history=False,
):
    """Damped-undamped Jacobi (omega=1) using the common stopping interface."""
    A = A.tocsr()
    b = np.asarray(b, dtype=float)
    x = np.zeros_like(b) if x0 is None else np.array(x0, dtype=float, copy=True)
    diagonal = A.diagonal()
    if np.any(diagonal == 0.0):
        raise ValueError("Jacobi requires a nonzero diagonal")
    scale = row_scaling(A)
    rhs_scale = max(float(np.max(np.abs(b) / scale)), np.finfo(float).tiny)
    start = time.perf_counter()
    history = []
    converged = False
    final_change = np.nan
    r = b - A @ x
    iterations = 0
    for k in range(maxiter + 1):
        rho = float(np.max(np.abs(r) / scale) / rhs_scale)
        if store_history and (k % 100 == 0 or rho <= tol):
            history.append((k, rho))
        if stop_criterion == "residual" and rho <= tol:
            converged = True
            iterations = k
            break
        if k == maxiter:
            iterations = k
            break
        old = x.copy()
        x = x + r / diagonal
        r = b - A @ x
        final_change = float(np.linalg.norm(x - old, 2))
        if stop_criterion == "iterate" and final_change <= iterate_tol:
            converged = True
            iterations = k + 1
            break
    elapsed = time.perf_counter() - start
    return SolveResult(
        x=x,
        method="jacobi",
        converged=converged,
        iterations=iterations,
        wall_time=elapsed,
        residuals=residual_metrics(A, b, x, row_kind),
        final_iterate_change=final_change,
        history=history,
        message="converged" if converged else "maximum iterations reached",
    )


def solve_red_black_sor(
    A,
    b,
    colours,
    omega=1.0,
    tol=1e-10,
    maxiter=100000,
    x0=None,
    row_kind=None,
    stop_criterion="residual",
    iterate_tol=1e-8,
    store_history=False,
):
    """Red-black Gauss-Seidel/SOR; one iteration is one red+black sweep.

    The five-point grid is bipartite.  Each half-sweep is vectorised, while the
    black half-sweep uses the newly updated red values.  omega=1 is red-black
    Gauss-Seidel.  Residual stopping is the default; iterate-difference stopping
    is available only as a documented diagnostic comparison.
    """
    if not (0.0 < omega < 2.0):
        raise ValueError("SOR requires 0<omega<2")
    A = A.tocsr()
    b = np.asarray(b, dtype=float)
    colours = _validate_red_black(A, colours)
    x = np.zeros_like(b) if x0 is None else np.array(x0, dtype=float, copy=True)
    diagonal = A.diagonal()
    if np.any(diagonal == 0.0):
        raise ValueError("Gauss-Seidel/SOR requires a nonzero diagonal")
    scale = row_scaling(A)
    rhs_scale = max(float(np.max(np.abs(b) / scale)), np.finfo(float).tiny)
    red = np.flatnonzero(colours == 0)
    black = np.flatnonzero(colours == 1)
    A_red = A[red, :].tocsr()
    A_black = A[black, :].tocsr()
    start = time.perf_counter()
    history = []
    converged = False
    final_change = np.nan
    iterations = 0
    for k in range(1, maxiter + 1):
        old = x.copy()
        # D^-1(b-offdiag*x)=x+D^-1(b-Ax) on each same-colour block.
        x[red] += omega * (b[red] - A_red @ x) / diagonal[red]
        x[black] += omega * (b[black] - A_black @ x) / diagonal[black]
        final_change = float(np.linalg.norm(x - old, 2))
        r, rho = _scaled_relative(A, b, x, scale, rhs_scale)
        if store_history and (k % 100 == 0 or rho <= tol):
            history.append((k, rho))
        if stop_criterion == "residual" and rho <= tol:
            converged = True
            iterations = k
            break
        if stop_criterion == "iterate" and final_change <= iterate_tol:
            converged = True
            iterations = k
            break
        iterations = k
    elapsed = time.perf_counter() - start
    name = "gauss_seidel_rb" if omega == 1.0 else "sor_rb"
    return SolveResult(
        x=x,
        method=name,
        converged=converged,
        iterations=iterations,
        wall_time=elapsed,
        residuals=residual_metrics(A, b, x, row_kind),
        omega=float(omega),
        final_iterate_change=final_change,
        history=history,
        message="converged" if converged else "maximum iterations reached",
    )


def solve_gauss_seidel(A, b, colours, **kwargs):
    return solve_red_black_sor(A, b, colours, omega=1.0, **kwargs)


def direct_condition_estimate(A):
    """Estimate ||A||_1 ||A^-1||_1 using a reusable SuperLU factorisation."""
    A = A.tocsc()
    lu = spla.splu(A)
    inverse = spla.LinearOperator(
        A.shape,
        matvec=lu.solve,
        rmatvec=lambda vector: lu.solve(vector, trans="T"),
        dtype=float,
    )
    return float(spla.onenormest(A) * spla.onenormest(inverse))


def direct_fill_statistics(A):
    lu = spla.splu(A.tocsc())
    total = int(lu.L.nnz + lu.U.nnz)
    return {"LU_nnz": total, "fill_nnz": total - A.nnz}


def gmres_cross_check(A, b, rtol=1e-12, maxiter=300):
    """Independent row-scaled GMRES solve for direct-solver forward-error audit."""
    scale = row_scaling(A)
    As = sp.diags(1.0 / scale, format="csr") @ A.tocsr()
    bs = b / scale
    count = {"iterations": 0}

    def callback(_):
        count["iterations"] += 1

    x, info = spla.gmres(
        As,
        bs,
        rtol=rtol,
        atol=0.0,
        restart=min(100, A.shape[0]),
        maxiter=maxiter,
        callback=callback,
        callback_type="legacy",
    )
    return x, info, count["iterations"]
