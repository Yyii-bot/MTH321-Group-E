"""Five-point Dirichlet solvers on square and aligned L-shaped grids."""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import scipy.sparse as sp


@dataclass(frozen=True)
class Grid2D:
    N: int
    h: float
    x: np.ndarray
    y: np.ndarray
    k_of: np.ndarray
    active: np.ndarray
    node_type: np.ndarray  # -1 inactive, 0 Dirichlet, 1 unknown
    name: str

    @property
    def n_unknowns(self):
        return int(np.max(self.k_of)) + 1

    def unknown_ij(self):
        """Return coordinates in ascending matrix-unknown order, not array order."""
        i, j = np.nonzero(self.k_of >= 0)
        order = np.argsort(self.k_of[i, j])
        return i[order], j[order]

    def unknown_coordinates(self):
        i, j = self.unknown_ij()
        return self.x[i], self.y[j]

    def coordinates(self):
        return np.meshgrid(self.x, self.y, indexing="ij")

    def full_solution(self, values, g):
        full = np.full(self.k_of.shape, np.nan, dtype=float)
        x, y = self.coordinates()
        dirichlet = self.node_type == 0
        full[dirichlet] = np.asarray(g(x, y))[dirichlet]
        i, j = self.unknown_ij()
        full[i, j] = values
        return full


@dataclass
class FDMSystem:
    A: sp.csr_matrix
    b: np.ndarray
    grid: Grid2D
    forcing: object
    g: object


def _make_grid(N, x, y, active, node_type, name):
    k_of = -np.ones((N + 1, N + 1), dtype=np.int32)
    unknown = np.argwhere(node_type == 1)
    # Row-major in j, with i fastest, matching the core index convention.
    unknown = sorted(unknown, key=lambda ij: (ij[1], ij[0]))
    for k, (i, j) in enumerate(unknown):
        k_of[i, j] = k
    return Grid2D(N, 1.0 / N if name == "square" else 2.0 / N,
                  x, y, k_of, active, node_type, name)


def square_grid(N):
    if N < 2:
        raise ValueError("N must be >=2")
    x = np.linspace(0.0, 1.0, N + 1)
    y = x.copy()
    active = np.ones((N + 1, N + 1), dtype=bool)
    node_type = np.ones((N + 1, N + 1), dtype=np.int8)
    node_type[0, :] = 0
    node_type[N, :] = 0
    node_type[:, 0] = 0
    node_type[:, N] = 0
    return _make_grid(N, x, y, active, node_type, "square")


def lshape_grid(N):
    """Aligned grid for (-1,1)^2 minus the upper-right open quadrant.

    N must be even so x=0 and y=0 are grid lines.  Re-entrant rays and the
    corner are Dirichlet nodes; nodes strictly inside the removed quadrant are
    inactive and never receive unknowns.
    """
    if N < 4 or N % 2:
        raise ValueError("The L-shaped grid requires an even N>=4")
    x = np.linspace(-1.0, 1.0, N + 1)
    y = x.copy()
    mid = N // 2
    active = np.ones((N + 1, N + 1), dtype=bool)
    node_type = np.ones((N + 1, N + 1), dtype=np.int8)
    for i in range(N + 1):
        for j in range(N + 1):
            removed = i > mid and i < N and j > mid and j < N
            outer = i in (0, N) or j in (0, N)
            reentrant = (i == mid and j >= mid) or (j == mid and i >= mid)
            if removed:
                active[i, j] = False
                node_type[i, j] = -1
            elif outer or reentrant:
                node_type[i, j] = 0
    return _make_grid(N, x, y, active, node_type, "L-shape")


def assemble_five_point_dirichlet(grid, forcing, g):
    """COO assembly with Dirichlet elimination and positive -Delta rows."""
    h = grid.h
    inv_h2 = 1.0 / h**2
    rows, cols, values = [], [], []
    b = np.zeros(grid.n_unknowns, dtype=float)
    for j in range(grid.N + 1):
        for i in range(grid.N + 1):
            k = grid.k_of[i, j]
            if k < 0:
                continue
            rows.append(k)
            cols.append(k)
            values.append(4.0 * inv_h2)
            b[k] = forcing(grid.x[i], grid.y[j])
            for ii, jj in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)):
                kk = grid.k_of[ii, jj]
                if kk >= 0:
                    rows.append(k)
                    cols.append(kk)
                    values.append(-inv_h2)
                elif grid.active[ii, jj]:
                    b[k] += inv_h2 * g(grid.x[ii], grid.y[jj])
                else:
                    raise RuntimeError("An unknown has an inactive immediate neighbour")
    A = sp.coo_matrix(
        (values, (rows, cols)), shape=(grid.n_unknowns, grid.n_unknowns)
    ).tocsr()
    A.sum_duplicates()
    return FDMSystem(A, b, grid, forcing, g)


def grid_error_metrics(system, values, exact):
    grid = system.grid
    i, j = grid.unknown_ij()
    u_exact = np.asarray(exact(grid.x[i], grid.y[j]))
    e = values - u_exact
    full_u = grid.full_solution(values, exact)
    full_exact = np.asarray(exact(*grid.coordinates()))
    full_error = np.zeros_like(full_u)
    full_error[grid.active] = full_u[grid.active] - full_exact[grid.active]
    from .metrics import discrete_energy_error

    return {
        "L2h": float(np.sqrt(grid.h**2 * np.sum(e**2))),
        "Linf": float(np.max(np.abs(e))),
        "H1h": discrete_energy_error(full_error, grid.active, grid.h),
        "solution_min": float(np.min(values)),
        "solution_max": float(np.max(values)),
    }
