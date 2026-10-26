"""Sparse finite-difference assembly for the mixed-boundary core problem."""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

from .data import CORE_DATA, CoreData


class NodeType(IntEnum):
    DIRICHLET = 0
    INTERIOR = 1
    RIGHT_BOUNDARY = 2
    RIGHT_CORNER = 3


@dataclass(frozen=True)
class MixedGrid:
    N: int
    h: float
    k_of: np.ndarray
    node_type: np.ndarray

    @property
    def n_unknowns(self):
        return self.N * (self.N - 1)

    def coordinates(self):
        i, j = np.indices((self.N + 1, self.N + 1))
        return i * self.h, j * self.h

    def unknown_ij(self):
        ks = np.arange(self.n_unknowns)
        # k=(j-1)N+(i-1): x index varies fastest within each row.
        i = ks % self.N + 1
        j = ks // self.N + 1
        return i, j

    def unknown_coordinates(self):
        i, j = self.unknown_ij()
        return i * self.h, j * self.h

    def red_black_colours(self):
        i, j = self.unknown_ij()
        return ((i + j) % 2).astype(np.int8)

    def full_solution(self, unknown_values, data: CoreData = CORE_DATA):
        x, y = self.coordinates()
        full = data.dirichlet(x, y)
        full = np.array(full, dtype=float, copy=True)
        i, j = self.unknown_ij()
        full[i, j] = unknown_values
        return full


@dataclass
class CoreSystem:
    A: sp.csr_matrix
    b: np.ndarray
    grid: MixedGrid
    variant: str
    closure: str
    row_kind: np.ndarray  # 0=PDE, 1=right-boundary closure
    data: CoreData

    @property
    def colours(self):
        return self.grid.red_black_colours()


def index_map(N: int):
    """Build the explicit (i,j)->k map and protect all Dirichlet corners.

    Unknowns are {(i,j): 1<=i<=N, 1<=j<=N-1}; therefore there are
    N(N-1) unknowns.  In particular, (N,0) and (N,N) can never receive a
    Neumann/Robin row.
    """
    if not isinstance(N, (int, np.integer)) or N < 2:
        raise ValueError("N must be an integer >= 2")
    k_of = -np.ones((N + 1, N + 1), dtype=np.int32)
    node_type = np.full((N + 1, N + 1), NodeType.DIRICHLET, dtype=np.int8)
    for j in range(1, N):
        for i in range(1, N + 1):
            k = (j - 1) * N + (i - 1)
            k_of[i, j] = k
            node_type[i, j] = (
                NodeType.RIGHT_BOUNDARY if i == N else NodeType.INTERIOR
            )
    node_type[N, 0] = NodeType.RIGHT_CORNER
    node_type[N, N] = NodeType.RIGHT_CORNER
    h = 1.0 / N
    return MixedGrid(N, h, k_of, node_type)


def _is_dirichlet(grid: MixedGrid, i: int, j: int):
    return grid.k_of[i, j] < 0


def assemble_core(
    N: int,
    bc_variant: str,
    closure: str,
    data: CoreData = CORE_DATA,
):
    """Assemble one prescribed mixed-boundary system in COO, then convert CSR.

    closure='one_sided' replaces the right-edge PDE row by a derivative row.
    closure='ghost' retains the boundary PDE row and eliminates U_{N+1,j}.
    Dirichlet values are eliminated, never inserted as unrelated identity rows.
    """
    variant = bc_variant.lower()
    closure = closure.lower()
    if variant not in {"neumann", "robin"}:
        raise ValueError("bc_variant must be 'neumann' or 'robin'")
    if closure not in {"one_sided", "ghost"}:
        raise ValueError("closure must be 'one_sided' or 'ghost'")

    grid = index_map(N)
    h = grid.h
    inv_h2 = 1.0 / h**2
    rows, cols, values = [], [], []
    b = np.zeros(grid.n_unknowns, dtype=float)
    row_kind = np.zeros(grid.n_unknowns, dtype=np.int8)

    def add(row, col, value):
        rows.append(row)
        cols.append(col)
        values.append(float(value))

    for j in range(1, N):
        y = j * h
        for i in range(1, N + 1):
            x = i * h
            k = grid.k_of[i, j]

            if i < N:
                # Positive -Delta_h: diagonal positive, neighbour entries <=0.
                add(k, k, 4.0 * inv_h2)
                b[k] = data.forcing(x, y)
                for ii, jj in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)):
                    kk = grid.k_of[ii, jj]
                    if kk >= 0:
                        add(k, kk, -inv_h2)
                    else:
                        b[k] += inv_h2 * data.dirichlet(ii * h, jj * h)
                continue

            # Right-edge unknown, 1<=j<=N-1.  Corners are structurally absent.
            row_kind[k] = 1
            k_west = grid.k_of[N - 1, j]
            if closure == "one_sided":
                add(k, k_west, -1.0 / h)
                if variant == "neumann":
                    add(k, k, 1.0 / h)
                    b[k] = data.neumann_right(y)
                else:
                    add(k, k, 1.0 / h + data.beta)
                    b[k] = data.robin_right(y)
            else:
                diagonal = 4.0 + (2.0 * h * data.beta if variant == "robin" else 0.0)
                add(k, k, diagonal * inv_h2)
                add(k, k_west, -2.0 * inv_h2)
                b[k] = data.forcing(x, y)
                b[k] += (
                    2.0 * data.neumann_right(y) / h
                    if variant == "neumann"
                    else 2.0 * data.robin_right(y) / h
                )
                # East and west are already accounted for by ghost elimination.
                for ii, jj in ((N, j - 1), (N, j + 1)):
                    kk = grid.k_of[ii, jj]
                    if kk >= 0:
                        add(k, kk, -inv_h2)
                    else:
                        b[k] += inv_h2 * data.dirichlet(ii * h, jj * h)

    A_coo = sp.coo_matrix(
        (values, (rows, cols)), shape=(grid.n_unknowns, grid.n_unknowns)
    )
    A = A_coo.tocsr()
    A.sum_duplicates()
    return CoreSystem(A, b, grid, variant, closure, row_kind, data)


def nonsymmetry_measure(A: sp.spmatrix):
    """Frobenius ||A-A^T||/||A|| without densifying the matrix."""
    A = A.tocsr()
    numerator = (A - A.T).power(2).sum()
    denominator = A.power(2).sum()
    return float(np.sqrt(numerator / denominator))


def matrix_audit(system: CoreSystem):
    """Return structural evidence, including sign pattern and row sums."""
    A = system.A
    diagonal = A.diagonal()
    off_diagonal = A - sp.diags(diagonal, format="csr")
    row_sums = np.asarray(A.sum(axis=1)).ravel()
    components = connected_components(A, directed=False, return_labels=False)
    return {
        "n": A.shape[0],
        "nnz": A.nnz,
        "min_diagonal": float(diagonal.min()),
        "max_off_diagonal": float(np.max(off_diagonal.data)) if off_diagonal.nnz else 0.0,
        "min_row_sum": float(row_sums.min()),
        "max_row_sum": float(row_sums.max()),
        "connected_components": int(components),
        "nonsymmetry": nonsymmetry_measure(A),
        "right_corner_unknowns": int(
            (system.grid.k_of[system.grid.N, 0] >= 0)
            + (system.grid.k_of[system.grid.N, system.grid.N] >= 0)
        ),
    }
