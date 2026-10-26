"""Independent manufactured-solution and row-by-row assembly audits."""

from __future__ import annotations
import numpy as np

from .core import CoreSystem


def exact_local_residuals(system: CoreSystem):
    """Return assembled and independently derived local residuals at u*.

    The residual is A*u*-b.  For ghost rows this is the *entire eliminated PDE
    row*, not merely the centred derivative condition.  Its leading unscaled
    term is O(h), even though the centred boundary derivative itself is O(h^2).
    """
    grid = system.grid
    h = grid.h
    x, y = grid.unknown_coordinates()
    u = system.data.exact(x, y)
    assembled = system.A @ u - system.b
    expected = np.empty_like(u)
    i, _ = grid.unknown_ij()
    pde = i < grid.N
    right = i == grid.N

    # For u=exp(x+y), each directional second difference contributes
    # 2u(cosh(h)-1)/h^2.  Thus -Delta_h(u)-f has the closed form below.
    pde_residual = 2.0 * u * (1.0 - 2.0 * (np.cosh(h) - 1.0) / h**2)
    expected[pde] = pde_residual[pde]

    if system.closure == "one_sided":
        # Robin cancels exactly because r=q+beta*u at the exact solution.
        expected[right] = u[right] * ((1.0 - np.exp(-h)) / h - 1.0)
    else:
        # Full ghost-eliminated boundary PDE row: interior truncation plus the
        # centred-normal-derivative defect multiplied by 2/h.
        expected[right] = pde_residual[right] + (
            2.0 * u[right] / h * (np.sinh(h) / h - 1.0)
        )
    return assembled, expected


def local_row_audit(system: CoreSystem):
    assembled, expected = exact_local_residuals(system)
    i, _ = system.grid.unknown_ij()
    right = i == system.grid.N
    pde = ~right
    return {
        "N": system.grid.N,
        "variant": system.variant,
        "closure": system.closure,
        "max_independent_formula_difference": float(
            np.max(np.abs(assembled - expected))
        ),
        "pde_local_residual_inf": float(np.max(np.abs(assembled[pde]))),
        "boundary_local_residual_inf": float(np.max(np.abs(assembled[right]))),
        "boundary_local_residual_expected_inf": float(
            np.max(np.abs(expected[right]))
        ),
    }


def dirichlet_corner_audit(system: CoreSystem):
    grid = system.grid
    corners = ((0, 0), (grid.N, 0), (0, grid.N), (grid.N, grid.N))
    return {
        "all_corners_dirichlet": all(grid.k_of[i, j] < 0 for i, j in corners),
        "right_corners_owned_by_dirichlet": bool(
            grid.node_type[grid.N, 0] == 3 and grid.node_type[grid.N, grid.N] == 3
        ),
        "unknown_count_matches_formula": grid.n_unknowns == grid.N * (grid.N - 1),
        "index_map_is_permutation": np.array_equal(
            np.sort(grid.k_of[grid.k_of >= 0]), np.arange(grid.n_unknowns)
        ),
    }
