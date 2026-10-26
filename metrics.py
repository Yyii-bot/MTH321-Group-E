"""Named algebraic residuals, discrete error norms, and observed-order fits."""

from __future__ import annotations
import numpy as np
import scipy.sparse as sp


def row_scaling(A: sp.spmatrix):
    """D_s=diag(sum_j |A_ij|), a nonzero scale for every assembled row."""
    scale = np.asarray(abs(A).sum(axis=1)).ravel()
    if np.any(scale == 0.0):
        raise ValueError("Cannot scale a matrix with an empty row")
    return scale


def residual_metrics(A, b, x, row_kind=None):
    """Separate algebraic residual measures from discretisation error.

    scaled_relative_inf is the project metric
        ||D_s^-1(b-Ax)||_inf / ||D_s^-1 b||_inf.
    componentwise_backward is max_i |r_i|/(|b_i|+(|A||x|)_i).
    """
    A = A.tocsr()
    b = np.asarray(b, dtype=float)
    x = np.asarray(x, dtype=float)
    r = b - A @ x
    scale = row_scaling(A)
    scaled_r = np.abs(r) / scale
    scaled_b = np.abs(b) / scale
    rhs_norm = max(float(np.max(scaled_b)), np.finfo(float).tiny)
    denominator = np.abs(b) + np.asarray(abs(A) @ np.abs(x)).ravel()
    component = np.zeros_like(r)
    np.divide(np.abs(r), denominator, out=component, where=denominator > 0.0)
    result = {
        "residual_l2": float(np.linalg.norm(r, 2)),
        "relative_residual_l2": float(
            np.linalg.norm(r, 2) / max(np.linalg.norm(b, 2), np.finfo(float).tiny)
        ),
        "scaled_residual_inf": float(np.max(scaled_r)),
        "scaled_relative_inf": float(np.max(scaled_r) / rhs_norm),
        "componentwise_backward": float(np.max(component)),
    }
    if row_kind is not None:
        row_kind = np.asarray(row_kind)
        for label, mask in (
            ("pde", row_kind == 0),
            ("boundary", row_kind == 1),
        ):
            result[f"scaled_{label}_inf"] = (
                float(np.max(scaled_r[mask])) if np.any(mask) else 0.0
            )
    return result


def error_metrics(unknown_values, exact_values, grid):
    """Discrete errors on the core unknown set, with interior/right splitting.

    Global and interior L2h norms use weight h^2 per included unknown node.
    The right-edge L2 norm additionally uses the line weight h and is named
    L2gamma to avoid presenting it as an area norm.
    """
    e = np.asarray(unknown_values) - np.asarray(exact_values)
    i, _ = grid.unknown_ij()
    interior = i <= grid.N - 1
    right = i == grid.N
    h2 = grid.h**2

    def l2h(mask):
        return float(np.sqrt(h2 * np.sum(e[mask] ** 2))) if np.any(mask) else 0.0

    def linf(mask):
        return float(np.max(np.abs(e[mask]))) if np.any(mask) else 0.0

    return {
        "L2h": l2h(np.ones(e.shape, dtype=bool)),
        "Linf": float(np.max(np.abs(e))),
        "L2h_interior": l2h(interior),
        "Linf_interior": linf(interior),
        "Linf_right": linf(right),
        "L2gamma_right": float(np.sqrt(grid.h * np.sum(e[right] ** 2))),
        "solution_min": float(np.min(unknown_values)),
        "solution_max": float(np.max(unknown_values)),
    }


def pairwise_orders(Ns, errors):
    """p_k=ln(E_{k-1}/E_k)/ln(N_k/N_{k-1}); undefined at zero error."""
    Ns = np.asarray(Ns, dtype=float)
    errors = np.asarray(errors, dtype=float)
    orders = np.full(errors.shape, np.nan)
    for k in range(1, len(errors)):
        if errors[k - 1] > 0.0 and errors[k] > 0.0:
            orders[k] = np.log(errors[k - 1] / errors[k]) / np.log(Ns[k] / Ns[k - 1])
    return orders


def asymptotic_fit(Ns, errors, last=3):
    """Fit E=C N^(-p) over the declared asymptotic range (last=3 by default)."""
    Ns = np.asarray(Ns, dtype=float)
    errors = np.asarray(errors, dtype=float)
    mask = np.isfinite(errors) & (errors > 0.0)
    indices = np.flatnonzero(mask)[-last:]
    if len(indices) < 2:
        return {"order": np.nan, "R2": np.nan, "max_log_deviation": np.nan, "fit_Ns": []}
    lx = np.log(Ns[indices])
    ly = np.log(errors[indices])
    slope, intercept = np.polyfit(lx, ly, 1)
    predicted = slope * lx + intercept
    ss_res = float(np.sum((ly - predicted) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    r2 = 1.0 if ss_tot == 0.0 and ss_res == 0.0 else (
        1.0 - ss_res / ss_tot if ss_tot > 0.0 else np.nan
    )
    return {
        "order": float(-slope),
        "R2": float(r2),
        "max_log_deviation": float(np.max(np.abs(ly - predicted))),
        "fit_Ns": [int(Ns[k]) for k in indices],
    }


def discrete_energy_error(full_error, active, h):
    """Discrete H1 seminorm on unique active east/north grid edges.

    For each edge, h^2*((e_Q-e_P)/h)^2=(e_Q-e_P)^2.  This uses
    differences of the numerical error, never the singular exact gradient.
    """
    full_error = np.asarray(full_error, dtype=float)
    active = np.asarray(active, dtype=bool)
    east = active[:-1, :] & active[1:, :]
    north = active[:, :-1] & active[:, 1:]
    east_difference = full_error[1:, :] - full_error[:-1, :]
    north_difference = full_error[:, 1:] - full_error[:, :-1]
    total = (float(np.sum(east_difference[east] ** 2))
             + float(np.sum(north_difference[north] ** 2)))
    return float(np.sqrt(total))
