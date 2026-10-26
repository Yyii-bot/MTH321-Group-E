"""E4: same aligned FDM solver, smooth versus re-entrant-singular L-shaped data."""

from __future__ import annotations
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..data import harmonic_exact, singular_exact
from ..dirichlet import lshape_grid, assemble_five_point_dirichlet, grid_error_metrics
from ..core import nonsymmetry_measure
from ..solvers import solve_direct
from .common import (save_table, save_figure, add_pairwise_orders, fit_summary,
                     fit_label)


def run_e4(Ns=(32, 64, 128, 256)):
    rows, representative = [], {}
    representative_N = max(Ns)
    zero = lambda x, y: np.zeros_like(np.asarray(x, dtype=float))
    for name, exact in (("smooth_harmonic", harmonic_exact),
                        ("singular_r_2_3", singular_exact)):
        for N in Ns:
            start = time.perf_counter()
            grid = lshape_grid(N)
            system = assemble_five_point_dirichlet(grid, zero, exact)
            assembly_time = time.perf_counter() - start
            direct = solve_direct(system.A, system.b)
            errors = grid_error_metrics(system, direct.x, exact)
            i, j = grid.unknown_ij()
            nodal_error = np.abs(direct.x - exact(grid.x[i], grid.y[j]))
            furthest = int(np.argmax(nodal_error))
            rows.append({
                "solution": name, "N": N, "h": grid.h,
                "n_dof": system.A.shape[0], "nnz": system.A.nnz,
                "assembly_time": assembly_time, "solve_time": direct.wall_time,
                **errors, **direct.residuals,
                "nonsymmetry": nonsymmetry_measure(system.A),
                "max_error_x": float(grid.x[i[furthest]]),
                "max_error_y": float(grid.y[j[furthest]]),
                "max_error_radius": float(np.hypot(grid.x[i[furthest]], grid.y[j[furthest]])),
            })
            if N == representative_N:
                representative[name] = (system, direct, exact)
            print(f"E4 completed {name}, N={N}", flush=True)
    frame = pd.DataFrame(rows)
    frame = add_pairwise_orders(frame, ["solution"], ["L2h", "Linf", "H1h"])
    save_table(frame, "E4_Lshape_smooth_vs_singular.csv")
    save_table(fit_summary(frame, ["solution"], ["L2h", "Linf", "H1h"]),
               "E4_asymptotic_fits.csv")
    _plots(frame, representative)
    return frame


def _plots(frame, representative):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    for solution, group in frame.groupby("solution"):
        group = group.sort_values("N")
        for ax, metric in zip(axes, ("L2h", "Linf", "H1h")):
            ax.loglog(group.N, group[metric],
                      label=fit_label(group, metric, solution))
    for ax, metric in zip(axes, ("L2h", "Linf", "H1h")):
        ax.set(xlabel="N (h=2/N)", ylabel=metric, title=f"{metric} error")
        ax.legend(fontsize=8)
    fig.suptitle("E4: measured FDM convergence on the L-shaped domain")
    fig.tight_layout()
    save_figure(fig, "E4_convergence.png")

    singular = frame[frame.solution == "singular_r_2_3"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    for ax, metric in zip(axes, ("H1h", "L2h", "Linf")):
        ax.loglog(singular.n_dof, singular[metric],
                  label=fit_label(singular, metric, "measured", x_column="n_dof"))
        ax.set(xlabel="number of unknowns", ylabel=metric,
               title=f"Singular {metric} versus DOFs")
    axes[0].loglog(singular.n_dof, singular.H1h.iloc[-1] *
                   (singular.n_dof / singular.n_dof.iloc[-1]) ** (-1 / 3),
                   marker="none", linestyle="--", label="$N_{dof}^{-1/3}$")
    axes[1].loglog(singular.n_dof, singular.L2h.iloc[-1] *
                   (singular.n_dof / singular.n_dof.iloc[-1]) ** (-2 / 3),
                   marker="none", linestyle="--", label="$N_{dof}^{-2/3}$")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, "E4_error_vs_dofs.png")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, (name, (system, direct, exact)) in zip(axes, representative.items()):
        grid = system.grid
        full = grid.full_solution(direct.x, exact)
        full_exact = exact(*grid.coordinates())
        error = np.abs(full - full_exact)
        x, y = grid.coordinates()
        image = ax.pcolormesh(x, y, np.maximum(error, 1e-18), shading="auto", cmap="magma")
        ax.plot([0], [0], marker="+", color="cyan", markersize=12, markeredgewidth=2)
        ax.set_aspect("equal")
        ax.set(xlabel="x", ylabel="y",
               title=f"{name}: nodal |error| (N={max(frame.N)})")
        fig.colorbar(image, ax=ax)
    fig.tight_layout()
    save_figure(fig, "E4_error_maps.png")
