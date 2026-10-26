"""E6 (optional): pure-Neumann compatibility defect, nullspace, and gauge."""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..fem import unit_square_triangulation, fem_error_integrals
from ..pure_neumann import fem_pure_neumann_solve, small_neumann_spectrum
from .common import (save_table, save_figure, add_pairwise_orders, fit_summary,
                     fit_label)


def run_e6(Ns=(4, 8, 16, 32)):
    rows = []
    for order in (2, 8):
        for N in Ns:
            mesh = unit_square_triangulation(N)
            result = fem_pure_neumann_solve(mesh, intorder=order)
            errors = fem_error_integrals(mesh, result["values"],
                                         lambda x, y: np.exp(x + y), intorder=8)
            rows.append({
                "N": N, "quadrature_order": order,
                "continuous_compatibility_defect": 0.0,
                **{key: value for key, value in result.items()
                   if key not in {"values", "A", "b", "saddle_residual"}},
                **{f"saddle_{key}": value for key, value in result["saddle_residual"].items()},
                **errors,
            })
            print(f"E6 completed quadrature={order}, N={N}", flush=True)
    frame = pd.DataFrame(rows)
    frame = add_pairwise_orders(frame, ["quadrature_order"],
                                ["L2", "Linf"])
    save_table(frame, "E6_pure_Neumann_audit.csv")
    save_table(fit_summary(frame, ["quadrature_order"], ["L2", "Linf"]),
               "E6_asymptotic_fits.csv")
    spectra = []
    for N in (4, 8):
        spectra.append({"N": N, **small_neumann_spectrum(unit_square_triangulation(N))})
    save_table(pd.DataFrame(spectra), "E6_small_mesh_nullspace.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for order, group in frame.groupby("quadrature_order"):
        group = group.sort_values("N")
        axes[0].loglog(group.N, np.maximum(np.abs(group.defect_before), 1e-18),
                       label=f"quadrature order {order}")
        axes[1].loglog(group.N, group.L2,
                       label=fit_label(group, "L2", f"quad {order}: $L^2$"))
        axes[1].loglog(group.N, group.Linf, linestyle="--",
                       label=fit_label(group, "Linf", f"quad {order}: $L_\\infty$"))
    axes[0].set(xlabel="N", ylabel="|1^T b| before correction",
                title="Discrete compatibility defect")
    axes[1].set(xlabel="N", ylabel="error", title="Mean-gauged solution error")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, "E6_compatibility_and_convergence.png")
    return frame
