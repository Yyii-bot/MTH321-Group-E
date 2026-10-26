"""E3: matched two-dimensional FDM versus P1 FEM on pure Dirichlet data."""

from __future__ import annotations
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..data import (square_exact, square_forcing, sine_exact, sine_forcing)
from ..dirichlet import square_grid, assemble_five_point_dirichlet
from ..fem import unit_square_triangulation, fem_dirichlet_solve, fem_error_integrals
from ..solvers import solve_direct
from .common import (save_table, save_figure, add_pairwise_orders, fit_summary,
                     fit_label)


def run_e3(Ns=(8, 16, 32, 64)):
    rows, stencil_rows = [], []
    data_sets = {
        "exponential": (square_exact, square_forcing),
        "sine": (sine_exact, sine_forcing),
    }
    for data_name, (exact, forcing) in data_sets.items():
        for N in Ns:
            h = 1.0 / N
            start = time.perf_counter()
            grid = square_grid(N)
            fd_system = assemble_five_point_dirichlet(grid, forcing, exact)
            fd_assembly = time.perf_counter() - start
            fd = solve_direct(fd_system.A, fd_system.b)
            i, j = grid.unknown_ij()
            common_exact = exact(grid.x[i], grid.y[j])
            fd_error = fd.x - common_exact

            mesh = unit_square_triangulation(N)
            fe = fem_dirichlet_solve(mesh, forcing, exact, intorder=6)
            fem_xy = mesh.p[:, fe.interior]
            if not (np.allclose(fem_xy[0], grid.x[i]) and np.allclose(fem_xy[1], grid.y[j])):
                raise RuntimeError("FDM/FEM common-node ordering does not match")
            if data_name == "exponential":
                # On these right-isosceles triangles the hypotenuse coupling is
                # zero: P1 stiffness and h^2 times the five-point stencil coincide.
                difference = fd_system.A - (N**2) * fe.A
                fd_pattern = fd_system.A.copy()
                fe_pattern = fe.A.copy()
                fd_pattern.data[:] = 1.0
                fe_pattern.data[:] = 1.0
                pattern_difference = fd_pattern - fe_pattern
                pattern_difference.eliminate_zeros()
                stencil_rows.append({
                    "N": N,
                    "n_dof": fd_system.A.shape[0],
                    "fdm_nnz": fd_system.A.nnz,
                    "fem_nnz": fe.A.nnz,
                    "pattern_mismatch_count": pattern_difference.nnz,
                    "relative_stencil_difference": float(np.sqrt(
                        difference.power(2).sum() / fd_system.A.power(2).sum()
                    )),
                })
            fe_error = fe.values[fe.interior] - common_exact
            integral = fem_error_integrals(mesh, fe.values, exact, intorder=8)

            for method, A, assembly, solve, residuals, error, nnz in (
                ("FDM", fd_system.A, fd_assembly, fd.wall_time, fd.residuals,
                 fd_error, fd_system.A.nnz),
                ("FEM_P1", fe.A, fe.assembly_time, fe.solve_time, fe.residuals,
                 fe_error, fe.A.nnz),
            ):
                rows.append({
                    "data": data_name, "method": method, "N": N, "h": h,
                    "n_dof": A.shape[0], "nnz": nnz,
                    "assembly_time": assembly, "solve_time": solve,
                    **residuals,
                    "L2h": float(np.sqrt(h**2 * np.sum(error**2))),
                    "Linf_common": float(np.max(np.abs(error))),
                    "FEM_integral_L2": integral["L2"] if method == "FEM_P1" else np.nan,
                })
            print(f"E3 completed {data_name}, N={N}", flush=True)
    frame = pd.DataFrame(rows)
    frame = add_pairwise_orders(frame, ["data", "method"], ["L2h", "Linf_common"])
    save_table(frame, "E3_matched_FDM_FEM.csv")
    save_table(fit_summary(frame, ["data", "method"], ["L2h", "Linf_common"]),
               "E3_asymptotic_fits.csv")
    save_table(pd.DataFrame(stencil_rows), "E3_stencil_equivalence.csv")
    _plots(frame)
    return frame


def _plots(frame):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, (data, subset) in zip(axes, frame.groupby("data")):
        for method, group in subset.groupby("method"):
            group = group.sort_values("N")
            ax.loglog(group.N, group.L2h,
                      label=fit_label(group, "L2h", f"{method} $L2_h$"))
            ax.loglog(group.N, group.Linf_common, linestyle="--",
                      label=fit_label(group, "Linf_common", f"{method} $L_\\infty$"))
        x = np.array(sorted(subset.N.unique()), dtype=float)
        fdm = subset[subset.method == "FDM"].sort_values("N")
        anchor = float(fdm.L2h.iloc[-1])
        ax.loglog(x, anchor * (x / x.max()) ** -2, marker="none", linestyle=":",
                  label="$N^{-2}$ reference anchored to finest FDM $L2_h$")
        ax.set(title=f"{data} manufactured data", xlabel="N", ylabel="error")
        ax.legend(fontsize=7)
    fig.tight_layout()
    save_figure(fig, "E3_error_comparison.png")

    subset = frame[frame.data == "exponential"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for method, group in subset.groupby("method"):
        group = group.sort_values("N")
        axes[0].loglog(group.N, group.assembly_time, label=method)
        axes[1].loglog(group.N, group.solve_time, label=method)
        axes[2].loglog(group.N, group.nnz, label=method)
    for ax, title, ylabel in zip(axes,
                                 ("Assembly time", "Solve time", "Matrix nonzeros"),
                                 ("seconds", "seconds", "nnz")):
        ax.set(title=title, xlabel="N", ylabel=ylabel)
        ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, "E3_cost_and_sparsity.png")
