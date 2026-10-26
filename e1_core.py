"""E1: two right-edge closures x Neumann/Robin, with full verification evidence."""

from __future__ import annotations
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from ..core import assemble_core, matrix_audit, nonsymmetry_measure
from ..solvers import (solve_direct, direct_condition_estimate,
                       direct_fill_statistics, gmres_cross_check)
from ..metrics import error_metrics
from ..verification import local_row_audit, dirichlet_corner_audit
from .common import (save_table, save_figure, add_pairwise_orders,
                     fit_summary, fit_label)


def run_e1(Ns=(8, 16, 32, 64)):
    rows, local_rows, matrix_rows, representative = [], [], [], {}
    for variant in ("neumann", "robin"):
        for closure in ("one_sided", "ghost"):
            for N in Ns:
                start = time.perf_counter()
                system = assemble_core(N, variant, closure)
                assembly_time = time.perf_counter() - start
                direct = solve_direct(system.A, system.b, system.row_kind)
                x, y = system.grid.unknown_coordinates()
                exact = system.data.exact(x, y)
                errors = error_metrics(direct.x, exact, system.grid)
                full = system.grid.full_solution(direct.x, system.data)
                full_exact = system.data.exact(*system.grid.coordinates())
                row = {
                    "experiment": "E1",
                    "variant": variant,
                    "closure": closure,
                    "N": N,
                    "h": system.grid.h,
                    "n_dof": system.A.shape[0],
                    "nnz": system.A.nnz,
                    "assembly_time": assembly_time,
                    "solve_time": direct.wall_time,
                    **errors,
                    **direct.residuals,
                    "nonsymmetry": nonsymmetry_measure(system.A),
                    "full_solution_min": float(np.nanmin(full)),
                    "full_solution_max": float(np.nanmax(full)),
                    "exact_grid_min": float(np.min(full_exact)),
                    "exact_grid_max": float(np.max(full_exact)),
                }
                if N == 32:
                    row["condition_1_estimate"] = direct_condition_estimate(system.A)
                else:
                    row["condition_1_estimate"] = np.nan
                rows.append(row)
                local_rows.append(local_row_audit(system))
                audit = matrix_audit(system)
                corners = dirichlet_corner_audit(system)
                matrix_rows.append({"variant": variant, "closure": closure, "N": N,
                                    **audit, **corners})
                if N == 32:
                    representative[(variant, closure)] = (system, direct)

    frame = pd.DataFrame(rows)
    metrics = ["L2h", "Linf", "L2h_interior", "Linf_interior", "Linf_right"]
    frame = add_pairwise_orders(frame, ["variant", "closure"], metrics)
    save_table(frame, "E1_core_convergence.csv")
    save_table(fit_summary(frame, ["variant", "closure"], metrics),
               "E1_asymptotic_fits.csv")

    local = pd.DataFrame(local_rows)
    local = add_pairwise_orders(local, ["variant", "closure"],
                                ["pde_local_residual_inf", "boundary_local_residual_inf"])
    save_table(local, "E1_local_row_consistency.csv")
    save_table(fit_summary(local, ["variant", "closure"],
                           ["pde_local_residual_inf", "boundary_local_residual_inf"]),
               "E1_local_residual_fits.csv")
    save_table(pd.DataFrame(matrix_rows), "E1_matrix_audits.csv")

    # Independent algebraic audit at a representative mesh.
    algebraic_rows = []
    for (variant, closure), (system, direct) in representative.items():
        x_gmres, info, iterations = gmres_cross_check(system.A, system.b)
        algebraic_rows.append({
            "variant": variant,
            "closure": closure,
            "N": 32,
            "condition_1_estimate": direct_condition_estimate(system.A),
            **direct_fill_statistics(system.A),
            **direct.residuals,
            "gmres_info": info,
            "gmres_iterations": iterations,
            "direct_gmres_difference_inf": float(np.max(np.abs(direct.x - x_gmres))),
            "direct_discretisation_Linf": float(frame.query(
                "variant == @variant and closure == @closure and N == 32"
            )["Linf"].iloc[0]),
        })
    save_table(pd.DataFrame(algebraic_rows), "E1_algebraic_error_audit.csv")

    _plot_matrices()
    _plot_convergence(frame)
    _plot_representative(representative)
    return frame


def _plot_matrices(N=16):
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    for ax, (variant, closure) in zip(axes.ravel(), (
        ("neumann", "one_sided"), ("neumann", "ghost"),
        ("robin", "one_sided"), ("robin", "ghost"),
    )):
        system = assemble_core(N, variant, closure)
        ax.spy(system.A, markersize=1.5, color="navy")
        ax.set_title(f"{variant.title()}, {closure.replace('_', '-')}\n$\\delta_A$={nonsymmetry_measure(system.A):.3e}")
        ax.set_xlabel("column k")
        ax.set_ylabel("row k")
    fig.suptitle("Sparse matrix structure (N=16); all four assemblies are nonsymmetric")
    fig.tight_layout()
    save_figure(fig, "E1_matrix_spy.png")


def _plot_convergence(frame):
    for variant in ("neumann", "robin"):
        subset = frame[frame.variant == variant]
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
        for closure, group in subset.groupby("closure"):
            group = group.sort_values("N")
            closure_label = closure.replace("_", "-")
            axes[0].loglog(group.N, group.L2h,
                           label=fit_label(group, "L2h", f"{closure_label} $L_h^2$"))
            axes[0].loglog(group.N, group.Linf, linestyle="--",
                           label=fit_label(group, "Linf", f"{closure_label} $L_\\infty$"))
            axes[1].loglog(group.N, group.Linf_right,
                           label=fit_label(group, "Linf_right", f"{closure_label} right"))
            axes[1].loglog(group.N, group.Linf_interior, linestyle="--",
                           label=fit_label(group, "Linf_interior", f"{closure_label} interior"))
        x = np.array(sorted(subset.N.unique()), dtype=float)
        # Anchor each reference to the finest measured value of its target curve,
        # not to an arbitrary coarse-grid maximum far above the error data.
        for order, (closure, metric) in ((1, ("one_sided", "Linf")),
                                         (2, ("ghost", "Linf"))):
            target = subset[subset.closure == closure].sort_values("N")
            anchor = float(target[metric].iloc[-1])
            axes[0].loglog(x, anchor * (x / x.max()) ** (-order), marker="none",
                           linestyle=":", label=f"$N^{{-{order}}}$ reference")
        axes[0].set_title("Global errors")
        axes[1].set_title("Interior versus right-boundary error")
        for ax in axes:
            ax.set_xlabel("N (h=1/N)")
            ax.set_ylabel("error")
            ax.set_xticks(x)
            ax.set_xticklabels([str(int(N)) for N in x])
            ax.xaxis.set_minor_formatter(NullFormatter())
        axes[0].legend(fontsize=8.5, loc="upper right")
        axes[1].legend(fontsize=8.5, loc="center right")
        fig.suptitle(f"E1 convergence — {variant.title()} right edge")
        fig.tight_layout()
        save_figure(fig, f"E1_{variant}_convergence.png")

    fig, ax = plt.subplots(figsize=(7, 4.8))
    for (variant, closure), group in frame.groupby(["variant", "closure"]):
        group = group.sort_values("N")
        label = f"{variant.title()}, {closure.replace('_', '-')}"
        ax.loglog(group.N, group.nonsymmetry, label=label)
    ax.set(xlabel="N", ylabel=r"$\|A-A^T\|_F/\|A\|_F$",
           title="Measured matrix nonsymmetry")
    x = np.array(sorted(frame.N.unique()), dtype=float)
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(N)) for N in x])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.legend(fontsize=9)
    fig.tight_layout()
    save_figure(fig, "E1_nonsymmetry.png")


def _plot_representative(representative):
    # One solution surface and four signed-error maps at exactly the same mesh.
    system, direct = representative[("neumann", "ghost")]
    full = system.grid.full_solution(direct.x, system.data)
    exact = system.data.exact(*system.grid.coordinates())
    x, y = system.grid.coordinates()
    fig = plt.figure(figsize=(11, 4.8))
    for index, (z, title) in enumerate(((full, "Numerical solution"), (exact, "Exact $e^{x+y}$")), 1):
        ax = fig.add_subplot(1, 2, index, projection="3d")
        ax.plot_surface(x, y, z, cmap="viridis", linewidth=0, antialiased=True)
        ax.set(xlabel="x", ylabel="y", zlabel="u", title=title)
    fig.suptitle("Representative ghost-Neumann solution (N=32)")
    save_figure(fig, "E1_solution_surface.png")

    fig, axes = plt.subplots(2, 2, figsize=(10, 8.5))
    for ax, ((variant, closure), (system, direct)) in zip(axes.ravel(), representative.items()):
        full = system.grid.full_solution(direct.x, system.data)
        exact = system.data.exact(*system.grid.coordinates())
        error = np.abs(full - exact)
        x, y = system.grid.coordinates()
        image = ax.pcolormesh(x, y, np.maximum(error, 1e-18), shading="auto", cmap="magma")
        ax.set_aspect("equal")
        ax.set_title(f"{variant.title()}, {closure.replace('_', '-')}; max={error.max():.2e}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(image, ax=ax, label="|error|")
    fig.suptitle("Nodal error distribution (N=32)")
    fig.tight_layout()
    save_figure(fig, "E1_error_maps.png")
