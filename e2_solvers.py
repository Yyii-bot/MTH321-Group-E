"""E2: direct versus Jacobi/Gauss-Seidel/SOR at matched scaled residuals."""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..core import assemble_core
from ..metrics import error_metrics
from ..solvers import (solve_direct, solve_jacobi, solve_gauss_seidel,
                       solve_red_black_sor, direct_fill_statistics)
from .common import save_table, save_figure


_METHOD_LABELS = {
    "direct": "SuperLU",
    "jacobi": "Jacobi",
    "gauss_seidel_rb": "red-black GS",
    "sor_rb": "red-black SOR",
}


def _right_midpoint_qoi(system, unknown_values):
    """Main physical quantity: numerical temperature u_h(1, 1/2)."""
    N = system.grid.N
    if N % 2:
        raise ValueError("The right-midpoint QoI requires an even N")
    k = system.grid.k_of[N, N // 2]
    if k < 0:
        raise RuntimeError("The right midpoint must be an unknown")
    exact = float(system.data.exact(1.0, 0.5))
    return float(unknown_values[k]), exact


def _iteration_matrix(A, colours, method, omega=1.0):
    """Dense red-black iteration matrix for small-mesh spectral audits only."""
    A = A.tocsr()
    red = np.flatnonzero(colours == 0)
    black = np.flatnonzero(colours == 1)
    perm = np.concatenate((red, black))
    Ap = A[perm][:, perm].toarray()
    nr = len(red)
    DR = np.diag(Ap[:nr, :nr])
    DB = np.diag(Ap[nr:, nr:])
    ARB = Ap[:nr, nr:]
    ABR = Ap[nr:, :nr]
    DRi, DBi = 1.0 / DR, 1.0 / DB
    n = A.shape[0]
    G = np.zeros((n, n))
    if method == "jacobi":
        return np.eye(n) - Ap / np.diag(Ap)[:, None]
    if method == "gauss_seidel":
        G[:nr, nr:] = -DRi[:, None] * ARB
        G[nr:, nr:] = DBi[:, None] * (ABR * DRi[None, :]) @ ARB
    elif method == "sor":
        w = omega
        G[:nr, :nr] = (1.0 - w) * np.eye(nr)
        G[:nr, nr:] = -w * DRi[:, None] * ARB
        G[nr:, :nr] = -w * (1.0 - w) * DBi[:, None] * ABR
        G[nr:, nr:] = ((1.0 - w) * np.eye(n - nr)
                        + w**2 * DBi[:, None] * (ABR * DRi[None, :]) @ ARB)
    else:
        raise ValueError(method)
    return G


def run_e2(Ns=(16, 32, 64, 128), extended_N=256, tolerances=(1e-10, 1e-6),
           accuracy_tolerances=(1e-10, 1e-11, 1e-12, 1e-13), accuracy_target=0.01):
    rows, omega_rows, accuracy_rows = [], [], []
    result_cache = {}
    all_Ns = tuple(Ns) + (() if extended_N is None else (extended_N,))
    accuracy_N = max(Ns)
    for N in all_Ns:
        system = assemble_core(N, "neumann", "ghost")
        x, y = system.grid.unknown_coordinates()
        exact = system.data.exact(x, y)
        direct = solve_direct(system.A, system.b, system.row_kind)
        direct_error = error_metrics(direct.x, exact, system.grid)
        direct_row = {
            "N": N, "method": "direct", "tolerance": np.nan,
            "iterations": 0, "converged": True,
            "wall_time": direct.wall_time,
            **direct.residuals,
            **direct_error,
            **direct_fill_statistics(system.A),
            "omega": np.nan,
            "difference_from_direct_inf": 0.0,
            "algebraic_to_discretisation": 0.0,
        }
        rows.append(direct_row)

        # The model formula is only an initial guess; every omega is measured.
        omega_star = 2.0 / (1.0 + np.sin(np.pi / N))
        candidates = np.linspace(max(0.05, omega_star - 0.05),
                                 min(1.999, omega_star + 0.03), 9)
        best = None
        for omega in candidates:
            trial = solve_red_black_sor(system.A, system.b, system.colours,
                                        omega=omega, tol=1e-8, maxiter=100000,
                                        row_kind=system.row_kind)
            omega_rows.append({"N": N, "omega": omega,
                               "iterations": trial.iterations,
                               "converged": trial.converged,
                               "wall_time": trial.wall_time})
            if trial.converged and (best is None or trial.iterations < best[0]):
                best = (trial.iterations, float(omega))
        omega = best[1] if best else omega_star

        for tolerance in tolerances:
            solvers = [
                ("jacobi", lambda: solve_jacobi(
                    system.A, system.b, tol=tolerance, maxiter=500000,
                    row_kind=system.row_kind)),
                ("gauss_seidel_rb", lambda: solve_gauss_seidel(
                    system.A, system.b, system.colours, tol=tolerance,
                    maxiter=500000, row_kind=system.row_kind)),
                ("sor_rb", lambda: solve_red_black_sor(
                    system.A, system.b, system.colours, omega=omega,
                    tol=tolerance, maxiter=100000, row_kind=system.row_kind)),
            ]
            # At N=256 the standard sweep extends direct/SOR; the two slowest
            # methods remain available through run_e2(extended_N=None).
            if extended_N is not None and N == extended_N and extended_N not in Ns:
                solvers = solvers[2:]
            for method, call in solvers:
                result = call()
                errors = error_metrics(result.x, exact, system.grid)
                difference = float(np.max(np.abs(result.x - direct.x)))
                ratio = difference / direct_error["Linf"]
                result_cache[(N, method, float(tolerance))] = (result, difference, ratio)
                rows.append({
                    "N": N, "method": method, "tolerance": tolerance,
                    "iterations": result.iterations,
                    "converged": result.converged,
                    "wall_time": result.wall_time,
                    **result.residuals, **errors,
                    "omega": result.omega,
                    "difference_from_direct_inf": difference,
                    "algebraic_to_discretisation": ratio,
                })
        if N == accuracy_N:
            accuracy_rows.extend(_accuracy_target_sweep(
                system, direct, direct_error, omega, accuracy_tolerances,
                target=accuracy_target, cache=result_cache,
            ))
        print(f"E2 completed N={N}", flush=True)

    frame = pd.DataFrame(rows)
    save_table(frame, "E2_solver_comparison.csv")
    omega_frame = pd.DataFrame(omega_rows)
    save_table(omega_frame, "E2_omega_scan.csv")
    _spectral_audit()
    _stopping_audit()
    _crossover_summary(frame)
    _plots(frame, omega_frame)
    accuracy_frame = pd.DataFrame(accuracy_rows)
    save_table(accuracy_frame, "E2_accuracy_target_sweep.csv")
    accuracy_summary = _accuracy_target_summary(accuracy_frame, frame)
    _accuracy_target_plot(accuracy_frame)
    _stagnation_plot(accuracy_frame)
    _matched_accuracy_cost_plot(accuracy_summary)
    return frame


def _accuracy_target_sweep(system, direct, direct_error, omega, tolerances,
                           target=0.01, cache=None):
    """Tighten tolerances until algebraic error, reported errors, and QoI stabilise.

    The main quantity of interest is u_h(1,1/2).  Successive changes in L2h,
    Linf, and this QoI are normalised by the corresponding direct-solve
    discretisation-error scales.  "Stable" therefore has a declared numerical
    meaning: every normalised change is at most ``target`` (1% by default).
    """
    tolerances = sorted((float(tol) for tol in tolerances), reverse=True)
    if len(tolerances) < 2:
        raise ValueError("At least two tolerances are needed for a stagnation audit")
    x, y = system.grid.unknown_coordinates()
    exact = system.data.exact(x, y)
    direct_qoi, qoi_exact = _right_midpoint_qoi(system, direct.x)
    direct_qoi_error = direct_qoi - qoi_exact
    qoi_scale = max(abs(direct_qoi_error),
                    np.finfo(float).eps * max(1.0, abs(qoi_exact)))
    rows = []
    previous = {}
    for tolerance in tolerances:
        solvers = [
            ("jacobi", lambda: solve_jacobi(
                system.A, system.b, tol=tolerance, maxiter=500000,
                row_kind=system.row_kind)),
            ("gauss_seidel_rb", lambda: solve_gauss_seidel(
                system.A, system.b, system.colours, tol=tolerance,
                maxiter=500000, row_kind=system.row_kind)),
            ("sor_rb", lambda: solve_red_black_sor(
                system.A, system.b, system.colours, omega=omega,
                tol=tolerance, maxiter=100000, row_kind=system.row_kind)),
        ]
        for method, call in solvers:
            key = (system.grid.N, method, float(tolerance))
            if cache is not None and key in cache:
                result, difference, ratio = cache[key]
            else:
                result = call()
                difference = float(np.max(np.abs(result.x - direct.x)))
                ratio = difference / direct_error["Linf"]
            errors = error_metrics(result.x, exact, system.grid)
            qoi, _ = _right_midpoint_qoi(system, result.x)
            meets_forward = bool(result.converged and ratio <= target)
            row = {
                "N": system.grid.N,
                "method": method,
                "tolerance": tolerance,
                "forward_error_target": target,
                "stagnation_target": target,
                "converged": result.converged,
                "iterations": result.iterations,
                "wall_time": result.wall_time,
                "omega": result.omega,
                **result.residuals,
                "L2h": errors["L2h"],
                "Linf": errors["Linf"],
                "Linf_right": errors["Linf_right"],
                "direct_L2h": direct_error["L2h"],
                "direct_Linf": direct_error["Linf"],
                "direct_qoi_right_midpoint": direct_qoi,
                "direct_qoi_error": direct_qoi_error,
                "qoi_right_midpoint": qoi,
                "qoi_exact": qoi_exact,
                "qoi_error": qoi - qoi_exact,
                "qoi_stability_scale": qoi_scale,
                "difference_from_direct_inf": difference,
                "algebraic_to_discretisation": ratio,
                "meets_forward_error_target": meets_forward,
                "stagnation_check_available": method in previous,
                "abs_L2h_change_from_previous": np.nan,
                "abs_Linf_change_from_previous": np.nan,
                "abs_qoi_change_from_previous": np.nan,
                "normalised_L2h_change": np.nan,
                "normalised_Linf_change": np.nan,
                "normalised_qoi_change": np.nan,
                "max_normalised_error_change": np.nan,
                "meets_stagnation_target": False,
                "meets_reporting_target": False,
            }
            if method in previous:
                old = previous[method]
                change_l2 = abs(errors["L2h"] - old["errors"]["L2h"])
                change_linf = abs(errors["Linf"] - old["errors"]["Linf"])
                change_qoi = abs(qoi - old["qoi"])
                normalised_l2 = change_l2 / direct_error["L2h"]
                normalised_linf = change_linf / direct_error["Linf"]
                normalised_qoi = change_qoi / qoi_scale
                max_error_change = max(normalised_l2, normalised_linf)
                meets_stagnation = bool(
                    result.converged
                    and max_error_change <= target
                    and normalised_qoi <= target
                )
                row.update({
                    "abs_L2h_change_from_previous": change_l2,
                    "abs_Linf_change_from_previous": change_linf,
                    "abs_qoi_change_from_previous": change_qoi,
                    "normalised_L2h_change": normalised_l2,
                    "normalised_Linf_change": normalised_linf,
                    "normalised_qoi_change": normalised_qoi,
                    "max_normalised_error_change": max_error_change,
                    "meets_stagnation_target": meets_stagnation,
                    "meets_reporting_target": bool(meets_forward and meets_stagnation),
                })
            rows.append(row)
            previous[method] = {"errors": errors, "qoi": qoi}
            print(f"E2 accuracy N={system.grid.N}, method={method}, "
                  f"tol={tolerance:g}, ratio={ratio:.3e}, "
                  f"reporting_target={row['meets_reporting_target']}", flush=True)
    return rows


def _accuracy_target_summary(frame, comparison):
    """Select the loosest tolerance passing both forward and stagnation gates."""
    rows = []
    direct_times = comparison[comparison.method == "direct"].set_index("N").wall_time
    for (N, method), group in frame.groupby(["N", "method"]):
        forward = group[group.meets_forward_error_target]
        passing = group[group.meets_reporting_target]
        target = float(group.forward_error_target.iloc[0])
        row = {
            "N": int(N),
            "method": method,
            "forward_error_target": target,
            "stagnation_target": float(group.stagnation_target.iloc[0]),
            "selection_rule": ("loosest tested tolerance with algebraic error and "
                               "successive L2h/Linf/QoI changes each <= target"),
            "least_stringent_forward_tolerance": np.nan,
            "forward_only_iterations": np.nan,
            "forward_only_ratio": np.nan,
            "selected_tolerance": np.nan,
            "selected_iterations": np.nan,
            "selected_wall_time": np.nan,
            "selected_scaled_residual": np.nan,
            "selected_algebraic_to_discretisation": np.nan,
            "selected_L2h": np.nan,
            "selected_Linf": np.nan,
            "selected_qoi_right_midpoint": np.nan,
            "selected_qoi_error": np.nan,
            "selected_max_normalised_error_change": np.nan,
            "selected_normalised_qoi_change": np.nan,
            "direct_wall_time": float(direct_times.loc[N]),
            "selected_to_direct_time_ratio": np.nan,
            "tightest_tested_tolerance": float(group.tolerance.min()),
            "tightest_tested_meets_reporting_target": bool(
                group.loc[group.tolerance.idxmin(), "meets_reporting_target"]
            ),
        }
        if not forward.empty:
            forward_selected = forward.loc[forward.tolerance.idxmax()]
            row.update({
                "least_stringent_forward_tolerance": float(forward_selected.tolerance),
                "forward_only_iterations": int(forward_selected.iterations),
                "forward_only_ratio": float(
                    forward_selected.algebraic_to_discretisation),
            })
        if not passing.empty:
            selected = passing.loc[passing.tolerance.idxmax()]
            row.update({
                "selected_tolerance": float(selected.tolerance),
                "selected_iterations": int(selected.iterations),
                "selected_wall_time": float(selected.wall_time),
                "selected_scaled_residual": float(selected.scaled_relative_inf),
                "selected_algebraic_to_discretisation": float(
                    selected.algebraic_to_discretisation),
                "selected_L2h": float(selected.L2h),
                "selected_Linf": float(selected.Linf),
                "selected_qoi_right_midpoint": float(selected.qoi_right_midpoint),
                "selected_qoi_error": float(selected.qoi_error),
                "selected_max_normalised_error_change": float(
                    selected.max_normalised_error_change),
                "selected_normalised_qoi_change": float(
                    selected.normalised_qoi_change),
                "selected_to_direct_time_ratio": float(
                    selected.wall_time / direct_times.loc[N]),
            })
        rows.append(row)
    summary = pd.DataFrame(rows)
    save_table(summary, "E2_accuracy_target_summary.csv")
    return summary


def _accuracy_target_plot(frame):
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for method, group in frame.groupby("method"):
        group = group.sort_values("tolerance", ascending=False)
        ax.loglog(group.tolerance, group.algebraic_to_discretisation,
                  label=_METHOD_LABELS[method])
    target = float(frame.forward_error_target.iloc[0])
    ax.axhline(target, color="black", linestyle="--", marker="none",
               label=f"{100 * target:g}% forward-error target")
    ax.invert_xaxis()
    N = int(frame.N.iloc[0])
    ax.set(title=f"Algebraic error versus residual tolerance (N={N})",
           xlabel="scaled-residual tolerance (tighter to the right)",
           ylabel=r"$\|x_{iter}-x_{direct}\|_\infty / E_{direct,\infty}$")
    ax.legend(fontsize=9)
    fig.tight_layout()
    save_figure(fig, "E2_accuracy_target_sweep.png")


def _stagnation_plot(frame):
    """Show successive changes in reported errors and the physical QoI."""
    methods = list(frame.method.unique())
    # A vertical layout keeps every curve and its fitted acceptance threshold
    # legible at the physical width of a printed report page.
    fig, axes = plt.subplots(len(methods), 1, figsize=(8.0, 9.6), squeeze=False)
    target = float(frame.stagnation_target.iloc[0])
    for ax, method in zip(axes.ravel(), methods):
        group = frame[frame.method == method].sort_values(
            "tolerance", ascending=False)
        ax.loglog(group.tolerance, group.max_normalised_error_change,
                  label=r"$\max(\Delta L_h^2, \Delta L_\infty)$")
        ax.loglog(group.tolerance, group.normalised_qoi_change,
                  label=r"$\Delta u_h(1,1/2)$")
        ax.axhline(target, color="black", linestyle="--", marker="none",
                   label=f"{100 * target:g}% stability target")
        ax.invert_xaxis()
        ax.set(title=_METHOD_LABELS[method],
               xlabel="residual tolerance (tighter to the right)",
               ylabel="successive change / direct error scale")
        ax.legend(fontsize=9)
    N = int(frame.N.iloc[0])
    fig.suptitle(f"Error and midpoint-QoI stagnation audit (N={N})")
    fig.tight_layout()
    save_figure(fig, "E2_tolerance_stagnation.png")


def _matched_accuracy_cost_plot(summary):
    """Cost comparison at tolerances selected by accuracy and stagnation."""
    if summary.selected_tolerance.isna().any():
        raise RuntimeError("A solver has no tolerance passing both reporting gates")
    N = int(summary.N.iloc[0])
    direct_time = float(summary.direct_wall_time.iloc[0])
    ordered = summary.sort_values("method")
    labels = [_METHOD_LABELS["direct"]] + [
        _METHOD_LABELS[method] for method in ordered.method]
    times = [direct_time] + ordered.selected_wall_time.tolist()
    tolerances = ["SuperLU"] + [f"tol={tol:g}" for tol in ordered.selected_tolerance]
    colours = ["black", "#4C78A8", "#F58518", "#54A24B"]
    fig, ax = plt.subplots(figsize=(8.2, 4.9))
    bars = ax.bar(labels, times, color=colours[:len(labels)])
    ax.set_yscale("log")
    ax.set(ylabel="wall time (seconds, log scale)",
           title=f"Matched accuracy-and-stagnation cost (N={N})")
    ax.tick_params(axis="x", rotation=12)
    for bar, seconds, tolerance in zip(bars, times, tolerances):
        ax.annotate(f"{seconds:.3g} s\n{tolerance}",
                    (bar.get_x() + bar.get_width() / 2, seconds),
                    xytext=(0, 4), textcoords="offset points",
                    ha="center", va="bottom", fontsize=8)
    ax.set_ylim(min(times) * 0.55, max(times) * 4.0)
    fig.tight_layout()
    save_figure(fig, "E2_matched_accuracy_cost.png")


def _spectral_audit(Ns=(8, 12, 16)):
    rows = []
    for N in Ns:
        system = assemble_core(N, "neumann", "ghost")
        omega = 2.0 / (1.0 + np.sin(np.pi / N))
        for method in ("jacobi", "gauss_seidel", "sor"):
            G = _iteration_matrix(system.A, system.colours, method, omega)
            rho = float(np.max(np.abs(np.linalg.eigvals(G))))
            rows.append({"N": N, "method": method, "measured_rho": rho,
                         "omega": omega if method == "sor" else np.nan,
                         "pure_Dirichlet_anchor": (
                             np.cos(np.pi / N) if method == "jacobi" else
                             np.cos(np.pi / N) ** 2 if method == "gauss_seidel" else
                             omega - 1.0)})
    save_table(pd.DataFrame(rows), "E2_spectral_radius_audit.csv")


def _stopping_audit(N=64):
    system = assemble_core(N, "neumann", "ghost")
    rows = []
    for name, omega in (("gauss_seidel_rb", 1.0), ("sor_rb", 1.85)):
        for criterion in ("residual", "iterate"):
            result = solve_red_black_sor(
                system.A, system.b, system.colours, omega=omega,
                tol=1e-9, iterate_tol=1e-8, maxiter=200000,
                stop_criterion=criterion, row_kind=system.row_kind,
            )
            rows.append({"N": N, "method": name, "criterion": criterion,
                         "iterations": result.iterations,
                         "final_iterate_change": result.final_iterate_change,
                         **result.residuals})
    save_table(pd.DataFrame(rows), "E2_stopping_criterion_audit.csv")


def _crossover_summary(frame):
    rows = []
    direct = frame[frame.method == "direct"].set_index("N")
    target_tolerance = frame.tolerance.dropna().min()
    for method in ("jacobi", "gauss_seidel_rb", "sor_rb"):
        subset = frame[(frame.method == method) & (frame.tolerance == target_tolerance)]
        ratios = []
        for _, row in subset.sort_values("N").iterrows():
            if row.N in direct.index:
                ratios.append((int(row.N), float(row.wall_time / direct.loc[row.N].wall_time)))
        rows.append({"method": method, "direct_faster_at_all_matched_N": all(r > 1 for _, r in ratios),
                     "observed_crossover": False if all(r > 1 for _, r in ratios) else True,
                     "time_ratio_by_N": str(ratios)})
    save_table(pd.DataFrame(rows), "E2_crossover_summary.csv")


def _plots(frame, omega_frame):
    target_tolerance = frame.tolerance.dropna().min()
    subset = frame[(frame.tolerance == target_tolerance) & frame.converged]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for method, group in subset.groupby("method"):
        group = group.sort_values("N")
        axes[0].loglog(group.N, group.iterations, label=_METHOD_LABELS[method])
    direct = frame[frame.method == "direct"]
    for method, group in pd.concat((direct, subset)).groupby("method"):
        group = group.sort_values("N")
        axes[1].loglog(group.N, group.wall_time, label=_METHOD_LABELS[method])
    axes[0].set(title="Residual-stopped iteration count", xlabel="N", ylabel="iterations")
    axes[1].set(title=f"Measured wall time (tol={target_tolerance:g})", xlabel="N", ylabel="seconds")
    mesh_levels = sorted(frame.N.unique())
    for ax in axes:
        ax.set_xticks(mesh_levels)
        ax.set_xticklabels([str(int(N)) for N in mesh_levels])
        ax.legend(fontsize=9)
    fig.tight_layout()
    save_figure(fig, "E2_cost_comparison.png")

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for N, group in omega_frame.groupby("N"):
        ax.semilogy(group.omega, group.iterations, label=f"N={N}")
    ax.set(xlabel="SOR relaxation omega", ylabel="iterations to scaled residual 1e-8",
           title="Measured omega sweep (not an assumed optimum)")
    ax.legend(fontsize=8)
    save_figure(fig, "E2_omega_scan.png")
