"""E5: uniform, red-graded, and radial-graded P1 FEM on the L-shaped domain."""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..data import singular_exact, singular_gradient
from ..fem import (lshape_triangulation, corner_graded_triangulation,
                   radial_graded_triangulation, fem_dirichlet_solve,
                   fem_error_integrals, mesh_size)
from ..metrics import asymptotic_fit
from .common import save_table, save_figure


def _mesh_quality(mesh):
    v = mesh.p[:, mesh.t].transpose(2, 1, 0)
    a = np.linalg.norm(v[:, 1] - v[:, 0], axis=1)
    b = np.linalg.norm(v[:, 2] - v[:, 1], axis=1)
    c = np.linalg.norm(v[:, 0] - v[:, 2], axis=1)
    area = 0.5 * np.abs((v[:, 1, 0] - v[:, 0, 0]) * (v[:, 2, 1] - v[:, 0, 1])
                        - (v[:, 2, 0] - v[:, 0, 0]) * (v[:, 1, 1] - v[:, 0, 1]))
    angles = []
    for opposite, x, y in ((a, b, c), (b, c, a), (c, a, b)):
        cosine = np.clip((x**2 + y**2 - opposite**2) / (2.0 * x * y), -1.0, 1.0)
        angles.append(np.arccos(cosine))
    min_angle = np.degrees(np.min(np.stack(angles)))
    altitude_min = 2.0 * area / np.maximum.reduce((a, b, c))
    incidence_count = np.sum(mesh.f2t >= 0, axis=0)
    interior_facet_ok = bool(
        np.all((incidence_count == 1) | (incidence_count == 2))
        and np.sum(incidence_count == 1) == len(mesh.boundary_facets())
    )
    return {
        "h_max": mesh_size(mesh),
        "min_angle_deg": float(min_angle),
        "max_edge_over_min_altitude": float(np.max(np.maximum.reduce((a, b, c)) / altitude_min)),
        "min_triangle_area": float(area.min()),
        "facet_incidence_valid": interior_facet_ok,
    }


def run_e5(Ns=(16, 32, 64, 128), kappa=6):
    rows = []
    cache = {}

    def compute(mesh, family, N, parameter):
        key = (family, N, parameter)
        if key in cache:
            return cache[key]
        zero = lambda x, y: np.zeros_like(np.asarray(x, dtype=float))
        result = fem_dirichlet_solve(mesh, zero, singular_exact, intorder=6)
        errors = fem_error_integrals(mesh, result.values, singular_exact,
                                     singular_gradient, singular_corner=True)
        row = {
            "family": family, "N": N, "parameter": parameter,
            "n_dof": result.A.shape[0], "nnz": result.A.nnz,
            "n_vertices": mesh.p.shape[1], "n_elements": mesh.t.shape[1],
            "assembly_time": result.assembly_time, "solve_time": result.solve_time,
            **result.residuals, **errors, **_mesh_quality(mesh),
        }
        rows.append(row)
        cache[key] = row
        print(f"E5 completed {family}, N={N}, parameter={parameter}", flush=True)
        return row

    for N in Ns:
        compute(lshape_triangulation(N), "uniform", N, 0.0)
        compute(corner_graded_triangulation(N, kappa=kappa), "red_graded", N, float(kappa))
    for local_kappa in (2, 3, 4, 5, 6):
        compute(corner_graded_triangulation(32, kappa=local_kappa),
                "kappa_scan", 32, float(local_kappa))
    for n in (16, 32):
        for mu in (1.0, 1.25, 1.5, 2.0, 3.0):
            compute(radial_graded_triangulation(n, mu), "radial_mu_scan", n, mu)

    frame = pd.DataFrame(rows)
    save_table(frame, "E5_graded_mesh_convergence.csv")
    fit_rows = []
    for family in ("uniform", "red_graded"):
        group = frame[frame.family == family].sort_values("n_dof")
        for metric in ("H1", "L2", "Linf"):
            fit = asymptotic_fit(group.n_dof.to_numpy(), group[metric].to_numpy(), last=3)
            fit_rows.append({"family": family, "metric": metric, **fit})
    save_table(pd.DataFrame(fit_rows), "E5_error_vs_dof_fits.csv")
    _plots(frame, pd.DataFrame(fit_rows))
    return frame


def _plots(frame, fits):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, metric in zip(axes, ("H1", "L2")):
        for family, group in frame[frame.family.isin(("uniform", "red_graded"))].groupby("family"):
            group = group.sort_values("n_dof")
            fit = fits.query("family == @family and metric == @metric").iloc[0]
            ax.loglog(group.n_dof, group[metric],
                      label=(f"{family}: slope {-fit.order:.3f}, "
                             f"$\\delta_{{\\log}}={fit.max_log_deviation:.2e}$"))
        uniform = frame[frame.family == "uniform"].sort_values("n_dof")
        reference_order = 1 / 3 if metric == "H1" else 2 / 3
        ax.loglog(uniform.n_dof,
                  uniform[metric].iloc[-1] * (uniform.n_dof / uniform.n_dof.iloc[-1]) ** (-reference_order),
                  marker="none", linestyle="--", label=f"uniform theory $N_dof^{{-{reference_order:.3f}}}$")
        graded = frame[frame.family == "red_graded"].sort_values("n_dof")
        target = 0.5 if metric == "H1" else 1.0
        ax.loglog(graded.n_dof,
                  graded[metric].iloc[-1] * (graded.n_dof / graded.n_dof.iloc[-1]) ** (-target),
                  marker="none", linestyle=":", label=f"graded target $N_dof^{{-{target}}}$")
        ax.set(xlabel="interior degrees of freedom", ylabel=f"integral {metric} error",
               title=f"{metric} error versus DOFs")
        ax.legend(fontsize=7)
    fig.tight_layout()
    save_figure(fig, "E5_error_vs_dofs.png")

    scan = frame[(frame.family == "kappa_scan") |
                 ((frame.family == "uniform") & (frame.N == 32))]
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for family, group in scan.groupby("family"):
        group = group.sort_values("n_dof")
        ax.loglog(group.n_dof, group.H1, label=family)
        for row in group.itertuples():
            label = "uniform" if family == "uniform" else rf"$\kappa={row.parameter:g}$"
            ax.annotate(label, (row.n_dof, row.H1), xytext=(5, 5),
                        textcoords="offset points", fontsize=8)
    ax.set(xlabel="interior DOFs", ylabel="integral H1 error",
           title="Local red-refinement scan on the same N=32 outer mesh")
    ax.legend()
    save_figure(fig, "E5_kappa_scan.png")

    radial = frame[frame.family == "radial_mu_scan"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for n, group in radial.groupby("N"):
        group = group.sort_values("parameter")
        axes[0].semilogy(group.parameter, group.H1, label=f"n={n}")
        axes[1].semilogy(group.parameter, group.L2, label=f"n={n}")
    axes[0].axvline(1.5, color="k", linestyle="--", label="theoretical threshold 3/2")
    axes[1].axvline(1.5, color="k", linestyle="--", label="theoretical threshold 3/2")
    for ax, metric in zip(axes, ("H1", "L2")):
        ax.set(xlabel="radial grading exponent mu", ylabel=f"integral {metric} error",
               title=f"Radial mu scan: {metric}")
        ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, "E5_mu_scan.png")

    _mesh_plots()


def _mesh_plots():
    meshes = [
        (lshape_triangulation(16), "Uniform N=16", (-1.05, 1.05, -1.05, 1.05)),
        (corner_graded_triangulation(16, 6), "Red-graded, kappa=6", (-1.05, 1.05, -1.05, 1.05)),
        (corner_graded_triangulation(16, 6), "Red-graded: corner zoom", (-0.35, 0.35, -0.35, 0.35)),
        (radial_graded_triangulation(12, 2.0), "Radial mu=2: corner zoom", (-0.35, 0.35, -0.35, 0.35)),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    for ax, (mesh, title, limits) in zip(axes.ravel(), meshes):
        ax.triplot(mesh.p[0], mesh.p[1], mesh.t.T, linewidth=0.35, color="navy")
        ax.set_aspect("equal")
        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.tight_layout()
    save_figure(fig, "E5_mesh_families.png")
