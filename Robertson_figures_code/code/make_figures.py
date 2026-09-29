"""make_figures.py — regenerate all report figures from results/."""
import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(FIG, exist_ok=True)

mpl.rcParams.update({
    "font.family": "serif", "font.size": 9.5, "axes.titlesize": 10,
    "axes.labelsize": 10, "axes.linewidth": 0.8, "axes.grid": True,
    "grid.alpha": 0.28, "grid.linewidth": 0.5, "legend.framealpha": 0.92,
    "legend.edgecolor": "0.7", "legend.fontsize": 8.5, "lines.linewidth": 1.5,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
    "ytick.right": True, "savefig.bbox": "tight",
    "axes.prop_cycle": mpl.cycler(color=["#8C2D2D", "#B07A3C", "#3E6B5B",
                                         "#4A4E69", "#7A7A52", "#5B3A4E"]),
})
C = dict(ee="#8C2D2D", rk4="#B07A3C", ie="#3E6B5B", radau="#4A4E69",
         bdf="#7A7A52", ref="0.35")


def savefig(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf"))
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300)
    plt.close(fig)
    print("wrote", name)


S = json.load(open(os.path.join(RES, "summary.json")))
ref = np.load(os.path.join(RES, "reference.npz"))
eig = np.load(os.path.join(RES, "eig.npz"))
t, Y = ref["t"], ref["y"]
LMAX40 = 3392.788124450614


# ---------------------------------------------------------------- fig 1
def fig_solution():
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    tp = t[1:]
    ax.semilogx(tp, Y[0, 1:], color=C["ie"], label=r"$y_1$ (A)")
    ax.semilogx(tp, Y[1, 1:], color=C["ee"], label=r"$y_2$ (B)")
    ax.semilogx(tp, Y[2, 1:], color=C["rk4"], label=r"$y_3$ (C)")
    ax.set_yscale("symlog", linthresh=1e-6)
    ax.set_xlabel(r"$t$"); ax.set_ylabel("concentration"); ax.set_xlim(1e-8, 40)
    ax.legend(loc="center left", bbox_to_anchor=(0.02, 0.62))
    ax2 = ax.inset_axes([0.44, 0.16, 0.30, 0.36])
    ax2.semilogx(tp, Y[1, 1:], color=C["ee"]); ax2.set_xlim(1e-8, 1e-1)
    ax2.grid(True, alpha=0.28); ax2.tick_params(labelsize=7.5)
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax2.set_title(r"$y_2$ early transient (peak $\approx 3.7\times10^{-5}$)",
                  fontsize=8)
    savefig(fig, "fig_solution")


# ---------------------------------------------------------------- fig 2
def fig_stiffness():
    lmax, lmin, Sr = eig["lmax"], eig["lmin"], eig["S"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 2.9))
    m = np.isfinite(lmax) & (t > 0)
    a1.loglog(t[m], lmax[m], color=C["ee"], label=r"$|\lambda_{\mathrm{fast}}(t)|$")
    a1.loglog(t[m], lmin[m], color=C["ie"], label=r"$|\lambda_{\mathrm{slow}}(t)|$")
    a1.set_xlabel(r"$t$"); a1.set_ylabel(r"$|\lambda|$")
    a1.legend(loc="lower left", fontsize=8)
    a1.set_title("Nonzero Jacobian eigenvalues", fontsize=9.5)
    a2.loglog(t[m], Sr[m], color=C["rk4"])
    a2.set_xlabel(r"$t$"); a2.set_ylabel(r"$S(t)$")
    a2.set_title("Stiffness ratio (two nonzero eigenvalues)", fontsize=9.5)
    for tv, sv in [(1e-4, 3.0e3), (1e-2, 5.4e3), (40, 1.584e5)]:
        a2.plot(tv, sv, "o", ms=4, color=C["ee"])
    a2.annotate(r"$1.58\times10^{5}$ at $t=40$", xy=(40, 1.584e5),
                xytext=(3e-3, 4e4), fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.7))
    fig.tight_layout(); savefig(fig, "fig_stiffness")


# ---------------------------------------------------------------- fig 3
def fig_stability():
    lmax = eig["lmax"]

    def R4(z):
        return 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24
    xs = np.linspace(-4.2, 2.2, 700); ys_ = np.linspace(-3.2, 3.2, 560)
    X, Yg = np.meshgrid(xs, ys_); Z = X + 1j * Yg
    fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.65), sharey=True)
    regions = [(np.abs(1 + Z), "Explicit Euler"),
               (np.abs(R4(Z)), "Classical RK4"),
               (np.abs(1 / (1 - Z)), "Implicit Euler")]
    for ax, (M, title) in zip(axes, regions):
        ax.contourf(X, Yg, M, levels=[0, 1], colors=["#DCE5DF"], alpha=0.9)
        ax.contour(X, Yg, M, levels=[1], colors=[C["ie"]], linewidths=1.4)
        ax.axhline(0, color="0.4", lw=0.6); ax.axvline(0, color="0.4", lw=0.6)
        ax.set_xlabel(r"$\mathrm{Re}(h\lambda)$"); ax.set_title(title, fontsize=9.5)
        ax.set_aspect("equal")
    axes[0].set_ylabel(r"$\mathrm{Im}(h\lambda)$")
    tmask = t >= 1e-6
    for ax, hh, mk, col, lab in [
            (axes[0], 5.0e-4, "o", C["rk4"], r"$h=5.0\times10^{-4}$"),
            (axes[0], 6.5e-4, "s", C["ee"], r"$h=6.5\times10^{-4}$"),
            (axes[1], 8.0e-4, "o", C["rk4"], r"$h=8.0\times10^{-4}$"),
            (axes[1], 9.5e-4, "s", C["ee"], r"$h=9.5\times10^{-4}$")]:
        zt = hh * lmax[tmask]
        sel = np.linspace(0, len(zt) - 1, 24).astype(int)
        ax.plot(-zt[sel], np.zeros_like(zt[sel]), mk, ms=3.2, color=col, label=lab)
    for ax in axes[:2]:
        ax.legend(fontsize=7.2, loc="lower left")
    axes[2].plot(-1e-4 * lmax[tmask], np.zeros(tmask.sum()), "o", ms=3,
                 color=C["rk4"], label=r"$h=10^{-4}$")
    axes[2].legend(fontsize=7.2, loc="lower right")
    axes[2].text(1.55, 2.5, "stable\noutside", fontsize=7.5, ha="center",
                 color=C["ie"])
    fig.tight_layout(); savefig(fig, "fig_stability")


# ---------------------------------------------------------------- fig 4
def fig_threshold():
    ee = S["threshold"]["ee"]; rk = S["threshold"]["rk4"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.8, 2.9))
    for rows, col, lab, bound in ((ee, C["ee"], "Explicit Euler", 2 / LMAX40),
                                  (rk, C["rk4"], "Classical RK4", 2.785293 / LMAX40)):
        h = [r["h"] for r in rows if r["err40"] is not None]
        e = [r["err40"] for r in rows if r["err40"] is not None]
        hb = [r["h"] for r in rows if r["err40"] is None]
        a1.loglog(h, e, "o-", color=col, ms=4, label=lab)
        for x in hb:
            a1.plot(x, 3e7, "^", color=col, ms=6)
        a1.axvline(bound, color=col, ls="--", lw=1.0, alpha=0.75)
    a1.set_xlabel(r"$h$"); a1.set_ylabel(r"weighted error at $t=40$")
    a1.set_title("Accuracy breakdown at the stability boundary", fontsize=9.5)
    a1.legend(fontsize=8, loc="lower left")
    a1.set_ylim(1e-8, 1e9); a1.set_xlim(1.5e-4, 1.3e-3)
    a1.set_xticks([2e-4, 4e-4, 6e-4, 8e-4, 1e-3])
    a1.set_xticklabels([r"$2$", r"$4$", r"$6$", r"$8$", r"$10$"])
    a1.text(2.1e-4, 1e-6, r"$h\;[\times10^{-4}]$", fontsize=8, color="0.35")
    a1.text(5.6e-4, 3e-5, "EE bound", rotation=90, fontsize=7, color=C["ee"])
    a1.text(8.6e-4, 3e-5, "RK4 bound", rotation=90, fontsize=7, color=C["rk4"])
    ts_ = np.load(os.path.join(RES, "traj_ee_stable.npz"))
    tm_ = np.load(os.path.join(RES, "traj_ee_marginal.npz"))
    m1 = ts_["t"] >= 30; m2 = tm_["t"] >= 30
    a2.semilogy(ts_["t"][m1], np.abs(ts_["y"][m1, 1]), color=C["ie"], lw=1.4,
                label=r"$h=5.9\times10^{-4}$ (stable)")
    a2.semilogy(tm_["t"][m2], np.abs(tm_["y"][m2, 1]), color=C["ee"], lw=1.4,
                label=r"$h=6.1\times10^{-4}$ (marginal)")
    a2.set_xlabel(r"$t$"); a2.set_ylabel(r"$|y_2(t)|$")
    a2.set_title("Explicit Euler: late-time view", fontsize=9.5)
    a2.legend(fontsize=8, loc="lower left")
    fig.tight_layout(); savefig(fig, "fig_threshold")


# ---------------------------------------------------------------- fig 5
def fig_convergence():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.8, 2.9))
    pr = S["pr_nonstiff"]
    h = [r["h"] for r in pr]
    for key, col, lab in (("ee", C["ee"], "Explicit Euler"),
                          ("ie", C["ie"], "Implicit Euler"),
                          ("rk4", C["rk4"], "Classical RK4")):
        a1.loglog(h, [r[key] for r in pr], "o-", color=col, ms=4, label=lab)
    a1.loglog(h, [0.08 * x for x in h], ":", color="0.5", lw=1)
    a1.loglog(h, [1.2e-4 * x**4 for x in h], "--", color="0.5", lw=1)
    a1.text(0.03, 5e-3, r"$\propto h$", fontsize=8, color="0.4")
    a1.text(0.035, 3e-9, r"$\propto h^4$", fontsize=8, color="0.4")
    a1.set_xlabel(r"$h$"); a1.set_ylabel(r"error at $t=1$")
    a1.set_title("Order verification (Prothero–Robinson, $\\lambda=-1$)",
                 fontsize=9.5)
    a1.legend(fontsize=8, loc="upper left")
    conv = S["convergence"]
    for name, col, mk in (("EE", C["ee"], "o"), ("IE", C["ie"], "s"),
                          ("RK4", C["rk4"], "^")):
        hs = [r["h"] for r in conv if r["method"] == name]
        es = [r["err40"] for r in conv if r["method"] == name]
        a2.loglog(hs, es, mk + "-", color=col, ms=4, label=name)
    a2.axvline(2 / LMAX40, color=C["ee"], ls="--", lw=1, alpha=0.7)
    a2.axvline(2.785293 / LMAX40, color=C["rk4"], ls="--", lw=1, alpha=0.7)
    hg = [0.008, 0.08]
    a2.loglog(hg, [93 * (x / 0.008) for x in hg], ":", color="0.5", lw=1)
    a2.text(0.011, 300, r"$\propto h$", fontsize=8, color="0.4")
    a2.axhline(1e-6, color="0.6", ls="-.", lw=0.8)
    a2.text(2.2e-5, 1.8e-6, "reference-solution noise floor", fontsize=7, color="0.4")
    a2.set_ylim(1e-7, 3e3)
    a2.set_xlabel(r"$h$"); a2.set_ylabel(r"weighted error at $t=40$")
    a2.set_title("Stiff problem: convergence and stability cut-offs", fontsize=9.5)
    a2.legend(fontsize=8, loc="upper right")
    fig.tight_layout(); savefig(fig, "fig_convergence")


# ---------------------------------------------------------------- fig 6
def fig_adaptive():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.8, 2.9))
    for rt, col in (("1e-03", C["ee"]), ("1e-04", C["rk4"]), ("1e-05", C["ie"])):
        d = np.load(os.path.join(RES, f"adaptive_{rt}.npz"))
        a1.loglog(d["t"][1:], d["h"], color=col, lw=1.2,
                  label=rf"rtol $=10^{{-{rt[-1]}}}$")
    a1.set_xlabel(r"$t$"); a1.set_ylabel(r"accepted step size $h$")
    a1.set_title("Adaptive IE: step-size history", fontsize=9.5)
    a1.legend(fontsize=8)
    d = np.load(os.path.join(RES, "adaptive_1e-04.npz"))
    a2.loglog(d["t"][1:], d["err"], color=C["ie"], lw=0.9)
    a2.axhline(1.0, color=C["ee"], ls="--", lw=1)
    a2.text(1.2e-6, 1.6, "acceptance threshold", fontsize=7.5, color=C["ee"])
    a2.set_xlabel(r"$t$"); a2.set_ylabel("estimated local error (weighted)")
    a2.set_title(r"Error estimate vs. tolerance (rtol $=10^{-4}$)", fontsize=9.5)
    fig.tight_layout(); savefig(fig, "fig_adaptive")


# ---------------------------------------------------------------- fig 7
def fig_diagnostics():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.8, 2.9))
    for name, col, lab in (("EE", C["ee"], "EE, $h=10^{-4}$"),
                           ("RK4", C["rk4"], "RK4, $h=10^{-4}$"),
                           ("IE", C["ie"], "IE, $h=10^{-4}$"),
                           ("IE_(h=0.05)", C["radau"], "IE, $h=0.05$")):
        d = np.load(os.path.join(RES, f"cons_{name}.npz"))
        a1.semilogy(d["t"], d["defect"] + 1e-18, lw=1.0, color=col, label=lab)
    a1.set_xlabel(r"$t$"); a1.set_ylabel(r"$|y_1+y_2+y_3-1|$")
    a1.set_title("Invariant defect stays at round-off", fontsize=9.5)
    a1.legend(fontsize=7.5)
    nt = S["newton_tol"]
    xs = [r["newton_tol"] for r in nt]
    res_ok = [r["err40"] if r["err40"] == r["err40"] else np.nan for r in nt]
    iters = [r["mean_newton_iter"] for r in nt]
    a2.loglog(xs, res_ok, "o-", color=C["ie"], ms=5,
              label="weighted error at $t=40$")
    a2.set_xscale("log"); a2.set_xlabel("Newton stopping tolerance")
    a2.set_ylabel("weighted error at $t=40$", color=C["ie"])
    a2b = a2.twinx()
    a2b.semilogx(xs, iters, "s--", color=C["rk4"], ms=4,
                 label="mean Newton iterations")
    a2b.set_ylabel("mean Newton iterations / step", color=C["rk4"])
    a2b.grid(False)
    a2.set_title("Inexact Newton: accuracy vs. cost", fontsize=9.5)
    a2.annotate("diverged (NaN)", xy=(1e-2, 5e2), fontsize=8, color=C["ee"],
                ha="center")
    a2.plot(1e-2, 5e2, "x", color=C["ee"], ms=7, mew=2)
    fig.tight_layout(); savefig(fig, "fig_diagnostics")


# ---------------------------------------------------------------- fig 8
def fig_workprecision():
    pts = S["work_precision"]
    fig, ax = plt.subplots(figsize=(6.4, 3.5))
    style = {"EE": (C["ee"], "o"), "RK4": (C["rk4"], "^"),
             "IE uniform": (C["ie"], "s"), "IE adaptive": (C["ie"], "d"),
             "Radau": (C["radau"], "v"), "BDF": (C["bdf"], ">")}
    for name, (col, mk) in style.items():
        xs = [p["cost"] for p in pts if p["method"] == name]
        ys = [p["err40"] for p in pts if p["method"] == name]
        ls = "-" if name in ("EE", "RK4", "IE uniform") else "--"
        ax.loglog(xs, ys, mk + ls, color=col, ms=5, lw=1.1,
                  mfc="white" if name == "IE adaptive" else col, label=name)
    ax.set_xlabel("wall-clock time [s]  (best of repeated runs)")
    ax.set_ylabel(r"weighted error at $t=40$")
    ax.legend(fontsize=8, ncol=2, loc="upper right")
    savefig(fig, "fig_workprecision")


if __name__ == "__main__":
    fig_solution(); fig_stiffness(); fig_stability(); fig_threshold()
    fig_convergence(); fig_adaptive(); fig_diagnostics(); fig_workprecision()
