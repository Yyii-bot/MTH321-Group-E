"""make_figures_supplementary.py — regenerate the two report figures that are
NOT produced by make_figures.py:

    fig_refcheck.pdf/.png  (report Figure 2) — Radau reference verification,
        recomputed deterministically from model.py (no stored data needed).
    fig_matchedcost.png    (report Figure 10) — matched-error EE/IE runtime
        benchmark.

Figure 10 plots the wall times measured on the report's benchmark machine
(Python 3.12.14, NumPy 2.2.6, SciPy 1.15.3; macOS 26.6.2, arm64). Because
timings are machine-dependent, the original values are stored below (they were
read back from the published figure; the four IE/EE ratios are 5.16-5.25 as
quoted in the report). Run

    python make_figures_supplementary.py --rerun

to re-measure the benchmark on the current machine (protocol: Appendix C.2 of
the report: one warm-up run, then five timed runs, minimum reported) and
replot Figure 10 with the fresh numbers.

Usage:  python make_figures_supplementary.py           # rebuild both figures
"""
import os
import sys
import time

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(FIG, exist_ok=True)

from model import rhs, jac, Y0, T_FINAL, output_grid
from methods import euler_explicit, euler_implicit

INK = "#253746"          # dark slate used for titles/labels
GRID = dict(alpha=0.28, linewidth=0.5)

mpl.rcParams.update({
    "font.family": "sans-serif", "font.size": 10,
    "axes.edgecolor": "0.25", "axes.linewidth": 0.9,
    "axes.grid": True, "grid.alpha": 0.28, "grid.linewidth": 0.5,
    "xtick.direction": "out", "ytick.direction": "out",
    "legend.framealpha": 0.95, "legend.edgecolor": "0.8",
})


def savefig(fig, name, dpi=200):
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- Figure 2
def fig_refcheck():
    """100-fold tolerance tightening of the Radau reference oracle.

    R0: Radau at (rtol, atol) = (1e-10, 1e-12);  R1: (1e-12, 1e-14).
    Both are sampled on the protocol output grid (t = 0 plus 241
    logarithmically spaced times in [1e-8, 40]); the absolute change
    |R0 - R1| is plotted. Deterministic: needs no stored data files.
    """
    t_out = output_grid()
    r0 = solve_ivp(rhs, (0, T_FINAL), Y0, method="Radau", jac=jac,
                   rtol=1e-10, atol=1e-12, t_eval=t_out)
    r1 = solve_ivp(rhs, (0, T_FINAL), Y0, method="Radau", jac=jac,
                   rtol=1e-12, atol=1e-14, t_eval=t_out)
    assert r0.success and r1.success
    d = np.abs(r0.y - r1.y)                      # shape (3, 242)
    endpoint = d[:, -1]                          # at t = 40
    gridmax = d[:, 1:].max(axis=1)               # over the output grid

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.8, 4.9),
                                   gridspec_kw=dict(width_ratios=[1.35, 1.0],
                                                    wspace=0.22))
    fig.suptitle("Robertson reference verification: 100-fold tolerance "
                 "tightening", fontsize=17, fontweight="bold", color=INK,
                 x=0.5, y=0.985)
    fig.text(0.5, 0.915,
             r"Radau: (rtol, atol) = $(10^{-10},\, 10^{-12})$  $\longrightarrow$"
             r"  $(10^{-12},\, 10^{-14})$",
             ha="center", fontsize=13, color=INK)

    # -- panel A: changes on the common output grid -------------------------
    tp = t_out[1:]
    style = [dict(color="#4977A4", ls="-",  lw=1.6),
             dict(color="#D8A250", ls="--", lw=1.6, marker="o", ms=4,
                  markevery=18),
             dict(color="#1F8779", ls=":",  lw=2.0, marker="h", ms=4.5,
                  markevery=18)]
    labels = [r"$y_1$", r"$y_2$", r"$y_3$"]
    for i in range(3):
        m = d[i, 1:] > 0           # exact-zero changes omitted from log axes
        axA.loglog(tp[m], d[i, 1:][m], label=labels[i], **style[i])
    axA.set_xlim(1e-8, 40)
    axA.set_ylim(3e-20, 1e-11)
    axA.set_xlabel(r"Time $t$", fontsize=12)
    axA.set_ylabel(r"Absolute change $|y_i^{R_0} - y_i^{R_1}|$", fontsize=12)
    axA.set_title("A   Changes on the common output grid", loc="left",
                  fontsize=13, fontweight="bold", color=INK)
    axA.legend(loc="upper left", ncol=3, fontsize=11, handlelength=2.6,
               columnspacing=1.4)

    # -- panel B: endpoint and trajectory checks ----------------------------
    x = np.arange(3)
    w = 0.38
    b1 = axB.bar(x - w / 2, endpoint, w, color="#4977A4", label=r"At $t = 40$")
    b2 = axB.bar(x + w / 2, gridmax, w, color="#D8A250",
                 label="Maximum on output grid")
    for bars in (b1, b2):
        for b in bars:
            axB.annotate(f"{b.get_height():.2e}",
                         (b.get_x() + b.get_width() / 2, b.get_height()),
                         xytext=(0, 3), textcoords="offset points",
                         ha="center", fontsize=9, color="0.25")
    axB.set_yscale("log")
    axB.set_ylim(3e-17, 1e-10)
    axB.set_xticks(x)
    axB.set_xticklabels([r"$y_1$", r"$y_2$", r"$y_3$"], fontsize=12)
    axB.set_ylabel("Absolute change", fontsize=12)
    axB.set_title("B   Endpoint and trajectory checks", loc="left",
                  fontsize=13, fontweight="bold", color=INK)
    axB.legend(loc="upper right", fontsize=10)
    axB.grid(axis="x", visible=False)

    fig.text(0.5, 0.052,
             "Output grid: $t = 0$ plus 241 logarithmically spaced times "
             "from 1e-8 to 40.",
             ha="center", fontsize=11, fontweight="bold", color=INK)
    fig.text(0.5, 0.016,
             "Exact-zero changes and $t = 0$ are omitted from log axes. "
             "Changes are consistency checks, not true-error bounds.",
             ha="center", fontsize=10.5, color="0.35")
    fig.subplots_adjust(left=0.065, right=0.99, top=0.83, bottom=0.185)
    savefig(fig, "fig_refcheck")


# ---------------------------------------------------------------- Figure 10
# Matched-error pairs (P1..P4): target = measured EE weighted endpoint error
# at N uniform steps; IE run at the same N is matched (error within 5 %).
MATCHED_N = [80_000, 160_000, 320_000, 640_000]
MATCHED_E = [5.938, 2.969, 1.485, 0.7423]      # weighted endpoint error E(40)
# Wall times [s] on the report's benchmark machine (minimum of 5 runs after
# one warm-up run; environment: Python 3.12.14, NumPy 2.2.6, SciPy 1.15.3,
# macOS 26.6.2 arm64). Values read back from the published figure; the four
# IE/EE ratios are 5.16-5.25 as quoted in the report.
MATCHED_EE_S = [0.121, 0.244, 0.490, 0.980]
MATCHED_IE_S = [0.627, 1.269, 2.530, 5.070]


def measure_matched(newton_tol=1e-13, repeats=5):
    """Re-run the matched-error benchmark of Appendix C.2 on this machine.

    EE and IE are run at the same four uniform step counts; each entry is the
    minimum wall time over `repeats` timed runs after one warm-up run.
    """
    def timed(fn, *args):
        fn(*args)                                # warm-up
        best = np.inf
        for _ in range(repeats):
            t0 = time.perf_counter()
            fn(*args)
            best = min(best, time.perf_counter() - t0)
        return best

    ee_s, ie_s, errs = [], [], []
    ref_end = solve_ivp(rhs, (0, T_FINAL), Y0, method="Radau", jac=jac,
                        rtol=1e-12, atol=1e-14).y[:, -1]
    w = 1e-8 + 1e-6 * np.abs(ref_end)
    for n in MATCHED_N:
        h = T_FINAL / n
        ee_s.append(timed(euler_explicit, rhs, (0, T_FINAL), Y0, h))
        ie_s.append(timed(euler_implicit, rhs, jac, (0, T_FINAL), Y0, h,
                          newton_tol))
        _, ys, _ = euler_implicit(rhs, jac, (0, T_FINAL), Y0, h, newton_tol)
        errs.append(float(np.max(np.abs(ys[-1] - ref_end) / w)))
        print(f"N={n:>7d}  EE {ee_s[-1]:.3f} s   IE {ie_s[-1]:.3f} s   "
              f"E_IE {errs[-1]:.4g}")
    return ee_s, ie_s, errs


def fig_matchedcost(ee_s=None, ie_s=None):
    ee_s = MATCHED_EE_S if ee_s is None else ee_s
    ie_s = MATCHED_IE_S if ie_s is None else ie_s
    E = MATCHED_E

    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    for k in range(4):                           # dotted pair connectors
        ax.plot([ee_s[k], ie_s[k]], [E[k], E[k]], ls=":", color="0.63",
                lw=1.4, zorder=1)
        xm = np.exp(0.5 * (np.log(ee_s[k]) + np.log(ie_s[k])))
        ax.annotate(f"P{k + 1}", (xm, E[k]), xytext=(-4, 7),
                    textcoords="offset points", fontsize=9, color="0.45")
    ax.loglog(ee_s, E, "-o", color="#B52438", lw=1.6, ms=6, label="EE",
              zorder=3)
    ax.loglog(ie_s, E, "-s", color="#006E66", lw=1.6, ms=6, label="IE",
              zorder=3)
    ax.set_xlim(0.06, 7)
    ax.set_ylim(0.45, 8)
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5])
    ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5"])
    ax.set_yticks([0.5, 1, 2, 5])
    ax.set_yticklabels(["0.5", "1", "2", "5"])
    ax.minorticks_off()
    ax.set_xlabel("Wall time [s] (minimum of 5 runs)", fontsize=10.5)
    ax.set_ylabel(r"Weighted endpoint error $E(40)$", fontsize=10.5)
    ax.set_title("A   Cost at matched endpoint accuracy", loc="left",
                 fontsize=11.5, fontweight="bold", color=INK)
    ax.legend(loc="upper right", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_matchedcost.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("wrote fig_matchedcost")


if __name__ == "__main__":
    fig_refcheck()
    if "--rerun" in sys.argv:
        ee_s, ie_s, _ = measure_matched()
        fig_matchedcost(ee_s, ie_s)
    else:
        fig_matchedcost()
