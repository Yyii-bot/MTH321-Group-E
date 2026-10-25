"""
final_check.py - Week 4 deliverable: pre-submission final check.

1. Key-number y2(40) tolerance scan: fixed-step implicit Euler with successively
   halved h, and adaptive implicit Euler with successively tightened tol; every
   value is compared against the team-built oracle and the true error is quoted
   - the report quotes only digits stable across this table (problem pack:
   the implicit method's y2(40) at several tolerances, cross-checked against
   our own high-accuracy reference).
2. Cost vs accuracy (brief learning outcome D / guide Rule 7):
   the cost measure is named first (right-hand-side evaluations nfev, Newton
   iterations, wall time; Jacobian assembly/factorisation, plotting and
   interpreter startup are NOT counted and are common to both variants), then
   fixed-step and step-doubling adaptive implicit Euler are compared at a
   matched accuracy target err@40 ~ 7e-5.
   Expected to be an honest negative result: step doubling solves three
   implicit steps per accepted step, so at matched error fixed stepping is
   cheaper here - the guide says a negative result is still a result: plot it
   and say so.
3. Report-number traceability table: maps every quoted number to a figure or
   a printed output line.

Run: python final_check.py
Author: Rong Yin (Role 5)
"""

from pathlib import Path
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from robertson_model import robertson_rhs, robertson_jacobian
from solvers import (solve_implicit_euler, adaptive_implicit_euler,
                     WorkCounter)
from reference import build_reference

_HERE = Path(__file__).resolve().parent
if _HERE.name == "validation":
    # team-repo layout (code/validation/*.py): figures -> <repo>/figures/validation
    FIGDIR = _HERE.parents[1] / "figures" / "validation"
else:
    # standalone layout (code/*.py): figures -> sibling figures/
    FIGDIR = _HERE.parent / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

# Problem-pack t=40 orientation value (final sanity comparison only; the
# report quotes the digits confirmed by our own runs)
PACK_Y40 = np.array([0.7158271, 9.1855e-6, 0.2841637])


def y2_tolerance_scan(ref):
    """Implicit-Euler y2(40) tolerance scan: fixed-step series + adaptive
    series."""
    y0 = np.array([1.0, 0.0, 0.0])
    print("[F1] fixed-step implicit Euler: successively halved h, y2(40) vs "
          "oracle")
    print(f"  {'h':>10} {'y2(40)':>14} {'err(y2)':>12} {'err(inf)':>12}")
    fixed_rows = []
    for N in [100, 200, 400, 800, 1600, 3200]:
        h = 40.0 / N
        t, y = solve_implicit_euler(robertson_rhs, robertson_jacobian,
                                    y0, 0.0, 40.0, h)
        err2 = abs(y[-1, 1] - ref["y40"][1])
        erri = np.max(np.abs(y[-1] - ref["y40"]))
        fixed_rows.append((h, y[-1, 1], err2, erri))
        print(f"  {h:>10.3e} {y[-1,1]:>14.8e} {err2:>12.3e} {erri:>12.3e}")

    print("[F1'] adaptive implicit Euler: successively tightened tol, y2(40) "
          "vs oracle")
    print(f"  {'tol':>9} {'accept':>8} {'y2(40)':>14} {'err(y2)':>12}")
    ad_rows = []
    for tol in [1e-2, 1e-4, 1e-6, 1e-8]:
        t, y, h_used = adaptive_implicit_euler(robertson_rhs,
                                               robertson_jacobian, y0,
                                               0.0, 40.0, 0.1, tol)
        err2 = abs(y[-1, 1] - ref["y40"][1])
        ad_rows.append((tol, len(h_used), y[-1, 1], err2))
        print(f"  {tol:>9.0e} {len(h_used):>8} {y[-1,1]:>14.8e} {err2:>12.3e}")

    # final sanity comparison against the problem-pack orientation value
    print("[F1''] comparison with problem-pack orientation value "
          "(0.7158271, 9.1855e-6, 0.2841637):")
    rel = np.abs(ref["y40"] - PACK_Y40) / PACK_Y40
    print(f"  oracle vs pack value, relative: {rel[0]:.2e}, {rel[1]:.2e}, "
          f"{rel[2]:.2e}")
    return fixed_rows, ad_rows


def cost_vs_accuracy(ref):
    """F2: cost comparison at a matched error target.

    Cost measure (named first): right-hand-side evaluations nfev, Newton
    iterations, wall time. Not counted (identical for both): Jacobian assembly
    and LU factorisation, plotting, interpreter startup.
    Error target: err@40 ~ 7e-5 (adaptive tol=1e-6 just reaches it, see fig05).
    """
    y0 = np.array([1.0, 0.0, 0.0])
    target = 7.0e-5

    # adaptive side (tol=1e-6 gives err@40 = 6.98e-5, from edge_cases E5)
    work_ad = WorkCounter()
    t0 = time.perf_counter()
    tA, yA, hA = adaptive_implicit_euler(robertson_rhs, robertson_jacobian,
                                         y0, 0.0, 40.0, 0.1, 1e-6, work=work_ad)
    t_ad = time.perf_counter() - t0
    err_ad = np.max(np.abs(yA[-1] - ref["y40"]))

    # fixed-step side: find the h reaching the same error (implicit Euler is
    # first order: err ~ C*h; C ~ 3.4e-3 from fig06 suggests h ~ 0.02, then
    # verified by measurement)
    best = None
    for N in [1600, 2000, 2400, 2800, 3200]:
        work_fx = WorkCounter()
        t0 = time.perf_counter()
        tF, yF = solve_implicit_euler(robertson_rhs, robertson_jacobian,
                                      y0, 0.0, 40.0, 40.0 / N, work=work_fx)
        t_fx = time.perf_counter() - t0
        err_fx = np.max(np.abs(yF[-1] - ref["y40"]))
        print(f"  fixed step N={N}: err@40={err_fx:.3e}  "
              f"nfev={work_fx.nfev}  nnewton={work_fx.nnewton}  t={t_fx:.3f}s")
        if err_fx <= target and best is None:
            best = (N, work_fx, t_fx, err_fx)

    print("[F2] cost vs accuracy (error target err@40 <= %.1e; cost measure: "
          "nfev / nnewton / wall time)" % target)
    print(f"  fixed step  N={best[0]}: err={best[3]:.3e}  "
          f"nfev={best[1].nfev}  nnewton={best[1].nnewton}  "
          f"wall={best[2]:.3f}s")
    print(f"  adaptive tol=1e-6:  err={err_ad:.3e}  nfev={work_ad.nfev}  "
          f"nnewton={work_ad.nnewton}  wall={t_ad:.3f}s")
    ratio_nfev = work_ad.nfev / best[1].nfev
    print(f"  => adaptive/fixed nfev ratio = {ratio_nfev:.2f}")
    print("  Negative result (reported honestly): the naive step-doubling "
          "controller solves three implicit steps per accepted step (one "
          "coarse + two half steps); that 3x per-step overhead is partially "
          "absorbed by large relaxed-tail steps, and the measured total nfev "
          "is about %.1f times the fixed-step run. Step doubling's value here "
          "is not saving work - it is not needing to know h in advance and "
          "crossing the transient automatically." % ratio_nfev)

    # figure: error vs nfev (several settings per method, connected curves)
    curve_fixed = []
    for N in [100, 200, 400, 800, 1600, best[0]]:
        work = WorkCounter()
        tF, yF = solve_implicit_euler(robertson_rhs, robertson_jacobian,
                                      y0, 0.0, 40.0, 40.0 / N, work=work)
        curve_fixed.append((work.nfev,
                            np.max(np.abs(yF[-1] - ref["y40"]))))
    curve_ad = []
    for tol in [1e-1, 1e-2, 1e-4, 1e-6, 1e-8]:
        work = WorkCounter()
        tA, yA, hA = adaptive_implicit_euler(robertson_rhs,
                                             robertson_jacobian, y0,
                                             0.0, 40.0, 0.1, tol, work=work)
        curve_ad.append((work.nfev, np.max(np.abs(yA[-1] - ref["y40"]))))
    cf = np.array(curve_fixed)
    ca = np.array(curve_ad)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(cf[:, 0], cf[:, 1], "o-", color="tab:blue",
              label="implicit Euler, fixed h")
    ax.loglog(ca[:, 0], ca[:, 1], "s-", color="tab:orange",
              label="implicit Euler, step-doubling adaptive")
    ax.axhline(target, color="k", ls="--", lw=1,
               label=f"matched accuracy target {target:.0e}")
    ax.set_xlabel("cost: right-hand-side evaluations nfev [---]")
    ax.set_ylabel("global error at t=40")
    ax.set_title("At matched error 7e-5, step-doubling costs ~1.6x the "
                 "f-evaluations\n(negative result reported honestly)")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig07_cost_vs_accuracy.png", dpi=150)
    plt.close(fig)
    return best, (work_ad, t_ad, err_ad)


def traceability_table():
    """F3: report-number traceability table (printed for the report)."""
    print("[F3] report-number traceability table")
    rows = [
        ("reference y(40) = (0.71582707, 9.1855e-6, 0.2841637)",
         "reference.py output [Level 1] (12/17/12 digits confirmed by the two "
         "tolerance runs)"),
        ("|lambda|max(40) = 3.393e3, S(40) = 1.58e5",
         "fig01_eigenvalue_history.png / edge_cases E1 output"),
        ("explicit Euler stability limit h < 5.9e-4 (frozen estimate, "
         "confirmed by measurement)",
         "fig02_euler_step_sweep.png (h >= 6e-4 unstable or negative)"),
        ("RK4 stability limit h < 8.2e-4 (frozen estimate, confirmed)",
         "fig03_rk4_step_sweep.png (error jumps 9 orders of magnitude at "
         "h >= 9e-4)"),
        ("observed orders 0.99±0.00 / 1.00±0.00 / 5.09±0.15 (>= 4)",
         "fig06_convergence_orders.png + convergence_study.py output"),
        ("conservation defect independent of Newton tolerance (algebraic "
         "preservation, e^T J = 0)",
         "fig04_newton_tol_sweep.png / E4 output (defect ~1e-16 for tol "
         "1e-1..1e-14, yet the solution is 28% wrong at tol >= 1e-2)"),
        ("adaptive tol >= 1e-4: controller inactive, err@40 = 1.7e-4",
         "fig05_adaptive_tol_sweep.png / E5 output"),
        ("at matched error 7e-5, adaptive nfev ~ 1.6x fixed step "
         "(measured 9396/6008)",
         "fig07_cost_vs_accuracy.png / F2 output"),
    ]
    for claim, source in rows:
        print(f"  {claim}\n      source: {source}")


def main():
    print("=" * 72)
    print("final_check.py - pre-submission final check (Week 4)")
    print("=" * 72)
    ref = build_reference()
    y2_tolerance_scan(ref)
    print("-" * 72)
    cost_vs_accuracy(ref)
    print("-" * 72)
    traceability_table()
    print("=" * 72)


if __name__ == "__main__":
    main()
