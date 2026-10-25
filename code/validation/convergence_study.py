"""
convergence_study.py - Week 3 deliverable: convergence-order verification
(>= 2 methods, with every visualisation-guide Rule 1 element).

Error-vs-h studies on the Robertson problem for three methods (full-state error
at t=40 against the team-built oracle):

  implicit Euler  h in {0.625, 0.3125, ..., 0.0195}   theory: order 1
  explicit Euler  h in {5e-4, 2.5e-4, 1.25e-4, 6.25e-5}
                  (stability limit 5.9e-4 forces steps this small)
  RK4             h in {8e-4 ... 4.5e-4}
                  (stability limit 8.2e-4; a 4th-order method hits the oracle
                  measurement floor almost immediately - the floor itself must
                  be annotated per the guide: it is a finding, not noise)

Every panel must contain (guide Rules 1/2/4): log-log axes, the theoretical
slope reference line, the polyfit(cov=True) slope +/- standard error, a 95%
confidence band, a round-off-floor annotation, and a conclusion-statement title.
Run: python convergence_study.py
Author: Rong Yin (Role 5)
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from robertson_model import robertson_rhs, robertson_jacobian
from solvers import solve_implicit_euler, solve_euler, solve_rk4
from reference import build_reference

_HERE = Path(__file__).resolve().parent
if _HERE.name == "validation":
    # team-repo layout (code/validation/*.py): figures -> <repo>/figures/validation
    FIGDIR = _HERE.parents[1] / "figures" / "validation"
else:
    # standalone layout (code/*.py): figures -> sibling figures/
    FIGDIR = _HERE.parent / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)


def fit_order(h_values, errors):
    """Log-log linear fit + slope standard error (the guide's polyfit recipe).

    Only points above the floor are used (the round-off floor would pollute
    the slope if all points were included); the floor is the first point where
    the error stops decreasing with h.
    """
    h_values = np.asarray(h_values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    floor_idx = int(np.argmin(errors))
    n_fit = max(3, floor_idx + 1) if floor_idx >= 2 else len(h_values)
    # no interior floor (minimum at the last point and still decreasing):
    # use all points
    if floor_idx == len(errors) - 1 and errors[-1] < errors[-2]:
        n_fit = len(errors)
    hs, es = h_values[:n_fit], errors[:n_fit]
    coef, cov = np.polyfit(np.log(hs), np.log(es), 1, cov=True)
    return coef[0], np.sqrt(cov[0, 0]), coef, cov


def plot_method(ax, name, h_values, errors, theory_order, color):
    """One convergence panel: data + theoretical slope + fit + 95% band +
    floor annotation."""
    h_values = np.asarray(h_values)
    errors = np.asarray(errors)
    order, order_err, coef, cov = fit_order(h_values, errors)

    ax.loglog(h_values, errors, "o-", color=color, label=name)

    # theoretical-slope reference line (guide: an order claim without the
    # reference line is a claim without a benchmark)
    h_lin = np.logspace(np.log10(h_values.min()), np.log10(h_values.max()), 50)
    c_ref = errors[0] / h_values[0] ** theory_order
    ax.loglog(h_lin, c_ref * h_lin ** theory_order, "k--",
              label=f"order-{theory_order} reference")

    # fitted line + 95% confidence band
    lx = np.linspace(np.log(h_values.min()), np.log(h_values.max()), 50)
    line = coef[0] * lx + coef[1]
    se = np.sqrt(cov[1, 1] + 2 * lx * cov[0, 1] + lx ** 2 * cov[0, 0])
    ax.loglog(np.exp(lx), np.exp(line), "--", color=color, alpha=0.7,
              label=f"fit {order:.2f} ± {order_err:.2f}")
    ax.fill_between(np.exp(lx), np.exp(line - 1.96 * se),
                    np.exp(line + 1.96 * se), color=color, alpha=0.15)

    # round-off-floor annotation (guide Rule 4: a floor is a finding - point
    # at it)
    floor_idx = int(np.argmin(errors))
    if floor_idx < len(errors) - 1 or errors[floor_idx] < errors[0] * 1e-6:
        ax.annotate("round-off floor",
                    xy=(h_values[floor_idx], errors[floor_idx]),
                    xytext=(h_values[floor_idx] * 8, errors[floor_idx] * 30),
                    arrowprops=dict(arrowstyle="->"), fontsize=8)

    verdict = ("consistent (>= design order)"
               if order >= theory_order - 2 * max(order_err, 0.05)
               else "MISMATCH")
    ax.set_title(f"{name}: order {order:.2f} ± {order_err:.2f} "
                 f"(theory {theory_order})", fontsize=11)
    ax.text(0.03, 0.06, verdict, transform=ax.transAxes, fontsize=9,
            color="darkgreen" if verdict.startswith("consistent") else "red")
    ax.set_xlabel("step size h [s]")
    ax.set_ylabel("global error at t=40")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(fontsize=7, loc="lower right")
    return order, order_err


def study_implicit_euler(ref):
    """Implicit Euler order study: successively halved h, unconditionally
    stable throughout."""
    y0 = np.array([1.0, 0.0, 0.0])
    N_values = [64, 128, 256, 512, 1024, 2048]
    h_values = [40.0 / N for N in N_values]
    errors = []
    for h in h_values:
        t, y = solve_implicit_euler(robertson_rhs, robertson_jacobian,
                                    y0, 0.0, 40.0, h)
        errors.append(np.max(np.abs(y[-1] - ref["y40"])))
    return h_values, np.array(errors)


def study_explicit_euler(ref):
    """Explicit Euler order study: every h must lie below the frozen stability
    limit 5.9e-4, otherwise the run goes unstable."""
    y0 = np.array([1.0, 0.0, 0.0])
    h_values = [5e-4, 2.5e-4, 1.25e-4, 6.25e-5]
    errors = []
    for h in h_values:
        t, y = solve_euler(robertson_rhs, y0, 0.0, 40.0, h)
        errors.append(np.max(np.abs(y[-1] - ref["y40"])))
    return h_values, np.array(errors)


def study_rk4(ref):
    """RK4 order study.

    The usable h-window is narrow: above, the stability limit 8.2e-4 (larger h
    goes unstable); below, the oracle measurement floor ~1e-13 (below it we
    measure the reference's own tolerance jitter). Points are taken inside the
    window. Note: within this window the effective slope is slightly steeper
    than 4 because the error is dominated by the stiff-transient resolution,
    so the verified claim is "at least the design order 4" - "exactly 4.00" is
    not identifiable on this problem and this limitation must be stated in the
    report.
    """
    y0 = np.array([1.0, 0.0, 0.0])
    h_values = [8e-4, 6e-4, 5e-4, 4.5e-4]
    errors = []
    for h in h_values:
        t, y = solve_rk4(robertson_rhs, y0, 0.0, 40.0, h)
        errors.append(np.max(np.abs(y[-1] - ref["y40"])))
    return h_values, np.array(errors)


def main():
    print("=" * 72)
    print("convergence_study.py - convergence-order verification (Week 3)")
    print("=" * 72)
    ref = build_reference()

    studies = [
        ("Implicit Euler", study_implicit_euler, 1, "tab:blue"),
        ("Explicit Euler", study_explicit_euler, 1, "tab:orange"),
        ("RK4", study_rk4, 4, "tab:green"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    summary = []
    for ax, (name, study, p_theory, color) in zip(axes, studies):
        h_values, errors = study(ref)
        order, order_err = plot_method(ax, name, h_values, errors,
                                       p_theory, color)
        summary.append((name, h_values, errors, order, order_err))
        print(f"[{name}] observed order = {order:.2f} ± {order_err:.2f} "
              f"(theory {p_theory})")
        for h, e in zip(h_values, errors):
            print(f"    h={h:.3e}   err@40={e:.3e}")

    fig.tight_layout()
    fig.savefig(FIGDIR / "fig06_convergence_orders.png", dpi=150)
    plt.close(fig)

    print("=" * 72)
    print(f"Figure written to {FIGDIR / 'fig06_convergence_orders.png'}")
    print("=" * 72)
    return summary


if __name__ == "__main__":
    main()
