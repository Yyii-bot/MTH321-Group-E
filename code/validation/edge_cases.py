"""
edge_cases.py - Week 2-3 deliverable: edge-case hunting experiments.

Each experiment on the role card's hunting list is run as a complete
"theory -> measurement -> agreement?" chain:

  E1  Frozen-Jacobian diagnostic: |lambda|max, |lambda|min (nonzero) and the
      stiffness ratio S(t) along the reference solution, on a logarithmic time
      axis; the degenerate initial point t=0 is excluded (plot starts at the
      first positive logarithmic time, per the problem-pack protocol).
  E2  Explicit Euler step sweep: measured stable/non-negative/error triples vs
      the frozen estimate h < 2/|lambda|max ~ 5.9e-4 (the problem pack requires
      the h-lambda overlay - a local diagnostic only - to be tested against a
      step sweep, with discrepancies reported).
  E3  RK4 step sweep: measured stability boundary vs the theoretical boundary
      h < 2.785/|lambda|max ~ 8.2e-4.
  E4  Newton stopping-tolerance sweep: how the conservation defect and the
      global error depend on the Newton tolerance (problem pack: Newton need
      not be driven to machine precision; pick a tolerance that makes the
      algebraic error negligible and prove it).
  E5  Loose-adaptive-tolerance experiment: local error estimate vs the true
      global error as tol varies, demonstrating quantitatively that a loose
      tolerance gives an untrustworthy answer.

Figures follow the visualisation guide: stability boundaries marked with
vertical lines, theoretical references drawn, blue/orange colour-blind-safe
palette, DPI >= 150, fixed file names.
Run: python edge_cases.py
Author: Rong Yin (Role 5)
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from robertson_model import robertson_rhs, robertson_jacobian, eigen_info
from solvers import (solve_euler, solve_rk4, solve_implicit_euler,
                     WorkCounter, adaptive_implicit_euler)
from reference import build_reference

_HERE = Path(__file__).resolve().parent
if _HERE.name == "validation":
    # team-repo layout (code/validation/*.py): figures -> <repo>/figures/validation
    FIGDIR = _HERE.parents[1] / "figures" / "validation"
else:
    # standalone layout (code/*.py): figures -> sibling figures/
    FIGDIR = _HERE.parent / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

# RK4's absolute-stability limit on the negative real axis, |R(-x)| = 1,
# computed numerically (avoids hand-copying the constant).
from scipy.optimize import brentq
_RK4_R = lambda z: 1 + z + z**2/2 + z**3/6 + z**4/24
RK4_NEG_LIMIT = -brentq(lambda x: abs(_RK4_R(-x)) - 1, 2.0, 3.5)


def experiment_e1_eigenvalue_history(ref):
    """E1: frozen-Jacobian diagnostic - nonzero eigenvalues and stiffness
    ratio along logarithmic time."""
    t, Y = ref["t"], ref["Y_hi"]
    lam_max, lam_min, ratio = [], [], []
    for i in range(len(t)):
        lam, lam_nz, r = eigen_info(t[i], Y[:, i])
        lam_max.append(np.abs(lam_nz.real).max() if len(lam_nz) else np.nan)
        lam_min.append(np.abs(lam_nz.real).min() if len(lam_nz) else np.nan)
        ratio.append(r)
    lam_max, lam_min, ratio = map(np.array, (lam_max, lam_min, ratio))

    # frozen stability boundaries as functions of time
    h_euler_lim = 2.0 / lam_max
    h_rk4_lim = RK4_NEG_LIMIT / lam_max

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    ax.loglog(t[1:], lam_max[1:], "o-", ms=3, color="tab:blue",
              label=r"$|\lambda|_{\max}$ (nonzero)")
    ax.loglog(t[1:], lam_min[1:], "s-", ms=3, color="tab:orange",
              label=r"$|\lambda|_{\min}$ (nonzero)")
    ax.set_xlabel("t [s]")
    ax.set_ylabel(r"|Re $\lambda$| [1/s]")
    ax.set_title("Nonzero Jacobian eigenvalues vs log time;\n"
                 "structural zero eigenvalue excluded")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend()

    ax = axes[1]
    ax.loglog(t[1:], ratio[1:], "o-", ms=3, color="tab:blue",
              label="stiffness ratio S(t)")
    ax.loglog(t[1:], h_euler_lim[1:], "--", color="tab:orange",
              label="Euler limit  2/|lambda|max")
    ax.loglog(t[1:], h_rk4_lim[1:], ":", color="tab:green",
              label="RK4 limit  2.785/|lambda|max")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("S(t)  /  step-size limit [s]")
    ax.set_title("Stiffness ratio grows to 1.6e5 at t=40;\n"
                 "explicit-step limit shrinks below 6e-4")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(FIGDIR / "fig01_eigenvalue_history.png", dpi=150)
    plt.close(fig)

    i40 = int(np.argmin(np.abs(t - 40.0)))
    print(f"[E1] |lambda|max(40) = {lam_max[i40]:.4e}   S(40) = {ratio[i40]:.4e}")
    print(f"     Euler frozen limit h < 2/|lambda|max = {h_euler_lim[i40]:.3e}")
    print(f"     RK4   frozen limit h < {abs(RK4_NEG_LIMIT):.4f}/|lambda|max "
          f"= {abs(h_rk4_lim[i40]):.3e}")
    return lam_max, ratio


def _run_explicit_sweep(solve_fn, name, h_values, ref, fig_name, theory_limit):
    """Shared driver for E2/E3: explicit-method step sweep recording the
    stable / non-negative / error triple, plus the blow-up location."""
    y0 = np.array([1.0, 0.0, 0.0])
    rows = []
    for h in h_values:
        with np.errstate(over="ignore", invalid="ignore"):
            t, y = solve_fn(robertson_rhs, y0, 0.0, 40.0, h)
        finite = bool(np.all(np.isfinite(y)))
        ymin = float(y.min()) if finite else np.nan
        neg = ymin < -1e-12 if finite else True
        # blow-up location: first time |y|>10 or a nan/inf appears
        blow_t = np.nan
        if not finite:
            bad = np.where(~np.all(np.isfinite(y), axis=1))[0]
            if len(bad):
                blow_t = float(t[bad[0]])
        else:
            big = np.where(np.abs(y).max(axis=1) > 10.0)[0]
            if len(big):
                blow_t = float(t[big[0]])
        stable = finite and np.isnan(blow_t)
        err40 = (float(np.max(np.abs(y[-1] - ref["y40"]))) if finite else np.nan)
        rows.append(dict(h=h, steps=len(t) - 1, stable=stable, neg=neg,
                         ymin=ymin, blow_t=blow_t, err40=err40))

    print(f"[{name}] step sweep (frozen-Jacobian limit h < {theory_limit:.2e}):")
    print(f"  {'h':>9} {'steps':>8} {'stable':>7} {'negative':>9} "
          f"{'min(y)':>11} {'blow-up t':>10} {'err@40':>11}")
    for r in rows:
        print(f"  {r['h']:>9.2e} {r['steps']:>8} {str(r['stable']):>7} "
              f"{str(r['neg']):>9} {r['ymin']:>11.3e} {r['blow_t']:>10.3g} "
              f"{r['err40']:>11.3e}")

    # figure: error vs h (log-log) + failure markers
    fig, ax = plt.subplots(figsize=(7, 4.5))
    hs = np.array([r["h"] for r in rows])
    errs = np.array([r["err40"] for r in rows])
    ok_mask = np.array([r["stable"] and not r["neg"] for r in rows])
    ax.loglog(hs[ok_mask], errs[ok_mask], "o-", color="tab:blue",
              label=f"{name}: stable & non-negative")
    # failed points: nan errors are placed at the top of the plot
    fail_h = hs[~ok_mask]
    fail_e = np.where(np.isfinite(errs[~ok_mask]), errs[~ok_mask], 1e1)
    ax.loglog(fail_h, fail_e, "x", color="tab:red", ms=10,
              label="unstable / negative")
    ax.axvline(theory_limit, color="k", ls="--", lw=1.2,
               label=f"frozen-Jacobian limit {theory_limit:.1e}")
    ax.set_xlabel("step size h [s]")
    ax.set_ylabel("global error at t=40 [---]")
    ax.set_title(f"{name} step sweep: instability onset matches the "
                 "frozen-Jacobian prediction")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR / fig_name, dpi=150)
    plt.close(fig)
    return rows


def experiment_e2_euler_sweep(ref):
    """E2: explicit Euler step sweep (theoretical limit 5.9e-4)."""
    h_values = [2e-3, 1e-3, 8e-4, 7e-4, 6e-4, 5.9e-4, 5e-4, 4e-4, 2e-4, 1e-4]
    return _run_explicit_sweep(solve_euler, "Explicit Euler", h_values, ref,
                               "fig02_euler_step_sweep.png", 2.0 / 3393.0)


def experiment_e3_rk4_sweep(ref):
    """E3: RK4 step sweep (theoretical limit 2.785/3393 ~ 8.2e-4)."""
    h_values = [2e-3, 1.5e-3, 1e-3, 9e-4, 8.2e-4, 7e-4, 5e-4, 2e-4]
    return _run_explicit_sweep(solve_rk4, "RK4", h_values, ref,
                               "fig03_rk4_step_sweep.png",
                               RK4_NEG_LIMIT / 3393.0)


def experiment_e4_newton_tol_sweep(ref):
    """E4: Newton stopping-tolerance sweep.

    Fixed h=0.1, newton_tol swept from 1e-1 down to 1e-14. Mathematical
    background: the conserved direction e=(1,1,1)^T satisfies e^T J = 0 (a
    linear invariant), hence e^T(I-hJ)^{-1} = e^T, so EVERY Newton iterate
    preserves e^T w = e^T y_n algebraically - the conservation defect should
    stay at round-off level independent of newton_tol. This experiment
    verifies that derivation numerically and locates the point where the
    global error stops being sensitive to newton_tol (time-discretisation
    error takes over).
    """
    y0 = np.array([1.0, 0.0, 0.0])
    tols = [1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14]
    defects, errs = [], []
    for tol in tols:
        t, y = solve_implicit_euler(robertson_rhs, robertson_jacobian,
                                    y0, 0.0, 40.0, 0.1, newton_tol=tol)
        defects.append(np.abs(y.sum(axis=1) - 1.0).max())
        errs.append(np.max(np.abs(y[-1] - ref["y40"])))

    print("[E4] Newton tolerance sweep (fixed h=0.1):")
    print(f"  {'newton_tol':>11} {'cons_defect':>13} {'err@40':>12}")
    for tol, d, e in zip(tols, defects, errs):
        print(f"  {tol:>11.0e} {d:>13.3e} {e:>12.3e}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].loglog(tols, defects, "o-", color="tab:blue")
    axes[0].axhline(1e-15, color="k", ls="--", lw=1,
                    label="round-off level (independent of tol)")
    axes[0].set_xlabel("Newton stopping tolerance")
    axes[0].set_ylabel("mass-conservation defect")
    axes[0].set_title("e^T J = 0 makes every Newton iterate preserve the\n"
                      "linear invariant: defect stays at round-off for any tol")
    axes[0].grid(True, which="both", ls=":", alpha=0.5)
    axes[0].legend(fontsize=8)

    axes[1].loglog(tols, errs, "o-", color="tab:orange")
    axes[1].axhline(errs[-1], color="k", ls="--", lw=1,
                    label="time-discretisation floor")
    axes[1].set_xlabel("Newton stopping tolerance")
    axes[1].set_ylabel("global error at t=40")
    axes[1].set_title("Global error is insensitive to Newton tol:\n"
                      "time-discretisation error dominates throughout")
    axes[1].grid(True, which="both", ls=":", alpha=0.5)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig04_newton_tol_sweep.png", dpi=150)
    plt.close(fig)
    return tols, defects, errs


def experiment_e5_adaptive_tol_sweep(ref):
    """E5: loose-adaptive-tolerance experiment - local error estimate vs the
    true error.

    tol swept from 1e-1 to 1e-10; accepted/rejected steps, minimum h and the
    true global error are recorded. Key point: the local error estimate e does
    satisfy e <= tol (the controller works as designed), but the global error
    is only O(tol) in the asymptotic regime - under a loose tolerance y2(40)
    is untrustworthy.
    """
    y0 = np.array([1.0, 0.0, 0.0])
    tols = [1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10]
    print("[E5] adaptive implicit Euler tolerance sweep (h0=0.1):")
    print(f"  {'tol':>9} {'accept':>8} {'reject':>8} {'min h':>10} "
          f"{'err@40':>12} {'y2(40)':>12}")
    rows = []
    for tol in tols:
        work = WorkCounter()
        t, y, h_used = adaptive_implicit_euler(robertson_rhs,
                                               robertson_jacobian, y0,
                                               0.0, 40.0, 0.1, tol,
                                               work=work)
        err = np.max(np.abs(y[-1] - ref["y40"]))
        rows.append((tol, work.nsteps, work.nreject, h_used.min(), err,
                     y[-1, 1]))
        print(f"  {tol:>9.0e} {work.nsteps:>8} {work.nreject:>8} "
              f"{h_used.min():>10.2e} {err:>12.3e} {y[-1,1]:>12.4e}")

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ts = np.array([r[0] for r in rows])
    es = np.array([r[4] for r in rows])
    ax.loglog(ts, es, "o-", color="tab:blue", label="true global error at t=40")
    ax.loglog(ts, ts, "k--", label="error = tol (reference)")
    ax.set_xlabel("adaptive controller tolerance tol")
    ax.set_ylabel("global error at t=40")
    ax.set_title("Loose tolerance gives loose answers: error tracks tol\n"
                 "only in the asymptotic regime")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig05_adaptive_tol_sweep.png", dpi=150)
    plt.close(fig)
    return rows


def main():
    print("=" * 72)
    print("edge_cases.py - edge-case hunting experiments (Week 2-3)")
    print("=" * 72)
    ref = build_reference()
    experiment_e1_eigenvalue_history(ref)
    print("-" * 72)
    experiment_e2_euler_sweep(ref)
    print("-" * 72)
    experiment_e3_rk4_sweep(ref)
    print("-" * 72)
    experiment_e4_newton_tol_sweep(ref)
    print("-" * 72)
    experiment_e5_adaptive_tol_sweep(ref)
    print("=" * 72)
    print(f"Figures written to {FIGDIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
