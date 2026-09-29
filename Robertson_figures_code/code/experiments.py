"""Experiment drivers for the Robertson stiff-ODE study.

Each function returns plain-data dicts/lists so run_all.py can serialise
them (JSON summary + .npz series) and plot the report figures.
"""
import time
import numpy as np
from scipy.integrate import solve_ivp

from model import (rhs, jac, jac_reduced, output_grid, prothero_robinson,
                   Y0, T_FINAL)
from methods import euler_explicit, rk4, euler_implicit
from adaptive import adaptive_ie


# ---------------------------------------------------------------- reference
def reference():
    t_out = output_grid()
    r1 = solve_ivp(rhs, (0, T_FINAL), Y0, method="Radau", jac=jac,
                   rtol=1e-12, atol=1e-14, t_eval=t_out)
    r2 = solve_ivp(rhs, (0, T_FINAL), Y0, method="Radau", jac=jac,
                   rtol=1e-13, atol=1e-15, t_eval=t_out)
    assert r1.success and r2.success
    return dict(t=r1.t, y=r1.y, tighten_diff=float(np.max(np.abs(r1.y - r2.y))),
                y40=r1.y[:, -1])


# ------------------------------------------------------------- eigenanalysis
def eigenanalysis(ref, zero_tol=1e-8):
    t, Y = ref["t"], ref["y"]
    n = t.size
    lmax = np.full(n, np.nan); lmin = np.full(n, np.nan); S = np.full(n, np.nan)
    max_imag = 0.0
    for i in range(n):
        e_full = np.linalg.eigvals(jac(t[i], Y[:, i]))
        e_red = np.linalg.eigvals(jac_reduced(Y[0, i], Y[1, i]))
        max_imag = max(max_imag, float(np.max(np.abs(np.imag(e_red)))))
        nz = np.abs(np.real(e_full)) > zero_tol
        if nz.sum() == 2:
            mags = np.sort(np.abs(np.real(e_full[nz])))
            lmax[i], lmin[i] = mags[-1], mags[0]
            S[i] = mags[-1] / mags[0]
    def at(tq):
        i = int(np.argmin(np.abs(t - tq)))
        return dict(t=float(t[i]), lmax=float(lmax[i]), lmin=float(lmin[i]),
                    S=float(S[i]))
    return dict(t=t, lmax=lmax, lmin=lmin, S=S, max_imag=max_imag,
                checks={str(k): at(k) for k in (1e-4, 1e-2, 1.0, 40.0)})


# ----------------------------------------------------------- order on PR test
def pr_order_study(lam, hs=(0.2, 0.1, 0.05, 0.025, 0.0125), T=1.0):
    f, jf, exact = prothero_robinson(lam)
    rows = []
    for h in hs:
        _, yE, _ = euler_explicit(f, (0, T), [1.0], h)
        _, yR, _ = rk4(f, (0, T), [1.0], h)
        _, yI, _ = euler_implicit(f, jf, (0, T), [1.0], h)
        e = exact(T)
        rows.append(dict(h=h, ee=float(abs(yE[-1, 0] - e)),
                         ie=float(abs(yI[-1, 0] - e)),
                         rk4=float(abs(yR[-1, 0] - e))))
    return rows


# ------------------------------------------------------- threshold sweeps
def weighted_err(y_end, ref_end):
    w = 1e-8 + 1e-6 * np.abs(ref_end)
    return float(np.max(np.abs(y_end - ref_end) / w))


def threshold_sweeps(ref_end):
    def run(method, h):
        with np.errstate(all="ignore"):
            if method == "ee":
                ts, ys, _ = euler_explicit(rhs, (0, T_FINAL), Y0, h)
            else:
                ts, ys, _ = rk4(rhs, (0, T_FINAL), Y0, h)
        fin = bool(np.isfinite(ys).all())
        return dict(h=h, finite=fin,
                    ymax=float(np.nanmax(np.abs(ys))) if fin else None,
                    ymin=float(np.nanmin(ys)) if fin else None,
                    err40=weighted_err(ys[-1], ref_end)
                    if fin and np.nanmax(np.abs(ys)) < 10 else None)
    ee = [run("ee", h) for h in
          (2e-4, 4e-4, 5e-4, 5.5e-4, 5.8e-4, 5.9e-4, 5.95e-4,
           6.0e-4, 6.1e-4, 6.2e-4, 6.3e-4, 6.5e-4)]
    rk = [run("rk4", h) for h in
          (2e-4, 4e-4, 6e-4, 7e-4, 8e-4, 8.2e-4, 8.5e-4, 9e-4, 1e-3)]
    return dict(ee=ee, rk4=rk)


def blowup_times(ref, h_ee=6.5e-4, h_rk=1e-3):
    """Observed instability onset vs frozen-Jacobian crossing prediction."""
    out = {}
    t, lmax = ref["t"], None
    for name, h, fn in (("ee", h_ee, euler_explicit), ("rk4", h_rk, rk4)):
        with np.errstate(all="ignore"):
            ts, ys, _ = fn(rhs, (0, T_FINAL), Y0, h)
        bad = np.where(~np.isfinite(ys).all(axis=1)
                       | (np.abs(ys) > 10).any(axis=1))[0]
        out[name] = dict(h=h, t_fail=float(ts[bad[0]]) if len(bad) else None)
    return out


# ------------------------------------------------ Robertson convergence sweep
def convergence_robertson(ref_end):
    rows = []
    for N in (160, 320, 640, 1280, 2560, 5120):
        _, ys, info = euler_implicit(rhs, jac, (0, T_FINAL), Y0, 40.0 / N)
        rows.append(dict(method="IE", h=40.0 / N, err40=weighted_err(ys[-1], ref_end),
                         ymin=float(np.min(ys)), nfev=info["nfev"]))
    for N in (80000, 160000, 320000, 640000):
        _, ys, _ = euler_explicit(rhs, (0, T_FINAL), Y0, 40.0 / N)
        rows.append(dict(method="EE", h=40.0 / N, err40=weighted_err(ys[-1], ref_end),
                         ymin=float(np.min(ys)), nfev=N))
        _, ys, _ = rk4(rhs, (0, T_FINAL), Y0, 40.0 / N)
        rows.append(dict(method="RK4", h=40.0 / N, err40=weighted_err(ys[-1], ref_end),
                         ymin=float(np.min(ys)), nfev=4 * N))
    return rows


def rk4_short_interval_order():
    """RK4 (and EE/IE) order on the Robertson problem over [0, 0.01],
    where the error is still above the reference-solution noise floor."""
    rS = solve_ivp(rhs, (0, 0.01), Y0, method="Radau", jac=jac,
                   rtol=1e-13, atol=1e-15, t_eval=[0.01])
    wS = 1e-10 + 1e-8 * np.abs(rS.y[:, -1])
    rows = []
    for h in (1e-3, 5e-4, 2.5e-4, 1.25e-4):
        _, yR, _ = rk4(rhs, (0, 0.01), Y0, h)
        _, yE, _ = euler_explicit(rhs, (0, 0.01), Y0, h)
        _, yI, _ = euler_implicit(rhs, jac, (0, 0.01), Y0, h)
        rows.append(dict(h=h,
                         rk4=float(np.max(np.abs(yR[-1] - rS.y[:, -1]) / wS)),
                         ee=float(np.max(np.abs(yE[-1] - rS.y[:, -1]) / wS)),
                         ie=float(np.max(np.abs(yI[-1] - rS.y[:, -1]) / wS))))
    return rows


# ------------------------------------------------------------- adaptivity
def adaptive_study(ref_end, rtols=(1e-3, 1e-4, 1e-5, 1e-6)):
    rows, series = [], {}
    for rt in rtols:
        t0 = time.perf_counter()
        ta, ya, ha, ea, info = adaptive_ie(rhs, jac, (0, T_FINAL), Y0, rtol=rt)
        el = time.perf_counter() - t0
        rows.append(dict(rtol=rt, nacc=info["nacc"], nrej=info["nrej"],
                         nfev=info["nfev"], wall=el,
                         hmin=float(ha.min()), hmax=float(ha.max()),
                         err40=weighted_err(ya[-1], ref_end),
                         y2_40=float(ya[-1, 1])))
        series[f"{rt:.0e}"] = dict(t=ta, h=ha, err=ea, y=ya)
    return dict(rows=rows, series=series)


# ------------------------------------------------- conservation & Newton tol
def newton_tol_sweep(ref_end, h=0.05,
                     tols=(1e-2, 1e-4, 1e-6, 1e-10, 1e-13)):
    rows = []
    for nt in tols:
        ts, ys, info = euler_implicit(rhs, jac, (0, T_FINAL), Y0, h, newton_tol=nt)
        rows.append(dict(newton_tol=nt,
                         max_defect=float(np.max(np.abs(ys.sum(axis=1) - 1.0))),
                         err40=weighted_err(ys[-1], ref_end),
                         mean_newton_iter=info["niter"] / (len(ts) - 1)))
    return rows


def conservation(ref_end, h=1e-4):
    out = {}
    for name, fn in (("EE", euler_explicit), ("RK4", rk4)):
        ts, ys, _ = fn(rhs, (0, T_FINAL), Y0, h)
        out[name] = dict(t=ts, defect=np.abs(ys.sum(axis=1) - 1.0),
                         ymin=float(np.min(ys)), err40=weighted_err(ys[-1], ref_end))
    for name, hh in (("IE", 1e-4), ("IE (h=0.05)", 0.05)):
        ts, ys, _ = euler_implicit(rhs, jac, (0, T_FINAL), Y0, hh)
        out[name] = dict(t=ts, defect=np.abs(ys.sum(axis=1) - 1.0),
                         ymin=float(np.min(ys)), err40=weighted_err(ys[-1], ref_end))
    return out


# ------------------------------------------------------- work vs accuracy
def work_precision(ref_end):
    def timeit(fn, repeats):
        best = np.inf; res = None
        for _ in range(repeats):
            t0 = time.perf_counter(); res = fn(); best = min(best, time.perf_counter() - t0)
        return best, res
    pts = []
    for N, rep in ((80000, 2), (160000, 2), (320000, 2), (640000, 1)):
        el, (_, ys, _) = timeit(lambda: euler_explicit(rhs, (0, T_FINAL), Y0, 40.0 / N), rep)
        pts.append(dict(method="EE", cost=el, err40=weighted_err(ys[-1], ref_end), steps=N))
        el, (_, ys, _) = timeit(lambda: rk4(rhs, (0, T_FINAL), Y0, 40.0 / N), rep)
        pts.append(dict(method="RK4", cost=el, err40=weighted_err(ys[-1], ref_end), steps=N))
    for N, rep in ((160, 5), (320, 5), (640, 5), (1280, 3), (2560, 2), (5120, 2)):
        el, (_, ys, _) = timeit(lambda: euler_implicit(rhs, jac, (0, T_FINAL), Y0, 40.0 / N), rep)
        pts.append(dict(method="IE uniform", cost=el, err40=weighted_err(ys[-1], ref_end), steps=N))
    for rt, rep in ((1e-3, 5), (1e-4, 5), (1e-5, 3), (1e-6, 2)):
        el, (_, ya, *_ ) = timeit(lambda: adaptive_ie(rhs, jac, (0, T_FINAL), Y0, rtol=rt)[:2], rep)
        pts.append(dict(method="IE adaptive", cost=el, err40=weighted_err(ya[-1], ref_end), steps=None))
    for name in ("Radau", "BDF"):
        for rt, rep in ((1e-4, 5), (1e-6, 5), (1e-8, 3), (1e-10, 2)):
            el, sol = timeit(lambda: solve_ivp(rhs, (0, T_FINAL), Y0, method=name,
                                               jac=jac, rtol=rt, atol=1e-12), rep)
            pts.append(dict(method=name, cost=el, err40=weighted_err(sol.y[:, -1], ref_end),
                            steps=sol.t.size - 1))
    return pts


# ------------------------------------------------------------- y2(40) table
def threshold_trajectories():
    """Down-sampled EE trajectories just below/above the predicted stability
    bound, for the late-time instability figure."""
    out = {}
    for h, name in ((5.9e-4, "ee_stable"), (6.1e-4, "ee_marginal")):
        with np.errstate(all="ignore"):
            ts, ys, _ = euler_explicit(rhs, (0, T_FINAL), Y0, h)
        sel = np.linspace(0, len(ts) - 1, 4000).astype(int)
        out[name] = dict(t=ts[sel], y=ys[sel])
    return out


def y2_table(ref_end):
    rows = []
    for rt in (1e-3, 1e-4, 1e-5, 1e-6):
        _, ya, _, _, _ = adaptive_ie(rhs, jac, (0, T_FINAL), Y0, rtol=rt)
        rows.append(dict(method=f"IE adaptive rtol={rt:.0e}", y2=float(ya[-1, 1])))
    for N in (160, 640, 2560):
        _, ys, _ = euler_implicit(rhs, jac, (0, T_FINAL), Y0, 40.0 / N)
        rows.append(dict(method=f"IE uniform N={N}", y2=float(ys[-1, 1])))
    for rt in (1e-6, 1e-9, 1e-12):
        sol = solve_ivp(rhs, (0, T_FINAL), Y0, method="Radau", jac=jac, rtol=rt, atol=1e-14)
        rows.append(dict(method=f"Radau rtol={rt:.0e}", y2=float(sol.y[1, -1])))
    rows.append(dict(method="reference (Radau 1e-12/1e-14)", y2=float(ref_end[1])))
    return rows
