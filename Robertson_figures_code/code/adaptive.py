"""Adaptive step-size control: implicit Euler with step doubling.

Local error estimate for a method of order p=1 from one full step y1 and
two half steps y2:  err ~ || y2 - y1 || / (2^p - 1) = ||y2 - y1||.
Weighted norm with (atol, rtol), solve_ivp style; accept if err <= 1.
Controller: integral (I) controller h *= 0.9 * err^(-1/(p+1)) with
limiters [0.3, 2.0].
"""
import numpy as np
from methods import newton_solve


def adaptive_ie(f, jacf, t_span, y0, rtol=1e-4, atol=1e-9, h0=1e-6,
                newton_tol=1e-10, maxsteps=3_000_000):
    t0, tf = t_span
    y = np.asarray(y0, float); t = t0; h = h0
    ts = [t0]; ys = [y.copy()]; hs = []; errs = []
    nacc = nrej = nfev = 0

    def ie_step(tt, yy, hh):
        F = lambda x: x - yy - hh * f(tt + hh, x)
        JF = lambda x: np.eye(len(yy)) - hh * jacf(tt + hh, x)
        x0 = yy + hh * f(tt, yy)          # explicit predictor as initial guess
        x, it = newton_solve(F, JF, x0, tol=newton_tol)
        return x, it + 1

    while t < tf - 1e-14 and (nacc + nrej) < maxsteps:
        h = min(h, tf - t)
        y1, c1 = ie_step(t, y, h)
        ym, c2 = ie_step(t, y, h / 2)
        y2, c3 = ie_step(t + h / 2, ym, h / 2)
        nfev += c1 + c2 + c3
        err = np.linalg.norm((y2 - y1) / (atol + rtol * np.abs(y2)), np.inf)
        if err <= 1.0:
            t += h; y = y2; nacc += 1
            ts.append(t); ys.append(y.copy()); hs.append(h); errs.append(err)
            fac = min(2.0, max(0.3, 0.9 * (1.0 / max(err, 1e-300)) ** 0.5))
        else:
            nrej += 1
            fac = max(0.3, 0.9 * (1.0 / err) ** 0.5)
        h *= fac
    stats = dict(nacc=nacc, nrej=nrej, nfev=nfev, t_final=t)
    return np.array(ts), np.array(ys), np.array(hs), np.array(errs), stats
