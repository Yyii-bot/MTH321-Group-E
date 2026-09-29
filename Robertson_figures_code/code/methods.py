"""Core time-stepping methods: explicit Euler, classical RK4, implicit Euler
with Newton iteration. All methods return (t, y) history arrays.

Design note: the core integrators are implemented directly (course
requirement); scipy.integrate.solve_ivp is used ONLY as an independent
reference oracle and as a production-solver benchmark (Radau, BDF).
"""
import numpy as np


def euler_explicit(f, t_span, y0, h):
    """Explicit (forward) Euler with uniform step size h."""
    t0, tf = t_span
    y = np.asarray(y0, float)
    ts = [t0]; ys = [y.copy()]; t = t0; nfev = 0
    while t < tf - 1e-15:
        hh = min(h, tf - t)
        y = y + hh * f(t, y); nfev += 1; t += hh
        ts.append(t); ys.append(y.copy())
    return np.array(ts), np.array(ys), nfev


def rk4(f, t_span, y0, h):
    """Classical 4-stage Runge-Kutta method with uniform step size h."""
    t0, tf = t_span
    y = np.asarray(y0, float)
    ts = [t0]; ys = [y.copy()]; t = t0; nfev = 0
    while t < tf - 1e-15:
        hh = min(h, tf - t)
        k1 = f(t, y)
        k2 = f(t + hh / 2, y + hh / 2 * k1)
        k3 = f(t + hh / 2, y + hh / 2 * k2)
        k4 = f(t + hh, y + hh * k3)
        y = y + hh / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        nfev += 4; t += hh
        ts.append(t); ys.append(y.copy())
    return np.array(ts), np.array(ys), nfev


def newton_solve(F, JF, x0, tol=1e-10, maxit=12):
    """Newton iteration for F(x)=0. Stops on a mixed absolute/relative
    correction criterion ||dx|_inf <= tol*(1+||x||_inf)."""
    x = x0.copy()
    for it in range(1, maxit + 1):
        d = np.linalg.solve(JF(x), -F(x))
        x = x + d
        if np.linalg.norm(d, np.inf) <= tol * (1.0 + np.linalg.norm(x, np.inf)):
            return x, it
    return x, maxit


def euler_implicit(f, jacf, t_span, y0, h, newton_tol=1e-10):
    """Implicit (backward) Euler with uniform step h; the nonlinear stage
    equation is solved by Newton's iteration with the analytic Jacobian.
    Initial guess: one explicit-Euler predictor step (extrapolation)."""
    t0, tf = t_span
    y = np.asarray(y0, float); n = y.size
    ts = [t0]; ys = [y.copy()]; t = t0
    nfev = njev = nlu = niter = 0
    while t < tf - 1e-15:
        hh = min(h, tf - t)
        F = lambda x: x - y - hh * f(t + hh, x)
        JF = lambda x: np.eye(n) - hh * jacf(t + hh, x)
        x0 = y + hh * f(t, y); nfev += 1
        y, it = newton_solve(F, JF, x0, tol=newton_tol)
        niter += it; njev += it; nlu += it; nfev += it
        t += hh; ts.append(t); ys.append(y.copy())
    stats = dict(nfev=nfev, njev=njev, nlu=nlu, niter=niter)
    return np.array(ts), np.array(ys), stats
