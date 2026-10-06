"""
solvers.py - Validation-owned solver collection (Role 5: Testing & Validation).

Note: the Algorithm Implementation role (Role 3) produces the team's production
solvers; for independent validation I re-implemented the same mathematical
definitions here. If two independent implementations agree in the fine-mesh
limit, that is exactly the self-consistency check required by the project brief.

Contents:
1. Explicit Euler (order 1)
2. Classical RK4 (order 4)
3. Implicit Euler + Newton iteration (Newton tolerance adjustable for the
   tolerance-sweep experiments)
4. Step-doubling adaptive implicit Euler (local-error-estimate-driven h)

Every solver records a "work ledger" (f evaluations, Newton iterations) for the
cost-vs-accuracy experiments.
Author: Rong Yin (Role 5)
"""

import numpy as np


# ---------------------------------------------------------------------------
# Work ledger: cost counters shared by all solvers (visualisation guide
# Rule 7 requires naming the cost measure before quoting any number)
# ---------------------------------------------------------------------------

class WorkCounter:
    """Lightweight cost counter.

    nfev    : right-hand-side f evaluations
    njeval  : analytic Jacobian evaluations
    nnewton : accepted Newton iterations (total)
    nsteps  : accepted time steps
    nreject : rejected time steps (adaptive methods only)
    """

    def __init__(self):
        self.nfev = 0
        self.njeval = 0
        self.nnewton = 0
        self.nsteps = 0
        self.nreject = 0


# ---------------------------------------------------------------------------
# Explicit Euler: y_{n+1} = y_n + h f(t_n, y_n)
# ---------------------------------------------------------------------------

def euler_step(f, t, y, h, work=None):
    """One explicit Euler step. Returns y_{n+1}."""
    if work is not None:
        work.nfev += 1
    return y + h * f(t, y)


def solve_euler(f, y0, t0, t1, h, work=None):
    """Fixed-step explicit Euler integration to t1 (final step shortened to
    land exactly on t1). Returns (t, Y): time grid (n+1,) and trajectory
    (n+1, d)."""
    n_steps = int(np.ceil((t1 - t0) / h))
    t = np.empty(n_steps + 1)
    t[0] = t0
    y = np.zeros((n_steps + 1,) + np.shape(y0))
    y[0] = y0
    for n in range(n_steps):
        step = min(h, t1 - t[n])
        t[n + 1] = t[n] + step
        y[n + 1] = euler_step(f, t[n], y[n], step, work)
    if work is not None:
        work.nsteps += n_steps
    return t, y


# ---------------------------------------------------------------------------
# Classical RK4 (order 4)
# ---------------------------------------------------------------------------

def rk4_step(f, t, y, h, work=None):
    """One classical four-stage fourth-order Runge-Kutta step. Returns y_{n+1}."""
    if work is not None:
        work.nfev += 4
    k1 = f(t, y)
    k2 = f(t + h / 2, y + h * k1 / 2)
    k3 = f(t + h / 2, y + h * k2 / 2)
    k4 = f(t + h, y + h * k3)
    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0


def solve_rk4(f, y0, t0, t1, h, work=None):
    """Fixed-step RK4 integration to t1. Returns (t, Y)."""
    n_steps = int(np.ceil((t1 - t0) / h))
    t = np.empty(n_steps + 1)
    t[0] = t0
    y = np.zeros((n_steps + 1,) + np.shape(y0))
    y[0] = y0
    for n in range(n_steps):
        step = min(h, t1 - t[n])
        t[n + 1] = t[n] + step
        y[n + 1] = rk4_step(f, t[n], y[n], step, work)
    if work is not None:
        work.nsteps += n_steps
    return t, y


# ---------------------------------------------------------------------------
# Implicit Euler + Newton iteration
#    y_{n+1} = y_n + h f(t_{n+1}, y_{n+1})
#    F(w) = w - y_n - h f(t_{n+1}, w) = 0,  J_F = I - h J_f
# ---------------------------------------------------------------------------

def newton_solve_ie(f, jac, t_next, y_n, h, tol=1e-12, max_iter=50,
                    work=None, trace=False):
    """Newton iteration for one implicit Euler step.

    Initial guess w0 = y_n (warm start). Returns (y_next, converged, iterates)
    where iterates is the list of residual norms per iteration (for analysing
    the residual trajectory in the "bad Newton initial guess" edge case).
    """
    w = y_n.copy()
    iterates = []
    for _ in range(max_iter):
        F = w - y_n - h * f(t_next, w)
        if work is not None:
            work.nfev += 1
        r = np.linalg.norm(F, np.inf)
        iterates.append(r)
        if trace:
            print(f"        ||F||_inf = {r:.3e}")
        if r < tol:
            return w, True, iterates
        JF = np.eye(len(w)) - h * jac(t_next, w)
        if work is not None:
            work.njeval += 1
        delta = np.linalg.solve(JF, F)
        w = w - delta          # full Newton step (no damping needed with warm start here)
        if work is not None:
            work.nnewton += 1
    return w, False, iterates


def implicit_euler_step(f, jac, t, y, h, newton_tol=1e-12, work=None,
                        trace=False):
    """One implicit Euler step (Newton solve for the nonlinear system)."""
    y_next, ok, it = newton_solve_ie(f, jac, t + h, y, h, tol=newton_tol,
                                     work=work, trace=trace)
    if not ok:
        raise RuntimeError(f"Newton failed to converge at t={t:.6g} (h={h:.3e})")
    return y_next


def solve_implicit_euler(f, jac, y0, t0, t1, h, newton_tol=1e-12, work=None):
    """Fixed-step implicit Euler integration to t1. Returns (t, Y)."""
    n_steps = int(np.ceil((t1 - t0) / h))
    t = np.empty(n_steps + 1)
    t[0] = t0
    y = np.zeros((n_steps + 1,) + np.shape(y0))
    y[0] = y0
    for n in range(n_steps):
        step = min(h, t1 - t[n])
        t[n + 1] = t[n] + step
        y[n + 1] = implicit_euler_step(f, jac, t[n], y[n], step,
                                       newton_tol=newton_tol, work=work)
    if work is not None:
        work.nsteps += n_steps
    return t, y


# ---------------------------------------------------------------------------
# Step-doubling adaptive implicit Euler
#    y_coarse = one step of size h; y_fine = two steps of size h/2;
#    e = ||y_fine - y_coarse|| (an O(h^2) local-error estimate)
#    accept if e <= tol (take y_fine); halve h if e > tol; double h (capped at
#    h0) if e < tol/10.
# ---------------------------------------------------------------------------

def adaptive_implicit_euler(f, jac, y0, t0, t1, h0, tol, newton_tol=1e-12,
                            work=None):
    """Adaptive implicit Euler. Returns (t, Y, h_used)."""
    t = [t0]
    y = [np.array(y0, dtype=float)]
    h_used = []
    h = h0
    while t[-1] < t1:
        h = min(h, t1 - t[-1])              # never overshoot the final time
        y_coarse = implicit_euler_step(f, jac, t[-1], y[-1], h,
                                       newton_tol=newton_tol, work=work)
        y_half = implicit_euler_step(f, jac, t[-1], y[-1], h / 2,
                                     newton_tol=newton_tol, work=work)
        y_fine = implicit_euler_step(f, jac, t[-1] + h / 2, y_half, h / 2,
                                     newton_tol=newton_tol, work=work)
        e = np.linalg.norm(y_fine - y_coarse, np.inf)
        if e <= tol:
            t.append(t[-1] + h)
            y.append(y_fine)                # accept the more accurate candidate
            h_used.append(h)
            if work is not None:
                work.nsteps += 1
            if e < tol / 10:
                h = min(2 * h, h0)          # well under budget: double the step
        else:
            h /= 2                          # over budget: halve and retry
            if work is not None:
                work.nreject += 1
    return np.array(t), np.array(y), np.array(h_used)
