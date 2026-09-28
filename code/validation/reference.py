"""
reference.py - Team-built high-accuracy reference solution(oracle) - Week 1.
Reference solution must be constructed and endorsed by ourselves. This module implements a three-level credibility chain:

Level 1  solve_ivp(method='Radau', rtol=1e-11, atol=1e-13) on the protocol grid;
         repeated with rtol=1e-13, atol=1e-15 (tolerances 100x tighter) and
         compared point by point.
Level 2  Cross-family check: LSODA and BDF (guards against a single solver's
         systematic bias).
Level 3  Sanity checks: y(40) vs the problem pack's orientation value; the mass
         invariant must hold at round-off level.

Only digits confirmed by both Level-1 runs are quoted (problem-pack protocol).

Author: Rong Yin (Role 5: Testing & Validation)
"""

import numpy as np
from scipy.integrate import solve_ivp

from robertson_model import robertson_rhs, protocol_grid


def build_reference(t_grid=None):
    """Build the three-level-validated reference solution.

    Parameters
    ----------
    t_grid : ndarray or None   output grid; None uses the protocol grid

    Returns
    -------
    ref : dict with keys
        t          : output grid
        Y_lo       : Radau rtol=1e-11 solution (row format, shape (3, n))
        Y_hi       : Radau rtol=1e-13 (100x tighter) solution
        diff       : pointwise absolute difference |Y_lo - Y_hi|
        max_diff   : maximum pointwise difference over the grid
        y40        : y(40) from the tighter run
        y40_digits : decimal digits of y(40) confirmed by BOTH runs, per component
    """
    if t_grid is None:
        t_grid = protocol_grid()

    y0 = [1.0, 0.0, 0.0]

    # Level 1: tight-tolerance Radau, two tolerance sets 100x apart
    # (problem-pack reproducible protocol; Radau's rtol floor is ~2.2e-14,
    # hence rtol = 1e-11 vs 1e-13)
    sol_lo = solve_ivp(robertson_rhs, [t_grid[0], t_grid[-1]], y0,
                       method="Radau", t_eval=t_grid,
                       rtol=1e-11, atol=1e-13)
    sol_hi = solve_ivp(robertson_rhs, [t_grid[0], t_grid[-1]], y0,
                       method="Radau", t_eval=t_grid,
                       rtol=1e-13, atol=1e-15)
    if not (sol_lo.success and sol_hi.success):
        raise RuntimeError("Radau reference solve failed")

    Y_lo, Y_hi = sol_lo.y, sol_hi.y
    diff = np.abs(Y_lo - Y_hi)
    max_diff = diff.max()

    # Confirmed digits of y(40): how many leading digits agree between the runs
    y40 = Y_hi[:, -1]
    y40_digits = []
    for i in range(3):
        d = diff[i, -1]
        digits = 0 if d <= 0 else max(0, int(-np.floor(np.log10(d))))
        y40_digits.append(digits)

    return dict(t=t_grid, Y_lo=Y_lo, Y_hi=Y_hi, diff=diff,
                max_diff=max_diff, y40=y40, y40_digits=y40_digits)


def cross_check_methods(t_grid=None):
    """Level 2: cross-family validation (LSODA / BDF vs the tight Radau run)."""
    if t_grid is None:
        t_grid = protocol_grid()
    y0 = [1.0, 0.0, 0.0]
    ref = build_reference(t_grid)

    results = {}
    for method in ("LSODA", "BDF"):
        sol = solve_ivp(robertson_rhs, [t_grid[0], t_grid[-1]], y0,
                        method=method, t_eval=t_grid,
                        rtol=1e-12, atol=1e-14)
        if not sol.success:
            raise RuntimeError(f"{method} cross-check solve failed")
        results[method] = np.abs(sol.y - ref["Y_hi"]).max()
    return ref, results


def reference_report():
    """Print the full reference-chain report (numbers for the report are copied
    from this output)."""
    print("=" * 68)
    print("reference.py - self-built reference solution, three-level chain")
    print("=" * 68)

    ref, cross = cross_check_methods()
    t, Y = ref["t"], ref["Y_hi"]

    print("[Level 1] Radau dual-tolerance comparison")
    print("  tolerance sets: rtol=1e-11/atol=1e-13  vs  rtol=1e-13/atol=1e-15 "
          "(100x apart)")
    print(f"  max pointwise component difference over grid = {ref['max_diff']:.3e}")
    print(f"  confirmed digits of y(40) per component: {ref['y40_digits']}")
    print(f"  y(40) = [{ref['y40'][0]:.10f}, {ref['y40'][1]:.10e}, "
          f"{ref['y40'][2]:.10f}]")
    # Problem-pack orientation value (sanity check only, NOT the answer)
    pack = np.array([0.7158271, 9.1855e-6, 0.2841637])
    rel = np.abs(ref["y40"] - pack) / pack
    print(f"  vs problem-pack orientation value (0.7158271, 9.1855e-6, "
          f"0.2841637), relative: {rel[0]:.2e}, {rel[1]:.2e}, {rel[2]:.2e}")

    print("[Level 2] cross-family validation (max pointwise difference vs Radau)")
    for m, d in cross.items():
        print(f"  {m:6s}: {d:.3e}")

    print("[Level 3] sanity checks")
    # Y is (3, n) row-format: sum along the component axis (axis=0)
    cons = np.abs(Y.sum(axis=0) - 1.0).max()
    print(f"  reference mass-conservation defect max|sum-1| = {cons:.3e} "
          "(expected at round-off level)")
    print(f"  reference minimum component min(Y) = {Y.min():.3e} (expected non-negative)")
    print("=" * 68)
    return ref


if __name__ == "__main__":
    reference_report()
