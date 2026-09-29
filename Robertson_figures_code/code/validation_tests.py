"""validation_tests.py — PASS/FAIL unit tests for the Robertson code base.

Run:  python validation_tests.py
Exit code 0 iff every test passes.
"""
import numpy as np

from model import rhs, jac, jac_reduced, prothero_robinson, Y0, T_FINAL
from methods import euler_explicit, rk4, newton_solve, euler_implicit
from adaptive import adaptive_ie
from scipy.integrate import solve_ivp

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}")


def main():
    # T1: invariant — vector field orthogonal to (1,1,1) everywhere
    rng = np.random.default_rng(0)
    for _ in range(50):
        y = rng.random(3)
        check_T1 = abs(np.sum(rhs(0.0, y))) < 1e-12
        if not check_T1:
            break
    check("T1 vector field conserves y1+y2+y3", check_T1)

    # T2: Jacobian matches complex-step differentiation of rhs
    y = np.array([0.6, 1e-5, 0.4]); J = jac(0.0, y)
    Jcs = np.column_stack([
        np.imag(rhs(0.0, y + 1e-30j * np.eye(3)[:, i])) / 1e-30
        for i in range(3)])
    check("T2 analytic Jacobian vs complex-step derivative",
          np.max(np.abs(J - Jcs)) < 1e-6 * max(1.0, np.max(np.abs(J))),
          f"max diff {np.max(np.abs(J - Jcs)):.2e}")

    # T3: Newton solver on a known nonlinear system (x^2+y^2=1, x+y=0.5)
    F = lambda x: np.array([x[0]**2 + x[1]**2 - 1.0, x[0] + x[1] - 0.5])
    JF = lambda x: np.array([[2 * x[0], 2 * x[1]], [1.0, 1.0]])
    x, it = newton_solve(F, JF, np.array([0.9, -0.3]))
    check("T3 Newton solves 2x2 nonlinear system", np.linalg.norm(F(x)) < 1e-12,
          f"residual {np.linalg.norm(F(x)):.2e}, {it} iterations")

    # T4: order verification on Prothero-Robinson (nonstiff)
    f, jf, ex = prothero_robinson(-1.0)
    errs = {}
    for h in (0.1, 0.05):
        errs[("EE", h)] = abs(euler_explicit(f, (0, 1), [1.0], h)[1][-1, 0] - ex(1.0))
        errs[("RK4", h)] = abs(rk4(f, (0, 1), [1.0], h)[1][-1, 0] - ex(1.0))
        errs[("IE", h)] = abs(euler_implicit(f, jf, (0, 1), [1.0], h)[1][-1, 0] - ex(1.0))
    oE = np.log(errs[("EE", 0.1)] / errs[("EE", 0.05)]) / np.log(2)
    oR = np.log(errs[("RK4", 0.1)] / errs[("RK4", 0.05)]) / np.log(2)
    oI = np.log(errs[("IE", 0.1)] / errs[("IE", 0.05)]) / np.log(2)
    check("T4 observed orders on PR test (EE~1, IE~1, RK4~4)",
          abs(oE - 1) < 0.1 and abs(oI - 1) < 0.1 and abs(oR - 4) < 0.2,
          f"orders {oE:.2f}/{oI:.2f}/{oR:.2f}")

    # T5: invariant preservation by all three methods (h=1e-4)
    for name, fn in (("EE", lambda: euler_explicit(rhs, (0, 1), Y0, 1e-4)),
                     ("RK4", lambda: rk4(rhs, (0, 1), Y0, 1e-4)),
                     ("IE", lambda: euler_implicit(rhs, jac, (0, 1), Y0, 1e-4)[:2])):
        out = fn()
        ys = out[1]
        check(f"T5 invariant defect ({name}) < 1e-10",
              np.max(np.abs(ys.sum(axis=1) - 1.0)) < 1e-10,
              f"max {np.max(np.abs(ys.sum(axis=1) - 1.0)):.2e}")

    # T6: short-run agreement with tight-tolerance Radau oracle
    ref = solve_ivp(rhs, (0, 1), Y0, method="Radau", jac=jac,
                    rtol=1e-12, atol=1e-14, t_eval=[1.0]).y[:, -1]
    _, ys, _ = euler_implicit(rhs, jac, (0, 1), Y0, 1e-4)
    check("T6 IE (h=1e-4) agrees with Radau oracle at t=1",
          np.max(np.abs(ys[-1] - ref)) < 1e-6,
          f"max diff {np.max(np.abs(ys[-1] - ref)):.2e}")

    # T7: stiffness-ratio check values from the problem statement
    ref40 = solve_ivp(rhs, (0, T_FINAL), Y0, method="Radau", jac=jac,
                      rtol=1e-12, atol=1e-14, t_eval=[40.0]).y[:, -1]
    S40 = None
    e = np.linalg.eigvals(jac_reduced(ref40[0], ref40[1]))
    mags = np.sort(np.abs(np.real(e)))
    S40 = mags[-1] / mags[0]
    check("T7 stiffness ratio S(40) ~ 1.58e5 and |lmax|(40) ~ 3.39e3",
          abs(S40 - 1.584e5) / 1.584e5 < 0.01 and abs(mags[-1] - 3393) / 3393 < 0.01,
          f"S(40)={S40:.4e}, |lmax|={mags[-1]:.2f}")

    # T8: adaptive controller reduces error under tolerance tightening
    _, ya1, *_ = adaptive_ie(rhs, jac, (0, T_FINAL), Y0, rtol=1e-3)
    _, ya2, *_ = adaptive_ie(rhs, jac, (0, T_FINAL), Y0, rtol=1e-5)
    refe = solve_ivp(rhs, (0, T_FINAL), Y0, method="Radau", jac=jac,
                     rtol=1e-12, atol=1e-14, t_eval=[T_FINAL]).y[:, -1]
    wv = 1e-8 + 1e-6 * np.abs(refe)
    e1 = np.max(np.abs(ya1[-1] - refe) / wv); e2 = np.max(np.abs(ya2[-1] - refe) / wv)
    check("T8 adaptive IE: tighter rtol => smaller endpoint error",
          e2 < e1 / 5, f"err {e1:.1e} -> {e2:.1e}")

    n_fail = sum(1 for _, ok in RESULTS if not ok)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} tests passed")
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
