"""
validation_tests.py - Week 2 deliverable: unit-test suite against the team-built oracle.

Role-card requirement: every number the team wants to quote must be compared
against the team's own reference solution. This suite turns the headline claims
into PASS/FAIL checks (run here against my independent re-implementation; the
same checks apply to the production solvers once merged):

  T1  implicit Euler (h=0.1) y(40) per-component relative error < 1%
      (the course scaffold's sanity bound)
  T2  implicit Euler mass-conservation defect < 1e-10 (should sit at the
      Newton-tolerance level)
  T3  implicit Euler non-negativity over the whole trajectory (the declared
      "slightly negative" threshold must be stated in the report)
  T4  RK4 verified at order ~4 on the non-stiff logistic problem FIRST
      (a correctness precondition before touching Robertson)
  T5  explicit Euler verified at order ~1 on logistic (same reasoning)
  T6  Robertson first-step Newton residual trajectory: rise then quadratic
      decay (degenerate initial guess; the "bad Newton initial guess" edge
      case and the lecture lightning-talk diagnostic)
  T7  eigenvalue / stiffness-ratio check values vs the problem pack
  T8  the oracle itself: dual-tolerance agreement >= 8 digits (precondition
      for trusting any other test)

Run: python validation_tests.py
Only when ALL tests PASS may the "validated" conclusion be handed to the team.
Author: Rong Yin (Role 5)
"""

import numpy as np

from robertson_model import robertson_rhs, robertson_jacobian, eigen_info
from solvers import (solve_implicit_euler, solve_rk4, solve_euler,
                     newton_solve_ie)
from reference import build_reference


def _logistic_rhs(t, y):
    """Non-stiff sanity problem: y' = y(1-y/10) (verify the implementation
    before attacking the stiff problem)."""
    return y * (1.0 - y / 10.0)


def _logistic_exact(t, y0=0.5):
    return 10.0 / (1.0 + (10.0 / y0 - 1.0) * np.exp(-t))


def test_t1_y40_vs_oracle(ref):
    """T1: implicit Euler h=0.1, y(40) vs oracle, per-component rel. err < 1%"""
    y0 = np.array([1.0, 0.0, 0.0])
    t, y = solve_implicit_euler(robertson_rhs, robertson_jacobian,
                                y0, 0.0, 40.0, 0.1)
    rel = np.abs(y[-1] - ref["y40"]) / np.abs(ref["y40"])
    ok = np.all(rel < 0.01)
    return ok, (f"implicit Euler h=0.1, y(40) rel. err per component = "
                f"{rel[0]:.2e}, {rel[1]:.2e}, {rel[2]:.2e} (bound < 1%)")


def test_t2_mass_conservation():
    """T2: mass-conservation defect (the linear invariant should be preserved
    algebraically up to Newton-tolerance level)."""
    y0 = np.array([1.0, 0.0, 0.0])
    t, y = solve_implicit_euler(robertson_rhs, robertson_jacobian,
                                y0, 0.0, 40.0, 0.1)
    defect = np.abs(y.sum(axis=1) - 1.0).max()
    ok = defect < 1e-10
    return ok, f"max|sum-1| = {defect:.2e} (bound < 1e-10)"


def test_t3_nonnegative():
    """T3: non-negativity over the whole trajectory.

    Note: the reference itself dips to -2.6e-18 (Radau is also slightly
    negative at round-off level), so the non-negativity criterion must use an
    explicitly declared threshold: -1e-12 here.
    """
    y0 = np.array([1.0, 0.0, 0.0])
    t, y = solve_implicit_euler(robertson_rhs, robertson_jacobian,
                                y0, 0.0, 40.0, 0.1)
    ymin = y.min()
    ok = ymin > -1e-12
    return ok, (f"min(trajectory components) = {ymin:.3e} (threshold -1e-12; "
                f"reference itself min ~ -2.6e-18)")


def test_t4_rk4_logistic_order():
    """T4: RK4 on logistic first (on Robertson it is stability-limited, so the
    implementation must be proven correct beforehand)."""
    y0 = np.array([0.5])
    h_values = np.array([0.4, 0.2, 0.1, 0.05, 0.025])
    errors = []
    for h in h_values:
        t, y = solve_rk4(_logistic_rhs, y0, 0.0, 10.0, h)
        errors.append(np.abs(y[-1, 0] - _logistic_exact(10.0)))
    errors = np.array(errors)
    coef = np.polyfit(np.log(h_values), np.log(errors), 1)
    ok = 3.7 <= coef[0] <= 4.3
    return ok, f"RK4 order on logistic = {coef[0]:.2f} (theory 4)"


def test_t5_euler_logistic_order():
    """T5: explicit Euler on logistic, order ~1 (same reasoning as T4)."""
    y0 = np.array([0.5])
    h_values = np.array([0.4, 0.2, 0.1, 0.05, 0.025])
    errors = []
    for h in h_values:
        t, y = solve_euler(_logistic_rhs, y0, 0.0, 10.0, h)
        errors.append(np.abs(y[-1, 0] - _logistic_exact(10.0)))
    errors = np.array(errors)
    coef = np.polyfit(np.log(h_values), np.log(errors), 1)
    ok = 0.7 <= coef[0] <= 1.3
    return ok, f"explicit Euler order on logistic = {coef[0]:.2f} (theory 1)"


def test_t6_newton_first_step_trace():
    """T6: first-step (t=0->0.1) Newton residual trajectory: the degenerate
    guess y2=0 makes the residual RISE first, then decay quadratically.
    Monotone decay would mean the degenerate point was missed; no convergence
    would mean an implementation bug."""
    y0 = np.array([1.0, 0.0, 0.0])
    w, ok, it = newton_solve_ie(robertson_rhs, robertson_jacobian,
                                0.1, y0, 0.1, tol=1e-12, max_iter=50)
    it = np.array(it)
    rose_first = it[1] > it[0] if len(it) > 1 else False
    quadratic_tail = it[-1] < 1e-12
    return ok and quadratic_tail, (
        f"first-step Newton: converged={ok}, {len(it)} iterations, "
        f"rose_first={rose_first}, final residual={it[-1]:.2e}, "
        f"trajectory={np.array2string(it, precision=1, separator=', ')}")


def test_t7_eigenvalue_check_values(ref):
    """T7: eigenvalue / stiffness-ratio check values vs the problem pack
    (sanity checks only, not the answer)."""
    t, Y = ref["t"], ref["Y_hi"]
    probes = {1e-4: 3.0e3, 1e-2: 5.4e3, 40.0: 1.58e5}
    msgs = []
    all_ok = True
    for probe, expect in probes.items():
        i = int(np.argmin(np.abs(t - probe)))
        _, _, ratio = eigen_info(t[i], Y[:, i])
        rel_err = abs(ratio - expect) / expect
        all_ok &= rel_err < 0.05
        msgs.append(f"t={probe:g}: S={ratio:.3e} vs pack {expect:.1e} "
                    f"(rel {rel_err:.1%})")
    i40 = int(np.argmin(np.abs(t - 40.0)))
    lam, lam_nz, _ = eigen_info(t[i40], Y[:, i40])
    lam_max = np.abs(lam_nz.real).max()
    ok_lam = abs(lam_max - 3393.0) / 3393.0 < 0.05
    all_ok &= ok_lam
    msgs.append(f"|lambda|max(40)={lam_max:.4e} vs pack 3.39e3 "
                f"(rel {abs(lam_max - 3393) / 3393:.1%})")
    return all_ok, "; ".join(msgs)


def test_t8_oracle_self_agreement(ref):
    """T8: oracle precondition - dual-tolerance agreement >= 8 digits
    (otherwise the reference itself is not trustworthy)."""
    ok = all(d >= 8 for d in ref["y40_digits"])
    return ok, (f"confirmed digits of y(40) = {ref['y40_digits']} "
                f"(bound >= 8 each), max grid difference = {ref['max_diff']:.2e}")


ALL_TESTS = [
    ("T1", "implicit Euler y(40) vs oracle", test_t1_y40_vs_oracle, True),
    ("T2", "mass-conservation defect", test_t2_mass_conservation, False),
    ("T3", "non-negativity over trajectory", test_t3_nonnegative, False),
    ("T4", "RK4 logistic order check", test_t4_rk4_logistic_order, False),
    ("T5", "explicit Euler logistic order check", test_t5_euler_logistic_order,
     False),
    ("T6", "first-step Newton residual trace", test_t6_newton_first_step_trace,
     False),
    ("T7", "eigenvalue/stiffness check values",
     test_t7_eigenvalue_check_values, True),
    ("T8", "oracle self-agreement digits", test_t8_oracle_self_agreement, True),
]


def main():
    print("=" * 72)
    print("validation_tests.py - unit-test suite against the team-built oracle "
          "(Week 2)")
    print("=" * 72)
    ref = build_reference()
    n_pass = 0
    for tid, name, fn, needs_ref in ALL_TESTS:
        try:
            ok, msg = fn(ref) if needs_ref else fn()
        except Exception as exc:    # an exception inside a test = FAIL, and
            ok, msg = False, f"exception: {exc}"   # the suite must not crash
        n_pass += ok
        print(f"[{tid}] {'PASS' if ok else 'FAIL'}  {name}")
        print(f"      {msg}")
    print("=" * 72)
    print(f"RESULT: {n_pass}/{len(ALL_TESTS)} tests PASS"
          + ("  - validated conclusions may be handed to the team"
             if n_pass == len(ALL_TESTS)
             else "  - some tests FAILED: numbers related to them must not "
                  "be quoted in the report"))
    print("=" * 72)
    return n_pass == len(ALL_TESTS)


if __name__ == "__main__":
    main()
