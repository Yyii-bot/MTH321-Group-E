"""Test real Robertson trajectories, orders, and failure behavior."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from robertson.model import robertson_jacobian, robertson_rhs
from robertson.solvers import solve_fixed
import robertson.solvers as solver_module


Y0 = np.array([1.0, 0.0, 0.0])


@pytest.fixture(scope="module")
def short_oracle():
    oracle = solve_ivp(
        robertson_rhs, (0.0, 0.001), Y0, method="Radau",
        jac=robertson_jacobian, rtol=1e-12, atol=1e-16, dense_output=True,
    )
    assert oracle.success
    return oracle.sol


@pytest.mark.parametrize("method", ["EE", "RK4", "IE"])
def test_fixed_solver_matches_independent_oracle(method, short_oracle):
    t, Y, stats = solve_fixed(
        method, 0.0, 0.001, Y0, 2.5e-5, {"newton_tol": 1e-15}
    )
    assert stats["converged"]
    assert stats["failure"] is None
    assert t[-1] == pytest.approx(0.001, abs=1e-15)
    assert np.all(np.diff(t) > 0)
    assert Y.shape == (len(t), 3)
    error = np.max(np.linalg.norm(Y - short_oracle(t).T, axis=1))
    assert error < (1e-11 if method == "RK4" else 5e-7)
    assert np.max(np.abs(Y.sum(axis=1) - 1.0)) < 1e-12
    assert stats["steps"] == len(t) - 1
    assert stats["accepted_steps"] == stats["steps"]
    assert stats["rejected_steps"] == 0
    if method in ("EE", "RK4"):
        expected_calls = {"EE": 1, "RK4": 4}[method] * stats["steps"]
        assert stats["rhs_evals"] == expected_calls
        assert stats["newton_iters"] == 0
    else:
        assert stats["rhs_evals"] >= stats["newton_iters"]
        assert stats["linear_solves"] == stats["newton_iters"]


@pytest.mark.parametrize("method,order", [("EE", 1.0), ("RK4", 4.0), ("IE", 1.0)])
def test_four_successive_halvings_show_expected_order(method, order, short_oracle):
    errors = []
    for h in (2e-4, 1e-4, 5e-5, 2.5e-5):
        t, Y, stats = solve_fixed(
            method, 0.0, 0.001, Y0, h, {"newton_tol": 1e-15}
        )
        assert stats["converged"]
        errors.append(np.max(np.linalg.norm(Y - short_oracle(t).T, axis=1)))
    orders = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
    assert np.all(np.diff(errors) < 0), errors
    # Use the two finer comparisons; coarse grids need not be asymptotic.
    assert np.all(np.abs(orders[-2:] - order) < 0.35), (errors, orders)


def test_last_fixed_step_is_shortened_to_endpoint():
    t, Y, stats = solve_fixed("EE", 0.0, 0.001, Y0, 3e-4)
    assert stats["converged"]
    assert len(t) == 5
    assert t[-1] == pytest.approx(0.001, abs=1e-15)
    assert np.diff(t)[-1] == pytest.approx(1e-4)


def test_newton_failure_is_flagged_without_accepting_failed_state():
    t, Y, stats = solve_fixed(
        "IE", 0.0, 1.0, Y0, 0.1,
        {"newton_tol": 1e-15, "newton_maxiter": 1},
    )
    assert not stats["converged"]
    assert "newton" in str(stats["failure"]).lower()
    assert stats["failure_time"] == pytest.approx(0.1)
    assert stats["accepted_steps"] == 0
    assert len(t) == 1
    np.testing.assert_array_equal(Y[0], Y0)
    assert stats["rhs_evals"] > 0


def test_large_explicit_step_reports_failure_without_crashing():
    t, Y, stats = solve_fixed(
        "EE", 0.0, 40.0, Y0, 0.01, {"blowup_threshold": 100.0}
    )
    assert not stats["converged"]
    assert stats["failure"]
    assert stats["failure_time"] is not None
    assert stats["failure_time"] < 40.0
    assert t[-1] < 40.0
    assert stats["first_negative_time"] is not None
    assert np.all(np.isfinite(Y))


def test_max_steps_produces_explicit_incomplete_status():
    t, _, stats = solve_fixed("EE", 0.0, 0.001, Y0, 1e-5, {"max_steps": 2})
    assert not stats["converged"]
    assert stats["failure"]
    assert stats["accepted_steps"] == 2
    assert t[-1] < 0.001


@pytest.mark.parametrize("h", [0.0, -0.1, np.nan, np.inf])
def test_invalid_step_size_is_rejected(h):
    with pytest.raises(ValueError):
        solve_fixed("EE", 0.0, 0.001, Y0, h)


@pytest.mark.parametrize("method", ["EE", "RK4", "IE"])
def test_rhs_counter_equals_instrumented_calls(method, monkeypatch):
    calls = []
    original = solver_module.robertson_rhs

    def counted_rhs(t, y):
        calls.append(t)
        return original(t, y)

    monkeypatch.setattr(solver_module, "robertson_rhs", counted_rhs)
    _, _, stats = solve_fixed(method, 0.0, 0.001, Y0, 1e-4)
    assert stats["converged"]
    assert stats["rhs_evals"] == len(calls)
    assert len(calls) > 0
