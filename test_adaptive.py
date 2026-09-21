"""Adaptive checks include rejected work, endpoint, and actual accuracy."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from robertson.adaptive import solve_adaptive
from robertson.model import robertson_jacobian, robertson_rhs
import robertson.solvers as solver_module


@pytest.mark.parametrize("method", ["RK4", "IE"])
def test_adaptive_rejects_large_trial_and_reaches_endpoint(method):
    t_end = 0.003
    tol = 1e-10
    y0 = np.array([1.0, 0.0, 0.0])
    t, Y, stats = solve_adaptive(
        method, 0.0, t_end, y0, t_end, tol,
        {"newton_tol": 1e-14},
    )
    assert stats["converged"]
    assert stats["failure"] is None
    assert t[-1] == pytest.approx(t_end, abs=1e-15)
    assert np.all(np.diff(t) > 0)
    assert stats["accepted_steps"] == len(t) - 1
    assert stats["rejected_steps"] > 0
    assert len(stats["h_history"]) == stats["accepted_steps"]
    np.testing.assert_allclose(stats["h_history"], np.diff(t), rtol=1e-10, atol=1e-15)
    assert np.max(stats["h_history"]) > 1.01 * np.min(stats["h_history"])
    if method == "RK4":
        # Three RK4 advances per trial; rejected trials must remain in total cost.
        assert stats["rhs_evals"] == 12 * (
            stats["accepted_steps"] + stats["rejected_steps"]
        )
    else:
        assert stats["newton_iters"] > stats["accepted_steps"]
        assert stats["linear_solves"] == stats["newton_iters"]
    oracle = solve_ivp(
        robertson_rhs, (0.0, t_end), y0, method="Radau",
        jac=robertson_jacobian, rtol=1e-12, atol=1e-15, t_eval=t,
    )
    assert oracle.success
    error = np.max(np.linalg.norm(Y - oracle.y.T, axis=1))
    # Local tolerance is not a global-error promise. Check independent accuracy.
    assert error < (1e-8 if method == "RK4" else 1e-7)
    assert np.max(np.abs(Y.sum(axis=1) - 1.0)) < 1e-11


def test_adaptive_counter_includes_instrumented_rejected_work(monkeypatch):
    calls = []
    original = solver_module.robertson_rhs

    def counted_rhs(t, y):
        calls.append(t)
        return original(t, y)

    monkeypatch.setattr(solver_module, "robertson_rhs", counted_rhs)
    _, _, stats = solve_adaptive(
        "RK4", 0.0, 0.003, [1.0, 0.0, 0.0], 0.003, 1e-12
    )
    assert stats["converged"]
    assert stats["rejected_steps"] > 0
    assert stats["rhs_evals"] == len(calls)
    assert len(calls) > 12 * stats["accepted_steps"]
