"""Model checks use algebraic identities and independent finite differences."""

import numpy as np
import pytest

from robertson.model import robertson_jacobian, robertson_rhs


def test_initial_derivative():
    np.testing.assert_array_equal(
        robertson_rhs(0.0, np.array([1.0, 0.0, 0.0])), [-0.04, 0.04, 0.0]
    )


@pytest.mark.parametrize(
    "y", [[1.0, 0.0, 0.0], [0.8, 2e-5, 0.19998], [0.2, 1e-5, 0.79999]]
)
def test_rhs_and_jacobian_respect_linear_invariant(y):
    state = np.array(y)
    assert abs(np.sum(robertson_rhs(0.2, state))) < 1e-14
    np.testing.assert_allclose(
        np.ones(3) @ robertson_jacobian(0.2, state), np.zeros(3), atol=1e-12
    )


@pytest.mark.parametrize("y", [[0.8, 2e-5, 0.19998], [0.2, 1e-5, 0.79999]])
def test_analytic_jacobian_against_central_differences(y):
    state = np.array(y)
    perturbations = np.array([1e-6, 1e-8, 1e-6])
    columns = []
    for column, epsilon in enumerate(perturbations):
        displacement = np.eye(3)[column] * epsilon
        columns.append(
            (robertson_rhs(0.3, state + displacement)
             - robertson_rhs(0.3, state - displacement)) / (2 * epsilon)
        )
    np.testing.assert_allclose(
        robertson_jacobian(0.3, state), np.column_stack(columns), rtol=2e-7, atol=2e-9
    )


def test_structural_zero_and_conservation_manifold_modes():
    # On y3=1-y1-y2, chain rule gives J_reduced = J[:2,:2] - J[:2,2,None].
    state = np.array([0.8, 2e-5, 0.19998])
    full = robertson_jacobian(1.0, state)
    reduced = full[:2, :2] - full[:2, 2, None]
    eigenvalues = np.linalg.eigvals(full)
    zero_index = np.argmin(np.abs(eigenvalues))
    assert abs(eigenvalues[zero_index]) < 1e-10
    np.testing.assert_allclose(
        np.sort_complex(np.delete(eigenvalues, zero_index)),
        np.sort_complex(np.linalg.eigvals(reduced)),
        rtol=1e-9,
        atol=1e-9,
    )
