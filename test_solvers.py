import unittest
import numpy as np

from mixed_bvp.core import assemble_core
from mixed_bvp.solvers import (solve_direct, solve_jacobi, solve_gauss_seidel,
                               solve_red_black_sor)
from mixed_bvp.metrics import error_metrics
from mixed_bvp.experiments.e2_solvers import (_iteration_matrix,
    _accuracy_target_sweep, _right_midpoint_qoi)


class SolverTests(unittest.TestCase):
    def setUp(self):
        self.system = assemble_core(6, "neumann", "ghost")
        self.direct = solve_direct(self.system.A, self.system.b)

    def test_jacobi_matches_direct_and_stops_on_residual(self):
        result = solve_jacobi(self.system.A, self.system.b, tol=1e-10,
                              maxiter=100000, row_kind=self.system.row_kind)
        self.assertTrue(result.converged)
        self.assertLessEqual(result.residuals["scaled_relative_inf"], 1.01e-10)
        np.testing.assert_allclose(result.x, self.direct.x, atol=2e-8, rtol=0)

    def test_red_black_gauss_seidel_and_sor_match_direct(self):
        for omega in (1.0, 1.4):
            result = solve_red_black_sor(self.system.A, self.system.b,
                                         self.system.colours, omega=omega,
                                         tol=1e-10, maxiter=100000,
                                         row_kind=self.system.row_kind)
            self.assertTrue(result.converged)
            self.assertLessEqual(result.residuals["scaled_relative_inf"], 1.01e-10)
            np.testing.assert_allclose(result.x, self.direct.x, atol=2e-8, rtol=0)

    def test_small_iteration_matrices_are_convergent(self):
        for method in ("jacobi", "gauss_seidel", "sor"):
            G = _iteration_matrix(self.system.A, self.system.colours, method, 1.4)
            self.assertLess(np.max(np.abs(np.linalg.eigvals(G))), 1.0)

    def test_invalid_red_black_colouring_is_rejected(self):
        with self.assertRaises(ValueError):
            solve_gauss_seidel(self.system.A, self.system.b,
                               np.zeros(self.system.A.shape[0], dtype=int))

    def test_tolerance_sweep_tracks_errors_and_midpoint_qoi(self):
        x, y = self.system.grid.unknown_coordinates()
        direct_error = error_metrics(
            self.direct.x, self.system.data.exact(x, y), self.system.grid)
        qoi, qoi_exact = _right_midpoint_qoi(self.system, self.direct.x)
        k = self.system.grid.k_of[self.system.grid.N, self.system.grid.N // 2]
        self.assertEqual(qoi, self.direct.x[k])
        self.assertAlmostEqual(qoi_exact, np.exp(1.5))

        rows = _accuracy_target_sweep(
            self.system, self.direct, direct_error, omega=1.2,
            tolerances=(1e-8, 1e-10), target=0.01)
        self.assertEqual(len(rows), 6)
        for method in ("jacobi", "gauss_seidel_rb", "sor_rb"):
            group = [row for row in rows if row["method"] == method]
            self.assertFalse(group[0]["stagnation_check_available"])
            self.assertTrue(group[1]["stagnation_check_available"])
            self.assertTrue(np.isfinite(group[1]["max_normalised_error_change"]))
            self.assertTrue(np.isfinite(group[1]["normalised_qoi_change"]))
            self.assertEqual(
                group[1]["meets_reporting_target"],
                group[1]["meets_forward_error_target"]
                and group[1]["meets_stagnation_target"])


if __name__ == "__main__":
    unittest.main()
