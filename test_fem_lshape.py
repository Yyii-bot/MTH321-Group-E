import unittest
import numpy as np

from mixed_bvp.data import (harmonic_exact, singular_exact, singular_gradient)
from mixed_bvp.dirichlet import (square_grid, lshape_grid,
                                 assemble_five_point_dirichlet, grid_error_metrics)
from mixed_bvp.fem import (unit_square_triangulation, lshape_triangulation,
                           radial_graded_triangulation, fem_dirichlet_solve,
                           fem_error_integrals)
from mixed_bvp.pure_neumann import fem_pure_neumann_solve, small_neumann_spectrum
from mixed_bvp.solvers import solve_direct


class FemAndLshapeTests(unittest.TestCase):
    def test_singular_branch_and_zero_trace_rays(self):
        x = np.linspace(0.05, 0.95, 12)
        y = np.linspace(0.05, 0.95, 12)
        np.testing.assert_allclose(singular_exact(x, 0 * x), 0, atol=2e-15)
        np.testing.assert_allclose(singular_exact(0 * y, y), 0, atol=2e-15)
        with self.assertRaises(ValueError):
            singular_gradient(0.0, 0.0)

    def test_lshape_grid_has_no_unknown_next_to_inactive_node(self):
        grid = lshape_grid(8)
        for i, j in zip(*grid.unknown_ij()):
            for ii, jj in ((i-1, j), (i+1, j), (i, j-1), (i, j+1)):
                self.assertTrue(grid.active[ii, jj])

    def test_harmonic_zero_forcing_discrete_maximum_principle(self):
        grid = lshape_grid(8)
        zero = lambda x, y: np.zeros_like(np.asarray(x))
        system = assemble_five_point_dirichlet(grid, zero, harmonic_exact)
        result = solve_direct(system.A, system.b)
        boundary_values = harmonic_exact(*grid.coordinates())[grid.node_type == 0]
        self.assertGreaterEqual(result.x.min(), boundary_values.min() - 1e-12)
        self.assertLessEqual(result.x.max(), boundary_values.max() + 1e-12)

    def test_p1_reproduces_linear_function(self):
        mesh = unit_square_triangulation(4)
        exact = lambda x, y: 1 + 2 * x - 3 * y
        gradient = lambda x, y: (2 + 0 * x, -3 + 0 * x)
        result = fem_dirichlet_solve(mesh, lambda x, y: 0 * x, exact)
        errors = fem_error_integrals(mesh, result.values, exact, gradient)
        self.assertLess(errors["L2"], 1e-12)
        self.assertLess(errors["H1"], 1e-12)
        self.assertLess(errors["Linf"], 1e-12)

    def test_fdm_fem_common_node_coordinates(self):
        N = 4
        grid = square_grid(N)
        mesh = unit_square_triangulation(N)
        result = fem_dirichlet_solve(mesh, lambda x, y: -2 * np.exp(x + y),
                                     lambda x, y: np.exp(x + y))
        i, j = grid.unknown_ij()
        np.testing.assert_allclose(mesh.p[0, result.interior], grid.x[i])
        np.testing.assert_allclose(mesh.p[1, result.interior], grid.y[j])

    def test_matched_p1_stiffness_equals_scaled_five_point_stencil(self):
        N = 4
        grid = square_grid(N)
        fd = assemble_five_point_dirichlet(
            grid, lambda x, y: -2 * np.exp(x + y), lambda x, y: np.exp(x + y)
        )
        fe = fem_dirichlet_solve(
            unit_square_triangulation(N), lambda x, y: -2 * np.exp(x + y),
            lambda x, y: np.exp(x + y)
        )
        difference = fd.A - (N**2) * fe.A
        relative = np.sqrt(difference.power(2).sum() / fd.A.power(2).sum())
        self.assertLess(relative, 1e-12)
        self.assertEqual(fd.A.nnz, fe.A.nnz)

    def test_corner_polar_error_integration_is_finite_and_accurate_l2(self):
        mesh = lshape_triangulation(8)
        result = fem_dirichlet_solve(mesh, lambda x, y: 0 * x, singular_exact)
        base = fem_error_integrals(mesh, result.values, singular_exact,
                                   singular_gradient, singular_corner=False)
        polar = fem_error_integrals(mesh, result.values, singular_exact,
                                    singular_gradient, singular_corner=True)
        self.assertTrue(np.isfinite(polar["H1"]))
        self.assertGreater(polar["H1"], 0)
        self.assertLess(abs(polar["L2"] - base["L2"]) / polar["L2"], 0.01)
        self.assertEqual(polar["quadrature_corner_elements"], 4)

    def test_radial_mesh_contains_outer_polygon_corners(self):
        mesh = radial_graded_triangulation(4, 2.0)
        for point in ((-1, 1), (-1, -1), (1, -1), (1, 0), (0, 1)):
            self.assertTrue(np.any(np.linalg.norm(mesh.p - np.array(point)[:, None], axis=0) < 1e-12))

    def test_pure_neumann_nullspace_compatibility_and_gauge(self):
        mesh = unit_square_triangulation(4)
        result = fem_pure_neumann_solve(mesh, intorder=2)
        self.assertLess(abs(result["defect_after"]), 1e-12)
        self.assertLess(result["null_vector_relative_residual"], 1e-12)
        self.assertAlmostEqual(result["mean"], result["exact_mean"], places=12)
        spectrum = small_neumann_spectrum(mesh)
        self.assertEqual(spectrum["zero_eigenvalue_count_at_1e-10"], 1)


if __name__ == "__main__":
    unittest.main()
