import unittest
import numpy as np

from mixed_bvp.core import assemble_core, matrix_audit
from mixed_bvp.data import CORE_DATA, square_exact, square_forcing, sine_exact, sine_forcing
from mixed_bvp.metrics import error_metrics
from mixed_bvp.solvers import solve_direct
from mixed_bvp.verification import local_row_audit, dirichlet_corner_audit


class CoreAssemblyTests(unittest.TestCase):
    def test_manufactured_forcing_identities(self):
        x = np.array([0.13, 0.41, 0.82])
        y = np.array([0.27, 0.63, 0.09])
        np.testing.assert_allclose(square_forcing(x, y), -2 * square_exact(x, y))
        np.testing.assert_allclose(sine_forcing(x, y), 2 * np.pi**2 * sine_exact(x, y))
        np.testing.assert_allclose(CORE_DATA.robin_right(y), 3 * np.exp(1 + y))

    def test_index_map_and_corner_protection(self):
        for variant in ("neumann", "robin"):
            for closure in ("one_sided", "ghost"):
                system = assemble_core(6, variant, closure)
                audit = dirichlet_corner_audit(system)
                self.assertTrue(all(audit.values()))
                self.assertEqual(system.A.shape, (30, 30))

    def test_independent_local_residual_formulas(self):
        for variant in ("neumann", "robin"):
            for closure in ("one_sided", "ghost"):
                system = assemble_core(5, variant, closure)
                audit = local_row_audit(system)
                self.assertLess(audit["max_independent_formula_difference"], 5e-10)

    def test_sign_pattern_connectivity_and_nonsymmetry(self):
        for variant in ("neumann", "robin"):
            for closure in ("one_sided", "ghost"):
                audit = matrix_audit(assemble_core(5, variant, closure))
                self.assertGreater(audit["min_diagonal"], 0)
                self.assertLessEqual(audit["max_off_diagonal"], 1e-14)
                self.assertEqual(audit["connected_components"], 1)
                self.assertGreater(audit["nonsymmetry"], 1e-6)
                self.assertGreaterEqual(audit["min_row_sum"], -1e-12)
                self.assertEqual(audit["right_corner_unknowns"], 0)

    def test_dirichlet_rhs_sign(self):
        N = 6
        system = assemble_core(N, "neumann", "ghost")
        h = system.grid.h
        k = system.grid.k_of[1, 1]
        expected = (CORE_DATA.forcing(h, h)
                    + (CORE_DATA.dirichlet(0, h) + CORE_DATA.dirichlet(h, 0)) / h**2)
        self.assertAlmostEqual(system.b[k], expected, places=11)

    def test_ghost_robin_row_coefficients(self):
        N = 6
        system = assemble_core(N, "robin", "ghost")
        h, beta = system.grid.h, CORE_DATA.beta
        j = 3
        k = system.grid.k_of[N, j]
        self.assertAlmostEqual(system.A[k, k], (4 + 2 * h * beta) / h**2)
        self.assertAlmostEqual(system.A[k, system.grid.k_of[N - 1, j]], -2 / h**2)
        self.assertAlmostEqual(system.A[k, system.grid.k_of[N, j - 1]], -1 / h**2)
        self.assertAlmostEqual(system.A[k, system.grid.k_of[N, j + 1]], -1 / h**2)
        expected_b = CORE_DATA.forcing(1, j * h) + 2 * CORE_DATA.robin_right(j * h) / h
        self.assertAlmostEqual(system.b[k], expected_b)

    def test_one_sided_row_replaces_pde(self):
        N = 6
        system = assemble_core(N, "neumann", "one_sided")
        h = system.grid.h
        k = system.grid.k_of[N, 3]
        self.assertEqual(system.A[k].nnz, 2)
        self.assertAlmostEqual(system.A[k, k], 1 / h)
        self.assertAlmostEqual(system.b[k], CORE_DATA.neumann_right(3 * h))

    def test_representative_convergence(self):
        for closure in ("one_sided", "ghost"):
            errors = []
            for N in (8, 16):
                system = assemble_core(N, "neumann", closure)
                direct = solve_direct(system.A, system.b)
                x, y = system.grid.unknown_coordinates()
                metric = "Linf_right" if closure == "one_sided" else "L2h"
                errors.append(error_metrics(direct.x, CORE_DATA.exact(x, y), system.grid)[metric])
            ratio = errors[0] / errors[1]
            self.assertGreater(ratio, 1.6 if closure == "one_sided" else 3.0)
            self.assertLess(ratio, 2.6 if closure == "one_sided" else 5.0)


if __name__ == "__main__":
    unittest.main()
