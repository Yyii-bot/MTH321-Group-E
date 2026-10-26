"""Pure-Neumann compatibility, constant nullspace, and mean-gauge demonstration."""

from __future__ import annotations
import numpy as np
import scipy.sparse as sp
from skfem import Basis, FacetBasis, ElementTriP1, LinearForm, asm
from skfem.models.poisson import laplace, mass

from .solvers import solve_direct
from .data import square_forcing


def fem_pure_neumann_solve(mesh, intorder=8):
    """Solve the compatible, mean-gauged P1 pure-Neumann problem.

    The exact data come from u=exp(x+y).  The discrete compatibility defect is
    measured before correction.  A uniform forcing correction b-=defect*M*1
    enforces 1^T b=0, and a saddle-point row fixes the continuous exact mean
    (e-1)^2 rather than silently choosing an unrelated pressure gauge.
    """
    basis = Basis(mesh, ElementTriP1(), intorder=intorder)
    K = asm(laplace, basis)
    M = asm(mass, basis)

    @LinearForm
    def volume_load(v, w):
        return square_forcing(w.x[0], w.x[1]) * v

    F = asm(volume_load, basis)
    fbasis = FacetBasis(mesh, ElementTriP1(), intorder=intorder,
                        facets=mesh.boundary_facets())

    @LinearForm
    def boundary_load(v, w):
        x, y = w.x
        q = np.select(
            (np.isclose(x, 0.0), np.isclose(x, 1.0),
             np.isclose(y, 0.0), np.isclose(y, 1.0)),
            (-np.exp(x + y), np.exp(x + y),
             -np.exp(x + y), np.exp(x + y)),
        )
        return q * v

    Q = asm(boundary_load, fbasis)
    b = F + Q
    ones = np.ones(mesh.p.shape[1])
    defect = float(ones @ b)
    mean_weight = M @ ones  # 1^T(M1)=|Omega|=1 on the unit square
    b_compatible = b - defect * mean_weight
    corrected_defect = float(ones @ b_compatible)
    exact_mean = float((np.e - 1.0) ** 2)
    saddle = sp.bmat(
        [[K.tocsr(), sp.csr_matrix(mean_weight[:, None])],
         [sp.csr_matrix(mean_weight[None, :]), None]],
        format="csr",
    )
    rhs = np.concatenate((b_compatible, [exact_mean]))
    direct = solve_direct(saddle, rhs)
    values = direct.x[:-1]
    multiplier = float(direct.x[-1])
    null_vector_residual = float(np.linalg.norm(K @ ones, 2))
    matrix_frobenius = float(np.sqrt(K.power(2).sum()))
    return {
        "values": values,
        "defect_before": defect,
        "defect_after": corrected_defect,
        "null_vector_relative_residual": null_vector_residual / matrix_frobenius,
        "mean": float(mean_weight @ values),
        "exact_mean": exact_mean,
        "gauge_multiplier": multiplier,
        "saddle_residual": direct.residuals,
        "A": K,
        "b": b_compatible,
    }


def small_neumann_spectrum(mesh):
    """Dense small-mesh nullspace diagnostic; never used as a proof on fine grids."""
    basis = Basis(mesh, ElementTriP1())
    K = asm(laplace, basis).toarray()
    eigenvalues = np.linalg.eigvalsh(K)
    return {
        "smallest_absolute_eigenvalue": float(np.min(np.abs(eigenvalues))),
        "second_smallest_eigenvalue": float(np.sort(eigenvalues)[1]),
        "zero_eigenvalue_count_at_1e-10": int(np.sum(np.abs(eigenvalues) < 1e-10)),
    }
