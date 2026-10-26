"""Two-dimensional P1 finite elements, mesh construction, and error integration."""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import scipy.sparse as sp
from skfem import MeshTri, Basis, ElementTriP1, LinearForm, asm
from skfem.models.poisson import laplace
from skfem.quadrature import get_quadrature
from skfem.refdom import RefTri

from .solvers import solve_direct


def _mesh_from_points_triangles(points, triangles):
    points = np.asarray(points, dtype=float).T
    triangles = np.asarray(triangles, dtype=np.int32)
    # Ensure positive orientation and retain only referenced vertices.
    p0 = points[:, triangles[:, 0]]
    p1 = points[:, triangles[:, 1]]
    p2 = points[:, triangles[:, 2]]
    det = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])
    triangles[det < 0.0] = triangles[det < 0.0][:, [0, 2, 1]]
    used, inverse = np.unique(triangles.ravel(), return_inverse=True)
    points = points[:, used]
    triangles = inverse.reshape((-1, 3)).astype(np.int32)
    return MeshTri(points, triangles.T)


def unit_square_triangulation(N):
    """N by N squares, all split along the (i,j)-(i+1,j+1) diagonal."""
    points = [(i / N, j / N) for j in range(N + 1) for i in range(N + 1)]
    triangles = []
    for j in range(N):
        for i in range(N):
            v00 = j * (N + 1) + i
            v10 = v00 + 1
            v01 = v00 + N + 1
            v11 = v01 + 1
            triangles.extend(((v00, v10, v11), (v00, v11, v01)))
    return _mesh_from_points_triangles(points, triangles)


def lshape_triangulation(N):
    """Conforming aligned P1 mesh of the standard three-quadrant L polygon."""
    if N % 2 or N < 4:
        raise ValueError("L-shaped triangulation requires even N>=4")
    points = [(-1.0 + 2.0 * i / N, -1.0 + 2.0 * j / N)
              for j in range(N + 1) for i in range(N + 1)]
    triangles = []
    mid = N // 2
    for j in range(N):
        for i in range(N):
            if i >= mid and j >= mid:
                continue
            v00 = j * (N + 1) + i
            v10 = v00 + 1
            v01 = v00 + N + 1
            v11 = v01 + 1
            triangles.extend(((v00, v10, v11), (v00, v11, v01)))
    return _mesh_from_points_triangles(points, triangles)


def corner_graded_triangulation(N, kappa=3, delta=0.5):
    """Red refinement with nested distance bands and library-managed closure.

    At pass ell, elements with centroid radius < delta/2^ell are marked.
    Consequently h near radius r scales approximately as h0*r/delta, while
    scikit-fem's adaptive refinement supplies the conforming closure.
    """
    if kappa < 0:
        raise ValueError("kappa must be nonnegative")
    mesh = lshape_triangulation(N)
    for level in range(kappa):
        centres = mesh.p[:, mesh.t].mean(axis=1)
        radius = np.hypot(centres[0], centres[1])
        marked = np.flatnonzero(radius < delta / (2.0**level))
        if len(marked) == 0:
            break
        mesh = mesh.refined(marked)
    return mesh


def radial_graded_triangulation(n=16, mu=1.0):
    """Polar radial grading on the exact L-shaped outer polygon.

    The radial law is r_i=R(theta)*(i/n)^mu.  Angular rays include all three
    outer square corners, so the polygon is not changed to a circle.
    """
    if n < 2 or mu < 1.0:
        raise ValueError("Use n>=2 and mu>=1")
    ntheta = 6 * n  # includes pi/4 increments and all outer polygon corners
    theta = np.linspace(np.pi / 2.0, 2.0 * np.pi, ntheta + 1)
    corners = np.array([3.0 * np.pi / 4.0, 5.0 * np.pi / 4.0, 7.0 * np.pi / 4.0])
    theta = np.unique(np.concatenate((theta, corners)))
    c, s = np.cos(theta), np.sin(theta)
    # Distance from the re-entrant corner to the square's outer sides.
    candidates = []
    for value, direction in ((1.0, c), (-1.0, c), (1.0, s), (-1.0, s)):
        with np.errstate(divide="ignore", invalid="ignore"):
            distance = value / direction
        distance[(direction == 0.0) | (distance <= 0.0)] = np.inf
        candidates.append(distance)
    R = np.min(np.vstack(candidates), axis=0)

    points = [(0.0, 0.0)]
    # ring index starts at one; node indices are retrieved from a table.
    ring = np.zeros((n, len(theta)), dtype=np.int32)
    for i in range(1, n + 1):
        radius = R * (i / n) ** mu
        for j in range(len(theta)):
            ring[i - 1, j] = len(points)
            points.append((radius[j] * c[j], radius[j] * s[j]))

    triangles = []
    for j in range(len(theta) - 1):
        triangles.append((0, ring[0, j], ring[0, j + 1]))
        for i in range(n - 1):
            inner_j = ring[i, j]
            outer_j = ring[i + 1, j]
            inner_next = ring[i, j + 1]
            outer_next = ring[i + 1, j + 1]
            triangles.extend(((inner_j, outer_j, outer_next),
                              (inner_j, outer_next, inner_next)))
    return _mesh_from_points_triangles(points, triangles)


@dataclass
class FEMResult:
    mesh: MeshTri
    values: np.ndarray
    interior: np.ndarray
    boundary: np.ndarray
    A: sp.csr_matrix
    b: np.ndarray
    assembly_time: float
    solve_time: float
    residuals: dict


def fem_dirichlet_solve(mesh, forcing, g, intorder=6):
    """Assemble K and F separately, eliminate Dirichlet DOFs, then solve."""
    import time
    start = time.perf_counter()
    basis = Basis(mesh, ElementTriP1(), intorder=intorder)
    K = asm(laplace, basis)

    @LinearForm
    def load(v, w):
        return forcing(w.x[0], w.x[1]) * v

    F = asm(load, basis)
    boundary = mesh.boundary_nodes()
    interior = np.setdiff1d(np.arange(mesh.p.shape[1]), boundary)
    A = K[interior][:, interior].tocsr()
    boundary_values = np.asarray(g(mesh.p[0, boundary], mesh.p[1, boundary]))
    b = F[interior] - K[interior][:, boundary] @ boundary_values
    assembly_time = time.perf_counter() - start
    direct = solve_direct(A, b)
    values = np.zeros(mesh.p.shape[1], dtype=float)
    values[boundary] = boundary_values
    values[interior] = direct.x
    return FEMResult(mesh, values, interior, boundary, A, b,
                     assembly_time, direct.wall_time, direct.residuals)


def _batch_error_integrals(vertices, original, nodal_values, X, W,
                           exact, gradient):
    """Quadrature on subtriangles, with uh affine on the original element."""
    v0 = original[:, 0, :]
    d1 = original[:, 1, :] - v0
    d2 = original[:, 2, :] - v0
    det = d1[:, 0] * d2[:, 1] - d2[:, 0] * d1[:, 1]
    q0 = vertices[:, 0, :]
    e1 = vertices[:, 1, :] - q0
    e2 = vertices[:, 2, :] - q0
    q = (q0[:, None, :]
         + e1[:, None, :] * X[0][None, :, None]
         + e2[:, None, :] * X[1][None, :, None])
    r = q - v0[:, None, :]
    lam1 = (d2[:, None, 1] * r[:, :, 0] - d2[:, None, 0] * r[:, :, 1]) / det[:, None]
    lam2 = (-d1[:, None, 1] * r[:, :, 0] + d1[:, None, 0] * r[:, :, 1]) / det[:, None]
    lam0 = 1.0 - lam1 - lam2
    uh = (nodal_values[:, 0, None] * lam0
          + nodal_values[:, 1, None] * lam1
          + nodal_values[:, 2, None] * lam2)
    ue = exact(q[:, :, 0], q[:, :, 1])
    weights = np.abs((e1[:, 0] * e2[:, 1] - e2[:, 0] * e1[:, 1]))[:, None] * W[None, :]
    l2 = np.sum(weights * (uh - ue) ** 2, axis=1)
    h1 = np.zeros_like(l2)
    if gradient is not None:
        x0, y0 = original[:, 0, 0], original[:, 0, 1]
        x1, y1 = original[:, 1, 0], original[:, 1, 1]
        x2, y2 = original[:, 2, 0], original[:, 2, 1]
        gp = np.stack((
            np.stack(((y1 - y2) / det, (x2 - x1) / det), axis=1),
            np.stack(((y2 - y0) / det, (x0 - x2) / det), axis=1),
            np.stack(((y0 - y1) / det, (x1 - x0) / det), axis=1),
        ), axis=1)
        guh = np.sum(nodal_values[:, :, None] * gp, axis=1)
        gue0, gue1 = gradient(q[:, :, 0], q[:, :, 1])
        h1 = np.sum(weights * ((guh[:, None, 0] - gue0) ** 2
                               + (guh[:, None, 1] - gue1) ** 2), axis=1)
    return np.stack((l2, h1), axis=1)


def _affine_gradient(vertices, nodal_values):
    v0, v1, v2 = vertices
    det = ((v1[0] - v0[0]) * (v2[1] - v0[1])
           - (v2[0] - v0[0]) * (v1[1] - v0[1]))
    shape_gradients = np.array((
        ((v1[1] - v2[1]) / det, (v2[0] - v1[0]) / det),
        ((v2[1] - v0[1]) / det, (v0[0] - v2[0]) / det),
        ((v0[1] - v1[1]) / det, (v1[0] - v0[0]) / det),
    ))
    return np.sum(nodal_values[:, None] * shape_gradients, axis=0)


def _singular_corner_integral(vertices, nodal_values, need_h1, epsabs, epsrel):
    """Exact radial integration for u=r^lambda sin(lambda*phi) on a corner element.

    The corner value is zero, hence uh=r*G.e_r there.  All radial integrals are
    analytic; only a smooth one-dimensional angular integral remains.  This
    avoids both pointwise singular-gradient evaluation and slow dyadic quadrature.
    """
    from scipy.integrate import quad
    from .data import singular_angle

    radii = np.linalg.norm(vertices, axis=1)
    corner = int(np.argmin(radii))
    if radii[corner] > 1e-12 or abs(nodal_values[corner]) > 1e-12:
        raise ValueError("Polar corner integration requires a zero corner value")
    other = [index for index in range(3) if index != corner]
    v1, v2 = vertices[other]
    theta_values = singular_angle(vertices[other, 0], vertices[other, 1])[0]
    if theta_values[0] > theta_values[1]:
        v1, v2 = v2, v1
        theta_values = theta_values[::-1]
    theta0, theta1 = map(float, theta_values)
    normal = np.array((v2[1] - v1[1], v1[0] - v2[0]))
    constant = float(normal @ v1)
    middle = 0.5 * (theta0 + theta1)
    if constant / (normal @ np.array((np.cos(middle), np.sin(middle)))) < 0.0:
        normal, constant = -normal, -constant

    G = _affine_gradient(vertices, nodal_values)
    lam = 2.0 / 3.0
    G2 = float(G @ G)

    def quantities(theta):
        c, s = np.cos(theta), np.sin(theta)
        er = np.array((c, s))
        et = np.array((-s, c))
        R = constant / (normal @ er)
        phi = theta - np.pi / 2.0
        sl, cl = np.sin(lam * phi), np.cos(lam * phi)
        return R, float(G @ er), float(G @ et), sl, cl

    def l2_integrand(theta):
        R, Ge, _, sl, _ = quantities(theta)
        return (Ge**2 * R**4 / 4.0
                - 2.0 * Ge * sl * R ** (lam + 3.0) / (lam + 3.0)
                + sl**2 * R ** (2.0 * lam + 2.0) / (2.0 * lam + 2.0))

    l2, error_l2 = quad(l2_integrand, theta0, theta1, epsabs=epsabs, epsrel=epsrel, limit=200)
    h1, error_h1 = 0.0, 0.0
    if need_h1:
        def h1_integrand(theta):
            R, Ge, Gt, sl, cl = quantities(theta)
            cross = Ge * sl + Gt * cl
            return (G2 * R**2 / 2.0
                    - 2.0 * lam * cross * R ** (lam + 1.0) / (lam + 1.0)
                    + lam**2 * R ** (2.0 * lam) / (2.0 * lam))
        h1, error_h1 = quad(h1_integrand, theta0, theta1, epsabs=epsabs, epsrel=epsrel, limit=200)
    return np.array((l2, h1)), error_l2 + error_h1


def fem_error_integrals(mesh, values, exact, gradient=None, intorder=8,
                        singular_corner=False, quad_epsabs=1e-13,
                        quad_epsrel=1e-11):
    """Return integral L2/H1 errors with a certified corner-quadrature route.

    Degree-8 triangle quadrature is the default away from the corner.  With
    singular_corner=True, elements having a vertex at (0,0) use analytic radial
    integration plus adaptive one-dimensional angular quadrature.  The exact
    gradient is never evaluated at the re-entrant corner.
    """
    X, W = get_quadrature(RefTri, intorder)
    original = mesh.p[:, mesh.t].transpose(2, 1, 0)
    nodal = values[mesh.t].T
    integrals = _batch_error_integrals(original, original, nodal, X, W, exact, gradient)
    angular_error = 0.0
    corner_elements = 0
    if singular_corner:
        vertex_radii = np.linalg.norm(original, axis=2)
        selected = np.flatnonzero(np.min(vertex_radii, axis=1) < 1e-12)
        corner_elements = len(selected)
        for element in selected:
            value, estimate = _singular_corner_integral(
                original[element], nodal[element], gradient is not None,
                quad_epsabs, quad_epsrel
            )
            integrals[element] = value
            angular_error += estimate
    return {
        "L2": float(np.sqrt(max(integrals[:, 0].sum(), 0.0))),
        "H1": float(np.sqrt(max(integrals[:, 1].sum(), 0.0))) if gradient is not None else np.nan,
        "Linf": float(np.max(np.abs(values - exact(mesh.p[0], mesh.p[1])))),
        "quadrature_angular_error": float(angular_error),
        "quadrature_corner_elements": int(corner_elements),
    }


def mesh_size(mesh):
    vertices = mesh.p[:, mesh.t]
    edges = np.stack((np.linalg.norm(vertices[:, 1] - vertices[:, 0], axis=0),
                      np.linalg.norm(vertices[:, 2] - vertices[:, 1], axis=0),
                      np.linalg.norm(vertices[:, 0] - vertices[:, 2], axis=0)))
    return float(edges.max())
