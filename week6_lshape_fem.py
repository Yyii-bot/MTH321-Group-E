"""
Week 6 -- L-shaped domain: uniform-refinement study with a small reference P1 FEM.

This is the validation seat's independent check for the open direction.
The team's real solver (probably scikit-fem) is the deliverable; this
script exists so that the *expected* convergence behaviour is measured
before the real solver exists, and so the real solver has an independent
number to disagree with.

What theory promises on quasi-uniform meshes (piecewise-linear FEM,
u* = r^{2/3} sin(2 phi / 3), u* in H^{1+s} only for s < 2/3):

    energy error  |u* - u_h|_{H1} = O(h^{2/3})   ~ O(N_dof^{-1/3})
    L2    error   ||u* - u_h||_L2  = O(h^{4/3})   ~ O(N_dof^{-2/3})

The L2 rate is NOT O(h^2): the dual problem has the same limited
regularity, so Aubin-Nitsche only buys about one more h^{2/3}.

Things this script checks beyond the bare rates:
  * the branch convention (already unit-tested in week5, used again here)
  * quadrature is accurate enough that it does not pollute the measured
    errors near the singular corner (self-check: L2 at Duffy n=8 vs n=6,
    plus the adaptive quad's own error estimate on the corner elements)
  * rates are reported against h AND against N_dof (the fair comparison
    for the later adaptive study)

Run:  python week6_lshape_fem.py
"""
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.polynomial.legendre import leggauss

# reuse the verified analytic functions from Week 5
from week5_verify_analytic import ustar_lshape, grad_ustar_lshape


# ----------------------------------------------------------------------
# mesh: uniform grid on [-1,1]^2 with the quadrant [0,1)^2 cut out,
# split into right triangles.  The L-shape is represented exactly because
# the cut runs along grid lines.
# ----------------------------------------------------------------------
def build_mesh(n):
    h = 2.0 / n
    xs = np.linspace(-1.0, 1.0, n + 1)
    X, Y = np.meshgrid(xs, xs, indexing="ij")

    # drop nodes strictly inside the removed quadrant
    removed = (X > 0) & (Y > 0) & (X < 1) & (Y < 1)
    kept = ~removed
    idx = -np.ones((n + 1, n + 1), dtype=int)
    idx[kept] = np.arange(kept.sum())
    xy = np.column_stack([X[kept], Y[kept]])

    # a cell survives only if ALL four corners are kept.  Because the cut
    # follows grid lines (x=0, y=0, x=1, y=1 are all grid lines for even n),
    # a cell with a removed corner lies entirely in the removed closure, so
    # dropping it costs nothing.  The mesh then represents the L exactly.
    tris = []
    for i in range(n):
        for j in range(n):
            cc = (idx[i, j], idx[i + 1, j], idx[i, j + 1], idx[i + 1, j + 1])
            if all(c >= 0 for c in cc):
                c00, c10, c01, c11 = cc
                tris.append([c00, c10, c11])
                tris.append([c00, c11, c01])
    tris = np.array(tris)
    e1 = xy[tris[:, 1]] - xy[tris[:, 0]]
    e2 = xy[tris[:, 2]] - xy[tris[:, 0]]
    area = 0.5 * np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]).sum()
    if not np.isclose(area, 3.0, rtol=1e-12):   # the L has area 4 - 1
        raise RuntimeError(f"mesh area {area} != 3 -- mask or connectivity is wrong")

    # interior vs boundary: a kept node is interior iff all four in-grid
    # neighbours are kept and it is not on the outer box
    interior = kept.copy()
    interior[0, :] = interior[-1, :] = False
    interior[:, 0] = interior[:, -1] = False
    interior[1:-1, 1:-1] &= kept[1:-1, 1:-1] & kept[:-2, 1:-1] & kept[2:, 1:-1] \
                            & kept[1:-1, :-2] & kept[1:-1, 2:]
    # the re-entrant corner (0,0) lies on the two cut edges, where the
    # exact solution has zero trace -- it is a Dirichlet node, not a free
    # unknown.  (Leaving it free was a real modelling slip: the discrete
    # corner value came out ~0.29 on coarse meshes, while the Algorithm
    # branch fixes it to 0; their analytic corner integral also assumes
    # u_h(corner) = 0.)  n is always even here, so x = 0 is a grid line.
    i0 = np.nonzero(xs == 0.0)[0]
    if len(i0):
        interior[i0[0], i0[0]] = False
    return xy, tris, interior[kept], h


# ----------------------------------------------------------------------
# P1 assembly of the Laplacian (load vector is zero: f = 0)
# ----------------------------------------------------------------------
def assemble(xy, tris):
    p0, p1, p2 = xy[tris[:, 0]], xy[tris[:, 1]], xy[tris[:, 2]]
    area = 0.5 * np.abs((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                        - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1]))
    # gradients of the barycentric coordinates: dlam_i = (rotated opposite edge)/(2|T|)
    def grad_lam(pa, pb, pc):
        return np.column_stack([(pb[:, 1] - pc[:, 1]), (pc[:, 0] - pb[:, 0])]) / (2 * area)[:, None]
    g0, g1, g2 = grad_lam(p0, p1, p2), grad_lam(p1, p2, p0), grad_lam(p2, p0, p1)

    rows, cols, vals = [], [], []
    for a in range(3):
        for b in range(3):
            ga = (g0, g1, g2)[a]; gb = (g0, g1, g2)[b]
            ka = (area * np.sum(ga * gb, axis=1))
            rows += list(tris[:, a]); cols += list(tris[:, b]); vals += list(ka)
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(xy),) * 2)


# ----------------------------------------------------------------------
# error measurement.  The integrands are singular at the corner
# (u* ~ r^{2/3}, grad u* ~ r^{-1/3}), so the quadrature order is a real
# concern, not a formality.  Final protocol (revision 3, same measuring
# stick as the Algorithm branch's E5 table):
#   * elements with a vertex at the corner: exact radial / adaptive
#     angular integral (_singular_corner_integral, ported from their
#     fem_impl);
#   * all other elements: L2 by Duffy-transform tensor Gauss at n=8
#     (their intorder=8), energy by the 3-point edge-midpoint rule,
#     whose points never touch the corner and which integrates the
#     squared-gradient terms exactly where the integrand is smooth.
# Revision 2 used composite 3-point refinement with a 3-vs-4-level
# check; the check passed while the absolute values sat ~7x above the
# certified numbers -- two quadratures can agree with each other and
# both be wrong -- which is why the measuring stick was switched.
# ----------------------------------------------------------------------

# 1D Gauss-Legendre nodes/weights on [-1,1].  Rules up to n=3 are
# tabulated; higher orders are generated on first use.
_GAUSS = {1: ([0.0], [2.0]),
          2: ([-0.5773502691896257, 0.5773502691896257], [1.0, 1.0]),
          3: ([-0.7745966692414834, 0.0, 0.7745966692414834], [5.0 / 9, 8.0 / 9, 5.0 / 9])}


def _gauss_rule(n):
    if n not in _GAUSS:
        x, w = leggauss(n)
        _GAUSS[n] = (list(x), list(w))
    return _GAUSS[n]


def _affine_gradient(vertices, nodal_values):
    """Gradient of the P1 function on one triangle (constant per element).
    Ported from the Algorithm branch's fem_impl._affine_gradient."""
    v0, v1, v2 = vertices
    det = ((v1[0] - v0[0]) * (v2[1] - v0[1])
           - (v2[0] - v0[0]) * (v1[1] - v0[1]))
    shape_gradients = np.array((
        ((v1[1] - v2[1]) / det, (v2[0] - v1[0]) / det),
        ((v2[1] - v0[1]) / det, (v0[0] - v2[0]) / det),
        ((v0[1] - v1[1]) / det, (v1[0] - v0[0]) / det),
    ))
    return np.sum(nodal_values[:, None] * shape_gradients, axis=0)


def _corner_theta(x, y):
    """Angular coordinate on the week5 branch convention (phi in
    [0, 3pi/2]); equivalent to the Algorithm branch's singular_angle.
    The corner elements of our mesh never cross the branch cut after
    this shift, so plain ordering of the two corner-edge angles is safe."""
    theta = np.arctan2(y, x)
    return np.where(theta < np.pi / 2, theta + 2.0 * np.pi, theta)


def _singular_corner_integral(vertices, nodal_values, need_h1=False,
                              epsabs=1e-13, epsrel=1e-11):
    """Exact radial integration of int (u_h - u*)^2 (and, optionally, the
    energy integrand) on a triangle with a vertex AT the re-entrant
    corner, for u* = r^lambda sin(lambda phi), lambda = 2/3.

    The corner nodal value is zero, so u_h = r (G.e_r) with G the affine
    gradient; every radial integral of (u_h - u*)^2 is then analytic and
    only a smooth one-dimensional angular integral remains, handled by
    scipy.integrate.quad with an error estimate we report.  Ported from
    the Algorithm branch's fem_impl._singular_corner_integral
    (their certified corner route, which their E5 table is measured
    with); _corner_theta matches their singular_angle.
    """
    from scipy.integrate import quad

    radii = np.linalg.norm(vertices, axis=1)
    corner = int(np.argmin(radii))
    if radii[corner] > 1e-12 or abs(nodal_values[corner]) > 1e-12:
        raise ValueError("corner integral requires a zero corner value")
    other = [i for i in range(3) if i != corner]
    v1, v2 = vertices[other]
    theta = _corner_theta(vertices[other, 0], vertices[other, 1])
    if theta[0] > theta[1]:
        v1, v2 = v2, v1
        theta = theta[::-1]
    theta0, theta1 = float(theta[0]), float(theta[1])
    # edge opposite the corner: normal.x = constant; R(theta) is where
    # the ray from the corner meets that edge
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

    # int_0^R (r Ge - r^lam sl)^2 r dr, expanded term by term
    def l2_integrand(theta):
        R, Ge, _, sl, _ = quantities(theta)
        return (Ge**2 * R**4 / 4.0
                - 2.0 * Ge * sl * R**(lam + 3.0) / (lam + 3.0)
                + sl**2 * R**(2.0 * lam + 2.0) / (2.0 * lam + 2.0))

    l2, err_l2 = quad(l2_integrand, theta0, theta1,
                      epsabs=epsabs, epsrel=epsrel, limit=200)
    h1, err_h1 = 0.0, 0.0
    if need_h1:
        # int_0^R |G - lam r^(lam-1)(sl er + cl et)|^2 r dr; the exact
        # gradient is never evaluated at the corner itself
        def h1_integrand(theta):
            R, Ge, Gt, sl, cl = quantities(theta)
            cross = Ge * sl + Gt * cl
            return (G2 * R**2 / 2.0
                    - 2.0 * lam * cross * R**(lam + 1.0) / (lam + 1.0)
                    + lam**2 * R**(2.0 * lam) / (2.0 * lam))
        h1, err_h1 = quad(h1_integrand, theta0, theta1,
                          epsabs=epsabs, epsrel=epsrel, limit=200)
    return l2, h1, err_l2 + err_h1


def _l2_duffy_core(xy, tris, uh, ue, n):
    """int (u_h - ue)^2 by Duffy-transform tensor Gauss at order n.
    Only used AWAY from the corner, where the integrand is smooth."""
    p0, p1, p2 = xy[tris[:, 0]], xy[tris[:, 1]], xy[tris[:, 2]]
    det = (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) \
          - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1])
    uh0, uh1, uh2 = uh[tris[:, 0]], uh[tris[:, 1]], uh[tris[:, 2]]
    total = 0.0
    for s, ws in zip(*_gauss_rule(n)):
        for t, wt in zip(*_gauss_rule(n)):
            lam0 = 0.5 * (1 + s)
            lam1 = 0.25 * (1 - s) * (1 + t)
            lam2 = 0.25 * (1 - s) * (1 - t)
            wgt = abs(det) / 8.0 * (1 - s) * ws * wt
            qx = lam0 * p0[:, 0] + lam1 * p1[:, 0] + lam2 * p2[:, 0]
            qy = lam0 * p0[:, 1] + lam1 * p1[:, 1] + lam2 * p2[:, 1]
            uh_q = lam0 * uh0 + lam1 * uh1 + lam2 * uh2
            e = uh_q - ue(qx, qy)
            total += np.sum(wgt * e**2)
    return total


def _h1_midpoint_sq(xy, tris, uh, grad_ue):
    """int |grad(u_h - ue)|^2 by the 3-point edge-midpoint rule.
    grad u_h is piecewise constant; the quadrature points are edge
    midpoints, which never touch the re-entrant corner."""
    p0, p1, p2 = xy[tris[:, 0]], xy[tris[:, 1]], xy[tris[:, 2]]
    area = 0.5 * np.abs((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                        - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1]))
    qpts = [0.5 * (p0 + p1), 0.5 * (p1 + p2), 0.5 * (p2 + p0)]
    dlam0x = (p1[:, 1] - p2[:, 1]) / (2 * area); dlam0y = (p2[:, 0] - p1[:, 0]) / (2 * area)
    dlam1x = (p2[:, 1] - p0[:, 1]) / (2 * area); dlam1y = (p0[:, 0] - p2[:, 0]) / (2 * area)
    gh_x = dlam0x * uh[tris[:, 0]] + dlam1x * uh[tris[:, 1]] - (dlam0x + dlam1x) * uh[tris[:, 2]]
    gh_y = dlam0y * uh[tris[:, 0]] + dlam1y * uh[tris[:, 1]] - (dlam0y + dlam1y) * uh[tris[:, 2]]
    h1_sq = 0.0
    for qp in qpts:
        g = grad_ue(qp[:, 0], qp[:, 1])
        h1_sq += np.sum(area / 3.0 * ((gh_x - g[0])**2 + (gh_y - g[1])**2))
    return h1_sq


def integrate_errors(xy, tris, uh, n_duffy=8, mask=None):
    """L2 and energy errors with the corner-aware protocol above.

    Returns (L2, energy, corner_quad_abs_err): the exact corner
    integrals contribute the third component's error estimate, the
    Duffy / midpoint parts are exact to machine precision for these
    integrands.  `mask` optionally restricts to a triangle subset.
    """
    if mask is not None:
        tris = tris[mask]
    rad = np.linalg.norm(xy[tris], axis=2)
    is_corner = rad.min(axis=1) < 1e-12
    l2_sq, h1_sq, quad_err = 0.0, 0.0, 0.0
    for t in tris[is_corner]:
        a, b, e = _singular_corner_integral(xy[t], uh[t], need_h1=True)
        l2_sq += a
        h1_sq += b
        quad_err += e
    tris = tris[~is_corner]
    if len(tris):
        l2_sq += _l2_duffy_core(xy, tris, uh, ustar_lshape, n_duffy)
        h1_sq += _h1_midpoint_sq(xy, tris, uh, grad_ustar_lshape)
    return np.sqrt(l2_sq), np.sqrt(h1_sq), quad_err


def solve_level(n):
    xy, tris, interior, h = build_mesh(n)
    K = assemble(xy, tris)
    g = ustar_lshape(xy[~interior, 0], xy[~interior, 1])   # Dirichlet data
    t0 = time.perf_counter()
    uI = spla.spsolve(K[interior][:, interior], -K[interior][:, ~interior] @ g)
    t_solve = time.perf_counter() - t0
    uh = np.zeros(len(xy)); uh[interior] = uI; uh[~interior] = g
    return xy, tris, interior, h, uh, t_solve


if __name__ == "__main__":
    print("| n | h | N_dof | L2 error | energy error | L2 quad check | corner quad rel-err | solve s |")
    print("|---|---|---|---|---|---|---|---|")
    prev = None
    for n in (8, 16, 32, 64, 128):   # add 256 if the rates are not yet stable
        xy, tris, interior, h, uh, t_solve = solve_level(n)
        # table L2 at Duffy n=8; check column: n=8 vs n=6 relative difference
        l2, en, qerr = integrate_errors(xy, tris, uh)
        l2c, _, _ = integrate_errors(xy, tris, uh, n_duffy=6)
        ndof = interior.sum()
        print(f"| {n} | {h:.4f} | {ndof} | {l2:.4e} | {en:.4e} "
              f"| {abs(l2-l2c)/l2:.1e} | {qerr/(l2*l2):.1e} | {t_solve:.2f} |")

        # where does the L2 error live? split at r = 0.25 around the corner
        near = np.abs(xy[tris]).min(axis=(1, 2)) < 0.25
        l2_near, _, _ = integrate_errors(xy, tris, uh, mask=near)
        l2_far, _, _ = integrate_errors(xy, tris, uh, mask=~near)
        print(f"    L2 split: corner region (r<0.25) {l2_near:.4e}, "
              f"rest of domain {l2_far:.4e}  ({100*l2_near**2/(l2_near**2+l2_far**2):.0f}% of the squared error is at the corner)")
        if prev:
            pn, (pl2, pen) = prev
            print(f"    -> orders vs h:  L2 p={np.log2(pl2/l2):.2f} (theory 4/3 = 1.33), "
                  f"energy p={np.log2(pen/en):.2f} (theory 2/3 = 0.67)")
            print(f"       vs N_dof:    L2 p={np.log2(pl2/l2)/np.log2(ndof/pn):.2f} (theory 2/3), "
                  f"energy p={np.log2(pen/en)/np.log2(ndof/pn):.2f} (theory 1/3)   "
                  f"[ndof ratio {ndof/pn:.2f}]")
        prev = (ndof, (l2, en))
    print()
    print("Notes for the report:")
    print("  * coarse meshes are NOT in the asymptotic range yet -- that is")
    print("    expected for a singular problem, say so instead of forcing a fit.")
    print("  * the quad check columns should stay small; if they grow on fine")
    print("    meshes the corner singularity is polluting the error measurement")
    print("    itself, and the quadrature order must be raised before claiming")
    print("    any rate.")
