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
    errors near the singular corner (composite-rule comparison)
  * rates are reported against h AND against N_dof (the fair comparison
    for the later adaptive study)

Run:  python week6_lshape_fem.py
"""
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

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
# concern, not a formality: a plain 3-point rule mis-measures the L2 norm
# by ~50% near the corner.  We therefore integrate on composite rules
# (each triangle split into 4^k subtriangles, 3-point on each) and check
# one level against the next before trusting a number.
# ----------------------------------------------------------------------
def refine_once(xy, tris, uh):
    """Split every triangle into 4 via edge midpoints.  Midpoint values of
    a P1 function are exact averages, so u_h is preserved."""
    mids, mid_vals, new_nodes = {}, [], []
    def mid(a, b):
        key = (min(a, b), max(a, b))
        if key not in mids:
            mids[key] = len(xy) + len(new_nodes)
            new_nodes.append(0.5 * (xy[a] + xy[b]))
            mid_vals.append(0.5 * (uh[a] + uh[b]))
        return mids[key]
    sub = []
    for t in tris:
        a, b, c = t
        ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
        sub += [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]
    return np.vstack([xy, np.array(new_nodes)]), np.array(sub), \
        np.concatenate([uh, np.array(mid_vals)])


def _integrate_3pt(xy, tris, uh):
    """Plain 3-point quadrature of the L2 and energy errors on one
    (possibly refined) triangle list."""
    l2_sq, h1_sq = 0.0, 0.0
    p0, p1, p2 = xy[tris[:, 0]], xy[tris[:, 1]], xy[tris[:, 2]]
    area = 0.5 * np.abs((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                        - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1]))
    # 3-point rule (degree 2 exact) on each triangle
    qpts = [0.5 * (p0 + p1), 0.5 * (p1 + p2), 0.5 * (p2 + p0)]
    # grad of u_h is piecewise constant per triangle
    dlam0x = (p1[:, 1] - p2[:, 1]) / (2 * area); dlam0y = (p2[:, 0] - p1[:, 0]) / (2 * area)
    dlam1x = (p2[:, 1] - p0[:, 1]) / (2 * area); dlam1y = (p0[:, 0] - p2[:, 0]) / (2 * area)
    gh_x = dlam0x * uh[tris[:, 0]] + dlam1x * uh[tris[:, 1]] - (dlam0x + dlam1x) * uh[tris[:, 2]]
    gh_y = dlam0y * uh[tris[:, 0]] + dlam1y * uh[tris[:, 1]] - (dlam0y + dlam1y) * uh[tris[:, 2]]
    for qp in qpts:
        lam0 = ((qp[:, 0] - p1[:, 0]) * (p2[:, 1] - p1[:, 1])
                - (p2[:, 0] - p1[:, 0]) * (qp[:, 1] - p1[:, 1])) / (2 * area)
        lam1 = ((qp[:, 0] - p2[:, 0]) * (p0[:, 1] - p2[:, 1])
                - (p0[:, 0] - p2[:, 0]) * (qp[:, 1] - p2[:, 1])) / (2 * area)
        lam2 = 1.0 - lam0 - lam1
        uh_q = lam0 * uh[tris[:, 0]] + lam1 * uh[tris[:, 1]] + lam2 * uh[tris[:, 2]]
        ue_q = ustar_lshape(qp[:, 0], qp[:, 1])
        l2_sq += np.sum(area / 3.0 * (uh_q - ue_q) ** 2)
        g = grad_ustar_lshape(qp[:, 0], qp[:, 1])
        h1_sq += np.sum(area / 3.0 * ((gh_x - g[0]) ** 2 + (gh_y - g[1]) ** 2))
    return np.sqrt(l2_sq), np.sqrt(h1_sq)


def integrate_errors(xy, tris, uh, base=1, corner=4, radius=0.5, mask=None):
    """Errors with corner-aware composite quadrature.

    The integrands (u*-u_h)^2 and |grad(u*-u_h)|^2 are singular at the
    re-entrant corner, and most of the L2 error lives there.  Triangles
    within `radius` of the corner get `corner` levels of 4-way subdivision,
    the rest get `base`.  Subsets are refined and integrated independently,
    which is fine for integration even where refinement levels differ
    across a shared edge.  `mask` optionally restricts to a triangle subset.
    """
    if mask is not None:
        tris = tris[mask]
    dmin = np.abs(xy[tris]).min(axis=(1, 2))      # nearest vertex to the corner
    levels = np.where(dmin < radius, corner, base)
    l2_sq, h1_sq = 0.0, 0.0
    for lv in np.unique(levels):
        xyl, tsl, uhl = xy, tris[levels == lv], uh
        for _ in range(int(lv)):
            xyl, tsl, uhl = refine_once(xyl, tsl, uhl)
        a, b = _integrate_3pt(xyl, tsl, uhl)
        l2_sq += a * a
        h1_sq += b * b
    return np.sqrt(l2_sq), np.sqrt(h1_sq)


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
    print("| n | h | N_dof | L2 error | energy error | L2 quad check | energy quad check | solve s |")
    print("|---|---|---|---|---|---|---|---|")
    prev = None
    for n in (8, 16, 32, 64, 128):   # add 256 if the rates are not yet stable
        xy, tris, interior, h, uh, t_solve = solve_level(n)
        # table values: 4 corner levels; check column: 3 vs 4 corner levels
        # everything within r < 0.5 of the corner is composite-refined;
        # outside that the integrand is smooth enough for plain 3-point
        l2, en = integrate_errors(xy, tris, uh, corner=4)
        l2c, enc = integrate_errors(xy, tris, uh, corner=3)
        ndof = interior.sum()
        print(f"| {n} | {h:.4f} | {ndof} | {l2:.4e} | {en:.4e} "
              f"| {abs(l2-l2c)/l2:.1e} | {abs(en-enc)/en:.1e} | {t_solve:.2f} |")

        # where does the L2 error live? split at r = 0.25 around the corner
        near = np.abs(xy[tris]).min(axis=(1, 2)) < 0.25
        l2_near, _ = integrate_errors(xy, tris, uh, corner=3, mask=near)
        l2_far, _ = integrate_errors(xy, tris, uh, mask=~near)
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
    print("  * the quad check column should stay small; if it grows on fine")
    print("    meshes the corner singularity is polluting the error measurement")
    print("    itself, and the quadrature order must be raised before claiming")
    print("    any rate.")
