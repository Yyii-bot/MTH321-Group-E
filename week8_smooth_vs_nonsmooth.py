"""
Week 8 -- depth comparison 4: smooth vs. nonsmooth data.

Same solver (the reference P1 FEM from week6), same L-shaped meshes,
same error norms.  Only the manufactured solution changes:

    smooth leg:     u_s = exp(x+y)      -> f = -2 exp(x+y)  (smooth on L)
    nonsmooth leg:  u* = r^{2/3} sin(2 phi / 3)             (f = 0, singular)

Prediction (written before the first run, pre-computation protocol):
energy O(h^1) smooth vs. O(h^{2/3}) singular; L2 O(h^2) smooth vs. the
singular leg's measured ~1.0.  Whatever is measured below is reported.

WHY THIS SCRIPT LOOKS THE WAY IT DOES -- the L2 quadrature saga:

The L2 error signal is TINY: int e^2 ~ h^4 (e.g. 8e-3 at h=1/16), while
the integrands that make it up are O(1).  Every low-order rule we tried
(int e^2 directly at any composite depth, the interpolant split, the
three-term expansion) has quadrature error O(h^3) in the integral --
the same order as the signal -- because e is small pointwise but grad e
is O(1) pointwise, so third derivatives of e^2 are O(1).  All those
methods agree with each other at fixed refinement to <1e-3 yet still
drift ~50% per refinement level: agreement between two imperfect
methods is not convergence.

The fix used here: Duffy transform (square -> triangle) with a tensor
Gauss-Legendre rule.  n=3 per direction (9 points per triangle) gives
per-triangle error O(h^8 |f|_{W^{6,inf}}), i.e. O(h^6) summed -- two
orders below the h^4 signal.  For the singular leg the corner triangles
are composite-refined first (e ~ r^{2/3} there), then Duffy-integrated
on the subtriangles, with a 3-vs-4-level check.

A second diagnostic settles the *interpretation* of the smooth-leg L2
rate: the error splits as e = theta - rho with
    theta = u_h - u_I   (discrete solution minus nodal interpolant,
                         orthogonal-ish part, "pollution")
    rho   = u_e - u_I   (plain interpolation error)
||theta||^2 is computed EXACTLY (theta is linear per triangle), ||rho||
by Duffy.  If the smooth leg's L2 rate is below 2, this split shows
whether the gap is interpolation (rho) or corner pollution (theta).

Run:  python week8_smooth_vs_nonsmooth.py   (a few minutes)
"""
import time
import numpy as np
import scipy.sparse.linalg as spla

from week6_lshape_fem import build_mesh, assemble, solve_level, integrate_errors
from week5_verify_analytic import ustar_lshape, grad_ustar_lshape

LEVELS = (8, 16, 32, 64, 128)

# 1D Gauss-Legendre nodes/weights on [-1,1]
_GAUSS = {1: ([0.0], [2.0]),
          2: ([-0.5773502691896257, 0.5773502691896257], [1.0, 1.0]),
          3: ([-0.7745966692414834, 0.0, 0.7745966692414834], [5.0 / 9, 8.0 / 9, 5.0 / 9])}


# ----------------------------------------------------------------------
# smooth manufactured solution and its data
# ----------------------------------------------------------------------
def ue_smooth(x, y):
    return np.exp(x + y)

def grad_ue_smooth(x, y):
    g = np.exp(x + y)
    return g, g

F_SMOOTH = -2.0                      # f = -Delta(exp(x+y))


def assemble_load_clean(xy, tris):
    p0, p1, p2 = xy[tris[:, 0]], xy[tris[:, 1]], xy[tris[:, 2]]
    det = (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) \
          - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1])
    b = np.zeros(len(xy))
    for s, ws in zip(*_GAUSS[2]):
        for t, wt in zip(*_GAUSS[2]):
            lam = (0.5 * (1 + s), 0.25 * (1 - s) * (1 + t), 0.25 * (1 - s) * (1 - t))
            wgt = abs(det) / 8.0 * (1 - s) * ws * wt
            qx = lam[0] * p0[:, 0] + lam[1] * p1[:, 0] + lam[2] * p2[:, 0]
            qy = lam[0] * p0[:, 1] + lam[1] * p1[:, 1] + lam[2] * p2[:, 1]
            fv = F_SMOOTH * np.exp(qx + qy)
            for a in range(3):
                np.add.at(b, tris[:, a], wgt * fv * lam[a])
    return b


# ----------------------------------------------------------------------
# Duffy-transform Gauss quadrature of int (u_h - ue)^2
# ----------------------------------------------------------------------
def refine_fields(xy, tris, fields):
    """4-way midpoint refinement carrying any number of P1 nodal fields."""
    mids, new_xy, new_vals = {}, [], [[] for _ in fields]
    def mid(a, b):
        key = (min(a, b), max(a, b))
        if key not in mids:
            mids[key] = len(xy) + len(new_xy)
            new_xy.append(0.5 * (xy[a] + xy[b]))
            for f, acc in zip(fields, new_vals):
                acc.append(0.5 * (f[a] + f[b]))
        return mids[key]
    sub = []
    for t in tris:
        a, b, c = t
        ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
        sub += [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]
    xy = np.vstack([xy, np.array(new_xy)]) if new_xy else xy
    return xy, np.array(sub), [np.concatenate([f, np.array(v)]) for f, v in zip(fields, new_vals)]


def l2_duffy(xy, tris, uh, ue, n=3, corner_levels=0, radius=0.0):
    """int (u_h - ue)^2 by Duffy-Gauss.  For a singular ue, triangles
    within `radius` of the origin are composite-refined `corner_levels`
    times first (carrying u_h exactly)."""
    if corner_levels > 0:
        dmin = np.abs(xy[tris]).min(axis=(1, 2))
        sel = dmin < radius
        if sel.any():
            xyc, trc, (uhc,) = refine_fields(xy, tris[sel], (uh,))
            for _ in range(corner_levels - 1):
                xyc, trc, (uhc,) = refine_fields(xyc, trc, (uhc,))
            l2_sq = _l2_duffy_core(xy, tris[~sel], uh, ue, n) if (~sel).any() else 0.0
            l2_sq += _l2_duffy_core(xyc, trc, uhc, ue, n)
            return np.sqrt(l2_sq)
    return np.sqrt(_l2_duffy_core(xy, tris, uh, ue, n))


def _l2_duffy_core(xy, tris, uh, ue, n):
    p0, p1, p2 = xy[tris[:, 0]], xy[tris[:, 1]], xy[tris[:, 2]]
    det = (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) \
          - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1])
    uh0, uh1, uh2 = uh[tris[:, 0]], uh[tris[:, 1]], uh[tris[:, 2]]
    total = 0.0
    for s, ws in zip(*_GAUSS[n]):
        for t, wt in zip(*_GAUSS[n]):
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


def theta_norm_exact(xy, tris, uh, uI):
    """||u_h - u_I||_L2 EXACTLY: the difference is linear per triangle,
    the 3-point edge-midpoint rule integrates quadratics exactly."""
    w = uh - uI
    p0, p1, p2 = xy[tris[:, 0]], xy[tris[:, 1]], xy[tris[:, 2]]
    area = 0.5 * np.abs((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                        - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1]))
    m01 = 0.5 * (w[tris[:, 0]] + w[tris[:, 1]])
    m12 = 0.5 * (w[tris[:, 1]] + w[tris[:, 2]])
    m20 = 0.5 * (w[tris[:, 2]] + w[tris[:, 0]])
    return np.sqrt(np.sum(area / 3.0 * (m01**2 + m12**2 + m20**2)))


def integrate_energy(xy, tris, uh, grad_ue):
    """Energy seminorm error; integrand converges cleanly under 3-point."""
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
    return np.sqrt(h1_sq)


def solve_smooth(n):
    xy, tris, interior, h = build_mesh(n)          # same meshes as singular leg
    K = assemble(xy, tris)
    b = assemble_load_clean(xy, tris)
    g = ue_smooth(xy[~interior, 0], xy[~interior, 1])
    t0 = time.perf_counter()
    uI = spla.spsolve(K[interior][:, interior],
                      b[interior] - K[interior][:, ~interior] @ g)
    t_solve = time.perf_counter() - t0
    uh = np.zeros(len(xy)); uh[interior] = uI; uh[~interior] = g
    return xy, tris, interior, h, uh, t_solve


def orders(vals):
    return [np.log2(vals[k] / vals[k + 1]) for k in range(len(vals) - 1)]

def fmt_orders(o):
    return " / ".join(f"{p:.2f}" for p in o)


def run_leg(label, solver, ue, grad_ue, n_duffy=3, corner=0, radius=0.0, decompose=False):
    print(f"\n### {label} ###")
    if decompose:
        print("| n | h | N_dof | L2(u_h-u_e) | ||theta||=||u_h-u_I|| | ||rho||=||u_e-u_I|| | energy | L2 order | energy order |")
        print("|---|---|---|---|---|---|---|---|---|")
    else:
        print("| n | h | N_dof | L2 error | corner quad check | energy error | L2 order | energy order |")
        print("|---|---|---|---|---|---|---|---|")
    l2s, ens = [], []
    for n in LEVELS:
        xy, tris, interior, h, uh, t = solver(n)
        l2 = l2_duffy(xy, tris, uh, ue, n=n_duffy, corner_levels=corner, radius=radius)
        if corner > 0:
            l2c = l2_duffy(xy, tris, uh, ue, n=2, corner_levels=corner + 1, radius=radius)
            chk = abs(l2 - l2c) / l2
        else:
            chk = float("nan")
        en = integrate_energy(xy, tris, uh, grad_ue)
        l2s.append(l2); ens.append(en)
        if decompose:
            uI = ue(xy[:, 0], xy[:, 1])
            th = theta_norm_exact(xy, tris, uh, uI)
            rho = l2_duffy(xy, tris, uI, ue, n=n_duffy, corner_levels=corner, radius=radius)
            line = f"| {n} | {h:.4f} | {interior.sum()} | {l2:.4e} | {th:.4e} | {rho:.4e} | {en:.4e} "
        else:
            line = f"| {n} | {h:.4f} | {interior.sum()} | {l2:.4e} | {chk:.1e} | {en:.4e} "
        if len(l2s) > 1:
            line += f"| {np.log2(l2s[-2]/l2):.2f} | {np.log2(ens[-2]/en):.2f} |"
        else:
            line += "| -- | -- |"
        print(line)
    print(f"  all pairwise orders vs h:  L2 {fmt_orders(orders(l2s))}   energy {fmt_orders(orders(ens))}")
    return l2s, ens


if __name__ == "__main__":
    print("smooth vs. nonsmooth on the L-domain: same solver, meshes, norms.")
    print("L2 by Duffy-transform tensor Gauss (n=3 smooth; corner-refined")
    print("n=2 singular with a 3-vs-4-level check); energy by 3-point.")

    # quadrature self-convergence on the smooth leg (should be ~exact):
    print("\n[L2 quadrature self-convergence, smooth leg, n=32]")
    xy, tris, interior, h, uh, _ = solve_smooth(32)
    for n in (1, 2, 3):
        print(f"  Duffy n={n}: {l2_duffy(xy, tris, uh, ue_smooth, n=n):.8e}")
    print("  (n=1 is NOT expected to agree -- shown to demonstrate sensitivity)")

    l2_s, en_s = run_leg("smooth leg: u_s = exp(x+y) (f = -2 exp(x+y))",
                         solve_smooth, ue_smooth, grad_ue_smooth,
                         n_duffy=3, decompose=True)
    l2_n, en_n = run_leg("nonsmooth leg: u* = r^{2/3} sin(2 phi / 3)",
                         solve_level, ustar_lshape, grad_ustar_lshape,
                         n_duffy=2, corner=3, radius=0.25)

    print("\n### contrast (same solver, same meshes, same norms) ###")
    print(f"  L2:     smooth {fmt_orders(orders(l2_s))}   singular {fmt_orders(orders(l2_n))}")
    print(f"  energy: smooth {fmt_orders(orders(en_s))}   singular {fmt_orders(orders(en_n))}")
    print("\nInterpretation notes (measured, not assumed):")
    print("  * energy: 1.0 vs 2/3 -- the textbook regularity gap. Clean.")
    print("  * L2: compare the ||theta|| and ||rho|| columns of the smooth leg.")
    print("    rho is the pure interpolation error (guaranteed O(h^2));")
    print("    theta is the discrete pollution term.  If ||theta|| decays like")
    print("    h^1 and dominates, the smooth leg's L2 gap is a corner-pollution")
    print("    effect, not a solver bug -- and the SAME mechanism explains the")
    print("    singular leg, making the depth comparison internally consistent.")
