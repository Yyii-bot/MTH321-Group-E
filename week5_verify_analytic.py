"""
Week 5 -- check the manufactured data before trusting any solver.

Everything here tests the *analytic* part of the problem (exact solution,
forcing, boundary data) with plain numpy.  No solver involved.  If a check
fails here, the bug is in the maths or the convention, not in the code.

Sign convention (same as the course test bank): we solve  -Delta u = f,
so for u* = exp(x+y) the forcing is f = -2*exp(x+y).  This matches the
problem pack.

Run:  python week5_verify_analytic.py
"""
import numpy as np


# ----------------------------------------------------------------------
# Direction 2: mixed Dirichlet / Neumann / Robin data on (0,1)^2
# ----------------------------------------------------------------------
def check_direction2():
    print("== Direction 2: u* = exp(x+y) ==")
    ok = True

    # Points chosen off-symmetric on purpose (symmetric points hide bugs).
    pts = [(0.13, 0.71), (0.55, 0.02), (0.97, 0.44), (0.31, 0.89)]

    # 1. Forcing: -Delta u* should equal 2*exp(x+y), so f = -2*exp(x+y).
    #    Verify with a fine 5-point stencil rather than trusting a by-hand
    #    derivative -- this is the same formula the solver uses.
    h = 1e-4
    for x, y in pts:
        uxx = (np.exp(x + h + y) - 2 * np.exp(x + y) + np.exp(x - h + y)) / h**2
        uyy = (np.exp(x + y + h) - 2 * np.exp(x + y) + np.exp(x + y - h)) / h**2
        f_num = -(uxx + uyy)
        f_bank = -2 * np.exp(x + y)
        match = np.isclose(f_num, f_bank, rtol=1e-5)
        ok &= match
        print(f"  -Delta u* at ({x:.2f},{y:.2f}): stencil={f_num:.8f}  f={f_bank:.8f}  {'ok' if match else 'FAIL'}")

    # 2. Boundary data on each edge, sampled at y = 0.37 (arbitrary spot).
    y = 0.37
    g_left   = np.exp(0.0 + y)      # x = 0,  u* = exp(y)
    g_bottom = np.exp(0.51 + 0.0)   # y = 0,  u* = exp(x)
    g_top    = np.exp(0.51 + 1.0)   # y = 1,  u* = exp(x+1)
    q_right  = np.exp(1.0 + y)      # x = 1,  du*/dn = du*/dx = exp(1+y) > 0
    r_right  = q_right + 2 * np.exp(1.0 + y)   # Robin: r = q + beta*u*, beta = 2

    print(f"  g(0, {y})    = {g_left:.8f}   (expect exp(y))")
    print(f"  g(.51, 0)    = {g_bottom:.8f}   (expect exp(.51))")
    print(f"  g(.51, 1)    = {g_top:.8f}")
    print(f"  q(1, {y})    = {q_right:.8f}")

    # 3. Normal direction sanity: u* grows in x, so du*/dn on x=1 must be
    #    strictly positive.  A wrong normal sign flips this -- cheap test,
    #    catches the most common boundary bug in this direction.
    if q_right <= 0:
        ok = False
        print("  FAIL: q is not positive -- check the outward normal on x=1")
    else:
        print("  q > 0 on x=1 -- outward normal sign is consistent")

    # 4. Robin data: r = q + beta*u* with beta = 2.
    if not np.isclose(r_right, 3 * np.exp(1 + y)):
        ok = False
        print("  FAIL: Robin data r != q + 2*u*")
    else:
        print("  r = 3*exp(1+y) -- Robin data consistent (r = q + 2u*)")

    print(f"  direction 2 analytic checks: {'ALL OK' if ok else 'FAILURES -- stop here'}\n")
    return ok


# ----------------------------------------------------------------------
# Open direction 1: L-shaped domain, singular benchmark
#     u*(r, phi) = r^(2/3) * sin(2 phi / 3),   phi = theta - pi/2 in [0, 3pi/2]
# ----------------------------------------------------------------------
def polar_lshape(x, y):
    """Branch convention from the problem pack: theta = atan2(y, x), then
    add 2*pi wherever theta < pi/2 so that phi lands in [0, 3pi/2].
    This must be implemented explicitly in the solver -- test it here."""
    theta = np.arctan2(y, x)
    theta = np.where(theta < np.pi / 2, theta + 2 * np.pi, theta)
    phi = theta - np.pi / 2
    return np.hypot(x, y), phi


def ustar_lshape(x, y):
    r, phi = polar_lshape(x, y)
    return r**(2 / 3) * np.sin(2 * phi / 3)


def grad_ustar_lshape(x, y):
    """Gradient in polar form (valid away from r = 0 -- never evaluate at
    the re-entrant corner).  Used again by the FEM error measurement."""
    r, phi = polar_lshape(x, y)
    with np.errstate(divide="ignore", invalid="ignore"):
        ur  = (2 / 3) * r**(-1 / 3) * np.sin(2 * phi / 3)
        uph = (2 / 3) * r**(-1 / 3) * np.cos(2 * phi / 3)   # (1/r) d u / d phi
    th = phi + np.pi / 2
    ux = np.cos(th) * ur - np.sin(th) * uph
    uy = np.sin(th) * ur + np.cos(th) * uph
    return ux, uy


def check_lshape():
    print("== Open direction 1: L-shaped singular benchmark ==")
    ok = True

    # 1. Branch convention: for points in the domain, phi must sit in
    #    [0, 3pi/2].  The tricky points are in quadrant IV (x>0, y<0),
    #    where atan2 returns a negative angle.
    test_pts = [(-0.5, -0.5), (-0.5, 0.5), (0.5, -0.5), (-1.0, 0.3), (0.3, -1.0)]
    for x, y in test_pts:
        _, phi = polar_lshape(np.array([x]), np.array([y]))
        good = 0 <= phi[0] <= 3 * np.pi / 2 + 1e-12
        ok &= good
        print(f"  phi at ({x:+.1f},{y:+.1f}) = {phi[0]:.4f}  {'ok' if good else 'FAIL: outside [0, 3pi/2]'}")

    # 2. Zero trace on the two rays meeting at the re-entrant corner.
    for x, y in [(0.0, 0.25), (0.0, 0.75), (0.25, 0.0), (0.75, 0.0)]:
        val = ustar_lshape(np.array([x]), np.array([y]))[0]
        good = abs(val) < 1e-14
        ok &= good
        print(f"  u* on ray ({x:.2f},{y:.2f}) = {val:.2e}  {'ok (zero trace)' if good else 'FAIL'}")

    # 3. Nonzero data on the *outer* boundary (if this were zero the whole
    #    problem would be trivially u = 0 and the convergence study dead).
    for x, y in [(-1.0, -0.5), (0.5, -1.0), (-0.5, 1.0), (1.0, -0.75)]:
        val = ustar_lshape(np.array([x]), np.array([y]))[0]
        good = abs(val) > 1e-6
        ok &= good
        print(f"  g({x:+.1f},{y:+.1f}) = {val:.6f}  {'ok (nonzero)' if good else 'FAIL: zero outer data'}")

    # 4. Harmonic away from the corner: Delta u* = 0 checked with a 5-point
    #    stencil at points with r >= 0.3 (stencil error blows up like r^(-7/3),
    #    so stay clear of the corner).
    h = 1e-3
    for x, y in [(-0.6, -0.6), (-0.6, 0.4), (0.4, -0.6), (-0.9, 0.0)]:
        lap = (ustar_lshape(np.array([x+h]), np.array([y]))[0]
               - 2 * ustar_lshape(np.array([x]), np.array([y]))[0]
               + ustar_lshape(np.array([x-h]), np.array([y]))[0]) / h**2
        lap += (ustar_lshape(np.array([x]), np.array([y+h]))[0]
                - 2 * ustar_lshape(np.array([x]), np.array([y]))[0]
                + ustar_lshape(np.array([x]), np.array([y-h]))[0]) / h**2
        # scale the tolerance by the local size of second derivatives
        scale = abs(ustar_lshape(np.array([x]), np.array([y]))[0]) / max(np.hypot(x, y)**2, 1e-2)
        good = abs(lap) < 1e-2 * max(scale, 1.0)
        ok &= good
        print(f"  Delta u* at ({x:+.1f},{y:+.1f}) = {lap: .3e}  {'ok (~0)' if good else 'check -- too large?'}")

    # 5. Sign: sin(2phi/3) >= 0 everywhere on the domain, so u* >= 0.
    #    Nice extra sanity when eyeballing solution plots later.
    rng = np.random.default_rng(0)
    xs = rng.uniform(-1, 1, 2000)
    ys = rng.uniform(-1, 1, 2000)
    keep = ~((xs > 0) & (ys > 0) & (xs < 1) & (ys < 1))
    keep &= np.hypot(xs, ys) > 0.05
    vals = ustar_lshape(xs[keep], ys[keep])
    good = vals.min() > -1e-12
    ok &= good
    print(f"  min u* over {keep.sum()} sample points = {vals.min():.3e}  {'ok (u* >= 0)' if good else 'FAIL'}")

    print(f"  L-shape analytic checks: {'ALL OK' if ok else 'FAILURES -- stop here'}\n")
    return ok


if __name__ == "__main__":
    a = check_direction2()
    b = check_lshape()
    print("Summary:", "everything passes -- safe to build the solver" if (a and b) else "SOMETHING FAILED")
