"""
robertson_model.py - Robertson chemical kinetics model (Topic 3).

Model (problem_pack.pdf, Section 3, authoritative definition):
    y1' = -0.04*y1 + 1e4*y2*y3
    y2' =  0.04*y1 - 1e4*y2*y3 - 3e7*y2^2
    y3' =                 3e7*y2^2
    y(0) = [1, 0, 0],  0 <= t <= 40

Mathematical properties used as the basis of the validation checks:
1. Linear invariant (mass conservation): (1,1,1).f(y) = 0, so the exact
   trajectory always satisfies y1+y2+y3 = 1.
2. The full Jacobian always has one structural zero eigenvalue (caused by the
   invariant), so the stiffness ratio must be defined on the two NONZERO
   eigenvalues only.
3. The initial state y=(1,0,0) is degenerate: y2=y3=0 makes the nonzero
   eigenvalues vanish too, so the stiffness ratio is undefined at t=0 and
   plots must start at the first positive logarithmic time.
4. Stiffness comes not from the rate constant 3e7 itself but from y2=O(1e-5)
   after the transient, which pushes |lambda|max to O(1e3-1e4).

Author: Rong Yin (Role 5: Testing & Validation)
Course: Numerical Analysis of ODEs and PDEs
"""

import numpy as np

# Reaction rate constants (given by the problem pack, do not change)
K1 = 0.04       # A -> B
K2 = 3.0e7      # B + B -> B + C (note: the term is exactly 3e7*y2^2, no extra 2)
K3 = 1.0e4      # B + C -> A + C


def robertson_rhs(t, y):
    """Right-hand side f(t, y) of the Robertson system.

    Parameters
    ----------
    t : float          current time (autonomous system, kept for API compatibility)
    y : ndarray(3,)    current state [y1, y2, y3]

    Returns
    -------
    f : ndarray(3,)    dy/dt
    """
    y1, y2, y3 = y
    return np.array([
        -K1 * y1 + K3 * y2 * y3,
        K1 * y1 - K3 * y2 * y3 - K2 * y2 * y2,
        K2 * y2 * y2,
    ])


def robertson_jacobian(t, y):
    """Analytic 3x3 Jacobian of the Robertson RHS, J[i,j] = d f_i / d y_j.

    Element-by-element differentiation:
        df1/dy1 = -0.04      df1/dy2 = 1e4*y3            df1/dy3 = 1e4*y2
        df2/dy1 =  0.04      df2/dy2 = -1e4*y3 - 6e7*y2  df2/dy3 = -1e4*y2
        df3/dy1 =  0         df3/dy2 =  6e7*y2           df3/dy3 = 0
    (note: d/dy2 of 3e7*y2^2 is 6e7*y2)
    """
    y1, y2, y3 = y
    return np.array([
        [-K1,        K3 * y3,               K3 * y2],
        [K1, -K3 * y3 - 2.0 * K2 * y2,     -K3 * y2],
        [0.0,        2.0 * K2 * y2,         0.0],
    ])


def reduced_jacobian(t, y):
    """Reduced 2x2 Jacobian after eliminating y3 = 1 - y1 - y2.

    Substituting y3 into the first two equations and differentiating with
    respect to (y1, y2). When both eigenvalues are nonzero, the stiffness
    ratio from this reduced Jacobian equals the nonzero-eigenvalue ratio of
    the full 3x3 Jacobian (the cleaner approach endorsed by the problem pack).
    """
    y1, y2, y3 = y
    return np.array([
        [-K1 - K3 * y2,  K3 * (y3 - y2)],
        [K1 - K3 * y2,  -K3 * y3 - 2.0 * K2 * y2 - K3 * y2],
    ])


def mass(y):
    """Linear invariant: mass y1+y2+y3 (exactly 1 on the exact trajectory)."""
    return np.sum(y, axis=-1)


def eigen_info(t, y, zero_threshold=1e-8):
    """Eigenvalue information and stiffness ratio at a given state, as defined
    by the problem pack.

    Parameters
    ----------
    zero_threshold : float   eigenvalues with |Re(lambda)| < threshold are
                             classified as the structural zero eigenvalue

    Returns
    -------
    lam : ndarray(3,)    full-Jacobian eigenvalues (sorted by real part)
    lam_nz : ndarray     nonzero eigenvalues (used for the stiffness ratio)
    ratio : float        S = max|Re lambda| / min|Re lambda| over the nonzero
                         spectrum; np.nan where fewer than 2 nonzero
                         eigenvalues exist (the degenerate initial state)
    """
    J = robertson_jacobian(t, y)
    lam = np.linalg.eigvals(J)
    lam = lam[np.argsort(lam.real)]
    lam_nz = lam[np.abs(lam.real) > zero_threshold]
    if len(lam_nz) < 2:
        return lam, lam_nz, np.nan
    re = np.abs(lam_nz.real)
    ratio = re.max() / re.min()
    return lam, lam_nz, ratio


# Protocol output grid (problem pack reproducible protocol):
# t=0 plus >=200 logarithmically spaced points (1e-8, 40]
def protocol_grid(n_log=200, t_end=40.0):
    """Standard output grid shared by every method and the reference
    (problem-pack reproducible protocol).

    Returns ndarray of shape (n_log+1,): t[0]=0, the rest log-spaced from
    1e-8 to t_end.
    """
    t = np.empty(n_log + 1)
    t[0] = 0.0
    t[1:] = np.logspace(-8, np.log10(t_end), n_log)
    return t
