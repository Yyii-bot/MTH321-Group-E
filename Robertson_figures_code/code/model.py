"""Robertson chemical kinetics model (1966 benchmark).

System (mass-action kinetics, rate constants k1=0.04, k2=3e7, k3=1e4):
    y1' = -0.04*y1 + 1e4*y2*y3
    y2' =  0.04*y1 - 1e4*y2*y3 - 3e7*y2^2
    y3' =  3e7*y2^2
    y(0) = (1, 0, 0),  0 <= t <= 40.

Invariant: y1+y2+y3 = 1 (the vector field is orthogonal to (1,1,1)).
Hence the full 3x3 Jacobian always has a zero eigenvalue; the two
nonzero eigenvalues coincide with those of the reduced 2x2 Jacobian
obtained by eliminating y3 = 1 - y1 - y2.
"""
import numpy as np

K1, K2, K3 = 0.04, 3.0e7, 1.0e4
Y0 = np.array([1.0, 0.0, 0.0])
T_FINAL = 40.0


def rhs(t, y):
    """Right-hand side f(t, y) of the Robertson system."""
    y1, y2, y3 = y
    return np.array([-K1 * y1 + K3 * y2 * y3,
                      K1 * y1 - K3 * y2 * y3 - K2 * y2**2,
                      K2 * y2**2])


def jac(t, y):
    """Full 3x3 Jacobian df/dy."""
    y1, y2, y3 = y
    return np.array([[-K1,  K3 * y3,            K3 * y2],
                     [ K1, -K3 * y3 - 2 * K2 * y2, -K3 * y2],
                     [0.0,  2 * K2 * y2,        0.0]])


def jac_reduced(y1, y2):
    """2x2 Jacobian of the reduced system with y3 = 1 - y1 - y2 eliminated.

    Its two eigenvalues are exactly the two nonzero eigenvalues of the
    full 3x3 Jacobian; used for the stiffness-ratio diagnostic.
    """
    return np.array([[-K1 - K3 * y2,  K3 * (1 - y1 - 2 * y2)],
                     [ K1 + K3 * y2, -K3 * (1 - y1 - 2 * y2) - 2 * K2 * y2]])


def output_grid(n_log=241):
    """Protocol output times: t=0 plus >=200 log-spaced times in [1e-8, 40]."""
    return np.concatenate([[0.0], np.logspace(-8, np.log10(T_FINAL), n_log)])


# ----- Prothero-Robinson test problem (order verification) -----
def prothero_robinson(lam):
    """y' = lam*(y - exp(-t)) - exp(-t),  y(0)=1,  exact solution y=exp(-t)."""
    f = lambda t, y: lam * (y - np.exp(-t)) - np.exp(-t)
    jf = lambda t, y: np.array([[lam]])
    exact = lambda t: np.exp(-t)
    return f, jf, exact
