"""Analytic data for all manufactured boundary-value problems.

Every function accepts NumPy arrays as well as scalars.  Keeping the physics in
this module prevents mesh code and experiment drivers from defining competing
versions of the same data.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class CoreData:
    """Data for -Delta(u)=f with a Neumann or Robin right edge."""

    beta: float = 2.0

    @staticmethod
    def exact(x, y):
        return np.exp(x + y)

    @staticmethod
    def forcing(x, y):
        return -2.0 * np.exp(x + y)

    @staticmethod
    def dirichlet(x, y):
        return np.exp(x + y)

    @staticmethod
    def neumann_right(y):
        # Outward normal on x=1 is (1,0), hence u_n=u_x.
        return np.exp(1.0 + y)

    def robin_right(self, y):
        # r=q+beta*u on x=1 and beta=2, hence r=3*exp(1+y).
        return self.neumann_right(y) + self.beta * self.exact(1.0, y)


CORE_DATA = CoreData()


def square_exact(x, y):
    return np.exp(x + y)


def square_forcing(x, y):
    return -2.0 * np.exp(x + y)


def sine_exact(x, y):
    return np.sin(np.pi * x) * np.sin(np.pi * y)


def sine_forcing(x, y):
    return 2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)


def harmonic_exact(x, y):
    """Smooth harmonic control on the L-shaped domain."""
    return np.exp(x) * np.cos(y)


def harmonic_gradient(x, y):
    return np.exp(x) * np.cos(y), -np.exp(x) * np.sin(y)


def singular_angle(x, y):
    """Return the project-prescribed branch angle theta and phi.

    theta=atan2(y,x), followed by theta += 2*pi whenever theta<pi/2.
    The added branch is essential on the positive-x re-entrant ray and on the
    analytic continuation to the outer boundary.  No gradient is evaluated at
    the re-entrant corner.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    theta = np.arctan2(y, x)
    theta = theta + np.where(theta < np.pi / 2.0, 2.0 * np.pi, 0.0)
    phi = theta - np.pi / 2.0
    return theta, phi


def singular_exact(x, y):
    """u=r^(2/3) sin(2 phi/3), including its specified angle branch."""
    x = np.asarray(x)
    y = np.asarray(y)
    r = np.hypot(x, y)
    _, phi = singular_angle(x, y)
    return r ** (2.0 / 3.0) * np.sin((2.0 / 3.0) * phi)


def singular_gradient(x, y):
    """Gradient away from r=0; raises if any evaluation point is the corner."""
    x = np.asarray(x)
    y = np.asarray(y)
    r = np.hypot(x, y)
    if np.any(r == 0.0):
        raise ValueError("The exact singular gradient is undefined at (0,0).")
    theta, phi = singular_angle(x, y)
    lam = 2.0 / 3.0
    radial = lam * r ** (lam - 1.0) * np.sin(lam * phi)
    angular = lam * r ** (lam - 1.0) * np.cos(lam * phi)
    erx, ery = np.cos(theta), np.sin(theta)
    etx, ety = -np.sin(theta), np.cos(theta)
    return radial * erx + angular * etx, radial * ery + angular * ety


def pure_neumann_data(x, y, side):
    """Outward-normal derivative of exp(x+y) on a unit-square side."""
    if side == "left":
        return -np.exp(x + y)
    if side == "right":
        return np.exp(x + y)
    if side == "bottom":
        return -np.exp(x + y)
    if side == "top":
        return np.exp(x + y)
    raise ValueError(f"Unknown side: {side}")
