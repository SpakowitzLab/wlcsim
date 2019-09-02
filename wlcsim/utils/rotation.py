import numpy as np
from numba import jit

@jit(nopython=True)
def cot_plane(u, unit=True, tol=1e-8):
    """Return an arbitrary basis to the cotangent plane given the tangent
    vector (assumed unit vector unless otherwise specified)."""
    if not unit:
        u = u/np.linalg.norm(u)
    n1 = np.array([0.0, 0.0, 1.0])
    n1 = n1 - (n1@u)*u
    mag = np.linalg.norm(n1)
    # if we accidentally chose vector parallel to u
    if np.abs(mag) < tol:
        n1 = np.array([0.0, 1.0, 0.0])
        n1 = n1 - (n1@u)*u
        mag = np.linalg.norm(n1)
    n1 = n1/mag
    # n2 = np.cross(n1, u)
    # for jit
    n2 = np.array([n1[1]*u[2] - n1[2]*u[1],
                   n1[2]*u[0] - n1[0]*u[2],
                   n1[0]*u[1] - n1[1]*u[0]])

    n2 = n2/np.linalg.norm(n2)
    return n1, n2

def Ry(theta):
    r"""Rotation matrix about the z axis with angle theta.

    Notes
    -----
    ..math::

        \frac{1}{\sqrt{3}}
        \begin{bmatrix}
            np.cos(theta) & 0 & -np.sin(theta) \\
                        0 & 1 &              0 \\
            np.sin(theta) & 0 &  np.cos(theta)
        \end{bmatrix}
    """
    return np.array([[np.cos(theta), 0, -np.sin(theta)],
                     [            0, 1,              0],
                     [np.sin(theta), 0,  np.cos(theta)]])


def Rz(theta):
    r"""Rotation matrix about the z axis with angle theta.

    Notes
    -----
    ..math::

        \frac{1}{\sqrt{3}}
        \begin{bmatrix}
            np.cos(theta) & -np.sin(theta) & 0 \\
            np.sin(theta) &  np.cos(theta) & 0 \\
                        0 &             0  & 1
        \end{bmatrix}
    """
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta),  np.cos(theta), 0],
                    [            0,             0,  1]])

