import numpy as np

from functools import reduce


def center_by_mass(x, particle_axis=0):
    """Subtract center of mass (unweighted) from a collection of vectors
    corresponding to particles's coordinates."""
    shape = x.shape
    # center of mass is the average position of the particles, so average over
    # the particles
    centers_of_mass = np.mean(x, particle_axis)
    # tile the centers of mass for clean elementwise division
    # is there a way to do this faster by relying on numpy's array projection
    # semantics?
    tile_shape = [shape[i] if i == particle_axis else 1
                  for i in range(len(shape))]
    centers_of_mass = np.tile(centers_of_mass, tile_shape)
    return x - centers_of_mass

from scipy import arange, pi, sqrt, zeros
from scipy.special import spherical_jn
from scipy.optimize import brentq
def Jn(r,n):
    return spherical_jn(n, r)

def spherical_jn_zeros(n, nt):
    """recursive method: computes zeros ranges of sphereical_jn(r,n) from zeros
    of sphereical_jj(r,n-1)

    pros : you are certain to find the right zeros values;
    cons : all zeros of the n-1 previous Jn have to be computed;
    note : Jn(r,0) = sin(r)/r
    """
    zerosj = zeros((n+1, nt), dtype=float)
    zerosj[0] = arange(1,nt+1)*pi
    points = arange(1,nt+n+1)*pi
    racines = zeros(nt+n, dtype=float)
    for i in range(1,n+1):
        for j in range(nt+n-i):
            foo = brentq(Jn, points[j], points[j+1], (i,))
            racines[j] = foo
        points = racines
        zerosj[i][:nt] = racines[:nt]
    return (zerosj)
