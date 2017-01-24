import numpy as np


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

