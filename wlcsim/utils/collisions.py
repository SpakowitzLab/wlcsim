"""Functions for keeping track of collisions between particles
faster-than-naive-ly.

We've got bin_* algorithms for bin-based collision detection that I'll
implement first, and sort_* algorithms for an insertion-sort based collision
detection method that's fast if particles dont swap places very often."""
#TODO compare bin_* method with method where particle position array itself is
# binned in the original simulation code, for cache-help
import numpy as np
from numba import jit
import warnings
import unittest

def bin_init(particle_positions, max_particles_per_bin,
             num_bins_per_edge, bounds):
    """Returns a structure tracking which paritcles are in which bins.
    if max_particles_per_bin is violated, undefined behavior.
    bounds should be a tuple (len == ndims) of tuples (len == 2 ea) holding the
    min and max values for the coordinates of each dimension.
    num_bins_per_edge should be how many bins you'll have in each dimension."""
    ndims = len(bounds)
    num_particles, ndim_x = particle_positions.shape
    if ndims != ndim_x:
        raise ValueError("Number of dimensions of bounds provided {ndims:d} "
                         "does not match number of dimensions inferred from "
                         "coordinates of each particle position provided "
                         "{ndim_x:d}.".format(ndims=ndims, ndim_x=ndim+x))
    if ndims*num_bins_per_edge > num_particles:
        warnings.warn("More boxes than particles is inefficient! Use a "
                      "tree-based approach for finding which boxes are "
                      "currently in use if this is truly required.",
                      RuntimeWarning)
    # each bin corresponds to a list of integer indices to particles that are
    # in that bin and to a number of particles in that bin
    count_dims = [num_bins_per_edge for _ in range(ndims)]
    bin_dims = count_dims + [max_particles_per_bin]
    binned_particles = np.zeros(bin_dims, dtype=np.int)
    binned_counts = np.zeros(count_dims, dtype=np.int)
    # bin_walls = [np.linspace(min, max, num_bins_per_edge+1)
    #              for min, max in bounds]
    # transpose to help cacheing
    particle_bin_idxs = np.zeros((ndims, num_particles), dtype=np.int)
    xs = particle_positions.T
    for i, (min,max) in enumerate(bounds):
        for j, x in enumerate(xs):
            particle_bin_idxs[i,j] = binx(x, num_bins_per_edge, min, max)
    # untranspose to help caching
    particle_bin_idxs = particle_bin_idxs.T
    for i in range(num_particles):
        # since our box is a cube, with same number of discretizations on each
        # dimension, the multiidx->linear idx conversion is just a polynomial
        # evaluation, if we wanted to do it
        # bin_idx = np.polyval(particle_bin_idxs[i,:], num_bins_per_edge)
        c_idx = particle_bin_idxs[i,:]
        b_idx = c_idx + (binned_counts[c_idx],)
        binned_particles[b_idx] = i
        binned_counts[c_idx] += 1


def binx(x, nbins, min, max):
    """Our binning algorithm, separated out for the purposes of testing before
    vectorization."""
    return np.floor(nbins*(x - min)/(max - min)).astype(np.int)


class TestBinning(unittest.TestCase):
    """Suppose we have 4 bins, with boundaries 0,1,2,3,4.
    Then we know that the following x -> binx(x) should hold:
        -1 -> anything < 0
        0 -> 0
        0.5 -> 0
        1.0 -> 1
        4.0 -> anything > 3
    """
    test_min = 0
    test_max = 4
    nbins = 4
    bin_test = lambda x: binx(x, TestBinning.nbins, TestBinning.test_min, TestBinning.test_max)

    def test_binx_below_min(self):
        self.assertLess(TestBinning.bin_test(-1), 0)

    def test_binx_equal_min(self):
        self.assertEqual(TestBinning.bin_test(0), 0)

    def test_binx_equal_internal_bdry(self):
        self.assertEqual(TestBinning.bin_test(1), 1)

    def test_binx_equal_max(self):
        self.assertGreater(TestBinning.bin_test(4), 3)
