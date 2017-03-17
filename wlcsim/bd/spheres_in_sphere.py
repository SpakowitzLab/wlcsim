"""This module will be for reproducing Tom Lampo's results on MSDs of colloidal
particles in confinement. Should ask Phillip Geissler's student about in what
cases in found normal diffusivity."""
import math
from ..mc import sphere_insertion

max_sphere_density = math.pi/(3*math.sqrt(2))
max_circle_density = math.pi/(2*math.sqrt(3))

def confinement_force(sphere1, containing_sphere):
    """Returns the force exerted on sphere1 by its proximity to the confinement
    sphere."""
    pass

def space_exclusion_force(sphere1, sphere2):
    """Returns the force exerted on each of two spheres due to their
    interactions with each other.

    After reading a lot about how people work with confinement, I've got the
    intuition that it's standard to arbitrarily choose among the following
    depending on what makes your particular problem more tractable:

    1) full lennard-jones
    2) lennard-jones, but only the (a/r)^12 bit without the attractive bit
    3) something "slightly harder" than gaussian, in order to increase time
    step, e.g.
    \frac{A}{3 b \delta^2}(h-\delta)^3, where \delta = b*(L/N)^{1/2}/2, A
    arbitrarily ~25k_BT"""
    pass

def bd_step(dt, spheres, containing_sphere):
    """Moves the spheres forward a time step.

    TODO: Check if returning the list of new spheres is faster."""
    pass

def bin_from_scratch():
    #TODO use histogramdd
    pass

def bd_sim(nstep, dt, target_density, sphere_radii, confinement_radius,
           bins_per_side, num_init_steps=10000):
    """Run a simulation of the spheres of sizes requested."""
    # the bins will tile the bounding box of the confining sphere, so the overall
    # bounding box has
    side_length = 2*confinement_radius
    # and each little bin has a side_length given by
    bin_length = side_length/bins_per_side
    bin_volume = bin_length**3
    # max_spheres_per_bin*volume_of_sphere = bin_volume*max_sphere_density
    sphere_volume = (4/3)*math.pi*(sphere_radii**3)
    max_spheres_in_bin = math.ceil(max_sphere_density*bin_volume/sphere_volume)
    num_bins = bins_per_side**3
    sphere_centers = np.zeros((num_bins, max_spheres_in_bin, 3))
    num_in_bin = np.zeros((num_bins,))
    # run num_init_steps MC steps to get spheres in reasonable positions
    init_spheres = mc.sphere_insertion.step_away_mc(num_init_steps, target_density, sphere_radii, confinement_radius)
    #TODO put spheres into bins
    for i in range(nstep):
        #TODO ipmlement bin-based algorithm for fast simulation

