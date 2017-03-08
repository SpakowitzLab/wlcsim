"""This module will be for reproducing Tom Lampo's results on MSDs of colloidal
particles in confinement. Should ask Phillip Geissler's student about in what
cases in found normal diffusivity."""

#TODO write whole module, lol

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

def bd_sim(nstep, dt, sphere_radii, confinement_radius):
    """Run a simulation of the spheres of sizes requested."""

