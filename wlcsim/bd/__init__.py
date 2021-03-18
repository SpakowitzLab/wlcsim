r"""
Modular Brownian Dynamics simulations of various types of polymers.

A Brownian dynamics simulation requires specification of three things:

    1. **The system state**. For us, this is :math:`\vec{r}_j(t_i) \in
       \mathbb{r}^{N\times 3}`, the positions of the :math:`N` discrete beads
       representing the polymer at each time :math:`t_i`.
    2. **The forces**. Functions taking as input the system state and returning
       as output the forces on each bead. The main forces in any polymer
       simulation are the chain forces, defined in `wlcsim.bd.rouse` and
       `wlcsim.bd.twlc` for their respective chains. For defining forces
       external to the polymer, there are many helpful functions in
       `wlcsim.bd.forces` covering confinement, spring-like bead tethers, and
       others.
    3. **The stochastic integrator**. The "main loop" of the simulation, this
       takes as input the initial system state and a force function. It returns
       as output the system state at some later times.

Using this module, running a simulation simply involves defining an initial
state, a force function, and then passing these to one of the provided
stochastic integrators. Either:

    1. `wlcsim.runge_kutta.rk4_thermal_lena`, which is a modified
       Euler-Maruyama method that instead uses the well-known deterministic
       fouth order Runge-Kutta scheme to compute the force at each time point.
       This method should work well for stiffer systems, where the
       deterministic forces dominate the stochastic forces.
    2. `wlcsim.runge_kutta.srk1_roberts`, which is a stochastic Runge-Kutta
       method of strong order 1. It approximates both the stochastic and
       deterministic forces to first order in time, and so will be faster for
       systems that are not very stiff, such as well-spaced beads on a Rouse
       polymer.

.. note::

    For more details on deterministic integrators, my favorite introductory
    resource is the book by Iserles, "A first course in the numerical analysis
    of differential equations". For more examples of stochastic integrators, a
    course on Ito calculus is required, after which the recent work by Andreas
    Rößler is an excellent starting point.


Example
-------

To simulate a Rouse chain in a soft, spherical confinement, for example, we can
simply combine the existing force functions for the Rouse chain and for a soft
elliptical confinement, and use the ``srk1`` integrator:

.. code-block:: python

    # we will use a "soft" confinement
    from wlcsim.bd.forces import f_conf_ellipse
    a = 1  # confinement radiu
    a_ex = 1  # confinement strength

    # we will use a Rouse chain
    from wlcsim.bd.forces import f_elas_linear_rouse
    # see wlcsim.bd.rouse.with_integrator for better parameterization
    k_over_xi = 3*Dhat/bhat**2

    def confined_rouse(x, t):
        '''The total force on each bead at each time step.'''
        return f_elas_linear_rouse(x, k) + f_conf_ellipse(x, a_ex, a, a, a)

    # start with all beads at the origin
    N = 100  # number of beads
    x0 = np.zeros((N, 3))

    # the better integrator for Rouse chains
    from wlcsim.bd.runge_kutta import srk1_roberts
    Nt = 1e5
    t = np.linspace(0, 1, Nt)
    # Use unit diffusivity. To compute diffusivity for Rouse chains, see
    # discussion in wlcsim.bd.rouse.measured_D_to_rouse.
    X = srk1_roberts(confined_rouse, D=1, t, x0)

After the above, ``X`` will hold the positions of each of the ``N`` beads at
each simulated time ``t`` (i.e. it will have shape ``(Nt, N, 3)``). For a more
realistic example, as used to generate data for publication, see the
`wlcsim.bd.homolog` documentation.
"""

from . import rouse
from . import runge_kutta
