.. _features:

Fortran Simulation Features
###########################

.. Why use wlcsim? Does it do what I need?

Semiflexable polymers
=====================

The feature that wlcsim was originally designed around and for which it is named
is it's ability to quickly and accurately simulation polymers with stiffness
ranging from quite rigid to very flexible.

For very stiff polymers, the usual wormlike chain is simulated.

For relatively more flexible polymers, the "stretchable, shearable" chain is
used.

For *VERY* stretchable polymers, a purely Gaussian chain is used.

Of these, the middle strategy is where this code base stands out. The
"stretchable, shearable wormlike chain" is described in `"Discretizing elastic
chains for coarse-grained polymer models"
<https://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.013304>`_ by E. F.
Koslover and A. J. Spakowitz (2013). This model allows for quick but still
accurate simulation of polymer statistics of semiflexible polymers of a wide
range of stiffnesses (i.e. persistence lengths).

The effective potential for the stretchable, searable wormlike chain is given in
terms of constants. The bending modulus, the stretch modulus, the shear modulus,
the bend-shear coupling, and the ground state segment length. These constants
are tabulated for different descritizations in `input/dssWLCparams`. The module
`MC_wlc` handles interpolating from this table and calculating the potential
between two beads. See :ref:`MC_wlc`.

Field theoretic interactions
============================

Monte-Carlo simulations of solutions and melts of polymers are greatly
accelerated by the use of a field-theoretic, binned-density approach for
calculating interbead interactions. This approach is described in
"Theoretically informed coarse grain simulations of polymeric systems." by Pike
et. al. (2009). With a further description of our implementation given in
"Field-theoretic simulations of random copolymers with structural rigidity" by
Mao et. al. (2017).

Under this approach the local volume fraction of each constituent is calculated
for each of a grid of bins. The code for calculating the volume fractions and
updating them after a Monte-Carlo move are

.. f:autosrcfile:: mc_int.f03

The process of calculating the change in volume fractions after a move is often
the computational bottle-neck in the simulation. Routines for doing this that
are optimized for different move types are documented here :ref:`other_int_fun`.

The density is interpolated linearly between bins as described by Pike. This
interpolation is performed in

.. f:autosrcfile:: mc_interp.f03

To run a simulation with the field theoretic interaction turned on set
`WLC_P__FIELD_INT_ON` to true. You will also need to specify the relevant
interaction parameters to your problem, for example `WLC_P__CHI`. The details
of the interaction depend on problem (i.e. `WLC_P__FIELDINTERACTIONTYPE`). What
each of these interactions are (or to add your own) see
`src/mc/mc_hamiltonain.f03`.

.. f:autosrcfile:: mc_hamiltonian.f03

If you edit this file be sure to edit both the initial calculation of volume
fraction and change in volume fraction options.

As a final note, for liquid crystal systems the simulation has additional volume
fraction fields that are turned on by `WLC_P__CHI_L2_ABLE`. There are five such
fields and they correspond to the `m` values for sphericals harmonics of `l=2`
between -2 and 2.

Nucleosome geometry
====================

For simulations of chromatin that require the geometry of nucleosomes to
determine the entry exit angles of the DNA linkers the module `nucleosome` has
the details.

:ref:`nucleosome`

Streamlined Energy Components
=============================

The are many different energy contributions. To keep track of these (or to add
more) you should look to the `energies` module.

:ref:`energies`

