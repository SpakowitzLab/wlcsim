# TODO

0) Use SPAG to beautify code
    0. agree on SPAG conventions - DONE!
    1. run SPAG - TODO
    2. double check output for sanity
1) Merge Quinn + Brad's codebases per the following notes.
merge input styles
put all parameters into a structure
combine code by cat'ing
MPI stuff (argggh?)

Quinn's simMod for structure of params, add Brad's and Bruno's parameters, remove
"dependent variables check", and split into two structures, one which should be
intent(in) everywhere except for in parallel tempering, and one which is changed
everywhere. All "energy constants" should be in the first structure, and the
"things the energy is multiply by" should be in the changing structure.

Bruno's getpara

in principle, we want
1) getInput - using quinn's param structure
2) validateInput - using quinn's idiot checks
3) getRenormalizedParamsFromInput - starting from Bruno's getpara

for linting:
4 spaces, google style, f95, no more labels when possible.
caps:
variables: camelcase, nameThemLikeThis
subroutines: caps camelcase, preferably starting with a verb NameThemLikeThis
modules: caps camelcase

for functions:
intent(in/out), use _dp from setPrecision, implicit none,
delete all unused variables

for documentation:
sphinx



CODE STUFF:

for global code, brad uses wlcsim_test
for MC, brad uses energy self chain instead of energy ponp

master code:
    call 1), 2), 3), then calls different wlcsim versions...

# define variable, mpi_id

para_in, para_inout = getpara....

mc_init() or load_from_file(para_inout)

if (simtype == 'replica couple integer parameters'):
    for i in range(numTimePoints):
        wlcsim(i, para_in, para_inout, output)
        #get_structure_quantities_brad(output, [1, 3, 4, 6])
        save_output(output)

elseif (simtype == 'bruno'):
    for i in range(numTimePoints):
        other_wlcsim(i, para, output)
        save_output(output)

OUTPUT:
r, u

REORGANIZATION:
1) elena's dssWLC should get own git repo
2) SIMcode should just have wlcsim, getpara, initcond
3) MISCcode purge
4) move collision code to misccode


STACK:
1. change mt199764123o49pu12r for mersenne_twister
2. finish contributors scratch information
3. add input: setType for initcond and verify parameters are passed correctly
4. change makefile to have FC=mpifort, and scan_wlcsim.py to run with mpirun if
   necessary
5. add number of mpi processes per simulation to use into scan_wlcsim.py
6. combine confinement branch with MC_confine.f95
6.5. use Brad's MC_self....
7. change call to MC_eelas
8. change call to MC_move
9. check set_type is set to 6 for rings
10. check set_type matches boundary conditions
11. inton -> field_interactions, self interactions -> INTERP_BEAD_LENNARD_JONES
12. ASK ANDY why both sides of IB1/IB2 are used to calculate energy change for
    WLC code in MC_elas. if IB1==IB2, isn't that double counting....
13. document that adapt.f90 and adaptCof.f90 are teh "crux" of optimizing the MC


