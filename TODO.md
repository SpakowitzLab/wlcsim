# TODO

0. Use SPAG to beautify code
    0. agree on SPAG conventions - DONE!
        for linting:
        4 spaces, google style, f95, no more labels when possible.
        caps:
        variables: camelcase, nameThemLikeThis
        subroutines: caps camelcase, preferably starting with a verb NameThemLikeThis
        modules: caps camelcase
        for functions:
        intent(in/out), use _dp from setPrecision, implicit none,
        delete all unused variables
    1. run SPAG - TODO
    2. double check output for sanity
1. Elena's dssWLC should have it's own git repo, be a submodule
    0. make git repo, convert dssWLC dir to submodule
    1. figure out how to "make" parameters from Elena's code
    2. add make target in Makefile to make dssWLC parameters
2. Merge Quinn + Brad's codebases per the following notes.
Stack of things to do:
    0. Profit!
3. Make FASTFLAGS actually fast flags instead of just O3.
4. Make constant parameters compile time constants.
5. Implement a binning based u(r) two body interaction in Monte-Carlo.
6. Implement an umbrella sampling procedure.
7. update files that use eps/epsapprox to use precision.mod


sitting in the src/wlcsim directory, and you've just run make -i in the top
level, the following now succeeds if you remove
references to mcsim from wlcsim_bruno:
gfortran wlcsim.f03 -I .. -I ../third_party/FLAP/exe/mod \
../third_party/FLAP/exe/obj/flap* ../third_party/FLAP/exe/obj/penf* \
../util/stop_if_err.o params.o ../util/inputparams.o wlcsim_bruno.o \
../third_party/mersenne_twister.o initcond.o get_derived_parameters.o \
../bd/BDsim.o ../bd/concalc.o ../bd/force_ponp.o ../bd/force_elas.o \
../util/colchecker.o ../util/colsort.o ../third_party/kdtree2.o ../bd/stress.o \
../bd/stressp.o ../third_party/mt19937.o ../bd/RKstep.o ../third_party/dgtsv.o


EELAS (1) = wlc
EELAS (2) = wlc
EELAS (3) = wlc
EELAS (4) = twist
