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
    00. Make sure that the following variable name changes don't break the code
        (should be able to grep for most as mc%${var_name} or md%${var_name}
        settype -> initCondType
        repSufix -> repSuffix
        G -> nBpM ??Quinn
        N -> nMpP (in quinn's code)
        N -> nB (in Andy's code)
        V -> beadVolume
        h_A -> hA
        ? -> tweak_param_defaults(wlcsim_p, wlcsim_d)
            => but instead, seems like MCAMP and SUCCESS and PHIT and WINDOW should
            be in wlcsim_p instead of _d, because they're being tweaked as parameters,
            but later on their used as proxy for current simulation state, so which
            are they?
        simType -> solType
        NULL -> simType ! wlc, sswlc, or gc
        NULL -> codeAuth ! whose code to run
        inton -> field_interactions
        self interactions -> intrapolymer_stick_crossing_enforced
        col_type -> fptColType
        ind -> now passed from wlcsim
        wlcsim_params_saveparameters -> save_parameters
        wlcsim_params_printEnergies -> printEnergies
        time_ind now refers to what time step you're on
        step_ind now refers to what mc step you're on


        notes to Quinn:
        in places where you write out simulation state, writing out paramters as
        well seems unecessary, but no need to 'fix':
        e.g. wlcsim_params_appendAdaptData
    1. Make suite of tests, with contributions from Andy's personal
       codebase (ask when you get here)
    2. attempt to compile, see if any bugs
    3. fix miscellaneous discrepancies that might turn into bugs
        00. check for TODOs in code
        0. for MC, brad uses energy self chain instead of energy ponp
        1. change all uses of mt19937 to mersenne_twister
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
        13. document that adapt.f90 and adaptCof.f90 are the "crux" of optimizing the MC
        15. Quinn should probably move shedule.f95 into wlcsim_quinn.f95, or similar
        16. Quinn should probably combine MCparrll_mpi.f90 and wlcsim_quinn.f95, to
            create a single function wlcsim_quinn, that wlcsim.f95 can call.
        17. Brad should probably change wlcsim_brad.f95 so that wlcsim.f95 can call it.
        18. Add hook to makefile to perform a recursive git submodule udpate if third
            party libraries not present, using ifeq directive, maybe checking
            if some $(wildcard ...) expression evaluates to the empty string?
        19. merge quinn and andy's move settings in params.f03
        20. guard writing via appendAdaptData and appendEnergyData to
            when the relevant variables will exist in
            save_simulation_state
        21. fix *_saveparameters functions to write out structures
            correctly
        22. when should repeatingBC be on?
        23. quinn, where did you get that seeding code?
        24. shouldn't use wlc_d%ind anymore
        25. we no longer have decom, from MCsim, should it be there? if so,
            make sure lhc (== fcom) has correct magnitude. seems unlikely that
            it's needed since it comes from an outdated input file type.
            appears to be used in force_ponp, and mc_self, maybe it has a new
            name now?
        26. same as lhc but for vhc
        27. add brad's temporary file things to .gitignore. looks like *~, ._*,
            #*# would do it
        28. make sure time_ind is incremented in bdsim
        29. make sure head node generates random seed sfor all workers
        30. change all case's in params to CAPS_WITH_UNDERSCORES
        31. only run VerifyEnegiesFromScratch if not manually changing mc
            parameters (Quinn)
        32. in CalculateEnergiesFromScratch, do Brad's eknot, self, etc
        33. finish code to initialize energies in params and to get out x_'s
            from CalculateEnergiesFromScratch

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
