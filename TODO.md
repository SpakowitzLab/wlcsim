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
        N -> nMpP
        V -> beadVolume
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

        notes to Quinn:
        in places where you write out simulation state, writing out paramters as
        well seems unecessary, but no need to 'fix': 
        e.g. wlcsim_params_appendAdaptData
    1. Make suite of tests, with contributions from Andy's personal
       codebase (ask when you get here)
    2. attempt to compile, see if any bugs
    3. fix miscellaneous discrepancies that might turn into bugs
        0. for MC, brad uses energy self chain instead of energy ponp
        1. change all uses of mt19973 to mersenne_twister
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
        13. document that adapt.f90 and adaptCof.f90 are the "crux" of optimizing the MC
        15. Quinn should probably move shedule.f95 into wlcsim_quinn.f95, or similar
        16. Quinn should probably combine MCparrll_mpi.f90 and wlcsim_quinn.f95, to
            create a single function wlcsim_quinn, that wlcsim.f95 can call.
        17. Brad should probably change wlcsim_brad.f95 so that wlcsim.f95 can call it.
        18. Add hook to makefile to perform a recursive git submodule udpate if third
            party libraries not present, using ifeq directive, maybe checking
            if some $(wildcard ...) expression evaluates to the empty string?
    4. merge input styles of all code
        1. modify wlcsim_bruno.f03 to use new calling convention
        2. remove "dependent variables check" from quinn's params code in
           favor of always specifyign the same paramters and usign
           get_derived_parameters by default
        2. modify get_derived_parameters to take wlc_p


