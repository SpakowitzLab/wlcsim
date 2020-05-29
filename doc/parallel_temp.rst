.. _parallel_temp:

About Parallel Tempering
########################

This code is set up to allow parallel tempering over multiple variables if desired.  The replicas will be along a path in parameter space parameterized by `s` which goes between 0 and `WLC_P__INITIAL_MAX_S`.

To determine which variables will be parallel tempered over set the corresponding switches.  For example, to parallel temper over chi set `WLC_P__PT_CHI` to true.  By default, the path is set to linearly go between 0 and `WLC_P__CHI` in this case.  However, you can define a different path in the function `cof_path_by_energy_type` in the file `src/wlcsim/wlcsim_quinn`.  For example, MU starts from -2.5 rather than zero.

How to run with parallel tempering
==================================

1. Install MPI fortran if not already installed.
2. In defines file, set `WLC_P__PTON=.TRUE.`, set the variable(s) you want to temper over, and set `WLC_P__CODENAME="quinn"`.
3. In file `Makefile` set ‘FC=mpifort`
4. Compile using `make`.
5. Run using `mpirun -np 16`.  This will run with 15 replicas and 1 head node.  You may need to provide the path to mpirun on your computer, mine is located at `usr/bin/mpirun`.

Adjustable spacing
==================

By default, replicas are spaced evenly between 0 and the maximum value.  You can adjust the spacing by writing a more complicated function into cof_path_by_energy_type` in the file `src/wlcsim/wlcsim_quinn`.  There is a built in way to adjust the spacing of the replicas to achieve an acceptable exchange rate between sequential replicas.  To do this set `WLC_P__INDSTARTREPADAT`, and `WLC_P__INDENDREPADAPT` to the save point range during which you would like the adaptation to occur.  For more details see the file `src/mc/adaptCof`.

