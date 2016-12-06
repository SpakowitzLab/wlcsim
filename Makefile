#
# Makefile
#
# wlcsim
#
# Author: Bruno Beltran <brunobeltran0@gmail.com>
#

# # replaced rm **.o with more portable find commands
# # and regex search =~ in distclean with manual check for Y/y
SHELL:=/bin/bash
# .SHELLFLAGS="-O extglob -c"

# script to automatically generate dependencies
MAKEDEPEND=./fort_depend.py

# $(DEP_FILE) is a .dep file generated by fort_depend.py
DEP_FILE = wlcsim.dep

# compiler
FC = mpifort
#TODO: fallback to gfortran gracefully, maybe with a dummy mpi.mod file?

# compile flags
INCLUDE_DIRS = -Isrc -Isrc/third_party/FLAP/exe/mod -Isrc/third_party
DEBUGFLAGS = -ggdb -Jsrc ${INCLUDE_DIRS} -cpp
FASTFLAGS = -O3 -Jsrc ${INCLUDE_DIRS} -cpp
PEDANTICFLAGS = -ggdb -Jsrc ${INCLUDE_DIRS} -cpp -fcheck=all -Wall -pedantic -fall-intrinsics -Wno-surprising # need instrincis because need sizeof, Wno-surprising to enforce Werror even though gfortran has a bug https://gcc.gnu.org/ml/fortran/2013-08/msg00050.html
FCFLAGS = ${FASTFLAGS}

# link flags
FLFLAGS =

# all non-legacy and non-test files should be compiled into wlcsim
SRC := $(shell find "src" -type f -name '*.f*' \
			    -not -path "src/legacy/*" \
				-not -path "src/tests/*" \
				-not -path "src/third_party/FLAP/*" \
				-not -path '*/\.*') # \
# -not -path 'src/wlcsim/wlcsim_brad.f03')
# -not -path 'src/wlcsim/wlcsim_quinn.f03' -not -path 'src/wlcsim/MCparrll_mpi.f90' -not -path 'src/wlcsim/restart_mpi.f90')

# takes each *.f* -> *.o
OBJ := $(addsuffix .o,$(basename $(SRC)))
TEST := $(shell find "src/tests" -type f -name '*.f*')

# program name
PROGRAM = wlcsim.exe

# by default, compile only
all: $(PROGRAM) flap depend

# a target to just run the main program
run: $(PROGRAM) dataclean
	./$(PROGRAM)

# target to build main program, depends on all object files, and on the
# constructed makefile output
$(PROGRAM): flap depend dummy_prog

# ugly line, needs to ask for FLAP objects at runtime, there's probably a better
# way to do this
dummy_prog: $(OBJ)
	$(FC) $(FCFLAGS) $(FLFLAGS) -o $@ $^ $(INCLUDE_DIRS) $(shell find "src/third_party/FLAP/exe/obj" -type f -not -path "src/third_party/FLAP/exe/obj/test_minimal.o")

# build third party dependencies that require "make" by hand
FLAP_DIR = src/third_party/FLAP
flap: flap_exists $(shell find ${FLAP_DIR} -type f -name '*.f*')
	make -C ${FLAP_DIR}

flap_exists: ${FLAP_DIR}
	git submodule update --init --recursive


# Make dependencies, easier to type
depend: $(DEP_FILE)

# this forces the "included" part of makefile itself to be rebuilt whenever the
# corresponding src files change, or when the Makefile itself changes
# While this is useful, it will force the dependency checker to run even if we
# only want to do e.g. make clean.
$(DEP_FILE): $(SRC) Makefile
	@echo "Making dependencies!"
	$(MAKEDEPEND) -w -o $(DEP_FILE) -f $(SRC) -c "$(FC) -c $(FCFLAGS) "

include $(DEP_FILE)

.PHONY: depend clean destroy dataclean

clean: dataclean
	find src \( -iname '*.o' -or -iname '*.mod' \) -delete
	rm -f ${PROGRAM} wlcsim.dep

dataclean:
	mkdir -p data trash
	touch "data/`date`"
	mv data/* trash/.

DEATH=rm -rf trash data savedata par-run-dir.*
distclean: clean
	@echo "About to destroy all simulation data, are you sure? ";
	read REPLY; \
	echo ""; \
	if [[ "$${REPLY:0:1}" == "Y" || "$${REPLY:0:1}" == "y" ]]; then \
		echo 'Running `${DEATH}`'; \
		${DEATH}; \
	else \
		echo 'Canceling data deletion!'; \
	fi

