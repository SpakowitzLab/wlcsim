SHELL:=/bin/bash
.SHELLFLAGS="-O extglob -c"
# compiler
FC = gfortran

# compile flags
FCFLAGS = -g -Jcode -cpp -Wall -pedantic

# link flags
FLFLAGS =

# source files, mt19937 first since everyone depends on it
# the `find` lists all fortran files in subdir 'code' except for mt19937 and
# those in the subdirectory code/MISCcode
SRCS := code/SIMcode/mt19937.f90 $(shell find code -type f \( -iname '*.f*' ! -name 'mt19937.f90' \) -not -path 'code/MISCcode/*')
#OBJS := $(addsuffix .o, $(basename $(SRCS)))
#HDRS := $(shell find code -name '*.h')
#MODS := $(patsubst %.h, %.mod, $(HDRS))

# program name
PROGRAM = wlcsim

# test:
# 	@echo $(value SRCS)

# by default, compile only
all: $(PROGRAM) Makefile

# a target to just run the main program
run: $(PROGRAM) dataclean
	$(PROGRAM)

# target to build main program
$(PROGRAM): $(SRCS)
	$(FC) $(FCFLAGS) $(FLFLAGS) -o $@ $^

# targets to build each object file
%.o: %.F
	$(FC) $(FCFLAGS) -o $@ $<
%.o: %.F77
	$(FC) $(FCFLAGS) -o $@ $<
%.o: %.F90
	$(FC) $(FCFLAGS) -o $@ $<
%.o: %.F95
	$(FC) $(FCFLAGS) -o $@ $<
%.o: %.f77
	$(FC) $(FCFLAGS) -o $@ $<
%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<
%.o: %.f95
	$(FC) $(FCFLAGS) -o $@ $<

.PHONY: clean destroy

clean: dataclean
	rm -f code/**.o code/**.mod code/**.MOD wlcsim

dataclean:
	mkdir -p data trash
	touch "data/`date`"
	mv data/* trash/.

destroy: clean
	rm -rf trash/* data/*
