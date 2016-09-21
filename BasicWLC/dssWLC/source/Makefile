OBJS = chainutil.o quatutil.o keys.o readkey.o genutil.o mt19937.o inputparams.o montecarlo.o inpututil.o browndyn.o redisc.o cylinder.o dbrent.o nrtype.o nrutils.o manychains.o dgtsl.o sampleutil.o besselI0.o

FC = gfortran
FLAGS = -O3 
DEBUGFLAGS =
LIBS = -lblas -llapack

# executable name
EXE = ../shearableWLC.exe
TESTEXE = ../testrun.exe

.SUFFIXES: .o .f90 

# general suffix rules
%.o: %.f90
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

default: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(EXE) $(LIBS)

test: $(OBJS) testmain.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testmain.o -o $(TESTEXE) $(LIBS)

clean: 
	rm *.o *.mod

veryclean: clean
	rm *~ $(EXE)

# ------------
# individual file dependencies
# ------------
chainutil.o: keys.o quatutil.o mt19937.o genutil.o
readkey.o: inputparams.o keys.o genutil.o mt19937.o
genutil.o: mt19937.o quatutil.o
browndyn.o: chainutil.o genutil.o dgtsl.o
redisc.o: chainutil.o cylinder.o
cylinder.o: dbrent.o 
dbrent.o: nrtype.o nrutils.o
montecarlo.o: chainutil.o redisc.o manychains.o browndyn.o
manychains.o: chainutil.o
sampleutil.o: dbrent.o besselI0.o