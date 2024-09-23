################################################################
# define how to compile Fortran 2018 code with coarray support #
################################################################

FC=caf
# use these flags for optimized build
FFLAGS=-O3 -march=native -fopenmp -std=f2018 -cpp -ffree-line-length-none
# and these for developing/debugging
FFLAGS=-O0 -g -fopenmp -std=f2018 -cpp -ffree-line-length-none -fcheck=all
LDFLAGS=

################################################################
# Define where to find the fortuno_coaray installation.        #
################################################################

FORTUNO_INCLUDE=-I/beegfs/apps/unsupported/fortuno/include
FORTUNO_LDFLAGS=-L/beegfs/apps/unsupported/fortuno/lib -lfortuno-coarray

######################################
# no system-specific settings below! #
######################################

# by default, build all examples:
all: main_hello.x main_dotprod.x main_benchmarks.f90 main_sorting.x

###############################################
# these should always be executed when called #
###############################################

.PHONY: clean test

##############################
# general rules              #
##############################

%.o: %.f90 Makefile
	${FC} ${FFLAGS} -c $<

%.x: %.f90 
	${FC} ${FFLAGS} -o $@ $^ ${LDFLAGS}

# note: for this target we first have to compile dotprod.f90 to produce dotprod.o and m_dotprod.mod
main_dotprod.x: main_dotprod.f90 dotprod.o

# this driver needs the m_benchmarks.mod module, which is produced when compiling benchmarks.f90
main_benchmarks.x: main_benchmarks.f90 benchmarks.o
	${FC} ${FFLAGS} -o $@ $^ ${LDFLAGS}

main_sorting.x: main_sorting.f90 benchmarks.o sorting.o
	${FC} ${FFLAGS} -o $@ $^ ${LDFLAGS}

# this module requires the Fortuno test framework with coarray support (fortuno-coarray)
unit_tests.o: unit_tests.f90 dotprod.o sorting.o Makefile
	${FC} ${FFLAGS} -c $< ${FORTUNO_INCLUDE}

main_tests.x: main_tests.f90 unit_tests.o benchmarks.o dotprod.o sorting.o
	${FC} ${FFLAGS} ${FORTUNO_INCLUDE} -o $@ $^ ${LDFLAGS} ${FORTUNO_LDFLAGS}

test: main_tests.x
	mpirun -np 1 ./main_tests.x
	mpirun -np 2 ./main_tests.x
	mpirun -np 3 ./main_tests.x
	mpirun -np 4 ./main_tests.x


clean:
	-rm *.x *.o *.mod

