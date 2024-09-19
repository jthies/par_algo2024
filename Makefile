################################################################
# define how to compile Fortran 2018 code with coarray support #
################################################################

FC=caf
FFLAGS=-O3 -march=native -fopenmp -std=f2018 -cpp -ffree-line-length-none
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
all: main_hello.x main_dotprod.x main_benchmarks.f90

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

clean:
	-rm *.x *.o *.mod

test: main_tests.x
	srun -n 4 --mem-per-cpu=1G ./main_tests.x
