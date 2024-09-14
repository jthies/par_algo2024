################################################################
# define how to compile Fortran 2018 code with coarray support #
################################################################

FC=gfortran
FFLAGS=-O3 -march=native -fcoarray=lib -fopenmp -std=f2018 -cpp -ffree-line-length-none
LDFLAGS=-L${OPENCOARRAYS_ROOT}/lib -lcaf_mpi

################################################################
# Define where to find the fortuno_coaray installation.        #
################################################################

FORTUNO_INCLUDE=-I/beegfs/apps/unsupported/fortuno/include
FORTUNO_LDFLAGS=-L/beegfs/apps/unsupported/fortuno/lib -lfortuno-coarray

######################################
# no system-specific settings below! #
######################################

# by default, build all examples:
all: main_hello.x main_dotprod.x main_benchmarks.f08

##############################
# general rules              #
##############################

%.o: %.f08 Makefile
	${FC} ${FFLAGS} -c $<

%.x: %.f08 
	${FC} ${FFLAGS} -o $@ $^ ${LDFLAGS}

# note: for this target we first have to compile dotprod.f08 to produce dotprod.o and m_dotprod.mod
main_dotprod.x: main_dotprod.f08 dotprod.o

# this driver needs the m_benchmarks.mod module, which is produced when compiling benchmarks.f08
main_benchmarks.x: main_benchmarks.f08 benchmarks.o
	${FC} ${FFLAGS} -o $@ $^ ${LDFLAGS}

clean:
	-rm *.x *.o *.mod

test: main_tests.x
	srun -n 4 --mem-per-cpu=1G ./main_tests.x
