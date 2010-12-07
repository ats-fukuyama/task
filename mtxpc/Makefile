include ../../task/make.header


#FFLAGS=$(OFLAGS) -I$(INCLUDE_MUMPS)
FFLAGS=$(DFLAGS) -I$(INCLUDE_MUMPS)

all: $(EXEC)

##### Without PETSc library #####
#FLINKER=$(FCFREE)
#INCLUDE_MUMPS=.
#EXEC=testbnd
#.SUFFIXES:
#.SUFFIXES: .f .f90 .mod .o .a
#.f90.o :
#	$(FCFREE) $(FFLAGS) -c $< -o $@
#clean :
#	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod

##### With PETSc library(complex) #####
include ${PETSC_DIR}/conf/base
INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS_4.8.4/include
INCLUDE_PETSC=${PETSC_DIR}/include
EXEC=testbndc testkspc testzmumps testmtxp 

LIBMPI=../mtxp/libmtxmpi.o
LIBNOMPI=../mtxp/libmtxnompi.o
LIBKSPC=$(LIBMPI) libmtxkspc.o
LIBZMUMPS=$(LIBMPI) libmtxzmumps.o
LIBBNDC=$(LIBNOMPI) libmtxbndc.o


all: $(EXEC)

libs: $(LIBNOMPI) $(LIBMPI) libmtxbndc.o libmtxkspc.o libmtxzmumps.o

testkspc: $(LIBKSPC) testmtxc.o
	-${FLINKER} -g -o testkspc testmtxc.o $(LIBKSPC) ${PETSC_KSP_LIB} $(GFLIBS) ../lib/tasklib.a

testzmumps: $(LIBZMUMPS) testmtxc.o
	-${FLINKER} -g -o testzmumps testmtxc.o $(LIBZMUMPS) ${PETSC_KSP_LIB} ${PETSC_DIR}/intel111/lib/libzmumps.a $(GFLIBS) ../lib/tasklib.a

testbndc: $(LIBBND) testmtxc.o 
	-${FLINKER} -g -o testbndc testmtxc.o $(LIBBNDC) ../../task/lib/libbnd.o

testmtxc.o: testmtxc.f90
../mtxp/libmtxmpi.o: ../mtxp/libmtxmpi.F90
../mtxp/libmtxnompi.o: ../mtxp/libmtxnompi.f90
libmtxbndc.o: libmtxbndc.f90
libmtxzmumps.o: libmtxzmumps.f90
libmtxkspc.o: libmtxkspc.F90
