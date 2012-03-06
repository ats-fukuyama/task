include ../../task/make.header


#FFLAGS=$(OFLAGS) -I$(INCLUDE_MUMPS) -I../mtxp
FFLAGS=$(DFLAGS) -I$(INCLUDE_MUMPS) -I../mtxp -I../../task/lib/mod

all: $(EXEC)

##### Without PETSc library #####
#FLINKER=$(FCFREE)
#INCLUDE_MUMPS=.
#EXEC=testbndc
#.SUFFIXES:
#.SUFFIXES: .f .f90 .mod .o .a
#.f90.o :
#	$(FCFREE) $(FFLAGS) -c $< -o $@
#clean :
#	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod

##### With PETSc library 3.0 (complex) #####
#include ${PETSC_DIR}/conf/base
#INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS_4.8.4/include
#INCLUDE_PETSC=${PETSC_DIR}/include
#EXEC=testbndc testkspc testzmumps

##### With PETSc library 3.2 (complex) #####
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS_4.10.0/include
INCLUDE_PETSC=${PETSC_DIR}/include
EXEC=testbndc testkspc testzmumps

LIBMPI=../mtxp/libmtxmpi.o
LIBNOMPI=../mtxp/libmtxnompi.o
LIBKSPC=$(LIBMPI) libmtxkspc.o
LIBZMUMPS=$(LIBMPI) libmtxzmumps.o
LIBBNDC=$(LIBNOMPI) libmtxbndc.o

all: $(EXEC)

libs: $(LIBNOMPI) $(LIBMPI) libmtxbndc.o libmtxkspc.o libmtxzmumps.o

testkspc: $(LIBKSPC) testmtxc.o
	-${FCFREE} -g -o testkspc testmtxc.o $(LIBKSPC) ${PETSC_KSP_LIB} ../lib/tasklib.a $(FLIBS)

testzmumps: $(LIBZMUMPS) testmtxc.o
	-${FCFREE} -g -o testzmumps testmtxc.o $(LIBZMUMPS) ${PETSC_KSP_LIB} ${PETSC_DIR}/${PETSC_ARCH}/lib/libzmumps.a ../lib/tasklib.a $(FLIBS)

testbndc: $(LIBBNDC) testmtxc.o 
	-${FCFREE} -g -o testbndc testmtxc.o $(LIBBNDC) ../lib/tasklib.a $(FLIBS) 

testmtxc.o: testmtxc.f90
../mtxp/libmtxmpi.o: ../mtxp/libmtxmpi.F90
../mtxp/libmtxnompi.o: ../mtxp/libmtxnompi.f90
libmtxbndc.o: libmtxbndc.f90
libmtxzmumps.o: libmtxzmumps.F90
libmtxkspc.o: libmtxkspc.F90
