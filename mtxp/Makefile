include ../../task/make.header


#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

all: $(EXEC)

##### Without PETSc library #####
#FLINKER=$(FCFREE)
#INCLUDE_MUMPS=.
#EXEC=testbnd testbndc
#.SUFFIXES:
#.SUFFIXES: .f .f90 .mod .o .a
#.f90.o :
#	$(FCFREE) $(FFLAGS) -c $< -o $@
#clean :
#	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod

##### With PETSc library 3.0 (real) #####
#include ${PETSC_DIR}/conf/base
#INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS_4.8.4/include
#INCLUDE_PETSC=${PETSC_DIR}/include
#EXEC=testbnd testbndc testksp testdmumps testzmumps

##### With PETSc library 3.2 (real) #####
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS_4.10.0/include
INCLUDE_PETSC=${PETSC_DIR}/include
EXEC=testbnd testbndc testksp testdmumps testzmumps

INCMOD=-I../lib/mod
LIBMOD=-L../lib -l libgrf
LIBMPI=libmtxmpi.o
LIBNOMPI=libmtxnompi.o
LIBKSP=$(LIBMPI) libmtxksp.o
LIBDMUMPS=$(LIBMPI) libmtxdmumps.o
LIBZMUMPS=$(LIBMPI) libmtxzmumps.o
LIBBND=$(LIBNOMPI) libmtxbnd.o
LIBBNDC=$(LIBNOMPI) libmtxbndc.o
LIBDMUMPSX=-L${PETSC_DIR}/${PETSC_ARCH}/lib -ldmumps -lmumps_common -lpord -lscalapack -lblacs -lblas -lparmetis -lmetis
LIBZMUMPSX=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lblas -lparmetis -lmetis

all: $(EXEC)

libs: $(LIBNOMPI) $(LIBMPI) libmtxbnd.o libmtxksp.o libmtxdmumps.o

testksp: $(LIBKSP) testmtx.o
	-${FLINKER} -g -o testksp testmtx.o $(LIBKSP) ${PETSC_KSP_LIB}

testdmumps: $(LIBDMUMPS) testmtx.o
	-${MF90} -g -o testdmumps testmtx.o $(LIBDMUMPS) $(LIBDMUMPSX)

testzmumps: $(LIBZMUMPS) testmtxc.o
	-${MF90} -g -o testzmumps testmtxc.o $(LIBZMUMPS) $(LIBZMUMPSX) \
	../lib/tasklib.a $(FLIBS)

testbnd: $(LIBBND) testmtx.o 
	-${FLINKER} -g -o testbnd testmtx.o $(LIBBND) ../lib/libbnd.o

testbndc: $(LIBBNDC) testmtxc.o 
	-${FLINKER} -g -o testbndc testmtxc.o $(LIBBNDC) \
	../lib/tasklib.a $(FLIBS)


testmtx.o: testmtx.f90
testmtxc.o: testmtxc.f90
	$(MF90) $(FFLAGS) -c testmtxc.f90 $(INCMOD)
libmtxmpi.o: libmtxmpi.f90
libmtxnompi.o: libmtxnompi.f90
libmtxbnd.o: libmtxbnd.f90
libmtxbndc.o: libmtxbndc.f90
libmtxdmumps.o: libmtxdmumps.f90
	$(MF90) $(FFLAGS) -c libmtxdmumps.f90 -I$(INCLUDE_MUMPS)
libmtxzmumps.o: libmtxzmumps.f90
	$(MF90) $(FFLAGS) -c libmtxzmumps.f90 -I$(INCLUDE_MUMPS)
libmtxksp.o: libmtxksp.F90
