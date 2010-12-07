include ../../task/make.header


FFLAGS=$(OFLAGS) -I$(INCLUDE_MUMPS)
#FFLAGS=$(DFLAGS) -I$(INCLUDE_MUMPS)

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

##### With PETSc library(real) #####
#include ${PETSC_DIR}/conf/base
#INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS_4.8.4/include
#INCLUDE_PETSC=${PETSC_DIR}/include
#EXEC=testbnd testksp testdmumps testzmumps testmtxp
#
#LIBMPI=libmtxmpi.o
#LIBNOMPI=libmtxnompi.o
#LIBKSP=$(LIBMPI) libmtxksp.o
#LIBDMUMPS=$(LIBMPI) libmtxdmumps.o
#LIBZMUMPS=$(LIBMPI) libmtxzmumps.o
#LIBBND=$(LIBNOMPI) libmtxbnd.o

##### With PETSc library(complex) #####
include ${PETSC_DIR}/conf/base
INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS_4.8.4/include
INCLUDE_PETSC=${PETSC_DIR}/include
EXEC=testbnd testkspc testdmumps testzmumps testmtxp 

LIBMPI=libmtxmpi.o
LIBNOMPI=libmtxnompi.o
LIBKSPC=$(LIBMPI) libmtxkspc.o
LIBDMUMPS=$(LIBMPI) libmtxdmumps.o
LIBZMUMPS=$(LIBMPI) libmtxzmumps.o
LIBBND=$(LIBNOMPI) libmtxbnd.o


all: $(EXEC)

libs: $(LIBNOMPI) $(LIBMPI) libmtxbnd.o libmtxksp.o libmtxdmumps.o libmtxzmumps.o

testksp: $(LIBKSP) testmtx.o
	-${FLINKER} -g -o testksp testmtx.o $(LIBKSP) ${PETSC_KSP_LIB}

testkspc: $(LIBKSPC) testmtxc.o
	-${FLINKER} -g -o testkspc testmtxc.o $(LIBKSPC) ${PETSC_KSP_LIB} $(GFLIBS) ../lib/tasklib.a

testdmumps: $(LIBDMUMPS) testmtx.o
	-${FLINKER} -g -o testdmumps testmtx.o $(LIBDMUMPS) ${PETSC_KSP_LIB}

testzmumps: $(LIBZMUMPS) testmtxc.o
	-${FLINKER} -g -o testzmumps testmtxc.o $(LIBZMUMPS) ${PETSC_KSP_LIB} ${PETSC_DIR}/intel111/lib/libzmumps.a $(GFLIBS) ../lib/tasklib.a

testbnd: $(LIBBND) testmtx.o 
	-${FLINKER} -g -o testbnd testmtx.o $(LIBBND) ../../task/lib/libbnd.o

testmtxp:  $(LIBDMUMPS) testmtxp.o
	-${FLINKER} -g -o testmtxp testmtxp.o $(LIBDMUMPS) ${PETSC_KSP_LIB}


testmtx.o: testmtx.f90
testmtxc.o: testmtxc.f90
libmtxmpi.o: libmtxmpi.F90
libmtxnompi.o: libmtxnompi.f90
libmtxbnd.o: libmtxbnd.f90
libmtxdmumps.o: libmtxdmumps.f90
libmtxzmumps.o: libmtxzmumps.f90
libmtxksp.o: libmtxksp.F90
libmtxkspc.o: libmtxkspc.F90