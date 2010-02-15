include ../../task/make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

FLINKER=$(FCFREE)
FFLAGS=$(DFLAGS)

#include ${PETSC_DIR}/conf/base
#INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS_4.8.4/include
#INCLUDE_PETSC=${PETSC_DIR}/include
#FFLAGS=-g -I${INCLUDE_MUMPS}

LIBMPI=libmtxmpi.o
LIBNOMPI=libmtxnompi.o
LIBKSP=$(LIBMPI) libmtxksp.o
LIBMUMPS=$(LIBMPI) libmtxdmumps.o
LIBBND=$(LIBNOMPI) libmtxbnd.o

#all: testbnd testksp testdmumps
all: testbnd

libs: $(LIBNOMPI) $(LIBMPI) libmtxbnd.o libmtxksp.o libmtxdmumps.o

testksp: $(LIBKSP) testmtx.o
	-${FLINKER} -g -o testksp testmtx.o $(LIBKSP) ${PETSC_KSP_LIB}

testdmumps: $(LIBMUMPS) testmtx.o
	-${FLINKER} -g -o testdmumps testmtx.o $(LIBMUMPS) ${PETSC_KSP_LIB}

testbnd: $(LIBBND) testmtx.o 
	-${FLINKER} -g -o testbnd testmtx.o $(LIBBND) ../../task/lib/libbnd.o

testmtx.o: testmtx.f90
libmtxmpi.o: libmtxmpi.F90
libmtxnompi.o: libmtxnompi.f90
libmtxbnd.o: libmtxbnd.f90
libmtxdmumps.o: libmtxdmumps.f90
libmtxksp.o: libmtxksp.F90
