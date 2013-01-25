include ../../task/make.header

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)


##### Without PETSc library (no MPI, real8 and complex8) #####
FLINKER=$(FCFREE)
INCLUDE_MUMPS=.
EXEC=testbnd testbndc
all: $(EXEC)
.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@
clean :
	-rm -f core a.out *.o *.a ./#* ./*~ ./*.mod

##### With MUMPS library (real8 and complex8)
#FLINKER=$(MF90)
#INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS_4.10.0/include
#EXEC=testbnd testbndc testdmumps testzmumps
#all: $(EXEC)
#.SUFFIXES:
#.SUFFIXES: .f .f90 .mod .o .a
#.f90.o :
#	$(FCFREE) $(FFLAGS) -c $< -o $@
#clean :
#	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod

##### With PETSc library 3.2 (real8) #####
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules
#INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS_4.10.0/include
#INCLUDE_PETSC=${PETSC_DIR}/include
#EXEC=testbnd testbndc testdmumps testzmumps testksp
#all: $(EXEC)

INCMOD=-I $(MOD) -I../lib/$(MOD)
LIBMOD=-L../lib -l libgrf
LIBMPI=libmtxmpi.o
LIBNOMPI=libmtxnompi.o
LIBBND=$(LIBNOMPI) libmtxbnd.o
LIBKSP=$(LIBMPI) libmtxksp.o
LIBMUMPS=$(LIBMPI) libmtxmumps.o
LIBMUMPSX=-L${PETSC_DIR}/${PETSC_ARCH}/lib -ldmumps -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lparmetis -lmetis -llapack -lblas

libs: $(LIBNOMPI) $(LIBMPI) libmtxbnd.o libmtxmumps.o libmtxksp.o 

testbnd: $(LIBBND) testmtx.o 
	-${FLINKER} -g -o testbnd testmtx.o $(LIBBND) ../lib/libbnd.o

testbndc: $(LIBBNDC) testmtxc.o 
	-${FLINKER} -g -o testbndc testmtxc.o $(LIBBND) \
	../lib/tasklib.a $(FLIBS)

testksp: $(LIBKSP) testmtx.o
	-${FLINKER} -g -o testksp testmtx.o $(LIBKSP) ${PETSC_KSP_LIB}

testdmumps: $(LIBMUMPS) testmtx.o
	-${MF90} -g -o testdmumps testmtx.o $(LIBMUMPS) $(LIBMUMPSX)

testzmumps: $(LIBMUMPS) testmtxc.o
	-${MF90} -g -o testzmumps testmtxc.o $(LIBMUMPS) $(LIBMUMPSX) \
	../lib/tasklib.a $(FLIBS)

testmtx.o: testmtx.f90
	$(FCFREE) $(FFLAGS) -c testmtx.f90 $(INCMOD)
testmtxc.o: testmtxc.f90
	$(FCFREE) $(FFLAGS) -c testmtxc.f90 $(INCMOD)
libmtxmpi.o: libmtxmpi.f90
	$(MF90) $(FFLAGS) -c libmtxmpi.f90 $(INCMOD)
libmtxnompi.o: libmtxnompi.f90
libmtxbnd.o: libmtxbnd.f90
libmtxmumps.o: libmtxmumps.f90
	$(MF90) $(FFLAGS) -c libmtxmumps.f90 -I$(INCLUDE_MUMPS)
libmtxksp.o: libmtxksp.F90
