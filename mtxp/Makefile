include ../../task/make.header

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

##### Without PETSc library (no MPI, real8 and complex8) #####
FLINKER=$(FCFREE)
INCLUDE_MUMPS=.
EXEC=testbnd testbndc testpcg
all: $(EXEC)
.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@
clean :
	-rm -f core a.out *.o *.a ./#* ./*~ ./*.mod

#### With MUMPS library (real8 and complex8)
#EXEC=testbnd testbndc testdmumps testzmumps
#all: $(EXEC)
#PETSc 3.4.2
#INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS/include
#PETSc 3.5.2
#INCLUDE_MUMPS=${PETSC_DIR}/${PETSC_ARCH}/externalpackages/MUMPS/include
#FLINKER=$(MF90)
#.SUFFIXES:
#.SUFFIXES: .f .f90 .mod .o .a
#.f90.o :
#	$(FCFREE) $(FFLAGS) -c $< -o $@
#clean :	
#	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod

##### With PETSc library 3.3 (real8) #####
#EXEC=testbnd testpcg testbndc testdmumps testzmumps testksp
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules
#PETSc 3.4.2
#INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS/include
#PETSc 3.5.2
#INCLUDE_MUMPS=${PETSC_DIR}/${PETSC_ARCH}/externalpackages/MUMPS/include
#INCLUDE_PETSC=${PETSC_DIR}/include
#INCLUDE_PETSC2=${PETSC_DIR}/${PETSC_ARCH}/include

INCMOD=-I $(MOD) -I../lib/$(MOD)
LIBMOD=-L../lib -l libgrf
#LIBMPI=libmtxmpi.o -L/usr/local/lib -lpmi
LIBMPI=libmtxmpi.o
LIBNOMPI=libmtxnompi.o
LIBBND=$(LIBNOMPI) libmtxbnd.o
LIBPCG=$(LIBNOMPI) libmtxpcg.o
LIBMUMPS=$(LIBMPI) libmtxmumps.o
LIBMUMPSX=-L${PETSC_DIR}/${PETSC_ARCH}/lib -ldmumps -lzmumps -lmumps_common -lpord -lscalapack -lparmetis -lmetis -llapack -lblas ../lib/tasklib.a
LIBKSP=$(LIBMPI) libmtxksp.o
LIBKSPX=-L${PETSC_LIB_DIR} -lpetsc -mkl -ldmumps -lzmumps -lmumps_common --lpord -lscalapack -lparmetis -lmetis -llapack -lblas

testbnd: $(LIBBND) testmtx.o 
	-${FCFREE} -g -o testbnd testmtx.o $(LIBBND) ../lib/libbnd.o

testbndc: $(LIBBND) testmtxc.o 
	-$(FCFREE) -g -o testbndc testmtxc.o $(LIBBND) \
	../lib/tasklib.a $(FLIBS)

testpcg: $(LIBPCG) testmtx.o 
	-$(FCFREE) $(FFLAGS) -o testpcg testmtx.o $(LIBPCG) ../lib/libpcgpme.o

testksp: $(LIBKSP) testmtx.o
	-${FLINKER} -g -o testksp testmtx.o $(LIBKSP) ${PETSC_KSP_LIB}  

testdmumps: $(LIBMUMPS) testmtx.o
	-${FLINKER} -g -o testdmumps testmtx.o $(LIBMUMPS) $(LIBMUMPSX) $(FLIBS)

testzmumps: $(LIBMUMPS) testmtxc.o
	-${FLINKER} -g -o testzmumps testmtxc.o $(LIBMUMPS) $(LIBMUMPSX) $(FLIBS)


testmtx.o: testmtx.f90
	$(FLINKER) $(FFLAGS) -c testmtx.f90 $(INCMOD)
testmtxc.o: testmtxc.f90
	$(FLINKER) $(FFLAGS) -c testmtxc.f90 $(INCMOD)
libmtxmpi.o: libmtxmpi.f90
	$(FLINKER) $(FFLAGS) -c libmtxmpi.f90 $(INCMOD)
libmtxnompi.o: libmtxnompi.f90
	$(FLINKER) $(FFLAGS) -c libmtxnompi.f90 $(INCMOD)
libmtxbnd.o: libmtxbnd.f90
libmtxpcg.o: libmtxpcg.f90
	$(FCFREE) $(FFLAGS) -c libmtxpcg.f90 $(INCMOD)
libmtxmumps.o: libmtxmumps.f90
	$(FLINKER) $(FFLAGS) -c libmtxmumps.f90 -I$(INCLUDE_MUMPS)
libmtxksp.o: libmtxksp.F90
	$(FLINKER) $(FFLAGS) -c libmtxksp.F90 $(INCMOD) -I$(INCLUDE_PETSC) -I$(INCLUDE_PETSC2) 
