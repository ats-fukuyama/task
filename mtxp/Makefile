include ../make.header
include ../mtxp/make.mtxp

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a
.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE=-I./$(MOD) -I../lib/$(MOD)

all: $(EXEC)

testbnd: $(LIB_MTX_BND) testmtx.o 
	-$(FLINKER) -g -o testbnd testmtx.o $(LIBX_MTX_BND) \
	../lib/tasklib.a $(FLIBS)

testbndc: $(LIB_MTX_BND) testmtxc.o 
	-$(FLINKER) -g -o testbndc testmtxc.o $(LIBX_MTX_BND) \
	../lib/tasklib.a $(FLIBS)

testpcg: $(LIB_MTX_PCG) testmtx.o 
	-$(FLINKER) $(FFLAGS) -o testpcg testmtx.o $(LIBX_MTX_PCG) \
	../lib/tasklib.a $(FLIBS)


testdmumps: $(LIB_MTX_MUMPS) testmtx.o
	-${FLINKER} -g -o testdmumps testmtx.o $(LIBX_MTX_MUMPS) \
	../lib/tasklib.a $(FLIBS)

testzmumps: $(LIB_MTX_MUMPS) testmtxc.o
	-${FLINKER} -g -o testzmumps testmtxc.o $(LIBX_MTX_MUMPS) \
	../lib/tasklib.a $(FLIBS)


testksp: $(LIB_MTX_KSP) testmtx.o
	-${FLINKER} -g -o testksp testmtx.o $(LIBX_MTX_KSP) \
	../lib/tasklib.a $(FLIBS)


testmtx.o: testmtx.f90
	$(FLINKER) $(FFLAGS) -c testmtx.f90 $(MODINCLUDE)
testmtxc.o: testmtxc.f90
	$(FLINKER) $(FFLAGS) -c testmtxc.f90 $(MODINCLUDE)
../mtxp/libmtxmpi.o: libmtxmpi.f90
	$(FLINKER) $(FFLAGS) -c libmtxmpi.f90 $(MODDIR) $(MODINCLUDE)
../mtxp/libmtxnompi.o: libmtxnompi.f90
	$(FLINKER) $(FFLAGS) -c libmtxnompi.f90 $(MODDIR) $(MODINCLUDE)
../mtxp/libmtxbnd.o: libmtxbnd.f90
	$(FLINKER) $(FFLAGS) -c libmtxbnd.f90 $(MODDIR) $(MODINCLUDE)
../mtxp/libmtxbndmpi.o: libmtxbndmpi.f90
	$(FLINKER) $(FFLAGS) -c libmtxbndmpi.f90 $(MODDIR) $(MODINCLUDE) $(INCLUDE_MTX_MPI)
../mtxp/libmtxpcg.o: libmtxpcg.f90
	$(FCFREE) $(FFLAGS) -c libmtxpcg.f90 $(MODDIR) $(MODINCLUDE)
../mtxp/libmtxpcgmpi.o: libmtxpcgmpi.f90
	$(FLINKER) $(FFLAGS) -c libmtxpcgmpi.f90 $(MODDIR) $(MODINCLUDE) $(INCLUDE_MTX_MPI)
../mtxp/libmtxmumps.o: libmtxmumps.f90
	$(FLINKER) $(FFLAGS) -c libmtxmumps.f90 $(MODDIR) $(MODINCLUDE) $(INCLUDE_MTX_MUMPS)
../mtxp/libmtxksp.o: libmtxksp.F90
	$(FLINKER) $(FFLAGS) -c libmtxksp.F90 $(MODDIR) $(MODINCLUDE) \
	$(INCLUDE_MTX_PETSC)

clean :
	-rm -f core a.out *.o *.a *.i90 ./#* ./*~ ./*.mod ./mod/*.mod

veryclean : clean
	-rm -f $(EXEC)
