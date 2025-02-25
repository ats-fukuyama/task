### Makefile for wm
include ../make.header
include ../mtxp/make.mtxp

LIB_MTX=$(LIB_MTX_BND)
LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_PCG)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIB_MTX=$(LIBXMTX_KSP)
#LIBX_MTX=$(LIBX_MTX_KSP)

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

##### Without PETSc library (no MPI, real8 and complex8) #####
FLINKER=$(FCFREE)
LIBNOMPI=../mtxp/libmtxnompi.o
LIBMTXP=$(LIBNOMPI) ../mtxp/libmtxbnd.o
LIBMUMPSX=.

##### With MUMPS library (MPI, real8 and complex8)
#FLINKER=$(MF90)
#LIBMPI=../mtxp/libmtxmpi.o
#LIBMTXP=$(LIBMPI) ../mtxp/libmtxmumps.o
#LIBMUMPSX=-L${PETSC_DIR}/${PETSC_ARCH}/lib -ldmumps -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lparmetis -lmetis -llapack -lblas

MODINCLUDE= -I$(MOD) -I../lib/$(MOD) -I../pl/$(MOD) -I../eq/$(MOD) \
            -I../dp/$(MOD) -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

SRCFIXED = wmfem_comm.f wmf_hhg.f ../eq/wmeqin.f wmmenu.f wminit.f wmeign.f \
       wmexec.f wmloop.f \
       wmxprf.f wmdprf.f \
       wmgout.f wmgsub.f wmdout.f wmfout.f wmtae.f \
       wmfem_allocate.f wmfem_pre.f wmfem_calc.f wmfem_main.f wmfem_post.f \
       wmfem_sub.f wmfem_setg.f wmfem_vmec.f\
       wmvout.f wmtpcs.f

SRCFREE = 

# wmtrfile.f

SRCS = $(SRC0)

OBJS = $(SRCFIXED:.f=.o) $(SRCFREE:.f90=.o)

LIBS = libwm.a ../dp/libdp.a ../eq/libeq.a ../pl/libpl.a ../dp/libdp.a \
       $(LIBMTXP) ../lib/libgrf.a ../lib/libtask.a ../../bpsd/libbpsd.a 

.f.o :
	$(FCFIXED)  $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE)   $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: wm

libwm.a: $(OBJS)
	$(LD) $(LDFLAGS) libwm.a $(OBJS)

../mtxp/libmtxmpi.o:
	(cd ../mtxp; make libmtxmpi.o)
../mtxp/libmtxnompi.o:
	(cd ../mtxp; make libmtxnompi.o)
../mtxp/libmtxksp.o:
	(cd ../mtxp; make libmtxksp.o)
../mtxp/libmtxmumps.o:
	(cd ../mtxp; make libmtxmumps.o)
../mtxp/libmtxbnd.o:
	(cd ../mtxp; make libmtxbnd.o)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../dp/libdp.a:
	(cd ../dp; make libdp.a)
libs:
	(cd ../lib; make libtask.a)
	(cd ../mtxp; make libmtxnompi.o)
	(cd ../mtxp; make libmtxbnd.o)
	(cd ../mtxp; make libmtxmpi.o)
	(cd ../mtxp; make libmtxmumps.o)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)
	(cd ../dp; make libdp.a)

wm : wmmain.o $(LIBS)
	$(FCFIXED) wmmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

clean:
	rm -f core *.o ./*~ ../mtx/*.o *.a $(MOD)/*.mod $(MOD)/*.f90

veryclean: clean
	rm -f wm

WMCOMM = wmcomm.inc wmcom0.inc wmcom1.inc ../pl/plcomm.f90 ../dp/dpcomm.f90

wmmain.o : wmmain.f $(WMCOMM)

wmmenu.o : wmmenu.f $(WMCOMM)
wminit.o : wminit.f $(WMCOMM)
wmeign.o : wmeign.f $(WMCOMM)
wmsetg.o : wmsetg.f $(WMCOMM)

wmloop.o : wmloop.f $(WMCOMM)
wmexec.o : wmexec.f $(WMCOMM)
wmprof.o : wmprof.f $(WMCOMM)
wmdisp.o : wmdisp.f $(WMCOMM)

wmsetm.o : wmsetm.f $(WMCOMM)
wmsetf.o : wmsetf.f $(WMCOMM)
wmxprf.o : wmxprf.f $(WMCOMM)
wmdprf.o : wmdprf.f $(WMCOMM)

wmpout.o : wmpout.f $(WMCOMM)
wmgout.o : wmgout.f $(WMCOMM)
wmgsub.o : wmgsub.f $(WMCOMM)
wmdout.o : wmdout.f $(WMCOMM)	
wmfout.o : wmfout.f $(WMCOMM)	

wmsetm0.o : wmsetm0.f $(WMCOMM)
wmpout0.o : wmpout0.f $(WMCOMM)
wmsetm1.o : wmsetm1.f $(WMCOMM)
wmpout1.o : wmpout1.f $(WMCOMM)

wmtae.o : wmtae.f $(WMCOMM)

wmfem.o : wmfem.f wmfem_comm.f
wmfem_comm.o : wmfem_comm.f
wmf_hermite.o : wmf_hermite.f
wmfem_calc.o : wmfem_calc.f wmfem_comm.f
wmfem_sub.o : wmfem_sub.f wmfem_comm.f $(WMCOMM)	
wmfem_setg.o : wmfem_setg.f wmfem_comm.f $(WMCOMM)
wmfem_vmec.o : wmfem_vmec.f wmfem_comm.f $(WMCOMM) vmcomm.inc

wmvmec.o : wmvmec.f $(WMCOMM) vmcomm.inc
wmbooz.o : wmbooz.f $(WMCOMM) vmcomm.inc
wmvout.o : wmvout.f $(WMCOMM) vmcomm.inc

wmsolv.o : wmsolv.f $(WMCOMM)
nowmsolv.o : nowmsolv.f $(WMCOMM)
wmvout.o : wmvout.f $(WMCOMM) vmcomm.inc
../eq/wmeqin.o: ../eq/wmeqin.f

wmtrfile.o : wmtrfile.f ../tr/trcomm.inc ../tr/trcom0.inc
wmtpcs.o : wmtpcs.f $(WMCOMM)
