### $Id$ ###
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

MODINCLUDE= -I$(MOD) -I../lib/$(MOD) -I../dp/$(MOD) \
            -I../pl/$(MOD) -I../eq/$(MOD) \
            -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

SRCFIXED1 = wmfem_comm.f  ../eq/wmeqin.f
SRCFIXED2 =        wmexec.f wmloop.f \
       wmxprf.f wmdprf.f \
       wmgout.f 
SRCFIXED3 = wmgsub.f wmdout.f wmfout.f wmtae.f \
       wmfem_allocate.f wmfem_pre.f wmfem_main.f wmfem_post.f \
       wmfem_sub.f wmfem_setg.f wmfem_vmec.f\
       wmvout.f wmtpcs.f \
       wmmenu.f wminit.f wmeign.f

SRCFREE1 = fem_hhg.f90
SRCFREE2 = wmfem_calc.f90

# wmtrfile.f

SRCS = $(SRC0)

OBJS = $(SRCFIXED1:.f=.o) $(SRCFREE1:.f90=.o) $(SRCFIXED2:.f=.o) $(SRCFREE2:.f90=.o) $(SRCFIXED3:.f=.o) 

LIB0 = ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a ../pl/libpl.a \
       ../eq/libeq.a ../dp/libdp.a

LIBS = libwm.a ../dp/libdp.a ../eq/libeq.a ../pl/libpl.a \
       $(LIBMTXP) ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a 

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

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
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../dp/libdp.a:
	(cd ../dp; make libdp.a)
libs:
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../mtxp; make libmtxnompi.o)
	(cd ../mtxp; make libmtxbnd.o)
	(cd ../mtxp; make libmtxmpi.o)
	(cd ../mtxp; make libmtxmumps.o)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)
	(cd ../dp; make libdp.a)

wm : wmmain.o $(LIBS)
	$(FLINKER) wmmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	rm -f core *.o ./*~ ../mtx/*.o *.a $(MOD)/*.mod

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
wmfem_calc.o : wmfem_calc.f90 wmfem_comm.f
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
