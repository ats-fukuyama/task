### $Id: Makefile,v 1.69 2012/07/25 07:53:52 nuga Exp $ ###
include ../make.header
include ../mtxp/make.mtxp
#LIB_MTX=$(LIB_MTX_BND)
#LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_PCG)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIB_MTX=$(LIBXMTX_KSP)
#LIBX_MTX=$(LIBX_MTX_KSP)

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../lib/$(MOD) -I../pl/$(MOD) -I../eq/$(MOD)  -I../dp/$(MOD) -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

SRCS = fpcomm.f90 eg_read.f90 fpinit.f90 fpparm.f90 \
       fpmpi.f90 fpbroadcast.f90 fpsub.f90 \
       fpexec.f90 fpsave.f90 \
       fpcaltp.f90 fpcalte.f90 fpcaldeff.f90 fpcalchieff.f90 \
       fpoutdata.f90 fpwmin.f90 fpwrin.f90 fpcalwm.f90 fpcalwr.f90 fpcalw.f90 \
       fpcalcn.f90 fpcalcnr.f90 fpcalc.f90 fpnfrr.f90 fpnflg.f90 \
       cdbmfp.f90 fpcalr.f90 sv_read.f90 fpcoef.f90 fpcale.f90 \
       fpbounce.f90 fpdisrupt.f90 fpprep.f90 fploop.f90 \
       fpcont.f90 fpfout.f90 fpgsub.f90 fpgout.f90 fpfile.f90 fpmenu.f90

OBJS = $(SRCS:.f90=.o)

LIBS = libfp.a ../pl/libpl.a ../eq/libeq.a ../dp/libdp.a \
        ../lib/libgrf.a ../lib/libtask.a ../../bpsd/libbpsd.a

all: fp

libfp.a: $(OBJS)
	$(LD) $(LDFLAGS) libfp.a $(OBJS)

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../dp/libeq.a:
	(cd ../dp; make libdp.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)
	(cd ../dp; make libdp.a)

fp : libs libfp.a fpmain.o 
	$(FLINKER) fpmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean :
	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod

veryclean: clean
	-rm -f fp

new:
	-mkdir ../fpnew
	cp -f Makefile ../fpnew
	cp -f *.f ../fpnew
	cp -f *.inc ../fpnew

fpmain.o : fpinit.f90 fpparm.f90 fpmain.f90 fpcomm.f90
fpmenu.o : fpwrin.f90 fpgout.f90 fpfout.f90 fploop.f90 fpprep.f90 fpparm.f90 fpfile.f90 fpmenu.f90 fpcomm.f90
fpinit.o : fpinit.f90 fpcomm.f90
fpparm.o : fpparm.f90 fpcomm.f90
fpprep.o : fpoutdata.f90 fpcale.f90 fpbounce.f90 fpnfrr.f90 fpwmin.f90 fpwrin.f90 fpsave.f90 fpcalw.f90 fpbroadcast.f90 fpparm.f90 fpcoef.f90 fpdisrupt.f90 eg_read.f90 sv_read.f90 fpprep.f90 fpcomm.f90
fploop.o : fpoutdata.f90 fpnfrr.f90 fpsave.f90 fpexec.f90 fpprep.f90 fpcoef.f90 fpdisrupt.f90 eg_read.f90 fploop.f90 fpcomm.f90 fpoutdata.f90
fpexec.o : fpexec.f90 fpcomm.f90
fpdrexec.o : fpdrexec.f90 fpcomm.f90
fpbroadcast.o : fpbroadcast.f90 fpcomm.f90
fpcoef.o : fpcalr.f90 fpnfrr.f90 fpcalwm.f90 fpcalwr.f90 fpcalw.f90 fpcalc.f90 fpcoef.f90 fpcomm.f90
fpcalc.o : fpcalcnr.f90 fpcalcn.f90 fpcalc.f90 fpcomm.f90
fpcalcn.o : fpcalcn.f90 fpcomm.f90
fpcalcnr.o : fpcoef.f90 fpcalcnr.f90 fpcomm.f90
fpcalw.o : fpcalw.f90 fpcalwr.f90 fpcalwm.f90 fpcomm.f90
fpcalwr.o : fpwrin.f90 fpcalwr.f90 fpcomm.f90
fpcalwm.o : fpwmin.f90 fpcalwm.f90 fpcomm.f90
fpcont.o : fpcont.f90 fpcomm.f90
fpfout.o : fpfout.f90 fpcomm.f90
fpgsub.o : fpgsub.f90 fpcomm.f90
fpgout.o : fpfout.f90 fpcont.f90 fpgsub.f90 fpgout.f90 fpcomm.f90
fpfile.o : fpprep.f90 eg_read.f90 sv_read.f90 fpfile.f90 fpcomm.f90
fpsave.o : fpexec.f90 eg_read.f90 fpsave.f90 fpcomm.f90
fpwrin.o : fpwrin.f90 fpcomm.f90
fpwmin.o : fpwrin.f90 fpcomm.f90
fpnfrr.o : fpparm.f90 fpnfrr.f90 fpcomm.f90
fpbounce.o : fpparm.f90 fpbounce.f90 fpcomm.f90
fpcale.o: fpcale.f90 fpmpi.f90
fpdisrupt.o: fpcale.f90 fpcoef.f90 fpdisrupt.f90 fpcomm.f90
fpcalr.o: fpcalr.f90 fpcomm.f90
eg_read.o: eg_read.f90 fpcomm.f90
sv_read.o: sv_read.f90 fpcomm.f90
fpoutdata.o: fpoutdata.f90 eg_read.f90 fpcomm.f90
fpcaltp.o: fpcaltp.f90 fpcomm.f90
fpcaldeff.o: fpcaldeff.f90 fpcomm.f90
