### Makefile for wm
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

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../pl/$(MOD) -I../eq/$(MOD) -I../tr/$(MOD) \
            -I../dp/$(MOD) -I../lib/$(MOD) -I../mtxp/$(MOD) \
            -I../../bpsd/$(MOD)

SRCS = wmcomm.f90 vmcomm.f90 wmlib.f90 wmsub.f90 \
       wminit.f90 wmparm.f90 wmview.f90 \
       wmxprf.f90 wmdprf.f90 wmdout.f90 \
       wmprof.f90 wmdisp.f90 \
       wmgsub.f90 wmgout.f90 \
       wmvmec.f90 wmbooz.f90 wmvout.f90 \
       wmsetg.f90 \
       wmsetf0.f90 wmsetf1.f90 wmsetf2.f90 \
       wmsetm0.f90 wmsetm1.f90 wmsetm2.f90 \
       wmsolv.f90 wmemfp.f90 wmpout.f90 \
       wmexec.f90 \
       wmfile.f90 \
       wmloop.f90 \
       wmeign.f90 wmtaem.f90 \
       wmmenu.f90
SRC7 = wmeqin.f

OBJS = $(SRC7:.f=.o) $(SRCS:.f90=.o)

LIBS = libwm.a ../dp/libdp.a ../eq/libeq.a ../pl/libpl.a ../tr/libtr.a \
       ../lib/libtask.a ../lib/libmds.a ../lib/libgrf.a ../../bpsd/libbpsd.a

all: libs wm

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../lib/libmds.a:
	(cd ../lib; make libmds.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../tr/libtr.a:
	(cd ../tr; make libtr.a)
../dp/libdp.a:
	(cd ../dp; make libdp.a)
$(LIB_MTX):
	(cd ../mtxp; make)
libwm.a: $(OBJS)
	$(LD) $(LDFLAGS) $@ $(OBJS)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../lib; make libmds.a)
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)
	(cd ../tr; make libtr.a)
	(cd ../dp; make libdp.a)

wm : $(LIB_MTX) $(LIBS) wmmain.o 
	$(FLINKER) wmmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean :
	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod

veryclean: clean
	-rm -f wm

wmcomm.o: wmcomm.f90
vmcomm.o: vmcomm.f90
wmlib.o: wmlib.f90 wmcomm.f90
wmsub.o: wmsub.f90 wmcomm.f90 
wminit.o: wminit.f90 wmcomm.f90
wmparm.o: wmparm.f90 wmlib.f90 ../dp/dpcomm.f90 ../pl/plcomm.f90 wmcomm.f90
wmview.o: wmview.f90 wmcomm.f90

wmxprf.o: wmxprf.f90 ../pl/plprof.f90 wmcomm.f90
wmdprf.o: wmdprf.f90 wmcomm.f90
wmdout.o: wmdout.f90 wmcomm.f90

wmgout.o: wmgout.f90 wmgsub.f90 wmcomm.f90
wmgsub.o: wmgsub.f90 wmprof.f90 wmcomm.f90

wmvmec.o: wmvmec.f90 wmgsub.f90 vmcomm.f90 wmcomm.f90
wmbooz.o: wmbooz.f90 vmcomm.f90 wmcomm.f90
wmvout.o: wmvout.f90 vmcomm.f90 wmcomm.f90

wmsetg.o: wmsetg.f90 wmdprf.f90 wmxprf.f90 ../pl/plprof.f90 \
          ../tr/trfile.f90 ../eq/wmeqin.f wmcomm.f90
wmprof.o: wmprof.f90 ../pl/plprof.f90 wmcomm.f90
wmdisp.o: wmdisp.f90 wmprof.f90 ../dp/dptnsr0.f90 ../pl/plprofw.f90 \
          ../pl/plprof.f90 wmcomm.f90

wmsetf0.o: wmsetf0.f90 wmprof.f90 wmdisp.f90 wmsub.f90 wmcomm.f90
wmsetf1.o: wmsetf1.f90 wmprof.f90 wmdisp.f90 wmsub.f90 wmcomm.f90
wmsetf2.o: wmsetf2.f90 wmprof.f90 wmdisp.f90 wmsub.f90 wmcomm.f90
wmsetm0.o: wmsetm0.f90 wmsetf0.f90 ../pl/plprof.f90 wmcomm.f90
wmsetm1.o: wmsetm1.f90 wmsetf1.f90 ../pl/plprof.f90 wmcomm.f90
wmsetm2.o: wmsetm2.f90 wmsetf2.f90 ../pl/plprof.f90 wmcomm.f90
wmsolv.o: wmsolv.f90 wmsetm0.f90 wmsetm2.f90 wmcomm.f90

wmemfp.o: wmemfp.f90 wmprof.f90 wmsub.f90 wmsetf0.f90 wmsetf2.f90 wmcomm.f90
wmpout.o: wmpout.f90 wmcomm.f90

wmexec.o: wmexec.f90 wmsetg.f90 wmsolv.f90 wmemfp.f90 wmpout.f90 \
          wmsub.f90 ../dp/dpparm.f90 wmcomm.f90

wmfile.o: wmfile.f90 wmsetg.f90 wmexec.f90 ../dp/dpparm.f90 wmcomm.f90
wmloop.o: wmloop.f90 wmsetg.f90 wmexec.f90 wmfile.f90 wmcomm.f90

wmeign.o: wmeign.f90 wmparm.f90 wmsetg.f90 wmsolv.f90 wmemfp.f90 wmcomm.f90
wmtaem.o: wmtaem.f90 wmprof.f90 ../pl/plprof.f90 wmsub.f90 wmcomm.f90

wmmenu.o: wmmenu.f90 wminit.f90 wmparm.f90 wmview.f90 wmexec.f90 \
          wmgout.f90 wmdout.f90 wmfile.f90 wmeign.f90 wmtaem.f90 wmcomm.f90
wmmain.o: wmmain.f90 ../pl/plinit.f90 ../eq/equnit.f ../dp/dpinit.f90 \
          ../dp/dpparm.f90 wminit.f90 wmparm.f90 wmmenu.f90 wmcomm.f90
