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

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../pl/$(MOD) -I../eq/$(MOD) -I../dp/$(MOD) \
            -I../lib/$(MOD) -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

SRCS = wmcomm.f90 vmcomm.f90 \
       wminit.f90 wmparm.f90 wmview.f90 \
       wmsetg.f90 wmprof.f90 wmdisp.f90 \
       wmsetm.f90 wmsetf.f90 \
       wmmenu.f90

SRCX = wmxprf.f90 wmdprf.f90 \
       wmpout.f90 wmgout.f90 wmgsub.f90 wmdout.f90 wmfout.f90 \
       wmsolv.f90 wmexec.f90 wmloop.f90 wmtaem.f90 \
       wmvmec.f90 wmbooz.f90 wmvout.f90 wmtpcs.f90 \
       wmeign.f90 \
       wmcont.f90 wmtest.f90

OBJS = $(SRCS:.f90=.o)

LIBS = libwm.a ../dp/libdp.a ../eq/libeq.a ../pl/libpl.a \
       ../lib/libtask.a ../../bpsd/libbpsd.a

all: libs wm

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../dp/libdp.a:
	(cd ../dp; make libdp.a)
$(LIB_MTX):
	(cd ../mtxp; make)
libwm.a: $(OBJS)
	$(LD) $(LDFLAGS) $@ $(OBJS)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)
	(cd ../dp; make libdp.a)

wm : $(LIB_MTX) $(LIBS) wmmain.o 
	$(FLINKER) wmmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean :
	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod

veryclean: clean
	-rm -f wm

wmcomm.o: wmcomm.f90
vmcomm.o: vmcomm.f90

wminit.o: wminit.f90 wmcomm.f90
wmparm.o: wmparm.f90 wmcomm.f90
wmview.o: wmview.f90 wmcomm.f90

wmmenu.0: wmmenu.f90 wmcomm.f90
wmmain.0: wmmain.f90 wmcomm.f90
