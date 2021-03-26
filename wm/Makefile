### $Id$ ###
include ../make.header
include ../mtxp/make.mtxp

#LIB_MTX=$(LIB_MTX_BND)
#LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../pl/$(MOD) -I../eq/$(MOD) -I../tr/$(MODE) \
            -I../dp/$(MOD) -I../lib/$(MOD) -I../mtxp/$(MOD) \
            -I../../bpsd/$(MOD)

SRCS = mpilib.f wmmenu.f wminit.f wmeign.f wmsetg.f \
       ../eq/wmeqin.f wmprof.f wmdisp.f \
       wmsetm.f wmsetf.f wmxprf.f wmdprf.f wmsolv.f \
       wmpout.f wmgout.f wmgsub.f wmdout.f wmfout.f \
       wmexec.f wmloop.f wmtae.f \
       wmvmec.f wmbooz.f wmvout.f wmtpcs.f

SRCX = wmcont.f90 wmtest.f90 
#wmdpfa.f90 wmdpfa2.f90 wmdpfaA.f90 wmdpfa3.f90

OBJS = $(SRCX:.f90=.o) $(SRCS:.f=.o)

LIBS = libwm.a ../dp/libdp.a ../pl/libpl.a  ../eq/libeq.a \
       ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a

all: wm

libwm.a: $(OBJS)
	$(LD) $(LDFLAGS) libwm.a $(OBJS)

../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../eq/wmeqin.o:
	(cd ../eq; make wmeqin.o)
../dp/libdp.a:
	(cd ../dp; make libdp.a)
libs:
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)
	(cd ../dp; make libdp.a)

wm : libs libwm.a wmmain.o
	$(FLINKER) wmmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	rm -f core *.o ./*~ *.a $(MOD)/*.mod

veryclean: clean
	rm -f wm

WMCOMM = wmcomm.inc wmcom0.inc wmcom1.inc ../pl/plcomm.f90 ../dp/dpcomm.f90

wmmain.o : wmmain.f $(WMCOMM)

wmmenu.o : wmmenu.f $(WMCOMM)
wminit.o : wminit.f $(WMCOMM)
wmeign.o : wmeign.f $(WMCOMM)
wmsetg.o : wmsetg.f $(WMCOMM)

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

wmsetm.o : wmsetm.f $(WMCOMM)
wmsetm.o : wmsetm.f $(WMCOMM)

wmtae.o : wmtae.f $(WMCOMM)

wmvmec.o : wmvmec.f $(WMCOMM) vmcomm.inc
wmbooz.o : wmbooz.f $(WMCOMM) vmcomm.inc
wmvout.o : wmvout.f $(WMCOMM) vmcomm.inc

wmsolv.o : wmsolv.f $(WMCOMM)
wmvout.o : wmvout.f $(WMCOMM) vmcomm.inc

../eq/wmeqin.o : ../eq/wmeqin.f
wmtrfile.o : wmtrfile.f ../tr/trcomm.inc ../tr/trcom0.inc
wmtpcs.o : wmtpcs.f $(WMCOMM)

wmloop.o : wmloop.f $(WMCOMM)
wmcont.o : wmcont.f90
wmtest.o : wmtest.f90
