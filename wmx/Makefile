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

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)

OBJDIR   = ./obj

#.SUFFIXES:
#.SUFFIXES: .f .f90 .mod .o .a

${OBJDIR}/%.o:%.f
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f90
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

OBJS = $(addprefix $(OBJDIR)/, $(SRC7:.f=.o) $(SRCS:.f90=.o))

LIBS = libwm.a ../dp/libdp.a ../eq/libeq.a ../pl/libpl.a ../tr/libtr2.a \
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
../tr/libtr2.a:
	(cd ../tr; make libtr2.a)
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
	(cd ../tr; make libtr2.a)
	(cd ../dp; make libdp.a)

wm : $(LIB_MTX) $(LIBS) $(OBJDIR)/wmmain.o 
	$(FLINKER) $(OBJDIR)/wmmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

test_nph : $(OBJDIR)/test_nph.o
	$(FLINKER) $(OBJDIR)/test_nph.o -o $@ $(FFLAGS) $(FLIBS)

clean :
	-rm -f core a.out $(OBJDIR)/*.o *.a ./#* ./*~ ./mod/*.mod

veryclean: clean
	-rm -f wm

$(OBJDIR)/wmcomm.o: wmcomm.f90
$(OBJDIR)/vmcomm.o: vmcomm.f90
$(OBJDIR)/wmlib.o: wmlib.f90 wmcomm.f90
$(OBJDIR)/wmsub.o: wmsub.f90 wmcomm.f90 
$(OBJDIR)/wminit.o: wminit.f90 wmcomm.f90
$(OBJDIR)/wmparm.o: wmparm.f90 wmlib.f90 ../dp/dpcomm.f90 ../pl/plcomm.f90 wmcomm.f90
$(OBJDIR)/wmview.o: wmview.f90 wmcomm.f90

$(OBJDIR)/wmxprf.o: wmxprf.f90 ../pl/plprof.f90 wmcomm.f90
$(OBJDIR)/wmdprf.o: wmdprf.f90 wmcomm.f90
$(OBJDIR)/wmdout.o: wmdout.f90 wmcomm.f90

$(OBJDIR)/wmgout.o: wmgout.f90 wmgsub.f90 wmcomm.f90
$(OBJDIR)/wmgsub.o: wmgsub.f90 wmprof.f90 wmcomm.f90

$(OBJDIR)/wmvmec.o: wmvmec.f90 wmgsub.f90 vmcomm.f90 wmcomm.f90
$(OBJDIR)/wmbooz.o: wmbooz.f90 vmcomm.f90 wmcomm.f90
$(OBJDIR)/wmvout.o: wmvout.f90 vmcomm.f90 wmcomm.f90

$(OBJDIR)/wmsetg.o: wmsetg.f90 wmdprf.f90 wmxprf.f90 ../pl/plprof.f90 \
          ../tr/trfile.f90 ../eq/wmeqin.f wmcomm.f90
$(OBJDIR)/wmprof.o: wmprof.f90 ../pl/plprof.f90 wmcomm.f90
$(OBJDIR)/wmdisp.o: wmdisp.f90 wmprof.f90 ../dp/dptnsr0.f90 ../pl/plprofw.f90 \
          ../pl/plprof.f90 wmcomm.f90

$(OBJDIR)/wmsetf0.o: wmsetf0.f90 wmprof.f90 wmdisp.f90 wmsub.f90 wmcomm.f90
$(OBJDIR)/wmsetf1.o: wmsetf1.f90 wmprof.f90 wmdisp.f90 wmsub.f90 wmcomm.f90
$(OBJDIR)/wmsetf2.o: wmsetf2.f90 wmprof.f90 wmdisp.f90 wmsub.f90 wmcomm.f90
$(OBJDIR)/wmsetm0.o: wmsetm0.f90 wmsetf0.f90 ../pl/plprof.f90 wmcomm.f90
$(OBJDIR)/wmsetm1.o: wmsetm1.f90 wmsetf1.f90 ../pl/plprof.f90 wmcomm.f90
$(OBJDIR)/wmsetm2.o: wmsetm2.f90 wmsetf2.f90 ../pl/plprof.f90 wmcomm.f90
$(OBJDIR)/wmsolv.o: wmsolv.f90 wmsetm0.f90 wmsetm2.f90 wmcomm.f90

$(OBJDIR)/wmemfp.o: wmemfp.f90 wmprof.f90 wmsub.f90 wmsetf0.f90 wmsetf2.f90 wmcomm.f90
$(OBJDIR)/wmpout.o: wmpout.f90 wmcomm.f90

$(OBJDIR)/wmexec.o: wmexec.f90 wmsetg.f90 wmsolv.f90 wmemfp.f90 wmpout.f90 \
          wmsub.f90 ../dp/dpparm.f90 wmcomm.f90

$(OBJDIR)/wmfile.o: wmfile.f90 wmsetg.f90 wmexec.f90 ../dp/dpparm.f90 wmcomm.f90
$(OBJDIR)/wmloop.o: wmloop.f90 wmsetg.f90 wmexec.f90 wmfile.f90 wmcomm.f90

$(OBJDIR)/wmeign.o: wmeign.f90 wmparm.f90 wmsetg.f90 wmsolv.f90 wmemfp.f90 wmcomm.f90
$(OBJDIR)/wmtaem.o: wmtaem.f90 wmprof.f90 ../pl/plprof.f90 wmsub.f90 wmcomm.f90

$(OBJDIR)/wmmenu.o: wmmenu.f90 wminit.f90 wmparm.f90 wmview.f90 wmexec.f90 \
          wmgout.f90 wmdout.f90 wmfile.f90 wmeign.f90 wmtaem.f90 wmcomm.f90
$(OBJDIR)/wmmain.o: wmmain.f90 ../pl/plinit.f90 ../eq/equnit.f ../dp/dpinit.f90 \
          ../dp/dpparm.f90 wminit.f90 wmparm.f90 wmmenu.f90 wmcomm.f90
