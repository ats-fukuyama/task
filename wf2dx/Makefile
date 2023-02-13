### $Id$ ###
include ../make.header
include ../mtxp/make.mtxp

#LIB_MTX=$(LIB_MTX_BND)
#LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

OBJDIR   = ./obj
${OBJDIR}/%.o:%.f
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I$(MOD) -I../pl/$(MOD) -I../eq/$(MOD) \
            -I../dp/$(MOD) -I../lib/$(MOD) -I../mtxp/$(MOD) \
            -I../../bpsd/$(MOD) \
            -I../lib/$(MOD)

SRCS = 	wfcomm.f90 femmeshprep.f90 feminterpolate.f90 \
        wfinit.f90 wfparm.f90 wfload.f90 wfprof.f90 \
        wfdiv.f90 wfsub.f90 \
        wfindx.f90 wfsort.f90 wfgsub.f90 wfant.f90 wfwave.f90 \
	wfsolv.f90 wfgout.f90 wftest.f90 wfmenu.f90

OBJS = $(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))

LIBS = libwf.a ../dp/libdp.a ../pl/libpl.a ../eq/libeq.a   \
        ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a

all: libs wf

libwf.a: $(OBJS)
	$(LD) $(LDFLAGS) libwf.a $(OBJS)

../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make liobeq.a)
../dp/libdp.a:
	(cd ../dp; make libdp.a)
$(LIB_MTX):
	(cd ../mtxp; make)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpd.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)

wf : libs $(LIBS) $(OBJDIR)/wfmain.o
	$(FLINKER) $(OBJDIR)/wfmain.o $(LIBS) -o $@ $(FFLAGS) \
$(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out $(OBJDIR)/*.o *.a ./#* ./*~ ./mod/*.mod

veryclean: clean
	rm -f wf

$(OBJDIR)/wfcomm.o: wfcomm.f90 ../pl/plcomm.f90
$(OBJDIR)/wfmain.o: wfmain.f90 wfcomm.f90
$(OBJDIR)/wfmenu.o: wfmenu.f90 wfcomm.f90
$(OBJDIR)/wfinit.o: wfinit.f90 wfcomm.f90
$(OBJDIR)/wfparm.o: wfparm.f90 wfcomm.f90
$(OBJDIR)/wfload.o: wfload.f90 wfcomm.f90
$(OBJDIR)/wfprof.o: wfprof.f90 wfcomm.f90
$(OBJDIR)/wfdiv.o: wfdiv.f90 wfcomm.f90
$(OBJDIR)/wfsub.o: wfsub.f90 wfcomm.f90
$(OBJDIR)/wfindx.o: wfindx.f90 wfcomm.f90
$(OBJDIR)/wfsort.o: wfsort.f90 wfcomm.f90
$(OBJDIR)/wfgsub.o: wfgsub.f90 wfcomm.f90
$(OBJDIR)/wfant.o: wfant.f90 wfcomm.f90
$(OBJDIR)/wfwave.o: wfwave.f90 wfcomm.f90
$(OBJDIR)/wfsolv.o: wfsolv.f90 wfcomm.f90
$(OBJDIR)/wfgout.o: wfgout.f90 wfcomm.f90
$(OBJDIR)/wftest.o: wftest.f90 wfcomm.f90
$(OBJDIR)/feminterpolate.o: feminterpolate.f90 wfcomm.f90
$(OBJDIR)/femmeshprep.o: femmeshprep.f90 wfcomm.f90
