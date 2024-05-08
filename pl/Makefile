### $Id$ ###
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

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

OBJDIR=./obj

MODINCLUDE=-I$(MOD) -I../eq/$(MOD) -I../mtxp/$(MOD) \
           -I../lib/$(MOD) -I../../bpsd/$(MOD)

${OBJDIR}/%.o:%.f
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

SRCS = plvmec.f90 plcomm.f90 plinit.f90 plparm.f90 \
       platmd.f90 plcoll.f90 \
       plload.f90 plprof_travis.f90 plprof2d.f90 \
       plprof.f90 plprofw.f90  \
       plgout.f90 plmenu.f90
SRCO = noeqlib.f noequlib.f

OBJS=$(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))
OBJO=$(addprefix $(OBJDIR)/, $(SRCO:.f90=.o))

LIBS = libpl.a ../eq/libeq.a ../lib/libgrf.a ../lib/libtask.a ../../bpsd/libbpsd.a

all : libs pl

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
libpl.a: $(OBJS)
	$(LD) $(LDFLAGS) libpl.a $(OBJS)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../mtxp; make)
	(cd ../eq; make libeq.a)

pl : libs $(LIB_MTX) libpl.a $(OBJDIR)/plmain.o
	$(FLINKER) $(OBJDIR)/plmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out $(OBJDIR)/*.o *.mod ./*~ ./#* *.a $(MOD)/*.mod

veryclean: clean
	-rm -f pl

$(OBJDIR)/plvmec.o: plvmec.f90 ../../bpsd/bpsd_kinds.f90 ../../bpsd/bpsd_constants.f90
$(OBJDIR)/plcomm.o: plcomm.f90 ../../bpsd/bpsd_kinds.f90 ../../bpsd/bpsd_constants.f90
$(OBJDIR)/plinit.o: plinit.f90 plcomm.f90
$(OBJDIR)/plparm.o: plparm.f90 plcomm.f90
$(OBJDIR)/platmd.o: platmd.f90 plcomm.f90
$(OBJDIR)/plcoll.o: plcoll.f90 plcomm.f90 platmd.f90
$(OBJDIR)/plload.o: plload.f90 plcomm.f90
$(OBJDIR)/plprof_travis.o: plprof_travis.f90 plcomm.f90
$(OBJDIR)/plprof2d.o: plprof2d.f90 plcoll.f90 plload.f90 plcomm.f90
$(OBJDIR)/plprof.o: plprof.f90 plcomm.f90 plload.f90 \
          plprof2d.f90 plcoll.f90 plprof_travis.f90
$(OBJDIR)/plprofw.o: plprofw.f90 plcomm.f90 plprof.f90
$(OBJDIR)/plgout.o: plgout.f90 plcomm.f90 plprof.f90
$(OBJDIR)/plmenu.o: plmenu.f90 plcomm.f90 plparm.f90 plgout.f90 plload.f90 plvmec.f90
$(OBJDIR)/plmain.o: plmain.f90 plcomm.f90 plinit.f90 plparm.f90 plmenu.f90
