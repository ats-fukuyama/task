### $Id: Makefile,v 1.69 2012/07/25 07:53:52 nuga Exp $ ###
include ../make.header
include ../mtxp/make.mtxp
# LIB_MTX=$(LIB_MTX_BND)
# LIBX_MTX=$(LIBX_MTX_BND)
# LIB_MTX=$(LIB_MTX_PCG)
# LIBX_MTX=$(LIBX_MTX_PCG)
LIB_MTX=$(LIB_MTX_MUMPS)
LIBX_MTX=$(LIBX_MTX_MUMPS)
# LIB_MTX=$(LIBXMTX_KSP)
# LIBX_MTX=$(LIBX_MTX_KSP)

FLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

OBJDIR   = ./obj

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	@[ -d $(OBJDIR) ]
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../lib/$(MOD) -I../pl/$(MOD) -I../eq/$(MOD) \
						-I../mtxp/$(MOD) -I../../bpsd/$(MOD) -I../ob/$(MOD)

SRCS = fpcomm.f90 eg_read.f90 fpinit.f90 fpmpi.f90 fpbroadcast.f90 fpsub.f90 \
       fpexec.f90 fpsave.f90 fpwrite.f90 fpdifcoef.f90 fprad.f90\
       fpoutdata.f90 fpwmin.f90 fpwrin.f90 fpcalwm.f90 fpcalwr.f90 fpcalw.f90 \
       fpcalcn.f90 fpcalcnr.f90 fpcalc.f90 fpnfrr.f90 fpnflg.f90 \
       cdbmfp.f90 fpcalr.f90 sv_read.f90 fpcoef.f90 fpcale.f90 \
       fpbounce.f90 fpdisrupt.f90 fpprep.f90 fploop.f90 \
       fpcont.f90 fpfout.f90 fpgout.f90 fpfile.f90 fpmenu.f90

SRCSFOW = fowcomm.f90 fowprep.f90  foworbitclassify.f90 foworbit.f90 fowbounceaverage.f90 fowdistribution.f90

# OBJS = $(SRCS:.f90=.o)
OBJS = $(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))
OBJSFOW = $(addprefix $(OBJDIR)/, $(SRCSFOW:.f90=.o))


LIBS = fplib.a ../pl/libpl.a ../eq/libeq.a  \
       ../lib/libgrf.a  ../lib/libtask.a ../../bpsd/libbpsd.a ../ob/libob.a

LIBFOW = fowlib.a

# all: libs fp
all: libs fow

fplib.a: $(OBJS)
	$(LD) $(LDFLAGS) fplib.a $(OBJS)
fowlib.a: $(OBJSFOW)
	$(LD) $(LDFLAGS) fowlib.a $(OBJSFOW)

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../ob/libob.a:
	(cd ../ob; make libob.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)
	(cd ../ob; make libob.a)

fp : $(LIB_MTX) $(LIBS) fpmain.o
	$(FLINKER) fpmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

fow : $(LIB_MTX) $(LIBS) $(LIBFOW) $(OBJDIR)/fowmain.o $(OBJDIR)
	$(FLINKER) $(OBJDIR)/fowmain.o $(LIBS) $(LIBFOW) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX) $(MODINCLUDE)


libregion.so : $(LIBEQ) $(LIB_MTX) $(LIBS) $(OBJDIR)/global_region.o $(OBJDIR)/dylibregion.o
	$(FCFREE) -shared $(OBJDIR)/dylibregion.o $(LIBEQ) $(LIBS) $(OBJDIR)/global_region.o -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX) $(MODINCLUDE)

clean :
	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod
	-rm -f ./obj/*.o
	-rm -f fp fow out
	-rm -rf ./*.dSYM
	-rm -f ./csv/*

veryclean: clean
	-rm -f fp

new:
	-mkdir ../fpnew
	cp -f Makefile ../fpnew
	cp -f *.f ../fpnew
	cp -f *.inc ../fpnew

fpmain.o : fpinit.f90 fpmain.f90 fpcomm.f90
fpmenu.o : fpwrin.f90 fpgout.f90 fpfout.f90 fploop.f90 fpprep.f90 fpinit.f90 fpfile.f90 fpmenu.f90 fpcomm.f90
fpinit.o : fpinit.f90 fpcomm.f90
fpprep.o : fpoutdata.f90 fpcale.f90 fpbounce.f90 fpnfrr.f90 fpwmin.f90 fpwrin.f90 fpsave.f90 fpcalw.f90 fpbroadcast.f90 fpinit.f90 fpcoef.f90 fpdisrupt.f90 eg_read.f90 sv_read.f90 fpprep.f90 fpcomm.f90
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
fpgout.o : fpfout.f90 fpcont.f90 fpgout.f90 fpcomm.f90
fpfile.o : fpprep.f90 eg_read.f90 sv_read.f90 fpfile.f90 fpcomm.f90
fpsave.o : fpexec.f90 eg_read.f90 fpsave.f90 fpcomm.f90
fpwrin.o : fpwrin.f90 fpcomm.f90
fpwmin.o : fpwrin.f90 fpcomm.f90
fpnfrr.o : fpinit.f90 fpnfrr.f90 fpcomm.f90
fpbounce.o : fpinit.f90 fpbounce.f90 fpcomm.f90
fpcale.o: fpcale.f90 fpmpi.f90
fpdisrupt.o: fpcale.f90 fpcoef.f90 fpdisrupt.f90 fpcomm.f90
fpcalr.o: fpcalr.f90 fpcomm.f90
eg_read.o: eg_read.f90 fpcomm.f90
sv_read.o: sv_read.f90 fpcomm.f90
fpoutdata.o: fpoutdata.f90 eg_read.f90 fpcomm.f90
fpdifcoef.o: fpdifcoef.f90 fpcomm.f90 fpwrite.f90
fprad.o: fprad.f90 fpcomm.f90
fpwrite.o : fpcomm.f90

fowmain.o : fowmain.f90 fowcomm.f90 fowprep.f90 foworbit.f90\
						fpprep.f90 fpinit.f90 fpcomm.f90 fpwrite.f90 foworbitclassify.f90
fowcomm.o : fowcomm.f90
fowprep.o : fowprep.f90 fowcomm.f90
foworbit.o : foworbit.f90 fowcomm.f90 fpcomm.f90 foworbitclassify.f90 \
			../ob/obcomm.f90 ../ob/obcalc.f90 ../ob/obprep.f90 ../obinit.f90 ../ob/obparm.f90
foworbitclassify.o : foworbitclassify.f90 fowcomm.f90 fpcomm.f90
fowdistribution.o : fowdistribution.f90 foworbit.f90 foworbitclassify.f90 fowcomm.f90 fpcomm.f90
fowbounceaverage.o : fowbounceaverage.f90 fowcomm.f90 fpcomm.f90

# foweq.o : foweq.f90 fpwrite.f90
# fowjacobi.o : fowjacobi.f90 fowcomm.f90 fpcomm.f90
