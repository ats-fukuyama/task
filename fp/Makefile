include ../../task/make.header
include ../mtxp/make.mtxp
#LIBX_MTX=$(LIBX_MTX_BND)
LIBX_MTX=$(LIBX_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_KSP)

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../../task/lib/$(MOD) \
            -I../../task/plx/$(MOD) -I../../task/eq/$(MOD) \
            -I../../bpsd/$(MOD) -I../../task/mtxp/$(MOD)

SRCS = fpcomm.f90 fpinit.f90 fpmpi.f90 fpbroadcast.f90 fpexec.f90 fpsave.f90 \
       fpwmin.f90 fpwrin.f90 fpcalw.f90 fpcalwm.f90 fpcalwr.f90 \
       fpcalcn.f90 fpcalcnr.f90 fpcalc.f90 fpnfrr.f90 fpcoef.f90 fpcale.f90 \
       fpbounce.f90 fpdisrupt.f90 fpprep.f90 fploop.f90 \
       fpcont.f90 fpfout.f90 fpgout.f90 fpfile.f90 fpmenu.f90

OBJS = $(SRCS:.f90=.o)

LIBS = fplib.a ../../task/plx/pllib.a ../../task/eq/eqlib.a   \
        ../../task/lib/tasklib.a ../../bpsd/bpsdlib.a

all: libs fp

fplib.a: $(OBJS)
	$(LD) $(LDFLAGS) fplib.a $(OBJS)

../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../../task/plx/pllib.a:
	(cd ../plx; make pllib.a)
../../task/eq/eqlib.a:
	(cd ../../task/eq; make eqlib.a)
../../task/lib/tasklib.a:
	(cd ../../task/lib; make tasklib.a)

libs:
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../../task/eq; make eqlib.a)
	(cd ../../task/mtxp; make)
	(cd ../../task/plx; make pllib.a)
	(cd ../../task/lib; make tasklib.a)

fp : $(LIBS) fpmain.o 
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

fpmain.o : fpmain.f90 fpcomm.f90
fpmenu.o : fpmenu.f90 fpcomm.f90
fpinit.o : fpinit.f90 fpcomm.f90
fpprep.o : fpprep.f90 fpcomm.f90
fploop.o : fploop.f90 fpcomm.f90
fpexec.o : fpexec.f90 fpcomm.f90
fpdrexec.o : fpdrexec.f90 fpcomm.f90
fpbroadcast.o : fpbroadcast.f90 fpcomm.f90
fpcoef.o : fpcoef.f90 fpcomm.f90
fpcalc.o : fpcalc.f90 fpcomm.f90
fpcalcn.o : fpcalcn.f90 fpcomm.f90
fpcalcnr.o : fpcalcnr.f90 fpcomm.f90
fpcalw.o : fpcalw.f90 fpcomm.f90
fpcalwr.o : fpcalwr.f90 fpcomm.f90
fpcalwm.o : fpcalwm.f90 fpcomm.f90
fpcont.o : fpcont.f90 fpcomm.f90
fpfout.o : fpfout.f90 fpcomm.f90
fpgout.o : fpgout.f90 fpcomm.f90
fpfile.o : fpfile.f90 fpcomm.f90
fpsave.o : fpsave.f90 fpcomm.f90
fpwrin.o : fpwrin.f90 fpcomm.f90
fpwmin.o : fpwrin.f90 fpcomm.f90
fpbounce.o : fpbounce.f90 fpcomm.f90
fpcale.o: fpmpi.f90
fpdisrupt.o: fpdisrupt.f90 fpcomm.f90
fpmpi.o: fpmpi.f90
	$(FLINKER) $(FFLAGS) -c fpmpi.f90 $(MODDIR) $(MODINCLUDE)
