### $Id$ ###
include ../make.header

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE= -I.$(MOD) -I../../bpsd$(MOD) -I../pl$(MOD)

SRCS = fpmenu.f fpinit.f fploop.f fpexec.f fpcoef.f fpcalc.f fpcalcn.f \
       fpcalw.f fpcalwr.f fpcalwm.f \
       fpgraf.f fpfile.f fpfout.f fpsave.f fpwrin.f fpwmin.f

OBJS = $(SRCS:.f=.o)

LIBS = fplib.a ../pl/pllib.a ../eq/eqlib.a ../../bpsd/bpsdlib.a \
       ../lib/tasklib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : fp

fplib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o fplib.a $(OBJS)

../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../pl/pllib.a:
	(cd ../pl; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../lib/tasklib.a:
	(cd ../lib; make tasklib.a)

libs:
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../pl; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../lib; make tasklib.a)

fp : fpmain.o $(LIBS)
	$(FCFIXED) fpmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MODDIR) $(MODINCLUDE)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean :
	-rm -f core a.out *.o *.a ./#* ./*~

veryclean: clean
	-rm -f fp

new:
	-mkdir ../fpnew
	cp -f Makefile ../fpnew
	cp -f *.f ../fpnew
	cp -f *.inc ../fpnew

fpmain.o : fpmain.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpmenu.o : fpmenu.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpinit.o : fpinit.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fploop.o : fploop.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpexec.o : fpexec.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpcoef.o : fpcoef.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpcalc.o : fpcalc.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpcalcn.o : fpcalcn.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpcalw.o : fpcalw.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpcalwr.o : fpcalwr.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc ../wr/wrcom1.inc fpcom2.inc
fpcalwm.o : fpcalwm.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpgraf.o : fpgraf.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpfout.o : fpfout.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpfile.o : fpfile.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpsave.o : fpsave.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpwrin.o : fpwrin.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc ../wr/wrcom1.inc fpcom2.inc
