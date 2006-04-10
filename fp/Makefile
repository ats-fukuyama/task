### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = fpmenu.f fpinit.f fploop.f fpexec.f fpcoef.f fpcalv.f fpcaln.f \
       fpsave.f fpgraf.f fpfile.f fpcalw.f fpcals.f fpwrin.f fpfout.f \
       fpwmin.f

OBJS = $(SRCS:.f=.o)

LIBS = fplib.a ../pl/pllib.a ../eq/eqlib.a ../lib/tasklib.a

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all : fp

fplib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o fplib.a $(OBJS)

../pl/pllib.a:
	(cd ../pl; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../lib/tasklib.a:
	(cd ../lib; make tasklib.a)

libs:
	(cd ../pl; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../lib; make tasklib.a)

fp : fpmain.o $(LIBS)
	$(FC) fpmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

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
fpcalv.o : fpcalv.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpcaln.o : fpcaln.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpsave.o : fpsave.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpgraf.o : fpgraf.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpfout.o : fpfout.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpfile.o : fpfile.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpcalw.o : fpcalw.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc
fpcals.o : fpcals.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc ../wr/wrcom1.inc fpcom2.inc
fpwrin.o : fpwrin.f fpcomm.inc fpcom1.inc ../pl/plcom1.inc ../wr/wrcom1.inc fpcom2.inc
