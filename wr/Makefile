### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

MODINCLUDE= -I./$(MOD) -I../../bpsd/$(MOD) -I../pl/$(MOD)

SRCS = wrmenu.f wrinit.f wrcalc.f wrbeam.f wrgout.f

OBJS = $(SRCS:.f=.o)

WMLIBS = ../wm/wmlib.a ../mtx/mtxlib.a ../eq/wmeqin.o ../mpi/nompilib.o
LIBS = wrlib.a ../pl/pllib.a ../eq/eqlib.a ../dp/dplib.a \
       ../../bpsd/bpsdlib.a ../lib/tasklib.a $(WMLIBS)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : wr

wrlib.a: $(OBJS)
	$(LD) $(LDFLAGS) wrlib.a $(OBJS)

../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../pl/pllib.a:
	(cd ../pl; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../dp/dplib.a:
	(cd ../dp; make dplib.a)
../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
libs:
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../pl; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../dp; make dplib.a)
	(cd ../lib; make tasklib.a)

wr : wrmain.o $(LIBS)
	$(FCFIXED) wrmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MODDIR) $(MODINCLUDE)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core a.out *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f wr

new:
	-mkdir ../wrnew
	cp Makefile ../wrnew
	cp *.f ../wrnew
	cp *.inc ../wrnew

wrmain.o: wrmain.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc
wrmenu.o: wrmenu.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc
wrinit.o: wrinit.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc
wrcalc.o: wrcalc.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc ../pl/plcom2.inc
wrbeam.o: wrbeam.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc ../pl/plcom2.inc
wrgout.o: wrgout.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc
