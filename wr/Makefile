### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

MODINCLUDE= -I./$(MOD) -I../../bpsd/$(MOD) -I../lib/$(MOD) -I../plx/$(MOD)

SRCS = wrmenu.f wrinit.f wrcalc.f wrbeam.f wrgout.f

OBJS = $(SRCS:.f=.o)

LIBS = wrlib.a ../plx/pllib.a ../eq/eqlib.a ../dp/dplib.a \
       ../../bpsd/bpsdlib.a ../lib/tasklib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : libs wr

wrlib.a: $(OBJS)
	$(LD) $(LDFLAGS) wrlib.a $(OBJS)

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../dp/dplib.a:
	(cd ../dp; make dplib.a)
libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../plx; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../dp; make dplib.a)

wr : libs wrmain.o
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

WRCOMM=wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../plx/plcomm.f90
wrmain.o: wrmain.f $(WRCOMM)
wrmenu.o: wrmenu.f $(WRCOMM)
wrinit.o: wrinit.f $(WRCOMM)
wrcalc.o: wrcalc.f $(WRCOMM) ../plx/plprof.f90
wrbeam.o: wrbeam.f $(WRCOMM) ../plx/plprof.f90
wrgout.o: wrgout.f $(WRCOMM)
