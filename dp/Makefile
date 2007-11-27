### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

SRCS = dpmenu.f dpinit.f dproot.f dpcont.f dpdisp.f dptnsr.f\
       dpfpin.f dphotf.f dphotr.f

OBJS = $(SRCS:.f=.o) 

LIBS = dplib.a ../pl/pllib.a  ../eq/eqlib.a ../lib/tasklib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@

all : dp

dplib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o dplib.a $(OBJS)

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

dp : dpmain.o $(LIBS)
	(cd ../pl; make pllib.a)
	(cd ../lib; make tasklib.a)
	$(FCFIXED) dpmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core a.out *.o *.a ./*~

veryclean: clean
	-rm -f dp

new:
	-mkdir ../dpnew
	cp -f Makefile ../dpnew
	cp -f *.f ../dpnew
	cp -f *.inc ../dpnew

DPCOMM=dpcomm.inc dpcom1.inc dpcom2.inc ../pl/plcom0.inc ../pl/plcom1.inc

dpmain.o: dpmain.f $(DPCOMM)
dpmenu.o: dpmenu.f $(DPCOMM)
dpinit.o: dpinit.f $(DPCOMM)
dproot.o: dproot.f $(DPCOMM)
dpcont.o: dpcont.f $(DPCOMM) ../pl/plcom2.inc
dpdisp.o: dpdisp.f $(DPCOMM) ../pl/plcom2.inc
dptnsr.o: dptnsr.f $(DPCOMM) ../pl/plcom2.inc
dpfpin.o: dpfpin.f $(DPCOMM) ../pl/plcom2.inc
dphotf.o: dphotf.f $(DPCOMM) ../pl/plcom2.inc
dphotr.o: dphotr.f $(DPCOMM) ../pl/plcom2.inc
