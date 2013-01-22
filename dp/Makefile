### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

MODINCLUDE= -I $(MOD) -I../plx/$(MOD) -I../lib/$(MOD) -I../../bpsd/$(MOD)

SRCS = dpmenu.f dpinit.f dproot.f dpcont.f dpdisp.f \
       dptens.f dptnsr1.f dptnsr2.f \
       dpfpin.f dphotf.f dphotr.f dpdkdt.f dpgsub.f 

OBJS = $(SRCS:.f=.o) 

LIBS = dplib.a ../eq/eqlib.a ../plx/pllib.a ../lib/tasklib.a \
       ../../bpsd/bpsdlib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : dp

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../plx/pllib.a:
	(cd ../pl; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
dplib.a: $(OBJS)
	$(LD) $(LDFLAGS) dplib.a $(OBJS)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../plx; make pllib.a)
	(cd ../eq; make eqlib.a)

dp : dpmain.o $(LIBS)
	(cd ../plx; make pllib.a)
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

DPCOMM=dpcomm.inc dpcom1.inc dpcom2.inc ../plx/plcomm.f90

dpmain.o: dpmain.f $(DPCOMM)
dpmenu.o: dpmenu.f $(DPCOMM)
dpinit.o: dpinit.f $(DPCOMM)
dproot.o: dproot.f $(DPCOMM)
dpcont.o: dpcont.f $(DPCOMM) ../plx/plprof.f90
dpdisp.o: dpdisp.f $(DPCOMM) ../plx/plprof.f90
dptens.o: dptens.f $(DPCOMM) ../plx/plprof.f90
dptnsr1.o: dptnsr1.f $(DPCOMM) ../plx/plprof.f90
dptnsr2.o: dptnsr2.f $(DPCOMM) ../plx/plprof.f90
dpfpin.o: dpfpin.f $(DPCOMM) ../plx/plprof.f90
dphotf.o: dphotf.f $(DPCOMM) ../plx/plprof.f90
dphotr.o: dphotr.f $(DPCOMM) ../plx/plprof.f90
dpdkdt.o: dpdkdt.f $(DPCOMM) ../plx/plprof.f90
dpgsub.o: dpgsub.f
