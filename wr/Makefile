### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = wrmain.f wrinit.f wrcalc.f wrbeam.f wrgout.f
SRCPL= ../pl/plinit.f ../pl/plprof.f
SRCEQ= ../eq/eqload.f ../eq/pleqin.f ../eq/eqcalq.f ../eq/eqdsk.f
SRCDP= ../dp/dpinit.f ../dp/dproot.f ../dp/dpcont.f ../dp/dpdisp.f \
       ../dp/dptnsr.f ../dp/dpfpin.f ../dp/dphotf.f ../dp/dphotr.f
SRCLIB=../lib/libdsp.f ../lib/libbes.f ../lib/libbrn.f \
       ../lib/libspl.f ../lib/libgut.f ../lib/libfp.f \
       ../lib/libode.f ../lib/librkf.f ../lib/librk.f \
       ../lib/libchar.f ../lib/libfio.f

OBJS = $(SRCS:.f=.o) $(SRCPL:.f=.o) $(SRCEQ:.f=.o) $(SRCDP:.f=.o) \
       $(SRCLIB:.f=.o)

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

wr : $(OBJS)
	$(FC) $(OBJS) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) $(SRCPL) $(SRCEQ) $(SRCDP) $(SRCLIB)| less

clean:
	-rm -f core a.out *.o ./*~

veryclean: clean
	-rm -f wr

new:
	-mkdir ../wrnew
	cp Makefile ../wrnew
	cp *.f ../wrnew
	cp *.inc ../wrnew

wrmain.o: wrmain.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc
wrinit.o: wrinit.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc
wrcalc.o: wrcalc.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc ../pl/plcom2.inc
wrbeam.o: wrbeam.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc ../pl/plcom2.inc
wrgout.o: wrgout.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc

../pl/plinit.o: ../pl/plinit.f ../pl/plcomm.inc ../pl/plcom1.inc
../pl/plprof.o: ../pl/plprof.f ../pl/plcomm.inc ../pl/plcom1.inc ../pl/plcom2.inc
../eq/eqload.o: ../eq/eqload.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../eq/eqcalq.o: ../eq/eqcalq.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../eq/eqdsk.o: ../eq/eqdsk.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../eq/pleqin.o: ../eq/pleqin.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../dp/dpinit.o: ../dp/dpinit.f ../dp/dpcomm.inc ../dp/dpcom1.inc ../pl/plcom1.inc
../dp/dproot.o: ../dp/dproot.f ../dp/dpcomm.inc ../dp/dpcom1.inc ../pl/plcom1.inc
../dp/dpcont.o: ../dp/dpcont.f ../dp/dpcomm.inc ../dp/dpcom1.inc ../pl/plcom1.inc
../dp/dpdisp.o: ../dp/dpdisp.f ../dp/dpcomm.inc ../dp/dpcom1.inc ../pl/plcom1.inc ../pl/plcom2.inc
../dp/dptnsr.o: ../dp/dptnsr.f ../dp/dpcomm.inc ../dp/dpcom1.inc ../pl/plcom1.inc ../pl/plcom2.inc
../dp/dpfpin.o: ../dp/dpfpin.f ../dp/dpcomm.inc ../dp/dpcom1.inc ../pl/plcom1.inc
../dp/dphotf.o: ../dp/dphotf.f ../dp/dpcomm.inc ../dp/dpcom1.inc ../pl/plcom1.inc
../dp/dphotr.o: ../dp/dphotr.f ../dp/dpcomm.inc ../dp/dpcom1.inc ../pl/plcom1.inc
