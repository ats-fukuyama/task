### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = dpmain.f dpinit.f dproot.f dpcont.f dpdisp.f dptnsr.f\
       dpfpin.f dphotf.f dphotr.f

OBJS = $(SRCS:.f=.o) ../pl/plinit.o ../pl/plprof.o \
       ../eq/eqload.o ../eq/pleqin.o ../eq/eqcalq.o ../eq/eqdsk.o \
       ../lib/libdsp.o ../lib/libbes.o ../lib/libbrn.o \
       ../lib/libspl.o ../lib/libgut.o ../lib/libchar.o ../lib/libfio.o

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

dp : $(OBJS)
	$(FC) $(OBJS) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core a.out *.o ./*~

veryclean: clean
	-rm -f dp

new:
	-mkdir ../dpnew
	cp -f Makefile ../dpnew
	cp -f *.f ../dpnew
	cp -f *.inc ../dpnew

dpmain.o: dpmain.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dpinit.o: dpinit.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dproot.o: dproot.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dpcont.o: dpcont.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dpdisp.o: dpdisp.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc ../pl/plcom2.inc
dptnsr.o: dptnsr.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc ../pl/plcom2.inc
dpfpin.o: dpfpin.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dphotf.o: dphotf.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dphotr.o: dphotr.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
../pl/plinit.o: ../pl/plinit.f ../pl/plcomm.inc ../pl/plcom1.inc
../pl/plprof.o: ../pl/plprof.f ../pl/plcomm.inc ../pl/plcom1.inc ../pl/plcom2.inc
../eq/eqload.o: ../eq/eqload.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../eq/eqcalq.o: ../eq/eqcalq.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../eq/eqdsk.o: ../eq/eqdsk.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../eq/pleqin.o: ../eq/pleqin.f ../eq/eqcom1.inc ../eq/eqcom3.inc
