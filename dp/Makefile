### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = dpmain.f dpinit.f dproot.f dpcont.f dpdisp.f dptnsr.f\
       dpfpin.f dphotf.f dphotr.f

OBJS = $(SRCS:.f=.o) ../pl/plinit.o ../pl/plprof.o \
       ../eq/eqload.o ../eq/pleqin.o ../eq/eqcalq.o \
       ../lib/libdsp.o ../lib/libbes.o ../lib/libbrn.o \
       ../lib/libspl.o ../lib/libgut.o

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
	cp -f *.h ../dpnew

dpmain.o: dpmain.f dpcomm.h dpcom1.h ../pl/plcom1.h
dpinit.o: dpinit.f dpcomm.h dpcom1.h ../pl/plcom1.h
dproot.o: dproot.f dpcomm.h dpcom1.h ../pl/plcom1.h
dpcont.o: dpcont.f dpcomm.h dpcom1.h ../pl/plcom1.h
dpdisp.o: dpdisp.f dpcomm.h dpcom1.h ../pl/plcom1.h ../pl/plcom2.h
dptnsr.o: dptnsr.f dpcomm.h dpcom1.h ../pl/plcom1.h ../pl/plcom2.h
dpfpin.o: dpfpin.f dpcomm.h dpcom1.h ../pl/plcom1.h
dphotf.o: dphotf.f dpcomm.h dpcom1.h ../pl/plcom1.h
dphotr.o: dphotr.f dpcomm.h dpcom1.h ../pl/plcom1.h
