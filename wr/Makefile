### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = wrmain.f wrinit.f wrcalc.f wrbeam.f wrgout.f

OBJS = $(SRCS:.f=.o) ../pl/plinit.o ../pl/plprof.o \
       ../eq/eqload.o ../eq/pleqin.o ../eq/eqcalq.o \
       ../dp/dpinit.o ../dp/dproot.o ../dp/dpcont.o ../dp/dpdisp.o \
       ../dp/dpfpin.o ../dp/dphotf.o ../dp/dphotr.o \
       ../lib/libdsp.o ../lib/libbes.o ../lib/libbrn.o \
       ../lib/libspl.o ../lib/libgut.o \
       ../lib/libode.o ../lib/librkf.o ../lib/librk.o

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

wr : $(OBJS)
	$(FC) $(OBJS) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core a.out *.o ./*~

veryclean: clean
	-rm -f wr

new:
	-mkdir ../wrnew
	cp Makefile ../wrnew
	cp *.f ../wrnew
	cp *.h ../wrnew

wrmain.o: wrmain.f wrcomm.h wrcom1.h ../dp/dpcom1.h ../pl/plcom1.h
wrinit.o: wrinit.f wrcomm.h wrcom1.h ../dp/dpcom1.h ../pl/plcom1.h
wrcalc.o: wrcalc.f wrcomm.h wrcom1.h ../dp/dpcom1.h ../pl/plcom1.h ../pl/plcom2.h
wrbeam.o: wrcalc.f wrcomm.h wrcom1.h ../dp/dpcom1.h ../pl/plcom1.h ../pl/plcom2.h
wrgout.o: wrgout.f wrcomm.h wrcom1.h ../dp/dpcom1.h ../pl/plcom1.h
