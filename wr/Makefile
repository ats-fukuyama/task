### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = wrmain.f wrinit.f wrcalc.f wrbeam.f wrgout.f
SRCPL= ../pl/plinit.f ../pl/plprof.f
SRCEQ= ../eq/eqload.f ../eq/pleqin.f ../eq/eqcalq.f
SRCDP= ../dp/dpinit.f ../dp/dproot.f ../dp/dpcont.f ../dp/dpdisp.f \
       ../dp/dptnsr.f ../dp/dpfpin.f ../dp/dphotf.f ../dp/dphotr.f
SRCLIB=../lib/libdsp.f ../lib/libbes.f ../lib/libbrn.f \
       ../lib/libspl.f ../lib/libgut.f \
       ../lib/libode.f ../lib/librkf.f ../lib/librk.f

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
	cp *.h ../wrnew

wrmain.o: wrmain.f wrcomm.h wrcom1.h ../dp/dpcom1.h ../pl/plcom1.h
wrinit.o: wrinit.f wrcomm.h wrcom1.h ../dp/dpcom1.h ../pl/plcom1.h
wrcalc.o: wrcalc.f wrcomm.h wrcom1.h ../dp/dpcom1.h ../pl/plcom1.h ../pl/plcom2.h
wrbeam.o: wrcalc.f wrcomm.h wrcom1.h ../dp/dpcom1.h ../pl/plcom1.h ../pl/plcom2.h
wrgout.o: wrgout.f wrcomm.h wrcom1.h ../dp/dpcom1.h ../pl/plcom1.h
