### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = plmain.f plinit.f plprof.f

OBJS = $(SRCS:.f=.o) ../eq/eqload.o ../eq/pleqin.o ../eq/eqcalq.o \
	../lib/libspl.o ../lib/libbrn.o ../lib/libgut.o

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

pl : $(OBJS)
	$(FC) $(OBJS) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core a.out *.o ./*~ ./#*

veryclean: clean
	-rm -f pl

new:
	-mkdir ../plnew
	cp -f Makefile ../plnew
	cp -f *.f ../plnew
	cp -f *.h ../plnew

plmain.o: plmain.f plcomm.h plcom1.h
plinit.o: plinit.f plcomm.h plcom1.h
plprof.o: plprof.f plcomm.h plcom1.h plcom2.h
