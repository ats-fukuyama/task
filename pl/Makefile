### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = plmain.f plinit.f plprof.f

OBJS = $(SRCS:.f=.o) ../eq/eqload.o ../eq/pleqin.o ../eq/eqcalq.o ../eq/eqdsk.o \
	../lib/libspl.o ../lib/libbrn.o ../lib/libgut.o ../lib/libchar.o ../lib/libfio.o

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
	cp -f *.inc ../plnew

plmain.o: plmain.f plcomm.inc plcom1.inc
plinit.o: plinit.f plcomm.inc plcom1.inc
plprof.o: plprof.f plcomm.inc plcom1.inc plcom2.inc
../eq/eqload.o: ../eq/eqload.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../eq/eqcalq.o: ../eq/eqcalq.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../eq/eqdsk.o: ../eq/eqdsk.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../eq/pleqin.o: ../eq/pleqin.f ../eq/eqcom1.inc ../eq/eqcom3.inc
