### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = plmenu.f plinit.f plprof.f pldata.f

OBJS = $(SRCS:.f=.o) 

LIBS = pllib.a ../eq/pleqin.o ../eq/eqlib.a ../lib/tasklib.a

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all : pl
pllib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o pllib.a $(OBJS)

../eq/pleqin.o:
	(cd ../eq; make pleqin.o)
../lib/tasklib.a:
	(cd ../lib; make)

pl : plmain.o $(LIBS)
	$(FC) plmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core a.out *.o ./*~ ./#* *.a

veryclean: clean
	-rm -f pl

new:
	-mkdir ../plnew
	cp -f Makefile ../plnew
	cp -f *.f ../plnew
	cp -f *.inc ../plnew

plmain.o: plmain.f plcomm.inc plcom1.inc
plmenu.o: plmenu.f plcomm.inc plcom1.inc
plinit.o: plinit.f plcomm.inc plcom1.inc
plprof.o: plprof.f plcomm.inc plcom1.inc plcom2.inc
pldata.o: pldata.f plcomd.inc plcom3.inc
../eq/pleqin.o: ../eq/pleqin.f ../eq/eqcom1.inc ../eq/eqcom3.inc
