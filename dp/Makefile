### $Id$ ###
include ../make.header

FFLAGS=$(oFLAGS)

SRCS = dpmenu.f dpinit.f dproot.f dpcont.f dpdisp.f dptnsr.f\
       dpfpin.f dphotf.f dphotr.f

OBJS = $(SRCS:.f=.o) 

LIBS = dplib.a ../pl/pllib.a  ../eq/pleqin.o ../eq/eqlib.a ../lib/tasklib.a

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all : dp

dplib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o dplib.a $(OBJS)

../pl/pllib.a:
	(cd ../pl; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../eq/pleqin.o:
	(cd ../eq; make pleqin.o)
../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
libs:
	(cd ../pl; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../eq; make pleqin.o)
	(cd ../lib; make tasklib.a)

dp : dpmain.o $(LIBS)
	(cd ../pl; make pllib.a)
	(cd ../eq; make eqlib.a pleqin.o)
	(cd ../lib; make tasklib.a)
	$(FC) dpmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

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

dpmain.o: dpmain.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dpmenu.o: dpmenu.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dpinit.o: dpinit.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dproot.o: dproot.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dpcont.o: dpcont.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dpdisp.o: dpdisp.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc ../pl/plcom2.inc
dptnsr.o: dptnsr.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc ../pl/plcom2.inc
dpfpin.o: dpfpin.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dphotf.o: dphotf.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
dphotr.o: dphotr.f dpcomm.inc dpcom1.inc ../pl/plcom1.inc
../eq/pleqin.o: ../eq/pleqin.f ../eq/eqcom1.inc ../eq/eqcom3.inc
