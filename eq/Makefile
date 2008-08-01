### $Id$ ###
include ../make.header

#FFLAGS=$(DFLAGS)
FFLAGS=$(OFLAGS)

MODINCLUDE=-I.$(MOD) -I../../bpsd$(MOD) -I../pl$(MOD) -I../lib$(MOD)

SRCS = eq_bpsd.f eqinit.f eqmenu.f eqcalc.f eqcalq.f eqcalv.f eqcalx.f \
       eqsub.f eqfunc.f eqintf.f eqsplf.f eqfile.f \
       eqgout.f eqgoutx.f eqgsub.f \
       newton.f invematrix.f equnit.f
SRCT = trtest.f treqin.f
SRCP = pltest.f
SRCA = aktest.f akeqin.f

OBJS = $(SRCS:.f=.o)
OBJT = $(SRCT:.f=.o)
OBJP = $(SRCP:.f=.o)
OBJA = $(SRCA:.f=.o)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : eq

eqlib.a: $(OBJS)
	$(LD) $(LDFLAGS) eqlib.a $(OBJS)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../pl; make pllib.a)

eq : eqmain.o eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a
	$(FCFIXED) eqmain.o eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a -o $@ $(FFLAGS) $(FLIBS)
tr : $(OBJT) eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a
	$(FCFIXED) $(OBJT) eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a -o $@ $(FFLAGS) $(FLIBS)
pl : $(OGJP) eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a
	$(FCFIXED) $(OBJP) eqlib.a ../lib/tasklib.a ../../bspd/bpsdlib.a ../pl/pllib.a -o $@ $(FFLAGS) $(FLIBS)
ak : $(OBJA) eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a
	$(FCFIXED) $(OBJA) eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less
checktr :
	ftnchek $(FCKFLAGS) $(SRCT) | less
checkpl :
	ftnchek $(FCKFLAGS) $(SRCP) | less
checkak :
	ftnchek $(FCKFLAGS) $(SRCA) | less

clean:
	-rm -f core a.out *.o *.mod ./*~ *.a .$(MOD)/*.mod
	-rm -f doc/*~

veryclean: clean
	rm -f eq pl ak

new:
	-mkdir ../eqnew
	cp -f Makefile ../eqnew
	cp -f *.f ../eqnew
	cp -f *.inc ../eqnew

EQCOMM = ../pl/plcom0.inc ../pl/plcom1.inc eqcom0.inc eqcom1.inc

eqmain.o: eqmain.f $(EQCOMM)
eqinit.o: eqinit.f $(EQCOMM)
eqcalc.o: eqcalc.f $(EQCOMM)
eqcalx.o: eqcalx.f $(EQCOMM) eqcom5.inc
eqcalq.o: eqcalq.f $(EQCOMM) eqcom3.inc
eqcalv.o: eqcalv.f $(EQCOMM) eqcom2.inc
eqsub.o: eqsub.f   $(EQCOMM) eqcom2.inc
eqfunc.o: eqfunc.f $(EQCOMM) eqcom2.inc eqcom4.inc
eqsplf.o: eqsplf.f $(EQCOMM) eqcom3.inc
eqfile.o: eqfile.f $(EQCOMM) eqcom2.inc
eqintf.o: eqintf.f $(EQCOMM) eqcom2.inc
eqgout.o: eqgout.f $(EQCOMM) eqcom2.inc eqcom3.inc
trtest.o: trtest.f
treqin.o: treqin.f $(EQCOMM) eqcom3.inc eqcom4.inc
wmeqin.o: wmeqin.f $(EQCOMM) eqcom3.inc
pltest.o: pltest.f
aktest.o: aktest.f
akeqin.o: akeqin.f $(EQCOMM) eqcom3.inc
newton.o: newton.f
invematrix.o: invematrix.f
