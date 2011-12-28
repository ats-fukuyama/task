### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

MODINCLUDE=-I$(MOD) -I../../bpsd/$(MOD) -I../pl/$(MOD) -I../lib/$(MOD)

SRCS = eq_bpsd.f eqinit.f eqmenu.f eqcalc.f eqcalq.f eqcalv.f \
       eqsub.f eqfunc.f eqintf.f eqsplf.f eqfile.f \
       eqgout.f eqgsub.f \
       newton.f invematrix.f equnit.f eqrppl.f
SRCT = trtest.f treqin.f
SRCP = pltest.f
SRCA = aktest.f akeqin.f
SRC2D= eqg2d.f
SRC3D= eqg3d.f

OBJS = $(SRCS:.f=.o)
OBJT = $(SRCT:.f=.o)
OBJP = $(SRCP:.f=.o)
OBJA = $(SRCA:.f=.o)
OBJ2D=$(SRC2D:.f=.o)
OBJ3D=$(SRC3D:.f=.o)

LIBS2 = eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a
LIBS3 = eq3lib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a

.SUFFIXES: .f90

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : libs eq
#all : eq3

eqlib.a: $(OBJS) $(OBJ2D)
	$(LD) $(LDFLAGS) eqlib.a $(OBJS) $(OBJ2D)
eq3lib.a: $(OBJS) $(OBJ3D)
	$(LD) $(LDFLAGS) eq3lib.a $(OBJS) $(OBJ3D)

../lib/tasklib.a:
	(cd ../lib; make)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make)
../pl/pllib.a:
	(cd ../pl; make)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../pl; make pllib.a)

eq: eqmain.o eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a
	$(FCFIXED) eqmain.o $(LIBS2) -o $@ $(FFLAGS) $(FLIBS)
eq3: eqmain.o eq3lib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a
	$(FCFIXED) eqmain.o $(LIBS3) -o $@ $(FFLAGS) $(FLIBS)
tr : $(OBJT) eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a
	$(FCFIXED) $(OBJT) $(LIBS) -o $@ $(FFLAGS) $(FLIBS)
pl : $(OGJP) eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a
	$(FCFIXED) $(OBJP) $(LIBS) -o $@ $(FFLAGS) $(FLIBS)
ak : $(OBJA) eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a ../pl/pllib.a
	$(FCFIXED) $(OBJA) $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

check : $(SRCS)
	ftnchek $(FCKFLAGS) $(SRCS) $(SRC3D) | $(PAGER)
checktr : $(SRCT)
	ftnchek $(FCKFLAGS) $(SRCT) | $(PAGER)
checkpl : $(SRCP)
	ftnchek $(FCKFLAGS) $(SRCP) | $(PAGER)
checkak : $(SRCA)
	ftnchek $(FCKFLAGS) $(SRCA) | $(PAGER)

clean:
	-rm -f core a.out *.o *.mod ./*~ *.a $(MOD)/*.mod
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
eqsub.o : eqsub.f  $(EQCOMM) eqcom2.inc
eqfunc.o: eqfunc.f $(EQCOMM) eqcom2.inc eqcom4.inc
eqsplf.o: eqsplf.f $(EQCOMM) eqcom3.inc
eqfile.o: eqfile.f $(EQCOMM) eqcom2.inc
eqintf.o: eqintf.f $(EQCOMM) eqcom2.inc
eqgout.o: eqgout.f $(EQCOMM) eqcom2.inc eqcom3.inc
trtest.o: trtest.f
treqin.o: treqin.f $(EQCOMM) eqcom3.inc eqcom4.inc
pltest.o: pltest.f
aktest.o: aktest.f
akeqin.o: akeqin.f $(EQCOMM) eqcom3.inc
newton.o: newton.f
invematrix.o: invematrix.f
