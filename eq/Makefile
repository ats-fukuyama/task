### $Id$ ###
include ../make.header

#FFLAGS=-g $(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE=-I$(MOD) -I../../bpsd/$(MOD) -I../plx/$(MOD) -I../lib/$(MOD)

SRCX = equcom.f90 equread.f90
SRCS = eqbpsd.f eqinit.f eqmenu.f eqcalc.f eqcalq.f eqcalv.f \
       eqsub.f eqfunc.f eqintf.f eqsplf.f \
       eqf-efit.f eqf-qst.f eqfile.f \
       eqgout.f eqgsub.f \
       newton.f invematrix.f equnit.f eqrppl.f
SRCT = trtest.f treqin.f
SRCP = pltest.f
SRCA = aktest.f akeqin.f
SRC2D= eqg2d.f
SRC3D= eqg3d.f

OBJS = $(SRCX:.f90=.o) $(SRCS:.f=.o)
OBJT = $(SRCT:.f=.o)
OBJP = $(SRCP:.f=.o)
OBJA = $(SRCA:.f=.o)
OBJ2D=$(SRC2D:.f=.o)
OBJ3D=$(SRC3D:.f=.o)

LIBS2 = ../plx/pllib.a eqlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a
LIBS3 = ../plx/pllib.a eq3lib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a

.SUFFIXES: .f90

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : eq


../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
eqlib.a: $(OBJS) $(OBJ2D)
	(cd ../plx; make plcomm.o)
	$(LD) $(LDFLAGS) eqlib.a $(OBJS) $(OBJ2D)
eq3lib.a: $(OBJS) $(OBJ3D)
	$(LD) $(LDFLAGS) eq3lib.a $(OBJS) $(OBJ3D)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../plx; make pllib.a)

eq: $(LIBS2) eqmain.o
	$(FCFIXED) eqmain.o $(LIBS2) -o $@ $(FFLAGS) $(FLIBS)
eq3: $(LIBS3) eqmain.o
	$(FCFIXED) eqmain.o $(LIBS3) -o $@ $(FFLAGS) $(FLIBS)
tr : $(LIBS2) $(OBJT)
	$(FCFIXED) $(OBJT) $(LIBS2) -o $@ $(FFLAGS) $(FLIBS)
pl : $(LIBS2) $(OGJP)
	$(FCFIXED) $(OBJP) $(LIBS2) -o $@ $(FFLAGS) $(FLIBS)
ak : $(LIBS2) $(OBJA)
	$(FCFIXED) $(OBJA) $(LIBS2) -o $@ $(FFLAGS) $(FLIBS)

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

EQCOMM = eqcom0.inc eqcom1.inc

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
eqbpsd.o: eqbpsd.f $(EQCOMM) eqcom3.inc
eqjft2m.o: eqjft2m.f $(EQCOMM) eqcom3.inc
equnit.o: equnit.f eqbpsd.f
eqgsub.o: eqgsub.f
eqg2d.o: eqg2d.f $(EQCOMM) eqcom3.inc
eqg3d.o: eqg3d.f $(EQCOMM) eqcom3.inc
equcom.o: equcom.f90
equread.o: equread.f90 equcom.f90
equintf.o: equintf.f $(EQCOMM) eqcom3.inc
eqrppl.o: eqrppl.f $(EQCOMM) eqcom3.inc
trtest.o: trtest.f
treqin.o: treqin.f $(EQCOMM) eqcom3.inc eqcom4.inc
pltest.o: pltest.f
aktest.o: aktest.f
akeqin.o: akeqin.f $(EQCOMM) eqcom3.inc
newton.o: newton.f
invematrix.o: invematrix.f
