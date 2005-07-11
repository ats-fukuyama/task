### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = eqinit.f eqmenu.f eqcalc.f eqcalq.f eqcalv.f eqsub.f eqfunc.f \
       eqsplf.f eqfile.f eqgout.f eqgsub.f eqlib.f
SRCT = trtest.f treqin.f
SRCP = pltest.f pleqin.f
SRCA = aktest.f akeqin.f

OBJS = $(SRCS:.f=.o)
OBJT = $(SRCT:.f=.o)
OBJP = $(SRCP:.f=.o)
OBJA = $(SRCA:.f=.o)

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all : eq

eqlib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o eqlib.a $(OBJS)

../lib/tasklib.a: 
	(cd ../lib; make tasklib.a)
libs:
	(cd ../lib; make tasklib.a)

eq : eqmain.o eqlib.a ../lib/tasklib.a
	$(FC) eqmain.o eqlib.a ../lib/tasklib.a -o $@ $(FFLAGS) $(FLIBS)
tr : $(OBJT) eqlib.a ../lib/tasklib.a
	$(FC) $(OBJT) eqlib.a ../lib/tasklib.a -o $@ $(FFLAGS) $(FLIBS)
pl : $(OGJP) eqlib.a ../lib/tasklib.a
	$(FC) $(OBJP) eqlib.a ../lib/tasklib.a -o $@ $(FFLAGS) $(FLIBS)
ak : $(OBJA) eqlib.a ../lib/tasklib.a
	$(FC) $(OBJA) eqlib.a ../lib/tasklib.a -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less
checktr :
	ftnchek $(FCKFLAGS) $(SRCT) | less
checkpl :
	ftnchek $(FCKFLAGS) $(SRCP) | less
checkak :
	ftnchek $(FCKFLAGS) $(SRCA) | less

clean:
	rm -f core a.out *.o ./*~ *.a
	rm -f doc/*~

veryclean: clean
	rm -f eq pl ak

new:
	-mkdir ../eqnew
	cp -f Makefile ../eqnew
	cp -f *.f ../eqnew
	cp -f *.inc ../eqnew

eqmain.o: eqmain.f eqcom1.inc
eqinit.o: eqinit.f eqcom1.inc
eqcalc.o: eqcalc.f eqcom1.inc eqcom2.inc
eqcalq.o: eqcalq.f eqcom1.inc eqcom3.inc
eqcalv.o: eqcalv.f eqcom1.inc eqcom2.inc
eqsub.o: eqsub.f eqcom1.inc eqcom3.inc
eqgout.o: eqgout.f eqcom1.inc eqcom2.inc eqcom3.inc
eqfunc.o: eqfunc.f eqcom1.inc eqcom2.inc eqcom4.inc
eqsplf.o: eqsplf.f eqcom1.inc eqcom2.inc eqcom4.inc
eqfile.o: eqfile.f eqcom1.inc eqcom2.inc
trtest.o: trtest.f
treqin.o: treqin.f eqcom1.inc eqcom3.inc eqcom4.inc
pltest.o: pltest.f
pleqin.o: pleqin.f eqcom1.inc eqcom3.inc
aktest.o: aktest.f
akeqin.o: akeqin.f eqcom1.inc eqcom3.inc
