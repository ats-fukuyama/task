### $Id$ ###
include ../make.header

LFLAGS=$(OFLAGS)
FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

SRCC = eqinit.f eqcalc.f eqfunc.f eqcalq.f eqload.f eqgout.f eqgsub.f eqdsk.f
SRCS = eqmain.f $(SRCC)
SRCT = trtest.f treqin.f $(SRCC)
SRCP = pltest.f pleqin.f eqload.f eqcalq.f eqdsk.f
SRCA = aktest.f akeqin.f eqload.f eqcalq.f eqdsk.f

SRCQL= ../lib/libspl.f ../lib/libbrn.f ../lib/libgut.f ../lib/libchar.f ../lib/libfio.f
SRCCL= ../lib/libbnd.f $(SRCQL)

OBJS = $(SRCS:.f=.o) $(SRCCL:.f=.o) 
OBJT = $(SRCT:.f=.o) $(SRCCL:.f=.o) 
OBJP = $(SRCP:.f=.o) $(SRCQL:.f=.o) 
OBJA = $(SRCA:.f=.o) $(SRCQL:.f=.o) 

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all : eq

eq : $(OBJS)
	$(FC) $(OBJS) -o $@ $(FFLAGS) $(FLIBS)
tr : $(OBJT)
	$(FC) $(OBJT) -o $@ $(FFLAGS) $(FLIBS)
pl : $(OGJP)
	$(FC) $(OBJP) -o $@ $(FFLAGS) $(FLIBS)
ak : $(OBJA)
	$(FC) $(OBJA) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) $(SRCCL) | less
checktr :
	ftnchek $(FCKFLAGS) $(SRCT) $(SRCCL) | less
checkpl :
	ftnchek $(FCKFLAGS) $(SRCP) $(SRCTL) | less
checkak :
	ftnchek $(FCKFLAGS) $(SRCA) $(SRCTL) | less

clean:
	rm -f core a.out *.o ./*~
	rm -f doc/*~

veryclean: clean
	rm -f eq pl ak

new:
	-mkdir ../eqnew
	cp -f Makefile ../eqnew
	cp -f *.f ../eqnew
	cp -f *.inc ../eqnew
	-mkdir ../eqnew/in
	cp -f in/* ../eqnew/in

eqmain.o: eqmain.f eqcom1.inc
eqinit.o: eqinit.f eqcom1.inc
eqcalc.o: eqcalc.f eqcom1.inc eqcom2.inc
eqcalq.o: eqcalq.f eqcom1.inc eqcom3.inc
eqgout.o: eqgout.f eqcom1.inc eqcom2.inc eqcom3.inc
eqfunc.o: eqfunc.f eqcom1.inc eqcom2.inc eqcom4.inc
trtest.o: trtest.f
treqin.o: treqin.f eqcom1.inc eqcom3.inc eqcom4.inc
pltest.o: pltest.f
pleqin.o: pleqin.f eqcom1.inc eqcom3.inc
aktest.o: aktest.f
akeqin.o: akeqin.f eqcom1.inc eqcom3.inc
eqload.o: eqload.f eqcom1.inc eqcom3.inc
eqdsk.o: eqdsk.f eqcom1.inc eqcom3.inc
../lib/libbnd.o: ../lib/libbnd.f
	$(FC) -o $@ $(LFLAGS) -c ../lib/libbnd.f
../lib/libspl.o: ../lib/libspl.f
	$(FC) -o $@ $(LFLAGS) -c ../lib/libspl.f
../lib/libbrn.o: ../lib/libbrn.f
	$(FC) -o $@ $(LFLAGS) -c ../lib/libbrn.f
../lib/libgut.o: ../lib/libgut.f
	$(FC) -o $@ $(LFLAGS) -c ../lib/libgut.f
../lib/libchar.o: ../lib/libchar.f
	$(FC) -o $@ $(LFLAGS) -c ../lib/libchar.f
../lib/libfio.o: ../lib/libfio.f
	$(FC) -o $@ $(LFLAGS) -c ../lib/libfio.f
