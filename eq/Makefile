### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCC = eqinit.f eqcalc.f eqfunc.f eqcalq.f eqload.f eqgout.f eqgsub.f
SRCS = eqmain.f $(SRCC)
SRCT = trtest.f treqin.f $(SRCC)
SRCP = pltest.f pleqin.f eqload.f eqcalq.f
SRCA = aktest.f akeqin.f eqload.f eqcalq.f

SRCCL= ../lib/libbnd.f ../lib/libspl.f ../lib/libbrn.f ../lib/libgut.f
SRCQL= ../lib/libspl.f ../lib/libbrn.f ../lib/libgut.f

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
pl : $(OBJP)
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
	cp -f *.h ../eqnew
	-mkdir ../eqnew/in
	cp -f in/* ../eqnew/in

eqmain.o: eqmain.f eqcom1.h
eqinit.o: eqinit.f eqcom1.h
eqcalc.o: eqcalc.f eqcom1.h eqcom2.h eqcom4.h
eqcalq.o: eqcalq.f eqcom1.h eqcom3.h
eqgout.o: eqgout.f eqcom1.h eqcom2.h eqcom3.h
eqfunc.o: eqfunc.f eqcom1.h eqcom2.h eqcom4.h
trtest.o: trtest.f
treqin.o: treqin.f eqcom1.h eqcom4.h
pltest.o: pltest.f eqcom1.h eqcom3.h
pleqin.o: pleqin.f eqcom1.h eqcom3.h
aktest.o: aktest.f eqcom1.h eqcom3.h
akeqin.o: akeqin.f eqcom1.h eqcom3.h
eqload.o: eqload.f eqcom1.h eqcom3.h
