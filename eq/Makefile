### $Id$ ###
include ../make.header

FFLAGS=$(DFLAGS)

SRCS = eqmain.f eqinit.f eqcalc.f eqfunc.f eqcalq.f eqload.f \
       eqgout.f eqgsub.f
SRCT = eqtest.f pleqin.f eqload.f eqcalq.f
SRCA = aktest.f akeqin.f eqload.f eqcalq.f

SRCSL= ../lib/libbnd.f ../lib/libspl.f ../lib/libbrn.f \
       ../lib/libgut.f
SRCTL= ../lib/libspl.f ../lib/libbrn.f ../lib/libgut.f

OBJS = $(SRCS:.f=.o) $(SRCSL:.f=.o) 
OBJT = $(SRCT:.f=.o) $(SRCTL:.f=.o) 
OBJA = $(SRCA:.f=.o) $(SRCTL:.f=.o) 

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all : eq eqtest

eq : $(OBJS)
	$(FC) $(OBJS) -o $@ $(FFLAGS) $(FLIBS)
eqtest : $(OBJT)
	$(FC) $(OBJT) -o $@ $(FFLAGS) $(FLIBS)
ak : $(OBJA)
	$(FC) $(OBJA) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) $(SRCSL) | less

checkt :
	ftnchek $(FCKFLAGS) $(SRCT) $(SRCTL) | less
checka :
	ftnchek $(FCKFLAGS) $(SRCA) $(SRCTL) | less

clean:
	rm -f core a.out *.o ./*~
	rm -f doc/*~

veryclean: clean
	rm -f eq eqtest

new:
	-mkdir ../eqnew
	cp -f Makefile ../eqnew
	cp -f *.f ../eqnew
	cp -f *.h ../eqnew
	-mkdir ../eqnew/in
	cp -f in/* ../eqnew/in

eqmain.o: eqmain.f eqcom1.h
eqinit.o: eqinit.f eqcom1.h
eqcalc.o: eqcalc.f eqcom1.h eqcom2.h
eqcalq.o: eqcalq.f eqcom1.h eqcom3.h
eqgout.o: eqgout.f eqcom1.h eqcom2.h
eqshow.o: eqshow.f eqcom1.h eqcom3.h
eqtest.o: eqtest.f eqcom1.h eqcom3.h
pleqin.o: pleqin.f eqcom1.h eqcom3.h
eqload.o: eqload.f eqcom1.h eqcom3.h
