### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = eqmain.f eqinit.f eqcalc.f eqcalq.f eqload.f eqgout.f eqshow.f eqgraf.f
SRCT = eqtest.f eqload.f eqcalq.f eqread.f

OBJS = $(SRCS:.f=.o) ../lib/libbnd.o ../lib/libspl.o ../lib/libbrn.o \
       ../lib/libgut.o
OBJT = $(SRCT:.f=.o) ../lib/libspl.o ../lib/libbrn.o ../lib/libgut.o

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all : eq eqtest

eq : $(OBJS)
	$(FC) $(OBJS) -o $@ $(FFLAGS) $(FLIBS)
eqtest : $(OBJT)
	$(FC) $(OBJT) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

checkt :
	ftnchek $(FCKFLAGS) $(SRCT) | less

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
eqread.o: eqread.f eqcom1.h eqcom3.h
eqload.o: eqload.f eqcom1.h eqcom3.h
