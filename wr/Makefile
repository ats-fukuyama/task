### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = wrmenu.f wrinit.f wrcalc.f wrbeam.f wrgout.f

OBJS = $(SRCS:.f=.o)

LIBS = wrlib.a ../dp/dplib.a ../pl/pllib.a ../eq/pleqin.o ../eq/eqlib.a \
       ../lib/tasklib.a

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all : wr

wrlib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o wrlib.a $(OBJS)

wr : wrmain.o $(LIBS)
	(cd ../pl; make pllib.a)
	(cd ../eq; make eqlib.a pleqin.o)
	(cd ../dp; make dplib.a)
	(cd ../lib; make tasklib.a)
	$(FC) wrmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core a.out *.o ./*~

veryclean: clean
	-rm -f wr

new:
	-mkdir ../wrnew
	cp Makefile ../wrnew
	cp *.f ../wrnew
	cp *.inc ../wrnew

wrmain.o: wrmain.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc
wrmenu.o: wrmenu.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc
wrinit.o: wrinit.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc
wrcalc.o: wrcalc.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc ../pl/plcom2.inc
wrbeam.o: wrbeam.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc ../pl/plcom2.inc
wrgout.o: wrgout.f wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../pl/plcom1.inc
