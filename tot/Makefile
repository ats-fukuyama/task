### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = main.f

OBJS = $(SRCS:.f=.o)

LIBS = ../eq/eqlib.a \
       ../eq/pleqin.o ../pl/pllib.a \
       ../eq/treqin.o ../tr/tr2lib.a \
       ../dp/dplib.a ../wr/wrlib.a \
       ../eq/wmeqin.o ../wm/wmlib.a ../fp/fplib.a \
       ../mtx/mtxlib.a ../lib/tasklib.a $(OBJLA)

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all : tot

tot : main.o $(LIBS)
	(cd ../eq; make eqlib.a pleqin.o treqin.o wmeqin.o)
	(cd ../pl; make pllib.a)
	(cd ../tr; make tr2lib.a)
	(cd ../dp; make dplib.a)
	(cd ../wr; make wrlib.a)
	(cd ../wm; make wmlib.a)
	(cd ../fp; make fplib.a)
	(cd ../mtx; make mtxlib.a)
	$(FC) main.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

force : 
	(cd ../eq; make eqlib.a pleqin.o treqin.o wmeqin.o)
	(cd ../pl; make pllib.a)
	(cd ../tr; make tr2lib.a)
	(cd ../dp; make dplib.a)
	(cd ../wr; make wrlib.a)
	(cd ../wm; make wmlib.a)
	(cd ../fp; make fplib.a)
	(cd ../mtx; make mtxlib.a)
	$(FC) main.o $(LIBS) -o tot $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core a.out *.o ./*~

veryclean: clean
	-rm -f tot

new:
	-mkdir ../totnew
	cp Makefile ../totnew
	cp *.f ../totnew

main.o: main.f ../wm/wmcomm.inc


