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

../pl/pllib.a:
	(cd ../pl; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../eq/pleqin.o:
	(cd ../eq; make pleqin.o)
../eq/treqin.o:
	(cd ../eq; make treqin.o)
../eq/wmeqin.o:
	(cd ../eq; make wmeqin.o)
../dp/tr2lib.a:
	(cd ../dp; make tr2lib.a)
../dp/dplib.a:
	(cd ../dp; make dplib.a)
../wr/wrlib.a:
	(cd ../wr; make wrlib.a)
../wm/wmlib.a:
	(cd ../wm; make wmlib.a)
../fp/fplib.a:
	(cd ../fp; make fplib.a)
../mtx/mtxlib.a:
	(cd ../mtx; make mtxlib.a)
../lib/tasklib.a:
	(cd ../lib; make)

tot : main.o $(LIBS)
	$(FC) main.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

force : main.o $(LIBS)
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


