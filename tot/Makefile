### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = main.f

OBJS = $(SRCS:.f=.o)

LIB0 = ../eq/eqlib.a \
       ../eq/pleqin.o ../pl/pllib.a \
       ../eq/treqin.o ../tr/tr2lib.a \
       ../dp/dplib.a ../wr/wrlib.a \
       ../eq/wmeqin.o ../fp/fplib.a \
       ../lib/tasklib.a ../lib/mdslib.a
LIBS = $(LIB0) ../wm/wmlib.a  ../mtx/mtxlib.a ../mpi/nompilib.o
LIBM = $(LIB0) ../wm/mpiwmlib.a  ../mtx/mpimtxlib.a ../mpi/mpilib.o

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
../tr/tr2lib.a:
	(cd ../tr; make tr2lib.a)
../dp/dplib.a:
	(cd ../dp; make dplib.a)
../wr/wrlib.a:
	(cd ../wr; make wrlib.a)
../wm/wmlib.a:
	(cd ../wm; make wmlib.a)
../wm/mpiwmlib.a:
	(cd ../wm; make mpiwmlib.a)
../fp/fplib.a:
	(cd ../fp; make fplib.a)
../mtx/mtxlib.a:
	(cd ../mtx; make mtxlib.a)
../mtx/mpimtxlib.a:
	(cd ../mtx; make mpimtxlib.a)
../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../lib/mdslib.a:
	(cd ../lib; make mdslib.a)
../mpi/nompilib.o:
	(cd ../mpi; make nompilib.o)
../mpi/mpilib.o:
	(cd ../mpi; make mpilib.o)

libs:
	(cd ../pl; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../eq; make pleqin.o)
	(cd ../eq; make treqin.o)
	(cd ../eq; make wmeqin.o)
	(cd ../tr; make tr2lib.a)
	(cd ../dp; make dplib.a)
	(cd ../wr; make wrlib.a)
	(cd ../wm; make wmlib.a)
	(cd ../wm; make mpiwmlib.a)
	(cd ../fp; make fplib.a)
	(cd ../mtx; make mtxlib.a)
	(cd ../mtx; make mpimtxlib.a)
	(cd ../lib; make tasklib.a)
	(cd ../mpi; make nompilib.o)
	(cd ../mpi; make mpilib.o)

tot : main.o $(LIBS)
	$(FC) main.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(OBJLA)

mpitot : main.o $(LIBM)
	$(MFC) main.o $(LIBM) -o $@ $(FFLAGS) $(FLIBS) $(OBJLA)

xtot :
	(cd ../pl; make pllib.a)
	(cd ../eq; make eqlib.a pleqin.o treqin.o wmeqin.o)
	(cd ../tr; make tr2lib.a)
	(cd ../dp; make dplib.a)
	(cd ../wr; make wrlib.a)
	(cd ../wm; make wmlib.a)
	(cd ../fp; make fplib.a)
	(cd ../mtx; make mtxlib.a)
	(cd ../lib; make)
	(cd ../mpi; make nompilib.o)
	$(FC) main.o $(LIBS) -o tot $(FFLAGS) $(FLIBS)

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


