### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

# without eq
#MODINCLUDE= -I.$(MOD) -I../pl$(MOD)
#EQLIBS= ../pl/noeqlib.o ../pl/noequlib.o 

# with task/eq
MODINCLUDE= -I.$(MOD) -I../eq$(MOD) -I../pl$(MOD)
EQLIBS= ../pl/noequlib.o ../eq/eqlib.a 

# with topics/equ
#MODINCLUDE= -I.$(MOD) -I../pl$(MOD) -I../../topics/src$(MOD)
#EQLIBS= ../../topics/src/tpxlib.a ../../topics/src/eqlib.a  ../pl/noeqlib.o 

# with task/eq and topics/equ
#MODINCLUDE= -I.$(MOD) -I../eq$(MOD) -I../pl$(MOD) -I../../topics/src$(MOD)
#EQLIBS= ../eq/eqlib.a ../../topics/src/tpxlib.a ../../topics/src/eqlib.a

SRCS = main.f

OBJS = $(SRCS:.f=.o)

LIB0 = ../eq/treqin.o ../tr/tr2lib.a \
       ../dp/dplib.a ../wr/wrlib.a ../fp/fplib.a \
       ../pl/pllib.a ../eq/wmeqin.o $(EQLIBS) \
       ../../bpsd/bpsdlib.a \
       ../lib/mdslib.a 
       
LIBS = $(LIB0) ../wm/wmlib.a  ../mtx/mtxlib.a \
       ../lib/tasklib.a ../mpi/nompilib.o
LIBM = $(LIB0) ../wm/mpiwmlib.a  ../mtx/mpimtxlib.a \
       ../lib/tasklib.a ../mpi/mpilib.o

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : tot

../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../pl/pllib.a:
	(cd ../pl; make pllib.a)
../pl/noeqlib.o:
	(cd ../pl; make noeqlib.o)
../pl/noequlib.o:
	(cd ../pl; make noequlib.o)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
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
	(cd ../lib; make tasklib.a)
	(cd ../pl; make pllib.a)
	(cd ../pl; make noeqlib.o)
	(cd ../pl; make noequlib.o)
	(cd ../eq; make eqlib.a)
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
	(cd ../mpi; make nompilib.o)
	(cd ../mpi; make mpilib.o)

tot : main.o $(LIBS)
	$(FCFIXED) main.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MODINCLUDE)

mpitot : main.o $(LIBM)
	$(MFC) main.o $(LIBM) -o $@ $(FFLAGS) $(FLIBS) $(MODINCLUDE)

xtot :
	(cd ../pl; make pllib.a noeqlib.o)
	(cd ../eq; make eqlib.a treqin.o wmeqin.o)
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
	-rm -f core a.out *.o *.mod ./*~ *.a .$(MOD)/*.mod

veryclean: clean
	-rm -f tot

new:
	-mkdir ../totnew
	cp Makefile ../totnew
	cp *.f ../totnew

main.o: main.f ../wm/wmcomm.inc

