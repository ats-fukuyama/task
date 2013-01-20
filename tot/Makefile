### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

# without eq
#MODINCLUDE= -I$(MOD) -I../../bpsd/$(MOD) -I../plx/$(MOD)  -I../fp/$(MOD)
#EQLIBS= ../plx/noeqlib.o ../plx/noequlib.o 

# with task/eq
MODINCLUDE= -I$(MOD) -I../../bpsd/$(MOD) -I../eq/$(MOD) -I../plx/$(MOD) -I../fp/$(MOD)
EQLIBS= ../plx/noequlib.o ../eq/eqlib.a 

# with topics/equ
#MODINCLUDE= -I$(MOD) -I../../bpsd/$(MOD) -I../plx/$(MOD) -I../fp/$(MOD) -I../../topics/src/$(MOD)
#EQLIBS= ../../topics/src/tpxlib.a ../../topics/src/eqlib.a  ../plx/noeqlib.o 

# with task/eq and topics/equ
#MODINCLUDE= -I./$(MOD) -I../../bpsd/$(MOD) -I../eq/$(MOD) -I../plx/$(MOD) -I../fp/$(MOD) -I../../topics/src/$(MOD)
#EQLIBS= ../eq/eqlib.a ../../topics/src/tpxlib.a ../../topics/src/eqlib.a

##### Without PETSc library #####
LIBMTX=../mtxp/libmtxbnd.o ../mtxp/libmtxnompi.o
PETSCLIBS=
##### With PETSc library #####
#include ${PETSC_DIR}/conf/variables
#LIBMTX=../mtxp/libmtxksp.o ../mtxp/libmtxmpi.o
#PETSCLIBS=${PETSC_KSP_LIB}

SRCS = main.f

OBJS = $(SRCS:.f=.o)

LIB0 =  ../eq/treqin.o ../tr/trcom2.o ../tr/trcom3.o ../tr/tr2lib.a \
       ../dp/dplib.a ../wr/wrlib.a \
       ../eq/wmeqin.o \
	../fp/fplib.a \
       ../../bpsd/bpsdlib.a \
       ../lib/mdslib.a $(LIBMTX) $(PETSCLIBS)

LIBS = $(LIB0) ../wm/wmlib.a ../plx/pllib.a \
       $(EQLIBS) \
       ../mpi/nompilib.o ../mtx/mtxlib.a \
       ../lib/tasklib.a

LIBM = $(LIB0) ../wm/mpiwmlib.a  ../plx/pllib.a \
       $(EQLIBS) \
       ../mpi/mpilib.o ../mtx/mpimtxlib.a \
       ../lib/tasklib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : libs tot

../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
../plx/noeqlib.o:
	(cd ../plx; make noeqlib.o)
../plx/noequlib.o:
	(cd ../plx; make noequlib.o)
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
	(cd ../plx; make pllib.a)
	(cd ../plx; make noeqlib.o)
	(cd ../plx; make noequlib.o)
	(cd ../eq; make eqlib.a)
	(cd ../eq; make treqin.o)
	(cd ../eq; make wmeqin.o)
	(cd ../tr; make trcom2.o)
	(cd ../tr; make trcom3.o)
	(cd ../tr; make tr2lib.a)
	(cd ../dp; make dplib.a)
	(cd ../wr; make wrlib.a)
	(cd ../wm; make wmlib.a)
#	(cd ../wm; make mpiwmlib.a)
	(cd ../fp; make fplib.a)
	(cd ../mtx; make mtxlib.a)
#	(cd ../mtx; make mpimtxlib.a)
	(cd ../mpi; make nompilib.o)
#	(cd ../mpi; make mpilib.o)

tot : main.o $(LIBS)
	$(FCFIXED) main.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MODINCLUDE)

mpitot : main.o $(LIBM)
	$(MFC) main.o $(LIBM) -o $@ $(FFLAGS) $(FLIBS) $(MODINCLUDE)

xtot :
	(cd ../plx; make pllib.a noeqlib.o)
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
	-rm -f core a.out *.o *.mod ./*~ *.a $(MOD)/*.mod

veryclean: clean
	-rm -f tot

new:
	-mkdir ../totnew
	cp Makefile ../totnew
	cp *.f ../totnew

main.o: main.f ../wm/wmcomm.inc

