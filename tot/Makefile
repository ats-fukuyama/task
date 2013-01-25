### $Id$ ###
include ../make.header

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE= -I$(MOD) -I../../bpsd/$(MOD) -I../eq/$(MOD) -I../plx/$(MOD) -I../fp/$(MOD) -I../mtxp
EQLIBS= ../plx/noequlib.o ../eq/eqlib.a 

##### Without PETSc library #####
LIBMTXP=../mtxp/libmtxbnd.o ../mtxp/libmtxnompi.o
PETSCLIBS=
##### With PETSc library #####
#include ${PETSC_DIR}/conf/variables
#LIBMTXP=../mtxp/libmtxksp.o ../mtxp/libmtxmpi.o
#PETSCLIBS=${PETSC_KSP_LIB}

SRCS = main.f

OBJS = $(SRCS:.f=.o)

LIBS = ../lib/mdslib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a \
       $(LIBMTXP) ../mtx/mtxlib.a ../mpi/mpilib.o \
       ../plx/pllib.a \
       ../tr/trcom2.o ../tr/trcom3.o $(EQLIBS) \
       ../tr/tr2lib.a ../dp/dplib.a ../wr/wrlib.a ../wm/wmlib.a \
       ../fp/fplib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : tot

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../lib/mdslib.a:
	(cd ../lib; make mdslib.a)
../mtxp/libmtxmpi.o:
	(cd ../mtxp; make libmtxmpi.o)
../mtxp/libmtxnompi.o:
	(cd ../mtxp; make libmtxnompi.o)
../mtxp/libmtxksp.o:
	(cd ../mtxp; make libmtxksp.o)
../mtxp/libmtxmumps.o:
	(cd ../mtxp; make libmtxmumps.o)
../mtxp/libmtxbnd.o:
	(cd ../mtxp; make libmtxbnd.o)
../mpi/mpilib.o:
	(cd ../mpi; make mpilib.o)
../mtx/mtxlib.a:
	(cd ../mtx; make mtxlib.a)
../mtx/mpimtxlib.a:
	(cd ../mtx; make mpimtxlib.a)
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
../tr/trcom2.o:
	(cd ../tr; make trcom2.o)
../tr/trcom3.o:
	(cd ../tr; make trcom3.o)
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

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../mtxp; make)
	(cd ../mtx; make mtxlib.a)
	(cd ../mpi; make mpilib.o)
	(cd ../plx; make pllib.a)
	(cd ../plx; make noeqlib.o)
	(cd ../eq; make eqlib.a)
	(cd ../tr; make trcom2.o)
	(cd ../tr; make trcom3.o)
	(cd ../tr; make tr2lib.a)
	(cd ../dp; make dplib.a)
	(cd ../wr; make wrlib.a)
	(cd ../wm; make wmlib.a)
	(cd ../fp; make fplib.a)

tot : libs main.o
	$(FCFIXED) main.o -o $@ $(FFLAGS) $(MODINCLUDE) $(LIBS) $(FLIBS) 

mpitot : $(LIBM) main.o
	$(MF90) main.o -o $@ $(FFLAGS) $(MODINCLUDE) $(LIB) $(FLIBS) $(PETSCLIBS) 

clean:
	-rm -f core a.out *.o *.mod ./*~ *.a $(MOD)/*.mod

veryclean: clean
	-rm -f tot

main.o: main.f ../wm/wmcomm.inc

