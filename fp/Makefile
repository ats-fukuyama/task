### $Id$ ###
include ../../task/make.header

##### Without PETSc library #####
#LIBMTX=../mtxp/libmtxbnd.o ../mtxp/libmtxnompi.o
#PETSCLIBS=
##### With MUMPS library #####
#include ${PETSC_DIR}/conf/variables
LIBMTX=../mtxp/libmtxdmumps.o ../mtxp/libmtxmpi.o
PETSCLIBS=-L${PETSC_DIR}/${PETSC_ARCH}/lib -ldmumps -lmumps_common -lpord -lscalapack -lblacs -lblas -lparmetis -lmetis
##### With PETSc library #####
#include ${PETSC_DIR}/conf/variables
#LIBMTX=../mtxp/libmtxksp.o ../mtxp/libmtxmpi.o
#PETSCLIBS=${PETSC_KSP_LIB}

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE= -I.$(MOD) -I../../task/lib/$(MOD) \
            -I../../task/plx/$(MOD) -I../../task/eq/$(MOD) \
            -I../../bpsd/$(MOD) -I../../task/mtxp

SRCS = fpcomm.f90 fpinit.f90 fpmpi.f90 fpbroadcast.f90 fpexec.f90 \
       fpdrexec.f90 fpsave.f90 \
       fpwmin.f90 fpwrin.f90 fpcalwm.f90 fpcalwr.f90 fpcalw.f90 \
       fpcalcn.f90 fpcalcnr.f90 fpcalc.f90 fpnfrr.f90 fpcoef.f90 \
       fpbounce.f90 fpprep.f90 fploop.f90 \
       fpcont.f90 fpfout.f90 fpgout.f90 fpfile.f90 fpmenu.f90

OBJS = $(SRCS:.f90=.o)

LIBS =  fplib.a ../../task/plx/pllib.a ../../task/eq/eqlib.a   \
        ../../task/lib/tasklib.a ../../bpsd/bpsdlib.a 

.f.o :
	$(MF90) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(MF90) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: libs fp

fplib.a: $(OBJS)
	$(LD) $(LDFLAGS) fplib.a $(OBJS)

../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../../task/plx/pllib.a:
	(cd ../plx; make pllib.a)
../../task/eq/eqlib.a:
	(cd ../../task/eq; make eqlib.a)
../../task/mtxp/libmtxksp.o:
	(cd ../../task/mtxp; make libmtxksp.o)
../../task/mtxp/libmtxbnd.o:
	(cd ../../task/mtxp; make libmtxbnd.o)
../../task/mtxp/libmtxmpi.o:
	(cd ../../task/mtxp; make libmtxmpi.o)
../../task/mtxp/libmtxnompi.o:
	(cd ../../task/mtxp; make libmtxnompi.o)
../../task/lib/tasklib.a:
	(cd ../../task/lib; make tasklib.a)

libs:
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../../task/lib; make tasklib.a)
	(cd ../../task/mtxp; make libmtxmpi.o libmtxnompi.o)
	(cd ../../task/mtxp; make libmtxbnd.o libmtxksp.o libmtxdmumps.o)
	(cd ../../task/eq; make eqlib.a)
	(cd ../../task/plx; make pllib.a)

fp : $(LIBS) fpmain.o 
	$(MF90) fpmain.o -o $@ $(FFLAGS) $(LIBS) $(LIBMTX) $(PETSCLIBS) $(FLIBS) 

clean :
	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod

veryclean: clean
	-rm -f fp

new:
	-mkdir ../fpnew
	cp -f Makefile ../fpnew
	cp -f *.f ../fpnew
	cp -f *.inc ../fpnew

fpmain.o : fpmain.f90 fpcomm.f90
fpmenu.o : fpmenu.f90 fpcomm.f90
fpinit.o : fpinit.f90 fpcomm.f90
fpprep.o : fpprep.f90 fpcomm.f90
fploop.o : fploop.f90 fpcomm.f90
fpexec.o : fpexec.f90 fpcomm.f90
fpdrexec.o : fpdrexec.f90 fpcomm.f90
fpbroadcast.o : fpbroadcast.f90 fpcomm.f90
fpcoef.o : fpcoef.f90 fpcomm.f90
fpcalc.o : fpcalc.f90 fpcomm.f90
fpcalcn.o : fpcalcn.f90 fpcomm.f90
fpcalcnr.o : fpcalcnr.f90 fpcomm.f90
fpcalw.o : fpcalw.f90 fpcomm.f90
fpcalwr.o : fpcalwr.f90 fpcomm.f90
fpcalwm.o : fpcalwm.f90 fpcomm.f90
fpcont.o : fpcont.f90 fpcomm.f90
fpfout.o : fpfout.f90 fpcomm.f90
fpgout.o : fpgout.f90 fpcomm.f90
fpfile.o : fpfile.f90 fpcomm.f90
fpsave.o : fpsave.f90 fpcomm.f90
fpwrin.o : fpwrin.f90 fpcomm.f90
fpwmin.o : fpwrin.f90 fpcomm.f90
fpbounce.o : fpbounce.f90 fpcomm.f90
