include ../make.header
#FCFREE=ifort 
#include ${PETSC_DIR}/conf/variables
#PETSCLIBS=${PETSC_KSP_LIB}
TARGETP = ./t2
##### Without PETSc library #####
LIBMTX=../mtxp/libmtxnompi.o ../mtxp/libmtxbnd.o
#LIBMTX=../mtxp/libmtxnompi.o ../mtxp/libmtxpcg.o
PETSCLIBS=
FLINKER=$(FCFREE)
##### With PETSc library #####
#include ${PETSC_DIR}/conf/variables
#LIBMTX= ../mtxp/libmtxmpi.o ../mtxp/libmtxksp.o
#PETSCLIBS=${PETSC_KSP_LIB}
#FLINKER=$(MF90)

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a 

#FFLAGS= $(OFLAGS)
FFLAGS= $(DFLAGS)

MODINCLUDE= -I.$(MOD) -I../mtxp -I../lib/$(MOD) -I../../bpsd/$(MOD) \
              -I../t2x/$(MOD) -I../plx/$(MOD) -I../eq/$(MOD)


SRCS = libt2.f90 t2cnst.f90 t2comm.f90 t2intg.f90 t2ngra.f90 t2vgra.f90\
       t2init.f90 t2parm.f90 t2div.f90 t2prof.f90 t2writ.f90 \
       t2prep.f90 t2calvs.f90 t2exec.f90 t2save.f90 \
       t2step.f90 t2loop.f90 t2gout.f90 t2menu.f90

#SRCS = libt2.f90 t2cnst.f90 t2comm.f90 t2intg.f90 t2ngra.f90 t2vgra.f90\
#       t2init.f90 t2parm.f90 t2div.f90 t2prof.f90 t2writ.f90 \
#       t2prep.f90 t2calv.f90 t2exec.f90 t2save.f90 \
#       t2step.f90 t2loop.f90 t2gout.f90 t2menu.f90

OBJS = $(SRCS:.f90=.o) 

LIBS= ../lib/tasklib.a ../../bpsd/bpsdlib.a $(LIBMTX) \
      t2xlib.a ../plx/pllib.a ../eq/eqlib.a

.f90.o:
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: t2

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
../eq/eqlib.a:
	(cd ../eq;  make eqlib.a)
t2xlib.a: $(OBJS)
	$(LD) $(LDFLAGS) t2xlib.a $(OBJS)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../mtxp; make libmtxmpi.o libmtxnompi.o)
	(cd ../mtxp; make libmtxbnd.o libmtxksp.o libmtxmumps.o)

t2:$(LIBS) t2main.o
	$(FLINKER)  t2main.o -o $@ $(FFLAGS) $(LIBS) $(PETSCLIBS) $(FLIBS) 
#cleaning
clean:
	rm -f *.o *.mod *.a ./mod/*.mod t2 *.vtk *.dat *~
cleand:
	rm -f *.vtk *.dat 
t2main.o: t2main.f90 t2comm.f90 t2cnst.f90
t2menu.o: t2menu.f90 t2comm.f90 t2cnst.f90
t2init.o: t2init.f90 t2comm.f90 t2cnst.f90
t2parm.o: t2parm.f90 t2comm.f90 t2cnst.f90
t2prep.o: t2prep.f90 t2comm.f90 t2cnst.f90
t2div.o:  t2div.f90  t2comm.f90 t2cnst.f90
t2intg.o: t2intg.f90 t2comm.f90 t2cnst.f90
t2ngra.o: t2ngra.f90 t2comm.f90 t2cnst.f90
t2vgra.o: t2vgra.f90 t2comm.f90 t2cnst.f90
t2wrap.o: t2wrap.f90 t2comm.f90 t2cnst.f90
t2mfcs.o: t2mfcs.f90 t2comm.f90 t2cnst.f90
t2calv.o: t2calv.f90 t2comm.f90 t2cnst.f90
t2exec.o: t2exec.f90 t2comm.f90 t2cnst.f90
t2step.o: t2step.f90 t2comm.f90 t2cnst.f90
t2writ.o: t2writ.f90 t2comm.f90 t2cnst.f90
t2loop.o: t2loop.f90 t2comm.f90 t2cnst.f90
t2term.o: t2term.f90 t2comm.f90 t2cnst.f90
libt2.o: libt2.f90
