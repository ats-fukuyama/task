### $Id$ ###
include ../make.header

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

##### Without PETSc library (no MPI, real8 and complex8) #####
FLINKER=$(FCFREE)
LIBNOMPI=../mtxp/libmtxnompi.o
LIBMTXP=$(LIBNOMPI) ../mtxp/libmtxbnd.o
LIBMUMPSX=.

##### With MUMPS library (MPI, real8 and complex8)
#FLINKER=$(MF90)
#LIBMPI=../mtxp/libmtxmpi.o
#LIBMTXP=$(LIBMPI) ../mtxp/libmtxmumps.o
#LIBMUMPSX=-L${PETSC_DIR}/${PETSC_ARCH}/lib -ldmumps -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lparmetis -lmetis -llapack -lblas

MODINCLUDE= -I$(MOD) -I../lib/$(MOD) -I../plx/$(MOD) -I../eq/$(MOD) \
            -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

SRC0 = wmmenu.f wminit.f wmeign.f wmsetg.f \
       ../eq/wmeqin.f wmprof.f wmdisp.f \
       wmsetm.f wmsetf.f wmxprf.f wmdprf.f \
       wmpout.f wmgout.f wmgsub.f wmdout.f wmfout.f \
       wmexec.f wmloop.f wmtae.f \
       wmvmec.f wmbooz.f wmvout.f wmtpcs.f

SRCX = wmcont.f90 wmtest.f90 wmdpfa.f90

# wmtrfile.f

SRCS = $(SRC0) nowmsolv.f 
SRCM = $(SRC0) wmsolv.f 

OBJS = $(SRCX:.f90=.o) $(SRCS:.f=.o)
OBJM = $(SRCX:.f90=.o) $(SRCM:.f=.o)

LIB0 = ../lib/tasklib.a ../../bpsd/bpsdlib.a ../plx/pllib.a \
       ../eq/eqlib.a ../dp/dplib.a

LIBS = $(LIBMTXP) ../mpi/mpilib.o ../mtx/mtxlib.a $(LIB0) wmlib.a  
LIBM = $(LIBMTXP) ../mpi/mpilib.o ../mtx/mpimtxlib.a $(LIB0) mpiwmlib.a 

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: wm

wmlib.a: $(OBJS)
	$(LD) $(LDFLAGS) wmlib.a $(OBJS)
mpiwmlib.a: $(OBJM)
	$(LD) $(LDFLAGS) mpiwmlib.a $(OBJM)

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
../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../eq/wmeqin.o:
	(cd ../eq; make wmeqin.o)
../dp/dplib.a:
	(cd ../dp; make dplib.a)
libs:
	(cd ../lib; make tasklib.a)
	(cd ../mtxp; make)
	(cd ../mpi; make mpilib.o)
	(cd ../mtx; make mtxlib.a)
	(cd ../mtx; make mpimtxlib.a)
	(cd ../plx; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../dp; make dplib.a)

wm : $(LIBS) wmmain.o
	$(FCFIXED) wmmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

mpiwm : wmmain.o $(LIBM)
	$(MFC) wmmain.o $(LIBM) -o $@ $(FFLAGS) $(FLIBS)

clean:
	rm -f core *.o ./*~ ../mtx/*.o *.a $(MOD)/*.mod

veryclean: clean
	rm -f wm

WMCOMM = wmcomm.inc wmcom0.inc wmcom1.inc ../plx/plcomm.f90 \
 ../dp/dpcom1.inc

wmmain.o : wmmain.f $(WMCOMM)

wmmenu.o : wmmenu.f $(WMCOMM)
wminit.o : wminit.f $(WMCOMM)
wmeign.o : wmeign.f $(WMCOMM)
wmsetg.o : wmsetg.f $(WMCOMM)

wmexec.o : wmexec.f $(WMCOMM)
wmprof.o : wmprof.f $(WMCOMM)
wmdisp.o : wmdisp.f $(WMCOMM)

wmsetm.o : wmsetm.f $(WMCOMM)
wmsetf.o : wmsetf.f $(WMCOMM)
wmxprf.o : wmxprf.f $(WMCOMM)
wmdprf.o : wmdprf.f $(WMCOMM)

wmpout.o : wmpout.f $(WMCOMM)
wmgout.o : wmgout.f $(WMCOMM)
wmgsub.o : wmgsub.f $(WMCOMM)
wmdout.o : wmdout.f $(WMCOMM)	
wmfout.o : wmfout.f $(WMCOMM)	

wmsetm.o : wmsetm.f $(WMCOMM)
wmsetm.o : wmsetm.f $(WMCOMM)

wmtae.o : wmtae.f $(WMCOMM)

wmvmec.o : wmvmec.f $(WMCOMM) vmcomm.inc
wmbooz.o : wmbooz.f $(WMCOMM) vmcomm.inc
wmvout.o : wmvout.f $(WMCOMM) vmcomm.inc

wmsolv.o : wmsolv.f $(WMCOMM)
nowmsolv.o : nowmsolv.f $(WMCOMM)
wmvout.o : wmvout.f $(WMCOMM) vmcomm.inc

../eq/wmeqin.o : ../eq/wmeqin.f
wmtrfile.o : wmtrfile.f ../tr/trcomm.inc ../tr/trcom0.inc
wmtpcs.o : wmtpcs.f $(WMCOMM)

wmloop.o : wmloop.f $(WMCOMM)
wmcont.o : wmcont.f90
wmtest.o : wmtest.f90
