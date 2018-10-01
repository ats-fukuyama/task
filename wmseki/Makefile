### Makefile for WM ###
include ../make.header
include ../mtxp/make.mtxp

# using mtxp complex solver
LIB_MTX=$(LIB_MTX_BND)
LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)

# select compile options
#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

MODINCLUDE= -I$(MOD) -I../lib/$(MOD) -I../plx/$(MOD) -I../eq/$(MOD) \
            -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

SRCS = wmmenu.f wminit.f wmeign.f wmsetg.f \
       ../eq/wmeqin.f wmprof.f wmdisp.f \
       wmsetm.f wmsetf.f wmxprf.f wmdprf.f \
       wmpout.f wmgout.f wmgsub.f wmdout.f wmfout.f \
       wmsolv.f wmexec.f wmloop.f wmtae.f \
       wmvmec.f wmbooz.f wmvout.f wmtpcs.f

SRCX = wmcont.f90 wmtest.f90

# wmtrfile.f

OBJS = $(SRCX:.f90=.o) $(SRCS:.f=.o)

LIBS =  wmlib.a ../plx/pllib.a ../eq/eqlib.a ../dpseki/dplib.a \
        ../lib/tasklib.a ../../bpsd/bpsdlib.a $(LIB_MTX)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: wm

wmlib.a: $(OBJS)
	$(LD) $(LDFLAGS) wmlib.a $(OBJS)

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
../dpseki/dplib.a:
	(cd ../dpseki; make dplib.a)
libs:
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../lib; make tasklib.a)
	(cd ../mtxp; make)
	(cd ../plx; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../dpseki; make dplib.a)

wm : $(LIBS) wmmain.o
	$(FLINKER) wmmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	rm -f core *.o ./*~ ../mtx/*.o *.a $(MOD)/*.mod

veryclean: clean
	rm -f wm

WMCOMM = wmcomm.inc wmcom0.inc wmcom1.inc ../plx/plcomm.f90 \
 ../dpseki/dpcom1.inc

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
