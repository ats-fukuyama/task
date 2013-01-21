### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)

MODINCLUDE= -I$(MOD) -I../lib/$(MOD) -I../plx/$(MOD) -I../eq/$(MOD) \
            -I../mtxp -I../../bpsd/$(MOD)

SRC0 = wmmenu.f wminit.f wmeign.f wmsetg.f \
       wmprof.f wmdisp.f \
       wmsetm.f wmsetf.f wmxprf.f wmdprf.f \
       wmpout.f wmgout.f wmgsub.f wmdout.f wmfout.f \
       wmexec.f wmloop.f wmtae.f \
       wmvmec.f wmbooz.f wmvout.f wmtpcs.f

SRCX = wmcont.f90 wmtest.f90

# wmtrfile.f

SRCS = $(SRC0) nowmsolv.f 
SRCM = $(SRC0) wmsolv.f 

OBJS = $(SRCX:.f90=.o) $(SRCS:.f=.o)
OBJM = $(SRCX:.f90=.o) $(SRCM:.f=.o)

LIB0 = ../eq/wmeqin.o ../eq/eqlib.a \
       ../plx/pllib.a ../dp/dplib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a 

LIBS = $(LIB0) wmlib.a ../mtx/mtxlib.a ../mpi/nompilib.o
LIBM = $(LIB0) mpiwmlib.a ../mtx/mpimtxlib.a ../mpi/mpilib.o ../mtxp/libmtxmpi.o

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: libs wmlib.a wm

wmlib.a: $(OBJS)
	$(LD) $(LDFLAGS) wmlib.a $(OBJS)
mpiwmlib.a: $(OBJM)
	$(LD) $(LDFLAGS) mpiwmlib.a $(OBJM)

../plx/pllib.a:
	(cd ../plx; make pllib.a)
../dp/dplib.a:
	(cd ../dp; make dplib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../eq/wmeqin.o:
	(cd ../eq; make wmeqin.o)
../mpi/nompilib.o:
	(cd ../mpi; make nompilib.o)
../mpi/mpilib.o:
	(cd ../mpi; make mpilib.o)
../mtx/mtxlib.a:
	(cd ../mtx; make mtxlib.a)
../mtx/mpimtxlib.a:
	(cd ../mtx; make mpimtxlib.a)
../mtxp/libmtxmpi.o:
	(cd ../mtxp; make libmtxmpi.o)
../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
libs:
	(cd ../lib; make tasklib.a)
	(cd ../mtx; make mtxlib.a)
	(cd ../mpi; make nompilib.o)
	(cd ../plx; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../eq; make wmeqin.o)
	(cd ../dp; make dplib.a)
#	(cd ../mpi; make mpilib.o)
#	(cd ../mtxp; make libmtxmpi.o)

wm : wmmain.o $(LIBS)
	$(FCFIXED) wmmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

mpiwm : wmmain.o $(LIBM)
	$(MFC) wmmain.o $(LIBM) -o $@ $(FFLAGS) $(FLIBS)

check : 
	ftnchek -array=0 -usage=no-com-block-unused $(SRCS) | less

new :
	mkdir ../wmnew
	cp -f Makefile ../wmnew
	cp -f *.f ../wmnew
	cp -f *.inc ../wmnew

clean:
	rm -f core *.o ./*~ ../mtx/*.o *.a $(MOD)/*.mod

veryclean: clean
	rm -f wm

WMCOMM = wmcomm.inc wmcom0.inc wmcom1.inc ../mpi/mpicom.inc ../plx/plcomm.f90 \
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
../eq/wmeqin.o: ../eq/wmeqin.f ../eq/eqcom1.inc ../eq/eqcom3.inc

wmtrfile.o : wmtrfile.f ../tr/trcomm.inc ../tr/trcom0.inc
wmtpcs.o : wmtpcs.f $(WMCOMM)

wmloop.o : wmloop.f $(WMCOMM)
wmcont.o : wmcont.f90
wmtest.o : wmtest.f90