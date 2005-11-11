### $Id$ ###
include ../make.header

FFLAGS = $(OFLAGS)

SRC0 = wmmenu.f wminit.f wmeign.f wmsetg.f wmvmec.f wmbooz.f \
       wmexec.f wmprof.f wmdisp.f \
       wmsetm.f wmsetf.f wmxprf.f wmdprf.f \
       wmpout.f wmgout.f wmvout.f wmgsub.f wmdout.f wmfout.f \
       wmtae.f wmtrfile.f

SRCS = $(SRC0) nowmsolv.f 
SRCM = $(SRC0) wmsolv.f 

OBJS = $(SRCS:.f=.o)
OBJM = $(SRCM:.f=.o)

LIB0 = ../eq/wmeqin.o ../eq/eqlib.a \
       ../pl/pllib.a ../dp/dplib.a ../lib/tasklib.a

LIBS = wmlib.a $(LIB0) ../mtx/mtxlib.a ../mpi/nompilib.o
LIBM = mpiwmlib.a $(LIB0) ../mtx/mpimtxlib.a ../mpi/mpilib.o

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all : wm

wmlib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o wmlib.a $(OBJS)
mpiwmlib.a: $(OBJM)
	$(LD) $(LDFLAGS) -o mpiwmlib.a $(OBJM)

../pl/pllib.a:
	(cd ../pl; make pllib.a)
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
../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
libs:
	(cd ../pl; make pllib.a)
	(cd ../dp; make dplib.a)
	(cd ../eq; make eqlib.a)
	(cd ../eq; make wmeqin.o)
	(cd ../pl; make pllib.a)
	(cd ../dp; make dplib.a)
	(cd ../mpi; make nompilib.o)
	(cd ../mpi; make mpilib.o)
	(cd ../mtx; make mtxlib.a)
	(cd ../mtx; make mpimtxlib.a)
	(cd ../lib; make tasklib.a)

wm : wmmain.o $(LIBS)
	$(FC) wmmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

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
	rm -f core *.o ./*~ ../mtx/*.o *.a

veryclean: clean
	rm -f wm

WMCOMM = wmcomm.inc wmcom0.inc wmcom1.inc ../mpi/mpicom.inc ../pl/plcom0.inc ../pl/plcom1.inc ../dp/dpcom1.inc

wmmain.o : wmmain.f $(WMCOMM)
wmmenu.o : wmmenu.f $(WMCOMM)
wmeign.o : wmeign.f $(WMCOMM)
wminit.o : wminit.f $(WMCOMM) ../pl/plcnst.inc
wmsetg.o : wmsetg.f $(WMCOMM)
wmvmec.o : wmvmec.f $(WMCOMM) vmcomm.inc
wmexec.o : wmexec.f $(WMCOMM)
wmprof.o : wmprof.f $(WMCOMM) wmxprf.inc
wmsolv.o : wmsolv.f $(WMCOMM)
nowmsolv.o : nowmsolv.f $(WMCOMM)
wmsetm.o : wmsetm.f $(WMCOMM)
wmsetf.o : wmsetf.f $(WMCOMM)
wmdisp.o : wmdisp.f $(WMCOMM)
wmtae.o : wmtae.f $(WMCOMM)
wmdprf.o : wmdprf.f $(WMCOMM)
wmxprf.o : wmxprf.f $(WMCOMM) wmxprf.inc
wmpout.o : wmpout.f $(WMCOMM)
wmgout.o : wmgout.f $(WMCOMM)
wmvout.o : wmvout.f $(WMCOMM) vmcomm.inc
wmgsub.o : wmgsub.f $(WMCOMM)
wmdout.o : wmdout.f $(WMCOMM)	
wmbooz.o : wmbooz.f $(WMCOMM) vmcomm.inc
../eq/wmeqin.o: ../eq/wmeqin.f ../eq/eqcom1.inc ../eq/eqcom3.inc
wmtrfile.o : wmtrfile.f ../tr/trcomm.inc ../tr/trcom0.inc
