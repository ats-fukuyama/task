### $Id$ ###
include ../make.header

FFLAGS = $(OFLAGS)

SRC0 = wmmenu.f wminit.f wmeign.f wmsetg.f wmvmec.f wmbooz.f \
       wmexec.f wmprof.f wmdisp.f \
       wmsetm.f wmsetf.f wmxprf.f wmdprf.f \
       wmpout.f wmgout.f wmvout.f wmgsub.f wmdout.f wmtae.f

SRCS = $(SRC0) nowmsolv.f 
SRCM = $(SRC0) wmsolv.f 

OBJS = $(SRCS:.f=.o)
OBJM = $(SRCM:.f=.o)

LIB0 = ../eq/wmeqin.o ../eq/eqlib.a ../eq/pleqin.o \
       ../pl/pllib.a ../lib/tasklib.a

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
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../eq/wmeqin.o:
	(cd ../eq; make wmeqin.o)
../eq/pleqin.o:
	(cd ../eq; make pleqin.o)
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
	(cd ../eq; make eqlib.a)
	(cd ../eq; make wmeqin.o)
	(cd ../eq; make pleqin.o)
	(cd ../mpi; make nompilib.o)
	(cd ../mpi; make mpilib.o)
	(cd ../mtx; make mtxlib.a)
	(cd ../mtx; make mpimtxlib.a)
	(cd ../lib; make tasklib.a)

wm : wmmain.o $(LIBS)
	$(FC) wmmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(OBJLA)

mpiwm : wmmain.o $(LIBM)
	$(MFC) wmmain.o $(LIBM) -o $@ $(FFLAGS) $(FLIBS) $(OBJLA)

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

wmmain.o : wmmain.f wmcomm.inc
wmmenu.o : wmmenu.f wmcomm.inc
wmeign.o : wmeign.f wmcomm.inc
wminit.o : wminit.f wmcomm.inc
wmsetg.o : wmsetg.f wmcomm.inc
wmvmec.o : wmvmec.f wmcomm.inc vmcomm.inc
wmexec.o : wmexec.f wmcomm.inc
wmprof.o : wmprof.f wmcomm.inc wmxprf.inc
wmsolv.o : wmsolv.f wmcomm.inc
nowmsolv.o : nowmsolv.f wmcomm.inc
wmsetm.o : wmsetm.f wmcomm.inc
wmsetf.o : wmsetf.f wmcomm.inc
wmdisp.o : wmdisp.f wmcomm.inc
wmtae.o : wmtae.f wmcomm.inc
wmdprf.o : wmdprf.f wmcomm.inc
wmxprf.o : wmxprf.f wmcomm.inc wmxprf.inc
wmpout.o : wmpout.f wmcomm.inc
wmgout.o : wmgout.f wmcomm.inc
wmvout.o : wmvout.f wmcomm.inc vmcomm.inc
wmgsub.o : wmgsub.f wmcomm.inc
wmdout.o : wmdout.f wmcomm.inc	
wmbooz.o : wmbooz.f wmcomm.inc vmcomm.inc
../eq/wmeqin.o: ../eq/wmeqin.f ../eq/eqcom1.inc ../eq/eqcom3.inc
