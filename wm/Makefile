### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)
FFLAGS = $(OFLAGS)

SRCX = wminit.f wmsetg.f wmvmec.f wmsetj.f wmprof.f \
       nowmsolv.f wmsetm.f wmsetf.f wmdisp.f wmtae.f wmxprf.f \
       wmpout.f wmgout.f wmvout.f wmgsub.f wmdout.f nompilib.f

SRCS = wmmain.f wmeign.f $(SRCX)

OBJM = ../mtx/bandcd.o  ../mtx/mtxset.o ../mtx/precndcb.o \
       ../mtx/bcgcdb.o ../mtx/cgscdb.o ../mtx/bstabcdb.o \
       ../mtx/bandcdb.o

OBJL = ../lib/libdsp.o ../lib/libfft.o ../lib/libinv.o \
       ../lib/libspl.o ../lib/libbrn.o ../lib/libbnd.o 

OBJE = ../eq/eqload.o ../eq/eqcalq.o ../eq/wmeqin.o

OBJS = $(SRCS:.f=.o) $(OBJM) $(OBJL) $(OBJE) $(OBJLA)

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all : wm

wm : $(OBJS) wmcomm.h
	$(FC) $(OBJS) -o $@ $(FFLAGS) $(FLIBS)

check : 
	ftnchek -array=0 -usage=no-com-block-unused $(SRCS) | less

new :
	mkdir ../wmnew
	cp -f Makefile ../wmnew
	cp -f Makefile.mpi ../wmnew
	cp -f Makefile.nompi ../wmnew
	cp -f *.f ../wmnew
	cp -f *.h ../wmnew

clean:
	rm -f core *.o ./*~ ../mtx/*.o

veryclean: clean
	rm -f wm

wmmain.o : wmmain.f wmcomm.h
wmeign.o : wmeign.f wmcomm.h
wminit.o : wminit.f wmcomm.h
wmsetg.o : wmsetg.f wmcomm.h
wmvmec.o : wmvmec.f wmcomm.h vmcomm.h
wmsetj.o : wmsetj.f wmcomm.h
wmprof.o : wmprof.f wmcomm.h wmxprf.h
wmsolv.o : wmsolv.f wmcomm.h
nowmsolv.o : nowmsolv.f wmcomm.h
wmsetm.o : wmsetm.f wmcomm.h
wmsetf.o : wmsetf.f wmcomm.h
wmdisp.o : wmdisp.f wmcomm.h
wmtae.o : wmtae.f wmcomm.h
wmxprf.o : wmxprf.f wmcomm.h wmxprf.h
wmpout.o : wmpout.f wmcomm.h
wmgout.o : wmgout.f wmcomm.h
wmvout.o : wmvout.f wmcomm.h vmcomm.h
wmgsub.o : wmgsub.f wmcomm.h
wmdout.o : wmdout.f wmcomm.h	
../eq/eqload.o: ../eq/eqload.f ../eq/eqcom1.h ../eq/eqcom3.h
../eq/eqcalq.o: ../eq/eqcalq.f ../eq/eqcom1.h ../eq/eqcom3.h
../eq/wmeqin.o: ../eq/wmeqin.f ../eq/eqcom1.h ../eq/eqcom3.h
