### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = fpmain.f fpinit.f fploop.f fpexec.f fpcoef.f fpcalv.f fpcaln.f \
       fpsave.f fpgraf.f fpfile.f fpcalw.f fpcals.f fpwrin.f fpfout.f

OBJS = $(SRCS:.f=.o) ../pl/plinit.o ../pl/plprof.o \
                     ../eq/eqload.o ../eq/pleqin.o ../eq/eqcalq.o \
                     ../lib/libspl.o ../lib/libspl2.o ../lib/libgut.o \
                     ../lib/libbrn.o \
                     ../lib/libde.o ../lib/libfp.o \
                     ../lib/libbes.o ../lib/libpcg.o 

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

fp : $(OBJS)
	$(FC) $(OBJS) -o fp $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean :
	-rm -f core a.out *.o ./#* ./*~

veryclean: clean
	-rm -f fp

new:
	-mkdir ../fpnew
	cp -f Makefile ../fpnew
	cp -f *.f ../fpnew
	cp -f *.h ../fpnew

fpmain.o : fpmain.f fpcomm.h fpcom1.h ../pl/plcom1.h
fpinit.o : fpinit.f fpcomm.h fpcom1.h ../pl/plcom1.h
fploop.o : fploop.f fpcomm.h fpcom1.h ../pl/plcom1.h
fpexec.o : fpexec.f fpcomm.h fpcom1.h ../pl/plcom1.h
fpcoef.o : fpcoef.f fpcomm.h fpcom1.h ../pl/plcom1.h
fpcalv.o : fpcalv.f fpcomm.h fpcom1.h ../pl/plcom1.h
fpcaln.o : fpcaln.f fpcomm.h fpcom1.h ../pl/plcom1.h

fpsave.o : fpsave.f fpcomm.h fpcom1.h ../pl/plcom1.h
fpgraf.o : fpgraf.f fpcomm.h fpcom1.h ../pl/plcom1.h
fpfout.o : fpfout.f fpcomm.h fpcom1.h ../pl/plcom1.h
fpfile.o : fpfile.f fpcomm.h fpcom1.h ../pl/plcom1.h
fpcalw.o : fpcalw.f fpcomm.h fpcom1.h ../pl/plcom1.h
fpcals.o : fpcals.f fpcomm.h fpcom1.h ../pl/plcom1.h ../wr/wrcom1.h fpcom2.h
fpwrin.o : fpwrin.f fpcomm.h fpcom1.h ../pl/plcom1.h ../wr/wrcom1.h fpcom2.h
../pl/plinit.o: ../pl/plinit.f ../pl/plcomm.h ../pl/plcom1.h
../pl/plprof.o: ../pl/plprof.f ../pl/plcomm.h ../pl/plcom1.h ../pl/plcom2.h
../eq/eqload.o: ../eq/eqload.f ../eq/eqcom1.h ../eq/eqcom3.h
../eq/pleqin.o: ../eq/pleqin.f ../eq/eqcom1.h ../eq/eqcom3.h
../eq/eqcalq.o: ../eq/eqcalq.f ../eq/eqcom1.h ../eq/eqcom3.h
