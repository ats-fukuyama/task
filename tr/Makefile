### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)
FFLAGS = $(OFLAGS)

SRCS=trmain.f trinit.f trfile.f trloop.f trrslt.f \
     trcalc.f trcoef.f trpnb.f trprf.f trpnf.f trpel.f \
     trgrar.f trgrat.f trgrap.f trgrae.f trgsub.f trhelp.f trgra2.f
SRCD=trmain.f trinit.f trfile.f trloop.f trrslt.f \
     trcalc.f trcoef.f trpnb.f trprf.f trpnf.f trpel.f \
     trgrar.f trgrat.f trgrap.f trgrae.f trgsub3.f trhelp.f trgrad.f
SRCE = ../eq/treqin.f ../eq/eqinit.f ../eq/eqcalc.f ../eq/eqfunc.f \
       ../eq/eqcalq.f ../eq/eqload.f ../eq/eqgout.f ../eq/eqgsub.f
SRCL = ../lib/libbnd.f ../lib/libspl.f ../lib/libbrn.f ../lib/libgut.f
OBJS = $(SRCS:.f=.o)
OBJD = $(SRCD:.f=.o)
OBJE = $(SRCE:.f=.o)
OBJL = $(SRCL:.f=.o)
PROG = trf2
PROG3 = trf

all: $(PROG3)

$(PROG3) : $(OBJD) $(OBJE) trcomm.h
	$(FC) -o $@ $(FFLAGS) $(OBJD) $(OBJE) $(OBJL) -L/usr/local/lib -lgsgl-ixc

2d : $(PROG)

$(PROG) : $(OBJS) $(OBJE) trcomm.h
	$(FC) -o $@ $(FFLAGS) $(OBJS) $(OBJE) $(OBJL) $(FLIBS)

check : $(SRCS)
	ftnchek $(FCKFLAGS) $(SRCS) $(SRCE) $(SRCL) $(FCKPIPE) | less

clean : 
	rm -f ./#* ./*~ *.o a.out core

veryclean : clean
	rm -f $(PROG) $(PROG3)

new :
	-mkdir ../trnew
	cp Makefile ../trnew
	-cp -f HISTORY ../trnew
	cp *.f ../trnew
	cp *.h ../trnew
	-mkdir ../trnew/in
	-cp -f in/* ../trnew/in
	-mkdir ../trnew/doc
	-cp -f doc/* ../trnew/doc

trmain.o : trmain.f trcomm.h
trinit.o : trinit.f trcomm.h
trfile.o : trfile.f trcomm.h
trloop.o : trloop.f trcomm.h
trrslt.o : trrslt.f trcomm.h

trcalc.o : trcalc.f trcomm.h
trcoef.o : trcoef.f trcomm.h
trpnb.o  : trpnb.f  trcomm.h
trprf.o  : trprf.f  trcomm.h
trpnf.o  : trpnf.f  trcomm.h
trpel.o  : trpel.f  trcomm.h

trgrar.o : trgrar.f trcomm.h
trgrat.o : trgrat.f trcomm.h
trgrap.o : trgrap.f trcomm.h
trgrae.o : trgrae.f trcomm.h
trgsub.o : trgsub.f trcomm.h
trgsub3.o : trgsub3.f trcomm.h
trhelp.o : trhelp.f trcomm.h
trgrad.o : trgrad.f trcomm.h
trgrp2.o : trgrp2.f trcomm.h
trmtrx.o : trmtrx.f
	$(FC) -c $(LFLAGS) trmtrx.f
../eq/eqinit.o: ../eq/eqinit.f ../eq/eqcom1.h
../eq/eqcalc.o: ../eq/eqcalc.f ../eq/eqcom1.h ../eq/eqcom2.h
../eq/eqcalq.o: ../eq/eqcalq.f ../eq/eqcom1.h ../eq/eqcom3.h
../eq/eqgout.o: ../eq/eqgout.f ../eq/eqcom1.h ../eq/eqcom2.h ../eq/eqcom3.h
../eq/eqfunc.o: ../eq/eqfunc.f ../eq/eqcom1.h ../eq/eqcom2.h ../eq/eqcom4.h
../eq/treqin.o: ../eq/treqin.f ../eq/eqcom1.h ../eq/eqcom3.h ../eq/eqcom4.h
../eq/eqload.o: ../eq/eqcalq.f ../eq/eqcom1.h ../eq/eqcom3.h


