### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)
FFLAGS = $(OFLAGS)

SRCS=trmain.f trinit.f trfile.f trloop.f trrslt.f \
     trcalc.f trcoef.f trpnb.f trprf.f trpnf.f trpel.f \
     trgrar.f trgrat.f trgrap.f trgsub.f trhelp.f trmtrx.f 
OBJS = $(SRCS:.f=.o)
PROG = trf

all: $(PROG)

$(PROG) : $(OBJS) trcomm.h
	$(FC) -o $@ $(FFLAGS) $(OBJS) $(FLIBS)

check : $(SRCS)
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean : 
	rm -f ./#* ./*~ *.o a.out core

veryclean : clean
	rm -f $(PROG)

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
trgsub.o : trgsub.f trcomm.h
trhelp.o : trhelp.f trcomm.h
trmtrx.o : trmtrx.f
	$(FC) -c $(LFLAGS) trmtrx.f


