### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)
FFLAGS = $(OFLAGS)

SRCS = txmain.f txinit.f txexec.f txcalv.f txcala.f \
       txfile.f txrslt.f txgraf.f txlib.f txband.f
OBJS = $(SRCS:.f=.o)
PROG = tx

all: $(PROG)

$(PROG): $(OBJS) txcomm.inc
	$(FC) -o $@ $(FFLAGS) $(OBJS) $(FLIBS)

check: $(SRCS)
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	rm -f *.o ./*~ ./#* a.out core

veryclean: clean
	rm -f $(PROG)

new :
	mkdir ../txnew
	cp -f Makefile ../txnew
	cp -f HISTORY ../txnew
	cp -f *.f ../txnew
	cp -f *.inc ../txnew
	-mkdir ../txnew/doc
	-cp -f doc/* ../txnew/doc
	-mkdir ../txnew/in
	-cp -f in/* ../txnew/in

txmain.o: txmain.f txcomm.inc
txinit.o: txinit.f txcomm.inc
txexec.o: txexec.f txcomm.inc
txcalv.o: txcalv.f txcomm.inc
txcala.o: txcala.f txcomm.inc
txfile.o: txfile.f txcomm.inc
txrslt.o: txrslt.f txcomm.inc
txgraf.o: txgraf.f txcomm.inc

txband.o: ${@:.o=.f}
	$(FC) $(LFLAGS) -c ${@:.o=.f}
txlib.o: ${@:.o=.f}
	$(FC) $(LFLAGS) -c ${@:.o=.f}
