### $Id$ ###
include ../make.header

FFLAGS = $(OFLAGS)

SRCS = txmenu.f txinit.f txexec.f txcalv.f txcala.f \
       txfile.f txrslt.f txgraf.f txlib.f
OBJS = $(SRCS:.f=.o)
PROG = tx
LIBS=txlib.a ../lib/tasklib.a

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all: $(PROG)

txlib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o txlib.a $(OBJS)

$(PROG): txmain.o $(LIBS)
	(cd ../lib; make tasklib.a)
#	$(FC) txmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)
	$(FC) txmain.o $(LIBS) -o $@ $(FFLAGS) $(GLLIBS)

check: $(SRCS)
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	rm -f *.o ./*~ ./#* a.out core *.a

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
txmenu.o: txmenu.f txcomm.inc
txinit.o: txinit.f txcomm.inc
txexec.o: txexec.f txcomm.inc
txcalv.o: txcalv.f txcomm.inc
txcala.o: txcala.f txcomm.inc
txfile.o: txfile.f txcomm.inc
txrslt.o: txrslt.f txcomm.inc
txgraf.o: txgraf.f txcomm.inc
txlib.o: txlib.f
