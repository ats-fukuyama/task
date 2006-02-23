### $Id$ ###
include ../make.header

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

SRCS = txcnst.f90 txlib.f90 txrslt.f90 txcalv.f90 txcala.f90 txexec.f90 \
       txinit.f90 txfile.f90 txmenu.f90 txmain.f90 txgraf.f90 
OBJS = $(SRCS:.f90=.o)
PROG = tx
LIBS = txlib.a ../lib/tasklib.a $(OBJLA)

.SUFFIXES: .f90

.f90.o :
	$(F95) $(FFLAGS) -c $< -o $@

all: $(PROG)

txlib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o txlib.a $(OBJS)

../lib/tasklib.a:
	(cd ../lib; make FC="$(FC)" FFLAGS="$(FFLAGS)" LFLAGS="$(FFLAGS)")

libs:
	(cd ../lib; make)

$(PROG): $(LIBS)
#	$(F95) $(LIBS) -o $@ $(FFLAGS) $(GFLIBS)
	$(F95) $(LIBS) -o $@ $(FFLAGS) $(GLLIBS)

check: $(SRCS)
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	rm -f *.o ./*~ ./#* a.out core *.a *.mod

veryclean: clean
	rm -f $(PROG)

distclean: clean
	cd ../lib && make clean

new :
	mkdir ../txnew
	cp -f Makefile ../txnew
	cp -f HISTORY ../txnew
	cp -f *.f90 ../txnew
	cp -f *.inc ../txnew
	-mkdir ../txnew/doc
	-cp -f doc/* ../txnew/doc
	-mkdir ../txnew/in
	-cp -f in/* ../txnew/in

txmain.o: txmain.f90 txcomm.inc
txmenu.o: txmenu.f90 txcomm.inc
txinit.o: txinit.f90 txcomm.inc
txexec.o: txexec.f90 txcomm.inc
txcalv.o: txcalv.f90 txcomm.inc
txcala.o: txcala.f90 txcomm.inc
txfile.o: txfile.f90 txcomm.inc
txrslt.o: txrslt.f90 txcomm.inc
txgraf.o: txgraf.f90 txcomm.inc
txlib.o: txlib.f90
