### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = plmenu.f plinit.f plprof.f pldata.f plintf.f pl_bpsd.f
SRCO = noeqlib.f noequlib.f

OBJS = $(SRCS:.f=.o) 
OBJO = $(SRCO:.f=.o) 

LIBS = pllib.a ../eq/eqlib.a ../lib/tasklib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@

all : pllib.a  noeqlib.o

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
libs:
	(cd ../lib; make tasklib.a)

pllib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o pllib.a $(OBJS)

pl : plmain.o $(LIBS)
	$(FCFIXED) plmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core a.out *.o ./*~ ./#* *.a

veryclean: clean
	-rm -f pl

new:
	-mkdir ../plnew
	cp -f Makefile ../plnew
	cp -f *.f ../plnew
	cp -f *.inc ../plnew

plmain.o: plmain.f plcomm.inc plcom0.inc plcom1.inc
plmenu.o: plmenu.f plcomm.inc plcom0.inc plcom1.inc
plinit.o: plinit.f plcomm.inc plcom0.inc plcom1.inc plcnst.inc
plintf.o: plintf.f plcomm.inc plcom0.inc plcom1.inc
plprof.o: plprof.f plcomm.inc plcom0.inc plcom1.inc plcom2.inc
pldata.o: pldata.f plcomd.inc plcom0.inc plcom3.inc
