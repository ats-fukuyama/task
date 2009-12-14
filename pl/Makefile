### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE=-I$(MOD) -I../../bpsd/$(MOD)

SRCSFIXED = plmenu.f plinit.f plprof.f pldata.f plintf.f
SRCSFREE = modules.f90 pl_vmec.f90 plcomm.f90 plunit.f90
SRCO = noeqlib.f noequlib.f

OBJS = $(SRCSFIXED:.f=.o) $(SRCSFREE:.f90=.o) 
OBJO = $(SRCO:.f=.o) 

LIBS = pllib.a ../eq/eqlib.a ../lib/tasklib.a ../wm/wmlib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : pllib.a  noeqlib.o noequlib.o

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)

libs:
	(cd ../lib; make tasklib.a)

pllib.a: $(OBJS)
	$(LD) $(LDFLAGS) pllib.a $(OBJS)
#	$(LD) $(LDFLAGS) -o pllib.a $(OBJS)

pl : plmain.o $(LIBS)
	$(FCFIXED) plmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCSFIXED) | less

clean:
	-rm -f core a.out *.o *.mod ./*~ ./#* *.a $(MOD)/*.mod

veryclean: clean
	-rm -f pl

new:
	-mkdir ../plnew
	cp -f Makefile ../plnew
	cp -f *.f ../plnew
	cp -f *.inc ../plnew

WMCOMM = ../wm/wmcomm.inc ../wm/wmcom0.inc ../wm/wmcom1.inc ../mpi/mpicom.inc ../dp/dpcom1.inc

plmain.o: plmain.f plcomm.f90 plcom0.f90 plcom1.inc
plmenu.o: plmenu.f plcomm.inc plcom0.inc plcom1.inc
plinit.o: plinit.f plcomm.inc plcom0.inc plcom1.inc plcnst.inc
plintf.o: plintf.f plcomm.inc plcom0.inc plcom1.inc
plprof.o: plprof.f plcomm.inc plcom0.inc plcom1.inc plcom2.inc plxprf.inc $(WMCOMM)
pldata.o: pldata.f plcomd.inc plcom0.inc plcom3.inc
