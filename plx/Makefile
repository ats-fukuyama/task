### $Id$ ###
include ../../task/make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

MODINCLUDE=-I$(MOD) -I../../bpsd/$(MOD)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

SRCSFIXED = 
#SRCSFREE = modules.f90 pl_vmec.f90 plcomm.f90 plxprf.f90 plinit.f90 \
#           plprof.f90 plmenu.f90
SRCSFREE = modules.f90 pl_vmec.f90 plcomm.f90 plinit.f90 \
           plprof.f90 plmenu.f90
SRCO = noeqlib.f noequlib.f

OBJS = $(SRCSFIXED:.f=.o) $(SRCSFREE:.f90=.o) 
OBJO = $(SRCO:.f=.o) 

LIBS = pllib.a ../../task/eq/eqlib.a ../../task/lib/tasklib.a ../../task/wm/wmlib.a

all : pllib.a  noeqlib.o noequlib.o

../../task/lib/tasklib.a:
	(cd ../../task/lib; make tasklib.a)

libs:
	(cd ../../task/lib; make tasklib.a)

pllib.a: $(OBJS)
	$(LD) $(LDFLAGS) pllib.a $(OBJS)

pl : plmain.o $(LIBS)
	$(FCFIXED) plmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MODDIR) $(MODINCLUDE)

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

plmain.o: plmain.f90 plcomm.f90
plmenu.o: plmenu.f90 plcomm.f90
plinit.o: plinit.f90 plcomm.f90
plintf.o: plintf.f90 plcomm.f90
pldata.o: pldata.f90 plcomm.f90
plprof.o: plprof.f90 plcomm.f90 plxprf.f90
plxprf.o: plxprf.f90
plcomm.o: plcomm.f90
