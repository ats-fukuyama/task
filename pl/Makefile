### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

MODINCLUDE=-I$(MOD) -I../lib/$(MOD) -I../../bpsd/$(MOD)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

SRCSFIXED = 
SRCSFREE = modules.f90 pl_vmec.f90 plcomm.f90 plinit.f90 plparm.f90 \
           plload.f90 plprof2d.f90 plprof.f90 plprofw.f90 \
           plgout.f90 plmenu.f90
SRCO = noeqlib.f noequlib.f

OBJS = $(SRCSFIXED:.f=.o) $(SRCSFREE:.f90=.o) 
OBJO = $(SRCO:.f=.o) 

LIBS = libpl.a ../eq/libeq.a ../lib/libtask.a ../../bpsd/libbpsd.a

all : libpl.a  noeqlib.o noequlib.o pl

../lib/libtask.a:
	(cd ../lib; make libtask.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
libpl.a: $(OBJS)
	$(LD) $(LDFLAGS) libpl.a $(OBJS)

libs:
	(cd ../lib; make libtask.a)
	(cd ../../bpsd; make libbpsd.a)
	(cd ../eq; make libeq.a)

pl : $(LIBS) plmain.o
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

plmain.o: plmain.f90 plcomm.f90
plmenu.o: plmenu.f90 plcomm.f90
plinit.o: plinit.f90 plcomm.f90
plparm.o: plparm.f90 plcomm.f90
plintf.o: plintf.f90 plcomm.f90
pldata.o: pldata.f90 plcomm.f90
plprof.o: plprof.f90 plcomm.f90 plload.f90
plprofw.o: plprofw.f90 plcomm.f90 plprof.f90
plgout.o: plgout.f90 plcomm.f90
plload.o: plload.f90 plcomm.f90
plcomm.o: plcomm.f90
