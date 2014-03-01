### $Id$###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

MODINCLUDE=-I $(MOD) -I../../bpsd/$(MOD) -I../lib/$(MOD)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

SRCSFIXED = 
SRCSFREE = wicomm.f90 wiparm.f90 \
	   wiunmag.f90 wifluid.f90 wiexec.f90 wiscan.f90 wigout.f90 \
           wimenu.f90

OBJS = $(SRCSFIXED:.f=.o) $(SRCSFREE:.f90=.o) 

LIBS = wilib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a

all : wi

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
wilib.a: $(OBJS)
	$(LD) $(LDFLAGS) wilib.a $(OBJS)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)

wi : $(LIBS) wimain.o
	$(FCFREE) wimain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) \
           $(MODDIR) $(MODINCLUDE)

wiunmag : wiunmag-old.f
	$(FCFIXED) wiunmag-old.f ../lib/tasklib.a -o $@ $(FFLAGS) $(FLIBS)

clean:
	-rm -f core a.out *.o *.mod ./*~ ./#* *.a $(MOD)/*.mod

veryclean: clean
	-rm -f wi wiunmag

wimain.o: wimain.f90 wicomm.f90
wimenu.o: wimenu.f90 wicomm.f90
wiparm.o: wiparm.f90 wicomm.f90
wiexec.o: wiexec.f90 wicomm.f90
wiunmag.o: wiunmag.f90 wicomm.f90
wifluid.o: wifluid.f90 wicomm.f90
wiscan.o: wiscan.f90 wicomm.f90
wigout.o: wigout.f90 wicomm.f90
