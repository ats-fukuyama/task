### $Id$###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE=-I $(MOD) -I../../bpsd/$(MOD) -I../lib/$(MOD)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

SRCSFIXED = 
SRCSFREE = wicomm.f90 wiparm.f90 wiprep.f90 \
	   wicold.f90 wiwarm.f90 wihot.f90 wiexec.f90 \
           wigout.f90 wiscan.f90 wimenu.f90

OBJS = $(SRCSFIXED:.f=.o) $(SRCSFREE:.f90=.o) 

LIBS = libwi.a ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a

all : wi plotwi plotwi2

../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
libwi.a: $(OBJS)
	$(LD) $(LDFLAGS) libwi.a $(OBJS)

libs:
	(cd ../lib; make libgrf.a)
	(cd ../lib; make libtask.a)
	(cd ../../bpsd; make libbpsd.a)

wi : $(LIBS) wimain.o
	$(FCFREE) wimain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) \
           $(MODDIR) $(MODINCLUDE)

wiunmag : wiunmag-old.f
	$(FCFIXED) wiunmag-old.f ../lib/libtask.a -o $@ $(FFLAGS) $(FLIBS)

clean:
	-rm -f core a.out *.o *.mod ./*~ ./#* *.a $(MOD)/*.mod

plotwi : plotwi.o
	$(FCFREE) plotwi.o ../lib/libgrf.a ../lib/libtask.a -o $@ $(FFLAGS) $(FLIBS)

plotwi2 : plotwi2.o
	$(FCFREE) plotwi2.o ../lib/libgrf.a ../lib/libtask.a -o $@ $(FFLAGS) $(FLIBS)

veryclean: clean
	-rm -f wi wiunmag

wimain.o: wimain.f90 wicomm.f90
wimenu.o: wimenu.f90 wicomm.f90
wiparm.o: wiparm.f90 wicomm.f90
wiexec.o: wiexec.f90 wicomm.f90
wiprep.o: wiscan.f90 wicomm.f90
wicold.o: wicold.f90 wicomm.f90
wiwarm.o: wiwarm.f90 wicomm.f90
wihot.o: wihot.f90 wicomm.f90
wiscan.o: wiscan.f90 wicomm.f90
wigout.o: wigout.f90 wicomm.f90
plotwi.o: plotwi.f90
plotwi2.o: plotwi2.f90
