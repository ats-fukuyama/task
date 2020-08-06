### $Id$ ###
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
SRCSFREE = wimcomm.f90 wiminit.f90 wimparm.f90 wimexec.f90 wimgout.f90 \
           wimmenu.f90

OBJS = $(SRCSFREE:.f90=.o) 

LIBS = libwim.a ../lib/libgrf.a ../lib/libtask.a ../../bpsd/libbpsd.a

all : wim

../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
libwim.a: $(OBJS)
	$(LD) $(LDFLAGS) libwim.a $(OBJS)

libs:
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../../bpsd; make libbpsd.a)

wim : $(LIBS) wimmain.o
	$(FCFIXED) wimmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) \
           $(MODDIR) $(MODINCLUDE)

clean:
	-rm -f core a.out *.o *.mod ./*~ ./#* *.a $(MOD)/*.mod

veryclean: clean
	-rm -f xx

wimmain.o: wimmain.f90 wimcomm.f90
wimmenu.o: wimmenu.f90 wimcomm.f90
wiminit.o: wiminit.f90 wimcomm.f90
wimparm.o: wimparm.f90 wimcomm.f90
wimexec.o: wimexec.f90 wimcomm.f90
wimgout.o: wimgout.f90 wimcomm.f90
