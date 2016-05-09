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

LIBS = wimlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a

all : wim

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
wimlib.a: $(OBJS)
	$(LD) $(LDFLAGS) wimlib.a $(OBJS)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)

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
