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
SRCSFREE = xxcomm.f90 xxexec.f90 xxgout.f90 \
           xxinit.f90 xxparm.f90 xxmenu.f90

#OBJS = $(SRCSFIXED:.f=.o) $(SRCSFREE:.f90=.o) 
OBJS = $(SRCSFREE:.f90=.o) 

LIBS = xxlib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a

all : xx

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
xxlib.a: $(OBJS)
	$(LD) $(LDFLAGS) xxlib.a $(OBJS)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)

xx : $(LIBS) xxmain.o
	$(FCFIXED) xxmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) \
           $(MODDIR) $(MODINCLUDE)

clean:
	-rm -f core a.out *.o *.mod ./*~ ./#* *.a $(MOD)/*.mod

veryclean: clean
	-rm -f xx

xxmain.o: xxmain.f90 xxcomm.f90
xxmenu.o: xxmenu.f90 xxcomm.f90
xxinit.o: xxinit.f90 xxcomm.f90
xxparm.o: xxparm.f90 xxcomm.f90
xxexec.o: xxexec.f90 xxcomm.f90
xxgout.o: xxgout.f90 xxcomm.f90
