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
SRCSFREE = xxcomm.f90 xxinit.f90 xxparm.f90 \
           xxprep.f90 xxexec.f90 xxgout.f90 \
           xxmenu.f90

OBJS = $(SRCSFIXED:.f=.o) $(SRCSFREE:.f90=.o) 

LIBS = libxx.a ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a

all : xx

../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
libxx.a: $(OBJS)
	$(LD) $(LDFLAGS) libxx.a $(OBJS)

libs:
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../../bpsd; make libbpsd.a)

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
xxprep.o: xxprep.f90 xxcomm.f90
xxexec.o: xxexec.f90 xxcomm.f90
xxgout.o: xxgout.f90 xxcomm.f90
