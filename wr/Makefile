### $Id$ ###
include ../make.header
include ../mtxp/make.mtxp
LIB_MTX=$(LIB_MTX_BND)
LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_PCG)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIB_MTX=$(LIBXMTX_KSP)
#LIBX_MTX=$(LIBX_MTX_KSP)

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE= -I./$(MOD) -I../plx/$(MOD) -I../lib/$(MOD) -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

SRCS = wrmenu.f wrinit.f wrcalc.f wrbeam.f wrgout.f wrfile.f

OBJS = $(SRCS:.f=.o)

LIBS = wrlib.a  ../dp/dplib.a ../eq/eqlib.a ../plx/pllib.a \
       $(LIB_MTX) ../lib/tasklib.a ../../bpsd/bpsdlib.a
.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : wr

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../dp/dplib.a:
	(cd ../dp; make dplib.a)
wrlib.a: $(OBJS)
	$(LD) $(LDFLAGS) wrlib.a $(OBJS)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../plx; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../dp; make dplib.a)

wr : $(LIBS) wrmain.o
	$(FLINKER) wrmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f wr

new:
	-mkdir ../wrnew
	cp Makefile ../wrnew
	cp *.f ../wrnew
	cp *.inc ../wrnew

WRCOMM=wrcomm.inc wrcom1.inc ../dp/dpcom1.inc ../plx/plcomm.f90
wrmain.o: wrmain.f $(WRCOMM)
wrmenu.o: wrmenu.f $(WRCOMM)
wrinit.o: wrinit.f $(WRCOMM)
wrcalc.o: wrcalc.f $(WRCOMM) ../plx/plprof.f90
wrbeam.o: wrbeam.f $(WRCOMM) ../plx/plprof.f90
wrgout.o: wrgout.f $(WRCOMM)
