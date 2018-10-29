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

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../plx/$(MOD) -I../dpx/$(MOD) -I../lib/$(MOD) \
            -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

SRCS = wrcomm.f90 wrinit.f90 wrparm.f90 \
       wrcalc.f90 wrbeam.f90 wrgout.f90 wrfile.f90 \
       wrmenu.f90

OBJS = $(SRCS:.f90=.o)

LIBS = wrlib.a  ../dpx/dplib.a ../eq/eqlib.a ../plx/pllib.a \
       ../lib/tasklib.a ../../bpsd/bpsdlib.a

all : wr

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../dpx/dplib.a:
	(cd ../dpx; make dplib.a)
wrlib.a: $(OBJS)
	$(LD) $(LDFLAGS) wrlib.a $(OBJS)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../mtxp; make)
	(cd ../plx; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../dpx; make dplib.a)

wr : $(LIBS) wrmain.o
	$(FLINKER) wrmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f wr

WRCOMM=wrcomm.f90 ../dpx/dpcomm.f90 ../plx/plcomm.f90
wrmain.o: wrmain.f90 $(WRCOMM)
wrmenu.o: wrmenu.f90 $(WRCOMM)
wrinit.o: wrinit.f90 $(WRCOMM)
wrparm.o: wrparm.f90 $(WRCOMM)
wrcalc.o: wrcalc.f90 $(WRCOMM) ../plx/plprof.f90
wrbeam.o: wrbeam.f90 $(WRCOMM) ../plx/plprof.f90
wrgout.o: wrgout.f90 $(WRCOMM)
wrfile.o: wrfile.f90 $(WRCOMM)
