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

#FFLAGS=-g $(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE=-I$(MOD) -I../../bpsd/$(MOD) -I../pl/$(MOD) -I../lib/$(MOD)

SRCX = equcom.f90 equread.f90
SRCS = eqbpsd.f eqinit.f eqmenu.f eqcalc.f eqcalq.f eqcalv.f \
       eqsub.f eqfunc.f eqintf.f eqsplf.f equintf.f \
       eq-eqdsk.f eq-qst.f eqfile.f \
       eqgout.f eqgsub.f \
       newton.f invematrix.f equnit.f eqrppl.f
SRCT = trtest.f treqin.f
SRCP = pltest.f
SRCA = aktest.f akeqin.f
SRC2D= eqg2d.f
SRC3D= eqg3d.f

OBJS = $(SRCX:.f90=.o) $(SRCS:.f=.o)
OBJT = $(SRCT:.f=.o)
OBJP = $(SRCP:.f=.o)
OBJA = $(SRCA:.f=.o)
OBJ2D=$(SRC2D:.f=.o)
OBJ3D=$(SRC3D:.f=.o)

LIBS2 = ../pl/libpl.a libeq.a ../lib/libgrf.a ../lib/libtask.a ../../bpsd/libbpsd.a
LIBS3 = ../pl/libpl.a libeq3.a ../lib/libgrf.a ../lib/libtask.a ../../bpsd/libbpsd.a

.SUFFIXES: .f90

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : eq


../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgr.a)
../pl/libpl.a:
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
libeq.a: $(OBJS) $(OBJ2D)
	(cd ../pl; make plcomm.o)
	$(LD) $(LDFLAGS) libeq.a $(OBJS) $(OBJ2D)
libeq3.a: $(OBJS) $(OBJ3D)
	$(LD) $(LDFLAGS) libeq3.a $(OBJS) $(OBJ3D)

libs:
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../../bpsd; make libbpsd.a)
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)

eq: libs libeq.a eqmain.o
	$(FLINKER) eqmain.o $(LIBS2) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

eq3: libs libeq3.a eqmain.o
	$(FLINKER) eqmain.o $(LIBS3) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

tr : libs libeq.a $(OBJT)
	$(FLINKER) $(OBJT) $(LIBS2) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

pl : libs libeq.a $(OGJP)
	$(FLINKER) $(OBJP) $(LIBS2) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

ak : libs libeq.a $(OBJA)
	$(FLINKER) $(OBJA) $(LIBS2) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out *.o *.mod ./*~ *.a $(MOD)/*.mod
	-rm -f doc/*~

veryclean: clean
	rm -f eq pl ak

EQCOMM = eqcom0.inc eqcom1.inc

eqmain.o: eqmain.f $(EQCOMM)
eqinit.o: eqinit.f $(EQCOMM)
eqcalc.o: eqcalc.f $(EQCOMM)
eqcalx.o: eqcalx.f $(EQCOMM) eqcom5.inc
eqcalq.o: eqcalq.f $(EQCOMM) eqcom3.inc
eqcalv.o: eqcalv.f $(EQCOMM) eqcom2.inc
eqsub.o : eqsub.f  $(EQCOMM) eqcom2.inc
eqfunc.o: eqfunc.f $(EQCOMM) eqcom2.inc eqcom4.inc
eqsplf.o: eqsplf.f $(EQCOMM) eqcom3.inc
eqfile.o: eqfile.f $(EQCOMM) eqcom2.inc
eq-eqdsk.o: eq-eqdsk.f $(EQCOMM) eqcom2.inc
eq-qst.o: eq-qst.f $(EQCOMM) eqcom2.inc
eqintf.o: eqintf.f $(EQCOMM) eqcom2.inc
eqgout.o: eqgout.f $(EQCOMM) eqcom2.inc eqcom3.inc
eqbpsd.o: eqbpsd.f $(EQCOMM) eqcom3.inc
eqjft2m.o: eqjft2m.f $(EQCOMM) eqcom3.inc
equnit.o: equnit.f eqbpsd.f
eqgsub.o: eqgsub.f
eqg2d.o: eqg2d.f $(EQCOMM) eqcom3.inc
eqg3d.o: eqg3d.f $(EQCOMM) eqcom3.inc
equcom.o: equcom.f90
equread.o: equread.f90 equcom.f90
equintf.o: equintf.f $(EQCOMM) eqcom3.inc
eqrppl.o: eqrppl.f $(EQCOMM) eqcom3.inc
trtest.o: trtest.f
treqin.o: treqin.f $(EQCOMM) eqcom3.inc eqcom4.inc
pltest.o: pltest.f
aktest.o: aktest.f
akeqin.o: akeqin.f $(EQCOMM) eqcom3.inc
newton.o: newton.f
invematrix.o: invematrix.f
