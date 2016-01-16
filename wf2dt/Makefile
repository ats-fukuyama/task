### $Id: Makefile,v 1.20 2012/03/05 06:29:02 maruyama Exp $
include ../make.header
include ../mtxp/make.mtxp
#LIBX_MTX=$(LIBX_MTX_BND)
#LIBX_MTX=$(LIBX_MTX_PCG)
LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_KSP)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(OFLAGS)
#LFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)
LFLAGS = $(DFLAGS)

MODINCLUDE= -I./$(MOD) -I../plx/$(MOD) -I../eq/$(MOD) \
            -I../../bpsd/$(MOD) -I../mtxp/$(MOD) \
            -I../lib/$(MOD)

SRCS = 	wfcomm.f90 wfmenu.f90 wfinit.f90 wfdiv.f90 wflib.f90 wfsub.f90 \
        wfindx.f90 wfsort.f90 wfgsub.f90 wfgutl.f90 wfant.f90 wfwave.f90 \
	wfsolv.f90 wfgout.f90 wfprof.f90

OBJS = $(SRCS:.f90=.o)

LIBS = wflib.a ../plx/pllib.a ../eq/eqlib.a   \
        ../lib/tasklib.a ../../bpsd/bpsdlib.a $(LIB_MTX)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: wf

wflib.a: $(OBJS)
	$(LD) $(LDFLAGS) wflib.a $(OBJS)

../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../lib/tasklib.a:
	(cd ../lib; make tasklib.a)

libs:
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../lib; make tasklib.a)
	(cd ../mtxp; make)
	(cd ../plx; make pllib.a)
	(cd ../eq; make eqlib.a)

wf : libs $(LIBS) wfmain.o
	$(FLINKER) wfmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

save :
	mkdir ../wfx.save
	cp -f Makefile ../wfx.save
	cp -f HISTORY ../wfx.save
	cp -f *.f90 ../wfx.save
	cp -f *.c ../wfx.save
	cp -f *.inc ../wfx.save
	mkdir ../wfx.save/in
	cp -f in/wf.in* ../wfx.save/in
	mkdir ../wfx.save/doc
	cp -f doc/*.txt doc/*.tex ../wfx.save/doc

clean:
	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod

veryclean: clean
	rm -f wf testin 

wfmain.o: wfmain.f90 wfcomm.f90
wfmenu.o: wfmenu.f90 wfcomm.f90
wfinit.o: wfinit.f90 wfcomm.f90
wfdiv.o: wfdiv.f90 wfcomm.f90
wfant.o: wfant.f90 wfcomm.f90
wfwave.o: wfwave.f90 wfcomm.f90
wfsolv.o: wfsolv.f90 wfcomm.f90
wfprof.o: wfprof.f90 wfcomm.f90
wfsub.o: wfsub.f90 wfcomm.f90
wfnas.o: wfnas.f90 wfcomm.f90
wfsort.o: wfsort.f90 wfcomm.f90
wffile.o: wffile.f90 wfcomm.f90
wfindx.o: wfindx.f90 wfcomm.f90
wfgout.o: wfgout.f90 wfcomm.f90
wfgsub.o: wfgsub.f90 wfcomm.f90
wfefwg.o: wfefwg.f90 wfcomm.f90
wfgdum.o: wfgdum.f90
	$(FCFREE) $(LFLAGS) -c $< 
wfgutl.o: wfgutl.f90
	$(FCFREE) $(LFLAGS) -c $< 
wflib.o: wflib.f90 wfcomm.f90
	$(FCFREE) $(LFLAGS) -c $< $(MODDIR) $(MODINCLUDE)
wfbes.o: wfbes.f90
	$(FCFREE) $(LFLAGS) -c $< 
