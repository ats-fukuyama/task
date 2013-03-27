### $Id$
include ../make.header

##### Without PETSc library #####
#FLINKER=$(FCFREE)
#LIBMTX=../../task/mtxp/libmtxbnd.o ../../task/mtxp/libmtxnompi.o
##### With MUMPS library #####
FLINKER=$(MF90)
LIBMTX=../../task/mtxp/libmtxmumps.o ../../task/mtxp/libmtxmpi.o
LIBMUMPSX=-L${PETSC_DIR}/${PETSC_ARCH}/lib -ldmumps -lzmumps -lmumps_common \
  -lpord -lscalapack -lblacs -lparmetis -lmetis -llapack -lblas

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE= -I../../task/plx/$(MOD) -I../../task/eq/$(MOD) \
            -I../../bpsd/$(MOD) -I../../task/mtxp -I $(MOD)

SRCS = wfcomm.f90 wfmenu.f90 wfinit.f90 wfdiv.f90 wfant.f90 wfwave.f90 \
       wfsolv.f90 wfprof.f90 wfsub.f90 wfnas.f90 wfsort.f90 wffile.f90 \
       wflib.f90 wfindx.f90 wfefwg.f90 wfbes.f90 wfgout.f90 wfmag.f90 \
       wfgsub.f90 wfgutl.f90 mbant.f90
OBJS = $(SRCS:.f90=.o)

LIBS = wflib.a ../../task/plx/pllib.a ../../task/eq/eqlib.a   \
        ../../task/lib/tasklib.a ../../bpsd/bpsdlib.a $(LIBMTX)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: wf

wflib.a: $(OBJS)
	$(LD) $(LDFLAGS) wflib.a $(OBJS)

../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../../task/plx/pllib.a:
	(cd ../../plx; make pllib.a)
../../task/eq/eqlib.a:
	(cd ../../task/eq; make eqlib.a)
../../task/mtxp/libmtxbnd.o:
	(cd ../../task/mtxp; make libmtxbnd.o)
../../task/mtxp/libmtxmumps.o:
	(cd ../../task/mtxp; make libmtxmumps.o)
../../task/mtxp/libmtxmpi.o:
	(cd ../../task/mtxp; make libmtxmpi.o)
../../task/mtxp/libmtxnompi.o:
	(cd ../../task/mtxp; make libmtxnompi.o)
../../task/lib/tasklib.a:
	(cd ../../task/lib; make tasklib.a)

libs:
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../../task/eq; make eqlib.a)
	(cd ../../task/mtxp; make libmtxbnd.o libmtxmumps.o)
	(cd ../../task/plx; make pllib.a)
	(cd ../../task/lib; make tasklib.a)

wf : $(LIBS) wfmain.o
	$(FLINKER) wfmain.o -o $@ $(FFLAGS) $(LIBS) $(FLIBS) $(LIBMUMPSX)

save :
	mkdir ../wfx.save
	cp -f Makefile ../wfx.save
	cp -f HISTORY ../wfx.save
	cp -f *.f ../wfx.save
	cp -f *.c ../wfx.save
	cp -f *.inc ../wfx.save
	mkdir ../wfx.save/in
	cp -f in/wf.in* ../wfx.save/in
	mkdir ../wfx.save/doc
	cp -f doc/*.txt doc/*.tex ../wfx.save/doc

clean:
	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod

veryclean: clean
	rm -f wf testin nascheck

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
wfmag.o: wfmag.f90
	$(FCFREE) $(LFLAGS) -c $< 
wfgdum.o: wfgdum.f90
	$(FCFREE) $(LFLAGS) -c $< 
wfgutl.o: wfgutl.f90
	$(FCFREE) $(LFLAGS) -c $< 
wflib.o: wflib.f90 wfcomm.f90
	$(FCFREE) $(LFLAGS) -c $< $(MODDIR) $(MODINCLUDE)
wfbes.o: wfbes.f90
	$(FCFREE) $(LFLAGS) -c $< 
