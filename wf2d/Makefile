### $Id$ ###
include ../make.header
include ../mtxp/make.mtxp

#LIB_MTX=$(LIB_MTX_BND)
#LIBX_MTX=$(LIBX_MTX_BND)
LIB_MTX=$(LIB_MTX_MUMPS)
LIBX_MTX=$(LIBX_MTX_MUMPS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(OFLAGS)
#LFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)
LFLAGS = $(DFLAGS)

MODINCLUDE= -I./$(MOD) -I../pl/$(MOD) -I../eq/$(MOD) \
            -I../../bpsd/$(MOD) -I../mtxp/$(MOD) \
            -I../lib/$(MOD)

SRCS = 	wfcomm.f90 wfinit.f90 wfparm.f90 wfload.f90 wfprof.f90 \
        wfdiv.f90 wflib.f90 wfsub.f90 \
        wfindx.f90 wfsort.f90 wfgsub.f90 wfgutl.f90 wfant.f90 wfwave.f90 \
	wfsolv.f90 wfgout.f90 wfmenu.f90

OBJS = $(SRCS:.f90=.o)

LIBS = libwf.a ../pl/libpl.a ../eq/libeq.a   \
        ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: wf

libwf.a: $(OBJS)
	$(LD) $(LDFLAGS) libwf.a $(OBJS)

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make loibbpsd.a)
../pl/liobpl.a:
	(cd ../pl; make loibpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libgrf.a)
	(cd ../lib; make libtask.a)
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)

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

wfcomm.o: wfcomm.f90 ../pl/plcomm.f90
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
wfload.o: wfload.f90 wfcomm.f90
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
