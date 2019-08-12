### $Id$
include ../make.header
include ../mtxp/make.mtxp
#LIB_MTX=$(LIB_MTX_BND)
#LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_PCG)
LIB_MTX=$(LIB_MTX_MUMPS)
LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIB_MTX=$(LIBXMTX_KSP)
#LIBX_MTX=$(LIBX_MTX_KSP)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE= -I../pl/$(MOD) -I../lib/$(MOD) -I../eq/$(MOD) \
            -I../../bpsd/$(MOD) -I../mtxp/$(MOD) -I $(MOD)

SRCS = wfcomm.f90 wfmenu.f90 wfinit.f90 wfdiv.f90 wfant.f90 wfwave.f90 \
       wfsolv.f90 wfprof.f90 wfsub.f90 wfnas.f90 wfsort.f90 wffile.f90 \
       wflib.f90 wfindx.f90 wfefwg.f90 wfbes.f90 wfgout.f90 wfmag.f90 \
       wfgsub.f90 wfgutl.f90 mbant.f90
OBJS = $(SRCS:.f90=.o)

LIBS = libwf.a ../pl/libpl.a ../eq/libeq.a   \
        ../lib/libtask.a ../../bpsd/libbpsd.a $(LIBMTX)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: wf

libwf.a: $(OBJS)
	$(LD) $(LDFLAGS) libwf.a $(OBJS)

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../mtxp/libmtxbnd.o:
	(cd ../mtxp; make libmtxbnd.o)
../mtxp/libmtxmumps.o:
	(cd ../mtxp; make libmtxmumps.o)
../mtxp/libmtxmpi.o:
	(cd ../mtxp; make libmtxmpi.o)
../mtxp/libmtxnompi.o:
	(cd ../mtxp; make libmtxnompi.o)
../pl/libpl.a:
	(cd ../../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)

wf : $(LIB_MTX) $(LIBS) wfmain.o
	$(FLINKER) wfmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

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
#	$(FCFREE) $(LFLAGS) -c $< 
wfgdum.o: wfgdum.f90
#	$(FCFREE) $(LFLAGS) -c $< 
wfgutl.o: wfgutl.f90
#	$(FCFREE) $(LFLAGS) -c $< 
wflib.o: wflib.f90 wfcomm.f90
#	$(FCFREE) $(LFLAGS) -c $< $(MODDIR) $(MODINCLUDE)
