### $Id$
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(DFLAGS)
FFLAGS = $(OFLAGS)

MODINCLUDE=-I ../lib/mod -I $(MOD)

SRCS = wfinit.f wfdiv.f wfant.f wfwave.f wfsolv.f \
       wfprof.f wfsub.f wfnas.f wfsort.f wffile.f wflib.f \
       wfindx.f wfefwg.f wfgout.f wfgsub.f

OBJS = $(SRCS:.f=.o)

LIBS=libwf.a ../lib/libtask.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: wf

libwf.a: $(OBJS)
	$(LD) $(LDFLAGS) libwf.a $(OBJS)

../lib/libtask.a:
	(cd ../lib; make libtask.a)

wf : wfcomm.inc $(LIBS) wfmain.o
	$(FCFIXED) wfmain.o -o $@ $(FFLAGS) $(LIBS) $(FLIBS)

testin: testin.f wfxnas.f wfcomm.inc
	$(FCFIXED) testin.f wfxnas.f -o $@ $(FFLAGS) $(FLIBS)

nascheck: nascheck.f
	$(FCFIXED) nascheck.f -o $@ $(FFLAGS) $(FLIBS)

clean:
	rm -f *.o ./#* ./*~ a.out core *.a *.mod $(MOD)/*.mod

veryclean: clean
	rm -f wf testin nascheck

wfmain.o: wfmain.f wfcomm.inc
wfinit.o: wfinit.f wfcomm.inc
wfdiv.o: wfdiv.f wfcomm.inc
wfant.o: wfant.f wfcomm.inc
wfwave.o: wfwave.f wfcomm.inc
wfsolv.o: wfsolv.f wfcomm.inc
wfprof.o: wfprof.f wfcomm.inc
wfsub.o: wfsub.f wfcomm.inc
wfnas.o: wfnas.f wfcomm.inc
wfsort.o: wfsort.f wfcomm.inc
wffile.o: wffile.f wfcomm.inc
wfindx.o: wfindx.f wfcomm.inc
wfgout.o: wfgout.f wfcomm.inc
wfgsub.o: wfgsub.f wfcomm.inc
wfefwg.o: wfefwg.f wfcomm.inc
wfgdum.o: wfgdum.f
wflib.o: wflib.f
