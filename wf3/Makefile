### $Id$
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(DFLAGS)
FFLAGS = $(OFLAGS)

SRCS = wfinit.f wfdiv.f wfant.f wfwave.f wfsolv.f \
       wfprof.f wfsub.f wfnas.f wfsort.f wffile.f wflib.f \
       wfindx.f wfefwg.f wfbes.f wfgout.f wfgsub.f wfgutl.f

OBJS = $(SRCS:.f=.o)

LIBS=wflib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: wf

wflib.a: $(OBJS)
	$(LD) $(LDFLAGS) wflib.a $(OBJS)

wf : wfcomm.inc $(LIBS) wfmain.o
	$(FCFIXED) wfmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

testin: testin.f wfxnas.f wfcomm.inc
	$(FCFIXED) testin.f wfxnas.f -o $@ $(FFLAGS) $(FLIBS)

nascheck: nascheck.f
	$(FCFIXED) nascheck.f -o $@ $(FFLAGS) $(FLIBS)

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
	$(FCFIXED) $(LFLAGS) -c $< 
wfgutl.o: wfgutl.f
	$(FCFIXED) $(LFLAGS) -c $< 
wflib.o: wflib.f
	$(FCFIXED) $(LFLAGS) -c $< 
wfbes.o: wfbes.f
	$(FCFIXED) $(LFLAGS) -c $< 
