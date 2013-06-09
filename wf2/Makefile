### $Id$ ###
include ../make.header
#include make.header

.SUFFIXES:
.SUFFIXES: .f90 .f .mod .o .a

MODINCLUDE=

.f90.o:
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f.o:
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

LFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

SRCS = wfmain.f wfinit.f wfdiv.f wfzone.f wfant.f \
       wfwave.f wfevol.f wfgout.f \
       libspl.f \
       wffreq.f wfprof.f wffile.f wfsub.f wfgsub.f wflib.f \
       wfatmd.f

OBJS = $(SRCS:.f=.o)

wf : $(OBJS) wfcomm.inc
	$(FCFIXED) $(OBJS) -o $@ $(FFLAGS) $(FLIBS)

check : 
	ftnchek -array=0 -usage=no-com-block-unused $(SRCS) | less
new :
	mkdir ../wfnew
	cp -f Makefile ../wfnew
	cp -f make.header.org ../wfnew
	cp -f *.f ../wfnew
	cp -f *.inc ../wfnew
	cp -f HISTORY ../wfnew
	mkdir ../wfnew/doc
	cp -f doc/*.tex ../wfnew/doc
	cp -f doc/*.sty ../wfnew/doc

clean:
	rm -f core *.o *.P ./*~

veryclean: clean
	rm -f wf

wfmain.o: wfmain.f wfcomm.inc
wfinit.o: wfinit.f wfcomm.inc
wfdiv.o: wfdiv.f wfcomm.inc
wfdivn.o: wfdivn.f wfcomm.inc
wfzone.o: wfzone.f wfcomm.inc
wfant.o: wfant.f wfcomm.inc
wfwave.o: wfwave.f wfcomm.inc
wfevol.o: wfevol.f wfcomm.inc
wfgout.o: wfgout.f wfcomm.inc
wffreq.o: wffreq.f wfcomm.inc
wfprof.o: wfprof.f wfcomm.inc
wffile.o: wffile.f wfcomm.inc
wfsub.o: wfsub.f wfcomm.inc
wfgsub.o: wfgsub.f wfcomm.inc
wflib.o: wflib.f
	$(FCFIXED) $(LFLAGS)  -c $< 
libspl.o: libspl.f
	$(FCFIXED) $(LFLAGS)  -c $< 
