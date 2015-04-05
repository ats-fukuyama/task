### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f90 .f .mod .o .a

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)

MODINCLUDE=-I ../lib/mod -I $(MOD)

SRCS = wfinit.f wfdiv.f wfzone.f wfant.f \
       wfwave.f wfevol.f wfgout.f \
       wffreq.f wfprof.f wffile.f wfsub.f wfgsub.f wflib.f \
       wfatmd.f

OBJS = $(SRCS:.f=.o)

LIBS = wflib.a ../lib/tasklib.a

.f90.o:
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f.o:
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: wf

wflib.a: $(OBJS)
	$(LD) $(LDFLAGS) wflib.a $(OBJS)

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)

wf : wfcomm.inc $(LIBS) wfmain.o
	$(FCFIXED) wfmain.o -o $@ $(FFLAGS) $(LIBS) $(FLIBS)

clean:
	rm -f *.o ./#* ./*~ a.out core *.a *.mod $(MOD)/*.mod

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
