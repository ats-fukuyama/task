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

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

MODINCLUDE= -I $(MOD) -I../plx/$(MOD) -I../lib/$(MOD) -I../../bpsd/$(MOD) -I../mtxp/$(MOD) 

SRCS = dpmenu.f dpinit.f dproot.f dpcont.f dpdisp.f \
       dptens.f dptnsr1.f dptnsr2.f \
       dpfpin.f dphotf.f dphotr.f dpdkdt.f dpgsub.f 

OBJS = $(SRCS:.f=.o) 

LIBS = $(LIBX_MTX) ../lib/tasklib.a ../../bpsd/bpsdlib.a ../plx/pllib.a \
       ../eq/eqlib.a dplib.a 

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : dp

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
dplib.a: $(OBJS)
	$(LD) $(LDFLAGS) dplib.a $(OBJS)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../plx; make pllib.a)
	(cd ../eq; make eqlib.a)

dp : $(LIBS) dpmain.o
	$(FLINKER) dpmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core a.out *.o *.a ./*~

veryclean: clean
	-rm -f dp

new:
	-mkdir ../dpnew
	cp -f Makefile ../dpnew
	cp -f *.f ../dpnew
	cp -f *.inc ../dpnew

DPCOMM=dpcomm.inc dpcom1.inc dpcom2.inc ../plx/plcomm.f90

dpmain.o: dpmain.f $(DPCOMM)
dpmenu.o: dpmenu.f $(DPCOMM)
dpinit.o: dpinit.f $(DPCOMM)
dproot.o: dproot.f $(DPCOMM)
dpcont.o: dpcont.f $(DPCOMM) ../plx/plprof.f90
dpdisp.o: dpdisp.f $(DPCOMM) ../plx/plprof.f90
dptens.o: dptens.f $(DPCOMM) ../plx/plprof.f90
dptnsr1.o: dptnsr1.f $(DPCOMM) ../plx/plprof.f90
dptnsr2.o: dptnsr2.f $(DPCOMM) ../plx/plprof.f90
dpfpin.o: dpfpin.f $(DPCOMM) ../plx/plprof.f90
dphotf.o: dphotf.f $(DPCOMM) ../plx/plprof.f90
dphotr.o: dphotr.f $(DPCOMM) ../plx/plprof.f90
dpdkdt.o: dpdkdt.f $(DPCOMM) ../plx/plprof.f90
dpgsub.o: dpgsub.f
