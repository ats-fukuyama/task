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

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I $(MOD) -I../plx/$(MOD) -I../lib/$(MOD) -I../../bpsd/$(MOD) -I../mtxp/$(MOD) 

SRCS = dpcomm.f90 dpinit.f90 dpfpin.f90 dpparm.f90 \
       dphotf.f90 dphotr.f90 \
       dptnsr1.f90 dptnsr2.f90 dptens.f90 dpdisp.f90 \
       dpgsub.f90 dproot.f90 dpcont.f90 \
       dpmenu.f90	
	
SRCSX = dpdkdt.f dptestk.f


OBJS = $(SRCS:.f90=.o)

LIBS = ../lib/tasklib.a ../../bpsd/bpsdlib.a ../plx/pllib.a \
       ../eq/eqlib.a dplib.a 

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
	$(FLINKER) dpmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f dp

DPCOMM=dpcomm.f90 ../plx/plcomm.f90

dpcomm.o: dpcomm.f90 ../plx/plcomm.f90
dpinit.o: dpinit.f90 $(DPCOMM)
dpparm.o: dpparm.f90 $(DPCOMM)
dproot.o: dproot.f90 $(DPCOMM)
dpcont.o: dpcont.f90 $(DPCOMM)
dpscan.o: dpscan.f $(DPCOMM)
dpdisp.o: dpdisp.f90 $(DPCOMM) ../plx/plprof.f90
dptens.o: dptens.f90 $(DPCOMM)
dptnsr1.o: dptnsr1.f90 $(DPCOMM)
dptnsr2.o: dptnsr2.f90 $(DPCOMM)
dpfpin.o: dpfpin.f90 $(DPCOMM) ../plx/plprof.f90
dphotf.o: dphotf.f90 $(DPCOMM)
dphotr.o: dphotr.f90 $(DPCOMM)
dpdkdt.o: dpdkdt.f90 $(DPCOMM) ../plx/plprof.f90
dpgsub.o: dpgsub.f90
dpmenu.o: dpmenu.f90 $(DPCOMM) ../plx/plprof.f90
dpmain.o: dpmain.f90 $(DPCOMM)
