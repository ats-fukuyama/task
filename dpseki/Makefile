### $Id: Makefile,v 1.48 2013/01/24 03:52:13 fukuyama Exp $ ###
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

MODINCLUDE= -I $(MOD) -I../pl/$(MOD) -I../lib/$(MOD) -I../../bpsd/$(MOD)

SRCS = dpmenu.f dpinit.f dproot.f dpcont.f dpdisp.f \
       dptens.f dptnsr1.f dptnsr2.f \
       dpfpin.f dphotf.f dphotr.f dpdkdt.f dpgsub.f 

OBJS = $(SRCS:.f=.o) 

LIBS = ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a ../pl/libpl.a \
       ../eq/libeq.a libdp.a 

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : dp

../lib/libtask.a:
	(cd ../lib; make libtask.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
libdp.a: $(OBJS)
	$(LD) $(LDFLAGS) libdp.a $(OBJS)

libs:
	(cd ../lib; make libtask.a)
	(cd ../../bpsd; make libbpsd.a)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)

dp : $(LIBS) dpmain.o
	$(FLINKER) dpmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

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

DPCOMM=dpcomm.inc dpcom1.inc dpcom2.inc ../pl/plcomm.f90

dpmain.o: dpmain.f $(DPCOMM)
dpmenu.o: dpmenu.f $(DPCOMM)
dpinit.o: dpinit.f $(DPCOMM)
dproot.o: dproot.f $(DPCOMM)
dpcont.o: dpcont.f $(DPCOMM) ../pl/plprof.f90
dpdisp.o: dpdisp.f $(DPCOMM) ../pl/plprof.f90
dptens.o: dptens.f $(DPCOMM) ../pl/plprof.f90
dptnsr1.o: dptnsr1.f $(DPCOMM) ../pl/plprof.f90
dptnsr2.o: dptnsr2.f $(DPCOMM) ../pl/plprof.f90
dpfpin.o: dpfpin.f $(DPCOMM) ../pl/plprof.f90
dphotf.o: dphotf.f $(DPCOMM) ../pl/plprof.f90
dphotr.o: dphotr.f $(DPCOMM) ../pl/plprof.f90
dpdkdt.o: dpdkdt.f $(DPCOMM) ../pl/plprof.f90
dpgsub.o: dpgsub.f
