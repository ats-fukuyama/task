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

MODINCLUDE= -I $(MOD) -I../pl/$(MOD) -I../lib/$(MOD) -I../../bpsd/$(MOD) -I../mtxp/$(MOD) 

SRCS = dpcomm.f90 dpinit.f90 dpfpin.f90 dpparm.f90 \
       dphotf.f90 dphotr.f90 \
       dptnsr1.f90 dptnsr2.f90 dptnsr3.f90 dptnsr4.f90 dptnsr0.f90 dpdisp.f90 \
       dpgsub.f90 dpglib.f90 dproot.f90 dpcont.f90 dpcont4.f90 \
       dpmenu.f90	
	
SRCSX = dpdkdt.f dptestk.f


OBJS = $(SRCS:.f90=.o)

LIBS = ../lib/libtask.a ../../bpsd/libbpsd.a ../pl/libpl.a \
       ../eq/libeq.a libdp.a 

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
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)

dp : $(LIBS) dpmain.o
	$(FLINKER) dpmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f dp

DPCOMM=dpcomm.f90
PLPROF=../pl/plprof.f90 ../pl/plprofw.f90

dpcomm.o: dpcomm.f90 ../pl/plcomm.f90
dpinit.o: dpinit.f90 $(DPCOMM)
dpparm.o: dpparm.f90 $(DPCOMM)
dproot.o: dproot.f90 $(DPCOMM)
dpcont.o: dpcont.f90 $(DPCOMM) $(PLPROF)
dpcont4.o: dpcont4.f90 $(DPCOMM) $(PLPROF)
dpscan.o: dpscan.f $(DPCOMM)
dpdisp.o: dpdisp.f90 $(DPCOMM) $(PLPROF)
dptnsr0.o: dptnsr0.f90 $(DPCOMM) $(PLPROF)
dptnsr1.o: dptnsr1.f90 $(DPCOMM) $(PLPROF)
dptnsr2.o: dptnsr2.f90 $(DPCOMM) $(PLPROF)
dptnsr3.o: dptnsr3.f90 $(DPCOMM) $(PLPROF)
dptnsr4.o: dptnsr4.f90 $(DPCOMM) $(PLPROF)
dpfpin.o: dpfpin.f90 $(DPCOMM) $(PLPROF)
dphotf.o: dphotf.f90 $(DPCOMM) $(PLPROF)
dphotr.o: dphotr.f90 $(DPCOMM) $(PLPROF)
dpdkdt.o: dpdkdt.f90 $(DPCOMM)
dpgsub.o: dpgsub.f90
dpmenu.o: dpmenu.f90 $(DPCOMM) $(PLPROF)
dpmain.o: dpmain.f90 $(DPCOMM)
