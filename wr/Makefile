### Madefile for wr
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

MODINCLUDE= -I./$(MOD) -I../pl/$(MOD) -I../dp/$(MOD) -I../lib/$(MOD) \
            -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

SRCS = wrcomm.f90 wrinit.f90 wrparm.f90 \
       wrsub.f90 wrexecr.f90 wrcalc.f90 wrbeam.f90 wrgout.f90 wrfile.f90 \
       wrmenu.f90

OBJS = $(SRCS:.f90=.o)

LIBS = libwr.a  ../dp/libdp.a ../eq/libeq.a ../pl/libpl.a \
       ../lib/libtask.a ../../bpsd/libbpsd.a

all : wr

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../dp/libdp.a:
	(cd ../dp; make libdp.a)
libwr.a: $(OBJS)
	$(LD) $(LDFLAGS) libwr.a $(OBJS)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)
	(cd ../dp; make libdp.a)

wr : $(LIBS) wrmain.o
	$(FLINKER) wrmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f wr

WRCOMM=wrcomm.f90 ../dp/dpcomm.f90 ../pl/plcomm.f90
wrmain.o: wrmain.f90 $(WRCOMM)
wrmenu.o: wrmenu.f90 $(WRCOMM)
wrinit.o: wrinit.f90 $(WRCOMM)
wrparm.o: wrparm.f90 $(WRCOMM)
wrcalc.o: wrcalc.f90 $(WRCOMM) ../pl/plprof.f90
wrexecr.o:wrexecr.f90 $(WRCOMM) ../pl/plprof.f90
wrbeam.o: wrbeam.f90 $(WRCOMM) ../pl/plprof.f90
wrsub.o:  wrsub.f90 $(WRCOMM)
wrgout.o: wrgout.f90 $(WRCOMM)
wrfile.o: wrfile.f90 $(WRCOMM)
