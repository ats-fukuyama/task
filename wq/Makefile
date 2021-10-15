### $Id: Makefile,v 1.5 2010/05/18 03:57:54 maruyama Exp $
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

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../lib/$(MOD) -I../mtxp/$(MOD) \
            -I../../bpsd/$(MOD)

SRCS = wqcomm.f90 \
       wqinit.f90 wqparm.f90 wqview.f90 \
       wqtens.f90 wqprep.f90 wqsolv0.f90 wqsolv1.f90 wqsolv2.f90 wqexec.f90 \
       wqgout.f90 \
       wqmenu.f90

OBJS = $(SRCS:.f90=.o)

LIBS = libwq.a ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a

all: wq

libwq.a: $(OBJS)
	$(LD) $(LDFLAGS) libwq.a $(OBJS)

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../task/lib/libtask.a:
	(cd ../lib; make libtask.a)
../task/lib/libgrf.a:
	(cd ../lib; make libgrf.a)
$(LIB_MTX):
	(cd ../mtxp; make)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../mtxp; make)

wq: $(LIB_MTX) $(LIBS) wqmain.o
	$(FLINKER) wqmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	rm -f core *.o ./*~ 

veryclean: clean
	rm -f wq

wqcomm.o: wqcomm.f90
wqinit.o: wqinit.f90 wqcomm.f90
wqparm.o: wqparm.f90 wqcomm.f90
wqview.o: wqview.f90 wqcomm.f90
wqtens.o: wqtens.f90 wqcomm.f90
wqprep.o: wqprep.f90 wqtens.f90 wqcomm.f90
wqsolv0.o: wqsolv0.f90 wqcomm.f90
wqsolv1.o: wqsolv1.f90 wqcomm.f90
wqsolv2.o: wqsolv2.f90 wqcomm.f90
wqexec.o: wqexec.f90 wqsolv0.f90 wqcomm.f90
wqgout.o: wqgout.f90 wqcomm.f90
wqmenu.o: wqmenu.f90 wqparm.f90 wqview.f90 wqexec.f90 wqgout.f90 wqcomm.f90
wqmain.o: wqmain.f90 wqinit.f90 wqparm.f90 wqmenu.f90 wqcomm.f90

