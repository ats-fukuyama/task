### Makefile for wm
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

MODINCLUDE= -I./$(MOD) -I../plx/$(MOD) -I../dpx/$(MOD) \
            -I../lib/$(MOD) -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

SRCS = wmcomm.f90 vmcomm.f90 \
       wminit.f90 wmparm.f90 wmview.f90 \
       wmmenu.f90

OBJS = $(SRCS:.f90=.o)

LIBS = libwm.a ../dpx/dplib.a ../eq/eqlib.a ../plx/pllib.a \
       ../lib/tasklib.a ../../bpsd/bpsdlib.a

all: libs wm

../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../dpx/dplib.a:
	(cd ../dpx; make dplib.a)
$(LIB_MTX):
	(cd ../mtxp; make)
libwm.a: $(OBJS)
	$(LD) $(LDFLAGS) $@ $(OBJS)

libs:
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../lib; make tasklib.a)
	(cd ../mtxp; make)
	(cd ../plx; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../dpx; make dplib.a)

wm : $(LIB_MTX) $(LIBS) wmmain.o 
	$(FLINKER) wmmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean :
	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod

veryclean: clean
	-rm -f wm

wmcomm.o: wmcomm.f90
vmcomm.o: vmcomm.f90

wminit.o: wminit.f90 wmcomm.f90
wmparm.o: wmparm.f90 wmcomm.f90
wmview.o: wmview.f90 wmcomm.f90

wmmenu.0: wmmenu.f90 wmcomm.f90
wmmain.0: wmmain.f90 wmcomm.f90
