###
### This Makefile should be located at task/ti/adpost
###
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

.SUFFIXES:
.SUFFIXES: .f .for .f90 .mod .o .a

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE= -I./$(MOD) -I../../bpsd/$(MOD) -I../mtxp/$(MOD) -I../lib/$(MOD)

SRCS = lib-adpost.f90

OBJS = $(SRCS:.f90=.o)

LIBS = lib-adpost.a ../mtxp/libmtxmpi.o ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a

.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : libs test-adpost

lib-adpost.a: $(OBJS)
	$(LD) $(LDFLAGS) $@ $(OBJS)

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../mtxp; make)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)

test-adpost : $(LIBS) test-adpost.o
	$(FLINKER) test-adpost.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MODDIR) $(MODINCLUDE)

clean:
	-rm -f core *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f test-adpost

lib-adpost.o: lib-adpost.f90
test-adpost.o: test-adpost.f90 lib-adpost.f90
