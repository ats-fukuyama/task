### Makefile for wq 2021-11-04 ###

include ../make.header
include ../mtxp/make.mtxp

#LIB_MTX=$(LIB_MTX_BND)
#LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_PCG)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIB_MTX=$(LIBXMTX_KSP)
#LIBX_MTX=$(LIBX_MTX_KSP)

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

OBJDIR   = ./obj

${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../lib/$(MOD) -I../mtxp/$(MOD) \
            -I../../bpsd/$(MOD)

SRCS = libpgkf.f90 \
       libpdkf.f90

OBJS = $(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))

LIBS = libpkf.a ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a

all: libs plotpkf

libpkf.a: $(OBJS)
	$(LD) $(LDFLAGS) $@ $(OBJS)

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

plotpkf: $(LIB_MTX) $(LIBS) $(OBJDIR)/plotpkf.o
	$(FLINKER) $(OBJDIR)/plotpkf.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out *.o *.a ./#* ./*~
	-rm -f ./mod/*.mod
	-rm -f ./obj/*.o
	-rm -rf ./*.dSYM

veryclean: clean
	-rm -f wq

$(OBJDIR)/plotpkf.o: plotpkf.f90 libpdkf.f90 libpgkf.f90
$(OBJDIR)/libpdkf.o: libpdkf.f90
$(OBJDIR)/libpgkf.o: libpgkf.f90

