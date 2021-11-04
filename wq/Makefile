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

SRCS = wqcomm.f90 \
       wqinit.f90 wqparm.f90 wqview.f90 \
       wqtens.f90 wqprep.f90 wqsolv1.f90 wqsolv2.f90 wqexec.f90 \
       wqgout.f90 \
       wqmenu.f90

OBJS = $(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))

LIBS = libwq.a ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a

all: wq

libwq.a: $(OBJS)
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

wq: $(LIB_MTX) $(LIBS) $(OBJDIR)/wqmain.o
	$(FLINKER) $(OBJDIR)/wqmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out *.o *.a ./#* ./*~
	-rm -f ./mod/*.mod
	-rm -f ./obj/*.o
	-rm -rf ./*.dSYM

veryclean: clean
	-rm -f wq

$(OBJDIR)/wqcomm.o: wqcomm.f90
$(OBJDIR)/wqinit.o: wqinit.f90 wqcomm.f90
$(OBJDIR)/wqparm.o: wqparm.f90 wqcomm.f90
$(OBJDIR)/wqview.o: wqview.f90 wqcomm.f90
$(OBJDIR)/wqtens.o: wqtens.f90 wqcomm.f90
$(OBJDIR)/wqprep.o: wqprep.f90 wqtens.f90 wqcomm.f90
$(OBJDIR)/wqsolv1.o: wqsolv1.f90 wqcomm.f90
$(OBJDIR)/wqsolv2.o: wqsolv2.f90 wqcomm.f90
$(OBJDIR)/wqexec.o: wqexec.f90 wqsolv1.f90 wqsolv2.f90 wqcomm.f90
$(OBJDIR)/wqgout.o: wqgout.f90 wqcomm.f90
$(OBJDIR)/wqmenu.o: wqmenu.f90 wqparm.f90 wqview.f90 wqexec.f90 wqgout.f90 wqcomm.f90
$(OBJDIR)/wqmain.o: wqmain.f90 wqinit.f90 wqparm.f90 wqmenu.f90 wqcomm.f90

