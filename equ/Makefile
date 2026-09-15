### Makefile

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

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS) -Wuninitialized

OBJDIR=./obj

${OBJDIR}/%.o:%.f
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../pl/$(MOD) \
            -I../lib/$(MOD) -I../mtxp/$(MOD) \
            -I../../bpsd/$(MOD)


SRCC =  equcomm.f90
SRCS =  equtpxssl.f90 equpfds.f90 equsub.f90 equinit.f90 \
        equfct.f90 equinp.f90 equset.f90 \
        equpl.f90 equtrn.f90 equfile.f90 \
	equgout.f90 equunit.f90 \
	equtrpl.f90 equtrset.f90 equtradv.f90 equtrunit.f90 \
        equ_eqdsk.f90 equ_convert.f90 \
        equmenu.f90

OBJC=$(addprefix $(OBJDIR)/, $(SRCC:.f90=.o))
OBJS=$(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))

LIBS = libequ.a ../pl/libpl.a ../eq/libeq.a ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a

all: libs equ

libequ.a: $(OBJS)
	$(LD) $(LDFLAGS) libequ.a $(OBJS)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib/; make libtask.a)
	(cd ../lib/; make libgrf.a)
	(cd ../mtxp/; make)
	(cd ../pl/; make libpl.a)

equ: $(OBJC) libs libequ.a $(OBJDIR)/equmain.o
	$(FLINKER) $(OBJDIR)/equmain.o $(OBJC) $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)
clean:
	-rm -f ./#* ./*~ $(OBJDIR)/*.o $(MOD)/*.mod core *.a
veryclean: clean
	-rm -f equ
#
$(OBJDIR)/equmain.o : equmain.f90
$(OBJDIR)/equmenu.o : equmenu.f90
$(OBJDIR)/equunit.o : equunit.f90 equcomm.f90
$(OBJDIR)/equfct.o : equfct.f90 equcomm.f90
$(OBJDIR)/equinp.o : equinp.f90 equcomm.f90
$(OBJDIR)/equinit.o : equinit.f90 equcomm.f90
$(OBJDIR)/equset.o : equset.f90 equcomm.f90
$(OBJDIR)/equsub.o : equsub.f90 equcomm.f90
$(OBJDIR)/equpfds.o : equpfds.f90 equcomm.f90
$(OBJDIR)/equintf.o : equintf.f90 equcomm.f90
$(OBJDIR)/equgout.o : equgout.f90 equcomm.f90
$(OBJDIR)/equfile.o : equgout.f90 equcomm.f90
$(OBJDIR)/equtrunit.o : equtrunit.f90  
$(OBJDIR)/equtrpl.o : equtrpl.f90  
$(OBJDIR)/equtrset.o : equtrset.f90  
$(OBJDIR)/equtradv.o : equtradv.f90 equtrpl.f90
$(OBJDIR)/equtrsvq.o : equtrsvq.f90 equtrpl.f90
$(OBJDIR)/equtpxssl.o : equtpxssl.f90
$(OBJDIR)/equ_eqdsk.o : equ_eqdsk.f90 equcomm.f90
$(OBJDIR)/equ_convert.o : equ_convert.f90 equcomm.f90
