### Makefile for task/imas

include ../../make.header
include ../../mtxp/make.mtxp

#LIB_MTX=$(LIB_MTX_NOMPI)
#LIBX_MTX=$(LIBX_MTX_NOMPI)
#LIB_MTX=$(LIB_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_PCG)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIB_MTX=$(LIBXMTX_KSP)
#LIBX_MTX=$(LIBX_MTX_KSP)

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MOD=mod
MODINCLUDE= -I./$(MOD) \
            -I../../pl/$(MOD) -I../../eq/$(MOD) \
            -I../../wrx/$(MOD) \
            -I../../lib/$(MOD) -I../../mtxp/$(MOD) -I../../../bpsd/$(MOD)

LIB_AL=`pkg-config --libs --cflags al-fortran`
F90INCMOD=$(LIB_AL) $(MODINCLUDE)

OBJDIR=./obj
${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(F90INCMOD)

# LIST OF FORTRAN FILES
SRCS = ids_equilibrium.f90 \
       ids_core_profile.f90 \
       ids_ec_launchers.f90 \
       ids_waves.f90

OBJS=$(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))

LIBS=../../wrx/libwr.a

all: libs

# LIBRARY COMPILATION

libids.a: $(OBJS)
	$(LD) $(LDFLAGS) libids.a $(OBJS)

libs: libids.a

# CLEAN DIRECTORY
clean:
	rm -f ./#* ./*~ $(MOD)/*.mod $(OBJDIR)/*.o *.a

$(OBJDIR)/ids_equilibrium.o : ids_equilibrium.f90
$(OBJDIR)/ids_core_profile.o : ids_core_profile.f90
$(OBJDIR)/ids_ec_launchers.o : ids_ec_launchers.f90
$(OBJDIR)/ids_waves.o : ids_waves.f90


