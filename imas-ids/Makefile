### Makefile for task/imas

include ../make.header
include ../mtxp/make.mtxp

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
            -I../pl/$(MOD) -I../eq/$(MOD) -I../dp/$(MOD) \
            -I../wrx/$(MOD) \
            -I../lib/$(MOD) -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

#LIB_AL=`pkg-config --libs --cflags al-fortran`
#-Wl,--defsym,DD_VER_3.41.0=0
#-Wl,--defsym,AL_VER_5.2.4=0

LIB_AL=-I/Users/fukuyama/imas/al-install/lib/pkgconfig/../../include/fortran -I/Users/fukuyama/imas/al-install/lib/pkgconfig/../../include -L/Users/fukuyama/imas/al-install/lib/pkgconfig/../../lib -lal-fortran-3.41.0  -Wl,-rpath,/Users/fukuyama/imas/al-install/lib/pkgconfig/../../lib -lal -Wl,-rpath,/Users/fukuyama/imas/al-install/lib/pkgconfig/../../lib


F90INCMOD=$(LIB_AL) $(MODINCLUDE)

OBJDIR=./obj
${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(F90INCMOD)

# LIST OF FORTRAN FILES
SRCS = task_ids_equilibrium.f90 \
       task_ids_core_profiles.f90 \
       task_ids_ec_launchers.f90 \
#       task_ids_waves.f90

OBJS=$(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))

LIBS=  libids.a ../wrx/libwr.a ../dp/libdp.a ../eq/libeq.a ../pl/libpl.a \
       ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a $(LIB_AL)

all: libs test

# LIBRARY COMPILATION

libids.a: $(OBJS)
	$(LD) $(LDFLAGS) libids.a $(OBJS)

libs: libids.a

# test ids

testids: libs $(OBJDIR)/testids.o
	echo $(FLINKER)
	echo $(LIBS)
	echo $(FFLAGS)
	echo $(FLIBS)
	echo $(LIBX_MTX)
	$(FLINKER) $(OBJDIR)/testids.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)


# CLEAN DIRECTORY
clean:
	rm -f ./#* ./*~ $(MOD)/*.mod $(OBJDIR)/*.o *.a

$(OBJDIR)/task_ids_equilibrium.o : task_ids_equilibrium.f90
$(OBJDIR)/task_ids_core_profiles.o : task_ids_core_profiles.f90
$(OBJDIR)/task_ids_ec_launchers.o : task_ids_ec_launchers.f90
$(OBJDIR)/task_ids_waves.o : task_ids_waves.f90


