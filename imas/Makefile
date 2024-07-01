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
            -I../pl/$(MOD) -I../eq/$(MOD) \
            -I../lib/$(MOD) -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

OBJDIR=./obj
${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

LIB_AL=`pkg-config --libs --cflags al-fortran`
F90INCMOD=$(LIB_AL) $(MODINCLUDE)

# LIST OF FORTRAN FILES
OBJS_ACTOR = $(OBJDIR)/wr_imas.o
OBJS_STDA  = $(OBJDIR)/wr_standalone.o


all: lib exe actor

# LIBRARY COMPILATION
lib: libwr_imas.a
libwr_imas.a: ${OBJS_ACTOR}
	ar -cr $@ $^

# STANDALONE COMPILATION
exe: wr_standalone.exe
wr_standalone.exe: libwr_imas.a ${OBJS_STDA}
	$(F90) ${F90FLAGS} -o wr_standalone.exe ${OBJS_STDA} $(F90LIB) -L. -lwr_imas

orig: wr_original.exe



# COMPILE THE ACTOR
task_wr_actor: libwr_imas_wr.a
	sed 's/__COMPILER__/${F90}/g' task_wr_template.xml > task_wr.xml
	$(FC2K)/bin/fc2k task_wr.xml -nokepler -pyworkspace $(ACTOR_FOLDER)

actor: libphysics_model_level1.a 
	sed 's/__COMPILER__/${F90}/g' task_wr_template.yaml > task_wr.yaml
	iwrap -f task_wr.yaml -i $(ACTOR_FOLDER)

# COMPILE THE ORIGINAL CODE
input_data.txt:
	echo "3.14156" >> $@
wr_original.exe: wr_original.f90 input_data.txt
	$(F90) -o $@ $<

# RULES FOR ALL OBJECT FILES
${OBJDIR}/%.o:%.f90
	$(F90) ${F90FLAGS} -c  $< -o $@ $(F90INCMOD)

# CLEAN DIRECTORY
clean:
	rm -f *~ *.exe *.a *.mod *.o task_wr.xml task_wr.exe input_data.txt output_data.txt
