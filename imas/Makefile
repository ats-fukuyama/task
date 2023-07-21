# DEFINE F90
include compiler.mk
# SET COMPILER OPTIONS
ifeq ($(F90), ifort)
# INTEL - LINKS TO THE IMAS LIBRARY AND INCLUDE DIRECTORY
  F90FLAGS=-fPIC -fpp -g # FPIC AND PREPROCESSING OPTIONS
  F90INC=-I. `pkg-config imas-ifort --cflags` -I../lib
  F90LIB=`pkg-config imas-ifort --libs`
else ifeq ($(F90), gfortran)
# GFORTRAN - LINKS TO THE IMAS LIBRARY AND INCLUDE DIRECTORY
  F90FLAGS=-fPIC -cpp # FPIC AND PREPROCESSING OPTIONS
  F90INC=-I. `pkg-config imas-gfortran --cflags`
  F90LIB=`pkg-config imas-gfortran --libs`
else
  $(warning Unsupported Fortran compiler $(F90), proceed with care...)
  F90FLAGS=-fPIC -cpp # FPIC AND PREPROCESSING OPTIONS
  F90INC=-I. `pkg-config imas-gfortran --cflags`
  F90LIB=`pkg-config imas-gfortran --libs`
endif

MOD=mod
MODINCLUDE = -I./$(MOD) -I../lib/$(MOD)
F90INCMOD=$(F90INC) $(MODINCLUDE)

OBJDIR=./obj

# LIST OF FORTRAN FILES
OBJS_ACTOR = $(OBJDIR)/wr_imas.o
OBJS_STDA  = $(OBJDIR)/wr_standalone.o

all: lib exe actor

lib: libwr_imas.a

exe: wr_standalone.exe

orig: wr_original.exe

# STANDALONE COMPILATION
wr_standalone.exe: libwr_imas.a ${OBJS_STDA}
	$(F90) ${F90FLAGS} -o wr_standalone.exe ${OBJS_STDA} $(F90LIB) -L. -lwr_imas

# LIBRARY COMPILATION
libwr_imas.a: ${OBJS_ACTOR}
	ar -cr $@ $^

# COMPILE THE ACTOR
task_wr_actor: libwr_imas_wr.a
	sed 's/__COMPILER__/${F90}/g'task_wr_template.xml > task_wr.xml
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
