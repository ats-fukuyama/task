#=======================================
# MMM - Secondary Makefile for libmmm7_1
#=======================================
#
# Lixiang Luo, 14-Oct-2010
# lixiang.luo@lehigh.edu
# Department of Physics, Lehigh University
#
# A component of the MMM7.1 package
#
# You need to change the MMMFC and FFLAGS according to your
# compiler requirements.
#
# Default target is libmmm.a, which doesn't include debug info.
# To build default target:
#
# $ make
#
# To clean up binaries:
#
# $ make clean
#
# To get a build with debug info:
#
# $ make clean
# $ make debug

lib_mmm := libmmm7_1.a # Final target
modules := w20mod.o
objects := \
   modmmm7_1.o

MMMFC  ?= gfortran # default to gfortran for maximum compatibility
FFLAGS := -O3 # Optimization
#FFLAGS := -O3 -Wall -fbounds-check -ffpe-trap=invalid,zero,overflow -finit-real=nan

all: $(lib_mmm)

.PHONY : debug
debug : FFLAGS := $(subst -O3,-g,$(FFLAGS))
debug : $(lib_mmm)

$(objects): $(modules) # Modules must be built first

# Generate the static-link library
$(lib_mmm): $(modules) $(objects)
	ar rcs $@ $^

.PHONY : clean
clean :
	rm -f *.MOD *.mod *.o *.a

%.o : %.f90
	$(MMMFC) -c $(FFLAGS) $< -o $@

%.o : %.f
	$(MMMFC) -c $(FFLAGS) $< -o $@

