#======================================
# MMM - Secondary Makefile for testmmm
#======================================
#
# Lixiang Luo, 9-April-2010
# lixiang.luo@lehigh.edu
# Department of Physics, Lehigh University
#
# A component of the MMM package
#
# You need to change the MMMFC and FFLAGS according to your
# compiler requirements.
#
# Default target is testmmm, which doesn't include debug info.
# To build default target (testmod):
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

lib_mmm   := ../libmmm7_1/libmmm7_1.a
libraries := $(lib_mmm)
f90module := -I../libmmm7_1

testmmm  := testmmm

# All object files are listed here
# Modules are always built first

objects = \
   testmmm.o

MMMFC  ?= gfortran
FFLAGS := -O3 $(f90module)
#FFLAGS := -O3 $(f90module) -Wall -fbounds-check -ffpe-trap=invalid,zero,overflow -finit-real=nan

all : $(testmmm)

.PHONY : debug
debug : FFLAGS := $(subst -O3,-g,$(FFLAGS))
debug : $(testmmm)

testmmm.o : ../libmmm7_1/modmmm7_1.o

$(testmmm) : $(objects) $(libraries)

$(testmmm) : $(objects) $(modules)
	$(MMMFC) -o testmmm $(objects) $(modules) $(libraries)

.PHONY : clean
clean :
	rm -f testmmm testmmm.exe *.MOD *.mod *.o

%.o : %.f
	$(MMMFC) -c $(FFLAGS) $< -o $@

%.o : %.f90
	$(MMMFC) -c $(FFLAGS) $< -o $@

