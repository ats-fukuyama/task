#===========================
# MMM7.1 - Master Makefile
#===========================
#
# Lixiang Luo, 14-Oct-2010
# lixiang.luo@lehigh.edu
# Department of Physics, Lehigh University
#
# A component of the MMM7.1 package
#
# Before building, change the environment variables MMMFC and FFLAGS
# according to your compiler requirements.
#
# Default target is testmmm, which doesn't include debug info.
# To build default target (testmmm):
#
#   $ make
#
# To clean up all binaries:
#
#   $ make clean
#
# To generate the debug-ready target:
#
#   $ make clean
#   $ make debug
#

# Default to gfortran for maximum compatibility
ifeq ($(origin MMMFC), undefined)
   MMMFC:=gfortran
endif
# Pass down the compiler selection to sublevel makefiles
export MMMFC

# Directory names of libraries
libraries := libmmm7_1

# Directory of executable programs
testmmm   := testmmm

.PHONY : all $(testmmm) $(libraries)
all : $(testmmm)

.PHONY : clean
clean : all
clean : TARGET=clean

.PHONY : debug
debug : all
debug : TARGET=debug

# The rule to make sublevel targets
$(testmmm) $(libraries):
	$(MAKE) --directory=$@ $(TARGET)

# Extra dependencies
$(testmmm) : $(libraries) # The driver depends on the module

