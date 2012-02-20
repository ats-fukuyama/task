#===========================
# MMM - Master Makefile
#===========================
#
# Lixiang Luo, 14-Oct-2010
# lixiang.luo@lehigh.edu
# Department of Physics, Lehigh University
#
# A component of the MMM package
#
# You need to set the environment variable MMMFC to your
# compiler. If not set, it defaults to gfortran.
#
# Default target is testmmm, which doesn't include debug info.
# To buile default target (testmmm):
#
#   $ make
#
# To clean up all binaries:
#
#   $ make clean
#
# To generate the debugging-ready target:
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

# Executable programs
testmmm   := testmmm

.PHONY : all $(testmmm) $(libraries)
all : $(testmmm)

.PHONY : clean
clean : all
clean : TARGET=clean

.PHONY : debug
debug : all
debug : TARGET=debug

$(testmmm) $(libraries):
	$(MAKE) --directory=$@ $(TARGET)

$(testmmm) : $(libraries) # The driver depends on the module

