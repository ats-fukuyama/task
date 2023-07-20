# Select the correct CXX and F90 depending on the COMPILER variable.
# COMPILER can be a C++ or a Fortran compiler,
# but it might also be a compiler suite, like gcc or intel.
# The rules are:
# 1. If COMPILER is defined and known, we set CXX, F90 accordingly.
# 2. If COMPILER is defined but unknown,
#    we set CXX, F90 to COMPILER only if they are not defined.
# 3. If COMPILER is not defined, we set CXX, F90
#     to a default (GCC) only if they are not defined.

# Logical OR: https://stackoverflow.com/questions/7656425/makefile-ifeq-logical-or
ifneq (, $(filter $(COMPILER), gfortran g++ gcc GCC))
	CXX  = g++
	F90  = gfortran
else ifneq (, $(filter $(COMPILER), ifort icc icpc intel))
	# should be icpc, but FC2K needs icc
	CXX  = icc
	F90  = ifort
# # continue here with similar tests for other compilers, like pgi or nvidia. Examples:
# else ifneq (, $(filter $(COMPILER), nvfortran nvc++ nvidia NVIDIA))
# 	CXX = nvcc
# 	F90 = nvfortran
# else ifneq (, $(filter $(COMPILER), pg90 pg95 pgfortran pgc++ pgi PGI))
# 	CXX = pgc++
# 	F90 = pgfortran

# Finally:
else ifdef COMPILER # unknown COMPILER
	CXX  ?= $(COMPILER)
	F90  ?= $(COMPILER)
else # undefined COMPILER
	CXX  ?= g++
	F90  ?= gfortran
endif

# Now we can make decisions based of CXX and F90...

