# 
# gmake makefile for MMM95 model module 
#
# NOTE Must use gmake not make.
# The user must set the environmental variable "CPU" to be
# one of the following:
#
# FUJITSU  for the Lahey Fujitsu F90/F95 compiler
# IBM      for the standard compiler on IBM/AIX workstations
# SOL      for the standard compiler on Sun Solaris workstations
# SGI      for the standard compiler on SGI IRIX workstations
# DEC      for the standard compiler on DEC alpha workstations
# HP       for the standard compiler on Hewlett Packard workstations
#
#       For example, in the Bourne shell on an SGI workstation, type:
#
#  export CPU=SGI
#
#       If you are using the C-shell on an SGI workstation, type
#
#  setenv CPU SGI
#
# Use `gmake clean' to remove all binaries.
#
# The MMM95 transport model has been used ONLY in DOUBLE PRECISION
# with the precision set by flags in the makefile
#
# Glenn Bateman and Arnold Kritz, Lehigh University
#
#
# environment variable CPU is used to set platform-dependent make variables.
#
ifeq ($(CPU),)
   # CPU is undefined; default to SGI.
   WARNcpu = 2
   CPU = SGI
endif
ifeq ($(filter $(CPU), FUJITSU, IBM SOL SGI DEC HP),)
   ERRcpu = 2
endif

#
#
# Compilation/linking flags
#
# Compilation/linking flags for LINUX
#
ifeq ($(CPU),FUJITSU)
   FC=lf95
   FFLAGS=--dbl -g --tpp --chk a,e,s,u
   COMPILE = $(FC) $(FFLAGS) -c
   LINK    = $(FC)
endif
#
# For IBM/AIX compilers, the compilation flags are (in order):
#   - Optimization (level 1 out of 3)
#   - Default real-variable size is 8 bytes
#   - Default integer-variable size is 4 bytes
#   - No limit on the executable's memory usage
#     (works with the -O flag to give greater optimization)
#   - Input files are written in fixed-format fortran, lines 72 chars long
#   - Compile only (linking is done in a seperate step)
#
ifeq ($(CPU),IBM)
   FFLAGS  = -O -qrealsize=8 -qintsize=4 -qmaxmem=-1
   FFLAGS2 = -qfixed=72 -c
   COMPILE = $(FC) $(FFLAGS) $(FFLAGS2)
   LINK    = $(FC) $(FFLAGS)
endif
#
# For Solaris compilers, these compilation flags are (in order):
#
ifeq ($(CPU),SOL)
   FFLAGS  = -O4 -Nx600 -Nn2000 -r8 -i4
   FFLAGS2 = -c
   COMPILE = $(FC) $(FFLAGS) $(FFLAGS2)
   LINK    = $(FC) $(FFLAGS)
endif
# For SGI compilers, these compilation flags are (in order):
#
ifeq ($(CPU),SGI)
   FFLAGS  = -default64
   FFLAGS2 = -c
   FC=f90
   COMPILE = $(FC) $(FFLAGS) $(FFLAGS2)
   LINK    = $(FC) $(FFLAGS)
endif
#
# For DEC-alpha workstations (dmc)
#
ifeq ($(CPU),DEC)
   FFLAGS  = -O -r8 -i4
   FFLAGS2 = -c
   COMPILE = $(FC) $(FFLAGS) $(FFLAGS2)
   LINK    = $(FC) $(FFLAGS)
endif
#
# For HP workstations (dmc)
#
ifeq ($(CPU),HP)
   FFLAGS = -O3 +autodblpad
   FFLAGS2 = -c
   FC=fort77
   COMPILE = $(FC) $(FFLAGS) $(FFLAGS2)
   LINK    = $(FC) $(FFLAGS)
endif

#
# Targets and Rules
#

xmmm95:	test_mmm95.o mmm95.o stripx.o weiland14.o tomsqz.o
	$(LINK) -o xmmm95 test_mmm95.o mmm95.o stripx.o \
                weiland14.o tomsqz.o

stripx.o:	stripx.f
	$(COMPILE) -o stripx.o stripx.f

tomsqz.o:	tomsqz.f
	$(COMPILE) -o tomsqz.o tomsqz.f

mmm95.o:	mmm95.f
	$(COMPILE) -o mmm95.o mmm95.f

weiland14.o:	weiland14.f
	$(COMPILE) -o weiland14.o weiland14.f

test_mmm95.o:	test_mmm95.f
	$(COMPILE) -o test_mmm95.o test_mmm95.f

cputst:
  ifeq ($(ERRcpu),2)
	@echo CPU=$(CPU) is not yet implemented in mmm95 makefile.
	@exit 2
  endif
  ifeq ($(WARNcpu),2)
	@echo CPU defaulted, CPU=IBM assumed.
  endif

clean:
	@rm *.o xmmm95

