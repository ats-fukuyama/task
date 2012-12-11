# Print FLAGS for debugging
# 01/27/00 CAL extracted from makeflags.mk
#              Don't want flagtest as a 1st rule in a Makefile
#
# to run:
# > make -f $CODESYSDIR/source/misc/flagtest.mk
#
#  1Dec2008  jim.conboy@jet.uk
#            Add target libroot - displays $(LIBROOT)
#----------------------------------------------------------------
#

ifdef CODESYSDIR
  include $(CODESYSDIR)/source/misc/makeflags.mk
else
  ifneq ("$(wildcard $(NTCC_ROOT)/etc/Make.flags)","")
    include $(NTCC_ROOT)/etc/Make.flags
  else
    ifneq ("$(wildcard /usr/ntcc/etc/Make.flags)","")
      include /usr/ntcc/etc/Make.flags
    else
      -include makeflags.mk
      -include Make.local
      -include Make.flags
    endif
  endif
endif

ifneq ("$(wildcard $(MODDIR)/*)","")
	EXMODFLAGS := $(MFLAGS)$(MODDIR)
endif

example:
	@echo " "
	@echo "   fortran example "
	@echo "---------------------"
	@echo $(FC) $(FFLAGS) $(OUT) f.o f.f
	@echo $(FC) $(LDFLAGS) -o prog f.o $(XFLIB)
	@echo " "
	@echo "   if you need ntcc libraries"
	@echo "--------------------------------"
	@echo "$(FC) $(FFLAGS) \\"
	@echo "      $(EXMODFLAGS) -I./ $(INCFLAGS) $(OUT) f.o f.f "
	@echo $(FC) $(LDFLAGS) -o prog f.o  -L$(LIBDIR) -lportlib  $(XFLIB)
	@echo " "
	@echo "  c++ example "
	@echo "----------------"
	@echo $(CXX) $(CXXFLAGS) -o c.o c.cpp
	@echo $(CXX) $(LDFLAGS) -o prog c.o 
	@echo " "
	@echo "   mixing c++ / fortran "
	@echo "--------------------------"
	@echo $(CXX) $(LDFLAGS) -o prog c.o f.o $(FORTLIBS)
	@echo $(FC) $(LDFLAGS) -o prog f.o c.o $(CLIBS)
	@echo " "
	@echo "   other libraries  "
	@echo "--------------------"
	@echo $(XLIBS) $(PPLLIB) $(MDSLIB)
	@echo " "

libroot:
   ifdef _64
	@echo " LIBROOT (64bit) : $(LIBROOT)"
   else
	@echo " LIBROOT         : $(LIBROOT)"
   endif

all: flagtest example

flagtest:
	@echo "------------------------------"
	@echo "flagtest results:"
  ifdef MKFUJITSU
	@echo MKFUJITSU set
  endif
  ifdef MKNAGFLAG
	@echo MKNAGFLAG set
  endif
  ifdef FORTRAN_TYPE
	@echo FORTRAN_TYPE = $(FORTRAN_TYPE)
  endif
  ifdef FORTRAN_VARIANT
	@echo FORTRAN_VARIANT = $(FORTRAN_VARIANT)
  endif
	@echo "...general:"
	@echo "  LIBROOT = $(LIBROOT)"
	@echo "  PREFIX  = $(PREFIX)"
	@echo "  LIBDIR  = $(LIBDIR)"
	@echo "  BINDIR  = $(BINDIR)"
	@echo "  ETCDIR  = $(ETCDIR)"
	@echo "  INCLDIR = $(INCLDIR)"
	@echo "  MODDIR  = $(MODDIR)"
	@echo "  MANDIR  = $(MANDIR)"
	@echo "  LLOC    = $(LLOC)"

	@echo "...fortran:"
	@echo "  FC        = $(FC)"
	@echo "  FC90      = $(FC90)"
	@echo "  FFLAGS    = $(FFLAGS)"
	@echo "  X_FFLAGS  = $(X_FFLAGS)"
	@echo "  F90FLAGS  = $(F90FLAGS)"
	@echo "  DFFLAGS   = $(DFFLAGS)"
	@echo "  DF90FLAGS = $(DF90FLAGS)"
	@echo "  STATIC    = $(STATIC)"
	@echo "  DYNAMIC   = $(DYNAMIC)"
	@echo "  DEFS      = $(DEFS)"
	@echo "  OUT       = $(OUT)"

	@echo "...C/C++:"
	@echo "  CC        = $(CC)"
	@echo "  CFLAGS    = $(CFLAGS)"
	@echo "  DCFLAGS   = $(DCFLAGS)"
	@echo "  OPT       = $(OPT)"
	@echo "  CPP       = $(CPP)"
	@echo "  CPPFLAGS  = $(CPPFLAGS)"
	@echo "  CXX       = $(CXX)"
	@echo "  CXXFLAGS  = $(CXXFLAGS)"
	@echo "  CDEFS     = $(CDEFS)"

	@echo "...MPI:"
	@echo "  MPI_FC       = $(MPI_FC)"
	@echo "  MPI_FFLAGS   = $(MPI_FFLAGS)"
	@echo "  MPI_DFFLAGS  = $(MPI_DFFLAGS)"
	@echo "  MPI_DEFS     = $(MPI_DEFS)"
	@echo "  MPI_CC       = $(MPI_CC)"
	@echo "  MPI_CFLAGS   = $(MPI_CFLAGS)"
	@echo "  MPI_DCFLAGS  = $(MPI_DCFLAGS)"
	@echo "  MPI_CXX      = $(MPI_CXX)"
	@echo "  MPI_CXXFLAGS = $(MPI_CXXFLAGS)"
	@echo "  MPI_CDEFS    = $(MPI_CDEFS)"
	@echo "  MPI_LDFLAGS  = $(MPI_LDFLAGS)"

	@echo "...Include/Module:"
	@echo "  IOPT    = $(IOPT)"
	@echo "  MFLAGS  = $(MFLAGS)"
	@echo "  MODEXT  = $(MODEXT)"
	@echo "  MFFLAGS = $(MFFLAGS)"
	@echo "  MODFLAGS= $(MODFLAGS)"

	@echo "...Linking, Libraries:"
	@echo "  FCEXE     = $(FCEXE)"
	@echo "  LDFLAGS   = $(LDFLAGS)"
	@echo "  LD        = $(LD)"
	@echo "  X_LDFLAGS = $(X_LDFLAGS)"
	@echo "  MVU       = $(MVU)"
	@echo "  NETCDF_DIR= $(NETCDF_DIR)"
	@echo "  NETCDF    = $(NETCDF)"
	@echo "  HDF       = $(HDF)"
	@echo "  BLAS      = $(BLAS)"
	@echo "  LAPACK    = $(LAPACK)"
	@echo "  USRLIB    = $(USRLIB)"
	@echo "  XLIBS     = $(XLIBS)"
	@echo "  XFLIB     = $(XFLIB)"
	@echo "  MDSLIB    = $(MDSLIB)"
	@echo "  FORTLIBS  = $(FORTLIBS)"
	@echo "  CLIBS     = $(CLIBS)"
	@echo "  EDITLIBS  = $(EDITLIBS)"

	@echo "...Miscellaneous:"
	@echo "  PERL      = $(PERL)"
	@echo "  INSTALL   = $(INSTALL)"
	@echo "  RANLIB    = $(RANLIB)"
	@echo "  OS        = $(OS)"
	@echo "  MACHINE   = $(MACHINE)"
	@echo "  TSTLNK    = $(TSTLNK)"
