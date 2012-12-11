# GNU Make file to execute module's Makefiles
# 10/12/99 C. Ludescher
#
# NOTE:
#======
# previously installed libraries (in $PREFIX/lib) are not re-built
# to re-build them do
# (cd <libname>; make all)
#
-include share/Make.local
-include share/Make.flags

# Look for required libraries
NEEDANY=$(wildcard */*_exe.link)
ifneq "$(strip $(NEEDANY))" ""
   NEEDALL := $(sort $(subst .a,, $(filter %.a, $(shell cat */*_exe.link))))
endif

# Look for libraries with public f90 modules
MMODS := $(sort $(notdir $(subst _lib.public,, $(wildcard */*_lib.public))))

#ifneq ("$(wildcard xplasma2/Makefile)","")
	MORDER :=  random ezcdf comput pspline adpak sigsub geqdsk_mds xplasma2 old_xplasma xplasma_debug plasma_state_kernel ps_xplasma2 plasma_state xstraln rfxqlo adas310_fortran_driver nubeam xdatmgr nubeam_comp_sub cppsub key_access
	MMIS := $(filter-out $(MMODS), $(MORDER))
	MEXTRA := $(filter-out $(MORDER), $(MMODS))
	MMODS := $(filter-out $(MMIS), $(MORDER)) $(MEXTRA)
#endif

# Look for Makefiles in sub-directories
MFILES = $(filter-out %~ ,$(wildcard */Makefile))
MDIRS  = $(subst /Makefile,, $(MFILES))
# Fix nun-standard names
MDIRSX := $(subst pest3.2,pest3_2, $(MDIRS))
# remove supplied libs from NEEDALL 
NEED := $(filter-out $(MDIRSX), $(NEEDALL))
# remove libs with public f90 modules
MSUBS := $(filter-out $(MMODS), $(filter-out sglib, $(MDIRS)))

# libs needed for test programs only are not installed
ifneq ("$(wildcard sglib/Make_install)","")
   TSTONLY= sgdummy cmdummy globdummy sddummy
else
   TSTONLY= sglib sgdummy cmdummy globdummy sddummy
endif
MEXEC0 = $(MDIRS)
MEXEC1 = $(filter-out $(TSTONLY), $(MEXEC0))

MEXEC    = $(sort $(MEXEC1) )
MINSTALL = $(sort $(MEXEC1) )


# sglib does not follow naming convention
SGLIB := $(findstring sglib, $(MDIRS))

# fgtok is also an exeption
FTOKEN := $(findstring ftoken, $(MDIRS))
ifneq "$(strip $(FTOKEN))" ""
   FGTOK = cpcheck fgtok idecl r8real
   MEXEC := $(filter-out $(FGTOK), $(MEXEC))
endif

# toric xfpprf exections
TORIC := $(findstring toric, $(MDIRS))
ifneq "$(strip $(TORIC))" ""
	TORICX := xfpprf toric4_main
endif

PGRMS := rplot tr_start xtctek util client $(FGTOK)
MDIRSX := $(filter-out $(PGRMS), $(MDIRSX))

# NEWLIB = current location
ifeq ("$(OBJ)","../$(MACHINE)")
     NEWLIB = $(MACHINE)/lib
     SUB = $(OBJ)/lib
else
     NEWLIB = $(OBJ)/lib
ifeq ("$(OBJ)",".")
     SUB = ../lib
else
     SUB = $(OBJ)/lib
endif
endif

.PHONY: clean realclean install all exec chklibs chklnklibs

all:    prep libs chklnklibs exec

prep:
ifneq ("$(wildcard prefpp/makefile.inf)","")
	(cd prefpp; $(MAKE) -f makefile.inf all)
endif
ifneq ("$(wildcard prefpp/Makefile)","")
	(cd prefpp; $(MAKE) all)
endif
# make libraries with public f90 modules first
# do not re-build existing libraries in $LIBDIR

libs:
# compile sglib
ifneq "$(SGLIB)" ""
	@echo --- CHECKING libsg.a
	@if test  -f $(LIBDIR)/libsg.a; then \
	  (cd sglib; $(MAKE) $(LIBDIR)/libsg.a OLDLIB=$(LIBDIR)/libsg.a);\
	  if test ! -f $(NEWLIB)/libsg.a; then \
	      echo "--- USING    $(LIBDIR)/libsg.a"; fi;\
          else echo "--- MAKING   sglib"; (cd sglib; $(MAKE));fi;
endif

ifneq "$(MMODS)" ""
	@echo Making f90 Modules
	@for m in $(MMODS); do \
          echo --- CHECKING lib$$m.a;\
	  if test  -f $(LIBDIR)/lib$$m.a; then \
	    (cd $$m; $(MAKE) $(LIBDIR)/lib$$m.a OLDLIB=$(LIBDIR)/lib$$m.a);\
	    if test ! -f $(NEWLIB)/lib$$m.a; then \
	      echo "--- USING    $(LIBDIR)/lib$$m.a"; fi;\
          else echo "--- MAKING   $$m"; (cd $$m; $(MAKE));fi; done
endif

# compile libraries, unless already installed
ifneq "$(MSUBS)" ""
	@echo Making Libraries 
	@for m in $(MSUBS); do \
          echo --- CHECKING lib$$m.a;\
	  if test  -f $(LIBDIR)/lib$$m.a; then \
	    (cd $$m; $(MAKE) $(LIBDIR)/lib$$m.a OLDLIB=$(LIBDIR)/lib$$m.a);\
	    if test ! -f $(NEWLIB)/lib$$m.a; then \
	      echo "--- USING    $(LIBDIR)/lib$$m.a"; fi;\
          else echo "--- MAKING   $$m"; (cd $$m; $(MAKE));fi; done
endif
	@echo Done with Libraries

chklnklibs:
# Check for libraries required to link test program(s)
ifneq "$(NEED)" ""
	@echo "Checking for Required Libraries"
	@for m in $(NEED); do \
	  if test ! -f $$m/Makefile; then \
	     if test  -f $(LIBDIR)/lib$$m.a; then \
	       echo "--- USING        $(LIBDIR)/lib$$m.a"; \
	     else if test  -f $(LIBROOT)/lib/lib$$m.a; then \
	       echo "--- USING        $(LIBROOT)/lib/lib$$m.a"; \
	     else if test  -f $(FLIBROOT)/lib/lib$$m.a; then \
	       echo "--- USING        $(FLIBROOT)/lib/lib$$m.a"; \
	     else if test  -f $(NETCDF_DIR)/lib/lib$$m.a; then \
	       echo "--- USING        $(NETCDF_DIR)/lib/lib$$m.a"; \
             else echo "--- NEED         $$m lib"; \
	       echo "    if it's installed define PREFIX"; fi; fi; fi; fi; fi; done
endif

# Make test/executable programs
# The FGTOK group is an exeption.
# These are utilities all dependent on libftoken.a

exec:
	@echo Making Test Programs
ifneq "$(MEXEC)" ""
	@for m in $(MEXEC); do \
	 l=$$m; \
	 if  test $$m = sglib; then l=sg; fi; \
	 if  test $$m = pest3.2; then l=pest3_2; fi; \
	 if  test $$m = tr_start; then l=mds_sub; fi; \
	 if  test $$m = nubeam_comp_exec; then l=nubeam_comp_sub; fi; \
	 if  test $$m = plasma_state_test; then l=plasma_state; fi; \
	 if  test $$m = cstate; then l=cstate_lib; fi; \
	 if  test $$m = update_state; then l=plasma_state; fi; \
	 if  test $$m = geq2ps; then l=cstate_lib; fi; \
	 if  test $$m = ps2geq; then l=cstate_lib; fi; \
	 if  test $$m = ps2xplasma; then l=plasma_state; fi; \
	 if  test $$m = jsocdf2ps; then l=plasma_state; fi; \
	 if  test $$m = ps_nbline_driver; then l=ps_nbline; fi; \
	 if  test $$m = ccps_test; then l=ccps_lib; fi; \
	 if test -f $(NEWLIB)/lib$$l.a; then \
	 echo --- make $$m programs; \
	 (cd $$m; $(MAKE) exec THISLIB=$(SUB)/lib$$l.a); fi; \
	 if  test $$m = xplasma; then (cd $$m; $(MAKE) exec ); fi; \
	 if  test $$m = nubeam_driver; then (cd $$m; $(MAKE) exec ); fi; \
	 if  test $$m = client; then (cd $$m; $(MAKE) exec ); fi; \
	 if  test $$m = rplot; then (cd $$m; $(MAKE) exec ); fi; \
	 if  test $$m = transpin; then (cd $$m; $(MAKE) exec ); fi; \
	 if  test $$m = util; then (cd $$m; $(MAKE) exec ); fi; done
endif

ifneq "$(TORICX)" ""
	@for m in $(TORICX); do \
	 echo --- make $$m programs; \
	(cd $$m; $(MAKE) exec ); done	
endif

ifneq "$(FGTOK)" ""
	@for m in $(FGTOK); do \
	 echo --- make $$m programs; \
	 (cd $$m; $(MAKE) exec THISLIB=$(SUB)/libftoken.a); done	
endif
	@echo Done with all Test Programs

# Install libraries and utilities
install:
ifneq "$(MINSTALL)" ""
ifdef PREFIX
	@if test ! -d $(PREFIX); \
	then echo "!!! PREFIX = $(PREFIX) does not exist !!!";\
	else \
	   for m in $(MINSTALL); do \
	   (cd $$m; $(MAKE) install); done; fi
else
	@echo "!!! You must define PREFIX !!!"
endif
endif

uninstall:
ifneq "$(MINSTALL)" ""
ifdef PREFIX
	@for m in $(MINSTALL); do \
	 (cd $$m; $(MAKE) uninstall); done
endif
endif

clean:
	@for m in $(MDIRS); do \
	 (cd $$m; $(MAKE) clean); done

realclean: clean
	@rm -f *.tar
	@rm -f *.zip
	@for m in $(MDIRS); do \
	  (cd $$m; $(MAKE) realclean); done

checklibs: chklnklibs 
	@for m in $(MDIRSX); do \
	  if test -f $(NEWLIB)/lib$$m.a; then \
	     echo "--- USING       $(NEWLIB)/lib$$m.a"; \
	  else if test -f $(LIBDIR)/lib$$m.a; then \
	     $(test_timestamp) \
	  else echo "--- WILL Make   $$m.a"; fi;fi; done   
	@echo done with libs


ifeq ($(OS),ALPHA)
  define test_timestamp
    if /bin/ksh -c "test $$m/timestamp.inf -nt $(LIBDIR)/lib$$m.a"; then\
    echo --- WILL Re-Make $(NEWLIB)/lib$$m.a;\
    else echo "--- USING        $(LIBDIR)/lib$$m.a"; fi
  endef
else 
ifeq ($(OS),SGI)
  define test_timestamp
    if /bin/ksh -c "test $$m/timestamp.inf -nt $(LIBDIR)/lib$$m.a"; then\
    echo --- WILL Re-Make $(NEWLIB)/lib$$m.a;\
    else echo "--- USING        $(LIBDIR)/lib$$m.a"; fi
  endef
else
  define test_timestamp
    if /usr/bin/test $$m/timestamp.inf -nt $(LIBDIR)/lib$$m.a; then \
    echo --- WILL Re-Make $(NEWLIB)/lib$$m.a;\
    else echo "--- USING        $(LIBDIR)/lib$$m.a"; fi
  endef
endif
endif

show_makeflags:
	cd share; $(MAKE) -f flagtest.mk flagtest
	@echo " "
	@echo " "
	@echo "----------------------------------------- "
	@echo " "
	@echo " "
	@echo "fortran preprocessing"
	@echo "=============="
	@echo "$(PREFPP) <FILE> $(OBJDIR)/<FILE>.f  $(INCFLAGS) $(DEFS) $(DPY)"
	@echo " "
	@echo "fortran compile"
	@echo "==============="
	@echo "$(FC) $(FFLAGS) $(MODFLAGS) -I./ $(INCFLAGS) $(FIXED) $(OUT) $(OBJDIR)/<FILE>.o $(OBJDIR)/<FILE>.f "
	@echo " "
	@echo " "
ifeq ($(MACHINE),CRAY)
	@echo "c processing"
	@echo "============"
	@echo "cat <FILE>.c | $(PERL) -ne \"
	@echo "'s/F77NAME(\([a-zA-Z0-9_ ]+)\)*/F77NAME\U$$1)/g; print;'"
	@echo " > $(OBJDIR)/<FILE>"
	@echo "$(CC) $(CFLAGS) -I./ $(INCFLAGS) $(CDEFS) -o $(OBJDIR)/<FILE>.o $(OBJDIR)/<FILE>"
	@echo "$(CXX) $(CXXFLAGS) -I$(OBJ)/obj $(INCFLAGS) $(CDEFS) -o $(OBJDIR)/<FILE>.o $(OBJDIR)/<FILE>"
else
	@echo "c compile"
	@echo "========="
	@echo "$(CC) $(CFLAGS) -I./ $(CINCL) $(CDEFS) -o $(OBJDIR)/<FILE>.o <FILE>.c"
	@echo "$(CXX) $(CXXFLAGS) -I./ $(CINCL) $(CDEFS) -o $(OBJDIR)/<FILE>.o <FILE>.cpp"
endif
