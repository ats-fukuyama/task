#GNU Makefile to build Library libglf.a 
#
#03/15/10 C. Ludescher-Furth
#
# define system dependent flags, etc.
-include  ../../share/Make.local
include ../../share/Make.flags

ifneq ($(MAKELEVEL),0)
# if OBJ was defined for main Makefile two levels down
  ifeq ("${OBJ}",".")
	override OBJ=../..
  else
	OBJ = ../../$(MACHINE)
  endif
else
  OBJ = ../../$(MACHINE)
endif
MODFLAGS = $(MFLAGS)$(MDIR)
PYDIR    = ../../fpreproc

ifdef DEBUG
	FFLAGS = $(DFFLAGS)
endif

ARC = $(OBJ)/lib/libglf.a
OBJDIR = $(OBJ)/obj/glf

# don't rebuild library if using installed one in $PREFIX 
ifeq ($(MAKELEVEL),0)
	THISLIB=$(ARC)
endif

#=============================================================
# f90 modules
Mnams = glf23_data_mod.mod
Mobjs = glf23_data_mod.o 
MODS = $(foreach file,$(Mobjs),$(ARC)($(file)))

ifeq ($(MODUP),Y)
 MODS0=$(foreach m,$(Mnams),$(shell  echo $(m) | tr 'a-z' 'A-Z'))
 MODULES=$(foreach m,$(MODS0),$(subst .MOD,.$(MODEXT),$(m)))
else
 MODULES = $(foreach m,$(Mnams),$(subst .mod,.$(MODEXT),$(m)))
endif 
 
# fortran archive members
F  = $(subst .f,.o, $(wildcard *.f)) $(subst .F,.o, $(wildcard *.F))
F90 = $(subst .f90,.o, $(wildcard *.f90)) $(subst .F90,.o, $(wildcard *.F90))

FALL = $(F) $(F90)  
MEM = $(foreach file,$(FALL),$(ARC)($(file)))

.PHONY: clean realclean checkdir all libs

libs: FORCE $(ARC)
 
all:  libs 
	@echo done  

$(OLDLIB): timestamp.inf
	@echo "--- DETECTED glf.a source update"
	@echo "--- Re-Making $(ARC)"
	@$(MAKE) libs

FORCE:
	@test -d $(OBJ)/lib || mkdir -p $(OBJ)/lib
	@test -d $(OBJ)/mod || mkdir -p $(OBJ)/mod
	@test -d $(OBJ)/obj/glf || mkdir -p $(OBJ)/obj/glf

 
# compile 
$(ARC)(%.o): %.f
	$(FC) $(FFLAGS) $(MODFLAGS) -I./ $(INCFLAGS) $(OUT)  $(OBJDIR)/$*.o $<
	$(make_ar)

# pre-process/compile fixed fortran code
$(ARC)(%.o): %.F
	$(F_proc)
	$(make_ar)

#  compile free form fortran 90 code
$(ARC)(%.o): %.f90
	$(FC90) $(FFLAGS) $(MODFLAGS) -I./ $(INCFLAGS) $(OUT) $(OBJDIR)/$*.o $<
	$(make_ar)

# pre-process/compile free form fortran 90 code
$(ARC)(%.o): %.F90
	$(F90_proc)
	$(make_ar)


$(ARC): $(MODS) $(MEM)
	@$(RANLIB) $@ > /dev/null

#
#---------------------------------------------------
chktest: 
	@test -d $(OBJ)/test || mkdir -p $(OBJ)/test


exec: chktest $(OBJS) $(OBJ)/test/testglf


clean: 
	@rm -f $(OBJDIR)/*

realclean: clean
	echo ARC = $(ARC)
	echo MDIR = $(MDIR)
	@rm -f $(ARC)
	@if test -d $(MDIR); then \
	  (cd $(MDIR); rm -f $(MODULES)); fi

install:
	@test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	@if test -f $(ARC); then \
	  umask 133; cp $(ARC) $(LIBDIR)/; \
	else test -f $(LIBDIR)/libglf.a || echo ">>> can not find $(ARC)";fi;  
#	@test -d $(MODDIR) || mkdir $(MODDIR)
#	@(cd $(MDIR); umask 133; cp $(MODULES) $(MODDIR)/)

uninstall:
	rm -f $(LIBDIR)/libglf.a
#	(cd $(MODDIR); rm -f $(MODULES))
