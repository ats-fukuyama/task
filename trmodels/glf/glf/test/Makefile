#GNU Makefile to build testglf 
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
 
OFILES = c9date.o zgeev.o blas_zgeev.o stripx.o
FALL = $(F) $(F90)  
OBJS = $(foreach file, $(OFILES),$(OBJDIR)/$(file))

LDLIBS = -L$(OBJ)/lib -lglf

.PHONY: clean realclean checkdir all exec

 
all:  exec 
	@echo done  


$(OBJDIR)/%.o: %.f
	$(FC) $(FFLAGS) $(MODFLAGS) -I./ $(INCFLAGS) $(OUT) $@ $<

#
#---------------------------------------------------
chktest: 
	@test -d $(OBJ)/test || mkdir -p $(OBJ)/test


exec: chktest $(OBJS) $(OBJ)/test/testglf $(OBJ)/test/xtverb


$(OBJ)/test/testglf: $(OBJDIR)/testglf.o $(OBJS)  $(THISLIB) 
	$(FC) $(LDFLAGS) -o $@  $<  $(OBJS) $(LDLIBS) $(LUSRLIB) $(XFLIB) $(XLIBS)
	cp $(OBJ)/test/testglf ../test/

$(OBJ)/test/xtverb: $(OBJDIR)/xtverb.o   
	$(FC) $(LDFLAGS) -o $@  $< 
	cp $(OBJ)/test/xtverb ../test/

$(THISLIB): $(ARC)
	(cd ../src; $(MAKE) libs)

clean: 
	@rm -f $(OBJDIR)/*
	rm -f ../test/testglf

realclean: clean
	rm -f $(OBJ)/test/testglf

