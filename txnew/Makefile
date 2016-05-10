### $Id: Makefile,v 1.66 2011/06/13 07:53:20 honda Exp $ ###
include ../make.header

#FOMP = -qopenmp
FOMP = 
FFLAGS_MAIN = $(OFLAGS)
#FFLAGS_MAIN = $(DFLAGS)
FFLAGS_LINK = $(FFLAGS_MAIN)
FFLAGS = $(FFLAGS_MAIN)
#FFLAGS = $(FFLAGS_MAIN) $(FOMP) #-fpe0
FFLAGS_NCLASS = $(OFLAGS)
#FFLAGS_NCLASS = $(DFLAGS)

#PREPRO = -fpp -Ddebug # debug option
#PREPRO = -fpp -Dlaself # when using self-compiled LAPACK
PREPRO = -fpp

SRCMI= mi/mat_inv.f90
SRCEQ= equ/eq_mod.f90 equ/eqread.f90
SRCS = txcomm.f90 spln.f90 txintf.f90 txlib.f90 txmisc.f90 txncls.f90 txrslt.f90 txrppl.f90 \
       txauxs.f90 txntv.f90 $(SRCMI) cdbm.f90 txcalv.f90 txcala.f90 $(SRCEQ) txgraf.f90 txequ.f90 \
       txinit.f90 txfile.f90 txexec.f90 txmenu.f90 txmain.f90 txmmm.f90 coulomb_log.f90
SRC2D= txg2d.f90
SRC3D= txg3d.f90
SRCNC= nclass/nclass_mod.f nclass/rarray_sum.f \
       nclass/rarray_copy.f nclass/rarray_zero.f nclass/u_erf.f \
       nclass/u_lu_backsub.f nclass/u_lu_decomp.f nclass/write_mod.f
INCNC= txncls.inc nclass/pamx_mi.inc nclass/pamx_ms.inc nclass/pamx_mz.inc
SRCMM= mmm/mmm95.f mmm/weiland14.f mmm/stripx.f mmm/tomsqz.f
OBJS = $(SRCS:.f90=.o)
OBJ2D= $(SRC2D:.f90=.o)
OBJ3D= $(SRC3D:.f90=.o)
OBJNC= $(SRCNC:.f=.o)
OBJMM= $(SRCMM:.f=.o)
PROG2= tx2
PROG3= tx
LIBS2= tx2lib.a ../lib/tasklib.a
LIBS3= tx3lib.a ../lib/tasklib.a
MODLIB= -I../lib/mod

.SUFFIXES: .f90

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODLA95) $(MODLIB)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODLA95) $(MODLIB)
%.mod: %.f90 %.o
	@:

#all: $(PROG2)
all: $(PROG3)

tx2lib.a: $(OBJS) $(OBJ2D) $(OBJNC) $(INCNC) $(OBJMM)
	$(LD) $(LDFLAGS) tx2lib.a $(OBJS) $(OBJ2D) $(OBJNC) $(OBJMM)

tx3lib.a: $(OBJS) $(OBJ3D) $(OBJNC) $(INCNC) $(OBJMM)
	$(LD) $(LDFLAGS) tx3lib.a $(OBJS) $(OBJ3D) $(OBJNC) $(OBJMM)

../lib/tasklib.a:
	(cd ../lib; make)

libs:
	(cd ../lib; make)

$(PROG2): $(LIBS2)
	$(FCFREE) $(FOMP) $(LIBS2) -o $@ $(FFLAGS_LINK) $(FLIBS)

$(PROG3): $(LIBS3)
	$(FCFREE) $(FOMP) $(LIBS3) -o $@ $(FFLAGS_LINK) $(GLIBS)

check: $(SRCS) $(SRC2D) $(SRC3D)
	ftnchek $(FCKFLAGS) $(SRCS) $(SRC2D) $(SRC3D) | $(PAGER)

clean:
	rm -f *.o ./*~ ./#* nclass/*.o mmm/*.o mi/*.o equ/*.o *.a $(MOD)/*.mod \
        *.optrpt nclass/*.optrpt mmm/*.optrpt mi/*.optrpt equ/*.optrpt

veryclean: clean
	rm -f $(PROG2) $(PROG3)

distclean: clean
	cd ../lib && make clean

new :
	mkdir ../txnew
	cp -f Makefile ../txnew
	cp -f HISTORY ../txnew
	cp -f *.f90 ../txnew
	cp -f *.inc ../txnew
	-mkdir ../txnew/doc
	-cp -f doc/* ../txnew/doc
	-mkdir ../txnew/in
	-cp -f in/* ../txnew/in

#
# Module dependencies
#
txauxs.o: txauxs.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_core_module.mod
txcala.o: txcala.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_core_module.mod
	$(FCFREE) $(FFLAGS) $(FOMP) -c $< -o $@ $(MODDIR) $(MODLA95) $(MODLIB)
txcalv.o: txcalv.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/cdbm_mod.mod \
	  $(MOD)/matrix_inversion.mod $(MOD)/aux_system.mod $(MOD)/sauter_mod.mod \
	  $(MOD)/tx_nclass_mod.mod $(MOD)/tx_ripple.mod $(MOD)/tx_core_module.mod
	$(FCFREE) $(FFLAGS) $(FOMP) -c $< -o $@ $(MODDIR) $(MODLA95) $(MODLIB)
txequ.o : txequ.f90  $(MOD)/tx_commons.mod $(MOD)/eqread_mod.mod $(MOD)/tx_interface.mod
txexec.o: txexec.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_coefficients.mod \
	  $(MOD)/tx_variables.mod $(MOD)/tx_graphic.mod $(MOD)/tx_ntv.mod
	$(FCFREE) $(FFLAGS) $(PREPRO) $(FOMP) -c $< -o $@ $(MODDIR) $(MODLA95) $(MODLIB)
txfile.o: txfile.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_coefficients.mod \
	  $(MOD)/tx_variables.mod $(MOD)/tx_parameter_control.mod $(MOD)/tx_graphic.mod \
	  $(MOD)/tx_core_module.mod
txgraf.o: txgraf.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_variables.mod \
	  $(MOD)/tx_ripple.mod $(MOD)/tx_core_module.mod $(MOD)/equ_params.mod
txinit.o: txinit.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_variables.mod \
	  $(MOD)/tx_graphic.mod $(MOD)/sauter_mod.mod $(MOD)/tx_ntv.mod $(MOD)/tx_core_module.mod \
	  $(MOD)/eqread_mod.mod
	$(FCFREE) $(FFLAGS) $(PREPRO) -c $< -o $@ $(MODDIR) $(MODLA95) $(MODLIB)
txlib.o : txlib.f90  $(MOD)/tx_commons.mod $(MOD)/mod_spln.mod
	$(FCFREE) $(FFLAGS) $(PREPRO) -c $< -o $@ $(MODDIR) $(MODLA95) $(MODLIB)
txmain.o: txmain.f90 $(MOD)/tx_commons.mod $(MOD)/tx_parameter_control.mod
txmenu.o: txmenu.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_main.mod \
	  $(MOD)/tx_graphic.mod $(MOD)/tx_variables.mod $(MOD)/tx_parameter_control.mod \
	  $(MOD)/tx_ripple.mod $(MOD)/tx_ntv.mod $(MOD)/eqread_mod.mod
txmisc.o: txmisc.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod
txmmm.o : txmmm.f90  $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod
txncls.o: txncls.f90 $(MOD)/tx_commons.mod
txntv.o : txntv.f90  $(MOD)/tx_commons.mod $(MOD)/tx_ripple.mod
txrppl.o: txrppl.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod
txrslt.o: txrslt.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_core_module.mod
mat_inv.o: mi/mat_inv.f90 $(MOD)/tx_commons.mod
nclass/nclass_mod.o : nclass/nclass_mod.f
	$(FCFIXED) $(FFLAGS_NCLASS) -c nclass/nclass_mod.f -o nclass/nclass_mod.o -Inclass
eqread.o: eqread.f90 $(MOD)/tx_interface.mod $(MOD)/tx_commons.mod $(MOD)/equ_params.mod \
	$(MOD)/mod_spln.mod $(MOD)/tx_core_module.mod
