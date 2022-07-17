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

#flag_GSAF = 0
flag_GSAF = 1

#PREPRO = -fpp
PREPRO = -cpp
#PREPRO += -Ddebug # debug option
PREPRO += -D_LAPACK=0 # 0: self-compiled, 1: LAPACK95, 2: LAPACK77

ifeq ($(flag_GSAF),0)
 PREPRO += -DnonGSAF
endif


SRCMI= mi/mat_inv.f90 mi/eqneo.f90
SRCTRC=rktrc/nrtype.f90 rktrc/nrutil.f90 rktrc/recipes.f90 rktrc/sub.f90 rktrc/trace.f90
SRCEQ= equ/eq_mod.f90 $(SRCTRC) equ/eqread.f90
SRCS = txcomm.f90 spln.f90 crosssec.f90 txintf.f90 txlib.f90 txrslt.f90 txmisc.f90 txncls.f90\
       txrppl.f90 txauxs.f90 txntv.f90 cdbm.f90 $(SRCEQ) $(SRCMI) txcalv.f90 txcala.f90 txgraf.f90 \
       txequ.f90 txinit.f90 txfile.f90 txexec.f90 txmenu.f90 txmain.f90 txmmm.f90 coulomb_log.f90
ifeq ($(flag_GSAF),1)
 SRC2D= txg2d.f90
 SRC3D= txg3d.f90
else
 SRC2D= 
 SRC3D= 
endif
SRCNC= nclass/nclass_mod.f nclass/rarray_sum.f \
       nclass/rarray_copy.f nclass/rarray_zero.f nclass/u_erf.f \
       nclass/u_lu_backsub.f nclass/u_lu_decomp.f nclass/write_mod.f
INCNC= txncls.inc nclass/pamx_mi.inc nclass/pamx_ms.inc nclass/pamx_mz.inc
SRCMM= mmm/mmm95.f mmm/weiland14.f mmm/stripx.f mmm/tomsqz.f

ifeq ($(flag_GSAF),1)
 SRCLIB=
 OBJLIB= $(SRCLIB:.f90=.o)
 LIBS2= libtx2.a ../lib/libgrf.a ../lib/libtask.a
 LIBS3= libtx3.a ../lib/libgrf.a ../lib/libtask.a
 MODLIB= -I../lib/mod -I/usr/local/include
else
 SRCLIB= lib/libbnd.f90 lib/libchar.f90 lib/libfio.f90 lib/libitp.f90 lib/libspl.f90 \
         lib/libspl2d.f90 lib/libell.f90 lib/libbes.f90 lib/libfft.f90 lib/libinv.f90 \
	 lib/libspf.f90
 SRCFIXLIB= lib/netlib-bes.f
 OBJLIB= $(SRCLIB:.f90=.o) $(SRCFIXLIB:.f=.o)
 LIBS2= libtx2.a
 LIBS3= libtx3.a
 MODLIB=
endif

OBJS = $(SRCS:.f90=.o)
OBJ2D= $(SRC2D:.f90=.o)
OBJ3D= $(SRC3D:.f90=.o)
OBJNC= $(SRCNC:.f=.o)
OBJMM= $(SRCMM:.f=.o)
PROG2= tx2
PROG3= tx

.SUFFIXES: .f90

.f.o :
	$(FCFIXED) $(PREPRO) $(FOMP) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODLA95) $(MODLIB)
.f90.o :
	$(FCFREE) $(PREPRO) $(FOMP) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODLA95) $(MODLIB)
%.mod: %.f90 %.o
	@:

all: $(PROG2)
#all: $(PROG3)

libtx2.a: $(OBJLIB) $(OBJS) $(OBJ2D) $(OBJNC) $(INCNC) $(OBJMM)
	$(LD) $(LDFLAGS) libtx2.a $(OBJS) $(OBJ2D) $(OBJNC) $(OBJMM) $(OBJLIB)

libtx3.a: $(OBJLIB) $(OBJS) $(OBJ3D) $(OBJNC) $(INCNC) $(OBJMM)
	$(LD) $(LDFLAGS) libtx3.a $(OBJS) $(OBJ3D) $(OBJNC) $(OBJMM) $(OBJLIB)

../lib/libtask.a:
	(cd ../lib; make)

libs:
	(cd ../lib; make)

ifeq ($(flag_GSAF),0)
 FLIBS = $(LIBLA)
 GLIBS = $(LIBLA)
endif

$(PROG2): $(LIBS2)
	$(FCFREE) $(FOMP) $(LIBS2) -o $@ $(FFLAGS_LINK) $(FLIBS)

$(PROG3): $(LIBS3)
	$(FCFREE) $(FOMP) $(LIBS3) -o $@ $(FFLAGS_LINK) $(GLIBS)

check: $(SRCS) $(SRC2D) $(SRC3D)
	ftnchek $(FCKFLAGS) $(SRCS) $(SRC2D) $(SRC3D) | $(PAGER)

clean:
	$(RM) ./*.o ./*.a equ/*.o lib/*.o mi/*.o mmm/*.o nclass/*.o rktrc/*.o \
              ./*~        equ/*~  lib/*~  mi/*~  mmm/*~  nclass/*~  rktrc/*~  \
        $(MOD)/*.mod

veryclean: clean
	$(RM) $(PROG2) $(PROG3)

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
txauxs.o: txauxs.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_core_module.mod \
          $(MOD)/mod_cross_section.mod 
txcala.o: txcala.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_core_module.mod
txcalv.o: txcalv.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/cdbm_mod.mod \
	  $(MOD)/matrix_inversion.mod $(MOD)/aux_system.mod $(MOD)/sauter_mod.mod \
	  $(MOD)/tx_nclass_mod.mod $(MOD)/tx_ripple.mod $(MOD)/tx_core_module.mod \
	  $(MOD)/mod_eqneo.mod $(MOD)/mod_cross_section.mod
txequ.o : txequ.f90  $(MOD)/tx_commons.mod $(MOD)/eqread_mod.mod $(MOD)/tx_interface.mod
txexec.o: txexec.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_coefficients.mod \
	  $(MOD)/tx_variables.mod $(MOD)/tx_graphic.mod $(MOD)/tx_ntv.mod $(MOD)/tx_glob.mod
txfile.o: txfile.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_coefficients.mod \
	  $(MOD)/tx_variables.mod $(MOD)/tx_parameter_control.mod $(MOD)/tx_graphic.mod \
	  $(MOD)/tx_core_module.mod $(MOD)/mod_eqneo.mod $(MOD)/mod_cross_section.mod \
          $(MOD)/tx_glob.mod
txgraf.o: txgraf.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_variables.mod \
	  $(MOD)/tx_ripple.mod $(MOD)/tx_core_module.mod $(MOD)/equ_params.mod
txinit.o: txinit.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_variables.mod \
	  $(MOD)/tx_graphic.mod $(MOD)/sauter_mod.mod $(MOD)/tx_ntv.mod $(MOD)/tx_core_module.mod \
	  $(MOD)/eqread_mod.mod	$(MOD)/mod_eqneo.mod $(MOD)/mod_cross_section.mod \
          $(MOD)/tx_glob.mod
txlib.o : txlib.f90  $(MOD)/tx_commons.mod $(MOD)/mod_spln.mod
txmain.o: txmain.f90 $(MOD)/tx_commons.mod $(MOD)/tx_parameter_control.mod
txmenu.o: txmenu.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_main.mod \
	  $(MOD)/tx_graphic.mod $(MOD)/tx_variables.mod $(MOD)/tx_parameter_control.mod \
	  $(MOD)/tx_ripple.mod $(MOD)/tx_ntv.mod $(MOD)/eqread_mod.mod $(MOD)/mod_cross_section.mod
txmisc.o: txmisc.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod
txmmm.o : txmmm.f90  $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod
txncls.o: txncls.f90 $(MOD)/tx_commons.mod
txntv.o : txntv.f90  $(MOD)/tx_commons.mod $(MOD)/tx_ripple.mod
txrppl.o: txrppl.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod
txrslt.o: txrslt.f90 $(MOD)/tx_commons.mod $(MOD)/tx_interface.mod $(MOD)/tx_core_module.mod
mi/mat_inv.o: mi/mat_inv.f90 $(MOD)/tx_commons.mod
mi/eqneo.o: mi/eqneo.f90 $(MOD)/tx_commons.mod $(MOD)/mod_spln.mod $(MOD)/equ_params.mod
nclass/nclass_mod.o : nclass/nclass_mod.f
	$(FCFIXED) $(FFLAGS_NCLASS) -c nclass/nclass_mod.f -o nclass/nclass_mod.o -Inclass
equ/eqread.o: equ/eqread.f90 $(MOD)/tx_interface.mod $(MOD)/tx_commons.mod $(MOD)/equ_params.mod \
	$(MOD)/mod_spln.mod $(MOD)/tx_core_module.mod
rktrc/nrutil.o: rktrc/nrutil.f90 $(MOD)/nrtype.mod
rktrc/trace.o: rktrc/trace.f90 $(MOD)/tx_commons.mod $(MOD)/equ_params.mod $(MOD)/subs.mod \
               $(MOD)/mod_num_recipe.mod
