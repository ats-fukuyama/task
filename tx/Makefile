### $Id$ ###
include ../make.header

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)

SRCS = txintf.f90 txcomm.f90 txlib.f90 txncls.f90 txrslt.f90 txcalv.f90 \
       txcala.f90 txinit.f90 txfile.f90 txgraf.f90 txexec.f90 txmenu.f90 \
       txmain.f90 txntv.f90 txrppl.f90
SRC2D= txg2d.f90
SRC3D= txg3d.f90
SRCNC= nclass/nclass_mod.f nclass/rarray_sum.f \
       nclass/rarray_copy.f nclass/rarray_zero.f nclass/u_erf.f \
       nclass/u_lu_backsub.f nclass/u_lu_decomp.f nclass/write_mod.f
INCNC= txncls.inc nclass/pamx_mi.inc nclass/pamx_ms.inc nclass/pamx_mz.inc
OBJS = $(SRCS:.f90=.o)
OBJ2D= $(SRC2D:.f90=.o)
OBJ3D= $(SRC3D:.f90=.o)
OBJNC= $(SRCNC:.f=.o)
PROG2= tx2
PROG3= tx
LIBS2= tx2lib.a ../lib/tasklib.a
LIBS3= tx3lib.a ../lib/tasklib.a

.SUFFIXES: .f90

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR)

#all: $(PROG2)
all: $(PROG3)

tx2lib.a: $(OBJS) $(OBJ2D) $(OBJNC) $(INCNC)
	$(LD) $(LDFLAGS) tx2lib.a $(OBJS) $(OBJ2D) $(OBJNC)

tx3lib.a: $(OBJS) $(OBJ3D) $(OBJNC) $(INCNC)
	$(LD) $(LDFLAGS) tx3lib.a $(OBJS) $(OBJ3D) $(OBJNC)

../lib/tasklib.a:
	(cd ../lib; make)

libs:
	(cd ../lib; make)

$(PROG2): $(LIBS2)
	$(FCFREE) $(LIBS2) $(LIBLA) -o $@ $(FFLAGS) $(FLIBS)

$(PROG3): $(LIBS3)
	$(FCFREE) $(LIBS3) $(LIBLA) -o $@ $(FFLAGS) $(FLIBS)

check: $(SRCS) $(SRC2D) $(SRC3D)
	ftnchek $(FCKFLAGS) $(SRCS) $(SRC2D) $(SRC3D) | less

clean:
	rm -f *.o ./*~ ./#* nclass/*.o *.a $(MOD)/*.mod

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

txcala.o: txcala.f90 txcomm.f90
txcalv.o: txcalv.f90 txcomm.f90
txexec.o: txexec.f90 txcomm.f90
txfile.o: txfile.f90 txcomm.f90
txgraf.o: txgraf.f90 txcomm.f90
txinit.o: txinit.f90 txcomm.f90
txlib.o : txlib.f90  txcomm.f90
txmain.o: txmain.f90 txcomm.f90
txmenu.o: txmenu.f90 txcomm.f90
trncls.o: txncls.f90 trncls.inc trcomm.inc trcom0.inc nclass/pamx_mi.inc \
          nclass/pamx_ms.inc nclass/pamx_mz.inc
txrslt.o: txrslt.f90 txcomm.f90
txrppl.o: txrppl.f90 txcomm.f90
nclass/nclass_mod.o : nclass/nclass_mod.f nclass/pamx_mi.inc \
                      nclass/pamx_ms.inc nclass/pamx_mz.inc
	$(FCFIXED) $(FFLAGS) -c nclass/nclass_mod.f -o nclass/nclass_mod.o -Inclass $(MODDIR)
