### $Id$ ###
include ../make.header

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)

SRC3D= txg3d.f90
SRCS = txcnst.f90 txcomm.f90 txlib.f90 txncls.f90 txrslt.f90 txcalv.f90 txcala.f90 \
       txinit.f90 txfile.f90 txgraf.f90 txexec.f90 txmenu.f90 txmain.f90
SRCNC= nclass/nclass_mod.f nclass/rarray_sum.f \
       nclass/rarray_copy.f nclass/rarray_zero.f nclass/u_erf.f \
       nclass/u_lu_backsub.f nclass/u_lu_decomp.f nclass/write_mod.f
INCNC=txncls.inc nclass/pamx_mi.inc nclass/pamx_ms.inc nclass/pamx_mz.inc
OBJ3D= $(SRC3D:.f90=.o)
OBJS = $(SRCS:.f90=.o)
OBJNC= $(SRCNC:.f=.o)
PROG = tx
LIBS = txlib.a ../lib/tasklib.a $(OBJLA)

.SUFFIXES: .f90

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ -module mod
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ -module mod

all: $(PROG)

txlib.a: $(OBJ3D) $(OBJS) $(OBJNC) $(INCNC)
	$(LD) $(LDFLAGS) -o txlib.a $(OBJ3D) $(OBJS) $(OBJNC) -Lmod

../lib/tasklib.a:
	(cd ../lib; make)

libs:
	(cd ../lib; make)

$(PROG): $(LIBS)
#	$(FCFREE) $(LIBS) -o $@ $(FFLAGS) $(GFLIBS) $(LIBLA)
	$(FCFREE) $(LIBS) -o $@ $(FFLAGS) $(GLLIBS) $(LIBLA)

check: $(SRCS)
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	rm -f *.o ./*~ ./#* nclass/*.o *.a mod/*.mod

veryclean: clean
	rm -f $(PROG)

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
nclass/nclass_mod.o : nclass/nclass_mod.f nclass/pamx_mi.inc \
                      nclass/pamx_ms.inc nclass/pamx_mz.inc
	$(FCFIXED) $(FFLAGS) -c nclass/nclass_mod.f -o nclass/nclass_mod.o -Inclass -module mod
