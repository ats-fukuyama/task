### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)
FFLAGS = $(OFLAGS)
#LFLAGS = $(DFLAGS)
#FFLAGS = $(DFLAGS)

SRCS=trmain.f trinit.f trfile.f trloop.f trrslt.f \
     trcalc.f trcoef.f trpnb.f trprf.f trpnf.f trpel.f \
     trgrar.f trgrat.f trgrap.f trgrae.f trgrad.f trgsub.f trhelp.f \
     callglf2d.f zgeev.f blas_zgeev.f ifspppl/ifspppl_model.f
SRCF=glf2d.F r8tomsqz.F
SRC2D=trg2d.f
SRC3D=trg3d.f
SRCE=../eq/treqin.f ../eq/eqinit.f ../eq/eqcalc.f ../eq/eqfunc.f \
     ../eq/eqcalq.f ../eq/eqload.f ../eq/eqgout.f ../eq/eqgsub.f \
     ../eq/eqdsk.f
SRCL=../lib/libbnd.f ../lib/libspl.f ../lib/libbrn.f ../lib/libgut.f \
     ../lib/libchar.f ../lib/libfio.f
SRCNC=trncls.f nclass/nclass_mod.f nclass/rarray_sum.f \
      nclass/rarray_copy.f nclass/rarray_zero.f nclass/u_erf.f \
      nclass/u_lu_backsub.f nclass/u_lu_decomp.f
INCNC=trncls.inc nclass/pamx_mi.inc nclass/pamx_ms.inc nclass/pamx_mz.inc
INCGLF=trglf.h glf.m

OBJS=$(SRCS:.f=.o) $(SRCF:.F=.o)
OBJ2D=$(SRC2D:.f=.o)
OBJ3D=$(SRC3D:.f=.o)
OBJE=$(SRCE:.f=.o)
OBJL=$(SRCL:.f=.o)
OBJNC=$(SRCNC:.f=.o)
PROG2=tr2
PROG3=trf

all: 3d 2d

3d: $(PROG3)

$(PROG3) : $(OBJS) $(OBJ3D) $(OBJE) $(OBJL) $(OBJNC) trcomm.h $(INCNC) $(INCGLF)
	-$(FC) -o $@ $(FFLAGS) $(OBJS) $(OBJ3D) $(OBJE) $(OBJL) $(OBJNC) $(GLLIBS)

2d : $(PROG2)

$(PROG2) : $(OBJS) $(OBJ2D) $(OBJE) $(OBJL) $(OBJNC) trcomm.h $(INCNC) $(INCGLF)
	$(FC) -o $@ $(FFLAGS) $(OBJS) $(OBJ2D) $(OBJE) $(OBJL) $(OBJNC) $(FLIBS)

check : $(SRCS)
	ftnchek $(FCKFLAGS) $(SRCS) $(SRCE) $(SRCL) $(FCKPIPE) | less

clean : 
	rm -f ./#* ./*~ *.o nclass/*.o a.out core

veryclean : clean
	rm -f $(PROG) $(PROG3)

new :
	-mkdir ../trnew
	cp Makefile ../trnew
	-cp -f HISTORY ../trnew
	cp *.f ../trnew
	cp *.h ../trnew
	-mkdir ../trnew/in
	-cp -f in/* ../trnew/in
	-mkdir ../trnew/doc
	-cp -f doc/* ../trnew/doc

trmain.o : trmain.f trcomm.h
trinit.o : trinit.f trcomm.h
trfile.o : trfile.f trcomm.h
trloop.o : trloop.f trcomm.h
trrslt.o : trrslt.f trcomm.h

trcalc.o : trcalc.f trcomm.h
trcoef.o : trcoef.f trcomm.h trglf.h
trpnb.o  : trpnb.f  trcomm.h
trprf.o  : trprf.f  trcomm.h
trpnf.o  : trpnf.f  trcomm.h
trpel.o  : trpel.f  trcomm.h

trgrar.o : trgrar.f trcomm.h
trgrat.o : trgrat.f trcomm.h
trgrap.o : trgrap.f trcomm.h
trgrae.o : trgrae.f trcomm.h
trgsub.o : trgsub.f trcomm.h
trhelp.o : trhelp.f trcomm.h
trgrad.o : trgrad.f trcomm.h
trgra2.o : trgra2.f trcomm.h
callglf2d.o : callglf2d.f glf.m
glf2d.o  : glf2d.F glf.m f77_dcomplx.h
	$(FC) $(FFLAGS) -c -o glf2d.o glf2d.F
r8tomsqz.o : r8tomsqz.F f77_dcomplx.h
	$(FC) $(FFLAGS) -c -o r8tomsqz.o r8tomsqz.F
zgeev.o  : zgeev.f
blas_zgeev.o : blas_zgeev.f
trncls.o : trncls.f trncls.inc trcomm.h nclass/pamx_mi.inc \
           nclass/pamx_ms.inc nclass/pamx_mz.inc
nclass/nclass_mod.o : nclass/nclass_mod.f nclass/pamx_mi.inc \
                      nclass/pamx_ms.inc nclass/pamx_mz.inc
ifspppl/ifspppl_model.o : ifspppl/ifspppl_model.f
../eq/eqinit.o: ../eq/eqinit.f ../eq/eqcom1.inc
../eq/eqcalc.o: ../eq/eqcalc.f ../eq/eqcom1.inc ../eq/eqcom2.inc
../eq/eqcalq.o: ../eq/eqcalq.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../eq/eqgout.o: ../eq/eqgout.f ../eq/eqcom1.inc ../eq/eqcom2.inc \
                ../eq/eqcom3.inc
../eq/eqfunc.o: ../eq/eqfunc.f ../eq/eqcom1.inc ../eq/eqcom2.inc \
                ../eq/eqcom4.inc
../eq/treqin.o: ../eq/treqin.f ../eq/eqcom1.inc ../eq/eqcom3.inc \
                ../eq/eqcom4.inc
../eq/eqload.o: ../eq/eqcalq.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../eq/eqdsk.o: ../eq/eqdsk.f ../eq/eqcom1.inc ../eq/eqcom3.inc
../lib/libbnd.o: ../lib/libbnd.f
	$(FC) -o $@ $(LFLAGS) -c ../lib/libbnd.f
../lib/libspl.o: ../lib/libspl.f
	$(FC) -o $@ $(LFLAGS) -c ../lib/libspl.f
../lib/libbrn.o: ../lib/libbrn.f
	$(FC) -o $@ $(LFLAGS) -c ../lib/libbrn.f
../lib/libgut.o: ../lib/libgut.f
	$(FC) -o $@ $(LFLAGS) -c ../lib/libgut.f
../lib/libchar.o: ../lib/libchar.f
	$(FC) -o $@ $(LFLAGS) -c ../lib/libchar.f
../lib/libfio.o: ../lib/libfio.f
	$(FC) -o $@ $(LFLAGS) -c ../lib/libfio.f
