### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)
FFLAGS = $(OFLAGS)
#LFLAGS = $(DFLAGS)
#FFLAGS = $(DFLAGS)

SRCS=trmain.f trinit.f trfile.f trloop.f trrslt.f \
     trcalc.f trcoef.f trpnb.f trprf.f trpnf.f trpel.f \
     trgrar.f trgrat.f trgrap.f trgrae.f trgrad.f trgsub.f trhelp.f \
     glf/callglf2d.f glf/zgeev.f glf/blas_zgeev.f itg/ifspppl_model.f \
     tritg.f itg/weiland.f
SRCF=glf/glf2d.F glf/r8tomsqz.F
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
INCGLF=trglf.inc glf/glf.m

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

$(PROG3) : $(OBJS) $(OBJ3D) $(OBJE) $(OBJL) $(OBJNC) trcomm.inc trcom0.inc $(INCNC) $(INCGLF)
	-$(FC) -o $@ $(FFLAGS) $(OBJS) $(OBJ3D) $(OBJE) $(OBJL) $(OBJNC) $(GLLIBS)

2d : $(PROG2)

$(PROG2) : $(OBJS) $(OBJ2D) $(OBJE) $(OBJL) $(OBJNC) trcomm.inc trcom0.inc $(INCNC) $(INCGLF)
	$(FC) -o $@ $(FFLAGS) $(OBJS) $(OBJ2D) $(OBJE) $(OBJL) $(OBJNC) $(FLIBS)

check : $(SRCS)
	ftnchek $(FCKFLAGS) $(SRCS) $(SRCE) $(SRCL) $(SRC2D) $(SRC3D) $(SRCNC) $(FCKPIPE) | less

clean : 
	rm -f ./#* ./*~ *.o nclass/*.o itg/*.o glf/*.o a.out core

veryclean : clean
	rm -f $(PROG2) $(PROG3)

new :
	-mkdir ../trnew
	cp Makefile ../trnew
	-cp -f HISTORY ../trnew
	cp *.f ../trnew
	cp *.inc ../trnew
	-mkdir ../trnew/in
	-cp -f in/* ../trnew/in
	-mkdir ../trnew/doc
	-cp -f doc/* ../trnew/doc

trmain.o : trmain.f trcomm.inc trcom0.inc 
trinit.o : trinit.f trcomm.inc trcom0.inc
trfile.o : trfile.f trcomm.inc trcom0.inc
trloop.o : trloop.f trcomm.inc trcom0.inc
trrslt.o : trrslt.f trcomm.inc trcom0.inc

trcalc.o : trcalc.f trcomm.inc trcom0.inc
trcoef.o : trcoef.f trcomm.inc trcom0.inc trglf.inc
trpnb.o  : trpnb.f  trcomm.inc trcom0.inc
trprf.o  : trprf.f  trcomm.inc trcom0.inc
trpnf.o  : trpnf.f  trcomm.inc trcom0.inc
trpel.o  : trpel.f  trcomm.inc trcom0.inc

trgrar.o : trgrar.f trcomm.inc trcom0.inc
trgrat.o : trgrat.f trcomm.inc trcom0.inc
trgrap.o : trgrap.f trcomm.inc trcom0.inc
trgrae.o : trgrae.f trcomm.inc trcom0.inc
trgrad.o : trgrad.f trcomm.inc trcom0.inc
trgsub.o : trgsub.f trcomm.inc trcom0.inc
trhelp.o : trhelp.f trcomm.inc trcom0.inc
trgrad.o : trgrad.f trcomm.inc trcom0.inc
trg2d.o  : trg2d.f  trcomm.inc trcom0.inc
trg3d.o  : trg3d.f  trcomm.inc trcom0.inc
tritg.o  : tritg.f  trcomm.inc trcom0.inc
glf/callglf2d.o : glf/callglf2d.f glf/glf.m
glf/glf2d.o  : glf/glf2d.F glf/glf.m glf/f77_dcomplx.h
	$(FC) $(FFLAGS) -c -o glf/glf2d.o glf/glf2d.F
glf/r8tomsqz.o : glf/r8tomsqz.F glf/f77_dcomplx.h
	$(FC) $(FFLAGS) -c -o glf/r8tomsqz.o glf/r8tomsqz.F
glf/zgeev.o  : glf/zgeev.f
glf/blas_zgeev.o : glf/blas_zgeev.f
trncls.o : trncls.f trncls.inc trcomm.inc trcom0.inc nclass/pamx_mi.inc \
           nclass/pamx_ms.inc nclass/pamx_mz.inc
nclass/nclass_mod.o : nclass/nclass_mod.f nclass/pamx_mi.inc \
                      nclass/pamx_ms.inc nclass/pamx_mz.inc
itg/ifspppl_model.o : itg/ifspppl_model.f
itg/weiland.o: itg/weiland.f
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
../eq/eqdsk.o : ../eq/eqdsk.f ../eq/eqcom1.inc ../eq/eqcom3.inc
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
