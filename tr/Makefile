## $Id$ ###
include ../make.header

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)

SRCS=trmenu.f trinit.f trgout.f trfile.f trloop.f trrslt.f trmdlt.f \
     trcalc.f trcoef.f trpnb.f trprf.f trpnf.f trpel.f \
     trgrar.f trgrat.f trgrap.f trgrae.f trgrad.f trgsub.f trhelp.f \
     glf/callglf2d.f glf/zgeev.f glf/blas_zgeev.f itg/ifspppl_model.f \
     tritg.f itg/weiland.f
SRCF=glf/glf2d.F glf/r8tomsqz.F
SRC2D=trg2d.f
SRC3D=trg3d.f
SRCNC=trncls.f nclass/nclass_mod.f nclass/rarray_sum.f \
      nclass/rarray_copy.f nclass/rarray_zero.f nclass/u_erf.f \
      nclass/u_lu_backsub.f nclass/u_lu_decomp.f nclass/write_mod.f
INCNC=trncls.inc nclass/pamx_mi.inc nclass/pamx_ms.inc nclass/pamx_mz.inc
INCGLF=trglf.inc glf/glf.m

OBJS=$(SRCS:.f=.o) $(SRCF:.F=.o)
OBJ2D=$(SRC2D:.f=.o)
OBJ3D=$(SRC3D:.f=.o)
OBJE=$(SRCE:.f=.o)
OBJL=$(SRCL:.f=.o)
OBJNC=$(SRCNC:.f=.o)

LIBS2=tr2lib.a ../pl/pllib.a ../eq/treqin.o ../eq/pleqin.o ../eq/eqlib.a ../lib/tasklib.a
LIBS3=tr3lib.a ../pl/pllib.a ../eq/treqin.o ../eq/pleqin.o ../eq/eqlib.a ../lib/tasklib.a

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

all: tr2 tr3

tr2lib.a: $(OBJS) $(OBJ2D) $(OBJE) $(OBJL) $(OBJNC) \
          trcomm.inc trcom0.inc $(INCNC) $(INCGLF)
	$(LD) $(LDFLAGS) -o tr2lib.a $(OBJS) $(OBJ2D) $(OBJE) $(OBJL) $(OBJNC)

tr3lib.a: $(OBJS) $(OBJ3D) $(OBJE) $(OBJL) $(OBJNC) \
          trcomm.inc trcom0.inc $(INCNC) $(INCGLF)
	$(LD) $(LDFLAGS) -o tr3lib.a $(OBJS) $(OBJ3D) $(OBJE) $(OBJL) $(OBJNC)

../pl/pllib.a:
	(cd ../pl; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../eq/treqin.o:
	(cd ../eq; make treqin.o)
../eq/pleqin.o:
	(cd ../eq; make pleqin.o)
../lib/tasklib.a:
	(cd ../lib; make)

libs:
	(cd ../pl; make pllib.a)
	(cd ../eq; make eqlib.a)
	(cd ../eq; make treqin.o)
	(cd ../eq; make pleqin.o)
	(cd ../lib; make)

tr2: trmain.o $(LIBS2)
	$(FC) trmain.o $(LIBS2) -o $@ $(FFLAGS) $(FLIBS) -lMdsLib
tr3: trmain.o $(LIBS3)
	$(FC) trmain.o $(LIBS3) -o $@ $(FFLAGS) $(GLLIBS) -lMdsLib

check : $(SRCS)
	ftnchek $(FCKFLAGS) $(SRCS) $(SRCE) $(SRCL) $(SRC2D) $(SRCNC) $(FCKPIPE) | less

clean : 
	rm -f ./#* ./*~ *.o nclass/*.o itg/*.o glf/*.o a.out core *.a

veryclean : clean
	rm -f tr2 tr3

new :
	-mkdir ../trnew
	cp Makefile ../trnew
	-cp -f HISTORY ../trnew
	cp *.f ../trnew
	cp *.inc ../trnew
	-mkdir ../trnew/itg
	-cp -f itg/*.f ../trnew/itg
	-mkdir ../trnew/glf
	-cp -f glf/*.f ../trnew/glf
	-cp -f glf/*.m ../trnew/glf
	-cp -f glf/*.F ../trnew/glf
	-cp -f glf/*.h ../trnew/glf
	-mkdir ../trnew/nclass
	-cp -f nclass/*.f ../trnew/nclass
	-cp -f nclass/*.inc ../trnew/nclass
	-mkdir ../trnew/doc
	-cp -f doc/* ../trnew/doc

trmain.o : trmain.f trcomm.inc trcom0.inc 
trmenu.o : trmenu.f trcomm.inc trcom0.inc 
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

trgout.o : trgout.f trcomm.inc trcom0.inc
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
trmdlt.o : trmdlt.f trcomm.inc trcom0.inc
glf/callglf2d.o : glf/callglf2d.f glf/glf.m
glf/glf2d.o  : glf/glf2d.F glf/glf.m glf/f77_dcomplx.h
	$(FC) $(FFLAGS) -c -o glf/glf2d.o -Iglf glf/glf2d.F
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
../eq/treqin.o: ../eq/treqin.f ../eq/eqcom1.inc ../eq/eqcom3.inc \
                ../eq/eqcom4.inc
