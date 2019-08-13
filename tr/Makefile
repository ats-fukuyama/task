### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = -g $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE= -I$(MOD) -I../../bpsd/$(MOD) -I../lib/$(MOD) -I../eq/$(MOD) \
	    -I../pl/$(MOD)
EQLIBS=../eq/libeq.a
PLLIBS=../../bpsd/libbpsd.a ../pl/libpl.a 

SRC90M=trcom0.f90 trcom1.f90 trcomm.f90 tr_bpsd.f90
SRC90=trinit.f90 trparm.f90 trprep.f90  \
      trprof.f90 trufsub.f90 trufile.f90 trfile.f90 trhelp.f90 \
      trexec.f90 trcalc.f90 \
      tradat.f90 cdbm.f90 trmodels.f90 trcoef.f90 tritg.f90 \
      trrslt.f90 trpnb.f90 trprf.f90 trpnf.f90 trpel.f90 trpsc.f90 trmdlt.f90 \
      trgout.f90 trgrar.f90 trgrat.f90 trgrap.f90 trgrae.f90 \
      trgrad.f90 trgram.f90 trgsub.f90 trfout.f90 \
      trunit.f90 trloop.f90 trmenu.f90

SRCITG= glf/callglf2d.f \
        itg/ifspppl_model.f itg/weiland.f
# glf/zgeev.f glf/blas_zgeev.f
SRCITGF= glf/glf2d.F glf/r8tomsqz.F
#SRCITG= noitg.f
#SRCITGF=
SRC2D=trg2d.f90
SRC3D=trg3d.f90
SRCNC=trncls.f nclass/nclass_mod.f nclass/rarray_sum.f \
      nclass/rarray_copy.f nclass/rarray_zero.f nclass/u_erf.f \
      nclass/u_lu_backsub.f nclass/u_lu_decomp.f nclass/write_mod.f
INCNC=trncls.inc nclass/pamx_mi.inc nclass/pamx_ms.inc nclass/pamx_mz.inc
INCGLF=trglf.inc glf/glf.m
SRCCYT=cytran/cytran_mod.f90 cytran/tr_cytran_mod.f90
SRCMMM95=mmm95/stripx.f mmm95/tomsqz.f mmm95/weiland14.f mmm95/mmm95.f
SRCMMM7=libmmm7_1/w20mod.f90 libmmm7_1/modmmm7_1.f90
SRCMBGB=mbgb/mixed_Bohm_gyro_Bohm.f

OBJM= $(SRC90M:.f90=.o)
OBJS= $(SRCITG:.f=.o) $(SRCITGF:.F=.o) $(SRC90:.f90=.o) 
OBJE=$(SRCE:.f=.o)
OBJL=$(SRCL:.f=.o)
OBJNC=$(SRCNC:.f=.o)
OBJCYT=$(SRCCYT:.f90=.o)
OBJMMM=$(SRCMMM95:.f=.o) $(SRCMMM7:.f90=.o)
OBJMBGB=$(SRCMBGB:.f=.o)
OBJTR=$(OBJM) $(OBJE) $(OBJL) $(OBJNC) $(OBJCYT) $(OBJMMM) $(OBJMBGB) $(OBJS)
OBJ2D=$(SRC2D:.f90=.o)
OBJ3D=$(SRC3D:.f90=.o)

LIBS2=../lib/libtask.a ../lib/libmds.a $(PLLIBS) $(EQLIBS) libtr2.a
LIBS3=../lib/libtask.a ../lib/libmds.a $(PLLIBS) $(EQLIBS) libtr3.a  

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

#all: libs tr2 tr3
all: libs tr2

libtr2.a: $(OBJTR) $(OBJ2D) $(INCNC) $(INCGLF)
	$(LD) $(LDFLAGS) libtr2.a $(OBJTR) $(OBJ2D)

libtr3.a: $(OBJTR) $(OBJ3D) $(INCNC) $(INCGLF)
	$(LD) $(LDFLAGS) libtr3.a $(OBJTR) $(OBJ3D)

../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libmds.a:
	(cd ../lib; make libmds.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make eqlib.a)

libs:
	(cd ../lib; make libtask.a)
	(cd ../../bpsd; make libbpsd.a)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)

tr2: $(LIBS2) trmain.o 
	$(FCFREE) trmain.o $(LIBS2) -o $@ $(FFLAGS) $(FLIBS) $(MODULE) $(MODO) $(MODINCLUDE)

tr3: $(LIBS3) trmain.o 
	$(FCFREE) trmain.o $(LIBS3) -o $@ $(FFLAGS) $(FLIBS) $(MODINCLUDE)

xtr: $(LIBS2) xtr.o 
	$(FCFREE) xtr.o $(LIBS2) -o $@ $(FFLAGS) $(FLIBS) $(MODULE) $(MODO) $(MODINCLUDE)

check : $(SRCS)
	ftnchek $(FCKFLAGS) $(SRCS) $(SRCE) $(SRCL) $(SRC2D) $(SRCNC) $(FCKPIPE) | less

clean : 
	-rm -f ./#* ./*~ *.o nclass/*.o itg/*.o glf/*.o cytran/*.o mmm95/*.o libmmm7_1/*.o mbgb/*.o a.out core *.a *.mod $(MOD)/*.mod

veryclean : clean
	-rm -f tr2 tr3

distclean: clean
	cd ../pl && make clean
	cd ../eq && make clean
	cd ../lib && make clean

trmain.o : trmain.f90 trcomm.f90  
trmenu.o : trmenu.f90 trcomm.f90 trcom1.f90 trunit.f90
trinit.o : trinit.f90 trcomm.f90 trcom1.f90
trparm.o : trparm.f90 trcomm.f90 trcom1.f90
trprep.o : trprep.f90 trcomm.f90 trcom1.f90
trprof.o : trprof.f90 trcomm.f90 trcom1.f90
trfile.o : trfile.f90 trcomm.f90 trcom1.f90
trufile.o : trufile.f90 trcomm.f90 trcom1.f90
trufsub.o : trufsub.f90 trcomm.f90 trcom1.f90
trloop.o : trloop.f90 trcomm.f90 trcom1.f90
trexec.o : trexec.f90 trcomm.f90 trcom1.f90
trpl.o : trpl.f90 trcomm.f90 trcom1.f90
trrslt.o : trrslt.f90 trcomm.f90
trcalc.o : trcalc.f90 trcomm.f90 trcom1.f90
trcoef.o : trcoef.f90 trcomm.f90 trglf.inc
trpnb.o  : trpnb.f90  trcomm.f90 
trprf.o  : trprf.f90  trcomm.f90 
trpnf.o  : trpnf.f90  trcomm.f90 
trpel.o  : trpel.f90  trcomm.f90 
trpsc.o  : trpsc.f90  trcomm.f90 

trgout.o : trgout.f90 trcomm.f90 
trgrar.o : trgrar.f90 trcomm.f90 
trgrat.o : trgrat.f90 trcomm.f90 trcom1.f90
trgrap.o : trgrap.f90 trcomm.f90 
trgrae.o : trgrae.f90 trcomm.f90 
trgrad.o : trgrad.f90 trcomm.f90 
trgram.o : trgram.f90 trcomm.f90 trcom1.f90
trgsub.o : trgsub.f90 trcomm.f90 
trhelp.o : trhelp.f90 trcomm.f90 
trg2d.o  : trg2d.f90  trcomm.f90 
trg3d.o  : trg3d.f90  trcomm.f90 
tritg.o  : tritg.f90  trcomm.f90 
trmodels.o : trmodels.f90  trcomm.f90 
trmdlt.o : trmdlt.f90 trcomm.f90 
glf/callglf2d.o : glf/callglf2d.f glf/glf.m
glf/glf2d.o  : glf/glf2d.F glf/glf.m glf/f77_dcomplx.h
	$(FCFIXED) $(FFLAGS) -c -o glf/glf2d.o -Iglf glf/glf2d.F
glf/r8tomsqz.o : glf/r8tomsqz.F glf/f77_dcomplx.h
	$(FCFIXED) $(FFLAGS) -c -o glf/r8tomsqz.o glf/r8tomsqz.F
glf/zgeev.o  : glf/zgeev.f
glf/blas_zgeev.o : glf/blas_zgeev.f
trncls.o : trncls.f90 trncls.inc trcomm.f90 nclass/pamx_mi.inc \
           nclass/pamx_ms.inc nclass/pamx_mz.inc
nclass/nclass_mod.o : nclass/nclass_mod.f nclass/pamx_mi.inc \
                      nclass/pamx_ms.inc nclass/pamx_mz.inc
itg/ifspppl_model.o : itg/ifspppl_model.f
itg/weiland.o: itg/weiland.f
../eq/treqin.o: ../eq/treqin.f90 ../eq/eqcom1.f90 ../eq/eqcom3.f90 \
                ../eq/eqcom4.f90
nclass/rarray_sum.o : nclass/rarray_sum.f
nclass/rarray_copy.o: nclass/rarray_copy.f
nclass/rarray_zero.o: nclass/rarray_zero.f
nclass/u_erf.o      : nclass/u_erf.f
nclass/u_lu_backsub.o: nclass/u_lu_backsub.f
nclass/u_lu_decomp.o : nclass/u_lu_decomp.f
nclass/write_mod.o   : nclass/write_mod.f
cytran/cytran_mod.o: cytran/cytran_mod.f90
cytran/tr_cytran_mod.o: cytran/tr_cytran_mod.f90
mmm95/stripx.o: mmm95/stripx.f
mmm95/tomsqz.o: mmm95/tomsqz.f
mmm95/weiland14.: mmm95/weiland14.f
mmm95/mmm95.o: mmm95/mmm95.f
libmmm7_1/w20mod.o: libmmm7_1/w20mod.f90
libmmm7_1/modmmm7_1.o: libmmm7_1/modmmm7_1.f90 libmmm7_1/w20mod.f90
mbgb/mixed_Bohm_gyro_Bohm.o: mbgb/mixed_Bohm_gyro_Bohm.f
