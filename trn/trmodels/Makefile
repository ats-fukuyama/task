TASKDIR=../../../task
BPSDDIR=../../../bpsd
include $(TASKDIR)/make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE=-I $(BPSDDIR)/mod -I $(TASKDIR)/plx/mod -I $(TASKDIR)/lib/mod \
           -I ../trcore/mod -I $(MOD)


SRC90=trcalnc.f90 trncls.f90 \
      trglf23.f90 trmbgb.f90 trmmm95.f90 trmmm7_1.f90 trcdbm.f90 \
      trsimple.f90 trcoeftb.f90 trcoefnc.f90
SRCMBGB=mbgb/mixed_Bohm_gyro_Bohm.f
SRCITG= glf/callglf2d.f  itg/ifspppl_model.f itg/weiland.f
SRCITGF= glf/glf2d.F glf/r8tomsqz.F
SRCCDBM=cdbm/cdbm.f90
SRCNC=nclass/nclass_mod.f nclass/rarray_sum.f \
      nclass/rarray_copy.f nclass/rarray_zero.f nclass/u_erf.f \
      nclass/u_lu_backsub.f nclass/u_lu_decomp.f nclass/write_mod.f
INCNC=trncls.inc nclass/pamx_mi.inc nclass/pamx_ms.inc nclass/pamx_mz.inc
INCGLF=glf/glf.m
SRCCYT=cytran/cytran_mod.f90 cytran/tr_cytran_mod.f90
SRCMMM=mmm95/stripx.f mmm95/tomsqz.f mmm95/weiland14.f mmm95/mmm95.f
SRCMMM7=mmm7_1/w20mod.f90 mmm7_1/modmmm7_1.f90

OBJS= $(SRCMMM7:.f90=.o) $(SRCMMM:.f=.o) $(SRCITGF:.F=.o) $(SRCITG:.f=.o) \
$(SRCMBGB:.f=.o) $(SRCCDBM:.f90=.o) $(SRC90:.f90=.o) $(SRCNC:.f=.o)
# $(SRCCYT:.f90=.o)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: trmodelslib.a

trmodelslib.a: $(OBJS) $(INCNC) $(INCGLF)
	$(LD) $(LDFLAGS) trmodelslib.a $(OBJS)

clean : 
	rm -f ./#* ./*f90~ *.o nclass/*.o itg/*.o glf/*.o cytran/*.o mmm7_1/*.o mmm95/*.o mbgb/*.o cdbm/*.o a.out core *.a ./mod/*.mod 

veryclean : 
	rm -f ./#* ./*~ *.o nclass/*.o itg/*.o glf/*.o cytran/*.o mmm7_1/*.o mmm95/*.o mbgb/*.o cdbm/*.o a.out core *.a *.mod

tritg.o  : tritg.f90
glf/callglf2d.o : glf/callglf2d.f glf/glf.m
glf/glf2d.o  : glf/glf2d.F glf/glf.m glf/f77_dcomplx.h
	$(FCFIXED) $(FFLAGS) -c -o glf/glf2d.o -Iglf glf/glf2d.F
glf/r8tomsqz.o : glf/r8tomsqz.F glf/f77_dcomplx.h
	$(FCFIXED) $(FFLAGS) -c -o glf/r8tomsqz.o glf/r8tomsqz.F
glf/zgeev.o  : glf/zgeev.f
glf/blas_zgeev.o : glf/blas_zgeev.f
trncls.o : trncls.f90 trncls.inc nclass/pamx_mi.inc \
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
	$(FCFIXED) -r8 $(FFLAGS) -c $? -o $@
mmm95/tomsqz.o: mmm95/tomsqz.f
	$(FCFIXED) -r8 $(FFLAGS) -c $? -o $@
mmm95/weiland14.o: mmm95/weiland14.f
	$(FCFIXED) -r8 $(FFLAGS) -c $? -o $@
mmm95/mmm95.o: mmm95/mmm95.f
	$(FCFIXED) -r8 $(FFLAGS) -c $? -o $@
mmm7_1/w20mod.o: mmm7_1/w20mod.f90
mmm7_1/modmmm7_1.o: mmm7_1/w20mod.f90 mmm7_1/modmmm7_1.f90
mbgb/mixed_Bohm_gyro_Bohm.o: mbgb/mixed_Bohm_gyro_Bohm.f
