### Makefile
include ../make.header
include ../mtxp/make.mtxp
#LIB_MTX=$(LIB_MTX_BND)
#LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_PCG)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIB_MTX=$(LIBXMTX_KSP)
#LIBX_MTX=$(LIBX_MTX_KSP)

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)

OBJDIR=./obj

${OBJDIR}/%.o:%.f
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) \
            -I../lib/$(MOD) -I../mtxp/$(MOD) \
            -I../../bpsd/$(MOD)

SRCITG= glf/callglf2d.f \
        itg/ifspppl_model.f itg/weiland.f
SRCITGF=glf/glf2d.F glf/r8tomsqz.F
SRCNC=nclass/nclass_mod.f nclass/rarray_sum.f \
      nclass/rarray_copy.f nclass/rarray_zero.f nclass/u_erf.f \
      nclass/u_lu_backsub.f nclass/u_lu_decomp.f nclass/write_mod.f
INCNC=nclass/pamx_mi.inc nclass/pamx_ms.inc nclass/pamx_mz.inc
INCGLF=glf/glf.m
SRCCYT=cytran/cytran_mod.f90
SRCMMM95=mmm95/stripx.f mmm95/tomsqz.f mmm95/weiland14.f mmm95/mmm95.f
SRCMMM7=libmmm7_1/w20mod.f90 libmmm7_1/modmmm7_1.f90
SRCMBGB=mbgb/mixed_Bohm_gyro_Bohm.f
SRCCDBM=trcdbm.f90

OBJNC=$(addprefix $(OBJDIR)/, $(SRCNC:.f=.o))
OBJCYT=$(addprefix $(OBJDIR)/, $(SRCCYT:.f90=.o))
OBJMMM=$(addprefix $(OBJDIR)/, $(SRCMMM95:.f=.o) $(SRCMMM7:.f90=.o))
OBJMBGB=$(addprefix $(OBJDIR)/, $(SRCMBGB:.f=.o))
OBJITG=$(addprefix $(OBJDIR)/, $(SRCITG:.f=.o) $(SRCITGF:.F=.o))
OBJCDBM=$(addprefix $(OBJDIR)/, $(SRCCDBM:.f90=.o))

OBJTR=$(OBJNC) $(OBJCYT) $(OBJMMM) $(OBJMBGB) $(OBJITG) $(OBJCDBM)

LIBS=../lib/libtask.a ../lib/libgrf.a  ../../bpsd/libbpsd.a

all: libs libtrm.a

libtrm.a: $(OBJTR) $(INCNC) $(INCGLF)
	$(LD) $(LDFLAGS) libtrm.a $(OBJTR)

../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)

clean : 
	-rm -f ./#* ./*~ $(OBJDIR)/*.o $(OBJDIR)/*/*.o $(MOD)/*.mod $(MOD)/*/*.mod core *.a

veryclean : clean
	-rm -f libtrm.a

$(OBJDIR)/glf/callglf2d.o : glf/callglf2d.f glf/glf.m
	$(FCFIXED) $(FFLAGS) -c -o $(OBJDIR)/glf/callglf2d.o -Iglf glf/callglf2d.f
$(OBJDIR)/glf/glf2d.o  : glf/glf2d.F glf/glf.m glf/f77_dcomplx.h
	$(FCFIXED) $(FFLAGS) -c -o $(OBJDIR)/glf/glf2d.o -Iglf glf/glf2d.F
$(OBJDIR)/glf/r8tomsqz.o : glf/r8tomsqz.F glf/f77_dcomplx.h
	$(FCFIXED) $(FFLAGS) -c -o $(OBJDIR)/glf/r8tomsqz.o glf/r8tomsqz.F
$(OBJDIR)/glf/zgeev.o  : glf/zgeev.f
$(OBJDIR)/glf/blas_zgeev.o : glf/blas_zgeev.f
$(OBJDIR)/nclass/nclass_mod.o : nclass/nclass_mod.f nclass/pamx_mi.inc \
                      nclass/pamx_ms.inc nclass/pamx_mz.inc
	$(FCFIXED) $(FFLAGS) -c -o $(OBJDIR)/nclass/nclass_mod.o -Inclass nclass/nclass_mod.f $(MODDIR) $(MODINCLUDE)
$(OBJDIR)/itg/ifspppl_model.o : itg/ifspppl_model.f
$(OBJDIR)/itg/weiland.o: itg/weiland.f
$(OBJDIR)/nclass/rarray_sum.o : nclass/rarray_sum.f
$(OBJDIR)/nclass/rarray_copy.o: nclass/rarray_copy.f
$(OBJDIR)/nclass/rarray_zero.o: nclass/rarray_zero.f
nclass/u_erf.o      : nclass/u_erf.f
$(OBJDIR)/nclass/u_lu_backsub.o: nclass/u_lu_backsub.f
$(OBJDIR)/nclass/u_lu_decomp.o : nclass/u_lu_decomp.f
$(OBJDIR)/nclass/write_mod.o   : nclass/write_mod.f
$(OBJDIR)/cytran/cytran_mod.o: cytran/cytran_mod.f90
$(OBJDIR)/mmm95/stripx.o: mmm95/stripx.f
$(OBJDIR)/mmm95/tomsqz.o: mmm95/tomsqz.f
$(OBJDIR)/mmm95/weiland14.: mmm95/weiland14.f
$(OBJDIR)/mmm95/mmm95.o: mmm95/mmm95.f
$(OBJDIR)/libmmm7_1/w20mod.o: libmmm7_1/w20mod.f90
$(OBJDIR)/libmmm7_1/modmmm7_1.o: libmmm7_1/modmmm7_1.f90 libmmm7_1/w20mod.f90
$(OBJDIR)/mbgb/mixed_Bohm_gyro_Bohm.o: mbgb/mixed_Bohm_gyro_Bohm.f
$(OBJDIR)/trcdbm.o: trcdbm.f90
