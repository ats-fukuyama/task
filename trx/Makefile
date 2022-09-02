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

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

OBJDIR=./obj

${OBJDIR}/%.o:%.f
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../pl/$(MOD) -I../eq/$(MOD) \
            -I../lib/$(MOD) -I../mtxp/$(MOD) \
            -I../../bpsd/$(MOD)

SRCM=trcom0.f90 trcom1.f90 trcomm.f90 trbpsd.f90
SRCS=trinit.f90 trparm.f90 \
     trmetric.f90 trprof.f90 trprep.f90 \
     trfile.f90 trhelp.f90 \
     trexec.f90 trcalc.f90 \
     tradat.f90 cdbm.f90 trmodels.f90 trcoef.f90 tritg.f90 \
     trrslt.f90 trpnb.f90 trprf.f90 trpnf.f90 trpel.f90 trpsc.f90 \
     trgout.f90 trgrar.f90 trgrat.f90 trgrap.f90 trgrae.f90 \
     trgrad.f90 trgram.f90 trgsub.f90 trfout.f90 \
     trloop.f90 trmenu.f90

SRC2D=trg2d.f90
SRC3D=trg3d.f90

SRCITG= glf/callglf2d.f \
        itg/ifspppl_model.f itg/weiland.f
SRCITGF=glf/glf2d.F glf/r8tomsqz.F
SRCNC=trncls.f nclass/nclass_mod.f nclass/rarray_sum.f \
      nclass/rarray_copy.f nclass/rarray_zero.f nclass/u_erf.f \
      nclass/u_lu_backsub.f nclass/u_lu_decomp.f nclass/write_mod.f
INCNC=trncls.inc nclass/pamx_mi.inc nclass/pamx_ms.inc nclass/pamx_mz.inc
INCGLF=trglf.inc glf/glf.m
SRCCYT=cytran/cytran_mod.f90 cytran/tr_cytran_mod.f90
SRCMMM95=mmm95/stripx.f mmm95/tomsqz.f mmm95/weiland14.f mmm95/mmm95.f
SRCMMM7=libmmm7_1/w20mod.f90 libmmm7_1/modmmm7_1.f90
SRCMBGB=mbgb/mixed_Bohm_gyro_Bohm.f

OBJM=$(addprefix $(OBJDIR)/, $(SRCM:.f90=.o))
OBJE=$(addprefix $(OBJDIR)/, $(SRCE:.f=.o))
OBJL=$(addprefix $(OBJDIR)/, $(SRCL:.f=.o))
OBJNC=$(addprefix $(OBJDIR)/, $(SRCNC:.f=.o))
OBJCYT=$(addprefix $(OBJDIR)/, $(SRCCYT:.f90=.o))
OBJMMM=$(addprefix $(OBJDIR)/, $(SRCMMM95:.f=.o) $(SRCMMM7:.f90=.o))
OBJMBGB=$(addprefix $(OBJDIR)/, $(SRCMBGB:.f=.o))
OBJITG=$(addprefix $(OBJDIR)/, $(SRCITG:.f=.o) $(SRCITGF:.F=.o))
OBJS=$(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))
OBJ2D=$(addprefix $(OBJDIR)/, $(SRC2D:.f90=.o))
OBJ3D=$(addprefix $(OBJDIR)/, $(SRC3D:.f90=.o))

OBJTR=$(OBJM) $(OBJE) $(OBJL) $(OBJNC) $(OBJCYT) $(OBJMMM) $(OBJMBGB) $(OBJITG) $(OBJS)

LIBS=../eq/libeq.a ../pl/libpl.a ../lib/libmds.a ../lib/libtask.a ../lib/libgrf.a  ../../bpsd/libbpsd.a
LIBS2=libtr2.a $(LIBS)
LIBS3=libtr3.a $(LIBS)

#all: libs tr2 tr3
all: libs tr2

libtr2.a: $(OBJTR) $(OBJ2D) $(INCNC) $(INCGLF)
	$(LD) $(LDFLAGS) libtr2.a $(OBJTR) $(OBJ2D)

libtr3.a: $(OBJTR) $(OBJ3D) $(INCNC) $(INCGLF)
	$(LD) $(LDFLAGS) libtr3.a $(OBJTR) $(OBJ3D)

../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libmds.a:
	(cd ../lib; make libmds.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)

tr2: $(LIBS2) $(OBJDIR)/trmain.o 
	$(FLINKER) $(OBJDIR)/trmain.o $(LIBS2) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

tr3: $(LIBS3) $(OBJDIR)/trmain.o 
	$(FLINKER) $(OBJDIR)/trmain.o $(LIBS3) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean : 
	-rm -f ./#* ./*~ $(OBJDIR)/*.o $(OBJDIR)/*/*.o $(MOD)/*.mod $(MOD)/*/*.mod

veryclean : clean
	-rm -f tr2 tr3

distclean: clean
	cd ../pl && make clean
	cd ../eq && make clean
	cd ../lib && make clean

$(OBJDIR)/trmain.o : trmain.f90 trcomm.f90  
$(OBJDIR)/trmenu.o : trmenu.f90 trcomm.f90 trcom1.f90
$(OBJDIR)/trinit.o : trinit.f90 trcomm.f90 trcom1.f90
$(OBJDIR)/trparm.o : trparm.f90 trcomm.f90 trcom1.f90
$(OBJDIR)/trprep.o : trprep.f90 trcomm.f90 trcom1.f90
$(OBJDIR)/trprof.o : trprof.f90 trcomm.f90 trcom1.f90
$(OBJDIR)/trmetric.o : trmetric.f90 trcomm.f90 trcom1.f90
$(OBJDIR)/trfile.o : trfile.f90 trcomm.f90 trcom1.f90
$(OBJDIR)/trbpsd.o : trbpsd.f90 trcomm.f90
$(OBJDIR)/tradat.o : tradat.f90 trcomm.f90
$(OBJDIR)/trloop.o : trloop.f90 trcomm.f90 trcom1.f90
$(OBJDIR)/trexec.o : trexec.f90 trcomm.f90 trcom1.f90
$(OBJDIR)/trpl.o : trpl.f90 trcomm.f90 trcom1.f90
$(OBJDIR)/trrslt.o : trrslt.f90 trcomm.f90
$(OBJDIR)/trcalc.o : trcalc.f90 trcomm.f90 trcom1.f90 cytran/tr_cytran_mod.f90
$(OBJDIR)/trcoef.o : trcoef.f90 trcomm.f90 trglf.inc
$(OBJDIR)/trpnb.o  : trpnb.f90  trcomm.f90 
$(OBJDIR)/trprf.o  : trprf.f90  trcomm.f90 
$(OBJDIR)/trpnf.o  : trpnf.f90  trcomm.f90 
$(OBJDIR)/trpel.o  : trpel.f90  trcomm.f90 
$(OBJDIR)/trpsc.o  : trpsc.f90  trcomm.f90 

$(OBJDIR)/trgout.o : trgout.f90 trcomm.f90 
$(OBJDIR)/trgrar.o : trgrar.f90 trcomm.f90 
$(OBJDIR)/trgrat.o : trgrat.f90 trcomm.f90 trcom1.f90
$(OBJDIR)/trgrap.o : trgrap.f90 trcomm.f90 
$(OBJDIR)/trgrae.o : trgrae.f90 trcomm.f90 
$(OBJDIR)/trgrad.o : trgrad.f90 trcomm.f90 
$(OBJDIR)/trgram.o : trgram.f90 trcomm.f90 trcom1.f90
$(OBJDIR)/trgsub.o : trgsub.f90 trcomm.f90 
$(OBJDIR)/trhelp.o : trhelp.f90 trcomm.f90 
$(OBJDIR)/trg2d.o  : trg2d.f90  trcomm.f90 
$(OBJDIR)/trg3d.o  : trg3d.f90  trcomm.f90 
$(OBJDIR)/tritg.o  : tritg.f90  trcomm.f90 
$(OBJDIR)/trmodels.o : trmodels.f90  trcomm.f90 
$(OBJDIR)/glf/callglf2d.o : glf/callglf2d.f glf/glf.m
	$(FCFIXED) $(FFLAGS) -c -o $(OBJDIR)/glf/callglf2d.o -Iglf glf/callglf2d.f
$(OBJDIR)/glf/glf2d.o  : glf/glf2d.F glf/glf.m glf/f77_dcomplx.h
	$(FCFIXED) $(FFLAGS) -c -o $(OBJDIR)/glf/glf2d.o -Iglf glf/glf2d.F
$(OBJDIR)/glf/r8tomsqz.o : glf/r8tomsqz.F glf/f77_dcomplx.h
	$(FCFIXED) $(FFLAGS) -c -o $(OBJDIR)/glf/r8tomsqz.o glf/r8tomsqz.F
$(OBJDIR)/glf/zgeev.o  : glf/zgeev.f
$(OBJDIR)/glf/blas_zgeev.o : glf/blas_zgeev.f
$(OBJDIR)/trncls.o : trncls.f90 trncls.inc trcomm.f90 nclass/pamx_mi.inc \
           nclass/pamx_ms.inc nclass/pamx_mz.inc
$(OBJDIR)/nclass/nclass_mod.o : nclass/nclass_mod.f nclass/pamx_mi.inc \
                      nclass/pamx_ms.inc nclass/pamx_mz.inc
	$(FCFIXED) $(FFLAGS) -c -o $(OBJDIR)/nclass/nclass_mod.o -Inclass nclass/nclass_mod.f $(MODDIR) $(MODINCLUDE)
$(OBJDIR)/itg/ifspppl_model.o : itg/ifspppl_model.f
$(OBJDIR)/itg/weiland.o: itg/weiland.f
$(OBJDIR)/nclass/rarray_sum.o : nclass/rarray_sum.f
$(OBJDIR)/nclass/rarray_copy.o: nclass/rarray_copy.f
$(OBJDIR)/nclass/rarray_zero.o: nclass/rarray_zero.f
$(OBJDIR)/nclass/u_erf.o      : nclass/u_erf.f
$(OBJDIR)/nclass/u_lu_backsub.o: nclass/u_lu_backsub.f
$(OBJDIR)/nclass/u_lu_decomp.o : nclass/u_lu_decomp.f
$(OBJDIR)/nclass/write_mod.o   : nclass/write_mod.f
$(OBJDIR)/cytran/cytran_mod.o: cytran/cytran_mod.f90
$(OBJDIR)/cytran/tr_cytran_mod.o: cytran/tr_cytran_mod.f90
$(OBJDIR)/mmm95/stripx.o: mmm95/stripx.f
$(OBJDIR)/mmm95/tomsqz.o: mmm95/tomsqz.f
$(OBJDIR)/mmm95/weiland14.: mmm95/weiland14.f
$(OBJDIR)/mmm95/mmm95.o: mmm95/mmm95.f
$(OBJDIR)/libmmm7_1/w20mod.o: libmmm7_1/w20mod.f90
$(OBJDIR)/libmmm7_1/modmmm7_1.o: libmmm7_1/modmmm7_1.f90 libmmm7_1/w20mod.f90
$(OBJDIR)/mbgb/mixed_Bohm_gyro_Bohm.o: mbgb/mixed_Bohm_gyro_Bohm.f
