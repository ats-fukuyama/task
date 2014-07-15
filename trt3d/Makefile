include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .c .mod .o .a

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)
CFLAGS = $(OCFLAGS)

# without eq
MODINCLUDE= -I.$(MOD) -I../../bpsd/$(MOD) -I../pl/$(MOD)
EQLIBS= ../pl/noeqlib.o ../pl/noequlib.o 
MAKE_TPXLIB=

# with task/eq
#MODINCLUDE= -I.$(MOD) -I../../bpsd$(MOD) -I../eq$(MOD) -I../pl$(MOD)
#EQLIBS= ../eq/eqlib.a ../pl/noequlib.o 
#MAKE_TPXLIB=

# with topics/equ
#MODINCLUDE= -I.$(MOD) -I../../bpsd$(MOD) -I../pl$(MOD) -I../../topics/src$(MOD)
#EQLIBS= ../../topics/src/tpxlib.a ../../topics/src/eqlib.a  ../pl/noeqlib.o 
#MAKE_TPXLIB= (cd ../../topics/src; make tpxlib.a eqlib.a)

# with task/eq and topics/equ
#MODINCLUDE= -I.$(MOD) -I../../bpsd$(MOD) -I../eq$(MOD) -I../pl$(MOD) -I../../topics/src$(MOD)
#EQLIBS= ../eq/eqlib.a ../../topics/src/tpxlib.a ../../topics/src/eqlib.a
#MAKE_TPXLIB= (cd ../../topics/src; make tpxlib.a eqlib.a)

PLLIBS=../../bpsd/bpsdlib.a ../pl/pl_vmec.o ../pl/modules.o

SRCC=t3d_NNW_DATA_custum_20121019.c t3d_intrfcDCOMNNW.c t3d_TR_FIT3D_CONNECT.c
#SRC03=t3d_FILE_io.f90 t3d_powerbalance.f90
SRC90=trcom0.f90 trcom1.f90 trcom2.f90 trcom3.f90 trcomm.f90 tr_bpsd.f90\
      trinit.f90 t3d_trprof.f90 t3d_trfile_20090422.f90 trhelp.f90 t3d_nogsaf.f90 \
	  t3d_nccoef_dgn_lhd.f90 \
	  t3d_er_param.f90 \
 	  t3d_tentative_flags.f90 \
          t3d_FILE_io.f90 \
          t3d_powerbalance.f90 \
	  t3d_fit3d_exec.f90 \
	  t3d_nccoef_shaing.f90 \
	  t3d_trexec.f90 \
	  t3d_er.f90  \
      t3d_trcalc.f90 \
      t3d_trcoef.f90 tritg.f90 trrslt.f90 \
      trpnb.f90 trprf.f90 trpnf.f90 trpel.f90 trmdlt.f90 \
      trgout.f90 trgrar.f90 trgrat.f90 trgrap.f90 trgrae.f90 \
      trgrad.f90 trgram.f90 trgsub.f90 trunit.f90 \
      t3d_vmec2tr.f90 t3d_input_prof.f90  \
	  t3d_trmenu.f90 \
	  t3d_trloop.f90 \
  	  depsum_100208.f90 \
	  BOZdata_translation.f90 \
	  t3d_parasurvey_loop.f90 \
	  fit_read20_100208.f90

#SRCITG= glf/callglf2d.f glf/zgeev.f glf/blas_zgeev.f itg/ifspppl_model.f itg/weiland.f
#SRCITGF= glf/glf2d.F glf/r8tomsqz.F
SRCITG= noitg.f 
SRCITGF=
SRC2D=trg2d.f90
SRC3D=trg3d.f90

OBJS= $(SRC03:.f90=.o) $(SRC90:.f90=.o) $(SRCITG:.f=.o) $(SRCITGF:.F=.o) $(SRCC:.c=.o)
OBJ2D=$(SRC2D:.f90=.o)
OBJ3D=$(SRC3D:.f90=.o)

LIBS2=$(PLLIBS) $(EQLIBS) tr2lib.a ../lib/tasklib.a ../lib/mdslib.a
LIBS3=$(PLLIBS) $(EQLIBS) tr3lib.a ../lib/tasklib.a ../lib/mdslib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -fpp -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE)  $(FFLAGS) -fpp -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.c.o :
	$(CC) $(CFLAGS) -c $< -o $@

#all: tr2 tr3
all: trt3d

tr2lib.a: $(OBJS) $(OBJ2D) $(OBJE) $(OBJL) $(OBJNC) $(INCNC) $(INCGLF)
	$(LD) $(LDFLAGS) tr2lib.a $(OBJS) $(OBJ2D) $(OBJE) $(OBJL) $(OBJNC)

tr3lib.a: $(OBJS) $(OBJ3D) $(OBJE) $(OBJL) $(OBJNC) $(INCNC) $(INCGLF)
	$(LD) $(LDFLAGS) tr3lib.a $(OBJS) $(OBJ3D) $(OBJE) $(OBJL) $(OBJNC)

../pl/pllib.a:
	(cd ../pl; make pllib.a)
../pl/noeqlib.o:
	(cd ../pl; make noeqlib.o)
../pl/noequlib.o:
	(cd ../pl; make noequlib.o)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../eq/treqin.o:
	(cd ../eq; make treqin.o)
../lib/tasklib.a:
	(cd ../lib; make)
../lib/mdslib.a:
	(cd ../lib; make)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../pl; make pllib.a)
	(cd ../eq; make eqlib.a)
	$(MAKE_TPXLIB)
	(cd ../eq; make treqin.o)

trt3d: $(LIBS2) t3d_trmain.o 
	$(FCFREE) t3d_trmain.o $(LIBS2) -o $@ $(FFLAGS) $(GFLIBS) $(MDSLIB) $(LIBLA) $(MODULE) $(MODO) $(MODINCLUDE)

tr2: $(LIBS2) t3d_trmain.o 
	$(FCFREE) t3d_trmain.o $(LIBS2) -o $@ $(FFLAGS) $(GFLIBS) $(MDSLIB) $(LIBLA) $(MODULE) $(MODO) $(MODINCLUDE)

tr3: $(LIBS3) t3d_trmain.o 
	$(FCFREE) t3d_trmain.o $(LIBS3) -o $@ $(FFLAGS) $(GLLIBS) $(MDSLIB) $(LIBLA) $(MODINCLUDE)

check : $(SRCS)
	ftnchek $(FCKFLAGS) $(SRCS) $(SRCE) $(SRCL) $(SRC2D) $(SRCNC) $(FCKPIPE) | less

clean : 
	rm -f ./#* ./*~ *.o nclass/*.o itg/*.o glf/*.o a.out core *.a *.mod .$(MOD)/*.mod

veryclean : clean
	rm -f tr2 tr3

distclean: clean
	cd ../pl && make clean
	cd ../eq && make clean
	cd ../lib && make clean

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

# as in the past 
trpl.o : trpl.f90 trcomm.f90 trcom1.f90
trrslt.o : trrslt.f90 trcomm.f90 trcom2.f90
trpnb.o  : trpnb.f90  trcomm.f90 
trprf.o  : trprf.f90  trcomm.f90 
trpnf.o  : trpnf.f90  trcomm.f90 
trpel.o  : trpel.f90  trcomm.f90 
trinit.o : trinit.f90 trcomm.f90 trcom1.f90 trcom3.f90
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
trmdlt.o : trmdlt.f90 trcomm.f90
depsum_100208.o :  depsum_100208.f90
# newly - TASK3D -
nomdsplus.o : nomdsplus.f
t3d_trmenu.o : t3d_trmenu.f90 trcomm.f90 trcom1.f90 t3d_FILE_io.f90 t3d_er_param.f90
t3d_nccoef_shaing.o : t3d_nccoef_shaing.f90 trcomm.f90 t3d_nccoef_dgn_lhd.f90
t3d_nccoef_dgn_lhd.o: t3d_nccoef_dgn_lhd.f90 t3d_nccoef_shaing.f90 trcomm.f90
t3d_FILE_io.o : t3d_FILE_io.f90 trcomm.f90 t3d_tentative_flags.f90
t3d_er_param.o : t3d_er_param.f90 trcomm.f90
t3d_tentative_flags.o : t3d_tentative_flags.f90
t3d_parasurvey_loop.o : t3d_parasurvey_loop.f90 trcomm.f90 trcom1.f90 t3d_FILE_io.f90 t3d_er_param.f90 t3d_tentative_flags.f90
t3d_fit3d_exec.o : t3d_fit3d_exec.f90 trcomm.f90 t3d_tentative_flags.f90
t3d_powerbalance.o : t3d_powerbalance.f90 trcomm.f90 t3d_tentative_flags.f90 t3d_FILE_io.f90 t3d_nccoef_dgn_lhd.f90
t3d_NNW_DATA_custum_20121019.o : t3d_NNW_DATA_custum_20121019.c
t3d_intrfcDCOMNNW.o : t3d_intrfcDCOMNNW.c
t3d_trmain.o : t3d_trmain.f90 trcomm.f90  
t3d_trcalc.o : t3d_trcalc.f90 trcomm.f90 trcom1.f90 t3d_er.f90 t3d_nccoef_dgn_lhd.f90 t3d_tentative_flags.f90 t3d_FILE_io.f90 t3d_er_param.f90
t3d_trcoef.o : t3d_trcoef.f90 trcomm.f90 trglf.inc t3d_er.f90 t3d_tentative_flags.f90 t3d_FILE_io.f90
t3d_trprof.o : t3d_trprof.f90 trcomm.f90 trcom1.f90
t3d_trfile_20090422.o : t3d_trfile_20090422.f90 trcomm.f90 trcom1.f90
t3d_trloop.o : t3d_trloop.f90 trcomm.f90 trcom1.f90 t3d_er_param.f90 t3d_tentative_flags.f90 t3d_tentative_flags.f90 t3d_FILE_io.f90 t3d_er.f90
t3d_trexec.o : t3d_trexec.f90 trcomm.f90 trcom1.f90 t3d_FILE_io.f90 t3d_er_param.f90 t3d_tentative_flags.f90
t3d_vmec2tr.o : t3d_vmec2tr.f90 trcomm.f90
t3d_input_prof.o : t3d_input_prof.f90 trcomm.f90
t3d_er.o : t3d_er.f90 trcomm.f90 t3d_nccoef_dgn_lhd.f90 t3d_nccoef_shaing.f90 t3d_FILE_io.f90 t3d_er_param.f90
# newly - FIT3D -
t3d_TR_FIT3D_CONNECT.o : t3d_TR_FIT3D_CONNECT.c
fit_read20_100208.o : fit_read20_100208.f90
depsum_100208.o : depsum_100208.f90
BOZdata_translation.o : BOZdata_translation.f90 t3d_tentative_flags.f90



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
