### $Id$ ###
include ../make.header
include ../mtxp/make.mtxp

LIB_MTX=$(LIB_MTX_BND)
LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_PCG)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIB_MTX=$(LIBXMTX_KSP)
#LIBX_MTX=$(LIBX_MTX_KSP)

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../open-adas/adf11/adf11-lib/$(MOD) \
            -I../adpost/$(MOD) -I../eq/$(MOD) -I../plx/$(MOD) \
            -I../mtxp/$(MOD) -I../lib/$(MOD) -I../../bpsd/$(MOD) 

SRCS=   ticomm.f90 tiadas.f90 \
	tiinit.f90 tiparm.f90 \
	ticoef.f90 tisource.f90 ticalc.f90 tirecord.f90 \
        tiprep.f90 tiexec.f90 \
        tigout.f90 timenu.f90

#      trprof.f90 trufsub.f90 trufile.f90 trfile.f90 trhelp.f90 \
#      trexec.f90 trcalc.f90 \
#      tradat.f90 cdbm.f90 trmodels.f90 trcoef.f90 tritg.f90 \
#      trrslt.f90 \
#      trpnb.f90 trprf.f90 trpnf.f90 trpel.f90 trpsc.f90 trmdlt.f90 \
#      trgout.f90 trgrar.f90 trgrat.f90 trgrap.f90 trgrae.f90 \
#      trgrad.f90 trgram.f90 trgsub.f90 trfout.f90 \
#      trunit.f90 trloop.f90 trmenu.f90 tr_bpsd.f90

OBJS=$(SRCS:.f90=.o)

LIBS=libti.a ../adpost/lib-adpost.a ../open-adas/adf11/adf11-lib/lib-adf11.a ../eq/eqlib.a ../plx/pllib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a


all: libs ti

libti.a: $(OBJS)
	$(LD) $(LDFLAGS) $@ $(OBJS)

../open-adas/adf11/adf11-lib/lib-adf11.a: 
	(cd ../open-adas/adf11/adf11-lib; make lib-adf11.a)
../adpost/lib-adpost.a: 
	(cd ../adpost; make lib-adpost.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
$(LIB_MTX):
	(cd ../mtxp; make)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)


libs:
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../lib; make tasklib.a)
	(cd ../mtxp; make)
	(cd ../open-adas/adf11/adf11-lib; make lib-adf11.a)
	(cd ../adpost; make lib-adpost.a)
	(cd ../plx; make pllib.a)
	(cd ../eq; make eqlib.a)

ti: $(LIB_MTX) $(LIBS) timain.o 
	$(FLINKER) timain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean : 
	rm -f ./#* ./*~ *.o nclass/*.o itg/*.o glf/*.o cytran/*.o mmm95/*.o mbgb/*.o a.out core *.a *.mod $(MOD)/*.mod

veryclean : clean
	rm -f ti

distclean: clean
	cd ../tr && make clean
	cd ../plx && make clean
	cd ../eq && make clean
	cd ../lib && make clean

ticomm.o : ticomm.f90
timain.o : timain.f90 ticomm.f90  timenu.f90
timenu.o : timenu.f90 ticomm.f90
tiinit.o : tiinit.f90 ticomm.f90
tiparm.o : tiparm.f90 ticomm.f90
tiprep.o : tiprep.f90 ticomm.f90 ../adpost/lib-adpost.f90 \
           ../open-adas/adf11/adf11-lib/lib-adf11.f90
tigout.o : tigout.f90 ticomm.f90
tiexec.o : tiexec.f90 ticomm.f90
ticoef.o : ticoef.f90 ticomm.f90
tisource.o : tisource.f90 ticomm.f90
ticalc.o : ticalc.f90 ticomm.f90
tirecord.o : tirecord.f90 ticomm.f90

nclass/rarray_sum.o : nclass/rarray_sum.f
nclass/rarray_copy.o: nclass/rarray_copy.f
nclass/rarray_zero.o: nclass/rarray_zero.f
nclass/u_erf.o      : nclass/u_erf.f
nclass/u_lu_backsub.o: nclass/u_lu_backsub.f
nclass/u_lu_decomp.o : nclass/u_lu_decomp.f
nclass/write_mod.o   : nclass/write_mod.f
