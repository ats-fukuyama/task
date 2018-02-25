### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE= -I./$(MOD) -I../../bpsd/$(MOD) -I../eq/$(MOD) -I../plx/$(MOD) -I../lib/$(MOD)

BPSDLIB=../../bpsd/bpsdlib.a
PLLIB=../plx/pllib.a 
EQLIB=../eq/eqlib.a

SRCS= ticomm.f90 libmtxbnd.f90 adpost/lib-adpost.f90 \
	tiinit.f90 tiparm.f90 tiprep.f90 \
	ticoef.f90 tisource.f90 ticalc.f90 tiloop.f90 tigout.f90 timenu.f90

#      trprof.f90 trufsub.f90 trufile.f90 trfile.f90 trhelp.f90 \
#      trexec.f90 trcalc.f90 \
#      tradat.f90 cdbm.f90 trmodels.f90 trcoef.f90 tritg.f90 \
#      trrslt.f90 \
#      trpnb.f90 trprf.f90 trpnf.f90 trpel.f90 trpsc.f90 trmdlt.f90 \
#      trgout.f90 trgrar.f90 trgrat.f90 trgrap.f90 trgrae.f90 \
#      trgrad.f90 trgram.f90 trgsub.f90 trfout.f90 \
#      trunit.f90 trloop.f90 trmenu.f90 tr_bpsd.f90

OBJS=$(SRCS:.f90=.o)
LIBS=tilib.a $(EQLIB) $(PLLIB) ../lib/tasklib.a $(BPSDLIB)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

#all: libs ti
all: ti

tilib.a: $(OBJS)
	$(LD) $(LDFLAGS) tilib.a $(OBJS)

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../lib/mdslib.a:
	(cd ../lib; make mdslib.a)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
	(cd ../tr; make tr3lib.a)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../plx; make pllib.a)
	(cd ../eq; make eqlib.a)

ti: $(LIBS) timain.o 
	$(FCFREE) timain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MODINCLUDE)

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
tiprep.o : tiprep.f90 ticomm.f90
tigout.o : tigout.f90 ticomm.f90
tiloop.o : tiloop.f90 ticomm.f90
ticoef.o : ticoef.f90 ticomm.f90
tisource.o : tisource.f90 ticomm.f90
ticalc.o : ticalc.f90 ticomm.f90

adpost/lib-adpost.o:adpost/lib-adpost.f90 

nclass/rarray_sum.o : nclass/rarray_sum.f
nclass/rarray_copy.o: nclass/rarray_copy.f
nclass/rarray_zero.o: nclass/rarray_zero.f
nclass/u_erf.o      : nclass/u_erf.f
nclass/u_lu_backsub.o: nclass/u_lu_backsub.f
nclass/u_lu_decomp.o : nclass/u_lu_decomp.f
nclass/write_mod.o   : nclass/write_mod.f
