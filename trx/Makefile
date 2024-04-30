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

MODINCLUDE= -I./$(MOD) -I../trmlib/$(MOD) -I../pl/$(MOD) -I../eq/$(MOD) \
            -I../equ/$(MOD) -I../lib/$(MOD) -I../mtxp/$(MOD) \
            -I../../bpsd/$(MOD)

${OBJDIR}/%.o:%.f
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

SRCM=trcomm.f90 trbpsd.f90
SRCS=trinit.f90 trparm.f90 trview.f90 \
     trmetric.f90 trfixed.f90  trprof.f90 trprep.f90 \
     trncls.f90 trcytran.f90 \
     trfile.f90 trhelp.f90 \
     trpnf.f90 \
     trexec.f90 trcalc.f90 \
     tradat.f90 trmodels.f90 trcoef.f90 tritg.f90 \
     trrslt.f90 trpnb.f90 trprf.f90 trpel.f90 trpsc.f90 \
     trgout.f90 trgrar.f90 trgrat.f90 trgrap.f90 trgrae.f90 \
     trgrad.f90 trgram.f90 trgsub.f90 trfout.f90 \
     trloop.f90 trmenu.f90

SRC2D=trg2d.f90
SRC3D=trg3d.f90

OBJM=$(addprefix $(OBJDIR)/, $(SRCM:.f90=.o))
OBJS=$(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))
OBJ2D=$(addprefix $(OBJDIR)/, $(SRC2D:.f90=.o))
OBJ3D=$(addprefix $(OBJDIR)/, $(SRC3D:.f90=.o))

OBJTR=$(OBJM) $(OBJS)

LIBS=../trmlib/libtrm.a ../eq/libeq.a ../pl/libpl.a ../lib/libmds.a \
     ../lib/libtask.a ../lib/libgrf.a  ../../bpsd/libbpsd.a
LIBS2=libtr2.a $(LIBS)
LIBS3=libtr3.a $(LIBS)

#all: libs tr2 tr3
all: libs tr2

libtr2.a: $(OBJTR) $(OBJ2D)
	$(LD) $(LDFLAGS) libtr2.a $(OBJTR) $(OBJ2D)

libtr3.a: $(OBJTR) $(OBJ3D)
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
../trmlib/libtrm.a:
	(cd ../trmlib; make libtrm.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)
	(cd ../trmlib; make libtrm.a)

tr2: $(LIBS2) $(OBJDIR)/trmain.o 
	$(FLINKER) $(OBJDIR)/trmain.o $(LIBS2) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

tr3: $(LIBS3) $(OBJDIR)/trmain.o 
	$(FLINKER) $(OBJDIR)/trmain.o $(LIBS3) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

trplot: $(OBJDIR)/trcomm.o $(OBJDIR)/trfixed.o $(OBJDIR)/trplot.o 
	$(FLINKER) $(OBJDIR)/trcomm.o $(OBJDIR)/trfixed.o $(OBJDIR)/trplot.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean : 
	-rm -f ./#* ./*~ $(OBJDIR)/*.o $(OBJDIR)/*/*.o $(MOD)/*.mod $(MOD)/*/*.mod core *.a

veryclean : clean
	-rm -f tr2 tr3

distclean: clean
	cd ../pl && make clean
	cd ../eq && make clean
	cd ../lib && make clean

$(OBJDIR)/trcomm.o : trcomm.f90  
$(OBJDIR)/trmain.o : trmain.f90 trcomm.f90  
$(OBJDIR)/trmenu.o : trmenu.f90 trcomm.f90
$(OBJDIR)/trinit.o : trinit.f90 trcomm.f90
$(OBJDIR)/trparm.o : trparm.f90 trcomm.f90
$(OBJDIR)/trview.o : trview.f90 trcomm.f90
$(OBJDIR)/trfixed.o : trfixed.f90 trcomm.f90
$(OBJDIR)/trprep.o : trprep.f90 trcomm.f90
$(OBJDIR)/trprof.o : trprof.f90 trcomm.f90
$(OBJDIR)/trmetric.o : trmetric.f90 trcomm.f90
$(OBJDIR)/trfile.o : trfile.f90 trcomm.f90
$(OBJDIR)/trbpsd.o : trbpsd.f90 trcomm.f90
$(OBJDIR)/tradat.o : tradat.f90 trcomm.f90
$(OBJDIR)/trufile.o : trufile.f90 trcomm.f90
$(OBJDIR)/tr_ufile_task.o : tr_ufile_task.f90 trcomm.f90
$(OBJDIR)/tr_ufile_topics.o : tr_ufile_topics.f90 trcomm.f90
$(OBJDIR)/trufsub.o : trufsub.f90 trcomm.f90
$(OBJDIR)/trloop.o : trloop.f90 trcomm.f90
$(OBJDIR)/trexec.o : trexec.f90 trcomm.f90 trfixed.f90
$(OBJDIR)/trpl.o : trpl.f90 trcomm.f90
$(OBJDIR)/trrslt.o : trrslt.f90 trcomm.f90
$(OBJDIR)/trcalc.o : trcalc.f90 trcomm.f90
$(OBJDIR)/trcoef.o : trcoef.f90 trcomm.f90 trglf.inc
$(OBJDIR)/trpnb.o  : trpnb.f90  trcomm.f90 
$(OBJDIR)/trprf.o  : trprf.f90  trcomm.f90 
$(OBJDIR)/trpnf.o  : trpnf.f90  trcomm.f90 
$(OBJDIR)/trpel.o  : trpel.f90  trcomm.f90 
$(OBJDIR)/trpsc.o  : trpsc.f90  trcomm.f90 

$(OBJDIR)/trgout.o : trgout.f90 trcomm.f90 
$(OBJDIR)/trgrar.o : trgrar.f90 trcomm.f90 
$(OBJDIR)/trgrat.o : trgrat.f90 trcomm.f90
$(OBJDIR)/trgrap.o : trgrap.f90 trcomm.f90 
$(OBJDIR)/trgrae.o : trgrae.f90 trcomm.f90 
$(OBJDIR)/trgrad.o : trgrad.f90 trcomm.f90 
$(OBJDIR)/trgram.o : trgram.f90 trcomm.f90
$(OBJDIR)/trgsub.o : trgsub.f90 trcomm.f90 
$(OBJDIR)/trhelp.o : trhelp.f90 trcomm.f90 
$(OBJDIR)/trg2d.o  : trg2d.f90  trcomm.f90 
$(OBJDIR)/trg3d.o  : trg3d.f90  trcomm.f90 
$(OBJDIR)/tritg.o  : tritg.f90  trcomm.f90 
$(OBJDIR)/trmodels.o : trmodels.f90  trcomm.f90 
$(OBJDIR)/trncls.o : trncls.f90 trcomm.f90 ../trmlib/nclass/pamx_mi.inc \
           ../trmlib/nclass/pamx_ms.inc ../trmlib/nclass/pamx_mz.inc
$(OBJDIR)/trcytran.o: trcytran.f90

$(OBJDIR)/trplot.o  : trplot.f90 trfixed.f90 trcomm.f90
