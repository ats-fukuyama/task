### $Id: Makefile,v 1.69 2012/07/25 07:53:52 nuga Exp $ ###
include ../make.header
include ../mtxp/make.mtxp
#LIB_MTX=$(LIB_MTX_BND)
#LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_PCG)
LIB_MTX=$(LIB_MTX_MUMPS)
LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIB_MTX=$(LIBXMTX_KSP)
#LIBX_MTX=$(LIBX_MTX_KSP)

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

OBJDIR   = ./obj

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	@[ -d $(OBJDIR) ]
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../lib/$(MOD) -I../pl/$(MOD) -I../eq/$(MOD) -I../mtxp/$(MOD) -I../../bpsd/$(MOD) -I../ob/$(MOD)

SRCS = fpcomm.f90 eg_read.f90 fpinit.f90 fpparm.f90 \
		fpmpi.f90 fpbroadcast.f90 fpsub.f90 \
		fpexec.f90 fpsave.f90 fpwrite.f90 \
		fpcaltp.f90 fpcalte.f90 fpcaldeff.f90 fpcalchieff.f90 \
		fpoutdata.f90 fpwmin.f90 fpwrin.f90 fpcalwm.f90 fpcalwr.f90 fpcalw.f90 \
		fpcalcn.f90 fpcalcnr.f90 fpcalc.f90 fpnfrr.f90 fpnflg.f90 \
		cdbmfp.f90 fpcalr.f90 sv_read.f90 fpcoef.f90 fpcale.f90 \
		fpbounce.f90 fpdisrupt.f90 fpprep.f90 fploop.f90 \
		fpcont.f90 fpfout.f90 fpgout.f90 fpfile.f90\
		fowcomm.f90 foworbit.f90 fowdistribution.f90 orbit_classify.f90 fowprep.f90 \
		fowsource.f90 \
		fowcoef.f90 fowexec.f90 check_neoclass.f90 fowloop.f90\
		fowmenu.f90 \
		fpmenu.f90

SRCSFOW =fowcomm.f90 foworbit.f90 fowdistribution.f90 orbit_classify.f90 fowprep.f90 \
		fowsource.f90 \
		fowcoef.f90 fowexec.f90 check_neoclass.f90 fowloop.f90\
		fowmenu.f90

OBJS = $(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))
OBJSFOW = $(addprefix $(OBJDIR)/, $(SRCSFOW:.f90=.o))

LIBS = libfp.a ../pl/libpl.a ../eq/libeq.a  \
        ../lib/libgrf.a ../lib/libtask.a ../../bpsd/libbpsd.a ../ob/libob.a

LIBFOW = libfow.a

# all: libs fp
all: libs fow

libfp.a: $(OBJS)
	$(LD) $(LDFLAGS) libfp.a $(OBJS)
libfow.a: $(OBJSFOW)
	$(LD) $(LDFLAGS) libfow.a $(OBJSFOW)

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../ob/libob.a:
	(cd ../ob; make libob.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)
	(cd ../ob; make libob.a)

fp : $(LIB_MTX) $(LIBS) $(OBJDIR)/fpmain.o
	$(FLINKER) $(OBJDIR)/fpmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

fow : $(LIB_MTX) $(LIBS) $(LIBFOW) $(OBJDIR)/fowmain.o
	$(FLINKER) $(OBJDIR)/fowmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean :
	-rm -f core a.out *.o *.a ./#* ./*~ ./mod/*.mod
	-rm -f ./obj/*.o
	-rm -f fp fow out
	-rm -rf ./*.dSYM
	-rm -f ./csv/*
	-rm -rf __pycache__

veryclean: clean
	-rm -f fp

new:
	-mkdir ../fpnew
	cp -f Makefile ../fpnew
	cp -f *.f ../fpnew
	cp -f *.inc ../fpnew

$(OBJDIR)/fpmain.o : fpinit.f90 fpparm.f90 fpmain.f90 fpcomm.f90
$(OBJDIR)/fpmenu.o : fpwrin.f90 fpgout.f90 fpfout.f90 fploop.f90 fpprep.f90 fpparm.f90 fpfile.f90 fpmenu.f90 fpcomm.f90
$(OBJDIR)/fpinit.o : fpinit.f90 fpcomm.f90
$(OBJDIR)/fpparm.o : fpparm.f90 fpcomm.f90
$(OBJDIR)/fpprep.o : fpoutdata.f90 fpcale.f90 fpbounce.f90 fpnfrr.f90 fpwmin.f90 fpwrin.f90 fpsave.f90 fpcalw.f90 fpbroadcast.f90 fpparm.f90 fpcoef.f90 fpdisrupt.f90 eg_read.f90 sv_read.f90 fpprep.f90 fpcomm.f90
fploop.o : fpoutdata.f90 fpnfrr.f90 fpsave.f90 fpexec.f90 fpprep.f90 fpcoef.f90 fpdisrupt.f90 eg_read.f90 fploop.f90 fpcomm.f90 fpoutdata.f90
$(OBJDIR)/fpexec.o : fpexec.f90 fpcomm.f90
$(OBJDIR)/fpdrexec.o : fpdrexec.f90 fpcomm.f90
$(OBJDIR)/fpbroadcast.o : fpbroadcast.f90 fpcomm.f90
$(OBJDIR)/fpcoef.o : fpcalr.f90 fpnfrr.f90 fpcalwm.f90 fpcalwr.f90 fpcalw.f90 fpcalc.f90 fpcoef.f90 fpcomm.f90
$(OBJDIR)/fpcalc.o : fpcalcnr.f90 fpcalcn.f90 fpcalc.f90 fpcomm.f90
$(OBJDIR)/fpcalcn.o : fpcalcn.f90 fpcomm.f90
$(OBJDIR)/fpcalcnr.o : fpcoef.f90 fpcalcnr.f90 fpcomm.f90
$(OBJDIR)/fpcalw.o : fpcalw.f90 fpcalwr.f90 fpcalwm.f90 fpcomm.f90
$(OBJDIR)/fpcalwr.o : fpwrin.f90 fpcalwr.f90 fpcomm.f90
$(OBJDIR)/fpcalwm.o : fpwmin.f90 fpcalwm.f90 fpcomm.f90
$(OBJDIR)/fpcont.o : fpcont.f90 fpcomm.f90
$(OBJDIR)/fpfout.o : fpfout.f90 fpcomm.f90
$(OBJDIR)/fpgout.o : fpfout.f90 fpcont.f90 fpgout.f90 fpcomm.f90
$(OBJDIR)/fpfile.o : fpprep.f90 eg_read.f90 sv_read.f90 fpfile.f90 fpcomm.f90
$(OBJDIR)/fpsave.o : fpexec.f90 eg_read.f90 fpsave.f90 fpcomm.f90
$(OBJDIR)/fpwrin.o : fpwrin.f90 fpcomm.f90
$(OBJDIR)/fpwmin.o : fpwrin.f90 fpcomm.f90
$(OBJDIR)/fpnfrr.o : fpparm.f90 fpnfrr.f90 fpcomm.f90
$(OBJDIR)/fpbounce.o : fpparm.f90 fpbounce.f90 fpcomm.f90
$(OBJDIR)/fpcale.o: fpcale.f90 fpmpi.f90
$(OBJDIR)/fpdisrupt.o: fpcale.f90 fpcoef.f90 fpdisrupt.f90 fpcomm.f90
$(OBJDIR)/fpcalr.o: fpcalr.f90 fpcomm.f90
$(OBJDIR)/eg_read.o: eg_read.f90 fpcomm.f90
$(OBJDIR)/sv_read.o: sv_read.f90 fpcomm.f90
$(OBJDIR)/fpoutdata.o: fpoutdata.f90 eg_read.f90 fpcomm.f90
$(OBJDIR)/fpcaltp.o: fpcaltp.f90 fpcomm.f90
$(OBJDIR)/fpcaldeff.o: fpcaldeff.f90 fpcomm.f90
$(OBJDIR)/fpwrite.o : fpwrite.f90 fpcomm.f90

$(OBJDIR)/fowmain.o : fowmain.f90 fowcomm.f90 fowprep.f90 foworbit.f90 \
            fpprep.f90 fpinit.f90 fpcomm.f90 fpwrite.f90
$(OBJDIR)/fowmenu.o : fowmain.f90 fowcomm.f90 fowprep.f90 foworbit.f90 \
            fpprep.f90 fpcomm.f90 fpwrite.f90
$(OBJDIR)/fowcomm.o : fowcomm.f90
$(OBJDIR)/fowexec.o : fowexec.f90 fowcomm.f90 fpcomm.f90
$(OBJDIR)/fowprep.o : fowprep.f90 fowcomm.f90 fowdistribution.f90 foworbit.f90 fowdistribution.f90 orbit_classify.f90
$(OBJDIR)/foworbit.o: foworbit.f90 fowcomm.f90 fpcomm.f90 \
											../ob/obcomm.f90 ../ob/obcalc.f90 ../ob/obprep.f90 ../ob/obinit.f90 \
											../ob/obparm.f90
$(OBJDIR)/fowdistribution.o : fowdistribution.f90 foworbit.f90 fowcomm.f90 fpcomm.f90
$(OBJDIR)/fowcoef.o : fowdistribution.f90 fowcoef.f90 fowcomm.f90 fpcomm.f90
$(OBJDIR)/fowsource.o : fowsource.f90 fowcomm.f90 fpcomm.f90 fowprep.f90 foworbit.f90 fowdistribution.f90
$(OBJDIR)/fowloop.o : fowexec.f90 fowdistribution.f90 fowcoef.f90 fowsource.f90 \
											fowcomm.f90 fpcomm.f90 fowprep.f90 fpwrite.f90 check_neoclass.f90 orbit_classify.f90
$(OBJDIR)/check_neoclass.o : check_neoclass.f90 fowcomm.f90 fpcomm.f90 fpwrite.f90 fowcoef.f90 fowexec.f90
$(OBJDIR)/orbit_classify.o : orbit_classify.f90 fowcomm.f90 fpcomm.f90 fpwrite.f90
