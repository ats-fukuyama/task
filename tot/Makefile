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

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I$(MOD) -I../open-adas/adf11/adf11-lib/$(MOD) \
            -I../adpost/$(MOD) -I../eq/$(MOD) -I../plx/$(MOD) \
            -I../fp/$(MOD) -I../ti/$(MOD) -I../mtxp/$(MOD) -I../lib/$(MOD) \
            -I../../bpsd/$(MOD) 

EQLIBS= ../plx/noequlib.o ../eq/eqlib.a 

SRCS = totmenu.f90

OBJS = $(SRCS:.f90=.o)

LIBS = ../lib/mdslib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a \
       ../adpost/lib-adpost.a ../open-adas/adf11/adf11-lib/lib-adf11.a \
       ../plx/pllib.a \
       ../tr/trcom2.o ../tr/trcom3.o $(EQLIBS) \
       ../tr/tr2lib.a ../dp/dplib.a ../wr/wrlib.a ../wm/wmlib.a \
       ../fp/fplib.a ../ti/libti.a $(OBJS)


all : libs tot

../open-adas/adf11/adf11-lib/lib-adf11.a:
	(cd ../open-adas/adf11/adf11-lib; make lib-adf11.a)
../adpost/lib-adpost.a:
	(cd ../adpost; make lib-adpost.a)
../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../lib/mdslib.a:
	(cd ../lib; make mdslib.a)
$(LIB_MTX):
	(cd ../mtxp; make)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
../plx/pllib.a:
	(cd ../plx; make pllib.a)
../plx/noeqlib.o:
	(cd ../plx; make noeqlib.o)
../plx/noequlib.o:
	(cd ../plx; make noequlib.o)
../eq/eqlib.a:
	(cd ../eq; make eqlib.a)
../tr/trcom2.o:
	(cd ../tr; make trcom2.o)
../tr/trcom3.o:
	(cd ../tr; make trcom3.o)
../tr/tr2lib.a:
	(cd ../tr; make tr2lib.a)
../dp/dplib.a:
	(cd ../dp; make dplib.a)
../wr/wrlib.a:
	(cd ../wr; make wrlib.a)
../wm/wmlib.a:
	(cd ../wm; make wmlib.a)
../wm/mpiwmlib.a:
	(cd ../wm; make mpiwmlib.a)
../fp/fplib.a:
	(cd ../fp; make fplib.a)
../ti/libti.a:
	(cd ../ti; make libti.a)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../lib; make mdslib.a)
	(cd ../../bpsd; make bpsdlib.a)
	(cd ../mtxp; make)
	(cd ../plx; make pllib.a)
	(cd ../plx; make noeqlib.o)
	(cd ../plx; make noequlib.o)
	(cd ../eq; make eqlib.a)
	(cd ../tr; make trcom2.o)
	(cd ../tr; make trcom3.o)
	(cd ../tr; make tr2lib.a)
	(cd ../dp; make dplib.a)
	(cd ../wr; make wrlib.a)
	(cd ../wm; make wmlib.a)
	(cd ../fp; make fplib.a)
	(cd ../ti; make libti.a)
	make totmenu.o

tot : $(LIB_MTX) $(LIBS) totmain.o
	$(FLINKER) totmain.o -o $@ $(FFLAGS) $(LIBS) $(FLIBS) $(LIBX_MTX)

mpitot : $(LIBM) totmain.o
	$(FLINKER) totmain.o -o $@ $(FFLAGS) $(LIBS) $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out *.o *.mod ./*~ *.a $(MOD)/*.mod

veryclean: clean
	-rm -f tot

totmain.o: totmain.f90
totmenu.o: totmenu.f90

