###### $Id$ ###
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

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I$(MOD) -I../open-adas/adf11/adf11-lib/$(MOD) \
            -I../adpost/$(MOD) -I../eq/$(MOD) -I../pl/$(MOD) \
            -I../dp/$(MOD) -I../wr/$(MOD) -I../wm/$(MOD) \
            -I../fp/$(MOD) -I../tr/$(MOD) -I../ti/$(MOD) \
            -I../mtxp/$(MOD) -I../lib/$(MOD) \
            -I../../bpsd/$(MOD) 

EQLIBS= ../pl/noequlib.o ../eq/libeq.a 

SRCS = totmenu.f90

OBJS = $(SRCS:.f90=.o)

LIBS = ../lib/libmds.a ../lib/libgrf.a ../lib/libtask.a ../../bpsd/libbpsd.a \
       ../adpost/lib-adpost.a ../open-adas/adf11/adf11-lib/lib-adf11.a \
       ../pl/libpl.a \
       $(EQLIBS) \
       ../tr/libtr.a ../dp/libdp.a ../wr/libwr.a ../wm/libwm.a \
       ../fp/libfp.a ../ti/libti.a


all : libs tot

../open-adas/adf11/adf11-lib/lib-adf11.a:
	(cd ../open-adas/adf11/adf11-lib; make lib-adf11.a)
../adpost/lib-adpost.a:
	(cd ../adpost; make lib-adpost.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libmds.a:
	(cd ../lib; make libmds.a)
$(LIB_MTX):
	(cd ../mtxp; make)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../pl/noeqlib.o:
	(cd ../pl; make noeqlib.o)
../pl/noequlib.o:
	(cd ../pl; make noequlib.o)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../tr/libtr.a:
	(cd ../tr; make libtr.a)
../dp/libdp.a:
	(cd ../dp; make libdp.a)
../wr/libwr.a:
	(cd ../wr; make libwr.a)
../wm/libwm.a:
	(cd ../wm; make libwm.a)
../fp/libfp.a:
	(cd ../fp; make libfp.a)
../ti/libti.a:
	(cd ../ti; make libti.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libmds.a)
	(cd ../mtxp; make)
	(cd ../adpost; make lib-adpost.a)
	(cd ../open-adas/adf11/adf11-lib; make lib-adf11.a)
	(cd ../pl; make libpl.a)
	(cd ../pl; make noeqlib.o)
	(cd ../pl; make noequlib.o)
	(cd ../eq; make libeq.a)
	(cd ../tr; make libtr.a)
	(cd ../dp; make libdp.a)
	(cd ../wr; make libwr.a)
	(cd ../wm; make libwm.a)
	(cd ../fp; make libfp.a)
	(cd ../ti; make libti.a)
	make totmenu.o

tot : $(LIB_MTX) $(LIBS) totmenu.o totmain.o
	$(FLINKER) totmenu.o totmain.o -o $@ $(FFLAGS) $(LIBS) $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out *.o *.mod ./*~ *.a $(MOD)/*.mod

veryclean: clean
	-rm -f tot

totmain.o: totmain.f90
totmenu.o: totmenu.f90

