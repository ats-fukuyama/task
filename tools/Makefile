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

#FFLAGS = $(DFLAGS)
FFLAGS = $(OFLAGS)

URSRCS=ufread.f ufsub.f
UR3SRCS=ufread3.f ufsub.f
GUISRCS=guiread.f ufsub.f
SPLSRC=usplit.f
UROBJS=$(URSRCS:.f=.o)
UR3OBJS=$(UR3SRCS:.f=.o)
GUIOBJS=$(GUISRCS:.f=.o)
SPLOBJ=$(SPLSRC:.f=.o)

LIBS=../../bpsd/libbpsd.a ../lib/libgrf.a ../lib/libtask.a ../lib/libmds.a ../tr/libtr.a ../eq/libeq.a ../pl/libpl.a 
MODINCLUDE=-I$(MOD) -I../../bpsd/$(MOD) -I../lib/$(MOD) -I../pl/$(MOD) -I../tr/$(MOD)


.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : ufread ufread3 guiread usplit
#all : ufread usplit

../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../tr/libtr.a:
	(cd ../tr; make libtr.a)

ufread : $(UROBJS) $(LIBS)
	-$(FLINKER) -o $@ $(LIBS) $(FFLAGS) $(UROBJS) $(FLIBS) $(LIBX_MTX)

ufread3 : $(UR3OBJS) $(LIBS)
	-$(FLINKER) -o $@ $(LIBS) $(FFLAGS) $(UR3OBJS) $(FLIBS) $(LIBX_MTX)

guiread : $(GUIOBJS) $(LIBS)
	-$(FLINKER) -o $@ $(LIBS) $(FFLAGS) $(GUIOBJS) $(FLIBS) $(LIBX_MTX)

usplit : $(SPLOBJ)
	-$(FLINKER) -o $@ $(LIBS) $(FFLAGS) $(SPLOBJ) $(FLIBS) $(LIBX_MTX)

check :
	ftnchek $(URSRCS) | less
	ftnchek $(UR3SRCS) | less
	ftnchek $(GUISRCS) | less
	ftnchek $(SPLSRC) | less

clean : 
	rm -f ./#*# ./*~ *.o core

veryclean : clean
	rm -f ufread ufread3 guiread usplit

distclean: clean
	cd ../lib && make clean

new:
	-mkdir ../toolsnew
	cp Makefile ../toolsnew
	cp *.f ../toolsnew

