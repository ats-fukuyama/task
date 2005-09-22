include ../make.header

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

LIBS=../lib/tasklib.a ../lib/mdslib.a $(OBJLA)

all : ufread ufread3 guiread usplit

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../lib/mdslib.a:
	(cd ../lib; make mdslib.a)

ufread : $(UROBJS) $(LIBS)
	$(FC) -o $@ $(LIBS) $(FFLAGS) $(UROBJS) $(GFLIBS) $(MDSLIB) $(OBJLA)

ufread3 : $(UR3OBJS) $(LIBS)
	-$(FC) -o $@ $(LIBS) $(FFLAGS) $(UR3OBJS) $(GLLIBS) $(MDSLIB) $(OBJLA)

guiread : $(GUIOBJS) $(LIBS)
	-$(FC) -o $@ $(LIBS) $(FFLAGS) $(GUIOBJS) $(GLLIBS) $(MDSLIB) $(OBJLA)

usplit : $(SPLOBJ)
	-$(FC) -o $@ $(FFLAGS) $(SPLOBJ)

check :
	ftnchek $(URSRCS) | less
	ftnchek $(UR3SRCS) | less
	ftnchek $(GUISRCS) | less
	ftnchek $(SPLSRC) | less

clean : 
	rm -f ./#*# ./*~ *.o core

veryclean : clean
	rm -f ufread ufread3 guiread usplit

new:
	-mkdir ../toolsnew
	cp Makefile ../toolsnew
	cp *.f ../toolsnew

