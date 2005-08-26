include ../make.header

#FFLAGS = $(DFLAGS)
FFLAGS = $(OFLAGS)

URSRCS=ufread.f ufsub.f
UR3SRCS=ufread3.f ufsub.f
GUISRCS=guiread.f ufsub.f
UROBJS=$(URSRCS:.f=.o)
UR3OBJS=$(UR3SRCS:.f=.o)
GUIOBJS=$(GUISRCS:.f=.o)

LIBS=../lib/tasklib.a ../lib/mdslib.a

all : ufread ufread3 guiread

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

check :
	ftnchek $(GUISRCS) | less

clean : 
	rm -f ./#*# ./*~ *.o core

veryclean : clean
	rm -f ufread ufread3 guiread

new:
	-mkdir ../toolsnew
	cp Makefile ../toolsnew
	cp *.f ../toolsnew

