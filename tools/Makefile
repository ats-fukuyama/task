include ../make.header

#FFLAGS = $(DFLAGS)
FFLAGS = $(OFLAGS)

MDSLIB=-lMdsLib

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
	$(FC) $(UROBJS) $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MDSLIB)

ufread3 : $(UR3OBJS) $(LIBS)
	-$(FC) $(UR3OBJS) $(LIBS) -o $@ $(FFLAGS) $(GLLIBS) $(MDSLIB)

guiread : $(GUIOBJS) $(LIBS)
	-$(FC) $(GUIOBJS) $(LIBS) -o $@ $(FFLAGS) $(GLLIBS) $(MDSLIB)

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

