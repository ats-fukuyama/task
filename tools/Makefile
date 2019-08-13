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

LIBS=../lib/libtask.a ../lib/libmds.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@

all : ufread ufread3 guiread usplit
#all : ufread usplit

../lib/liobtask.a:
	(cd ../lib; make libtask.a)
../lib/liobmds.a:
	(cd ../lib; make libmds.a)

ufread : $(UROBJS) $(LIBS)
	-$(FCFIXED) -o $@ $(LIBS) $(FFLAGS) $(UROBJS) $(GFLIBS) $(MDSLIB) $(LIBLA)

ufread3 : $(UR3OBJS) $(LIBS)
	-$(FCFIXED) -o $@ $(LIBS) $(FFLAGS) $(UR3OBJS) $(GLLIBS) $(MDSLIB) $(LIBLA)

guiread : $(GUIOBJS) $(LIBS)
	-$(FCFIXED) -o $@ $(LIBS) $(FFLAGS) $(GUIOBJS) $(GLLIBS) $(MDSLIB) $(LIBLA)

usplit : $(SPLOBJ)
	-$(FCFIXED) -o $@ $(FFLAGS) $(SPLOBJ)

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

