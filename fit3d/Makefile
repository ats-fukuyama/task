### $Id$ ###
include ../make.header

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

MODINCLUDE= -I./$(MOD) -Ihfreya2/$(MOD) -Imcnbi/$(MOD) -I../../bpsd/$(MOD) -I../pl/$(MOD)

SRCSFREE = mcnbi_all.f90
SRCSFIXED = 

OBJS = $(SRCSFREE:.f90=.o) $(SRCSFIXED:.f=.o)
LIBS = hfreya2/libhfreya2.a mcnbi/libmcnbi.a fit/libfit.a ../lib/libtask.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : fit3d

hfreya2/libhfreya2.a:
	(cd hfreya2; make)
mcnbi/libmcnbi.a:
	(cd mcnbi; make)
fit/libfit.a:
	(cd fit; make)

libs:
	(cd hfreya2; make)
	(cd mcnbi; make)
	(cd fit; make)

fit3d:  $(LIBS) $(OBJS)
	$(FCFIXED) $(OBJS) $(LIBS) -o fit3d $(FFLAGS) $(FLIBS)

clean :
	(cd hfreya2; make clean)
	(cd mcnbi; make clean)
	(cd fit; make clean)
	-rm -f *.o *.a ./*~ .$(MOD)/*.mod

veryclean: clean
	-rm -f fit3d

