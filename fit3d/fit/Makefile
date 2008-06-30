### $Id$ ###
include ../../make.header

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

MODINCLUDE= -I.$(MOD) -I../../../bpsd$(MOD) -I../../task/pl$(MOD)

SRCSFREE = fit_read20_050323.f90
SRCSFIXED = depsum_070416.f

OBJS = $(SRCSFREE:.f90=.o) $(SRCSFIXED:.f=.o)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : fitlib.a

fitlib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS)

clean :
	-rm -f *.o *.a ./*~ .$(MOD)/*.mod

veryclean: clean

