### $Id$ ###
include ../../make.header

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

# for sxf90
#FFLAGS=-Chopt -Wf"-init heap=zero"

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

MODINCLUDE= -I$(MOD) -I../../../bpsd/$(MOD) -I../../lib/$(MOD) \
            -I../../pl/$(MOD)

SRCS = mcnmod.f90 \
       nbidb_040527pb.f90  rkhn2.f90  iodisk41.f90 dinspl.f90 ranu2.f90

SRCX = fort_random.f

OBJS = $(SRCS:.f90=.o) $(SRCX:.f=.o)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : libmcnbi.a

libmcnbi.a: $(OBJS)
	$(LD) $(LDFLAGS) $@ $(OBJS)

clean :
	-rm -f *.o *.a ./*~ .$(MOD)/*.mod

veryclean: clean

