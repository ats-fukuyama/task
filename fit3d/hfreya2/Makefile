### $Id$ ###
include ../../make.header

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

# for sxf90
#FFLAGS=-Chopt -Wf"-A dbl4 -O extendreorder"

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

MODINCLUDE= -I.$(MOD) -I../../../bpsd$(MOD) -I../../task/pl$(MOD)

SRCS = hfrmod.f90 \
       hfreya2z_020502Tnlim_gauss2_rdv01.f90 \
       drive_rdv01.f90 iodisk3_hf.f90 hfrfnc.f90 
SRCT = chs_c05pcf.f

OBJS = $(SRCS:.f90=.o) $(SRCT:.f=.o)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : hfreya2lib.a

hfreya2lib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS)

clean :
	-rm -f *.o *.a ./*~ .$(MOD)/*.mod

veryclean: clean

