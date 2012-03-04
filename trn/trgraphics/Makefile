TASKDIR=../../../task
BPSDDIR=../../../bpsd
include $(TASKDIR)/make.header
MOD = ../mod
MODDIR = -module ../mod

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE=-I $(BPSDDIR)/mod -I $(TASKDIR)/plx/mod -I $(TASKDIR)/lib/mod \
           -I $(MOD) 

SRC90= trgcom.f90 trgtmp.f90 trgrad.f90 trgdgn.f90 trgout.f90 

OBJS= $(SRC90:.f90=.o) 

.f.o :
	$(FCFIXED) -r8 $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: trgraphicslib.a

trgraphicslib.a: $(OBJS)
	$(LD) $(LDFLAGS) trgraphicslib.a $(OBJS)

clean : 
	rm -f ./#* ./*f90~ *.o 

veryclean : 
	rm -f ./#* ./*~ *.o 
