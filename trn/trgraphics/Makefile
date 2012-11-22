TASKDIR=../../../task
BPSDDIR=../../../bpsd
include $(TASKDIR)/make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE=-I $(BPSDDIR)/mod -I $(TASKDIR)/plx/mod -I $(TASKDIR)/lib/mod \
           -I ../trcore/mod -I ../trufile/mod -I $(MOD)

SRC90= trgsub.f90 trgcom.f90 trgtmp.f90 trgrad.f90 trgexp.f90 trgdgn.f90 \
       trgcmp.f90 trgequ.f90 trgout.f90 

OBJS= $(SRC90:.f90=.o) 

.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: trgraphicslib.a

trgraphicslib.a: $(OBJS)
	$(LD) $(LDFLAGS) trgraphicslib.a $(OBJS)

clean : 
	rm -f ./#* ./*f90~ *.o 

veryclean : 
	rm -f ./#* ./*~ *.o 
