TASKDIR=../../../task
BPSDDIR=../../../bpsd
include $(TASKDIR)/make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS) -warn -g

MODINCLUDE=-I $(BPSDDIR)/mod -I $(TASKDIR)/plx/mod -I $(TASKDIR)/lib/mod \
           -I ../trcore/mod -I $(MOD)

SRC90= trufsub.f90 trufcalc.f90 trufile.f90

OBJS= $(SRC90:.f90=.o)

.f.o :
	$(FCFIXED) -r8 $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: trufilelib.a

trufilelib.a: $(OBJS)
	$(LD) $(LDFLAGS) trufilelib.a $(OBJS)

clean :
	rm -f ./#* ./*f90~ *.o

veryclean :
	rm -f ./#* ./*~ *.o