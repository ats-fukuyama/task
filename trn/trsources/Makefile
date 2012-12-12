TASKDIR=../../../task
BPSDDIR=../../../bpsd
include $(TASKDIR)/make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

MODINCLUDE=-I $(BPSDDIR)/mod -I $(TASKDIR)/plx/mod -I $(TASKDIR)/lib/mod \
           -I ../trufile/mod -I ../trcore/mod -I $(MOD)

SRC90= trsource.f90

OBJS= $(SRC90:.f90=.o)

.f.o :
	$(FCFIXED) -r8 $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: trsourceslib.a

trsourceslib.a: $(OBJS)
	$(LD) $(LDFLAGS) trsourceslib.a $(OBJS)

clean :
	rm -f ./#* ./*f90~ *.o

veryclean :
	rm -f ./#* ./*~ *.o