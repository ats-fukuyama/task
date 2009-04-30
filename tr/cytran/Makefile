include ../../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE= -I$(MOD) -I../../../bpsd/$(MOD) -I../../eq/$(MOD) -I../../pl/$(MOD)

SRC90=cytran_mod.f90 tr_cytran_mod.f90

OBJS= $(SRC90:.f90=.o)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: $(OBJS)

clean : 
	rm -f ./#* ./*~ *.o *.a *.mod $(MOD)/*.mod

veryclean : clean

tr_cytran.o : tr_cytran.f90 ../../../bpsd/bpsd_kinds.f90
cytran_mod.o : cytran_mod.f90 ../../../bpsd/bpsd_kinds.f90
