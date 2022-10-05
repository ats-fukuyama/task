### Makefile
include ../../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = -g $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE= -I../../../bpsd/$(MOD) -I../$(MOD) -I$(MOD)

SRCNC=nclass_mod.f rarray_sum.f \
      rarray_copy.f rarray_zero.f u_erf.f \
      u_lu_backsub.f u_lu_decomp.f write_mod.f
INCNC=pamx_mi.inc pamx_ms.inc pamx_mz.inc
OBJNC=$(SRCNC:.f=.o)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: libtr_nc.a

libtr_nc.a: $(OBJNC)
	$(LD) $(LDFLAGS) libtr_nc.a $(OBJNC)

clean : 
	-rm -f ./#* ./*~ *.o

veryclean : clean

nclass_mod.o  : nclass_mod.f pamx_mi.inc pamx_ms.inc pamx_mz.inc
rarray_sum.o  : rarray_sum.f
rarray_copy.o : rarray_copy.f
rarray_zero.o : rarray_zero.f
u_erf.o       : u_erf.f
u_lu_backsub.o: u_lu_backsub.f
u_lu_decomp.o : u_lu_decomp.f
write_mod.o   : write_mod.f
