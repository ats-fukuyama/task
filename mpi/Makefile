### $Id$ ###
include ../make.header

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

SRCN = nompilib.f
SRCM = mpilib.f

OBJN = $(SRCN:.f=.o)
OBJM = $(SRCM:.f=.o)

MODINCLUDE= -I./$(MOD) -I../mtxp

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : $(OBJN) $(OBJM)

check : 
	ftnchek -array=0 -usage=no-com-block-unused $(SRCN) $(SRCM) | less

new :
	mkdir ../mpinew
	cp -f Makefile ../mpinew
	cp -f *.f ../mpinew
	cp -f *.inc ../mpinew

clean:
	rm -f core *.a *.o ./*~

veryclean: clean

nompilib.o : nompilib.f mpicom.inc
mpilib.o : mpilib.f mpicom.inc mpilib.inc mpif.inc
