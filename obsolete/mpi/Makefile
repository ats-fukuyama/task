### $Id$ ###
include ../make.header

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE=-I../mtxp/mod
SRCS = mpilib.f

OBJS = $(SRCS:.f=.o)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : $(OBJS)

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

mpilib.o : mpilib.f
