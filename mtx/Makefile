### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)

SRCM = bandcdm.f bandcdbm.f \
       bcgcdbm.f bcgcdbma.f \
       cgscdbm.f cgscdbma.f \
       bstabcdbm.f bstabcdbma.f \
       precndcbm.f mpisrrc.f

SRCN = bandcd.f  mtxset.f precndcb.f \
       bcgcdb.f cgscdb.f bstabcdb.f \
       bandcdb.f

OBJM = $(SRCM:.f=.o)
OBJN = $(SRCN:.f=.o)

.f.o :
	$(FC) $(LFLAGS) -c $< -o $@

all : $(OBJM) $(ONJN)
mpi : $(OBJM)
nompi : $(OBJN)

check : 
	ftnchek -array=0 -usage=no-com-block-unused $(SRCS) | less

new:
	mkdir ../mtxnew
	cp -f Makefile ../mtxnew
	cp -f *.f ../mtxnew

clean:
	rm -f core *.o ./*~
veryclean: clean