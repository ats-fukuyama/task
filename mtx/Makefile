### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)

SRCN = bandcd.f  mtxset.f precndcb.f \
       bcgcdb.f cgscdb.f bstabcdb.f \
       bandcdb.f

SRCM = bandcdm.f bandcdbm.f \
       bcgcdbm.f bcgcdbma.f \
       cgscdbm.f cgscdbma.f \
       bstabcdbm.f bstabcdbma.f \
       precndcbm.f mpisrrc.f

OBJN = $(SRCN:.f=.o)
OBJM = $(SRCM:.f=.o)

.f.o :
	$(FC) $(LFLAGS) -c $< -o $@

all : mtxlib.a mpimtxlib.a
mtxlib.a : $(OBJN)
	$(LD) $(LDFLAGS) -o mtxlib.a $(OBJN)
mpimtxlib.a : $(OBJM)
	$(LD) $(LDFLAGS) -o mpimtxlib.a $(OBJM)

check : 
	ftnchek -array=0 -usage=no-com-block-unused $(SRCS) | less

new:
	mkdir ../mtxnew
	cp -f Makefile ../mtxnew
	cp -f *.f ../mtxnew

clean:
	rm -f core *.o ./*~
veryclean: clean