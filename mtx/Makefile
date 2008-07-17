### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)

SRCN = bandcd.f  mtxset.f precndcb.f \
       bstabcdb.f \
       bandcdb.f

SRCM = bstabcdbm.f \
       precndcbm.f mpisrrc.f
#      bandcdm.f bandcdbm.f \

OBJN = $(SRCN:.f=.o)
OBJM = $(SRCM:.f=.o)

.f.o :
	$(FCFIXED) $(LFLAGS) -c $< -o $@
.f90.o :
	$(FCFREE) $(LFLAGS) -c $< -o $@

all : mtxlib.a mpimtxlib.a
mtxlib.a : $(OBJN)
	$(LD) $(LDFLAGS) mtxlib.a $(OBJN)
mpimtxlib.a : $(OBJN) $(OBJM)
	$(LD) $(LDFLAGS) mpimtxlib.a $(OBJN) $(OBJM)

check : 
	ftnchek -array=0 -usage=no-com-block-unused $(SRCS) | less

new:
	mkdir ../mtxnew
	cp -f Makefile ../mtxnew
	cp -f *.f ../mtxnew

clean:
	rm -f core *.o *.a ./*~

veryclean: clean
