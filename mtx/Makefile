### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)

SRCS = bandcd.f bandcdm.f bandcdbm.f \
       bcgcdb.f bcgcdbm.f bcgcdbma.f \
       cgscdb.f cgscdbm.f cgscdbma.f \
       bstabcdb.f bstabcdbm.f bstabcdbma.f \
       bandcdb.f\
       precndcb.f precndcbm.f mpisrrc.f mtxset.f

OBJS = $(SRCS:.f=.o)

.f.o :
	$(MFC) $(LFLAGS) -c $< -o $@

all : $(OBJS)

check : 
	ftnchek -array=0 -usage=no-com-block-unused $(SRCS) | less

new:
	mkdir ../mtxnew
	cp -f Makefile ../mtxnew
	cp -f *.f ../mtxnew

clean:
	rm -f core *.o ./*~
veryclean: clean