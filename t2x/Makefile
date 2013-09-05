### $Id$

include ../make.header

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

MODINCLUDE= -I./$(MOD) -I../../bpsd/$(MOD) -I../lib/$(MOD)

SRCS = libt2.f90

OBJS = $(SRCS:.f90=.o)

LIBS = ../lib/tasklib.a ../../bpsd/bpsdlib.a t2lib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : testt2

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)

t2lib.a: $(OBJS)
	$(LD) $(LDFLAGS) t2lib.a $(OBJS)

testt2 : $(LIBS) testt2.o
	$(FCFIXED) testt2.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MODDIR) $(MODINCLUDE)

clean:
	-rm -f core a.out *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f testt2

testt2.o: testt2.f90
libt2.o: libt2.f90
