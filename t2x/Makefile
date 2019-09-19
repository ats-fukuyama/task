### $Id$

include ../make.header

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

MODINCLUDE= -I./$(MOD) -I../../bpsd/$(MOD) -I../lib/$(MOD)

SRCS = libt2.f90

OBJS = $(SRCS:.f90=.o)

LIBS = ../lib/libtask.a ../../bpsd/libbpsd.a libt2.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : testt2

../lib/libtask.a:
	(cd ../lib; make libtask.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)

libt2.a: $(OBJS)
	$(LD) $(LDFLAGS) libt2.a $(OBJS)

testt2 : $(LIBS) testt2.o
	$(FCFIXED) testt2.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MODDIR) $(MODINCLUDE)

clean:
	-rm -f core a.out *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f testt2

testt2.o: testt2.f90
libt2.o: libt2.f90
