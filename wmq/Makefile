### $Id: Makefile,v 1.5 2010/05/18 03:57:54 maruyama Exp $
include ../make.header

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

FFLAGS = $(DFLAGS)
#FFLAGS = $(OFLAGS)

MODINCLUDE= -I./$(MOD) -I../lib/$(MOD) -I../../bpsd/$(MOD)

SRCS = wmqtens.f90 wmqsolv.f90 wmqgout.f90

OBJS = $(SRCS:.f90=.o)

LIBS = libwmq.a ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: wmq

libwmq.a: $(OBJS)
	$(LD) $(LDFLAGS) libwmq.a $(OBJS)

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../task/lib/libtask.a:
	(cd ../lib; make libtask.a)
../task/lib/libgrf.a:
	(cd ../lib; make libgrf.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)

wmq: $(LIBS) wmqmain.o
	$(FCFREE) wmqmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(GFLIBS)

clean:
	rm -f core *.o ./*~ 

veryclean: clean
	rm -f wmq

wmqmain.o: wmqmain.f90
wmqtens.o: wmqtens.f90
wmqsolv.o: wmqsolv.f90
wmqgout.o: wmqgout.f90

