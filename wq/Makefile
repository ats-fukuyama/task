### $Id: Makefile,v 1.5 2010/05/18 03:57:54 maruyama Exp $
include ../make.header

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

FFLAGS = $(DFLAGS)
#FFLAGS = $(OFLAGS)

MODINCLUDE= -I./$(MOD) -I../lib/$(MOD) -I../../bpsd/$(MOD)

SRCS = wqtens.f90 wqsolv.f90 wqgout.f90

OBJS = $(SRCS:.f90=.o)

LIBS = libwq.a ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: wq

libwq.a: $(OBJS)
	$(LD) $(LDFLAGS) libwq.a $(OBJS)

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

wq: $(LIBS) wqmain.o
	$(FCFREE) wqmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(GFLIBS)

clean:
	rm -f core *.o ./*~ 

veryclean: clean
	rm -f wq

wqmain.o: wqmain.f90
wqtens.o: wqtens.f90
wqsolv.o: wqsolv.f90
wqgout.o: wqgout.f90

