###
### This Makefile should be located at task/ti/adpost
###
include ../make.header
include ../mtxp/make.mtxp

.SUFFIXES:
.SUFFIXES: .f .for .f90 .mod .o .a

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE= -I./$(MOD) -I../../bpsd/$(MOD) -I../mtxp/$(MOD) -I../lib/$(MOD)

SRCS = lib-adpost.f90

OBJS = $(SRCS:.f90=.o)

LIBS = lib-adpost.a ../mtxp/libmtxmpi.o ../lib/tasklib.a ../../bpsd/bpsdlib.a

.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : test-adpost

lib-adpost.a: $(OBJS)
	$(LD) $(LDFLAGS) $@ $(OBJS)

libs:
	(cd ../lib; make)

test-adpost : $(LIBS) test-adpost.o
	$(FLINKER) test-adpost.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MODDIR) $(MODINCLUDE)

clean:
	-rm -f core *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f test-adpost

lib-adpost.o: lib-adpost.f90
test-adpost.o: test-adpost.f90 lib-adpost.f90
