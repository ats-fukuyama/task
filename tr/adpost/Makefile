###
### This Makefile should be located at task/open-adas/adf11/adf11-lib
### Data file should be located at task/open-adas, e.g. task/open-adas/acd50
###
include ../../make.header

.SUFFIXES:
.SUFFIXES: .f .for .f90 .mod .o .a

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE= -I./$(MOD) -I../../../bpsd/mod  -I../../lib/mod

SRCS = lib-adpost.f90

OBJS = $(SRCS:.f90=.o)

LIBS = adpost-lib.a ../../lib/tasklib.a ../../../bpsd/bpsdlib.a

.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : test-adpost

adpost-lib.a: $(OBJS)
	$(LD) $(LDFLAGS) adpost-lib.a $(OBJS)

libs:
	(cd ../../lib; make)

test-adpost : $(LIBS) test-adpost.o
	$(FCFIXED) test-adpost.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MODDIR) $(MODINCLUDE)

clean:
	-rm -f core *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f test-adpost

lib-adpost.o: lib-adpost.f90
test-adpost.o: test-adpost.f90 lib-adpost.f90
