###
### This Makefile should be located at task/open-adas/adf11/adf11-lib
### Data file should be located at task/open-adas, e.g. task/open-adas/acd50
###
include ../../../make.header

.SUFFIXES:
.SUFFIXES: .f .for .f90 .mod .o .a

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE= -I./$(MOD) -I../../../lib/$(MOD)

SRCS = \
xxdata_11.for \
xxrptn.for \
i4unit.for \
i4fctn.for \
xxword.for \
xxcase.for \
xfelem.for \
xxslen.for 

SRCX = \
lib-adf11.f90 \

OBJS = $(SRCS:.for=.o) $(SRCX:.f90=.o) 

LIBS = adf11-lib.a ../../../lib/tasklib.a

.for.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : conv-adf11 test-adf11

adf11-lib.a: $(OBJS)
	$(LD) $(LDFLAGS) adf11-lib.a $(OBJS)

libs:
	(cd ../../../lib; make)

conv-adf11 : $(LIBS) conv-adf11.o
	$(FCFIXED) conv-adf11.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

test-adf11 : $(LIBS) test-adf11.o
	$(FCFIXED) test-adf11.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

clean:
	-rm -f core conv-adf11 test-adf11 *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f test-xxdata_11

xxdata_11.o: xxdata_11.for
xxrptn.o: xxrptn.for
i4unit.o: i4unit.for
i4fctn.o: i4fctn.for
xxword.o: xxword.for
xxcase.o: xxcase.for
xfelem.o: xfelem.for
xxslen.o: xxslen.for 
