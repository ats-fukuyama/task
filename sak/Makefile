### Madefile for cv
include ../make.header
include ../mtxp/make.mtxp

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../lib/$(MOD) -I../../bpsd/$(MOD)

SRCS = sakcomm.f90 saksub.f90 sak1.f90 sak2.f90 sak3.f90 sak4.f90

OBJS = $(SRCS:.f90=.o)
FLINKER=$(FCFREE)

LIBS = libsak.a  ../lib/libgrf.a ../lib/libtask.a

all : libs sak

libsak.a: $(OBJS)
	$(LD) $(LDFLAGS) libsak.a $(OBJS)

libs:
	(cd ../lib; make libgrf.a)
	(cd ../lib; make libtask.a)

sak : libs libsak.a sakmain.o
	$(FLINKER) sakmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

clean:
	-rm -f core a.out *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f sak

SAKCOMM=sakcomm.f90
saksub.o: $(SAKCOMM) saksub.f90
sak1.o: $(SAKCOMM) sak1.f90
sak2.o: $(SAKCOMM) saksub.f90 sak2.f90
sak3.o: $(SAKCOMM) saksub.f90 sak3.f90
sak4.o: $(SAKCOMM) saksub.f90 sak3.f90
sakmain.o: $(SAKCOMM) sak1.f90 sak2.f90 sak3.f90 sak4.f90 sakmain.f90
