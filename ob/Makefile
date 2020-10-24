### Madefile for ob
include ../make.header
include ../mtxp/make.mtxp

LIB_MTX=$(LIB_MTX_BND)
LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_PCG)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIB_MTX=$(LIBXMTX_KSP)
#LIBX_MTX=$(LIBX_MTX_KSP)

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../pl/$(MOD) -I../eq/$(MOD) -I../lib/$(MOD) \
            -I../mtxp/$(MOD) -I../../bpsd/$(MOD)

SRCS = obcomm.f90 obinit.f90 obparm.f90 obview.f90 \
       obprep.f90 obexec.f90 obcalc.f90 obgout.f90 \
       obmenu.f90

SRCX = ../eq/wmeqin.f

OBJS = $(SRCS:.f90=.o) $(SRCX:.f=.o)

LIBS = libob.a  ../eq/libeq.a ../pl/libpl.a \
       ../lib/libgrf.a ../lib/libtask.a ../../bpsd/libbpsd.a

all : libs ob

../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq; make libeq.a)
../eq/wmeqin.o:
	(cd ../eq; make wmeqin.o)
libob.a: $(OBJS)
	$(LD) $(LDFLAGS) libob.a $(OBJS)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libgrf.a)
	(cd ../lib; make libtask.a)
	(cd ../mtxp; make)
	(cd ../pl; make libpl.a)
	(cd ../eq; make libeq.a)
	(cd ../eq; make wmeqin.o)

testf: libs libob.a test.o
	$(FLINKER) test.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

ob : libs libob.a obmain.o
	$(FLINKER) obmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f core a.out *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f ob

OBCOMM=obcomm.f90 ../pl/plcomm.f90
obcomm.o: obcomm.f90
obinit.o: obinit.f90 $(OBCOMM)
obparm.o: obparm.f90 $(OBCOMM)
obview.o: obview.f90 $(OBCOMM)
obprep.o: obprep.f90 $(OBCOMM)
obcalc.o: obcalc.f90 obprep.f90 obexec.f90 $(OBCOMM)
obexec.o: obexec.f90 obprep.f90 $(OBCOMM)
obgout.o: obgout.f90 obprep.f90 ../pl/plgout.f90 ../pl/plprof.f90 $(OBCOMM)
../eq/wmeqin.o: ../eq/wmeqin.f
obmenu.o: obmenu.f90 obparm.f90 obview.f90 obprep.f90 obcalc.f90 obgout.f90  $(OBCOMM)
obmain.o: obmain.f90 ../pl/plinit.f90 ../pl/plparm.f90 obinit.f90 obparm.f90 obmenu.f90 $(OBCOMM)
