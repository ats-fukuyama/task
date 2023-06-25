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

TARGETP = ./t2
FFLAGS= $(OFLAGS)
#FFLAGS= $(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a 

#### ####

MODINCLUDE= -I$(MOD) -I../mtxp/$(MOD) -I../lib/$(MOD) -I../../bpsd/$(MOD) \
              -I../t2x/$(MOD) -I../pl/$(MOD) -I../eq/$(MOD)

SRCS = libt2.f90  t2cnst.f90 t2comm.f90 t2intg.f90 t2ngra.f90 t2vgra.f90\
       t2init.f90 t2parm.f90 t2div.f90  t2vout.f90 t2gout.f90 t2cout.f90\
       t2prep.f90 t2calv.f90 t2coef.f90 t2exec.f90 t2conv.f90 t2save.f90 \
       t2prec.f90 t2step.f90 t2loop.f90 t2menu.f90

OBJS = $(SRCS:.f90=.o) 

LIBS= ../lib/libgrf.a ../lib/libtask.a ../../bpsd/libbpsd.a $(LIBMTX) \
      libt2x.a ../pl/libpl.a ../eq/libeq.a
.f90.o:
	$(FLINKER) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: libs t2

../bpsd/libbpsd.a:
	(cd ../bpsd; make libbpsd.a)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../pl/libpl.a:
	(cd ../pl; make libpl.a)
../eq/libeq.a:
	(cd ../eq;  make libeq.a)
libt2x.a: $(OBJS)
	$(LD) $(LDFLAGS) libt2x.a $(OBJS)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../mtxp; make)

t2:$(LIBS) t2main.o
	$(FLINKER)  t2main.o -o $@ $(FFLAGS) $(LIBS)  $(FLIBS) $(LIBX_MTX)

clean:
	-rm -f *.o *.mod *.a ./mod/*.mod t2 *.vtk *.dat *~
cleand:
	-rm -f *.vtk *.dat 
veryclean: clean cleand
	-rm -f t2
t2main.o: t2main.f90 t2comm.f90 t2cnst.f90
t2menu.o: t2menu.f90 t2comm.f90 t2cnst.f90
t2init.o: t2init.f90 t2comm.f90 t2cnst.f90
t2parm.o: t2parm.f90 t2comm.f90 t2cnst.f90
t2prep.o: t2prep.f90 t2comm.f90 t2cnst.f90
t2div.o:  t2div.f90  t2comm.f90 t2cnst.f90
t2intg.o: t2intg.f90 t2comm.f90 t2cnst.f90
t2ngra.o: t2ngra.f90 t2comm.f90 t2cnst.f90
t2vgra.o: t2vgra.f90 t2comm.f90 t2cnst.f90
t2wrap.o: t2wrap.f90 t2comm.f90 t2cnst.f90
t2mfcs.o: t2mfcs.f90 t2comm.f90 t2cnst.f90
t2calv.o: t2calv.f90 t2comm.f90 t2cnst.f90
t2coef.o: t2coef.f90 t2comm.f90 t2cnst.f90
t2exec.o: t2exec.f90 t2comm.f90 t2cnst.f90
t2step.o: t2step.f90 t2comm.f90 t2cnst.f90
t2writ.o: t2writ.f90 t2comm.f90 t2cnst.f90
t2loop.o: t2loop.f90 t2comm.f90 t2cnst.f90
t2term.o: t2term.f90 t2comm.f90 t2cnst.f90
libt2.o: libt2.f90
