include ../make.header
include ../mtxp/make.mtxp
#LIBX_MTX=$(LIBX_MTX_BND)
LIBX_MTX=$(LIBX_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_KSP)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS=$(OFLAGS) -openmp
FFLAGS=$(DFLAGS) -fopenmp
#FFLAGS=$(OFLAGS) -fopenmp #-Wl, -stack_size,0x80000000,-stack_addr,0xf0000000

MODINCLUDE=-I $(MOD) -I../../bpsd/$(MOD) -I../mtxp/$(MOD) -I../lib/$(MOD)
FINCLUDE=-I /usr/local/mpich332-gfortran-clang/include
#FINCLUDE=-I /usr/local/mpich32-ifort-gcc5/include
#FINCLUDE=-I /usr/local/mpich3-intel140/include
#FINCLUDE=-I /home/hoge/lib/openmpi/include

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE) $(FINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE) $(FINCLUDE)

SRCSFREE = libgrf/plot2.f90 libgrf/grf2dpexec.f90 libgrf/grfpconvert.f90 \
           libgrf/grdpconvert.f90 libgrf/grf2dp.f90 libgrf/grd2dp.f90 \
           piclib.f90 piccomm.f90 picinit.f90 picparm.f90 \
           picsub.f90 picprep.f90 picexec.f90 \
           picgout.f90 picmenu.f90

OBJS = $(SRCSFREE:.f90=.o)

LIBS = ../../bpsd/libbpsd.a ../lib/libgrf.a ../lib/libtask.a libpic.a

all : pic

#disp: disp

../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
libpic.a: $(OBJS)
	$(LD) $(LDFLAGS) libpic.a $(OBJS)

libs:
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)
	(cd ../../bpsd; make libbpsd.a)
	(cd ../mtxp; make)

pic : $(LIBS) $(LIBX_MTX) picmain.o
	$(FLINKER) picmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX) -lfftw3

disp: disp.o
	$(FCFREE) disp.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

testplot2 : $(LIBS) libgrf/testplot2.o
	$(MF90) libgrf/testplot2.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) \
           $(MODDIR) $(MODINCLUDE)

clean:
	-rm -f core a.out *.o *.mod ./*~ ./#* *.a $(MOD)/*.mod libgrf/*.o

veryclean: clean
	-rm -f pic

picmain.o: picmain.f90 piccomm.f90 picinit.f90 picparm.f90 picmenu.f90
picmenu.o: picmenu.f90 piccomm.f90 picparm.f90 picprep.f90 picexec.f90 \
           picgout.f90
picinit.o: picinit.f90 piccomm.f90
picparm.o: picparm.f90 piccomm.f90
picprep.o: picprep.f90 piccomm.f90 picsub.f90 piclib.f90
picexec.o: picexec.f90 piccomm.f90 picsub.f90 piclib.f90
picgout.o: picgout.f90 piccomm.f90 picparm.f90 libgrf/grd2dp.f90
picsub.o: picsub.f90
piclib.o: piclib.f90
#disp.o: disp.f90
libgrf/grd2dp.o: libgrf/grd2dp.f90 libgrf/grdpconvert.f90 \
                 libgrf/grf2dpexec.f90
libgrf/grf2dp.o: libgrf/grf2dp.f90 libgrf/grfpconvert.f90 \
                 libgrf/grf2dpexec.f90
libgrf/grdpconvert.o: libgrf/grdpconvert.f90
libgrf/grfpconvert.o: libgrf/grfpconvert.f90
libgrf/grf2dpexec.o: libgrf/grf2dpexec.f90 libgrf/plot2.f
