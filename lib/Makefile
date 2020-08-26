### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

LFLAGS = $(OFLAGS)
#LFLAGS = $(DFLAGS)
#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE= -I $(MOD) -I ../../bpsd/$(MOD)

SRCSFREE  =  libbrn.f90 libinv.f90 libbnd.f90 \
             libpcgpme.f90 \
             librk.f90 libsympl.f90 \
             libspl.f90 libspl2d.f90 libspl3d.f90 \
             libeq.f90 libgut.f90 \
             libchar.f90 libkio.f90 \
             libfio.f90 ufile.f90 libitp.f90 libfem.f90 \
             libde.f90 libdsp.f90 libell.f90 libspf.f90 \
	     libfft.f90 \
             libqsort.f90 \
             libufile/ufinit.f90 libufile/uflist.f90 libufile/ufread.f90 \
             libbes.f90  netlib-bes.f90

SRCSGRF   =  libgrf/grftype.f90 libgrf/grfutils.f90 \
	     libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 \
	     libgrf/grfxyexec.f90 \
	     libgrf/grfconvert.f90 libgrf/grdconvert.f90 \
             libgrf/grf1d.f90 libgrf/grd1d.f90 \
	     libgrf/grf2d.f90 libgrf/grd2d.f90 \
             libgrf/grd1dframe.f90 libgrf/grd2dframe.f90 \
	     libgrf/grdxy.f90 \
	     libgrf.f90

#SRCSFIXED =  netlib-bes.f $(LAPACK)
SRCSFIXED =  $(LAPACK)

OBJSFREE  = $(SRCSFREE:.f90=.o)
OBJSFIXED = $(SRCSFIXED:.f=.o)
OBJSGRF = $(SRCSGRF:.f90=.o)
OBJSMDS = $(MDSPLUS:.f=.o)

.f.o :
	$(FCFIXED) $(LFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(LFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : libtask.a libgrf.a libmds.a

libtask.a: $(OBJSFREE) $(OBJSFIXED)
	$(LD) $(LDFLAGS) libtask.a $(OBJSFREE) $(OBJSFIXED)

libgrf.a: $(OBJSGRF)
	$(LD) $(LDFLAGS) libgrf.a $(OBJSGRF)

libmds.a: $(OBJSMDS)
	$(LD) $(LDFLAGS) libmds.a $(OBJSMDS)

check :
	ftnchek $(SRCSFREE) $(SRCSFIXED) $(MDSPLUS) | less

clean:
	-rm -f core* a.out *.o libgrf/*.o libufile/*.o *.mod *.a ./*~ ./#* $(MOD)/*.mod

veryclean: clean

testfem: libgrf.a libtask.a libtestfem.o fem_calc.o fem_calc_15.o fem_calc_16.o testfem.o 
	$(FCFREE) libtestfem.o fem_calc.o fem_calc_15.o fem_calc_16.o testfem.o libtask.a $(FFLAGS) -o $@ $(FLIBS) $(MODINCLUDE)

fem_calc.o:	libtestfem.f90 fem_calc.f90
fem_calc_15.o:	libtestfem.f90 fem_calc_15.f90
fem_calc_16.o:	libtestfem.f90 fem_calc_16.f90

testfem_func: libgrf.a libtask.a testfem_func.o 
	$(FCFREE) testfem_func.o libgrf.a libtask.a $(FFLAGS) -o $@ $(FLIBS) $(MODINCLUDE)

testbes: libgrf.a libtask.a testbes.o 
	$(FCFIXED) testbes.o libgrf.a libtask.a $(FFLAGS) -o $@ $(FLIBS)
testlamb: libgrf.a libtask.a testlamb.o 
	$(FCFREE) testlamb.o libgrf.a libtask.a $(FFLAGS) -o $@ $(FLIBS)

testell: libgrf.a libtask.a testell.o 
	$(FCFIXED) testell.o libgrf.a libtask.a $(FFLAGS) -o $@ $(FLIBS)

testspl: libgrf.a libtask.a testspl.o 
	$(FCFIXED) testspl.o libgrf.a libtask.a $(FFLAGS) -o $@ $(FLIBS)

testspl3d: libgrf.a libtask.a testspl3d.o 
	$(FCFIXED) testspl3d.o libgrf.a libtask.a $(FFLAGS) -o $@ $(FLIBS)

mmconvert: mmconvert.o
	$(FCFIXED) mmconvert.o $(FFLAGS) -o $@

testgrf: libgrf.a libtask.a testgrf.o 
	$(FCFIXED) testgrf.o libgrf.a libtask.a $(FFLAGS) -o $@ $(FLIBS)
testqsort: libtask.a testqsort.o 
	$(FCFREE) testqsort.o libgrf.a libtask.a $(FFLAGS) -o $@ $(FLIBS)
testgrd: libtask.a testgrd.o 
	$(FCFIXED) testgrd.o libgrf.a libtask.a $(FFLAGS) -o $@ $(FLIBS)

libgrf/grftype.o: libgrf/grftype.f90
libgrf/grfutils.o: libgrf/grfutils.f90 libgrf/grftype.f90
libgrf/grfconvert.o: libgrf/grfconvert.f90 libgrf/grfutils.f90 \
                     libgrf/grftype.f90
libgrf/grdconvert.o: libgrf/grdconvert.f90 libgrf/grfutils.f90 \
                     libgrf/grftype.f90
libgrf/grf1dexec.o: libgrf/grf1dexec.f90 libgrf/grfutils.f90 \
                    libgrf/grftype.f90
libgrf/grf2dexec.o: libgrf/grf2dexec.f90 libgrf/grfutils.f90 \
                    libgrf/grftype.f90
libgrf/grf1d.o: libgrf/grf1d.f90 libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 \
                libgrf/grfconvert.f90 libgrf/grftype.f90
libgrf/grd1d.o: libgrf/grd1d.f90 libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 \
                libgrf/grdconvert.f90 libgrf/grftype.f90
libgrf/grf2d.o: libgrf/grf2d.f90 libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 \
                libgrf/grfconvert.f90 libgrf/grftype.f90
libgrf/grd2d.o: libgrf/grd2d.f90 libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 \
                libgrf/grdconvert.f90 libgrf/grftype.f90
libgrf/grd1dframe.o: libgrf/grd1dframe.f90 \
                libgrf/grdconvert.f90 libgrf/grftype.f90
libgrf/grd2dframe.o: libgrf/grd2dframe.f90 \
                libgrf/grdconvert.f90 libgrf/grftype.f90

libufile/ufinit.o: libufile/ufinit.f90
libufile/uflist.o: libufile/uflist.f90 libufile/ufinit.f90
libufile/ufread.o: libufile/ufread.f90 libufile/ufinit.f90 libufile/uflist.f90
