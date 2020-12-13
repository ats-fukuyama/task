### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

LFLAGS = $(OFLAGS)
#LFLAGS = $(DFLAGS)
#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE= -I $(MOD)

SRCSFREE  =  task_kinds.f90 \
             libbnd.f90 libbrn.f90 libde.f90 libdsp.f90 \
             libell.f90 libeq.f90  libfem.f90 libfft.f90 \
             libgut.f90 libinv.f90 libitp.f90 libpcgpme.f90 \
             libqsort.f90 librk.f90 libspf.f90 libsympl.f90 \
             libspl.f90 libspl2d.f90 libspl3d.f90 \
             netlib-bes.f90 libbes.f90 \
             libchar.f90 libkio.f90 libfio.f90 \
             libufile/ufinit.f90 libufile/uflist.f90 libufile/ufread.f90 \
             ufile.f90

SRCSGRF   =  libgrf/grftype.f90 libgrf/grfutils.f90 \
	     libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 \
	     libgrf/grfxyexec.f90 \
	     libgrf/grfconvert.f90 libgrf/grdconvert.f90 \
             libgrf/grf1d.f90 libgrf/grd1d.f90 \
	     libgrf/grf2d.f90 libgrf/grd2d.f90 \
	     libgrf/grdxy.f90 \
             libgrf/grd1dframe.f90 libgrf/grd2dframe.f90 \
	     libgrf.f90

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
	-rm testbes testell testfem testfem_func testgrd testgrf testlamb testqsort testspl testspl3d

veryclean: clean

testfem: libgrf.a libtask.a libtestfem.o fem_calc.o fem_calc_15.o fem_calc_16.o testfem-sub.o testfem.o 
	$(FCFREE) libtestfem.o fem_calc.o fem_calc_15.o fem_calc_16.o testfem-sub.o testfem.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS) $(MODINCLUDE)

testfem-sub.o: $(TASKCOMM) testfem-sub.f90
libtestfem.o:  $(TASKCOMM) libtestfem.f90
fem_calc.o:    $(TASKCOMM) libtestfem.f90 fem_calc.f90
fem_calc_15.o: $(TASKCOMM) libtestfem.f90 fem_calc_15.f90
fem_calc_16.o: $(TASKCOMM) libtestfem.f90 fem_calc_16.f90


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

testgrf: libgrf.a libtask.a testgrf.o 
	$(FCFIXED) testgrf.o libgrf.a libtask.a $(FFLAGS) -o $@ $(FLIBS)
testgrd: libtask.a testgrd.o 
	$(FCFIXED) testgrd.o libgrf.a libtask.a $(FFLAGS) -o $@ $(FLIBS)

testqsort: libtask.a testqsort.o 
	$(FCFREE) testqsort.o libgrf.a libtask.a $(FFLAGS) -o $@ $(FLIBS)

TASKCOMM=task_kinds.f90

libbnd.o: $(TASKCOMM) libbnd.f90
libbrn.o: $(TASKCOMM) libbrn.f90
libde.o: $(TASKCOMM) libde.f90
libdsp.o: $(TASKCOMM) libdsp.f90
libell.o: $(TASKCOMM) libell.f90
libeq.o: $(TASKCOMM) libeq.f90
libfem.o: $(TASKCOMM) libfem.f90
libfft.o: $(TASKCOMM) libfft.f90
libgut.o: $(TASKCOMM) libgut.f90
libinv.o: $(TASKCOMM) libinv.f90
libitp.o: $(TASKCOMM) libitp.f90
libpcgpme.o: $(TASKCOMM) libpcgpme.f90
libqsort.o: $(TASKCOMM) libqsort.f90
librk.o: $(TASKCOMM) librk.f90
libspf.o: $(TASKCOMM) libspf.f90
libspl.o: $(TASKCOMM) libspl.f90
libspl2d.o: $(TASKCOMM) libspl2d.f90
libspl3d.o: $(TASKCOMM) libspl3d.f90
libsympl.o: $(TASKCOMM) libsympl.f90

netlib-bes.o: libspf.f90 netlib-bes.f90
libbes.o: $(TASKCOMM) netlib-bes.f90 libbes.f90

libchar.o: libchar.f90
libkio.o: libkio.f90
libfio.o: libfio.f90

libufile/ufinit.o: libufile/ufinit.f90
libufile/uflist.o: libufile/uflist.f90 libufile/ufinit.f90
libufile/ufread.o: libufile/ufread.f90 libufile/ufinit.f90 libufile/uflist.f90
ufile.o: ufile.f90

libgrf/grftype.o: libgrf/grftype.f90
libgrf/grfutils.o: libgrf/grftype.f90 libgrf/grfutils.f90
libgrf/grfconvert.o: libgrf/grftype.f90  libgrf/grfutils.f90 \
                     libgrf/grfconvert.f90
libgrf/grdconvert.o: $(TASKCOMM) libgrf/grftype.f90  libgrf/grfutils.f90 \
                     libgrf/grdconvert.f90
libgrf/grf1dexec.o: libgrf/grftype.f90 libgrf/grfutils.f90 \
                    libgrf/grf1dexec.f90
libgrf/grf2dexec.o: libgrf/grftype.f90 libgrf/grfutils.f90 \
                    libgrf/grf2dexec.f90
libgrf/grfxyexec.o: libgrf/grftype.f90 libgrf/grfutils.f90 \
                    libgrf/grfxyexec.f90
libgrf/grf1d.o: libgrf/grftype.f90 libgrf/grfconvert.f90 \
                libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 libgrf/grf1d.f90
libgrf/grd1d.o: $(TASKCOMM) libgrf/grftype.f90 libgrf/grdconvert.f90 \
                libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 libgrf/grd1d.f90
libgrf/grf2d.o: libgrf/grftype.f90 libgrf/grfconvert.f90 \
                libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 libgrf/grf2d.f90
libgrf/grd2d.o: $(TASKCOMM) libgrf/grftype.f90 libgrf/grdconvert.f90 \
                libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 libgrf/grd2d.f90
libgrf/grdxy.o: $(TASKCOMM) libgrf/grftype.f90 libgrf/grdconvert.f90 \
                libgrf/grfxyexec.f90 libgrf/grdxy.f90
libgrf/grd1dframe.o: $(TASKCOMM) libgrf/grftype.f90 libgrf/grdconvert.f90 \
                libgrf/grd1dframe.f90
libgrf/grd2dframe.o: $(TASKCOMM) libgrf/grftype.f90 libgrf/grdconvert.f90 \
                libgrf/grd2dframe.f90
libgrf.o: libgrf/grftype.f90 libgrf/grfutils.f90 libgrf/grf1d.f90 \
          libgrf/grd1d.f90 libgrf/grf2d.f90 libgrf/grd2d.f90 libgrf/grdxy.f90 \
          libgrf/grd2dframe.f90 libgrf/grd1dframe.f90 libgrf.f90 


