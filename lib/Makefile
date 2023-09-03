### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

LFLAGS = $(OFLAGS)
#LFLAGS = $(DFLAGS)
FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)

MODINCLUDE= -I $(MOD)
OBJDIR=./obj

${OBJDIR}/%.o:%.f
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

SRCSFREE  =  task_kinds.f90 task_constants.f90 \
             libspl1d.f90 libspl2d.f90 libspl3d.f90 \
             libderiv.f90 libpol.f90 \
             libbnd.f90 libbrent.f90 libde.f90 libdsp.f90 \
             libell.f90 libfem.f90 libfft.f90 \
             libinv.f90 libpcgpme.f90 \
             libqsort.f90 librk.f90 librkf.f90 \
             libspf.f90 libspf2.f90 libsympl.f90 \
             netlib-bes.f90 libbes.f90 \
             libchar.f90 libkio.f90 libfio.f90 \
             libitp.f90 libsmooth.f90 libplog.f90 \
             libufile/ufinit.f90 libufile/uflist.f90 libufile/ufread.f90 \
             ufile.f90

SRCSGRF   =  libgrf/grftype.f90 libgrf/grfutils.f90 libgrf/grdutils.f90 \
	     libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 \
	     libgrf/grfxyexec.f90 \
	     libgrf/grfconvert.f90 libgrf/grdconvert.f90 \
             libgrf/grf1d.f90 libgrf/grd1d.f90 \
	     libgrf/grf2d.f90 libgrf/grd2d.f90 \
	     libgrf/grdxy.f90 \
             libgrf/grd1dframe.f90 libgrf/grd2dframe.f90 \
	     libgrf.f90

SRCSFIXED =  $(LAPACK)

OBJSFREE  = $(addprefix $(OBJDIR)/,$(SRCSFREE:.f90=.o))
OBJSFIXED = $(addprefix $(OBJDIR)/,$(SRCSFIXED:.f=.o))
OBJSGRF = $(addprefix $(OBJDIR)/,$(SRCSGRF:.f90=.o))
OBJSMDS = $(addprefix $(OBJDIR)/,$(MDSPLUS:.f=.o))


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
	-rm -f core* a.out $(OBJDIR)/*.o $(OBJDIR)/*/*.o $(MOD)/*.mod $(MOD)/*/*.mod *.a ./*~ ./#* $(MOD)/*.mod
	-rm -f testairy testbes testell testfem testfem_func testgrd testgrf testqsort testspl testspl3d testpol testduo testrk testtri

veryclean: clean

testfem: libtask.a libgrf.a $(OBJDIR)/libtestfem.o $(OBJDIR)/fem_calc.o $(OBJDIR)/fem_calc_15.o $(OBJDIR)/fem_calc_16.o $(OBJDIR)/testfem-sub.o $(OBJDIR)/testfem.o 
	$(FCFREE) $(OBJDIR)/libtestfem.o $(OBJDIR)/fem_calc.o $(OBJDIR)/fem_calc_15.o $(OBJDIR)/fem_calc_16.o $(OBJDIR)/testfem-sub.o $(OBJDIR)/testfem.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS) $(MODINCLUDE)

$(OBJDIR)/testfem-sub.o: $(TASKCOMM) testfem-sub.f90
$(OBJDIR)/libtestfem.o:  $(TASKCOMM) libtestfem.f90
$(OBJDIR)/fem_calc.o:    $(TASKCOMM) libtestfem.f90 fem_calc.f90
$(OBJDIR)/fem_calc_15.o: $(TASKCOMM) libtestfem.f90 fem_calc_15.f90
$(OBJDIR)/fem_calc_16.o: $(TASKCOMM) libtestfem.f90 fem_calc_16.f90

testtri: libtask.a $(OBJDIR)/testtri.o 
	$(FCFREE) $(OBJDIR)/testtri.o libtask.a $(FFLAGS) -o $@ $(FLIBS) $(MODINCLUDE)
testduo: libtask.a $(OBJDIR)/testduo.o 
	$(FCFREE) $(OBJDIR)/testduo.o libtask.a $(FFLAGS) -o $@ $(FLIBS) $(MODINCLUDE)

testfem_func: libtask.a libgrf.a $(OBJDIR)/testfem_func.o 
	$(FCFREE) $(OBJDIR)/testfem_func.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS) $(MODINCLUDE)

testairy: libtask.a libgrf.a $(OBJDIR)/testairy.o
	$(FCFREE) $(OBJDIR)/testairy.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS)

testbes: libtask.a libgrf.a $(OBJDIR)/testbes.o
	$(FCFREE) $(OBJDIR)/testbes.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS)

testell: libtask.a libgrf.a $(OBJDIR)/testell.o 
	$(FCFIXED) $(OBJDIR)/testell.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS)

testspl: libtask.a libgrf.a $(OBJDIR)/testspl.o 
	$(FCFIXED) $(OBJDIR)/testspl.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS)

testspl3d: libtask.a libgrf.a $(OBJDIR)/testspl3d.o 
	$(FCFIXED) $(OBJDIR)/testspl3d.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS)

testgrf: libtask.a libgrf.a $(OBJDIR)/testgrf.o 
	$(FCFIXED) $(OBJDIR)/testgrf.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS)
testgrd: libtask.a libgrf.a $(OBJDIR)/testgrd.o 
	$(FCFIXED) $(OBJDIR)/testgrd.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS)

testqsort: libtask.a $(OBJDIR)/testqsort.o 
	$(FCFREE) $(OBJDIR)/testqsort.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS)

testrk: libtask.a $(OBJDIR)/testrk.o
	$(FCFREE) $(OBJDIR)/testrk.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS)
testpol: libtask.a $(OBJDIR)/testpol.o 
	$(FCFREE) $(OBJDIR)/testpol.o libtask.a libgrf.a $(FFLAGS) -o $@ $(FLIBS)

TASKCOMM=task_kinds.f90

$(OBJDIR)/task_kinds.o: task_kinds.f90
$(OBJDIR)/task_constants.o: $(TASKCOMM) task_constants.f90
$(OBJDIR)/libderiv.o: $(TASKCOMM) libderiv.f90
$(OBJDIR)/libpol.o: $(TASKCOMM) libpol.f90
$(OBJDIR)/libbnd.o: $(TASKCOMM) libbnd.f90
$(OBJDIR)/libbrent.o: $(TASKCOMM) libbrent.f90
$(OBJDIR)/libde.o: $(TASKCOMM) libde.f90
$(OBJDIR)/libdsp.o: $(TASKCOMM) libdsp.f90
$(OBJDIR)/libell.o: $(TASKCOMM) libell.f90
$(OBJDIR)/libfem.o: $(TASKCOMM) libfem.f90
$(OBJDIR)/libfft.o: $(TASKCOMM) libfft.f90
$(OBJDIR)/libinv.o: $(TASKCOMM) libinv.f90
$(OBJDIR)/libpcgpme.o: $(TASKCOMM) libpcgpme.f90
$(OBJDIR)/libqsort.o: $(TASKCOMM) libqsort.f90
$(OBJDIR)/librk.o: $(TASKCOMM) librk.f90
$(OBJDIR)/librkf.o: $(TASKCOMM) librkf.f90
$(OBJDIR)/libspf.o: $(TASKCOMM) libspf.f90
$(OBJDIR)/libspf2.o: $(TASKCOMM) libspf2.f90 special_functions.f90
$(OBJDIR)/libspl1d.o: $(TASKCOMM) libspl1d.f90
$(OBJDIR)/libspl2d.o: $(TASKCOMM) libspl2d.f90
$(OBJDIR)/libspl3d.o: $(TASKCOMM) libspl3d.f90
$(OBJDIR)/libsympl.o: $(TASKCOMM) libsympl.f90
$(OBJDIR)/libitp.o: $(TASKCOMM) libitp.f90 
$(OBJDIR)/libsmooth.o: $(TASKCOMM) libsmooth.f90 

$(OBJDIR)/netlib-bes.o: libspf.f90 netlib-bes.f90
$(OBJDIR)/libbes.o: $(TASKCOMM) netlib-bes.f90 libbes.f90

$(OBJDIR)/libchar.o: libchar.f90
$(OBJDIR)/libkio.o: libkio.f90 libchar.f90
$(OBJDIR)/libfio.o: libfio.f90

$(OBJDIR)/libufile/ufinit.o: libufile/ufinit.f90
$(OBJDIR)/libufile/uflist.o: libufile/uflist.f90 libufile/ufinit.f90
$(OBJDIR)/libufile/ufread.o: libufile/ufread.f90 libufile/ufinit.f90 \
                             libufile/uflist.f90
$(OBJDIR)/ufile.o: ufile.f90

$(OBJDIR)/libgrf/grftype.o: libgrf/grftype.f90
$(OBJDIR)/libgrf/grfutils.o: libgrf/grftype.f90 libgrf/grfutils.f90
$(OBJDIR)/libgrf/grdutils.o: libgrf/grftype.f90 libgrf/grdutils.f90
$(OBJDIR)/libgrf/grfconvert.o: libgrf/grftype.f90  libgrf/grfutils.f90 \
                     libgrf/grfconvert.f90
$(OBJDIR)/libgrf/grdconvert.o: $(TASKCOMM) libgrf/grftype.f90 \
                     libgrf/grfutils.f90 \
                     libgrf/grdconvert.f90
$(OBJDIR)/libgrf/grf1dexec.o: libgrf/grftype.f90 libgrf/grfutils.f90 \
                    libgrf/grf1dexec.f90
$(OBJDIR)/libgrf/grf2dexec.o: libgrf/grftype.f90 libgrf/grfutils.f90 \
                    libgrf/grf2dexec.f90
$(OBJDIR)/libgrf/grfxyexec.o: libgrf/grftype.f90 libgrf/grfutils.f90 \
                    libgrf/grfxyexec.f90
$(OBJDIR)/libgrf/grf1d.o: libgrf/grftype.f90 libgrf/grfconvert.f90 \
                libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 libgrf/grf1d.f90
$(OBJDIR)/libgrf/grd1d.o: $(TASKCOMM) libgrf/grftype.f90 \
                libgrf/grdconvert.f90 \
                libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 libgrf/grd1d.f90
$(OBJDIR)/libgrf/grf2d.o: libgrf/grftype.f90 libgrf/grfconvert.f90 \
                libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 libgrf/grf2d.f90
$(OBJDIR)/libgrf/grd2d.o: $(TASKCOMM) libgrf/grftype.f90 \
                libgrf/grdconvert.f90 \
                libgrf/grf1dexec.f90 libgrf/grf2dexec.f90 libgrf/grd2d.f90
$(OBJDIR)/libgrf/grdxy.o: $(TASKCOMM) libgrf/grftype.f90 \
                libgrf/grdconvert.f90 \
                libgrf/grfxyexec.f90 libgrf/grdxy.f90
$(OBJDIR)/libgrf/grd1dframe.o: $(TASKCOMM) libgrf/grftype.f90 \
                libgrf/grdconvert.f90 \
                libgrf/grd1dframe.f90
$(OBJDIR)/libgrf/grd2dframe.o: $(TASKCOMM) libgrf/grftype.f90 \
                libgrf/grdconvert.f90 \
                libgrf/grd2dframe.f90
$(OBJDIR)/libgrf.o: $(TASKCOMM) libgrf/grftype.f90 libgrf/grfutils.f90 \
          libgrf/grf1d.f90 \
          libgrf/grd1d.f90 libgrf/grf2d.f90 libgrf/grd2d.f90 libgrf/grdxy.f90 \
          libgrf/grd2dframe.f90 libgrf/grd1dframe.f90 libgrf.f90 \
          libgrf/grdutils.f90


