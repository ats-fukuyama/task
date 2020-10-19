### Madefile for cv
include ../make.header
include ../mtxp/make.mtxp

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../lib/$(MOD)

SRCS = cvcomm.f90 cvlib.f90 libcharx.f90 \
       cvinit.f90 cvparm.f90 cvview.f90 cvfile.f90 \
       cvpopulation.f90 cvrank.f90 \
       cvgout.f90 cvregion.f90 cvselect.f90 \
       cvmenu.f90

OBJS = $(SRCS:.f90=.o)
FLINKER=$(FCFREE)

LIBS = libcv.a  ../lib/libgrf.a ../lib/libtask.a

all : libs cv

libcv.a: $(OBJS)
	$(LD) $(LDFLAGS) libcv.a $(OBJS)

libs:
	(cd ../lib; make libgrf.a)
	(cd ../lib; make libtask.a)

cv : libs libcv.a cvmain.o
	$(FLINKER) cvmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

clean:
	-rm -f core a.out *.o *.a ./*~ $(MOD)/*.mod

veryclean: clean
	-rm -f cv

CVCOMM=cvcomm.f90
cvmain.o: cvmain.f90 cvinit.f90 cvparm.f90 cvmenu.f90 $(CVCOMM)
cvmenu.o: cvmenu.f90 cvparm.f90 cvview.f90 cvgout.f90 cvregion.f90 \
          cvselect.f90 cvpopulation.f90 cvlib.f90 $(CVCOMM)
cvinit.o: cvinit.f90 $(CVCOMM)
cvparm.o: cvparm.f90 $(CVCOMM)
cvview.o: cvview.f90 $(CVCOMM)
cvfile.o: cvfile.f90 $(CVCOMM)
cvgout.o: cvgout.f90 cvrank.f90 $(CVCOMM)
cvregion.o: cvregion.f90 $(CVCOMM)
cvselect.o: cvselect.f90 $(CVCOMM)
cvpopulation.o: cvpopulation.f90 $(CVCOMM)
cvrank.o: cvrank.f90 $(CVCOMM)
cvlib.o:  cvlib.f90
