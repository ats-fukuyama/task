### $Id$ ###
include ../make.header

LFLAGS = $(DFLAGS)

SRCSFREE  =  libbrn.f90 libinv.f90 libbnd.f90 \
             librk.f90 libsympl.f90 \
             libspl.f90 libeq.f90 \
             libgut.f90 \
             libchar.f90 libkio.f90 \
             libfio.f90 ufile.f90 libitp.f90
SRCSFIXED =  libspf.f libbes.f libell.f \
             libde.f libfft.f libdsp.f \
             libfp.f libgrf.f $(LAPACK)

OBJSFREE  = $(SRCSFREE:.f90=.o)
OBJSFIXED = $(SRCSFIXED:.f=.o)
OBJSMDS = $(MDSPLUS:.f=.o)

.SUFFIXES:
.SUFFIXES: .f .f90 .o .a

.f.o :
	$(FCFIXED) $(LFLAGS) -c $< -o $@
.f90.o :
	$(FCFREE) $(LFLAGS) -c $< -o $@

all : tasklib.a mdslib.a

tasklib.a: $(OBJSFREE) $(OBJSFIXED)
	$(LD) $(LDFLAGS) -o tasklib.a $(OBJSFREE) $(OBJSFIXED)

mdslib.a: $(OBJSMDS)
	$(LD) $(LDFLAGS) -o mdslib.a $(OBJSMDS)

new :
	mkdir ../libnew
	cp -f Makefile ../libnew
	cp -f *.f90 ../libnew
	cp -f *.f ../libnew

check :
	ftnchek $(SRCSFREE) $(SRCSFIXED) $(MDSPLUS) | less

clean:
	-rm -f core* a.out *.o *.mod *.a ./*~ ./#*

veryclean: clean

libspl.o   : libspl.f90
libsympl.o : libsympl.f90
librk.o    : librk.f90
libbrn.o   : libbrn.f90
libfio.o   : libfio.f90
libchar.o  : libchar.f90
libkio.o   : libkio.f90
libgut.o   : libgut.f90
libbnd.o   : libbnd.f90
libeq.o    : libeq.f90
libitp.o   : libitp.f90
ufile.o    : ufile.f90
