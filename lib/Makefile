### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

LFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE= -I.

SRCSFREE  =  libbrn.f90 libinv.f90 libbnd.f90 \
             librk.f90 libsympl.f90 \
             libspl.f90 libeq.f90 \
             libgut.f90 \
             libchar.f90 libkio.f90 \
             libfio.f90 ufile.f90 libitp.f90 libfem.f90
SRCSFIXED =  libspf.f libbes.f libell.f \
             libde.f libfft.f libdsp.f \
             libfp.f libgrf.f $(LAPACK)

OBJSFREE  = $(SRCSFREE:.f90=.o)
OBJSFIXED = $(SRCSFIXED:.f=.o)
OBJSMDS = $(MDSPLUS:.f=.o)

.f.o :
	$(FCFIXED) $(LFLAGS) -c $< -o $@ $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(LFLAGS) -c $< -o $@ $(MODINCLUDE)

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

femtest: tasklib.a femtest.o 
	$(FCFIXED) femtest.o tasklib.a $(FFLAGS) -o $@ $(FLIBS) $(MODINCLUDE)

mmconvert: mmconvert.o
	$(FCFIXED) mmconvert.o $(FFLAGS) -o $@

#libspl.o: libspl.f90
#	$(FCFREE) -O3 -tp p7 -c libspl.f90 -o libspl.o
