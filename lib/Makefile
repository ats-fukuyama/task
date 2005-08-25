### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)
MDSPLUS = mdsplus.f

SRCS = libdsp.f libsympl.f libspl.f libbrn.f libgut.f libbnd.f \
       libpcg.f libde.f libfp.f libbes.f librk.f libkio.f \
       libfft.f libinv.f libchar.f libfio.f libspl2.f libgrf.f ufile.f

OBJS = $(SRCS:.f=.o)
OBJMDS = $(MDSPLUS:.f=.o)

.f.o :
	$(FC) $(LFLAGS) -c $< -o $@

all : tasklib.a mdslib.a

tasklib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o tasklib.a $(OBJS)

mdslib.a: $(OBJMDS) mdslib.inc
	$(LD) $(LDFLAGS) -o mdslib.a $(OBJMDS)

testspl : testspl.o libspl.o ../grf/libgrf.o
	$(FC) -o $@ testspl.o $(FFLAGS) libspl.o ../grf/libgrf.o $(FLIBS)

testgrf : testgrf.o $(OBJS)
	$(FC) -o $@ testgrf.o $(FFLAGS) $(OBJS) $(FLIBS)

new :
	mkdir ../libnew
	cp -f Makefile ../libnew
	cp -f *.f ../libnew

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core a.out *.o *.a ./*~ ./#*

veryclean: clean

mdsplus.o  : mdsplus.f mdslib.inc
