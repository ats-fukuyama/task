### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)
FFLAGS = $(OFLAGS)

SRCS = libspf.f libbes.f libell.f \
       libbrn.f \
       libinv.f libbnd.f \
       librk.f libsympl.f libde.f \
       libfft.f libspl.f \
       libdsp.f \
       libfp.f libeq.f \
       libgut.f libgrf.f \
       libchar.f libkio.f \
       libfio.f ufile.f \
       libitp.f \
       $(LAPACK)

OBJS = $(SRCS:.f=.o)
OBJMDS = $(MDSPLUS:.f=.o)

.f.o :
	$(FCFIXED) $(LFLAGS) -c $< -o $@
.f90.o :
	$(FCFREE) $(LFLAGS) -c $< -o $@

all : tasklib.a mdslib.a

tasklib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o tasklib.a $(OBJS)

mdslib.a: $(OBJMDS) mdslib.inc
	$(LD) $(LDFLAGS) -o mdslib.a $(OBJMDS)

testspl : testspl.o libspl.o ../grf/libgrf.o
	$(FCFIXED) -o $@ testspl.o $(FFLAGS) libspl.o ../grf/libgrf.o $(FLIBS)

testgrf : testgrf.o $(OBJS)
	$(FCFIXED) -o $@ testgrf.o $(FFLAGS) $(OBJS) $(FLIBS)

testbes : testbes.o $(OBJS)
	$(FCFIXED) -o $@ testbes.o $(FFLAGS) $(OBJS) $(FLIBS)

testfft : testfft.o $(OBJS)
	$(FCFIXED) -o $@ testfft.o $(FFLAGS) $(OBJS) $(FLIBS)

testbrn : testbrn.o $(OBJS)
	$(FCFIXED) -o $@ testbrn.o $(FFLAGS) $(OBJS) $(FLIBS)

testell : testell.o $(OBJS)
	$(FCFIXED) -o $@ testell.o $(FFLAGS) $(OBJS) $(FLIBS)

new :
	mkdir ../libnew
	cp -f Makefile ../libnew
	cp -f *.f ../libnew

check :
	ftnchek $(SRCS) $(MDSPLUS) | less

clean:
	-rm -f core a.out *.o *.a ./*~ ./#*

veryclean: clean

mdsplus.o  : mdsplus.f mdslib.inc
