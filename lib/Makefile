### $Id$ ###
include ../make.header

LFLAGS = $(OFLAGS)

SRCS = libdsp.f libode.f librkf.f libspl.f libbrn.f libgut.f libbnd.f \
       libpcg.f libde.f libfp.f libbes.f librk.f \
       libfft.f libinv.f

OBJS = $(SRCS:.f=.o)

.f.o :
	$(FC) $(LFLAGS) -c $< -o $@

all : $(OBJS)

new :
	mkdir ../libnew
	cp -f Makefile ../libnew
	cp -f *.f ../libnew

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core a.out *.o ./*~ ./#*

veryclean: clean

