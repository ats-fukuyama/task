### $Id$ ###
include ../make.header

FFLAGS=$(OFLAGS)

SRCS = plcnst.f90 plcom0.f90 plcom1.f90 plcom2.f90 plcom3.f90 \
       plmenu.f90 plinit.f90 plprof.f90 pldata.f90

OBJS = $(SRCS:.f90=.o) 

LIBS = ../lib/tasklib.a ../eq/eqlib.a pllib.a

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a
.f90.o:
	$(FCFREE) $(FFLAGS) -c $< -o $@ -I.

all : pl

../eq/eqlib.a:
	(cd ../lib; make eqlib.a)

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)

pllib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o pllib.a $(OBJS)

pl : $(LIBS) plmain.o
	$(FCFREE) plmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

check :
	ftnchek $(FCKFLAGS) $(SRCS) | less

clean:
	-rm -f core* a.out *.o *.mod ./*~ ./#* *.a

veryclean: clean
	-rm -f pl

new:
	-mkdir ../plnew
	cp -f Makefile ../plnew
	cp -f *.f90 ../plnew

plmain.o: plmain.f90 plcom0.f90 plcom1.f90
plmenu.o: plmenu.f90 plcom0.f90 plcom1.f90
plinit.o: plinit.f90 plcom0.f90 plcom1.f90 plcnst.f90
plprof.o: plprof.f90 plcom0.f90 plcom1.f90 plcom2.f90 plcnst.f90
pldata.o: pldata.f90 plcom0.f90 plcom3.f90
plcom0.o: plcom0.f90
plcom1.o: plcom1.f90
plcom2.o: plcom2.f90
plcom3.o: plcom3.f90
plcnst.o: plcnst.f90
