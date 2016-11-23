### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE=-I $(MOD) -I../../bpsd/$(MOD) -I../lib/$(MOD)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

SRCSFIXED=w1init.f w1prep.f w1prof.f w1bcnd.f w1rslt.f w1out.f \
          w1grf1.f w1grf2.f w1disp.f w1mlm.f w1fem.f w1intg.f w1lib.f
SRCSFREE = 

OBJS = $(SRCSFIXED:.f=.o) $(SRCSFREE:.f90=.o) 
LIBS = w1lib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a

all: w1

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
w1lib.a: $(OBJS)
	$(LD) $(LDFLAGS) w1lib.a $(OBJS)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)

w1 : $(LIBS) w1main.o
	$(FCFREE) w1main.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) \
           $(MODDIR) $(MODINCLUDE)

clean:
	-rm -f core ./*.o *.mod ./*~ ./\#* *.a $(MOD)/*.mod
veryclean: clean
	-rm -f w1 core

w1main.o: w1comm.f w1main.f
w1init.o: w1comm.f w1init.f
w1prep.o: w1comm.f w1prep.f
w1prof.o: w1comm.f w1prof.f
w1bcnd.o: w1comm.f w1bcnd.f
w1rslt.o: w1comm.f w1rslt.f
w1out.o:  w1comm.f w1out.f
w1grf1.o: w1comm.f w1grf1.f
w1grf2.o: w1comm.f w1grf2.f
w1disp.o: w1comm.f w1disp.f
w1mlm.o:  w1comm.f w1mlm.f
w1fem.o:  w1comm.f w1fem.f
w1intg.o: w1comm.f w1intg.f
w1lib.o: w1lib.f
