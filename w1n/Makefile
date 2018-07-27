### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f90 .f .mod .o .a

FFLAGS = $(OFLAGS)
#FFLAGS = $(DFLAGS)

MODINCLUDE=-I $(MOD) -I../../bpsd/$(MOD) -I../lib/$(MOD)

.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

SRCSFREE =w1comm.f90 w1lib.f90 w1init.f90 w1parm.f90 \
          w1prep.f90 w1prof.f90 w1bcnd.f90 w1rslt.f90 w1out.f90 \
          w1alfa.f90 w1disp.f90 w1grf1.f90 w1grf2.f90 \
          w1mlm.f90 w1fem.f90 w1intg.f90 w1gdsp.f90\
          w1exec.f90 w1gout.f90 \
          w1menu.f90 w1main.f90


#OBJS = $(SRCSFIXED:.f=.o) $(SRCSFREE:.f90=.o) 
OBJS = $(SRCSFREE:.f90=.o) $(SRCSFIXED:.f=.o) 
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

w1main.o: w1comm.f90 w1main.f90
w1init.o: w1comm.f90 w1init.f90
w1parm.o: w1comm.f90 w1parm.f90
w1exec.o: w1comm.f90 w1exec.f90
w1prep.o: w1comm.f90 w1prep.f90
w1prof.o: w1comm.f90 w1prof.f90
w1bcnd.o: w1comm.f90 w1bcnd.f90
w1rslt.o: w1comm.f90 w1rslt.f90
w1out.o:  w1comm.f90 w1out.f90
w1grf1.o: w1comm.f90 w1grf1.f90
w1grf2.o: w1comm.f90 w1grf2.f90
w1disp.o: w1comm.f90 w1disp.f90
w1gdsp.o: w1comm.f90 w1gdsp.f90
w1mlm.o:  w1comm.f90 w1mlm.f90
w1fem.o:  w1comm.f90 w1fem.f90
w1intg.o: w1comm.f90 w1intg.f90
w1lib.o: w1lib.f90
