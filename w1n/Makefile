### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FLFLAGS = $(OFLAGS)
FLFLAGS = $(DFLAGS)
#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE= -I $(MOD) -I../../bpsd/$(MOD) -I../lib/$(MOD)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

SRCSFIXED= w1lib.f w1prep.f w1prof.f \
           w1disp.f w1mlm.f w1fem.f w1intg.f w1bcnd.f \
           w1rslt.f w1out.f w1grf1.f w1grf2.f
SRCSFREE1=  w1comm.f90
SRCSFREE2=  w1gout.f90 w1gdsp.f90 w1init.f90 w1exec.f90 w1parm.f90 \
            w1init.f90 w1menu.f90 

OBJS = $(SRCSFREE1:.f90=.o) $(SRCSFIXED:.f=.o) $(SRCSFREE2:.f90=.o) 

LIBS = w1lib.a ../lib/tasklib.a ../../bpsd/bpsdlib.a

all : w1

../lib/tasklib.a:
	(cd ../lib; make tasklib.a)
../../bpsd/bpsdlib.a:
	(cd ../../bpsd; make bpsdlib.a)
w1lib.a: $(OBJS)
	$(LD) $(LDFLAGS) w1lib.a $(OBJS)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)

w1:$(LIBS) w1main.o
	$(FCFREE) w1main.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) \
           $(MODDIR) $(MODINCLUDE)

clean:
	-rm -f core a.out *.o *.mod ./*~ ./#* *.a $(MOD)/*.mod
veryclean: clean
	-rm -f w1

w1main.o: w1comm.f90 w1main.f90
w1init.o: w1comm.f90 w1init.f90
w1menu.o: w1comm.f90 w1menu.f90
w1parm.o: w1comm.f90 w1parm.f90
w1exec.o: w1comm.f90 w1exec.f90
w1gout.o: w1comm.f90 w1gout.f90
w1prep.o: w1comm.f90 w1prep.f
w1prof.o: w1comm.f90 w1prof.f
w1bcnd.o: w1comm.f90 w1bcnd.f
w1rslt.o: w1comm.f90 w1rslt.f
w1out.o:  w1comm.f90 w1out.f
w1grf1.o: w1comm.f90 w1grf1.f
w1grf2.o: w1comm.f90 w1grf2.f
w1disp.o: w1comm.f90 w1disp.f
w1mlm.o:  w1comm.f90 w1mlm.f
w1fem.o:  w1comm.f90 w1fem.f
w1intg.o: w1comm.f90 w1intg.f
w1lib.o: w1lib.f
