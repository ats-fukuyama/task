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
          w1sub.f90 w1fflr.f90 w1prof.f90 w1out.f90 \
	  w1prep.f90 w1prepost.f90 w1prepostx.f90 \
          w1bcnd.f90 w1rslt.f90  \
          w1alfa.f90 w1disp.f90 w1grf1.f90 w1grf2.f90 \
          w1mlm.f90 w1fem.f90 w1intg.f90 w1gdsp.f90 \
          w1exec6.f90 w1exec7.f90 w1exec8.f90 \
          w1exec9.f90 w1exec10.f90 w1exec11.f90 \
          w1exec.f90 w1gout.f90 \
          w1menu.f90 w1main.f90


#OBJS = $(SRCSFIXED:.f=.o) $(SRCSFREE:.f90=.o) 
OBJS = $(SRCSFREE:.f90=.o) $(SRCSFIXED:.f=.o) 
LIBS = libw1.a ../lib/libtask.a ../../bpsd/libbpsd.a

all: w1

../lib/libtask.a:
	(cd ../lib; make libtask.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)
libw1.a: $(OBJS)
	$(LD) $(LDFLAGS) libw1.a $(OBJS)

libs: ../../bpsd/libbpsd.a ../lib/libtask.a libw1.a
 
w1 : libs w1main.o
	$(FCFREE) w1main.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) \
           $(MODDIR) $(MODINCLUDE)

clean:
	-rm -f core ./*.o *.mod ./*~ ./\#* *.a $(MOD)/*.mod
veryclean: clean
	-rm -f w1 core

w1main.o: w1comm.f90 w1main.f90
w1menu.o: w1comm.f90 w1menu.f90
w1init.o: w1comm.f90 w1init.f90
w1parm.o: w1comm.f90 w1parm.f90
w1exec.o: w1comm.f90 w1exec.f90
w1exec6.o: w1comm.f90 w1exec6.f90
w1exec7.o: w1comm.f90 w1exec7.f90
w1exec8.o: w1comm.f90 w1exec8.f90
w1exec9.o: w1comm.f90 w1exec9.f90
w1exec10.o: w1comm.f90 w1exec10.f90
w1exec11.o: w1comm.f90 w1exec11.f90
w1prepost.o: w1comm.f90 w1prepost.f90
w1prepostx.o: w1comm.f90 w1prepostx.f90
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
w1qtbl.o: w1comm.f90 w1qtbl.f90
w1sub.o: w1comm.f90 w1sub.f90
w1lib.o: w1lib.f90
