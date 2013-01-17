### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FLFLAGS = $(OFLAGS)
FLFLAGS = $(DFLAGS)
#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

MODINCLUDE= -I./$(MOD) -I../lib/$(MOD)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

SOURCES= w1main.f w1init.f w1prep.f w1prof.f w1bcnd.f w1rslt.f w1out.f \
         w1grf1.f w1grf2.f w1disp.f w1mlm.f w1fem.f w1intg.f w1lib.f

OBJS=$(SOURCES:.f=.o)
LIBS=../lib/tasklib.a

w1:$(OBJS)
	$(FCFIXED) $(OBJS) $(LIBS) $(FFLAGS) $(FLIBS) -o $@

check:
	ftnchek $(SOURCES) | less

clean:
	-rm -f core ./*.o ./*~ ./\#*
veryclean:
	-rm -f w1 core ./*.o ./*~ ./\#*

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
