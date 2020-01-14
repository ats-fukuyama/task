include ../make.header

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)
#FFLAGS = $(DFLAGS) -warn all


.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

SRCCOM =  eqcomm.f
SRCS1 =  tpxssl.f eqpfds.f eqsub.f eqinit.f \
         eqfct.f  eqinp.f eqset.f \
         eqpl.f eqtrn.f eqfile.f \
	 eqgout.f equunit.f \
	 trpl.f trset.f tradv.f trunit.f
SRCSX =  equ-eqdsk.f90
SRCS2 = eqmenu.f


BPSDLIB= ../../bpsd/libbpsd.a
OBJCOM = $(SRCCOM:.f=.o)
OBJS = $(SRCS1:.f=.o) $(SRCSX:.f90=.o) $(SRCS2:.f=.o)
LIBS = libequ.a ../pl/libpl.a ../lib/libtask.a $(BPSDLIB) 

MODINCLUDE= -I../../bpsd/$(MOD) -I../lib/$(MOD) -I../pl/$(MOD) -I$(MOD)


.f.o:
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o:
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: equ

libequ.a: $(OBJS)
	$(LD) $(LDFLAGS) libequ.a $(OBJS)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib/; make libtask.a)
	(cd ../pl/; make libpl.a)

equ: $(OBJCOM) libs libequ.a eqmain.o
	$(FCFIXED) eqmain.o $(OBJCOM) $(LIBS) -o $@ $(FFLAGS) $(FLIBS) 
clean:
	-rm -f *.a *.o $(MOD)/*.mod *~
veryclean: clean
	-rm -f equ
#
eqmain.o : eqmain.f
eqmenu.o : eqmenu.f
equunit.o : equunit.f eqcomm.f
eqfct.o : eqfct.f eqcomm.f
eqinp.o : eqinp.f eqcomm.f
eqinit.o : eqinit.f eqcomm.f
eqset.o : eqset.f eqcomm.f
eqsub.o : eqsub.f eqcomm.f
eqpfds.o : eqpfds.f eqcomm.f
eqintf.o : eqintf.f eqcomm.f
eqgout.o : eqgout.f eqcomm.f
eqfile.o : eqgout.f eqcomm.f
trunit.o : trunit.f  
trpl.o : trpl.f  
trset.o : trset.f  
tradv.o : tradv.f trpl.f
tpxssl.o : tpxssl.f
equ-eqdsk.o : equ-eqdsk.f90 eqcomm.f
