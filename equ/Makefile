include ../make.header

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)
#FFLAGS = $(DFLAGS) -warn all


.SUFFIXES:
.SUFFIXES: .f .mod .o .a

SRCCOM =  eqcomm.f

SRCTPX = tpxssl.f 

SRCEQ = eqpfds.f eqsub.f eqinit.f \
        eqfct.f  eqinp.f eqset.f \
        eqpl.f eqtrn.f eqfile.f \
	eqgout.f equunit.f

SRCTR = trpl.f trset.f tradv.f trunit.f

#
BPSDLIB= ../../bpsd/libbpsd.a
OBJCOM = $(SRCCOM:.f=.o)
OBJTPX= $(SRCTPX:.f=.o)
OBJEQ= $(SRCEQ:.f=.o)
OBJTR= $(SRCTR:.f=.o)
LIBS = ../lib/libtask.a $(BPSDLIB) $(OBJCOM)  libtpx.a libeq.a libtr.a
#
MODINCLUDE=-I$(MOD) -I../../bpsd/$(MOD) -I../lib/$(MOD)
#
#
.f.o:
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o:
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all: $(LIBS) equ

$(LIBBPSD):
	(cd ../../bpsd; make libbpsd.a)

libtpx.a: $(OBJTPX)
	$(LD) $(LDFLAGS) $@ $(OBJTPX)
libeq.a: $(OBJEQ)
	$(LD) $(LDFLAGS) $@ $(OBJEQ)
libtr.a: $(OBJTR)
	$(LD) $(LDFLAGS) $@ $(OBJTR)
equ: $(LIBS) eqmenu.o eqmain.o
	$(FCFIXED) eqmain.o eqmenu.o $(LIBS) $(FLIBS) -o $@ $(FFLAGS)
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
