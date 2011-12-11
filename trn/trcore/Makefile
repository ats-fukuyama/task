TASKDIR=../../../task
BPSDDIR=../../../bpsd
include $(TASKDIR)/make.header
MOD = ../mod
MODDIR = -module ../mod

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a

FLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS) -warn

MODINCLUDE=-I $(BPSDDIR)/mod -I $(TASKDIR)/plx/mod -I $(TASKDIR)/lib/mod \
           -I ../trmodels/mod -I $(MOD) 


TRLIBS=	../trmodels/trmodelslib.a \
	$(TASKDIR)/plx/pllib.a \
        $(TASKDIR)/eq/eqlib.a \
        $(TASKDIR)/lib/tasklib.a \
	$(BPSDDIR)/bpsdlib.a

SRCS = trcomm.f90 trcalc.f90 \
       trexec.f90 trstep.f90 trresult.f90 trloop.f90 \
       trinit.f90 trsetup.f90 trgout.f90 \
       trmenu.f90 trmain.f90

OBJS = $(SRCS:.f90=.o)

.f90.o:
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

../trn: $(OBJS) $(TRLIBS)
	ifort $(FFLAGS) -o ../trn $(OBJS) $(TRLIBS) $(FLIBS)

$(TASKDIR)/plx/pllib.a:
	(cd $(TASKDIR)/plx; make pllib.a)
$(TASKDIR)/eq/eqlib.a:
	(cd $(TASKDIR)/eq; make eqlib.a)
$(TASKDIR)/lib/tasklib.a:
	(cd $(TASKDIR)/lib; make tasklib.a)
$(BPSDDIR)/bpsdlib.a:
	(cd $(BPSDDIR)/lib; make bpsdlib.a)

clean:
	rm -f *.o *~ *.mod
veryclean: clean
	rm -f trn

COMMON=trcomm.f90 $(TASKDIR)/plx/plcomm.f90
trcomm.o: trcomm.f90 $(TASKDIR)/plx/plcomm.f90
trinit.o: trinit.f90 $(COMMON)
trcalc.o: trcalc.f90 $(COMMON)
trexec.o: trexec.f90 $(COMMON)
trstep.o: trstep.f90 trexec.f90 trcalc.f90 $(COMMON)
trloop.o: trloop.f90 trstep.f90 $(COMMON)
trsetup.o: trsetup.f90 $(COMMON)
trgout.o: trgout.f90 $(COMMON)
trmenu.o: trmenu.f90 $(COMMON)
trmain.o: trmain.f90 $(COMMON)
