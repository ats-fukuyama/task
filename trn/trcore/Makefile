TASKDIR=../../../task
BPSDDIR=../../../bpsd
include $(TASKDIR)/make.header

.SUFFIXES:
.SUFFIXES: .f90 .mod .o .a


#FFLAGS=$(OFLAGS)
FLAGS=$(DFLAGS) -g


MODINCLUDE=-I $(BPSDDIR)/mod -I $(TASKDIR)/plx/mod -I $(TASKDIR)/lib/mod \
	   -I $(TASKDIR)/eq/mod -I ../trgraphics/mod -I ../trmodels/mod  \
	   -I ../trsources/mod  -I ../trufile/mod -I $(MOD) 

TRLIBS=	../trmodels/trmodelslib.a \
	../trsources/trsourceslib.a \
	../trgraphics/trgraphicslib.a \
	../trufile/trufilelib.a \
	$(TASKDIR)/plx/pllib.a \
        $(TASKDIR)/eq/eqlib.a \
        $(TASKDIR)/lib/tasklib.a \
        $(TASKDIR)/lib/mdslib.a \
	$(BPSDDIR)/bpsdlib.a

SRCS = trcomm.f90 trbpsd.f90  trlib.f90   trcalv.f90   trcalc2.f90 \
       trexec.f90 trstep.f90  trcalc1.f90 trresult.f90 trwrite.f90 \
       trloop.f90 trinit.f90 trsetup.f90 trmenu.f90  trmain.f90 

OBJS = $(SRCS:.f90=.o)


.f90.o:
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)


../trn: $(OBJS) $(TRLIBS)
	$(FCFREE) $(FFLAGS) -o ../trn $(OBJS) $(TRLIBS) $(FLIBS)


$(TASKDIR)/plx/pllib.a:
	(cd $(TASKDIR)/plx; make pllib.a)
$(TASKDIR)/eq/eqlib.a:
	(cd $(TASKDIR)/eq; make eqlib.a)
$(TASKDIR)/lib/tasklib.a:
	(cd $(TASKDIR)/lib; make tasklib.a)
$(TASKDIR)/lib/mdslib.a:
	(cd $(TASKDIR)/lib; make mdslib.a)
$(BPSDDIR)/bpsdlib.a:
	(cd $(BPSDDIR)/lib; make bpsdlib.a)

clean:
	rm -f *.o *.f90~ *.mod $(MOD)/*.mod

veryclean:
	rm -f *.o *.~ *.mod $(MOD)/*.mod

COMMON=trcomm.f90 $(TASKDIR)/plx/plcomm.f90
trcomm.o: trcomm.f90 $(TASKDIR)/plx/plcomm.f90
trbpsd.o: trbpsd.f90 $(COMMON)
trlib.o: trlib.f90 $(COMMON)
trinit.o: trinit.f90 $(COMMON)
trcalv.o: trcalv.f90 $(COMMON)
trexec.o: trexec.f90 $(COMMON)
trstep.o: trstep.f90 trexec.f90 trcalc2.f90 $(COMMON)
trloop.o: trloop.f90 trcalc1.f90 trstep.f90 trbpsd.f90 $(COMMON)
trsetup.o: trsetup.f90 trbpsd.f90 $(COMMON)
trmenu.o: trmenu.f90 $(COMMON)
trmain.o: trmain.f90 $(COMMON)
