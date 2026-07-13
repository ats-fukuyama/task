### $Id$ ###
include ../make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

OBJDIR=./obj

MODINCLUDE=-I $(MOD) -I../../bpsd/$(MOD) -I../lib/$(MOD)

${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

SRCS = wimcomm.f90 wiminit.f90 wimparm.f90 libpdkf.f90 \
       wimcalc.f90 wimexec.f90 wimgout.f90 \
       wimmenu.f90

OBJS = $(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))

LIBS = libwim.a ../lib/libgrf.a ../lib/libtask.a ../../bpsd/libbpsd.a

all : libs wim

libwim.a: $(OBJS)
	$(LD) $(LDFLAGS) libwim.a $(OBJS)
../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libtask.a)
	(cd ../lib; make libgrf.a)

wim : $(LIBS) $(OBJDIR)/wimmain.o
	$(FCFREE) $(OBJDIR)/wimmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

wimgconv: $(LIBS) wimgconv.o
	$(FCREE)/wimgconv.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

clean:
	-rm -f core ./#* ./*~ $(OBJDIR)/*.o $(MOD)/*.mod coew *.a

veryclean: clean
	-rm -f wim

$(OBJDIR)/wimmain.o: wimmain.f90 wimcomm.f90
$(OBJDIR)/wimmenu.o: wimmenu.f90 wimcomm.f90
$(OBJDIR)/wiminit.o: wiminit.f90 wimcomm.f90
$(OBJDIR)/wimparm.o: wimparm.f90 wimcomm.f90
$(OBJDIR)/libpdkf.o: libpdkf.f90
$(OBJDIR)/wimcalc.o: wimcalc.f90 wimcomm.f90
$(OBJDIR)/wimexec.o: wimexec.f90 wimcomm.f90
$(OBJDIR)/wimgout.o: wimgout.f90 wimcomm.f90
