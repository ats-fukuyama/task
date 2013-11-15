### $Id$###
include ../make.header

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

MODINCLUDE= -I $(MOD) -I../lib/$(MOD)

SRCS = wiunmag.f

OBJS = $(SRCS:.f=.o) 

LIBS = ../lib/tasklib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : wiunmag

#dplib.a: $(OBJS)
#	$(LD) $(LDFLAGS) dplib.a $(OBJS)

libs:
	(cd ../lib; make tasklib.a)
	(cd ../../bpsd; make bpsdlib.a)

wiunmag : wiunmag.o
	$(FCFIXED) wiunmag.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS)

clean:
	-rm -f core a.out *.o *.a ./*~

veryclean: clean
	-rm -f wiunmag

wiunmag.o: wiunmag.f
