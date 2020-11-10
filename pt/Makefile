### plot/Makefile ###
include ../make.header
include ../mtxp/make.mtxp

LIB_MTX=$(LIB_MTX_BND)
LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_PCG)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIB_MTX=$(LIBXMTX_KSP)
#LIBX_MTX=$(LIBX_MTX_KSP)

#FFLAGS = $(OFLAGS)
FFLAGS = $(DFLAGS)

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) \
            -I../mtxp/$(MOD) -I../lib/$(MOD) -I../../bpsd/$(MOD) 

SRCS=   ptcomm.f90 ptinit.f90 ptparm.f90 ptview.f90 \
        ptexec1.f90 ptmenu.f90 

OBJS=$(SRCS:.f90=.o)

LIBS=libpt.a ../lib/libtask.a ../lib/libgrf.a ../../bpsd/libbpsd.a


all: pt

libpt.a: $(OBJS)
	$(LD) $(LDFLAGS) $@ $(OBJS)

../lib/libtask.a:
	(cd ../lib; make libtask.a)
../lib/libgrf.a:
	(cd ../lib; make libgrf.a)
../../bpsd/libbpsd.a:
	(cd ../../bpsd; make libbpsd.a)

libs:
	(cd ../../bpsd; make libbpsd.a)
	(cd ../lib; make libgrf.a)
	(cd ../lib; make libtask.a)
	(cd ../mtxp; make)

pt: libs libpt.a ptmain.o 
	$(FLINKER) ptmain.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(LIBX_MTX)

clean : 
	rm -f ./#* ./*~ *.o core *.a *.mod $(MOD)/*.mod

veryclean : clean
	rm -f pt

ptcomm.o : ptcomm.f90
ptexec1.o: ptcomm.f90 ptexec1.f90
ptview.o : ptcomm.f90 ptview.f90
ptparm.o : ptcomm.f90 ptparm.f90
ptinit.o : ptcomm.f90 ptinit.f90
ptmenu.o : ptcomm.f90 ptexec1.f90 ptparm.f90 ptview.f90 ptmenu.f90
ptmain.o : ptcomm.f90 ptinit.f90 ptparm.f90 ptmenu.f90 ptmain.f90
