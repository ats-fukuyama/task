FFLAGS=-g -I${INCLUDE_MUMPS}

include ${PETSC_DIR}/conf/base
INCLUDE_MUMPS=${PETSC_DIR}/externalpackages/MUMPS_4.8.4/include

all: testbnd testksp testdmumps
libs: libmtxbnd.o libmtxksp.o

testksp: libmtxksp.o testmtx.o 
	-${FLINKER} -g -o testksp testmtx.o libmtxksp.o ${PETSC_KSP_LIB}

testdmumps: libmtxdmumps.o testmtx.o 
	-${FLINKER} -g -o testdmumps testmtx.o libmtxdmumps.o ${PETSC_KSP_LIB}

testbnd: libmtxbnd.o testmtx.o 
	-${FLINKER} -g -o testbnd testmtx.o libmtxbnd.o ../../task/lib/libbnd.o

testmtx.o: testmtx.f90
libmtxbnd.o: libmtxbnd.f90
libmtxdmumps.o: libmtxdmumps.F
libmtxksp.o: libmtxksp.F
