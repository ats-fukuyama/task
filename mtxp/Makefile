FFLAGS=-g

include ${PETSC_DIR}/conf/base

testksp: libmtxksp.o testmtx.o 
	-${FLINKER} -g -o testksp testmtx.o libmtxksp.o ${PETSC_KSP_LIB}

testbnd: libmtxbnd.o testmtx.o 
	-${FLINKER} -g -o testbnd testmtx.o libmtxbnd.o ../lib/libbnd.o

testmtx.o: testmtx.f90
libmtxksp.o: libmtxksp.F
libmtxbnd.o: libmtxbnd.f90
