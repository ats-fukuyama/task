FFLAGS=-g

include ${PETSC_DIR}/conf/base

all: testbnd testksp testbndpara testksppara
libs: libmtxbnd.o libmtxksp.o libmtxksppara.o

testksp: libmtxksp.o testmtx.o 
	-${FLINKER} -g -o testksp testmtx.o libmtxksp.o ${PETSC_KSP_LIB}

testbnd: libmtxbnd.o testmtx.o 
	-${FLINKER} -g -o testbnd testmtx.o libmtxbnd.o ../lib/libbnd.o

testksppara: libmtxksppara.o testmtxpara.o 
	-${FLINKER} -g -o testksppara testmtxpara.o libmtxksppara.o \
         ${PETSC_KSP_LIB}

testbndpara: libmtxbnd.o testmtxpara.o 
	-${FLINKER} -g -o testbndpara testmtxpara.o libmtxbnd.o ../lib/libbnd.o

testmtx.o: testmtx.f90
testmtxpara.o: testmtxpara.f90
libmtxbnd.o: libmtxbnd.f90
libmtxksp.o: libmtxksp.F
libmtxksppara.o: libmtxksppara.F
