
include ${PETSC_DIR}/conf/base

all: testksp

testksp: libmtxksp.o testksp.o 
	-${FLINKER} -g -o testksp testksp.o libmtxksp.o ${PETSC_KSP_LIB}

veryclean:
	rm *.o *.mod testksp

libmtxksp.o: libmtxksp.F
testksp.o: testksp.f libmtxksp.F
