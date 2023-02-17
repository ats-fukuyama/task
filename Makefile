### $Id$ ###
include make.header

all: message

message:
	@echo '## This is TASK main directory'
	@echo '     change working directory to'
	@echo '         pl, eq, tr, dp, wm, wr, fp or tot'
	@echo '     and make it'

clean:
	(cd adpost; make clean)
	(cd dp; make clean)
	(cd eq; make clean)
	(cd equ; make clean)
	(cd fit3d; make clean)
	(cd fp; make clean)
	(cd lib; make clean)
	(cd mtxp; make clean)
	(cd ob; make clean)
	(cd open-adas/adf11/adf11-lib; make clean)
	(cd pic; make clean)
	(cd pl; make clean)
	(cd ti; make clean)
	(cd tools; make clean)
	(cd tot; make clean)
	(cd tr; make clean)
	(cd tx; make clean)
	(cd w1; make clean)
	(cd wf2; make clean)
	(cd wf2d; make clean)
	(cd wf2dt; make clean)
	(cd wf3; make clean)
	(cd wf3d; make clean)
	(cd wi; make clean)
	(cd wim; make clean)
	(cd wm; make clean)
	(cd wmf; make clean)
	(cd wmfn; make clean)
	(cd wmx; make clean)
	(cd wr; make clean)
	(cd wq; make clean)

	rm -f core a.out *.o ./*~

veryclean: clean
	(cd adpost; make veryclean)
	(cd dp; make veryclean)
	(cd eq; make veryclean)
	(cd equ; make veryclean)
	(cd fit3d; make veryclean)
	(cd fp; make veryclean)
	(cd lib; make veryclean)
	(cd mtxp; make veryclean)
	(cd ob; make veryclean)
	(cd open-adas/adf11/adf11-lib; make veryclean)
	(cd pic; make veryclean)
	(cd pl; make veryclean)
	(cd ti; make veryclean)
	(cd tools; make veryclean)
	(cd tot; make veryclean)
	(cd tr; make veryclean)
	(cd tx; make veryclean)
	(cd w1; make veryclean)
	(cd wf2; make veryclean)
	(cd wf2d; make veryclean)
	(cd wf2dt; make veryclean)
	(cd wf3; make veryclean)
	(cd wf3d; make veryclean)
	(cd wi; make veryclean)
	(cd wim; make veryclean)
	(cd wm; make veryclean)
	(cd wmf; make veryclean)
	(cd wmfn; make veryclean)
	(cd wmx; make veryclean)
	(cd wr; make veryclean)
	(cd wq; make veryclean)
