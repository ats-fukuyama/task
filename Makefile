### $Id$ ###
include make.header

all: message

message:
	@echo '## This is TASK main directory'
	@echo '     change working directory to'
	@echo '         pl, eq, tr, dp, wm, wr, fp or tot'
	@echo '     and make it'

clean:
	(cd lib; make clean)
	(cd mtx; make clean)
	(cd mtxp; make clean)
	(cd mpi; make clean)
	(cd eq; make clean)
	(cd equ; make clean)
	(cd pl; make clean)
	(cd plx; make clean)
	(cd dp; make clean)
	(cd dpx; make clean)
	(cd wr; make clean)
	(cd fp; make clean)
	(cd t2; make clean)
	(cd t2x; make clean)
	(cd tr; make clean)
	(cd trn; make clean)
	(cd tx; make clean)
	(cd w1; make clean)
	(cd w1n; make clean)
	(cd wf2; make clean)
	(cd wf2d; make clean)
	(cd wf2dt; make clean)
	(cd wf3; make clean)
	(cd wf3d; make clean)
	(cd wi; make clean)
	(cd wm; make clean)
	(cd wmf; make clean)
	(cd wmfn; make clean)
	(cd tot; make clean)
	(cd tools; make clean)
	rm -f core a.out *.o ./*~

veryclean: clean
	(cd lib; make veryclean)
	(cd mtx; make veryclean)
	(cd mtxp; make veryclean)
	(cd mpi; make veryclean)
	(cd eq; make veryclean)
	(cd equ; make veryclean)
	(cd pl; make veryclean)
	(cd plx; make veryclean)
	(cd dp; make veryclean)
	(cd dpx; make veryclean)
	(cd wr; make veryclean)
	(cd fp; make veryclean)
	(cd t2; make veryclean)
	(cd t2x; make veryclean)
	(cd tr; make veryclean)
	(cd trn; make veryclean)
	(cd tx; make veryclean)
	(cd w1; make veryclean)
	(cd w1n; make veryclean)
	(cd wf2; make veryclean)
	(cd wf2d; make veryclean)
	(cd wf2dt; make veryclean)
	(cd wf3; make veryclean)
	(cd wf3d; make veryclean)
	(cd wi; make veryclean)
	(cd wm; make veryclean)
	(cd wmf; make veryclean)
	(cd wmfn; make veryclean)
	(cd tot; make veryclean)
	(cd tools; make veryclean)

new:
	-mkdir ../tasknew
	cp Makefile ../tasknew
	cp make.header.org ../tasknew
	(cd lib; make new)
	mv libnew ../tasknew/lib
	(cd mtx; make new)
	mv mtxnew ../tasknew/mtx
	(cd mpi; make new)
	mv mpinew ../tasknew/mpi
	(cd eq; make new)
	mv eqnew ../tasknew/eq
	(cd pl; make new)
	mv plnew ../tasknew/pl
	(cd dp; make new)
	mv dpnew ../tasknew/dp
	(cd wr; make new)
	mv wrnew ../tasknew/wr
	(cd fp; make new)
	mv fpnew ../tasknew/fp
	(cd wm; make new)
	mv wmnew ../tasknew/wm
	(cd tr; make new)
	mv trnew ../tasknew/tr
	(cd tx; make new)
	mv txnew ../tasknew/tx
	(cd tot; make new)
	mv totnew ../tasknew/tot
	(cd tools; make new)
	mv toolsnew ../tasknew/tools
