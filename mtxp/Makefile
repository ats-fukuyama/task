include ${PETSC_DIR}/conf/base

testksp: libmtxksp.o testmtx.o 
	-${FLINKER} -g -o testksp testmtx.o libmtxksp.o ${PETSC_KSP_LIB}

testbnd: libmtxbnd.o testmtx.o 
	-${FLINKER} -g -o testbnd testmtx.o libmtxbnd.o ../task/lib/libbnd.o

ex1f: ex1f.o  chkopts
	-${FLINKER} -o ex1f ex1f.o  ${PETSC_KSP_LIB}
	${RM} ex1f.o

ex2f: ex2f.o  chkopts
	-${FLINKER} -o ex2f ex2f.o  ${PETSC_KSP_LIB}
	${RM} ex2f.o

runex1f:
	-@${MPIEXEC} -n 1 ./ex1f -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex1f_1.tmp 2>&1;	  \
	   if (${DIFF} output/ex1f_1.out ex1f_1.tmp) then true; \
	   else echo "Possible problem with ex1f_1, diffs above"; fi; \
	   ${RM} -f ex1f_1.tmp
runex2f:
	-@${MPIEXEC} -n 2 ./ex2f -pc_type jacobi -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex2f_1.tmp 2>&1; \
	   if (${DIFF} output/ex2f_1.out ex2f_1.tmp) then true; \
	   else echo "Possible problem with ex2f_1, diffs above"; fi; \
	   ${RM} -f ex2f_1.tmp
