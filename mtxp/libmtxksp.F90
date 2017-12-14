!     $Id$

      MODULE libmtx

!  for PETSc V3.8
#include "petsc/finclude/petsc.h"
      use petsc
!
      use libmpi
      use commpi
      PRIVATE

      PUBLIC mtx_initialize
      PUBLIC mtx_finalize

      PUBLIC mtx_setup
      PUBLIC mtx_set_matrix
      PUBLIC mtx_set_source
      PUBLIC mtx_set_vector
      PUBLIC mtx_solve
      PUBLIC mtx_get_vector_j
      PUBLIC mtx_get_vector
      PUBLIC mtx_gather_vector
      PUBLIC mtx_cleanup

      PUBLIC mtxc_setup
      PUBLIC mtxc_set_matrix
      PUBLIC mtxc_set_source
      PUBLIC mtxc_set_vector
      PUBLIC mtxc_solve
      PUBLIC mtxc_get_vector_j
      PUBLIC mtxc_get_vector
      PUBLIC mtxc_gather_vector
      PUBLIC mtxc_cleanup

      TYPE(mtx_mpi_type):: mtx_global

!
!  Description: Solves a linear system in parallel with KSP (Fortran code).
!
! -----------------------------------------------------------------------


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                    Include files
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  This program uses CPP for preprocessing, as indicated by the use of
!  PETSc include files in the directory petsc/include/finclude.  This
!  convention enables use of the CPP preprocessor, which allows the use
!  of the #include statements that define PETSc objects and variables.
!
!  Use of the conventional Fortran include statements is also supported
!  In this case, the PETsc include files are located in the directory
!  petsc/include/foldinclude.
!         
!  Since one must be very careful to include each file no more than once
!  in a Fortran routine, application programmers must exlicitly list
!  each file needed for the various PETSc components within their
!  program (unlike the C/C++ interface).
!
!  See the Fortran section of the PETSc users manual for details.
!
!  The following include statements are required for KSP Fortran programs:
!     petsc.h       - base PETSc routines
!     petscvec.h    - vectors
!     petscmat.h    - matrices
!     petscpc.h     - preconditioners
!     petscksp.h    - Krylov subspace methods
!  Include the following to use PETSc random numbers:
!     petscsys.h    - system routines
!  Additional include statements may be needed if using additional
!  PETSc routines in a Fortran program, e.g.,
!     petscviewer.h - viewers
!     petscis.h     - index sets
!
!  for PETSc V3.8: remove #include 
!#include "finclude/petsc.h"


!#include "finclude/petscvec.h"
!#include "finclude/petscmat.h"
!#include "finclude/petscpc.h"
!#include "finclude/petscksp.h"
!#include "finclude/petscsys.h"
!#include "finclude/petsckspdef.h"

!  for PETSc V3.8: remove #include 
!#include "finclude/petscvec.h90"
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                   Variable declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     ksp      - linear solver context
!     ksp      - Krylov subspace method context
!     pc       - preconditioner context
!     x, b, u  - approx solution, right-hand-side, exact solution vectors
!     A        - matrix that defines linear system
!     its      - iterations for convergence
!     norm     - norm of error in solution
!     rctx     - random number generator context
!
!  Note that vectors are declared as PETSc "Vec" objects.  These vectors
!  are mathematical objects that contain more than just an array of
!  double precision numbers. I.e., vectors in PETSc are not just
!        double precision x(*).
!  However, local vector data can be easily accessed via VecGetArray().
!  See the Fortran section of the PETSc users manual for details.
!  
      PetscInt       istart,iend,irange,imax
      Vec            x,b
      Mat            A 
      KSP            ksp
      PC             pc
      INTEGER,PARAMETER:: ione=1
      REAL(8),PARAMETER:: one=1.D0
      INTEGER,DIMENSION(:),POINTER:: istartx,iendx,isizex
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      CONTAINS

!-----

      SUBROUTINE mtx_initialize
      IMPLICIT NONE
      INTEGER:: ierr

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_initialize: PetscInitialize: ierr=',ierr
      ncomm=PETSC_COMM_WORLD
      CALL mtx_set_communicator_global(ncomm)
      mtx_global%comm=ncomm
      mtx_global%rank=nrank
      mtx_global%size=nsize
      mtx_global%rankg=0
      mtx_global%sizeg=1
      mtx_global%rankl=nrank
      mtx_global%sizel=nsize
      return
      END SUBROUTINE mtx_initialize

!-----

      SUBROUTINE mtx_finalize
      IMPLICIT NONE
      INTEGER:: ierr

      call PetscFinalize(ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_finalize: PetscFinalize: ierr=',ierr
      END SUBROUTINE mtx_finalize

!-----

      SUBROUTINE mtx_setup(imax_,istart_,iend_,jwidth,nzmax,idebug)

      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(OUT):: istart_,iend_  ! allocated range of lines 
      INTEGER,OPTIONAL,INTENT(IN):: jwidth ! band matrix width, not used in KSP
      INTEGER,OPTIONAL,INTENT(IN):: nzmax  ! number of nonzero components
      INTEGER,OPTIONAL,INTENT(IN):: idebug ! debug level
      INTEGER:: i,ierr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!      Compute the matrix and right-hand-side vector that define
!      the linear system, Ax = b.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!  Create parallel matrix, specifying only its global dimensions.
!  When using MatCreate(), the matrix format can be specified at
!  runtime. Also, the parallel partitioning of the matrix is
!  determined by PETSc at runtime.

      IF(PRESENT(idebug)) THEN
         idebug_save=idebug
      ELSE
         idebug_save=0
      END IF

      imax=imax_
      call MatCreate(ncomm,A,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: MatCreate: ierr=',ierr
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,imax,imax,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: MatSetSizes: ierr=',ierr
      call MatSetFromOptions(A,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: MatSetFromOptions: ierr=',ierr
      call MatSetUp(A,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: MatSetUp: ierr=',ierr

      call MatGetOwnershipRange(A,istart,iend,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: MatGetOwnershipRange: ierr=',ierr
      istart_=istart+1
      iend_=iend
      irange=iend-istart

      call VecCreate(ncomm,b,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: VecCreate: ierr=',ierr
      call VecSetSizes(b,PETSC_DECIDE,imax,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: VecSetSizes: ierr=',ierr
      call VecSetFromOptions(b,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: VecSetFromOptions: ierr=',ierr
      call VecDuplicate(b,x,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: VecDuplicate: ierr=',ierr

!  Currently, all PETSc parallel matrix formats are partitioned by
!  contiguous chunks of rows across the processors.  Determine which
!  rows of the matrix are locally owned. 

      ALLOCATE(istartx(0:nsize-1),iendx(0:nsize-1),isizex(0:nsize-1))

      call mtx_allgather1_integer(istart,istartx)
!      write(6,*) nrank,istart,istartx(nrank)
      call mtx_allgather1_integer(iend,iendx)

      do i=0,nsize-1
         isizex(i)=iendx(i)-istartx(i)
!         write(6,*) istartx(i),iendx(i),isizex(i)
      enddo

!      if(nrank.eq.0) then
!         write(6,'(A)') '# mtx_setup: '
!         do i=0,nsize-1
!            write(6,'(A,4I10)') '#  nrank,istart,iend,isize=', &
!                  i,istartx(i)+1,iendx(i),isizex(i)
!         enddo
!      endif

!  Create linear solver context

      call KSPCreate(ncomm,ksp,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: KSPCreate: ierr=',ierr

      return
      END SUBROUTINE mtx_setup

      SUBROUTINE mtx_set_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      REAL(8),INTENT(IN):: v    ! value to be inserted
      INTEGER:: ierr

      IF(i.GE.istart+1.AND.i.LE.iend) THEN
         call MatSetValues(A,ione,i-1,ione,j-1,v,INSERT_VALUES,ierr)
         IF(ierr.NE.0) WRITE(6,*) &
              'XX mtx_set_matrix: MatSetValues: ierr=',ierr
      ELSE
         write(6,'(A)') &
              'XX libmtxksp:mtx_set_matrix: i : out of range'
         write(6,'(A,3I10)') '   ,istart,iend,i=',istart+1,iend,i
      ENDIF

      return
      END SUBROUTINE mtx_set_matrix
      
      SUBROUTINE mtx_set_source(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted
      INTEGER:: ierr

      IF(j.GE.istart.AND.j.LE.iend) THEN
         call VecSetValues(b,ione,j-1,v,INSERT_VALUES,ierr)
         IF(ierr.NE.0) WRITE(6,*) &
              'XX mtx_set_source: VecSetValues: ierr=',ierr
      ELSE
         write(6,'(A)') &
              'XX libmtxksp:mtx_set_source: j : out of range'
         write(6,'(A,4I10)') '   nrank,istart,iend,j=',nrank,istart+1,iend,j
      ENDIF
      return
      END SUBROUTINE mtx_set_source

      SUBROUTINE mtx_set_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted
      INTEGER:: ierr

      IF(j.GE.istart.AND.j.LE.iend) THEN
         call VecSetValues(x,ione,j-1,v,INSERT_VALUES,ierr)
         IF(ierr.NE.0) WRITE(6,*) &
              'XX mtx_set_vector: VecSetValues: ierr=',ierr
      ELSE
         write(6,'(A)') &
              'XX libmtxksp:mtx_set_vector: j : out of range'
         write(6,'(A,4I10)') '   nrank,istart,iend,j=',nrank,istart+1,iend,j
      ENDIF
      RETURN
      END SUBROUTINE mtx_set_vector
      
      SUBROUTINE mtx_split_operation
      INTEGER:: ierr

      call VecAssemblyBegin(b,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_split_operation: VecAssemblyBegin: b: ierr=',ierr
      call VecAssemblyEnd(b,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_split_operation: VecAssemblyEnd: b: ierr=',ierr
      call VecAssemblyBegin(x,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_split_operation: VecAssemblyBegin: x: ierr=',ierr
      call VecAssemblyEnd(x,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_split_operation: VecAssemblyEnd: x: ierr=',ierr
      return
      END SUBROUTINE mtx_split_operation

      SUBROUTINE mtx_solve(itype,tolerance,its, &
           methodKSP,methodPC,damping_factor,emax,emin,max_steps)

      INTEGER,INTENT(IN):: itype     ! for future use
      REAL(8),INTENT(IN):: tolerance
      INTEGER,INTENT(OUT):: its
      INTEGER,OPTIONAL:: methodKSP,methodPC,max_steps
      REAL(8),OPTIONAL:: damping_factor,emax,emin
      INTEGER:: methodKSP_,methodPC_,max_steps_
      REAL(8):: damping_factor_,emax_,emin_
      INTEGER:: ierr

!      Options:
!        itype: preconditioning matrix pattern, initial guess
!        itype=0: DIFFERENT_NONZERO_PATTERN, initial guess zero
!        itype=1: DIFFERENT_NONZERO_PATTERN, initial guess nonzero
!        itype=2: SAME_NONZERO_PATTERN, initial guess zero
!        itype=3: SAME_NONZERO_PATTERN, initial guess nonzero
!        itype=4: SAME_PRECONDITIONER, initial guess zero
!        itype=5: SAME_PRECONDITIONER, initial guess nonzero
!
!        methodKSP: type of KSP solver (default=4)
!        methodKSP= 0: Richardson (optional damping factor = 1.0)
!        methodKSP= 1: Chebychev (optional emin=0.01 emax=100.)
!        methodKSP= 2: Conjugate Gradient
!        methodKSP= 3: BiConjugate Gradient
!        methodKSP= 4..7: Generalized Gradient Residual (optional max_steps=30)
!           methodKSP= 4:       (Clasical orthogonalization, no refine) 
!           methodKSP= 5:       (Clasical orthogonalization, refine if needed)
!           methodKSP= 6:       (Clasical orthogonalization, refine always)
!           methodKSP= 7:       (Modified orthogonalization) XX not supported
!        methodKSP= 8: BiCGSTAB
!        methodKSP= 9: Conjugate Gradient Squared (CGS)
!        methodKSP=10: Transpose-Free Quasi Minimal Residual (1)
!        methodKSP=11: Transpose-Free Quasi Minimal Residual (2)
!        methodKSP=12: Conjugate Residual Method
!        methodKSP=13: Least Square Method
!        methodKSP=14: Preconditioning only
!
!        methodPC : type of Pre-Conditioner (default=5 or 1)
!        mdthodPC=  0: Jacobi
!        mdthodPC=  1: Block Jacobi
!        mdthodPC=  2: SOR
!        mdthodPC=  3: SOR with Eisenstat trick
!        mdthodPC=  4: Incomplete Cholesky
!        mdthodPC=  5: Incomplete LU
!        mdthodPC=  6: Additive Schwarz
!        mdthodPC=  7: Linear solver
!        mdthodPC=  8: Combination of preconditioners
!        mdthodPC=  9: LU
!        mdthodPC= 10: Cholesky
!        mdthodPC= 11: No preconditioning
!        mdthodPC= 12: Shell for user-defined PC

!  Assemble matrix, using the 2-step process:
!       MatAssemblyBegin(), MatAssemblyEnd()
!  Computations can be done while messages are in transition,
!  by placing code between these two statements.

      IF(PRESENT(methodKSP)) THEN
         methodKSP_=methodKSP
      ELSE
         methodKSP_=4
      ENDIF
      IF(PRESENT(methodPC)) THEN
         methodPC_=methodPC
      ELSE
         IF(nsize.EQ.1) THEN
            methodPC_=5
         ELSE
            methodPC_=1
         ENDIF
      ENDIF
      IF(PRESENT(max_steps)) THEN
         max_steps_=max_steps
      ELSE
         max_steps_=30
      ENDIF
      IF(PRESENT(damping_factor)) THEN
         damping_factor_=damping_factor
      ELSE
         damping_factor_=1.D0
      ENDIF
      IF(PRESENT(emin)) THEN
         emin_=emin
      ELSE
         emin_=0.01D0
      ENDIF
      IF(PRESENT(emax)) THEN
         emax_=emax
      ELSE
         emax_=100.D0
      ENDIF

      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: MatAssemblyBegin: ierr=',ierr
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: MatAssemblyEnd: ierr=',ierr
      call VecAssemblyBegin(b,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: VecAssemblyBegin: b: ierr=',ierr
      call VecAssemblyEnd(b,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: VecAssemblyEnd: b: ierr=',ierr
      call VecAssemblyBegin(x,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: VecAssemblyBegin: x: ierr=',ierr
      call VecAssemblyEnd(x,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: VecAssemblyEnd: x: ierr=',ierr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!         Create the linear solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!  Set operators. Here the matrix that defines the linear system
!  also serves as the preconditioning matrix.
!
!  for PETSc V3.8
!
!      SELECT CASE(itype/2)
!      CASE(0)
!         CALL KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN,ierr)
!      CASE(1)
!         CALL KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN,ierr)
!      CASE(2)
!         CALL KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER,ierr)
!      END SELECT

      CALL KSPSetOperators(ksp,A,A,ierr)

      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: KSPSetOperators: ierr=',ierr

      IF(MOD(itype,2).eq.1) THEN
         CALL KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
         IF(ierr.NE.0) WRITE(6,*) &
              'XX mtx_solve: KSPSetInitialGuessNonzero: ierr=',ierr
      END IF

!  Set linear solver defaults for this problem (optional).
!   - By extracting the KSP and PC contexts from the KSP context,
!     we can then directly directly call any KSP and PC routines
!     to set various options.
!   - The following four statements are optional; all of these
!     parameters could alternatively be specified at runtime via
!     KSPSetFromOptions(). All of these defaults can be
!     overridden at runtime, as indicated below.

!     We comment out this section of code since the Jacobi
!     preconditioner is not a good general default.

!      call KSPGetPC(ksp,pc,ierr)
!      ptype = PCJACOBI
!      call PCSetType(pc,ptype,ierr)
!      tol = 1.e-7

      call KSPSetup(ksp,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: KSPSetup: ierr=',ierr

      call KSPGetPC(ksp,pc,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: KSPGetPC: ierr=',ierr

      SELECT CASE(methodKSP_)
      CASE(0) 
         CALL KSPSetType(ksp,KSPRICHARDSON,ierr)
         CALL KSPRichardsonSetScale(ksp,damping_factor_,ierr)
      CASE(1)
         CALL KSPSetType(ksp,KSPCHEBYCHEV,ierr)
         CALL KSPChebyshevSetEigenvalues(ksp,emax_,emin_,ierr)
      CASE(2)
         CALL KSPSetType(ksp,KSPCG,ierr)
      CASE(3)
         CALL KSPSetType(ksp,KSPBICG,ierr)
      CASE(4:7)
         CALL KSPSetType(ksp,KSPGMRES,ierr)
         CALL KSPGMRESSetRestart(ksp,max_steps_,ierr)
         SELECT CASE(methodKSP)
         CASE(5)
            CALL KSPGMRESSetCGSRefinementType(ksp, &
                 KSP_GMRES_CGS_REFINE_IFNEEDED,ierr)
         CASE(6)
            CALL KSPGMRESSetCGSRefinementType(ksp, &
                 KSP_GMRES_CGS_REFINE_ALWAYS,ierr)
         CASE(7)
            !    CALL KSPGMRESSetOrthogonalization(ksp, &
            !    KSPGMRESModifiedGramSchmidtOrthogonalization,ierr)
         END SELECT
      CASE(8)
         CALL KSPSetType(ksp,KSPBCGS,ierr)
      CASE(9)
         CALL KSPSetType(ksp,KSPCGS,ierr)
      CASE(10)
         CALL KSPSetType(ksp,KSPTFQMR,ierr)
      CASE(11)
         CALL KSPSetType(ksp,KSPTCQMR,ierr)
      CASE(12)
         CALL KSPSetType(ksp,KSPCR,ierr)
      CASE(13)
         CALL KSPSetType(ksp,KSPLSQR,ierr)
      CASE(14)
         CALL KSPSetType(ksp,KSPPREONLY,ierr)
      END SELECT
	
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: KSPSetType: methodKSP,ierr=',methodKSP_,ierr

      SELECT CASE(methodPC_)
      CASE(0) 
         CALL PCSetType(pc,PCJACOBI,ierr)
      CASE(1) 
         CALL PCSetType(pc,PCBJACOBI,ierr)
      CASE(2) 
         CALL PCSetType(pc,PCSOR,ierr)
      CASE(3) 
         CALL PCSetType(pc,PCEISENSTAT,ierr)
      CASE(4) 
         CALL PCSetType(pc,PCICC,ierr)
      CASE(5) 
         CALL PCSetType(pc,PCILU,ierr)
      CASE(6) 
         CALL PCSetType(pc,PCASM,ierr)
      CASE(7) 
         CALL PCSetType(pc,PCKSP,ierr)
      CASE(8) 
         CALL PCSetType(pc,PCCOMPOSITE,ierr)
      CASE(9) 
         CALL PCSetType(pc,PCLU,ierr)
      CASE(10) 
         CALL PCSetType(pc,PCCHOLESKY,ierr)
      CASE(11) 
         CALL PCSetType(pc,PCNONE,ierr)
      CASE(12) 
         CALL PCSetType(pc,PCSHELL,ierr)
      END SELECT

      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: PC: methodPC,ierr=',methodPC_,ierr

!      call KSPSetTolerances(ksp,tolerance, &
!           PETSC_DEFAULT_DOUBLE_PRECISION, &
!           PETSC_DEFAULT_DOUBLE_PRECISION, &
!           PETSC_DEFAULT_INTEGER,ierr)
      call KSPSetTolerances(ksp,tolerance, &
           PETSC_DEFAULT_REAL, &
           PETSC_DEFAULT_REAL, &
           PETSC_DEFAULT_INTEGER,ierr)

      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: KSPSetTolerances: ierr=',ierr

!  Set user-defined monitoring routine if desired

!      call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-my_ksp_monitor',  &
!                               flg,ierr)
!      if (flg) then
!        call KSPMonitorSet(ksp,MyKSPMonitor,PETSC_NULL_OBJECT,          &
!                           PETSC_NULL_FUNCTION,ierr)
!      endif


!  Set runtime options, e.g.,
!      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
!  These options will override those specified above as long as
!  KSPSetFromOptions() is called _after_ any other customization
!  routines.

      call KSPSetFromOptions(ksp,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: KSPSetFromOptions: ierr=',ierr

!  Set convergence test routine if desired

!      call PetscOptionsHasName(PETSC_NULL_CHARACTER,                    &
!           '-my_ksp_convergence',flg,ierr)
!      if (flg) then
!        call KSPSetConvergenceTest(ksp,MyKSPConverged,                  &
!                PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,ierr)
!      endif
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                      Solve the linear system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      call KSPSolve(ksp,b,x,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: KSPSolve: ierr=',ierr
      call KSPGetIterationNumber(ksp,its,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: KSPGetIterationNumber: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_solve

      SUBROUTINE mtx_get_vector_j(j,v)

      INTEGER,INTENT(IN):: j
      REAL(8),INTENT(OUT):: v
      PetscScalar,pointer:: x_value(:)
      INTEGER:: ierr

      call VecGetArrayF90(x,x_value,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_get_vector: VecGetArrayF90: ierr=',ierr
      v=x_value(j-Istart)
      call VecRestoreArrayF90(x,x_value,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_get_vector: VecRestoreArrayF90: ierr=',ierr

      RETURN
      END SUBROUTINE mtx_get_vector_j

      SUBROUTINE mtx_get_vector(v)

      REAL(8),DIMENSION(irange),intent(out):: v
      PetscScalar,pointer:: x_value(:)
      INTEGER:: j,ierr,imax_

      call VecGetArrayF90(x,x_value,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_vector: VecGetArrayF90: ierr=',ierr
      do j=1,irange
         v(j)=x_value(j)
      enddo
      call VecRestoreArrayF90(x,x_value,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_vector: VecRestoreArrayF90: ierr=',ierr

      RETURN
      END SUBROUTINE mtx_get_vector

      SUBROUTINE mtx_gather_vector(x_)

      REAL(8),DIMENSION(imax),INTENT(OUT):: x_
      REAL(8),DIMENSION(irange):: v
      PetscScalar,pointer:: x_value(:)
      INTEGER:: j,ierr,imax_

      call VecGetArrayF90(x,x_value,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_vector: VecGetArrayF90: ierr=',ierr
      do j=1,irange
         v(j)=x_value(j)
      enddo
      call VecRestoreArrayF90(x,x_value,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_vector: VecRestoreArrayF90: ierr=',ierr

      imax_=imax      

      call mtx_allgatherv_real8(v,iend-istart,x_,imax_,isizex,istartx)

      RETURN
      END SUBROUTINE mtx_gather_vector

      SUBROUTINE mtx_cleanup
      INTEGER:: ierr

!  Free work space.  All PETSc objects should be destroyed when they
!  are no longer needed.

      call KSPDestroy(ksp,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_cleanup: KSPDestroy: ierr=',ierr

      DEALLOCATE(istartx,iendx,isizex)

      call VecDestroy(x,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_cleanup: VecDestroy: x: ierr=',ierr
      call VecDestroy(b,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_cleanup: VecDestroy: b: ierr=',ierr
      call MatDestroy(A,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_cleanup: MatDestroy: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_cleanup

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Complex inteface for compatibility
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      SUBROUTINE mtxc_setup(imax_,istart_,iend_,jwidth,nzmax,idebug)

      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(OUT):: istart_,iend_  ! allocated range of lines 
      INTEGER,OPTIONAL,INTENT(IN):: jwidth ! band matrix width
      INTEGER,OPTIONAL,INTENT(IN):: nzmax  ! number of nonzero components
      INTEGER,OPTIONAL,INTENT(IN):: idebug ! debug level

      istart_=0
      iend_=0
      RETURN
      END SUBROUTINE mtxc_setup
      
!-----

      SUBROUTINE mtxc_set_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      COMPLEX(8),INTENT(IN):: v    ! value to be inserted

      RETURN
      END SUBROUTINE mtxc_set_matrix
      
!-----

      SUBROUTINE mtxc_set_source(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      COMPLEX(8),INTENT(IN):: v ! value to be inserted

      RETURN
      END SUBROUTINE mtxc_set_source
      
!-----

      SUBROUTINE mtxc_set_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      RETURN
      END SUBROUTINE mtxc_set_vector
      
!-----

      SUBROUTINE mtxc_solve(itype,tolerance,its, &
           methodKSP,methodPC,damping_factor,emax,emin,max_steps)
      INTEGER,INTENT(IN):: itype     ! info level
      REAL(8),INTENT(IN):: tolerance
      INTEGER,INTENT(OUT):: its
      INTEGER,OPTIONAL:: methodKSP,methodPC,max_steps
      REAL(8),OPTIONAL:: damping_factor,emax,emin

      its=0
      RETURN
      END SUBROUTINE mtxc_solve

!-----

      SUBROUTINE mtxc_get_vector_j(j,v)
        INTEGER,INTENT(IN):: j
        COMPLEX(8),INTENT(OUT):: v

        v=0.D0
        RETURN
      END SUBROUTINE mtxc_get_vector_j

!-----

      SUBROUTINE mtxc_get_vector(v)
        COMPLEX(8),DIMENSION(irange),INTENT(OUT):: v

        v(1:irange)=0.D0
        RETURN
      END SUBROUTINE mtxc_get_vector

!-----

      SUBROUTINE mtxc_gather_vector(v)
        COMPLEX(8),DIMENSION(imax),INTENT(OUT):: v

        v(1:imax)=0.D0
        RETURN
      END SUBROUTINE mtxc_gather_vector

!-----

      SUBROUTINE mtxc_cleanup

        RETURN
      END SUBROUTINE mtxc_cleanup

      END MODULE libmtx
