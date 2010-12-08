!     $Id$

      MODULE libmtxc

      use libmpi

      IMPLICIT NONE

      PUBLIC mtx_initialize
      PUBLIC mtx_finalize
      PUBLIC mtx_barrier
      PUBLIC mtx_broadcast_character
      PUBLIC mtx_broadcast_integer
      PUBLIC mtx_broadcast_real8
      PUBLIC mtx_broadcast_complex8
      PUBLIC mtx_gather_integer
      PUBLIC mtx_gather_real8
      PUBLIC mtx_allgather_integer
      PUBLIC mtx_gatherv_real8
      PUBLIC mtx_allgatherv_real8
      PUBLIC mtx_allgatherv_complex8

      PUBLIC mtx_setup
      PUBLIC mtx_set_matrix
      PUBLIC mtx_set_source
      PUBLIC mtx_set_vector
      PUBLIC mtx_solve
      PUBLIC mtx_get_vector
      PUBLIC mtx_gather_vector
      PUBLIC mtx_cleanup

      PRIVATE
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
!#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscpc.h"
#include "finclude/petscksp.h"
#include "finclude/petscsys.h"
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
      PetscInt       istart,iend,ilen,imax
      Vec            x,b
      Mat            A 
      KSP            ksp
      INTEGER,PARAMETER:: ione=1
      REAL(8),PARAMETER:: one=1.D0
      INTEGER,DIMENSION(:),POINTER:: istartx,iendx,isiz

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      CONTAINS

      SUBROUTINE mtx_setup(imax_,istart_,iend_,jwidth_)

      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(IN):: jwidth_         ! band matrix width
      INTEGER,INTENT(OUT):: istart_,iend_   ! allocated range of lines 
      INTEGER:: i,ierr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!      Compute the matrix and right-hand-side vector that define
!      the linear system, Ax = b.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!  Create parallel matrix, specifying only its global dimensions.
!  When using MatCreate(), the matrix format can be specified at
!  runtime. Also, the parallel partitioning of the matrix is
!  determined by PETSc at runtime.

      imax=imax_
      call MatCreate(PETSC_COMM_WORLD,A,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: MatCreate: ierr=',ierr
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,imax,imax,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: MatSetSizes: ierr=',ierr
      call MatSetFromOptions(A,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: MatSetFromOptions: ierr=',ierr

      call MatGetOwnershipRange(A,istart,iend,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: MatGetOwnershipRange: ierr=',ierr
      istart_=istart+1
      iend_=iend
      ilen=iend-istart

      call VecCreate(PETSC_COMM_WORLD,b,ierr)
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

      ALLOCATE(istartx(0:nsize-1),iendx(0:nsize-1),isiz(0:nsize-1))

      call mtx_allgather_integer(istart,istartx)
      call mtx_allgather_integer(iend,iendx)

      do i=0,nsize-1
         isiz(i)=iendx(i)-istartx(i)
      enddo

!      if(rank.eq.0) then
!         write(6,'(A)') '# mtx_setup: '
!         do i=0,nsize-1
!            write(6,'(A,4I10)') '#  rank,istart,iend,isiz=', &
!                  i,istartx(i)+1,iendx(i),isiz(i)
!         enddo
!      endif

!  Create linear solver context

      call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_setup: KSPCreate: ierr=',ierr

      return
      END SUBROUTINE mtx_setup

      SUBROUTINE mtx_set_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      COMPLEX(8),INTENT(IN):: v    ! value to be inserted
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
      COMPLEX(8),INTENT(IN):: v ! value to be inserted
      INTEGER:: ierr

      IF(j.GE.istart.AND.j.LE.iend) THEN
         call VecSetValues(b,ione,j-1,v,INSERT_VALUES,ierr)
         IF(ierr.NE.0) WRITE(6,*) &
              'XX mtx_set_source: VecSetValues: ierr=',ierr
      ELSE
         write(6,'(A)') &
              'XX libmtxksp:mtx_set_source: j : out of range'
         write(6,'(A,4I10)') '   rank,istart,iend,j=',rank,istart+1,iend,j
      ENDIF
      return
      END SUBROUTINE mtx_set_source
      
      SUBROUTINE mtx_set_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      COMPLEX(8),INTENT(IN):: v ! value to be inserted
      INTEGER:: ierr

      IF(j.GE.istart.AND.j.LE.iend) THEN
         call VecSetValues(x,ione,j-1,v,INSERT_VALUES,ierr)
         IF(ierr.NE.0) WRITE(6,*) &
              'XX mtx_set_vector: VecSetValues: ierr=',ierr
      ELSE
         write(6,'(A)') &
              'XX libmtxksp:mtx_set_vector: j : out of range'
         write(6,'(A,4I10)') '   rank,istart,iend,j=',rank,istart+1,iend,j
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

      SUBROUTINE mtx_solve(itype,tolerance,its)
      INTEGER,INTENT(IN):: itype     ! for futuer use
      REAL(8),INTENT(IN):: tolerance
      INTEGER,INTENT(OUT):: its
      INTEGER:: ierr

!  Assemble matrix, using the 2-step process:
!       MatAssemblyBegin(), MatAssemblyEnd()
!  Computations can be done while messages are in transition,
!  by placing code between these two statements.

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

      CALL KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_solve: KSPSetOperators: ierr=',ierr

      IF(itype.eq.1) THEN
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

      call KSPSetTolerances(ksp,tolerance, &
           PETSC_DEFAULT_DOUBLE_PRECISION, &
           PETSC_DEFAULT_DOUBLE_PRECISION, &
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

      SUBROUTINE mtx_get_vector(j,v)

      INTEGER,INTENT(IN):: j
      COMPLEX(8),INTENT(OUT):: v
      PetscScalar:: x_value(1)
      PetscOffset:: x_offset
      INTEGER:: ierr

      call VecGetArray(x,x_value,x_offset,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_get_vector: VecGetArray: ierr=',ierr
      v=x_value(x_offset+j-Istart)
      call VecRestoreArray(x,x_value,x_offset,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_get_vector: VecRestoreArray: ierr=',ierr

      RETURN
      END SUBROUTINE mtx_get_vector

      SUBROUTINE mtx_gather_vector(x_)

      COMPLEX(8),DIMENSION(imax),INTENT(OUT):: x_
      COMPLEX(8),DIMENSION(ilen):: v
      INTEGER:: j,ierr
      PetscScalar:: x_value(1)
      PetscOffset:: x_offset

      call VecGetArray(x,x_value,x_offset,ierr)

      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_vector: VecGetArray: ierr=',ierr
      do j=1,ilen
         v(j)=x_value(x_offset+j)
      enddo
      call VecRestoreArray(x,x_value,x_offset,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_vector: VecRestoreArray: ierr=',ierr

      call mtx_allgatherv_complex8(v,iend-istart,x_,imax,isiz,istartx)
      RETURN
      END SUBROUTINE mtx_gather_vector

      SUBROUTINE mtx_cleanup
      INTEGER:: ierr

!  Free work space.  All PETSc objects should be destroyed when they
!  are no longer needed.

      call KSPDestroy(ksp,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_cleanup: KSPDestroy: ierr=',ierr

      DEALLOCATE(istartx,iendx,isiz)

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

      END MODULE libmtxc
