!     $Id$

      MODULE libmtx

      IMPLICIT NONE
      PUBLIC mtx_initialize
      PUBLIC mtx_setup
      PUBLIC mtx_set_matrix
      PUBLIC mtx_add_matrix
      PUBLIC mtx_set_source
      PUBLIC mtx_add_source
      PUBLIC mtx_set_vector
      PUBLIC mtx_add_vector
      PUBLIC mtx_split_operation
      PUBLIC mtx_solve
      PUBLIC mtx_get_vector
      PUBLIC mtx_gather_vector
      PUBLIC mtx_cleanup
      PUBLIC mtx_finalize
      PUBLIC mtx_barrier
      PUBLIC mtx_broadcast_character
      PUBLIC mtx_broadcast_integer
      PUBLIC mtx_broadcast_real8
      PUBLIC mtx_broadcast_complex8
      PUBLIC mtx_gather_integer
      PUBLIC mtx_allgather_integer
      PUBLIC mtx_gatherv_real8
      PUBLIC mtx_allgatherv_real8
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
#include "finclude/petsc.h"
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
      PetscMPIInt    rank,size
      PetscInt       istart,iend,ilen,imax
      PetscErrorCode ierr
      Vec            x,b,xx
      Mat            A 
      KSP            ksp
      INTEGER,PARAMETER:: ione=1
      REAL(8),PARAMETER:: one=1.D0
      INTEGER,DIMENSION(:),POINTER:: istartx,iendx,isiz

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      CONTAINS

      SUBROUTINE mtx_initialize(rank_,size_)
      INTEGER,INTENT(OUT):: rank_,size_

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
      call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
      rank_=rank
      size_=size
!      if(rank.eq.0) write(6,'(A,I10)') '# mtx_initialize: size=',size
      return

      END SUBROUTINE mtx_initialize

      SUBROUTINE mtx_setup(imax_,istart_,iend_,jwidth_)

      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(IN):: jwidth_         ! band matrix width
      INTEGER,INTENT(OUT):: istart_,iend_   ! allocated range of lines 
      INTEGER:: i

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
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,imax,imax,ierr)
      call MatSetFromOptions(A,ierr)

      call MatGetOwnershipRange(A,istart,iend,ierr)
      istart_=istart+1
      iend_=iend
      ilen=iend-istart

      call VecCreate(PETSC_COMM_WORLD,b,ierr)
      call VecSetSizes(b,PETSC_DECIDE,imax,ierr)
      call VecSetFromOptions(b,ierr)
      call VecDuplicate(b,x,ierr)

!  Currently, all PETSc parallel matrix formats are partitioned by
!  contiguous chunks of rows across the processors.  Determine which
!  rows of the matrix are locally owned. 

      ALLOCATE(istartx(0:size-1),iendx(0:size-1),isiz(0:size-1))

      call mtx_allgather_integer(istart,istartx,size)
      call mtx_allgather_integer(iend,iendx,size)
!      call MPI_ALLGATHER(istart,1,MPI_INTEGER,
!     &                   istartx,1,MPI_INTEGER,
!     &                   PETSC_COMM_WORLD,ierr)
!      call MPI_ALLGATHER(iend,1,MPI_INTEGER,
!     &                   iendx,1,MPI_INTEGER,
!     &                   PETSC_COMM_WORLD,ierr)

      do i=0,size-1
         isiz(i)=iendx(i)-istartx(i)
      enddo
!      if(rank.eq.0) then
!         write(6,'(A)') '# mtx_setup: '
!         do i=0,size-1
!            write(6,'(A,4I10)') '#  rank,istart,iend,isiz=',
!     &           i,istartx(i)+1,iendx(i),isiz(i)
!         enddo
!      endif

!  Create linear solver context

      call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)

      return
      END SUBROUTINE mtx_setup

      SUBROUTINE mtx_set_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      REAL(8),INTENT(IN):: v    ! value to be inserted

      call MatSetValues(A,ione,i-1,ione,j-1,v,INSERT_VALUES,ierr)
      return
      END SUBROUTINE mtx_set_matrix
      
      SUBROUTINE mtx_add_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      REAL(8),INTENT(IN):: v    ! value to be inserted

      call MatSetValues(A,ione,i-1,ione,j-1,v,ADD_VALUES,ierr)
      return
      END SUBROUTINE mtx_add_matrix
      
      SUBROUTINE mtx_set_source(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      call VecSetValues(b,ione,j-1,v,INSERT_VALUES,ierr)
      return
      END SUBROUTINE mtx_set_source
      
      SUBROUTINE mtx_add_source(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      call VecSetValues(b,ione,j-1,v,ADD_VALUES,ierr)
      call VecAssemblyBegin(b,ierr)
      call VecAssemblyEnd(b,ierr)
      return
      END SUBROUTINE mtx_add_source
      
      SUBROUTINE mtx_set_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      call VecSetValues(x,ione,j-1,v,INSERT_VALUES,ierr)
      return
      END SUBROUTINE mtx_set_vector
      
      SUBROUTINE mtx_add_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      call VecAssemblyBegin(x,ierr)
      return
      END SUBROUTINE mtx_add_vector

      SUBROUTINE mtx_split_operation
      call VecAssemblyBegin(b,ierr)
      call VecAssemblyEnd(b,ierr)
      call VecAssemblyBegin(x,ierr)
      call VecAssemblyEnd(x,ierr)
      return
      END SUBROUTINE mtx_split_operation
      SUBROUTINE mtx_solve(itype,tolerance,its)
      INTEGER,INTENT(IN):: itype     ! for futuer use
      REAL(8),INTENT(IN):: tolerance
      INTEGER,INTENT(OUT):: its

!  Assemble matrix, using the 2-step process:
!       MatAssemblyBegin(), MatAssemblyEnd()
!  Computations can be done while messages are in transition,
!  by placing code between these two statements.

      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
      call VecAssemblyBegin(b,ierr)
      call VecAssemblyEnd(b,ierr)
      call VecAssemblyBegin(x,ierr)
      call VecAssemblyEnd(x,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!         Create the linear solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!  Set operators. Here the matrix that defines the linear system
!  also serves as the preconditioning matrix.

      call KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN,ierr)

      if(itype.eq.1) call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)

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

      call KSPSetTolerances(ksp,tolerance,
     &     PETSC_DEFAULT_DOUBLE_PRECISION,
     &     PETSC_DEFAULT_DOUBLE_PRECISION,
     &     PETSC_DEFAULT_INTEGER,ierr)

!  Set user-defined monitoring routine if desired

!      call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-my_ksp_monitor',  &
!     &                    flg,ierr)
!      if (flg) then
!        call KSPMonitorSet(ksp,MyKSPMonitor,PETSC_NULL_OBJECT,          &
!     &                     PETSC_NULL_FUNCTION,ierr)
!      endif


!  Set runtime options, e.g.,
!      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
!  These options will override those specified above as long as
!  KSPSetFromOptions() is called _after_ any other customization
!  routines.

      call KSPSetFromOptions(ksp,ierr)

!  Set convergence test routine if desired

!      call PetscOptionsHasName(PETSC_NULL_CHARACTER,                    &
!     &     '-my_ksp_convergence',flg,ierr)
!      if (flg) then
!        call KSPSetConvergenceTest(ksp,MyKSPConverged,                  &
!     &          PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,ierr)
!      endif
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                      Solve the linear system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      call KSPSolve(ksp,b,x,ierr)
      call KSPGetIterationNumber(ksp,its,ierr)
      RETURN
      END SUBROUTINE mtx_solve

      SUBROUTINE mtx_get_vector(j,v)

      INTEGER,INTENT(IN):: j
      REAL(8),INTENT(OUT):: v
      PetscScalar:: x_value(1)
      PetscOffset:: x_offset

      call VecGetArray(x,x_value,x_offset,ierr)
      v=x_value(x_offset+j-Istart)
      call VecRestoreArray(x,x_value,x_offset,ierr)

      RETURN
      END SUBROUTINE mtx_get_vector

      SUBROUTINE mtx_gather_vector(x_)

      REAL(8),DIMENSION(imax),INTENT(OUT):: x_
      REAL(8),DIMENSION(ilen):: v
      INTEGER:: j
      PetscScalar:: x_value(1)
      PetscOffset:: x_offset

      call VecGetArray(x,x_value,x_offset,ierr)
      do j=1,ilen
         v(j)=x_value(x_offset+j)
      enddo
      call VecRestoreArray(x,x_value,x_offset,ierr)

      call mtx_allgatherv_real8(v,iend-istart,x_,imax,isiz,istartx)
!      call mtx_gatherv_real8(v,iend-istart,x_,imax,isiz,istartx)
!      call MPI_GATHERV(v,isiz(rank),MPI_DOUBLE_PRECISION,
!     &                 x_,isiz,istartx,MPI_DOUBLE_PRECISION,
!     &                 0,PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_gather_vector

      SUBROUTINE mtx_cleanup

!  Free work space.  All PETSc objects should be destroyed when they
!  are no longer needed.

      call KSPDestroy(ksp,ierr)

      DEALLOCATE(istartx,iendx,isiz)

      call VecDestroy(x,ierr)
      call VecDestroy(b,ierr)
      call MatDestroy(A,ierr)
      RETURN
      END SUBROUTINE mtx_cleanup

      SUBROUTINE mtx_finalize

!  Always call PetscFinalize() before exiting a program.  This routine
!    - finalizes the PETSc libraries as well as MPI
!    - provides summary and diagnostic information if certain runtime
!      options are chosen (e.g., -log_summary).  See PetscFinalize()
!      manpage for more information.

      call PetscFinalize(ierr)
      END SUBROUTINE mtx_finalize

      SUBROUTINE mtx_barrier
      call MPI_Barrier(PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_barrier

      SUBROUTINE mtx_broadcast_character(kdata,n)

      CHARACTER(LEN=n),INTENT(INOUT):: kdata
      INTEGER,INTENT(IN):: n
      
      call mpi_bcast(kdata,n,MPI_CHARACTER,
     &               0,PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_broadcast_character

      SUBROUTINE mtx_broadcast_integer(idata,n)

      INTEGER,DIMENSION(n),INTENT(INOUT):: idata
      INTEGER,INTENT(IN):: n
      
      call mpi_bcast(idata,n,MPI_INTEGER,
     &               0,PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_broadcast_integer

      SUBROUTINE mtx_broadcast_real8(vdata,n)

      REAL(8),DIMENSION(n),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n
      
      call mpi_bcast(vdata,n,MPI_DOUBLE_PRECISION,
     &               0,PETSC_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE mtx_broadcast_real8

      SUBROUTINE mtx_broadcast_complex8(vdata,n)

      COMPLEX(8),DIMENSION(n),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n
      
      call mpi_bcast(vdata,n,MPI_DOUBLE_COMPLEX,
     &               0,PETSC_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE mtx_broadcast_complex8

      SUBROUTINE mtx_gather_integer(idata,itot,ntot)

      INTEGER,INTENT(IN):: idata
      INTEGER,INTENT(INOUT):: ntot
      INTEGER,DIMENSION(ntot),INTENT(OUT):: itot

      call MPI_GATHER(idata,1,MPI_INTEGER,
     &                itot,1,MPI_INTEGER,
     &                0,PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_gather_integer

      SUBROUTINE mtx_allgather_integer(idata,itot,ntot)

      INTEGER,INTENT(IN):: idata
      INTEGER,INTENT(INOUT):: ntot
      INTEGER,DIMENSION(ntot),INTENT(OUT):: itot

      call MPI_ALLGATHER(idata,1,MPI_INTEGER,
     &                   itot,1,MPI_INTEGER,
     &                   PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_allgather_integer

      SUBROUTINE mtx_gatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa)

      INTEGER,INTENT(IN):: ndata
      INTEGER,INTENT(INOUT):: ntot
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,DIMENSION(size):: ilena,iposa

      call MPI_GATHERV(vdata,ndata,MPI_DOUBLE_PRECISION,
     &                 vtot,ilena,iposa,MPI_DOUBLE_PRECISION,
     &                 0,PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_gatherv_real8

      SUBROUTINE mtx_allgatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa)

      INTEGER,INTENT(IN):: ndata
      INTEGER,INTENT(INOUT):: ntot
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,DIMENSION(size):: ilena,iposa

      call MPI_ALLGATHERV(vdata,ndata,MPI_DOUBLE_PRECISION,
     &                    vtot,ilena,iposa,MPI_DOUBLE_PRECISION,
     &                    PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_allgatherv_real8

      END MODULE libmtx
