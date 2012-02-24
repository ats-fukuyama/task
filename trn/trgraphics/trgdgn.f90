MODULE trgdgn

  USE trcomm, ONLY: ikind,rkind,nrmax
  USE libgrf, ONLY: grd1d
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_gr_diagnostic

  CHARACTER(LEN=30) :: label

  REAL(rkind),DIMENSION(:),ALLOCATABLE :: rhomg ! (1:nrmax)
  REAL(rkind),DIMENSION(:),ALLOCATABLE :: &     ! (0:nrmax)
       nrd1g,nrd2g,nrd3g,nrd4g
  REAL(rkind),DIMENSION(:),ALLOCATABLE :: &     ! (1:nrmax)
       nrd1mg,nrd2mg,nrd3mg,nrd4mg

CONTAINS

  SUBROUTINE tr_gr_diagnostic(k2)
    USE trcomm, ONLY: rhom,nrd1,nrd2,nrd3,nrd4

    CHARACTER(LEN=1),INTENT(IN) :: k2
    INTEGER(ikind) :: iosts,i2

    CALL tr_gr_diagnostic_alloc

    ! set axis
    rhomg(1:nrmax) = rhom(1:nrmax)

    READ(k2,'(I1)',IOSTAT=iosts) i2

    !--- for diagnostic array
    nrd1mg(1:nrmax) = nrd1(1:nrmax)
    nrd2mg(1:nrmax) = nrd2(1:nrmax)
    nrd3mg(1:nrmax) = nrd3(1:nrmax)
    nrd4mg(1:nrmax) = nrd4(1:nrmax)

    nrd1g(0:nrmax) = nrd1(0:nrmax)
    nrd2g(0:nrmax) = nrd2(0:nrmax)
    nrd3g(0:nrmax) = nrd3(0:nrmax)
    nrd4g(0:nrmax) = nrd4(0:nrmax)

    CALL PAGES
    LABEL = '/diagnostic1 vs rho/'
    CALL GRD1D(1,rhomg,nrd1mg, nrmax, nrmax, 1, LABEL, 0)
    LABEL = '/diagnostic2 vs rho/'
    CALL GRD1D(2,rhomg,nrd2mg, nrmax, nrmax, 1, LABEL, 0)
    LABEL = '/diagnostic3 vs rho/'
    CALL GRD1D(3,rhomg,nrd3mg, NRMAX+1, NRMAX+1, 1, LABEL, 0)
    LABEL = '/diagnostic4 vs rho/'
    CALL GRD1D(4,rhomg,nrd4mg, NRMAX+1, NRMAX+1, 1, LABEL, 0)
    CALL PAGEE    

  END SUBROUTINE tr_gr_diagnostic

  SUBROUTINE tr_gr_diagnostic_alloc

    INTEGER(ikind),SAVE :: nrmax_save
    INTEGER(ikind)      :: ierr

    IF(nrmax /= nrmax_save)THEN

       IF(nrmax_save /= 0) CALL tr_gr_diagnostic_dealloc

       DO
          ALLOCATE(rhomg(1:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT

          ALLOCATE(nrd1g(0:nrmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd2g(0:nrmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd3g(0:nrmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd4g(0:nrmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd1mg(1:nrmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd2mg(1:nrmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd3mg(1:nrmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd4mg(1:nrmax),STAT=ierr); IF(ierr /=0) EXIT

          nrmax_save   = nrmax
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_diagnostic_alloc: allocation error: ierr=',ierr

    END IF
    RETURN
  END SUBROUTINE tr_gr_diagnostic_alloc

  SUBROUTINE tr_gr_diagnostic_dealloc

    IF(ALLOCATED(rhomg)) DEALLOCATE(rhomg)
    IF(ALLOCATED(nrd1g)) DEALLOCATE(nrd1g)
    IF(ALLOCATED(nrd2g)) DEALLOCATE(nrd2g)
    IF(ALLOCATED(nrd3g)) DEALLOCATE(nrd3g)
    IF(ALLOCATED(nrd4g)) DEALLOCATE(nrd4g)
    IF(ALLOCATED(nrd1mg)) DEALLOCATE(nrd1mg)
    IF(ALLOCATED(nrd2mg)) DEALLOCATE(nrd2mg)
    IF(ALLOCATED(nrd3mg)) DEALLOCATE(nrd3mg)
    IF(ALLOCATED(nrd4mg)) DEALLOCATE(nrd4mg)

    RETURN
  END SUBROUTINE tr_gr_diagnostic_dealloc

END MODULE trgdgn
