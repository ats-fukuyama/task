MODULE trgtmp
! **************************************************************************
!           Temporal evolution outputs
! **************************************************************************

  USE trcomm,ONLY: ikind,rkind,nsamax,neqmax,neqrmax,neq_neqr,nsa_neq,ngt
  USE libgrf,ONLY: grd1d
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_gr_temporal

  CHARACTER(LEN=30) :: label

  INTEGER(ikind) :: nr,nsa,neq,neqr

  REAL(rkind),DIMENSION(:),ALLOCATABLE :: gt  !(0:ngt)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: & !(0:ngt,nsamax)
       gt1,gt2,gt3,gt4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: & !(0:ngt,5)
       gti1

CONTAINS

  SUBROUTINE tr_gr_temporal(k2)
! --------------------------------------------------------------------------
!            Control routine of temporal evolution outputs
! --------------------------------------------------------------------------
    USE trcomm, ONLY: gvt

    CHARACTER(LEN=1),INTENT(IN) :: k2
    INTEGER(ikind) :: i2,iosts

    CALL tr_gr_temp_alloc

    ! set axis
    gt(0:ngt) = gvt(0:ngt,0)

    ! control pages
    READ(k2,'(I1)',IOSTAT=iosts) i2
    IF(iosts /= 0)THEN
       WRITE(6,*) ' ERROR : Unsupported graph ID'
       RETURN
    END IF

    SELECT CASE(i2)
    CASE(1)
       CALL tr_gr_temp1 ! n(0),u(0),t(0),q(0),q(a)
    CASE(2)
       CALL tr_gr_temp2 ! I_pl,etc...
    END SELECT
    
    RETURN
  END SUBROUTINE tr_gr_temporal

! *************************************************************************

  SUBROUTINE tr_gr_temp1
  ! ----- time evolution of (n, u, T, q)-----
    USE trcomm,ONLY: gvt,gvts

    gt(0:ngt)=gvt(0:ngt,0)
    DO nsa=1,nsamax
       gt1(0:ngt,nsa)=gvts(0:ngt,nsa,1)
       gt2(0:ngt,nsa)=gvts(0:ngt,nsa,2)
       gt3(0:ngt,nsa)=gvts(0:ngt,nsa,3)
    END DO

    gt4(0:ngt,1)=gvt(0:ngt,1)
    gt4(0:ngt,2)=gvt(0:ngt,2)

    CALL PAGES
    LABEL = '/n(0) [10^20/m^3] vs t/'
    CALL GRD1D(1,gt,gt1,ngt+1,ngt+1,nsamax,label,0)
    LABEL = '/u(0) vs t/'
    CALL GRD1D(2,gt,gt2,ngt+1,ngt+1,nsamax,label,0)
    LABEL = '/T(0) [keV] vs t/'
    CALL GRD1D(3,gt,gt3,ngt+1,ngt+1,nsamax,label,0)
    LABEL = '/q(0),q(a) vs t/'
    CALL GRD1D(4,gt,gt4,ngt+1,ngt+1,2,label,0)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_temp1

! *************************************************************************

  SUBROUTINE tr_gr_temp2
    ! ----- time evolution of (j, q) -----
    USE trcomm, ONLY: gvt,gvti

    gt(0:ngt)=gvt(0:ngt,0)

    gti1(0:ngt,1) = gvti(0:ngt,1)
!    gti1(0:ngt,2) = gvti(0:ngt,2)
!    gti1(0:ngt,3) = gvti(0:ngt,3)

    gt4(0:ngt,1)=gvt(0:ngt,1)
    gt4(0:ngt,2)=gvt(0:ngt,2)

    CALL PAGES
    LABEL = '/Ipl [MA] vs t/'
    CALL GRD1D(1,gt,gti1,ngt+1,ngt+1,5,label,0)
    LABEL = '/q(0),q(a) vs t/'
    CALL GRD1D(4,gt,gt4,ngt+1,ngt+1,2,label,0)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_temp2


! *************************************************************************
! *************************************************************************
! *************************************************************************

  SUBROUTINE tr_gr_temp_alloc

    INTEGER(ikind),SAVE :: ngt_save, nsamax_save
    INTEGER(ikind)      :: ierr

    IF(ngt /= ngt_save .OR. nsamax /= nsamax_save)THEN

       IF(ngt_save /= 0) CALL tr_gr_temp_dealloc

       DO
          ALLOCATE(gt(0:ngt),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gt1(0:ngt,nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gt2(0:ngt,nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gt3(0:ngt,nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gt4(0:ngt,nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gti1(0:ngt,5),STAT=ierr); IF(ierr /= 0) EXIT

          ngt_save    = ngt
          nsamax_save = nsamax
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_temp_alloc: allocation error: ierr=',ierr

    END IF
    RETURN
  END SUBROUTINE tr_gr_temp_alloc

  SUBROUTINE tr_gr_temp_dealloc

    IF(ALLOCATED(gt)) DEALLOCATE(gt)
    IF(ALLOCATED(gt1)) DEALLOCATE(gt1)
    IF(ALLOCATED(gt2)) DEALLOCATE(gt2)
    IF(ALLOCATED(gt3)) DEALLOCATE(gt3)
    IF(ALLOCATED(gt4)) DEALLOCATE(gt4)
    IF(ALLOCATED(gti1)) DEALLOCATE(gti1)

  END SUBROUTINE tr_gr_temp_dealloc

END MODULE trgtmp
