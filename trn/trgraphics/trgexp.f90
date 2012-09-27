MODULE trgexp
! *************************************************************************
!       graphic output for experimental data
! *************************************************************************

  USE trcomm,ONLY: ikind,rkind,ntum,nrum,nsum,nrmax,rhog, &
                    mdluf,ntxmax,tlmax,time_slc
  USE trgsub,ONLY: tr_gr_time
  USE libgrf,ONLY: grd1d

  IMPLICIT NONE
  PUBLIC tr_gr_exp

  CHARACTER(LEN=50) :: label
  INTEGER(ikind) :: nr

  REAL(rkind),DIMENSION(:),ALLOCATABLE :: rhomg !(1:nrmax)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:nrmax,1:nsum)
       vgu1,vgu2,vgu3,vgu4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:nrmax,1:nsum)
       vmu1,vmu2,vmu3,vmu4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:nrmax,5)
       vgxu1,vgxu2,vgxu3,vgxu4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:nrmax,5)
       vmxu1,vmxu2,vmxu3,vmxu4

CONTAINS

  SUBROUTINE tr_gr_exp(k2,k3)
! -----------------------------------------------------------------------
!         Control routine of experimental data outputs
! -----------------------------------------------------------------------
    USE trcomm, ONLY: rhom
    CHARACTER(LEN=1),INTENT(IN) :: k2,k3
    INTEGER(ikind) :: i2,i3,iosts

    CALL tr_gr_exp_alloc

    ! set axis
    rhomg(1:nrmax) = rhom(1:nrmax)

    READ(k2,'(I1)',IOSTAT=iosts) i2
    READ(k3,'(I1)',IOSTAT=iosts) i3
    IF(iosts /= 0)THEN
       WRITE(6,*) ' ERROR : Unsupported graph ID'
       RETURN
    END IF

    IF(k3 .EQ. ' ')THEN
       SELECT CASE(i2)
       CASE(1)
          CALL tr_gr_exp1
       CASE(2)
          CALL tr_gr_exp2
       END SELECT
    ELSE IF(i2 == 1)THEN
       SELECT CASE(i3)
       CASE(1)
       CASE(2)
       END SELECT
    END IF   

    RETURN
  END SUBROUTINE tr_gr_exp

! ************************************************************************
  SUBROUTINE tr_gr_exp1

    USE trcomm, ONLY:rnu,rtu
    INTEGER(ikind) :: nsu

    vgu1(0:nrmax,1:nsum) = 0.d0
    vgu2(0:nrmax,1:nsum) = 0.d0
    vgu3(0:nrmax,1:nsum) = 0.d0
    vgu4(0:nrmax,1:nsum) = 0.d0

    DO nsu = 1, nsum
       vgu1(0:nrmax,nsu) = rnu(nsu,ntxmax,1:nrmax+1)
       vgu2(0:nrmax,nsu) = rtu(nsu,ntxmax,1:nrmax+1)
    END DO


    CALL PAGES
    label = '/n(exp) [10$+20$=/m$+3$=] vs rho/'
    CALL GRD1D(1,rhog,vgu1,nrmax+1,nrmax+1,nsum,label,0)
    label = '/T(exp) [keV] vs rho/'
    CALL GRD1D(2,rhog,vgu2,nrmax+1,nrmax+1,nsum,label,0)
    CALL tr_gr_time
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_exp1

! ************************************************************************
  SUBROUTINE tr_gr_exp2

!    USE trcomm, ONLY:

    RETURN
  END SUBROUTINE tr_gr_exp2

! ************************************************************************
! ************************************************************************
! ************************************************************************

  SUBROUTINE tr_gr_exp_alloc

    INTEGER(ikind),SAVE :: nrmax_save=0
    INTEGER(ikind)      :: ierr

    IF(nrmax /= nrmax_save)THEN
       
       IF(nrmax_save /= 0 ) CALL tr_gr_exp_dealloc

       DO
          ALLOCATE(rhomg(1:nrmax),STAT=ierr); IF(ierr /= 0) EXIT

          ALLOCATE(vgu1(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgu2(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgu3(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgu4(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmu1(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmu2(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmu3(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmu4(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgxu1(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgxu2(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgxu3(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgxu4(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmxu1(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmxu2(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmxu3(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmxu4(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT

          nrmax_save = nrmax
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_exp_alloc: allocation error: ierr= ', ierr

    END IF
    RETURN
  END SUBROUTINE tr_gr_exp_alloc

  SUBROUTINE tr_gr_exp_dealloc

    IF(ALLOCATED(rhomg)) DEALLOCATE(rhomg)
    IF(ALLOCATED(vgu1)) DEALLOCATE(vgu1)
    IF(ALLOCATED(vgu2)) DEALLOCATE(vgu2)
    IF(ALLOCATED(vgu3)) DEALLOCATE(vgu3)
    IF(ALLOCATED(vgu4)) DEALLOCATE(vgu4)
    IF(ALLOCATED(vmu1)) DEALLOCATE(vmu1)
    IF(ALLOCATED(vmu2)) DEALLOCATE(vmu2)
    IF(ALLOCATED(vmu3)) DEALLOCATE(vmu3)
    IF(ALLOCATED(vmu4)) DEALLOCATE(vmu4)
    IF(ALLOCATED(vgxu1)) DEALLOCATE(vgxu1)
    IF(ALLOCATED(vgxu2)) DEALLOCATE(vgxu2)
    IF(ALLOCATED(vgxu3)) DEALLOCATE(vgxu3)
    IF(ALLOCATED(vgxu4)) DEALLOCATE(vgxu4)
    IF(ALLOCATED(vmxu1)) DEALLOCATE(vmxu1)
    IF(ALLOCATED(vmxu2)) DEALLOCATE(vmxu2)
    IF(ALLOCATED(vmxu3)) DEALLOCATE(vmxu3)
    IF(ALLOCATED(vmxu4)) DEALLOCATE(vmxu4)

    RETURN
  END SUBROUTINE tr_gr_exp_dealloc


END MODULE trgexp
