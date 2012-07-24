MODULE trgcom

  USE trcomm, ONLY: ikind,rkind,nrmax,nsamax,neqrmax,lmaxtr,nitmax, &
       rhog
  USE libgrf, ONLY: grd1d
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_gr_comp

  CHARACTER(LEN=30) :: label
  INTEGER(ikind),PARAMETER :: nggmax=10
  INTEGER(ikind) :: nr,nit,nsa,ngg,ngg_interval
  
  REAL(rkind),DIMENSION(:),ALLOCATABLE :: rhomg      !(1:nrmax)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &          !(0:nrmax,1:nsamax)
       vga1,vga2
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &          !(0:nrmax,1:nsamax)
       vma1,vma2
  REAL(rkind),DIMENSION(:),ALLOCATABLE :: ig, err_ig !(1:lmaxtr)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: lt       !(0:nrmax,nsamax)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &
       vgap1,vgap2,vgap3,vgap4                   !(0:nrmax,0:nggmax)

CONTAINS

  SUBROUTINE tr_gr_comp(k2,k3)
! -------------------------------------------------------------------------
!        Control routine of computational parameters outputs
! -------------------------------------------------------------------------
    USE trcomm, ONLY: rhom

    CHARACTER(LEN=1),INTENT(IN) :: k2,k3
    INTEGER(ikind) :: i2,i3,iosts

    CALL tr_gr_comp_alloc

    ! set axis
    rhomg(1:nrmax) = rhom(1:nrmax)
    DO nit = 1, lmaxtr
       ig(nit) = dble(nit)
    END Do

    ! control pages
    READ(k2,'(I1)',IOSTAT=iosts) i2
    READ(k3,'(I1)',IOSTAT=iosts) i3
    IF(iosts /= 0)THEN
       WRITE(6,*) ' ERROR : Unsupported graph ID'
       RETURN
    END IF

    IF(k3 .EQ. ' ')THEN
       SELECT CASE(i2)
       CASE(1)
          CALL tr_gr_comp1 ! convergence of non-linear itearation
       CASE(2)
          CALL tr_gr_comp2 ! Outputs for Pereverzev method: dtr_prv, vtr_prv
       END SELECT
    ELSE IF(i2 == 1)THEN ! history of radial profile
       SELECT CASE(i3)
       CASE(2)
          CALL tr_gr_comp12 ! history of prof. of add. quantities (Pereverzev)
       END SELECT
    END IF

    RETURN
  END SUBROUTINE tr_gr_comp

! *************************************************************************

  SUBROUTINE tr_gr_comp1
  !--- convergence of non-linear iteration---
    USE trcomm, ONLY: error_it, rt_ecl, rt_icl

    lt = 0.d0

    err_ig(1:lmaxtr) = log10(error_it(1:lmaxtr))
    lt(0:nrmax,1) = - rt_ecl(0:nrmax)
    lt(0:nrmax,2) = - rt_icl(0:nrmax)

    ! GRD1D: MODE = 2 ; X:LINEAR  Y:LOG
    CALL PAGES
    label = '/convergence vs NIT/'
    CALL GRD1D(1,ig,  err_ig,lmaxtr, lmaxtr, 1,      label,2)
    label = '/Temp. scale length vs rho/'
    CALL GRD1D(2,rhog,lt,    nrmax+1,nrmax+1,nsamax, label,0)
    CALL PAGEE    

  END SUBROUTINE tr_gr_comp1

! *************************************************************************

  SUBROUTINE tr_gr_comp2
  !--- Outputs for Pereverzev method (Numerical stabilazation method) ---
    USE trcomm, ONLY: mdltr_prv,rt,dtr_prv,vtr_prv

    IF(mdltr_prv == 0)THEN
       vma1 = 0.d0
       vma2 = 0.d0
       DO nsa = 1, nsamax
          vga1(0:nrmax,nsa) = rt(nsa,0:nrmax)
       END DO
    ELSE
       DO nsa = 1, nsamax
          vga1(0:nrmax,nsa) = rt(nsa,0:nrmax)
          vma1(1:nrmax,nsa) = dtr_prv(1+3*nsa,1:nrmax)
          vma2(1:nrmax,nsa) = vtr_prv(1+3*nsa,1:nrmax)
       END DO
    END IF

    CALL PAGES
    label = '/T vs rho/'
    CALL GRD1D(1,rhog, vga1, nrmax+1,nrmax+1,nsamax, label, 0)
    label = '/add_Diff(chi) vs rho/'
    CALL GRD1D(2,rhomg,vma1, nrmax,  nrmax,  nsamax, label, 0)
    label = '/add_Conv(vel) vs rho/'
    CALL GRD1D(3,rhomg,vma2, nrmax,  nrmax,  nsamax, label, 0)
    CALL PAGEE
    
  END SUBROUTINE tr_gr_comp2

! *************************************************************************

  SUBROUTINE tr_gr_comp12
  !--- history of profiles of numerically additional quantities
  !                                               ( for Pereverzev method )
  !    * numerically net additional quantities in each nodal equation *
  !    * These quantities are relative value to total dtr.            *
    USE trcomm, ONLY: ngt,gvrts

    IF(ngt > 0)THEN
       ngg_interval = ngt/(MOD(ngt-1,nggmax)+1)
    ELSE IF(ngt <= 0 )THEN
       ngg_interval = 1
    END IF
    DO ngg = 0, nggmax
       ! thermal diffusivity
       vgap1(0:nrmax,ngg) = gvrts(0:nrmax,ngg*ngg_interval,1,6)
       vgap2(0:nrmax,ngg) = gvrts(0:nrmax,ngg*ngg_interval,2,6)
!!$       ! particle diffusivity
!!$       gparg3(0:nrmax,ngg) = gparts(0:nrmax,ngg*ngg_interval,1,1)
!!$       gparg4(0:nrmax,ngg) = gparts(0:nrmax,ngg*ngg_interval,2,1)
    END DO

    CALL PAGES
    label = '/add_Net(nT(1)) vs rho'
    CALL GRD1D(1,rhog,vgap1, nrmax+1,nrmax+1,nggmax+1, label,0)
    label = '/add_Net(nT(2)) vs rho'
    CALL GRD1D(2,rhog,vgap2, nrmax+1,nrmax+1,nggmax+1, label,0)
    CALL PAGEE
    
  END SUBROUTINE tr_gr_comp12

! *************************************************************************
! *************************************************************************
! *************************************************************************

  SUBROUTINE tr_gr_comp_alloc

    INTEGER(ikind),SAVE :: nrmax_save, lmaxtr_save
    INTEGER(ikind)      :: ierr

    IF(nrmax /= nrmax_save .OR. lmaxtr /= lmaxtr_save)THEN

       IF(nrmax_save /= 0) CALL tr_gr_comp_dealloc

       DO
          ALLOCATE(rhomg(1:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vga1(0:nrmax,nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vga2(0:nrmax,nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vma1(1:nrmax,neqrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vma2(1:nrmax,neqrmax),STAT=ierr); IF(ierr /= 0) EXIT

          ALLOCATE(ig(1:lmaxtr),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(err_ig(1:lmaxtr),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(lt(0:nrmax,nsamax),STAT=ierr); IF(ierr /= 0) EXIT

          ALLOCATE(vgap1(0:nrmax,0:nggmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgap2(0:nrmax,0:nggmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgap3(0:nrmax,0:nggmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgap4(0:nrmax,0:nggmax),STAT=ierr); IF(ierr /= 0) EXIT

          nrmax_save  = nrmax
          lmaxtr_save = lmaxtr
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_comp_alloc: allocation error: ierr=',ierr

    END IF
    RETURN
  END SUBROUTINE tr_gr_comp_alloc

  SUBROUTINE tr_gr_comp_dealloc

    IF(ALLOCATED(rhomg)) DEALLOCATE(rhomg)
    IF(ALLOCATED(vga1)) DEALLOCATE(vga1)
    IF(ALLOCATED(vga2)) DEALLOCATE(vga2)
    IF(ALLOCATED(vma1)) DEALLOCATE(vma1)
    IF(ALLOCATED(vma2)) DEALLOCATE(vma2)

    IF(ALLOCATED(ig)) DEALLOCATE(ig)
    IF(ALLOCATED(err_ig)) DEALLOCATE(err_ig)
    IF(ALLOCATED(lt)) DEALLOCATE(lt)

    IF(ALLOCATED(vgap1)) DEALLOCATE(vgap1)
    IF(ALLOCATED(vgap2)) DEALLOCATE(vgap2)
    IF(ALLOCATED(vgap3)) DEALLOCATE(vgap3)
    IF(ALLOCATED(vgap4)) DEALLOCATE(vgap4)

    RETURN
  END SUBROUTINE tr_gr_comp_dealloc

END MODULE trgcom
