MODULE trgrad
! **************************************************************************
!          Snap shot and histroy of radial profile
! **************************************************************************

  USE trcomm, ONLY:ikind,rkind,nrmax,nsamax,neqmax,neqrmax,neq_neqr,nsa_neq,&
       rhog
  USE libgrf,ONLY: grd1d
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_gr_radial

  CHARACTER(LEN=30) :: label
  INTEGER(ikind),PARAMETER :: nggmax=10
  INTEGER(ikind) :: nr,nsa,neq,neqr,ngg,ngg_interval

  REAL(rkind),DIMENSION(:),ALLOCATABLE :: rhomg !(1:nrmax)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:nrmax,neqrmax)
       vg1,vg2,vg3,vg4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:nrmax,neqrmax)
       vm1,vm2,vm3,vm4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:nrmax,0:nggmax)
       gg1,gg2,gg3,gg4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:nrmax,0:nggmax)
       gm1,gm2,gm3,gm4
  
CONTAINS

  SUBROUTINE tr_gr_radial(k2,k3)
! -------------------------------------------------------------------------
!          Control routine of radial profile outputs
! -------------------------------------------------------------------------
    USE trcomm, ONLY: rhom
    CHARACTER(LEN=1),INTENT(IN) :: k2,k3
    INTEGER(ikind) :: i2,i3,iosts

    CALL tr_gr_rad_alloc

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
          CALL tr_gr_rad1 ! rn,ru,rt,qp
       CASE(2)
          CALL tr_gr_rad2 ! dtr,vtr,str,htr
       END SELECT
    ELSE IF(i2 == 1)THEN ! history of radial profile
       SELECT CASE(i3)
       CASE(1)
          CALL tr_gr_rad11 ! rn,rt(e,D),qp
       END SELECT
    END IF

    RETURN
  END SUBROUTINE tr_gr_radial

! **************************************************************************

  SUBROUTINE tr_gr_rad1
  ! ----- current radial profile -----
    USE trcomm, ONLY: rn,ru,rt,dpdrho,qp

    vg1(0:nrmax,1:neqrmax) = 0.d0
    vg2(0:nrmax,1:neqrmax) = 0.d0
    vg3(0:nrmax,1:neqrmax) = 0.d0
    vg4(0:nrmax,1:neqrmax) = 0.d0

    DO neq = 1, neqmax
       nsa = nsa_neq(neq)
       IF(nsa /= 0)THEN
          vg1(0:nrmax,nsa)=rn(nsa,0:nrmax)
          vg2(0:nrmax,nsa)=ru(nsa,0:nrmax)
          vg3(0:nrmax,nsa)=rt(nsa,0:nrmax)
       END IF
    END DO
       vg4(0:nrmax,1)=dpdrho(0:nrmax)
       vg4(0:nrmax,1)=qp(0:nrmax)

    CALL PAGES
    label = '/n vs rho/'
    CALL GRD1D(1,rhog,vg1,nrmax+1,nrmax+1,nsamax,label,0)
    label = '/u vs rho/'
    CALL GRD1D(2,rhog,vg2,nrmax+1,nrmax+1,nsamax,label,0)
    label = '/T vs rho/'
    CALL GRD1D(3,rhog,vg3,nrmax+1,nrmax+1,nsamax,label,0)
    label = '/d psi d rho vs rho/'
    CALL GRD1D(4,rhog,vg4,nrmax+1,nrmax+1,     1,label,0)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_rad1

! **************************************************************************

  SUBROUTINE tr_gr_rad2
  ! ----- diffusion coefficients -----
    USE trcomm, ONLY: rt,dtr,vtr,str,htr,dtr_prv,vtr_prv,mdltr_prv

    REAL(rkind),DIMENSION(nsamax,1:nrmax) :: dtrg,vtrg

    vg1(0:nrmax,1:neqmax) = 0.d0
    vg2(0:nrmax,1:neqmax) = 0.d0
    vm1(1:nrmax,1:neqmax) = 0.d0
    vm2(1:nrmax,1:neqmax) = 0.d0

    DO neqr = 1, neqrmax
       neq = neq_neqr(neqr)
       nsa = nsa_neq(neq)
       IF(nsa /= 0)THEN
       IF(mdltr_prv /= 0)THEN
          ! for Pereverzev method
          dtrg(nsa,1:nrmax) = dtr(neq,neq,1:nrmax) - dtr_prv(neq-1,1:nrmax)
          vtrg(nsa,1:nrmax) = vtr(neq,neq,1:nrmax) - vtr_prv(neq-1,1:nrmax)

          vg1(0:nrmax,nsa)=rt(nsa,0:nrmax)
          vm1(1:nrmax,nsa)=MIN(dtrg(nsa,1:nrmax),20.D0)
          vg2(0:nrmax,nsa)=str(neq,0:nrmax)
          vm2(1:nrmax,nsa)=vtrg(nsa,1:nrmax)
       ELSE
          vg1(0:nrmax,nsa)=rt(nsa,0:nrmax)
          vm1(1:nrmax,nsa)=MIN(dtr(neq,neq,1:nrmax),20.D0)
          vg2(0:nrmax,nsa)=str(neq,0:nrmax)
          vm2(1:nrmax,nsa)=vtr(neq,neq,1:nrmax)
       END IF
       END IF
    END DO
    
    CALL PAGES
    label = '/T vs rho/'
    CALL GRD1D(1,rhog, vg1, nrmax+1,nrmax+1,nsamax,label, 0)
    label = '/Diffusion vs rho/'
    CALL GRD1D(2,rhomg,vm1, nrmax,  nrmax,  neqrmax,label, 0)
    label = '/Heat_pw vs rho/'
    CALL GRD1D(3,rhog, vg2, nrmax+1,nrmax+1,nsamax,label, 0)
    label = '/Convection vs rho/'
    CALL GRD1D(4,rhomg,vm2, nrmax,  nrmax,  nsamax,label, 0)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_rad2

! **************************************************************************

  SUBROUTINE tr_gr_rad11
  ! ----- history of radial profile -----
    USE trcomm, ONLY: ngt,gvrt,gvrts

    ngg_interval = ngt/(MOD(ngt-1,nggmax)+1)
    DO ngg = 0, nggmax
       gg1(0:nrmax,ngg) = gvrts(0:nrmax, ngg*ngg_interval, 1,3)
       gg2(0:nrmax,ngg) = gvrts(0:nrmax, ngg*ngg_interval, 2,3)
       gg3(0:nrmax,ngg) = gvrts(0:nrmax, ngg*ngg_interval, 1,1)
       gg4(0:nrmax,ngg) =  gvrt(0:nrmax, ngg*ngg_interval, 1)
    END DO

    CALL PAGES
    label = '/T1(t) vs rho/'
    CALL GRD1D(1,rhog,gg1,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '/T2(t) vs rho/'
    CALL GRD1D(2,rhog,gg2,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '/n1(t) vs rho/'
    CALL GRD1D(3,rhog,gg3,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '/qp(t) vs rho/'
    CALL GRD1D(4,rhog,gg4,nrmax+1,nrmax+1,nggmax+1,label,0)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_rad11

! **************************************************************************
! **************************************************************************
! **************************************************************************

  SUBROUTINE tr_gr_rad_alloc
    
    INTEGER(ikind),SAVE :: nrmax_save, neqmax_save
    INTEGER(ikind)      :: ierr

    IF(nrmax /= nrmax_save .OR. neqmax /= neqmax_save)THEN

       IF(nrmax_save /= 0) CALL tr_gr_rad_dealloc

       DO
          ALLOCATE(rhomg(1:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT

          ALLOCATE(vg1(0:nrmax,neqmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vg2(0:nrmax,neqmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vg3(0:nrmax,neqmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vg4(0:nrmax,neqmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vm1(1:nrmax,neqmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vm2(1:nrmax,neqmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vm3(1:nrmax,neqmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vm4(1:nrmax,neqmax),STAT=ierr); IF(ierr /=0) EXIT

          ALLOCATE(gg1(0:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gg2(0:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gg3(0:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gg4(0:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gm1(1:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gm2(1:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gm3(1:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gm4(1:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT

          nrmax_save  = nrmax
          neqmax_save = neqmax
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_rad_alloc: allocation error: ierr=',ierr

    END IF
    RETURN
  END SUBROUTINE tr_gr_rad_alloc

  SUBROUTINE tr_gr_rad_dealloc
    
    IF(ALLOCATED(rhomg)) DEALLOCATE(rhomg)

    IF(ALLOCATED(vg1)) DEALLOCATE(vg1)
    IF(ALLOCATED(vg2)) DEALLOCATE(vg2)
    IF(ALLOCATED(vg3)) DEALLOCATE(vg3)
    IF(ALLOCATED(vg4)) DEALLOCATE(vg4)
    IF(ALLOCATED(vm1)) DEALLOCATE(vm1)
    IF(ALLOCATED(vm2)) DEALLOCATE(vm2)
    IF(ALLOCATED(vm3)) DEALLOCATE(vm3)
    IF(ALLOCATED(vm4)) DEALLOCATE(vm4)

    IF(ALLOCATED(gg1)) DEALLOCATE(gg1)
    IF(ALLOCATED(gg2)) DEALLOCATE(gg2)
    IF(ALLOCATED(gg3)) DEALLOCATE(gg3)
    IF(ALLOCATED(gg4)) DEALLOCATE(gg4)
    IF(ALLOCATED(gm1)) DEALLOCATE(gm1)
    IF(ALLOCATED(gm2)) DEALLOCATE(gm2)
    IF(ALLOCATED(gm3)) DEALLOCATE(gm3)
    IF(ALLOCATED(gm4)) DEALLOCATE(gm4)

    RETURN
  END SUBROUTINE tr_gr_rad_dealloc

END MODULE trgrad
