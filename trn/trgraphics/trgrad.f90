MODULE trgrad
! **************************************************************************
!          Snap shot and histroy of radial profile
! **************************************************************************

  USE trcomm, ONLY:ikind,rkind,nrmax,nsamax,neqmax,neqrmax, &
       neq_neqr,nsa_neq,nva_neq,rhog
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
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:nrmax,5)
       vgx1,vgx2,vgx3,vgx4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:nrmax,5)
       vmx1,vmx2,vmx3,vmx4
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
          CALL tr_gr_rad2 ! dtr for particle equation
       CASE(3)
          CALL tr_gr_rad3 ! dtr for energy equation
       CASE(4)
          CALL tr_gr_rad4 ! vtr for particle and energy equation
       CASE(5)
          CALL tr_gr_rad5 ! jtot,jtor,qp
       CASE(6)
          CALL tr_gr_rad6 ! heating profile
       CASE(7)
          CALL tr_gr_rad7 ! rotational velocity, Vtor,Vpol,Vpar,Vprp
       END SELECT
    ELSE IF(i2 == 1)THEN ! history of radial profile
       SELECT CASE(i3)
       CASE(1)
          CALL tr_gr_rad11 ! rn,rt(e,D),qp
       CASE(5)
          CALL tr_gr_rad15 ! jtot,joh,etc... ,qp
       END SELECT
    END IF

    RETURN
  END SUBROUTINE tr_gr_radial

! **************************************************************************
  SUBROUTINE tr_gr_rad1
  ! ----- current radial profile of (n, u, T, q)-----
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
!       vg4(0:nrmax,1)=dpdrho(0:nrmax)
       vg4(0:nrmax,1)=qp(0:nrmax)

    CALL PAGES
    label = '/n [10^20/m^3] vs rho/'
    CALL GRD1D(1,rhog,vg1,nrmax+1,nrmax+1,nsamax,label,0)
    label = '/u vs rho/'
    CALL GRD1D(2,rhog,vg2,nrmax+1,nrmax+1,nsamax,label,0)
    label = '/T [keV] vs rho/'
    CALL GRD1D(3,rhog,vg3,nrmax+1,nrmax+1,nsamax,label,0)
    label = '/q vs rho/'
    CALL GRD1D(4,rhog,vg4,nrmax+1,nrmax+1,     1,label,0)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_rad1

! **************************************************************************
  SUBROUTINE tr_gr_rad2
  ! ----- particle diffusion coefficients of each species -----
    USE trcomm, ONLY: dtr,dtr_nc,dtr_tb,dtr_prv, &
         vtr,vtr_nc,vtr_tb,vtr_prv,mdltr_prv

    REAL(rkind),DIMENSION(nsamax,1:nrmax) :: dtrg_tb,dtrg_nc
    REAL(rkind),DIMENSION(3,1:nrmax) :: dtrg_s1,dtrg_s2
    INTEGER(ikind) :: neq,nsa,nva,nk

    vm1(1:nrmax,1:neqmax) = 0.d0
    vm2(1:nrmax,1:neqmax) = 0.d0

    DO neq = 1, neqmax
       nsa = nsa_neq(neq)
       nva = nva_neq(neq)
       IF(nva == 1)THEN ! particle
       IF(mdltr_prv /= 0)THEN ! for Pereverzev method
          dtrg_nc(nsa,1:nrmax) = dtr_nc(neq,neq,1:nrmax)

          dtrg_tb(nsa,1:nrmax) = dtr_tb(neq,neq,1:nrmax) &
                               - dtr_prv(neq,1:nrmax)

          IF(nsa == 1)THEN 
             dtrg_s1(1,1:nrmax) = dtr(neq,neq,1:nrmax)    &
                                - dtr_prv(neq,1:nrmax)
             dtrg_s1(2,1:nrmax) = dtr_nc(neq,neq,1:nrmax)
             dtrg_s1(3,1:nrmax) = dtr_tb(neq,neq,1:nrmax) &
                                - dtr_prv(neq,1:nrmax)
          ELSE IF(nsa == 2)THEN
             dtrg_s2(1,1:nrmax) = dtr(neq,neq,1:nrmax)    &
                                - dtr_prv(neq,1:nrmax)
             dtrg_s2(2,1:nrmax) = dtr_nc(neq,neq,1:nrmax)
             dtrg_s2(3,1:nrmax) = dtr_tb(neq,neq,1:nrmax) &
                                - dtr_prv(neq,1:nrmax)
          END IF
          vm1(1:nrmax,nsa) = dtrg_nc(nsa,1:nrmax)
          vm2(1:nrmax,nsa) = MIN(dtrg_tb(nsa,1:nrmax),50.d0)
       ELSE
          IF(nsa == 1)THEN 
             dtrg_s1(1,1:nrmax) = dtr(neq,neq,1:nrmax)
             dtrg_s1(2,1:nrmax) = dtr_nc(neq,neq,1:nrmax)
             dtrg_s1(3,1:nrmax) = dtr_tb(neq,neq,1:nrmax)
          ELSE IF(nsa == 2)THEN
             dtrg_s2(1,1:nrmax) = dtr(neq,neq,1:nrmax)
             dtrg_s2(2,1:nrmax) = dtr_nc(neq,neq,1:nrmax)
             dtrg_s2(3,1:nrmax) = dtr_tb(neq,neq,1:nrmax)
          END IF
          vm1(1:nrmax,nsa) = dtr_nc(neq,neq,1:nrmax)
          vm2(1:nrmax,nsa) = MIN(dtr_tb(neq,neq,1:nrmax),50.d0)
       END IF
       END IF
    END DO
    
    DO nk = 1, 3
       vmx1(1:nrmax,nk) = MIN(dtrg_s1(nk,1:nrmax),60.d0)
       vmx2(1:nrmax,nk) = MIN(dtrg_s2(nk,1:nrmax),60.d0)
    END DO
    
    CALL PAGES
    label = '/Dnc_s [m^2/s] vs rho/'
    CALL GRD1D(1,rhomg,vm1, nrmax, nrmax, nsamax,label, 0)
    label = '/Dtb_s [m^2/s] vs rho/'
    CALL GRD1D(2,rhomg,vm2, nrmax, nrmax, nsamax,label, 0)
    label = '/D(1) tot,nc,tb/'
    CALL GRD1D(3,rhomg,vmx1,nrmax, nrmax, 3,label,0)
    label = '/D(2) tot,nc,tb/'
    CALL GRD1D(4,rhomg,vmx2,nrmax, nrmax, 3,label,0)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_rad2

! **************************************************************************
  SUBROUTINE tr_gr_rad3
  ! ----- heat diffusion coefficients of each species -----
    USE trcomm, ONLY: dtr,dtr_nc,dtr_tb,dtr_prv,mdltr_prv

    REAL(rkind),DIMENSION(nsamax,1:nrmax) :: dtrg_tb,dtrg_nc
    REAL(rkind),DIMENSION(3,1:nrmax) :: dtrg_s1,dtrg_s2
    INTEGER(ikind) :: neq,nsa,nva,nk

    vm1(1:nrmax,1:neqmax) = 0.d0
    vm2(1:nrmax,1:neqmax) = 0.d0

    DO neq = 1, neqmax
       nsa = nsa_neq(neq)
       nva = nva_neq(neq)
       IF(nva == 3)THEN ! heat
       IF(mdltr_prv /= 0)THEN ! for Pereverzev method
          dtrg_nc(nsa,1:nrmax) = dtr_nc(neq,neq,1:nrmax)

          dtrg_tb(nsa,1:nrmax) = dtr_tb(neq,neq,1:nrmax) &
                               - dtr_prv(neq,1:nrmax)

          IF(nsa == 1)THEN 
             dtrg_s1(1,1:nrmax) = dtr(neq,neq,1:nrmax)    &
                                - dtr_prv(neq,1:nrmax)
             dtrg_s1(2,1:nrmax) = dtr_nc(neq,neq,1:nrmax)
             dtrg_s1(3,1:nrmax) = dtr_tb(neq,neq,1:nrmax) &
                                - dtr_prv(neq,1:nrmax)
          ELSE IF(nsa == 2)THEN
             dtrg_s2(1,1:nrmax) = dtr(neq,neq,1:nrmax)    &
                                - dtr_prv(neq,1:nrmax)
             dtrg_s2(2,1:nrmax) = dtr_nc(neq,neq,1:nrmax)
             dtrg_s2(3,1:nrmax) = dtr_tb(neq,neq,1:nrmax) &
                                - dtr_prv(neq,1:nrmax)
          END IF
          vm1(1:nrmax,nsa) = dtrg_nc(nsa,1:nrmax)
          vm2(1:nrmax,nsa) = MIN(dtrg_tb(nsa,1:nrmax),50.d0)
       ELSE
          IF(nsa == 1)THEN 
             dtrg_s1(1,1:nrmax) = dtr(neq,neq,1:nrmax)
             dtrg_s1(2,1:nrmax) = dtr_nc(neq,neq,1:nrmax)
             dtrg_s1(3,1:nrmax) = dtr_tb(neq,neq,1:nrmax)
          ELSE IF(nsa == 2)THEN
             dtrg_s2(1,1:nrmax) = dtr(neq,neq,1:nrmax)
             dtrg_s2(2,1:nrmax) = dtr_nc(neq,neq,1:nrmax)
             dtrg_s2(3,1:nrmax) = dtr_tb(neq,neq,1:nrmax)
          END IF
          vm1(1:nrmax,nsa) = dtr_nc(neq,neq,1:nrmax)
          vm2(1:nrmax,nsa) = MIN(dtr_tb(neq,neq,1:nrmax),50.d0)
       END IF
       END IF
    END DO
    
    DO nk = 1, 3
       vmx1(1:nrmax,nk) = MIN(dtrg_s1(nk,1:nrmax),60.d0)
       vmx2(1:nrmax,nk) = MIN(dtrg_s2(nk,1:nrmax),60.d0)
    END DO
    
    CALL PAGES
    label = '/chi_nc_s [m^2/s] vs rho/'
    CALL GRD1D(1,rhomg,vm1, nrmax, nrmax, nsamax,label, 0)
    label = '/chi_tb_s [m^2/s] vs rho/'
    CALL GRD1D(2,rhomg,vm2, nrmax, nrmax, nsamax,label, 0)
    label = '/chi(1)tot,nc,tb [m^2/s] vs rho/'
    CALL GRD1D(3,rhomg,vmx1,nrmax, nrmax, 3,label,0)
    label = '/chi(2)tot,nc,tb [m^2/s] vs rho/'
    CALL GRD1D(4,rhomg,vmx2,nrmax, nrmax, 3,label,0)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_rad3

! **************************************************************************
  SUBROUTINE tr_gr_rad4
  ! ----- pinch velocity profile for particle and energy equation -----
    USE trcomm, ONLY: vtr,vtr_nc,vtr_tb,vtr_prv,mdltr_prv,eta

    REAL(rkind),DIMENSION(1:nsamax,1:nrmax) :: vtrnc_d,vtrnc_chi
    INTEGER(ikind) :: neq,nsa,nva

    DO neq = 1, neqmax
       nsa = nsa_neq(neq)
       nva = nva_neq(neq)
       IF(mdltr_prv /= 0)THEN ! for Pereverzev method
       IF(nva == 1)THEN ! particle
          vmx1(1:nrmax,nsa) = vtr(neq,neq,1:nrmax)-vtr_prv(neq,1:nrmax)
       ELSE IF(nva == 3)THEN ! energy
          vmx2(1:nrmax,nsa) = vtr(neq,neq,1:nrmax)-vtr_prv(neq,1:nrmax)
       END IF

       ELSE IF(mdltr_prv == 0)THEN
       IF(nva == 1)THEN ! particle
          vmx1(1:nrmax,nsa) = vtr(neq,neq,1:nrmax)
       ELSE IF(nva == 3)THEN ! energy
          vmx2(1:nrmax,nsa) = vtr(neq,neq,1:nrmax)
       END IF
       END IF
    END DO
    vgx1(0:nrmax,1) = LOG10(eta(0:nrmax))

    CALL PAGES
    label = '/V_D [m/s] vs rho/'
    CALL GRD1D(1,rhomg,vmx1,nrmax,nrmax,nsamax,label,0)
    label = '/V_chi [m/s] vs rho/'
    CALL GRD1D(2,rhomg,vmx2,nrmax,nrmax,nsamax,label,0)
    label = '/eta_par [ohm m] vs rho/'
    CALL GRD1D(3,rhog,vgx1,nrmax+1,nrmax+1,1,label,2)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_rad4

! **************************************************************************
  SUBROUTINE tr_gr_rad5
  ! ----- current density profile -----
    USE trcomm, ONLY: jtot,joh,jtor,eta,qp,htr

    vgx1(0:nrmax,1) = 1.d-6*jtot(0:nrmax) + 1.d-6*htr(1,0:nrmax)
    vgx1(0:nrmax,2) = 1.d-6*joh(0:nrmax)
    vgx1(0:nrmax,3) = 1.d-6*htr(1,0:nrmax)

    vgx2(0:nrmax,1) = 1.d-6*jtor(0:nrmax)
    vgx3(0:nrmax,1) = 1.d-6*htr(1,0:nrmax)
    vgx4(0:nrmax,1) = LOG10(eta(0:nrmax))
    
!    vgx4(0:nrmax,1) = qp(0:nrmax)

    CALL PAGES
    label = '/j(tot,oh,bs) [MA/m^2] vs rho'
    CALL GRD1D(1,rhog,vgx1,nrmax+1,nrmax+1,5,label,0)
    label = '/j(tor) [MA/m^2] vs rho'
    CALL GRD1D(2,rhog,vgx2,nrmax+1,nrmax+1,5,label,0)
    label = '/j(bs,ex) [MA/m^2] vs rho'
    CALL GRD1D(3,rhog,vgx3,nrmax+1,nrmax+1,1,label,0)
    label = '/eta(para) [ohm m] vs rho'
    CALL GRD1D(4,rhog,vgx4,nrmax+1,nrmax+1,1,label,2)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_rad5

! **************************************************************************
  SUBROUTINE tr_gr_rad6
  ! ----- heating profile-----
    USE trcomm, ONLY: str,poh

    vgx1 = 0.d0

!    write(*,*) '*** Unavailable now ***'
    vgx1(0:nrmax,1) = poh(0:nrmax)*1.d-6

    CALL PAGES
    label = '/Pin [MW/m^3] vs rho/'
    CALL GRD1D(1,rhog,vgx1,nrmax+1,nrmax+1,5,label,0)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_rad6

! **************************************************************************
  SUBROUTINE tr_gr_rad7
  ! rotation velocity profile
    USE trcomm, ONLY: vtor,vpol,vpar,vprp

    vgx1(0:nrmax,1) = vtor(0:nrmax)
    vgx2(0:nrmax,1) = vpol(0:nrmax)
    vgx3(0:nrmax,1) = vpar(0:nrmax)
    vgx4(0:nrmax,1) = vprp(0:nrmax)

    CALL PAGES
    label = '/Vtor[m/s] vs rho'
    CALL GRD1D(1,rhog,vgx1,nrmax+1,nrmax+1,1,label,0)
    label = '/Vpol[m/s] vs rho'
    CALL GRD1D(2,rhog,vgx2,nrmax+1,nrmax+1,1,label,0)
    label = '/Vpar[m/s] vs rho'
    CALL GRD1D(3,rhog,vgx3,nrmax+1,nrmax+1,1,label,0)
    label = '/Vprp[m/s] vs rho'
    CALL GRD1D(4,rhog,vgx4,nrmax+1,nrmax+1,1,label,0)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_rad7

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
    label = '/T1(t) [keV] vs rho/'
    CALL GRD1D(1,rhog,gg1,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '/T2(t) [keV] vs rho/'
    CALL GRD1D(2,rhog,gg2,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '/n1(t) [10^20/m^3] vs rho/'
    CALL GRD1D(3,rhog,gg3,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '/qp(t) vs rho/'
    CALL GRD1D(4,rhog,gg4,nrmax+1,nrmax+1,nggmax+1,label,0)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_rad11

! **************************************************************************
  SUBROUTINE tr_gr_rad15
  ! ----- history of radial profile -----
    USE trcomm, ONLY: ngt,gvrt,gvrtj

    ngg_interval = ngt/(MOD(ngt-1,nggmax)+1)
    DO ngg = 0, nggmax
       gg1(0:nrmax,ngg) = 1.d-6*gvrtj(0:nrmax, ngg*ngg_interval, 1)
       gg2(0:nrmax,ngg) = 1.d-6*gvrtj(0:nrmax, ngg*ngg_interval, 2)
!       gg3(0:nrmax,ngg) = gvrtj(0:nrmax, ngg*ngg_interval, 3)
       ! history of qp profile
       gg4(0:nrmax,ngg) =  gvrt(0:nrmax, ngg*ngg_interval, 1)
    END DO

    CALL PAGES
    label = '/j_tot(t) [MA/m^2] vs rho/'
    CALL GRD1D(1,rhog,gg1,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '/j_oh(t) [MA/m^2] vs rho/'
    CALL GRD1D(2,rhog,gg2,nrmax+1,nrmax+1,nggmax+1,label,0)
!    label = '/j_ex(t) vs rho/'
!    CALL GRD1D(3,rhog,gg3,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '/qp(t) vs rho/'
    CALL GRD1D(4,rhog,gg4,nrmax+1,nrmax+1,nggmax+1,label,0)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_rad15

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
          ALLOCATE(vgx1(0:nrmax,5),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vgx2(0:nrmax,5),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vgx3(0:nrmax,5),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vgx4(0:nrmax,5),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vmx1(1:nrmax,5),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vmx2(1:nrmax,5),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vmx3(1:nrmax,5),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vmx4(1:nrmax,5),STAT=ierr); IF(ierr /=0) EXIT
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
    IF(ALLOCATED(vgx1)) DEALLOCATE(vgx1)
    IF(ALLOCATED(vgx2)) DEALLOCATE(vgx2)
    IF(ALLOCATED(vgx3)) DEALLOCATE(vgx3)
    IF(ALLOCATED(vmx4)) DEALLOCATE(vmx4)
    IF(ALLOCATED(vmx1)) DEALLOCATE(vmx1)
    IF(ALLOCATED(vmx2)) DEALLOCATE(vmx2)
    IF(ALLOCATED(vmx3)) DEALLOCATE(vmx3)
    IF(ALLOCATED(vmx4)) DEALLOCATE(vmx4)
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
