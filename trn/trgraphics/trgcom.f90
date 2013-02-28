MODULE trgcom

  USE trcomm, ONLY: ikind,rkind,nrmax,nsamax,neqmax,lmaxtr,nitmax, &
       ngt,rhog
  USE trgsub,ONLY: tr_gr1d_rad, tr_gr_time,  &
       tr_gr_vnr_alloc, tr_gr_vnrt_alloc,    &
       tr_gr_lmt_alloc,tr_gr_vnt_alloc,                             &
       tr_gr_init_vg,tr_gr_init_vm,tr_gr_init_lg,                   &
       tr_gr_init_gt, tr_gr_init_gti,  tr_gr_init_gg,               &
       vg1,vg2,vg3,vg4, vm1,vm2,vm3,vm4,   gg1,gg2,gg3,gg4, lg,lg1, &
       gt, gt1,gt2,gt3,gt4, gti1,gti2,gti3,gti4, nggmax,rhomg
  USE libgrf, ONLY: grd1d
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_gr_computation


  CHARACTER(LEN=40) :: label
  INTEGER(ikind)    :: nr,nit,nsa,nva,neq,ngg,ngg_interval,idexp

CONTAINS

  SUBROUTINE tr_gr_computation(k2,k3)
! -------------------------------------------------------------------------
!        Control routine of computational parameters outputs
! -------------------------------------------------------------------------
    USE trcomm, ONLY: rhom, gvt

    CHARACTER(LEN=1),INTENT(IN) :: k2,k3
    INTEGER(ikind) :: i2,i3,ierr,iosts

    CALL tr_gr_vnr_alloc(ierr)
    IF(ierr /= 0) RETURN
    CALL tr_gr_vnrt_alloc(ierr)
    IF(ierr /= 0) RETURN
    CALL tr_gr_lmt_alloc(ierr)
    IF(ierr /= 0) RETURN
    CALL tr_gr_vnt_alloc(ierr)
    IF(ierr /= 0) RETURN

    ! set axis
    rhomg(1:nrmax) = rhom(1:nrmax)
    DO nit = 1, lmaxtr
       lg(nit) = dble(nit)
    END DO

    gt(0:ngt) = gvt(0:ngt,0)

    ! control pages
    READ(k2,'(I1)',IOSTAT=iosts) i2
    READ(k3,'(I1)',IOSTAT=iosts) i3
    IF(iosts /= 0)THEN
       WRITE(6,*) ' ERROR : Unsupported graph ID'
       RETURN
    END IF

    idexp = 0 ! print simulation time on every GSAF page

    IF(k3 .EQ. ' ')THEN
       SELECT CASE(i2)
       CASE(1)
          CALL tr_gr_comp1 ! convergence of non-linear itearation
       CASE(2)
          CALL tr_gr_comp2 ! Outputs for Pereverzev method: dtr_prv, vtr_prv
       CASE(3)
          CALL tr_gr_comp3
       END SELECT
    ELSE IF(i2 == 1)THEN ! history of radial profile
       SELECT CASE(i3)
       CASE(2)
          CALL tr_gr_comp12 ! history of prof. of add. quantities (Pereverzev)
       END SELECT
    ELSE IF(i2 == 2)THEN
       SELECT CASE(i3)
       CASE(1)
       CASE(2)
          CALL tr_gr_comp22
       END SELECT
    END IF

    RETURN
  END SUBROUTINE tr_gr_computation

! *************************************************************************

  SUBROUTINE tr_gr_comp1
  !--- convergence of non-linear iteration---
    USE trcomm, ONLY: error_it, rt_ecl, rt_icl

    CALL tr_gr_init_vg
    CALL tr_gr_init_lg

    FORALL(nit=1:lmaxtr,error_it(nit) > 0)
       lg1(nit) = log10(error_it(nit))
    END FORALL
    vg1(0:nrmax,1) = - rt_ecl(0:nrmax)
    vg1(0:nrmax,2) = - rt_icl(0:nrmax)

    ! tr_gr1d_rad: MODE = 2 ; X:LINEAR  Y:LOG
    CALL PAGES
    label = '@Convergence vs NIT@'
    CALL tr_gr1d_rad(1,lg,lg1,lmaxtr,1,label,2)
    label = '@Temp. scale length vs rho@'
    CALL tr_gr1d_rad(2,rhog,vg1,nrmax+1,2,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE    

  END SUBROUTINE tr_gr_comp1

! *************************************************************************

  SUBROUTINE tr_gr_comp2
  !--- Outputs for Pereverzev method (Numerical stabilazation method) ---
    USE trcomm, ONLY: mdltr_prv,rt,dtr_prv,vtr_prv

    CALL tr_gr_init_vm

    IF(mdltr_prv == 0)THEN
       vm1 = 0.d0
       vm2 = 0.d0
       DO nsa = 1, nsamax
          vg1(0:nrmax,nsa) = rt(nsa,0:nrmax)
       END DO
    ELSE
       DO nsa = 1, nsamax
          vg1(0:nrmax,nsa) = rt(nsa,0:nrmax)
          vm1(1:nrmax,nsa) = dtr_prv(1+3*nsa,1:nrmax)
          vm2(1:nrmax,nsa) = vtr_prv(1+3*nsa,1:nrmax)
       END DO
    END IF

    CALL PAGES
    label = '@T vs rho@'
    CALL tr_gr1d_rad(1,rhog, vg1,nrmax+1,nsamax, label, 0,FMIN0=0.d0)
    label = '@add_Diff(chi) vs rho@'
    CALL tr_gr1d_rad(2,rhomg,vm1,nrmax,nsamax, label, 0,FMIN0=0.d0)
    label = '@add_Conv(vel) vs rho@'
    CALL tr_gr1d_rad(3,rhomg,vm2,nrmax,nsamax,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE
    
  END SUBROUTINE tr_gr_comp2


  SUBROUTINE tr_gr_comp3
    USE trcomm, ONLY: rkev,fluxtb,grdpf,dtr_nl,dtr_prv,nsa_neq,nva_neq

    CALL tr_gr_init_vm

    DO neq = 1, neqmax
       nva = nva_neq(neq)
       nsa = nsa_neq(neq)
       IF(nva == 3)THEN
          vm1(1:nrmax,nsa) = fluxtb(neq,1:nrmax) *1.d-6
          vm2(1:nrmax,nsa) = grdpf(neq,1:nrmax)  *1.d-6
          vm3(1:nrmax,nsa) = dtr_nl(neq,1:nrmax)
          vm4(1:nrmax,nsa) = dtr_prv(neq,1:nrmax)
       END IF
    END DO

    CALL PAGES
    label = '@energy flux [MJ/(m$+2$=s)] vs rho@'
    CALL tr_gr1d_rad(1,rhomg,vm1,nrmax,nsamax,label,0)
    label = '@grad p [MPa/m] vs rho@'
    CALL tr_gr1d_rad(2,rhomg,vm2,nrmax,nsamax,label,0,FMIN0=0.d0)
    label = '@D_nl [m$+2$=/s] vs rho@'
    CALL tr_gr1d_rad(3,rhomg,vm3,nrmax,nsamax,label,0)
    label = '@D_prv [m$+2$=/s] vs rho@'
    CALL tr_gr1d_rad(4,rhomg,vm4,nrmax,nsamax,label,0,FMIN0=0.d0)
    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_comp3

! *************************************************************************

  SUBROUTINE tr_gr_comp12
! -------------------------------------------------------------------------
!   history of profiles of numerically additional quantities
!                                               ( for Pereverzev method )
!    * numerically net additional quantities in each nodal equation *
!    * These quantities are relative value to total dtr.            *
! -------------------------------------------------------------------------
    USE trcomm, ONLY: ngt,gvrts

    CALL tr_gr_init_gg

    ngg_interval = ngt / nggmax
    IF(MOD(ngt,nggmax) /= 0) ngg_interval = ngg_interval + 1

    DO ngg = 0, nggmax-1
       ! thermal diffusivity
       gg1(0:nrmax,ngg) = gvrts(0:nrmax,ngg*ngg_interval,1,7)
       gg2(0:nrmax,ngg) = gvrts(0:nrmax,ngg*ngg_interval,2,7)
       ! particle diffusivity
       gg3(0:nrmax,ngg) = gvrts(0:nrmax,ngg*ngg_interval,1,5)
       gg4(0:nrmax,ngg) = gvrts(0:nrmax,ngg*ngg_interval,2,5)
    END DO
       ! thermal diffusivity
       gg1(0:nrmax,nggmax) = gvrts(0:nrmax,ngt,1,7)
       gg2(0:nrmax,nggmax) = gvrts(0:nrmax,ngt,2,7)
       ! particle diffusivity
       gg3(0:nrmax,nggmax) = gvrts(0:nrmax,ngt,1,5)
       gg4(0:nrmax,nggmax) = gvrts(0:nrmax,ngt,2,5)

    CALL PAGES
    label = '@add_Net(nT(1)) vs rho@'
    CALL tr_gr1d_rad(1,rhog,gg1,nrmax+1,nggmax+1, label,0)
    label = '@add_Net(nT(2)) vs rho@'
    CALL tr_gr1d_rad(2,rhog,gg2,nrmax+1,nggmax+1, label,0)
    label = '@add_Net(n(1)) vs rho@'
    CALL tr_gr1d_rad(3,rhog,gg3, nrmax+1,nggmax+1, label,0)
    label = '@add_Net(n(2)) vs rho@'
    CALL tr_gr1d_rad(4,rhog,gg4, nrmax+1,nggmax+1, label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE
    
  END SUBROUTINE tr_gr_comp12

! ************************************************************************
  
  SUBROUTINE tr_gr_comp22
    ! numerically additional term in nodal equation (relative value)
    USE trcomm, ONLY: gvts

    CALL tr_gr_init_gt

    DO nsa = 1, nsamax
       gt1(0:ngt,nsa) = gvts(0:ngt,nsa,10) ! add_prv for density
       gt2(0:ngt,nsa) = gvts(0:ngt,nsa,11) !         for velocity
       gt3(0:ngt,nsa) = gvts(0:ngt,nsa,12) !         for energy
    END DO

    CALL PAGES
    label = '@MAX_ABS(prv_add_n)  vs t@'
    CALL tr_gr1d_rad(1,gt,gt1,ngt+1,nsamax,label,0)
    label = '@MAX_ABS(prv_add_u) vs t@'
    CALL tr_gr1d_rad(2,gt,gt2,ngt+1,nsamax,label,0)
    label = '@MAX_ABS(prv_add_nT) vs t@'
    CALL tr_gr1d_rad(3,gt,gt3,ngt+1,nsamax,label,0)
    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_comp22

END MODULE trgcom
