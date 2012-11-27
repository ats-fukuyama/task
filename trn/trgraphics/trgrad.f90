MODULE trgrad
! **************************************************************************
!          Snap shot and histroy of radial profile
! **************************************************************************
  USE trcomm, ONLY:ikind,rkind,nrmax,nsamax,neqmax,neqrmax,        &
       nsabmax,nsafmax,neq_neqr,nsa_neq,nva_neq,nsab_nsaf,rhog
  USE trgsub, ONLY: tr_gr_time, tr_gr_vnr_alloc,tr_gr_vnrt_alloc,    &
       tr_gr_init_vg, tr_gr_init_vm, tr_gr_init_vgx, tr_gr_init_vmx, &
       tr_gr_init_gg,tr_gr_init_gm,                                  &
       vg1,vg2,vg3,vg4, vm1,vm2,vm3,vm4, vgx1,vgx2,vgx3,vgx4,        &
       vmx1,vmx2,vmx3,vmx4, gg1,gg2,gg3,gg4, gm1,gm2,gm3,gm4,        &
       nggmax,rhomg
  USE libgrf,ONLY: grd1d
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_gr_radial

  CHARACTER(LEN=64) :: label
  INTEGER(ikind)    :: nr,nsa,nsaf,neq,neqr,ngg,ngg_interval, idexp
  
CONTAINS

  SUBROUTINE tr_gr_radial(k2,k3)
! -------------------------------------------------------------------------
!          Control routine of radial profile outputs
! -------------------------------------------------------------------------
    USE trcomm, ONLY: rhom
    CHARACTER(LEN=1),INTENT(IN) :: k2,k3
    INTEGER(ikind) :: i2,i3,ierr,iosts

    ierr = 0
    CALL tr_gr_vnr_alloc(ierr)
    IF(ierr /= 0) RETURN
    CALL tr_gr_vnrt_alloc(ierr)
    IF(ierr /= 0) RETURN

    ! set axis
    rhomg(1:nrmax) = rhom(1:nrmax)

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
          CALL tr_gr_rad1 ! rn,(ru,) rp,rt,qp
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
       CASE(8)
          CALL tr_gr_rad8 ! er
       CASE(9)
          CALL tr_gr_rad9 ! fast ion particle
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
    USE trcomm, ONLY: rn,ru,rt,rp,dpdrho,qp

    CALL tr_gr_init_vg

    DO neq = 1, neqmax
       nsa = nsa_neq(neq)
       IF(nsa /= 0)THEN
          IF(nsab_nsaf(nsa) == 0)THEN! excluding fast ion species
             vg1(0:nrmax,nsa)=rn(nsa,0:nrmax)
             vg2(0:nrmax,nsa)=rp(nsa,0:nrmax) * 1.d-6
             vg3(0:nrmax,nsa)=rt(nsa,0:nrmax)
          END IF
       END IF
    END DO
!    write(6,*) rt(2,0:nrmax)
!       vg4(0:nrmax,1)=dpdrho(0:nrmax)
       vg4(0:nrmax,1)=qp(0:nrmax)

    CALL PAGES
    label = '@n [10$+20$=/m$+3$=] vs rho@'
    CALL GRD1D(1,rhog,vg1,nrmax+1,nrmax+1,nsamax,label,0)
    label = '@p [MPa] vs rho@'
    CALL GRD1D(2,rhog,vg2,nrmax+1,nrmax+1,nsamax,label,0)
    label = '@T [keV] vs rho@'
    CALL GRD1D(3,rhog,vg3,nrmax+1,nrmax+1,nsamax,label,0)
    label = '@q vs rho@'
    CALL GRD1D(4,rhog,vg4,nrmax+1,nrmax+1,     1,label,0)

    CALL tr_gr_time(idexp)
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

    CALL tr_gr_init_vm
    CALL tr_gr_init_vmx

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
    label = '@Dnc_s [m$+2$=/s] vs rho@'
    CALL GRD1D(1,rhomg,vm1, nrmax, nrmax, nsamax,label, 0)
    label = '@Dtb_s [m$+2$=/s] vs rho@'
    CALL GRD1D(2,rhomg,vm2, nrmax, nrmax, nsamax,label, 0)
    label = '@D(1) tot,nc,tb [m$+2$=/s] vs rho@'
    CALL GRD1D(3,rhomg,vmx1,nrmax, nrmax, 3,label,0)
    label = '@D(2) tot,nc,tb [m$+2$=/s] vs rho@'
    CALL GRD1D(4,rhomg,vmx2,nrmax, nrmax, 3,label,0)

    CALL tr_gr_time(idexp)
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

    CALL tr_gr_init_vm
    CALL tr_gr_init_vmx

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
    label = '@chi_nc_s [m$+2$=/s] vs rho@'
    CALL GRD1D(1,rhomg,vm1, nrmax, nrmax, nsamax,label, 0)
    label = '@chi_tb_s [m$+2$=/s] vs rho@'
    CALL GRD1D(2,rhomg,vm2, nrmax, nrmax, nsamax,label, 0)
    label = '@chi(1)tot,nc,tb [m$+2$=/s] vs rho@'
    CALL GRD1D(3,rhomg,vmx1,nrmax, nrmax, 3,label,0)
    label = '@chi(2)tot,nc,tb [m$+2$=/s] vs rho@'
    CALL GRD1D(4,rhomg,vmx2,nrmax, nrmax, 3,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_rad3

! **************************************************************************
  SUBROUTINE tr_gr_rad4
  ! ----- pinch velocity profile for particle and energy equation -----
    USE trcomm, ONLY: vtr,vtr_nc,vtr_tb,vtr_prv,mdltr_prv,eta

    REAL(rkind),DIMENSION(1:nsamax,1:nrmax) :: vtrnc_d,vtrnc_chi
    INTEGER(ikind) :: neq,nsa,nva,gsid

    CALL tr_gr_init_vmx
    CALL tr_gr_init_vgx

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
    gsid = 2
    DO nr = 0, nrmax
       IF(eta(nr) <= 0)THEN
          gsid = 0
       END IF
    END DO
    IF(gsid == 2)THEN
       vgx1(0:nrmax,1) = LOG10(eta(0:nrmax))
    ELSE
       WRITE(6,*) 'XX non-positive eta values. Graph scale is set to linear.'
    END IF

    CALL PAGES
    label = '@V_D [m/s] vs rho@'
    CALL GRD1D(1,rhomg,vmx1,nrmax,nrmax,nsamax,label,0)
    label = '@V_chi [m/s] vs rho@'
    CALL GRD1D(2,rhomg,vmx2,nrmax,nrmax,nsamax,label,0)
    label = '@eta_par [ohm m] vs rho@'
    CALL GRD1D(3,rhog,vgx1,nrmax+1,nrmax+1,1,label,gsid)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_rad4

! **************************************************************************
  SUBROUTINE tr_gr_rad5
  ! ----- current density profile -----
    USE trcomm, ONLY: jtot,joh,jcd_nb,jtor,jbs_nc,eta,qp,dpdrho
    INTEGER(ikind) :: gsid

    CALL tr_gr_init_vgx

    vgx1(0:nrmax,1) = 1.d-6*jtot(0:nrmax)
    vgx1(0:nrmax,2) = 1.d-6*joh(0:nrmax)
    vgx1(0:nrmax,3) = 1.d-6*jcd_nb(0:nrmax)
    vgx1(0:nrmax,4) = 1.d-6*jbs_nc(0:nrmax)

    vgx2(0:nrmax,1) = dpdrho(0:nrmax)
    vgx3(0:nrmax,1) = qp(0:nrmax)

    gsid = 2
    DO nr = 0, nrmax
       IF(eta(nr) <= 0)THEN
          gsid = 0
       END IF
    END DO
    IF(gsid == 2)THEN
       vgx4(0:nrmax,1) = LOG10(eta(0:nrmax))
    ELSE
       WRITE(6,*) 'XX non-positive eta values. Graph scale is set to linear.'
    END IF
    
    CALL PAGES
    label = '@j(tot,oh,nb,bs) [MA/m$+2$=] vs rho@'
    CALL GRD1D(1,rhog,vgx1,nrmax+1,nrmax+1,4,label,0)
    label = '@d Psi/d rho vs rho@'
    CALL GRD1D(2,rhog,vgx2,nrmax+1,nrmax+1,1,label,0)
    label = '@qp vs rho@'
    CALL GRD1D(3,rhog,vgx3,nrmax+1,nrmax+1,1,label,0)
    label = '@eta(para) [ohm m] vs rho@'
    CALL GRD1D(4,rhog,vgx4,nrmax+1,nrmax+1,1,label,gsid)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_rad5

! **************************************************************************
  SUBROUTINE tr_gr_rad6
  ! ----- heating profile-----
    USE trcomm, ONLY: str,poh,pnb,pec,pic,plh,pnf,prl

    CALL tr_gr_init_vgx

    vgx1(0:nrmax,1) = poh(1,0:nrmax)*1.d-6
    vgx1(0:nrmax,2) = (pnb(1,0:nrmax)+pnb(2,0:nrmax))*1.d-6
    vgx1(0:nrmax,3) = (pec(1,0:nrmax)+pec(2,0:nrmax))*1.d-6
    vgx1(0:nrmax,4) = (pic(1,0:nrmax)+pic(2,0:nrmax))*1.d-6
    vgx1(0:nrmax,5) = (plh(1,0:nrmax)+plh(2,0:nrmax))*1.d-6

    vgx2(0:nrmax,1) = poh(1,0:nrmax)*1.d-6
    vgx2(0:nrmax,2) = (prl(1,0:nrmax)+prl(2,0:nrmax))*1.d-6
    vgx2(0:nrmax,3) = (pnf(1,0:nrmax)+pnf(2,0:nrmax))*1.d-6


    CALL PAGES
    label = '@Pin(oh,nb,ec,ic,lh) [MW/m$+3$=] vs rho@'
    CALL GRD1D(1,rhog,vgx1,nrmax+1,nrmax+1,5,label,0)
    label = '@Pin(oh,rl,nf) [MW/m$+3$=] vs rho@'
    CALL GRD1D(2,rhog,vgx2,nrmax+1,nrmax+1,3,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_rad6

! **************************************************************************
  SUBROUTINE tr_gr_rad7
  ! rotation velocity profile
    USE trcomm, ONLY: vtor,vpol,vpar,vprp

    CALL tr_gr_init_vgx

    vgx1(0:nrmax,1) = vtor(0:nrmax)
    vgx2(0:nrmax,1) = vpol(0:nrmax)
    vgx3(0:nrmax,1) = vpar(0:nrmax)
    vgx4(0:nrmax,1) = vprp(0:nrmax)

    CALL PAGES
    label = '@Vtor[m/s] vs rho@'
    CALL GRD1D(1,rhog,vgx1,nrmax+1,nrmax+1,1,label,0)
    label = '@Vpol[m/s] vs rho@'
    CALL GRD1D(2,rhog,vgx2,nrmax+1,nrmax+1,1,label,0)
    label = '@Vpar[m/s] vs rho@'
    CALL GRD1D(3,rhog,vgx3,nrmax+1,nrmax+1,1,label,0)
    label = '@Vprp[m/s] vs rho@'
    CALL GRD1D(4,rhog,vgx4,nrmax+1,nrmax+1,1,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_rad7

! **************************************************************************
  SUBROUTINE tr_gr_rad8
    USE trcomm, ONLY: er,vexbp,dvexbpdr,wexbp

    CALL tr_gr_init_vmx

    vmx1(1:nrmax,1) = er(1:nrmax)
    vmx2(1:nrmax,1) = vexbp(1:nrmax)
    vmx3(1:nrmax,1) = dvexbpdr(1:nrmax)
    vmx4(1:nrmax,1) = wexbp(1:nrmax)

    CALL PAGES
    label = '@Er [V/m] vs rho@'
    CALL GRD1D(1,rhomg,vmx1,nrmax,nrmax,1,label,0)
    label = '@Vexb [1/s] vs rho@'
    CALL GRD1D(2,rhomg,vmx2,nrmax,nrmax,1,label,0)
    label = '@dpvexbp vs rho@'
    CALL GRD1D(3,rhomg,vmx3,nrmax,nrmax,1,label,0)
    label = '@Wexb [1/s] vs rho@'
    CALL GRD1D(4,rhomg,vmx4,nrmax,nrmax,1,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_rad8

! **************************************************************************
  SUBROUTINE tr_gr_rad9
    ! fast ion particle profile
    USE trcomm, ONLY: rt,rn,pnb,nsafmax
    INTEGER(ikind) :: nsafmaxg

    CALL tr_gr_init_vg
    CALL tr_gr_init_vgx

    IF(nsafmax > 0)THEN
       nsafmaxg = nsafmax
       DO neq = 1, neqmax
          nsa = nsa_neq(neq)
          IF(nsa /= 0)THEN
             nsaf = nsab_nsaf(nsa)
             IF(nsaf /= 0)THEN
                vg1(0:nrmax,nsaf) = rt(nsa,0:nrmax)
                vg2(0:nrmax,nsaf) = rn(nsa,0:nrmax)
             END IF
          END IF
       END DO
    ELSE
       nsafmaxg = 1
    END IF
       
    vgx1(0:nrmax,1) = (pnb(1,0:nrmax)+pnb(2,0:nrmax))*1.d-6

    CALL PAGES
    label = '@Ti(fast) [keV] vs rho@'
    CALL GRD1D(1,rhog,vg1,nrmax+1,nrmax+1,nsafmaxg,label,0)
    label = '@ni(fast) [10$+20$=/m$+3$=] vs rho@'
    CALL GRD1D(2,rhog,vg2,nrmax+1,nrmax+1,nsafmaxg,label,0)
    label = '@Pin(fast) [MW/m$+3$=] vs rho@'
    CALL GRD1D(3,rhog,vgx1,nrmax+1,nrmax+1,nsafmaxg,label,0)
    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_rad9

! **************************************************************************
  SUBROUTINE tr_gr_rad11
  ! ----- history of radial profile -----
    USE trcomm, ONLY: ngt,gvrt,gvrts

    CALL tr_gr_init_gg

    IF(ngt > 0)THEN
       ngg_interval = ngt/(MOD(ngt-1,nggmax)+1)
    ELSE IF(ngt <= 0 )THEN
       ngg_interval = 1
    END IF
    DO ngg = 0, nggmax
       gg1(0:nrmax,ngg) = gvrts(0:nrmax, ngg*ngg_interval, 1,3)
       gg2(0:nrmax,ngg) = gvrts(0:nrmax, ngg*ngg_interval, 2,3)
       gg3(0:nrmax,ngg) = gvrts(0:nrmax, ngg*ngg_interval, 1,1)
       gg4(0:nrmax,ngg) =  gvrt(0:nrmax, ngg*ngg_interval, 1)
    END DO

    CALL PAGES
    label = '@T1(t) [keV] vs rho@'
    CALL GRD1D(1,rhog,gg1,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '@T2(t) [keV] vs rho@'
    CALL GRD1D(2,rhog,gg2,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '@n1(t) [10$+20$=/m$+3$=] vs rho@'
    CALL GRD1D(3,rhog,gg3,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '@qp(t) vs rho@'
    CALL GRD1D(4,rhog,gg4,nrmax+1,nrmax+1,nggmax+1,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_rad11

! **************************************************************************
  SUBROUTINE tr_gr_rad15
  ! ----- history of radial profile -----
    USE trcomm, ONLY: ngt,gvrt

    CALL tr_gr_init_gg

    IF(ngt > 0)THEN
       ngg_interval = ngt/(MOD(ngt-1,nggmax)+1)
    ELSE IF(ngt <= 0 )THEN
       ngg_interval = 1
    END IF
    DO ngg = 0, nggmax
       gg1(0:nrmax,ngg) = 1.d-6*gvrt(0:nrmax, ngg*ngg_interval, 2)
       gg2(0:nrmax,ngg) = 1.d-6*gvrt(0:nrmax, ngg*ngg_interval, 3)
!       gg3(0:nrmax,ngg) = gvrt(0:nrmax, ngg*ngg_interval, 4)
       ! history of qp profile
       gg4(0:nrmax,ngg) = gvrt(0:nrmax, ngg*ngg_interval, 1)
    END DO

    CALL PAGES
    label = '@j_tot(t) [MA/m$+2$=] vs rho@'
    CALL GRD1D(1,rhog,gg1,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '@j_oh(t) [MA/m$+2$=] vs rho@'
    CALL GRD1D(2,rhog,gg2,nrmax+1,nrmax+1,nggmax+1,label,0)
!    label = '@j_ex(t) vs rho@'
!    CALL GRD1D(3,rhog,gg3,nrmax+1,nrmax+1,nggmax+1,label,0)
    label = '@qp(t) vs rho@'
    CALL GRD1D(4,rhog,gg4,nrmax+1,nrmax+1,nggmax+1,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_rad15

END MODULE trgrad
