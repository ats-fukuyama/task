MODULE trresult

  USE trcomm, ONLY: rkind,ikind

  PUBLIC tr_calc_global,tr_status,tr_save_ngt
  PRIVATE

CONTAINS

! ***** simple status report *****

  SUBROUTINE tr_status

    USE trcomm, ONLY : nsamax,ns_nsa,t,qp,rt,kidns,nitmax, &
         wp_t,taue1,taue2,taue3
    IMPLICIT NONE
    INTEGER(ikind):: nsa
    CHARACTER(LEN=64) :: fmt_1,fmt_2,fmt_3

    fmt_1 = '(A5,F7.3,A12,F7.2,A5,A7,F7.3,A12,F7.3)'
    fmt_2 = '(A3,A1,A1)'
    fmt_3 = '(F7.3,A7)'

    WRITE(6,fmt_1) ' # T:',t,'(s)      WP:',wp_t,'(MJ)', &
                   '  TAUE:',taue3,'(s)     Q0:',qp(0)

    DO nsa = 1, nsamax
       WRITE(6,fmt_2,ADVANCE='NO') '  T',kidns(ns_nsa(nsa)),': '
       WRITE(6,fmt_3,ADVANCE='NO') rt(nsa,0),'(keV)  '
    END DO
    WRITE(6,*) ! breaking line

    WRITE(6,'(A14,I3)') '  Iterations: ',nitmax
    nitmax = 0
    
    RETURN
  END SUBROUTINE tr_status

! ***** calculate global values *****

  SUBROUTINE tr_calc_global
    USE trcomm, ONLY: pi,rmu0,rkev,RR,ra,BB,nrmax,nsamax,nsabmax,ns_nsa, &
         rkap,pa,pz,rhog,bp,dvrho,rip,t,t_prev,dt,rn,rt,rn_prev,rt_prev, &
         beta,beta_va,betap,betap_va,betaq,betan,                        &
         poh,pnb,pec,pic,plh,pibw,pnf,prl,                               &
         pin_t,poh_t,pnb_t,pec_t,pic_t,plh_t,pibw_t,prl_t,pnf_t,         &
         rns_va,rts_va,ws_t,wp_t,wp_th,taue1,taue2,taue3,                &
         taue89,taue98,h89,h98y2,                                        &
         mdluf,mdleqn,mdlequ,mdleqt,mdleqm,mdlsrc,mdlgmt,mdlglb

    IMPLICIT NONE
    REAL(rkind) :: FCTR
    REAL(rkind) :: pvol,sumnt,sumntm,sumntq,rnsum,rntsum
    REAL(rkind) :: rns_vasum,rn_la,rne_la
    REAL(rkind) :: wp_prev,dwpdt,ami
    REAL(rkind),DIMENSION(0:nrmax) :: dvrhom,drhog,dsrho,dvm

    INTEGER(ikind) :: nr,ns,nsa

! ------------------------------------------------------------------------
!   * Plasma volume
!   * Local beta
!      - beta     : toroidal beta
!      - beta_va  : volume-averaged toroidal beta
!      - betap    : poloidal beta
!      - betap_va : volume-averaged poloidal beta
!      - betaq    : toroidal beta for reaction rate
!                    (ref: Tokamaks 3rd, p115)
    pvol   = 0.d0
    sumntm = 0.d0
    sumntq = 0.d0
    DO nr = 1, nrmax
       drhog(nr)  = rhog(nr)-rhog(nr-1)
       dvrhom(nr) = 0.5d0*(dvrho(nr-1)+dvrho(nr))
       
       sumnt = rkev*1.d20 &
               *SUM(0.5d0*(rn(1:nsamax,nr-1)*rt(1:nsamax,nr-1)   &
                          +rn(1:nsamax,nr  )*rt(1:nsamax,nr  )))
       
       pvol = pvol + dvrhom(nr)*drhog(nr)
       sumntm = sumntm + sumnt    * dvrhom(nr)*drhog(nr)
       sumntq = sumntq + sumnt**2 * dvrhom(nr)*drhog(nr)

       beta(nr)     = 2.d0*rmu0*sumntm / (pvol*BB**2)
       beta_va(nr)  = 2.d0*rmu0*sumnt  / (     BB**2)
       betap(nr)    = 2.d0*rmu0*sumntm / (pvol*BP(nrmax)**2)
       betap_va(nr) = 2.d0*rmu0*sumnt  / (     BP(nrmax)**2)
       betaq(nr)    = 2.d0*rmu0*SQRT(sumntq)/(SQRT(pvol)*BB**2)
    END DO
    beta(0)  = FCTR(rhog(1),rhog(2),beta(1),beta(2))
    betap(0) = FCTR(rhog(1),rhog(2),betap(1),betap(2))
    betaq(0) = FCTR(rhog(1),rhog(2),betaq(1),betaq(2))

    !   * Global beta
    !      - beta_a  : toroidal beta at separatrix
    !      - betap_a : poloidal beta at separatrix
    !      - betan   : normalized toroidal beta (Troyon beta)
    betan = beta(nrmax)*1.d2 / (rip/(ra*BB))

! ------------------------------------------------------------------------

!   * Volume-averaged density and temperature

!   * Line-averaged densities
    rn_la = 0.d0
    DO nsa = 1, nsamax
       rn_la = rn_la &
              +SUM(0.5d0*(rn(nsa,0:nrmax-1)+rn(nsa,1:nrmax))*drhog(1:nrmax))
    END DO
    ! only electron
    rne_la = SUM(0.5d0*(rn(1,0:nrmax-1)+rn(1,1:nrmax))*drhog(1:nrmax))

! ------------------------------------------------------------------------
!   * Ohmic, NBI and fusion powers
    DO
       dvm(1:nrmax) = dvrhom(1:nrmax)*drhog(1:nrmax)*1.d-6

       poh_t  = 0.d0
       DO ns = 1, 2
          poh_t = poh_t &
                + SUM(0.5d0*(poh(ns,0:nrmax-1)+poh(ns,1:nrmax))*dvm(1:nrmax))
       END DO

       IF(mdluf==2 .AND. mdlsrc==7) EXIT ! experimental data are already set
       pnb_t  = 0.d0
       pec_t  = 0.d0
       pic_t  = 0.d0
       plh_t  = 0.d0
       pibw_t = 0.d0
       prl_t  = 0.d0
       pnf_t  = 0.d0    
       DO ns = 1, 2 ! electron and ion
          pnb_t = pnb_t &
                + SUM(0.5d0*(pnb(ns,0:nrmax-1)+pnb(ns,1:nrmax))*dvm(1:nrmax))
          pec_t = pec_t &
                + SUM(0.5d0*(pec(ns,0:nrmax-1)+pec(ns,1:nrmax))*dvm(1:nrmax))
          pic_t = pic_t &
                + SUM(0.5d0*(pic(ns,0:nrmax-1)+pic(ns,1:nrmax))*dvm(1:nrmax))
          plh_t = plh_t &
                + SUM(0.5d0*(plh(ns,0:nrmax-1)+plh(ns,1:nrmax))*dvm(1:nrmax))
          pibw_t = pibw_t &
               + SUM(0.5d0*(pibw(ns,0:nrmax-1)+pibw(ns,1:nrmax))*dvm(1:nrmax))
          prl_t = prl_t &
               + SUM(0.5d0*(prl(ns,0:nrmax-1)+prl(ns,1:nrmax))*dvm(1:nrmax))
          pnf_t = pnf_t &
               + SUM(0.5d0*(pnf(ns,0:nrmax-1)+pnf(ns,1:nrmax))*dvm(1:nrmax))
       END DO
       EXIT
    END DO

!   * Radiation, charge exchage and ionization losses

! ------------------------------------------------------------------------
!   * Stored energy
    DO nsa = 1, nsamax
       rnsum  = 0.d0
       rntsum = 0.d0
       DO nr = 1, nrmax
          rnsum  = rnsum  + 0.5d0*(rn(nsa,nr-1)+rn(nsa,nr)) &
                                 *dvrhom(nr)*drhog(nr)
          rntsum = rntsum + 0.5d0*(rn(nsa,nr-1)*rt(nsa,nr-1)  &
                                  +rn(nsa,nr  )*rt(nsa,nr  )) &
                                 *dvrhom(nr)*drhog(nr)
       END DO

       rns_va(nsa) = rnsum / pvol
       IF(rnsum == 0.d0)THEN
          rts_va(nsa) = 0.d0
       ELSE
          rts_va(nsa) = rntsum / rnsum
       END IF

       ws_t(nsa) = 1.5d0*rntsum*rkev*1.d14
    END DO

!   * Fusion production rate

!   * Output source and power

    !   * Input and output sources and powers
    wp_t  = SUM(ws_t(1:nsamax))  ! including fast ions
    wp_th = SUM(ws_t(1:nsabmax)) ! excluding fast ions
    pin_t  = poh_t + pnb_t + pec_t + pic_t + plh_t + pnf_t - prl_t
!    pout_t =
!    sin_t  = 
!    sout_t =

!   * Ionization, fusion and NBI fuelling
!   * Pellet injection fuelling

! ------------------------------------------------------------------------
!   * Energy confinement times
!      - taue1 : steady state
!      - taue2 : transient
!      - taue3 : thermal energy confinement time (transient)
    rntsum = 0.d0
    DO nsa = 1, nsamax
       DO nr = 1, nrmax
          rntsum = rntsum + 0.5d0*(rn_prev(nsa,nr-1)*rt_prev(nsa,nr-1)  &
                                  +rn_prev(nsa,nr  )*rt_prev(nsa,nr  )) &
                                 *dvrhom(nr)*drhog(nr)
       END DO
    END DO
    wp_prev = 1.5d0*rntsum*rkev*1.d14

    IF(ABS(t-t_prev) <= 1.d-70)THEN
       dwpdt = 0.d0
    ELSE
       dwpdt = (wp_t-wp_prev) / (t-t_prev)
!       write(*,*) wp_t, wp_prev, rn(1,10), rn_prev(1,10)
    ENDIF
    t_prev  = t

    IF(t_prev == 0) THEN
       taue1 = wp_t / pin_t
       taue2 = taue1
       taue3 = taue1
    ELSE
       taue1 = wp_t / pin_t
       taue2 = wp_t  / (pin_t-dwpdt)
       taue3 = wp_th / (pin_t-dwpdt)
    END IF

!   * Confinement scaling
!    taue89: ITER89-P L-mode scaling
!      ** P.N. Yushmanov et al 1990 Nucl. Fusion 30 1999
!    taue98: IPB89(y,2) H-mode scaling with ELMs
!      ** ITER Physics Expert Group on Confinement and Transport et al 
!          1999 Nucl. Fusion 39 2175
!    h89   : confinement enhancement factor for ITER89-P
!    h98y2 : confinement enhancement factor for IPB98(y,2)
    
    rns_vasum = 0.d0
    ami       = 0.d0
    DO nsa = 1, nsamax
       ns = ns_nsa(nsa)
       IF(pz(ns) > 0)THEN
          ami = ami + pa(ns)*rns_va(nsa)
          rns_vasum = rns_vasum + rns_va(nsa)
       END IF
    END DO
    ami = ami / rns_vasum

!    [ ITER89-P ]
    taue89 = 4.8d-2 * (ABS(rip)**0.85d0) * (RR**1.2d0) * (ra**0.3d0) &
             * (rkap**0.5d0) * (rn_la**0.1d0)                        &
             * (ABS(BB)**0.2d0) * (ami**0.5d0) * (pin_t**(-0.5d0))

    h89 = taue3 / taue89

!    [ IPB98(y,2) ]
    taue98 = 5.62d-2 * (ABS(rip)**0.93d0) * (RR**1.39d0) * (ra**0.58d0) &
             * (rkap**0.78d0) * ((rne_la*10.d0)**0.41d0)                &
             * (ABS(BB)**0.15d0) * (ami**0.19d0) * (pin_t**(-0.69d0))
!     *** Tokamaks 3rd p.184 ***
!    taue98 = 0.145 * (ABS(rip)**0.93d0) * (RR**1.39d0) * (ra**0.58d0) &
!             * (rkap**0.78d0) * (rn_la**0.41d0)                       &
!             * (ABS(BB)**0.15d0) * (ami**0.19d0) * (pin_t**(-0.69d0))

    h98y2 = taue3 / taue98

! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
!   * Inductance and loop voltage
    !   * Current densities
    dsrho(0:nrmax) = dvrho(0:nrmax) / (2.d0*pi*RR)

!   * Distance of q=1 surface from magnetic axis
!   * Effective charge number at axis
    

    RETURN
  END SUBROUTINE tr_calc_global

! ***** save data for time history *****

  SUBROUTINE tr_save_ngt

    USE trcomm, ONLY : &
         nrmax,nsamax,ngtmax,neqmax,ngt,gvt,gvts,gvrt,gvrts,        &
         rkev,t,rn,ru,rt,rp,qp,jtot,joh,rip,wp_t,ws_t,              &
         jcd_nb,jcd_ec,jcd_ic,jcd_lh,                               &
         pin_t,poh_t,pnb_t,pec_t,pic_t,plh_t,pibw_t,pnf_t,          &
         beta,betap,betan,taue1,taue3,taue89,taue98,h89,h98y2
    USE trcoeftb, ONLY: Pereverzev_check
    IMPLICIT NONE
    INTEGER(ikind):: nsa,nr,wstmax
    REAL(rkind),DIMENSION(neqmax,0:nrmax):: add_prv

    IF(ngt >= ngtmax) RETURN

    ngt=ngt+1
    gvt(ngt, 0) = t

    ! *****   0D value --- values on the axis or the edge   *****
    gvt(ngt, 1) = qp(0)
    gvt(ngt, 2) = qp(nrmax)

    ! plasma curent
    gvt(ngt, 3) = rip

    ! stored energy
    gvt(ngt,8) = wp_t
    
    IF(nsamax >  4) wstmax = 4
    IF(nsamax <= 4) wstmax = nsamax
    DO nsa = 1, wstmax
       gvt(ngt,8+nsa) = ws_t(nsa)
    END DO

    ! beta
    gvt(ngt,13) = beta(0)
    gvt(ngt,14) = beta(nrmax)
    gvt(ngt,15) = betap(0)
    gvt(ngt,16) = betap(nrmax)
    gvt(ngt,17) = betan

    gvt(ngt,18) = pin_t
    gvt(ngt,19) = poh_t
    gvt(ngt,20) = pnb_t
    gvt(ngt,21) = pec_t
    gvt(ngt,22) = pic_t
    gvt(ngt,23) = plh_t
    gvt(ngt,24) = pnf_t

    ! confinement time
    gvt(ngt,25) = taue3
    gvt(ngt,26) = taue89
    gvt(ngt,27) = taue98
    gvt(ngt,28) = h89
    gvt(ngt,29) = h98y2

! ------------------------------------------------------------------------

    DO nsa=1,nsamax
       gvts(ngt,nsa, 1) = rn(nsa,0)
       gvts(ngt,nsa, 2) = ru(nsa,0)
       gvts(ngt,nsa, 3) = rt(nsa,0)
       gvts(ngt,nsa, 4) = rp(nsa,0)
    END DO

! ------------------------------------------------------------------------

    gvrt(0:nrmax,ngt,1) = qp(0:nrmax)
    gvrt(0:nrmax,ngt,2) = jtot(0:nrmax)
    gvrt(0:nrmax,ngt,3) = joh(0:nrmax)
    gvrt(0:nrmax,ngt,4) = jcd_nb(0:nrmax)
    gvrt(0:nrmax,ngt,5) = jcd_ec(0:nrmax)
    gvrt(0:nrmax,ngt,6) = jcd_ic(0:nrmax)
    gvrt(0:nrmax,ngt,7) = jcd_lh(0:nrmax)

    ! *****   radial profile   *****
    DO nsa = 1, nsamax
          gvrts(0:nrmax,ngt,nsa, 1) = rn(nsa,0:nrmax)
          gvrts(0:nrmax,ngt,nsa, 2) = ru(nsa,0:nrmax)
          gvrts(0:nrmax,ngt,nsa, 3) = rt(nsa,0:nrmax)
          gvrts(0:nrmax,ngt,nsa, 4) = rp(nsa,0:nrmax)
    END DO
!!$    ! experimental data
!!$    DO nsu = 1, nsum
!!$       IF(idnm(nsu))THEN
!!$          nsa = nsa_nsu(nsu)
!!$          IF(nsa /= 0)THEN
!!$             gvrtu(0:nrmax,ngt,nsa,1) = rnug(nsu,1:nrmax+1)
!!$!             gvrtu(0:nrmax,ngt,nsa,2) = ruug(nsu,1:nrmax+1)
!!$             gvrtu(0:nrmax,ngt,nsa,3) = rtug(nsu,1:nrmax+1)
!!$             gvrtu(0:nrmax,ngt,nsa,4) = &
!!$                  (rnug(nsu,1:nrmax+1)+rtug(nsu,1:nrmax+1))*rkev*1.d20
!!$          END IF
!!$       END IF
!!$    END DO
!!$    DO nsfu = 1, nsum
!!$       IF(idnfast(nsu))THEN
!!$          nsa = nsa_nsu(nsu)
!!$          IF(nsa /= 0)THEN
!!$             gvrtu(0:nrmax,ngt,nsa,1) = rnfug(nsu,1:nrmax+1)
!!$!             gvrtu(0:nrmax,ngt,nsa,3) = rtfug(nsu,1:nrmax+1)
!!$!             gvrtu(0:nrmax,ngt,nsa,4) = &
!!$!                  (rnfug(nsu,1:nrmax+1)+rtfug(nsu,1:nrmax+1))*rkev*1.d20
!!$          END IF
!!$       END IF
!!$    END DO

    ! for Pereverzev method
    ! numerically addtional term in nodal equation (relative value)
    IF(t == 0)THEN
       gvrts(0:nrmax,ngt,1:nsamax,5) = 0.d0
       gvrts(0:nrmax,ngt,1:nsamax,6) = 0.d0
       gvrts(0:nrmax,ngt,1:nsamax,7) = 0.d0
    ELSE
       CALL Pereverzev_check(add_prv)
       DO nsa=1,nsamax
          DO nr=0,nrmax
             gvrts(nr,ngt,nsa,5) = add_prv(1+3*nsa-2,nr)
             gvrts(nr,ngt,nsa,6) = add_prv(1+3*nsa-1,nr)
             gvrts(nr,ngt,nsa,7) = add_prv(1+3*nsa  ,nr)
          END DO
       END DO
    END IF

    RETURN
  END SUBROUTINE tr_save_ngt

END MODULE trresult
