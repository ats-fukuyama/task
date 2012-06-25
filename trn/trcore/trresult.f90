MODULE trresult

  USE trcomm, ONLY: rkind,ikind

  PUBLIC tr_calc_global,tr_status,tr_save_ngt
  PRIVATE

CONTAINS

! ***** calculate global values *****

  SUBROUTINE tr_calc_global
    USE trcomm, ONLY: pi,rmu0,rkev,RR,ra,BB,nrmax,nsamax, &
         rhog,bp,dvrho,rip,rn,rt,                         &
         beta,beta_va,betap,betap_va,betaq,betan,         &
         rns_va,rts_va,ws_t,wp_t,taue1,taue2

    IMPLICIT NONE
    REAL(rkind) :: FCTR
    REAL(rkind) :: pvol,sumnt,sumntm,sumntq,rnsum,rntsum
    REAL(rkind),DIMENSION(0:nrmax) :: dvrhom,drhog,dsrho

    INTEGER(ikind) :: nr,nsa

    !   * Plasma volume
    !   * Local beta
    !      - beta     : toroidal beta
    !      - beta_va  : volume-averaged toroidal beta
    !      - betap    : poloidal beta
    !      - betap_va : colume-averaged poloidal beta
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


    !   * Volume-averaged density and temperature
    !   * Stored energy
    rnsum  = 0.d0
    rntsum = 0.d0
    DO nsa = 1, nsamax
       DO nr = 1, nrmax
          rnsum  = rnsum  + 0.5d0*(rn(nsa,nr-1)+rn(nsa,nr)) &
                                 *dvrhom(nr)*drhog(nr)
          rntsum = rntsum + 0.5d0*(rn(nsa,nr-1)*rt(nsa,nr-1)  &
                                  +rn(nsa,nr  )*rt(nsa,nr  )) &
                                 *dvrhom(nr)*drhog(nr)
       END DO

       rns_va(nsa) = rnsum / pvol
       rts_va(nsa) = rntsum / rnsum

       ws_t(nsa) = 1.5d0*rntsum*rkev*1.d14
    END DO
       
!   * Line-averaged densities ???

!   * Ohmic, NBI and fusion powers
!   * External power typically for NBI from exp. data
!   * RF power
!   * Total NBI power distributed on electrons and bulk ions
!   * Total RF  power distributed on electrons and bulk ions
!   * Radiation, charge exchage and ionization losses

!   * Inductance and loop voltage
    !   * Current densities
       dsrho(0:nrmax) = dvrho(0:nrmax) / (2.d0*pi*RR)

!   * Fusion production rate
!   * Output source and power

    !   * Input and output sources and powers
    wp_t = SUM(ws_t(1:nsamax))

!   * Ionization, fusion and NBI fuelling
!   * Pellet injection fuelling

    !   * Energy confinement times
    !      - taue1 : steady state
    !      - taue2 : transient

!   * Confinement scaling (Tokamaks 3rd p183,184)



!   * Distance of q=1 surface from magnetic axis
!   * Effective charge number at axis
    

    RETURN
  END SUBROUTINE tr_calc_global

! ***** simple status report *****

  SUBROUTINE tr_status

    USE trcomm, ONLY : nsamax,t,qp,rt,kidnsa,nitmax, &
         wp_t
    IMPLICIT NONE
    REAL(rkind):: taue
    INTEGER(ikind):: nsa

    taue=0.d0

    WRITE(6,601) t,wp_t,taue,qp(0)
601 FORMAT(' # T: ',F7.3,'(s)     WP:',F7.2,'(MJ)  ', &
           '  TAUE:',F7.3,'(s)   Q0:',F7.3)
    WRITE(6,602) (kidnsa(nsa),rt(nsa,0),nsa=1,nsamax)
602 FORMAT(4(' T',A1,':',F7.3,'(keV)  ':))

    WRITE(6,'(A13,I4)') 'Iterations: ',nitmax
    nitmax = 0
    
    RETURN
  END SUBROUTINE tr_status

! ***** save data for time history *****

  SUBROUTINE tr_save_ngt

    USE trcomm, ONLY : &
         nrmax,nsamax,ngtmax,neqmax,ngt,gvt,gvts,gvrt,gvrts,  &
         t,rn,ru,rt,qp,jtot,joh,htr,rip,wp_t
    USE trcoeftb, ONLY: Pereverzev_check
    IMPLICIT NONE
    INTEGER(ikind):: nsa,nr
    REAL(rkind),DIMENSION(neqmax,0:nrmax):: add_prv

    IF(ngt >= ngtmax) RETURN

    ngt=ngt+1
    gvt(ngt, 0) = t

    ! ----- values on the axis or the edge -----
    gvt(ngt, 1) = qp(0)
    gvt(ngt, 2) = qp(nrmax)

    ! plasma curent
    gvt(ngt, 3) = rip

    ! stored energy
    gvt(ngt,8) = wp_t

    DO nsa=1,nsamax
       gvts(ngt,nsa, 1) = rn(nsa,0)
       gvts(ngt,nsa, 2) = ru(nsa,0)
       gvts(ngt,nsa, 3) = rt(nsa,0)
    END DO

    ! ----- radial profile -----
    DO nr=0,nrmax
       gvrt(nr,ngt, 1) = qp(nr)
    END DO

    DO nsa=1,nsamax
          gvrts(0:nrmax,ngt,nsa, 1) = rn(nsa,0:nrmax)
          gvrts(0:nrmax,ngt,nsa, 2) = ru(nsa,0:nrmax)
          gvrts(0:nrmax,ngt,nsa, 3) = rt(nsa,0:nrmax)
    END DO

    gvrt(0:nrmax,ngt,2) = jtot(0:nrmax) + htr(1,0:nrmax)
    gvrt(0:nrmax,ngt,3) = joh(0:nrmax)
    gvrt(0:nrmax,ngt,4) = htr(1,0:nrmax)


    ! for Pereverzev method
    ! numerically addtional term in nodal equation (relative value)
    IF(t == 0)THEN
       gvrts(0:nrmax,ngt,1:nsamax,4) = 0.d0
       gvrts(0:nrmax,ngt,1:nsamax,5) = 0.d0
       gvrts(0:nrmax,ngt,1:nsamax,6) = 0.d0
    ELSE
       CALL Pereverzev_check(add_prv)
       DO nsa=1,nsamax
          DO nr=0,nrmax
             gvrts(nr,ngt,nsa,4) = add_prv(1+3*nsa-2,nr)
             gvrts(nr,ngt,nsa,5) = add_prv(1+3*nsa-1,nr)
             gvrts(nr,ngt,nsa,6) = add_prv(1+3*nsa  ,nr)
          END DO
       END DO
    END IF

    RETURN
  END SUBROUTINE tr_save_ngt

END MODULE trresult
