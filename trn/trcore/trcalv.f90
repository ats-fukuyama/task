MODULE trcalv
!---------------------------------------------------------------------------
!
!   This module calculates variables to evaluate transport coefficients
!    and other effects. 
!
!---------------------------------------------------------------------------
  USE trcomm, ONLY: rkind,ikind,nrmax,nsamax

  PRIVATE
  PUBLIC tr_calc_variables,tr_calc_zeff,tr_calc_clseta,coulog

CONTAINS

  SUBROUTINE tr_calc_variables
!--------------------------------------------------------------------------
!   Calculate variables and some effects used in calculate transport 
!--------------------------------------------------------------------------
    USE trcomm, ONLY: nrmax,mdltr_nc,nrd1,nrd2,nrd3,nrd4

    IMPLICIT NONE
!    INTEGER(ikind) :: nr
!    REAL(rkind) :: deriv3 ! the 2rd order accuracy derivative in TASK/lib
!    REAL(rkind) :: deriv4 ! the 3rd order accuracy derivative in TASK/lib

    CALL tr_calc_grad

    CALL tr_calc_er

    CALL tr_calc_exb

    CALL tr_calc_zeff

    IF(mdltr_nc==0) CALL tr_calc_clseta

       ! Doppler shift
       ! AGMP(NR) = QP(NR)/EPS*WEXB(NR)

       ! sound speed for electron
!       v_se  = 0.d0
       ! pressure gradient for MHD instability

!    nrd1(0:nrmax) = er(0:nrmax)
!    nrd2(0:nrmax) = vexbp(0:nrmax)
!    nrd3(0:nrmax) = dvexbpdr(0:nrmax)
!    nrd4(0:nrmax) = bp(0:nrmax)

    RETURN
  END SUBROUTINE tr_calc_variables

!==========================================================================

  SUBROUTINE tr_calc_grad
!-------------------------------------------------------------------------
!   The derivatives are calculated on the each grid.
!-------------------------------------------------------------------------
    USE trcomm, ONLY: rkev,rmu0,pa,pz,pz0,idnsa,ns_nsa, &
         RR,BB,rhog,ar1rho,rmnrho,rt,rn,bp,abb1rho,qp,  &
         rn_i,rn_id,rn_icl,rn_e,rn_ed,rn_ecl,     &
         rt_i,rt_id,rt_icl,rt_e,rt_ed,rt_ecl,     &                     
         rp,rp_d,rp_tot,rp_totd,qp_d,             &
         mshear,mcurv,ai_ave,alpha
!         nrd1,nrd2,nrd3

    IMPLICIT NONE

    ! --- control and internal variables
    INTEGER(ikind) :: nr,ns,nsa
    REAL(rkind) :: rt_isum,pansum,dr

!    *****
!    ! This part will be implemented in introducing the NBI heating modules.
!    *****
!    additional pressure due to NBI
!
!    IF(SUMPBM.EQ.0.D0) THEN
!       PADD(1:NRMAX)=0.D0
!    ELSE
!       ! density [10^20 m^-3] * temperature [keV]
!       PADD(1:NRMAX)=PBM(1:NRMAX)*1.D-20/RKEV-RNF(1:NRMAX,1)*RT(1:NRMAX,2)
!                                                             RTF(nr, 1)
!       ! padd = experimental - calculating result (beam pressure)
!    ENDIF


    ! initilization (zero at NR = 0) ----
    rn_i(0:nrmax)    = 0.d0

    rp_d(1:nsamax,0:nrmax)   = 0.d0
    rp_totd(0:nrmax) = 0.d0
    rn_ed(0:nrmax)  = 0.d0
    rn_ecl(0:nrmax) = 0.d0
    rn_id(0:nrmax)  = 0.d0
    rn_icl(0:nrmax) = 0.d0
    rt_ed(0:nrmax)  = 0.d0
    rt_ecl(0:nrmax) = 0.d0
    rt_id(0:nrmax)  = 0.d0
    rt_icl(0:nrmax) = 0.d0
    qp_d(0:nrmax)   = 0.d0
    mshear(0:nrmax) = 0.d0
    mcurv(0:nrmax)  = 0.d0

    ! '_i' means 'hydrogenic ions'. 
    ! for now, in following, 'nsa' is used as the loop counter.
    DO nr = 0, nrmax
       rt_isum = 0.d0
       
       rt_e(nr) = rt(1,nr)
       rn_e(nr) = rn(1,nr)
       DO nsa = 1, nsamax
          IF(idnsa(nsa) == 1) THEN ! ion
             rn_i(nr) = rn_i(nr) + rn(nsa,nr)
             rt_isum  = rt_isum + rn(nsa,nr)*rt(nsa,nr)
          END IF            
       END DO         
       rt_i(nr)    = rt_isum / rn_i(nr)
    END DO


    ! derivatives with respect to 'r' on half grids
    !  ( d XX/dr = <|grad rho|>*d XX/d rho)
    ! need consideration for PADD
    DO nr = 1, nrmax
       dr = rmnrho(nr) - rmnrho(nr-1)
       DO nsa = 1, nsamax
          rp_d(nsa,nr) = (rp(nsa,nr)-rp(nsa,nr-1)) / dr
       END DO

       rn_ed(nr)   = (rn_e(nr)   -   rn_e(nr-1)) / dr
       rn_id(nr)   = (rn_i(nr)   -   rn_i(nr-1)) / dr
       rt_ed(nr)   = (rt_e(nr)   -   rt_e(nr-1)) / dr
       rt_id(nr)   = (rt_i(nr)   -   rt_i(nr-1)) / dr
       rp_totd(nr) = (rp_tot(nr) - rp_tot(nr-1)) / dr

       ! scale length ( (d XX/d rho)/XX ) on half grids
       rt_ecl(nr)  = rt_ed(nr) / (0.5d0*(rt_e(nr)+rt_e(nr-1)))
       rt_icl(nr)  = rt_id(nr) / (0.5d0*(rt_i(nr)+rt_i(nr-1)))
       rn_ecl(nr)  = rn_ed(nr) / (0.5d0*(rn_e(nr)+rn_e(nr-1)))
       rn_icl(nr)  = rn_id(nr) / (0.5d0*(rn_i(nr)+rn_i(nr-1)))

       ! safety factor
       qp_d(nr) = (qp(nr)-qp(nr-1)) / dr

       ! magnetic shear
       mshear(nr) = qp_d(nr) * rmnrho(nr)/(0.5d0*(qp(nr)+qp(nr-1)))

       ! magnetic curvature
       mcurv = 0.d0

       ! mean atomic mass of thermal ions [AMU] (on half grids)
       pansum = 0.d0
       DO nsa = 1, nsamax
          IF(idnsa(nsa)==1)THEN
             ns     = ns_nsa(nsa)
             pansum = pansum + pa(ns)*0.5d0*(rn(nsa,nr)+rn(nsa,nr-1))
          END IF
       END DO
       ai_ave(nr)   = pansum / (0.5d0*(rn_i(nr)+rn_i(nr-1)))

       ! MHD alpha
       alpha(nr) = -2.d0*rmu0*(0.5d0*(qp(nr)+qp(nr-1)))**2*rp_totd(nr) &
                   /((0.5d0*(abb1rho(nr)+abb1rho(nr-1)))**2            &
                    +(0.5d0*(bp(nr)+bp(nr-1)))**2)

       ! impurity and hydrogen density, density weighted charge/atomic number

    END DO

!    nrd1(1:nrmax) = rp_totd(1:nrmax)
!    nrd2(1:nrmax) = rt_id(1:nrmax)
!    nrd3(0:nrmax) = rmnrho(0:nrmax)

  END SUBROUTINE tr_calc_grad


  SUBROUTINE tr_calc_er
!---------------------------------------------------------------------------
!        Radial Electric Field (on grid)
!---------------------------------------------------------------------------
    USE trcomm, ONLY: nrmax,aee,bb,bp,rg,rt,rn,rg,rm,er, &
         pa,pz,mdler,vtor,vpol,rp_d

    IMPLICIT NONE
    INTEGER(ikind) :: nr
    REAl(rkind) :: term_dp

    DO nr = 1, nrmax

!    *****
!    ! This part will be implemented in introducing the NBI heating modules.
!    *****
!       DRL=RJCB(NR)/DR
!       IF(NR.EQ.NRMAX) THEN
!          ! second speices (main ion) only
!          DPD = DERIV3P(PNSS(2)*PTS(2), &
!                        RN(NR  ,2)*RT(NR  ,2)-PADD(NR  ), &
!                        RN(NR-1,2)*RT(NR-1,2)-PADD(NR-1), &
!                        RHOG(NR),RHOM(NR),RHOM(NR-1))
!          TERM_DP = DPD*RKEV/(PZ(2)*AEE*PNSS(2))
!       ELSE
       ! pressure gradient of bulk D(main ion)
!       DPD =(RN(NR+1,2)*RT(NR+1,2)-PADD(NR+1)-(RN(NR  ,2)*RT(NR  ,2)-PADD(NR  )))*DRL
!          TERM_DP = DPD*RKEV/(PZ(2)*AEE*0.5D0*(RN(NR+1,2)+RN(NR,2)))
!       ENDIF

! DPD -> rp_d(nsa,nr)

       ! second species (main ion: D(H)) only
       term_dp = rp_d(2,nr) / (pz(2)*aee*0.5d0*(rn(2,nr)+rn(2,nr-1))*1.d20)

       ! toroidal and poloidal rotation velocity <- from experiments for now
!       vtor = 0.d0
!       vpol = 0.d0

       IF(mdler.EQ.0) THEN
          ! pressure gradient only (nabla p)
          er(nr) = term_dp
       ELSE IF(mdler.EQ.1) THEN
          ! nabla p + toroidal rotation (V_tor)
          er(nr) = term_dp+vtor(nr)*bp(nr)
       ELSE IF(mdler.EQ.2) THEN
          ! nabla p + V_tor + poloidal rotation (V_pol) *** typical ER ***
          er(nr) = term_dp+vtor(nr)*bp(nr)-vpol(nr)*bb
!       ELSEIF(MDLER.EQ.3) THEN
!          !     Waltz definition
!          EPS = EPSRHO(NR)
!          F_UNTRAP = 1.D0-1.46D0*SQRT(EPS)+0.46D0*EPS**1.5D0
!          ALPHA_NEO = 1.D0-0.8839D0*F_UNTRAP/(0.3477D0+0.4058D0*F_UNTRAP)
!          IF(NR.EQ.NRMAX) THEN
!             TEL = PTS(1)
!             TIL = PTS(2)
!             RLNI = -DERIV3P(PNSS(2),RN(NR,2),RN(NR-1,2),RHOG(NR),RHOM(NR),R\
!             HOM(NR-1))/PNSS(2)
!             RLTI = -DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2),RHOG(NR),RHOM(NR),RH\
!             OM(NR-1))/PTS(2)
!          ELSE
!             TEL = 0.5D0*(RT(NR,1)+RT(NR+1,1))
!             TIL = 0.5D0*(RT(NR,2)+RT(NR+1,2))
!             RLNI = -(LOG(RN(NR+1,2))-LOG(RN(NR,2)))*DRL
!             RLTI = -(LOG(RT(NR+1,2))-LOG(RT(NR,2)))*DRL
!          ENDIF
!          CS = SQRT(TEL*RKEV/(PA(2)*AMM))
!          RHO_S = CS*PA(2)*AMM/(PZ(2)*AEE*BB)
!          ER(NR) =-BB*( (TIL/TEL)*RHO_S*CS*(RLNI+ALPHA_NEO*RLTI)-EPS/QP(NR)*\
!          VTOR(NR))
       ENDIF

    ENDDO

    RETURN
  END SUBROUTINE tr_calc_er


  SUBROUTINE tr_calc_exb
! --------------------------------------------------------------------------
!       E x B velocity and velocity shear
! --------------------------------------------------------------------------
    USE trcomm, ONLY: nrmax,RR,rhog,rhom,BB,dpdrho,bp,er, &
         vexbp,dvexbpdr,wexbp, &
         nrd1,nrd2,nrd3,nrd4

    IMPLICIT NONE
    INTEGER(ikind) :: nr
    REAL(rkind) :: bp_m
    REAL(rkind),DIMENSION(0:nrmax) :: dpvexbp,drhog

    vexbp(0:nrmax)    = 0.d0
    dvexbpdr(0:nrmax) = 0.d0
    wexbp(0:nrmax)    = 0.d0

    drhog(1:nrmax) = 0.5d0*(rhog(0:nrmax-1)-rhog(1:nrmax))

    ! on half grids

    ! ExB velocity (toroidal angular speed) [1/s]
    vexbp(1:nrmax) = -er(1:nrmax)/(RR*0.5d0*(bp(0:nrmax-1)+bp(1:nrmax)))

    dpvexbp(1:nrmax) = vexbp(1:nrmax) &
                      /(0.5d0*(dpdrho(0:nrmax-1)+dpdrho(1:nrmax)))

    DO nr = 1, nrmax
       bp_m  = 0.5d0*(bp(nr)+bp(nr-1))

       !:the 2nd order accuracy derivative of V_exb
!       dvexbpdr(nr) = deriv3(nr,rhom,dpvexbp,nrmax,1)
       IF(nr==1)THEN
          dvexbpdr(nr) = (dpvexbp(nr+1)-dpvexbp(nr)) &
                        /(0.5d0*(drhog(1)+drhog(2)))
       ELSE IF(nr == nrmax)THEN
          dvexbpdr(nr) = (dpvexbp(nr)-dpvexbp(nr-1)) &
                        /(0.5d0*(drhog(nrmax-1)+drhog(nrmax)))
       ELSE
          dvexbpdr(nr) = (dpvexbp(nr+1)-dpvexbp(nr-1)) &
                        /(0.5d0*(drhog(nr-1)*drhog(nr+1))+drhog(nr))
       END IF
!       dvexbpdr(nr) = 0.d0


!      ExB Rotation shear
!       "Effects of {ExB} velocity shear and magnetic shear
!           on turbulence and transport in magnetic confinement devices"
!       [Phys. of Plasmas, 4, 1499 (1997)]          
       wexbp(nr) = (RR*bp_m)**2/sqrt(BB**2+bp_m**2) &
                   *dvexbpdr(nr)
       ! wexbp(nr) = RR*bp(nr)*dvexbpdr(nr)/BB

    END DO

    nrd3(0:nrmax) = bp(0:nrmax)
    nrd4(1:nrmax) = dpvexbp(1:nrmax)

    RETURN
  END SUBROUTINE tr_calc_exb


  SUBROUTINE tr_calc_zeff
!---------------------------------------------------------------------------
!   Caluculate Z_eff (effective charge) on grids assuming quasi-neutrality
!    
!    Z_eff = (n_h + Z_imp^2 n_imp + Z_fi^2 n_fi) / n_e
!---------------------------------------------------------------------------
!!$    USE TRCOMM, ONLY : ANC, ANFE, MDLEQN, MDLUF, NRMAX, PZ, PZC, PZFE, RN, R\
!!$    NF, RT, ZEFF
    USE trcomm, ONLY: idnsa,ns_nsa,nrmax,pz,rn,mdluf,z_eff
    IMPLICIT NONE
    INTEGER(ikind) :: nr,nsa,ns

    z_eff(0:nrmax)   = 0.d0

    ! not include impurities
!    IF(mdluf == 0)THEN
       DO nr = 0, nrmax
          DO nsa = 1, nsamax
             IF(idnsa(nsa) == 1)THEN
                ns = ns_nsa(nsa)
                z_eff(nr) = z_eff(nr) + pz(ns)**2 *rn(nsa,nr)                
             END IF
          END DO
          z_eff(nr) = z_eff(nr)/rn(1,nr)
       END DO
!    END IF

!!$    IF(MDLUF.EQ.0) THEN
!!$       DO NR=1,NRMAX
!!$          TE =RT(NR,1)
!!$          ZEFF(NR) =(PZ(2)  *PZ(2)  *RN(NR,2) +PZ(3)  *PZ(3)  *RN(NR,3) +PZ(\
!!$4         )  *PZ(4)  *RN(NR,4) &
!!$               &     &                +TRZEC(TE)**2   *ANC (NR) +TRZEFE(TE)**2  *ANFE(NR))/RN(\
!!$          NR,1)
!!$       ENDDO
!!$    ELSE
!!$       IF(MDLEQN.EQ.0) THEN ! fixed density
!!$          DO NR=1,NRMAX
!!$             ZEFF(NR) =(PZ(2)  *PZ(2)  *RN(NR,2) +PZ(3)  *PZ(3)  *RN(NR,3) +\
!!$             PZ(2)  *PZ(2)  *RNF(NR,1))/RN(NR,1)
!!$          ENDDO
!!$       ENDIF
!!$    ENDIF
!!$
!!$    DO NR=1,NRMAX
!!$       TE=RT(NR,1)
!!$       PZC(NR)=TRZEC(TE)
!!$       PZFE(NR)=TRZEFE(TE)
!!$    ENDDO

    RETURN
  END SUBROUTINE tr_calc_zeff


  SUBROUTINE tr_calc_clseta
! ----------------------------------------------------------------------
!                classical resitivity
! ----------------------------------------------------------------------
    USE trcomm, ONLY: aee,ame,nrmax,rn,rt,eta,z_eff

    IMPLICIT NONE
    INTEGER(ikind) :: nr
    REAL(rkind) :: taue

      DO nr = 0, nrmax
!        ****** CLASSICAL RESISTIVITY (Spitzer) from JAERI Report ******

         ! electron collision time with ions
         taue = FTAUE(rn(1,nr),rn(2,nr),rt(1,nr),z_eff(nr))

         eta(nr) = ame/(rn(1,nr)*1.d20*aee**2*taue) &
                      *(0.29d0+0.46d0/(1.08d0+z_eff(nr)))

!         eta(nr) = 1.d-7
      END DO

    RETURN
  END SUBROUTINE tr_calc_clseta

! ----------------------------------------------------------------------
!           COULOMB LOGARITHM
! ----------------------------------------------------------------------
      REAL(8) FUNCTION COULOG(NS1,NS2,ANEL,TL)

!     ANEL : electron density [10^20 /m^3]
!     TL   : electron or ion temperature [keV]
!            in case of ion-ion collision, TL becomes ion temp.

      IMPLICIT NONE
      INTEGER(4):: NS1,NS2
      REAL(8)   :: ANEL,TL

      ! *** ref: 'Tokamaks 3rd Edition' p727 Coulomb logarithm ***
      IF(NS1.EQ.1.AND.NS2.EQ.1) THEN
         COULOG=14.9D0-0.5D0*LOG(ANEL)+LOG(TL)
      ELSE
         IF(NS1.EQ.1.OR.NS2.EQ.1) THEN
            COULOG=15.2D0-0.5D0*LOG(ANEL)+LOG(TL)
         ELSE
            COULOG=17.3D0-0.5D0*LOG(ANEL)+1.5D0*LOG(TL)
         ENDIF
      ENDIF

      RETURN
      END FUNCTION COULOG

! ----------------------------------------------------------------------
!           COLLISION TIME
! ----------------------------------------------------------------------
!     between electrons and ions
      REAL(8) FUNCTION FTAUE(ANEL,ANIL,TEL,ZL)

!     ANEL : electron density [10^20 /m^3]
!     ANIL : ion density [10^20 /m^3]
!     TEL  : electron temperature [kev]
!     ZL   : ion charge number

      USE TRCOMM, ONLY : AEE, AME, EPS0, PI, PZ, RKEV
      IMPLICIT NONE
      REAL(8) :: ANEL, ANIL, TEL, ZL
      REAL(8) :: COEF

      COEF = 6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)/(AEE**4*1.D20)
      IF(ZL-PZ(2).LE.1.D-7) THEN
         FTAUE = COEF*(TEL*RKEV)**1.5D0/(ANIL*ZL**2*COULOG(1,2,ANEL,TEL))
      ELSE
!     If the plasma contains impurities, we need to consider the
!     effective charge number instead of ion charge number.
!     From the definition of Zeff=sum(n_iZ_i^2)/n_e,
!     n_iZ_i^2 is replaced by n_eZ_eff at the denominator of tau_e.
         FTAUE = COEF*(TEL*RKEV)**1.5D0/(ANEL*ZL*COULOG(1,2,ANEL,TEL))
      ENDIF

      RETURN
      END FUNCTION FTAUE

!     between ions and ions
      REAL(8) FUNCTION FTAUI(ANEL,ANIL,TIL,ZL,PAL)

!     ANEL : electron density [10^20 /m^3]
!     ANIL : ion density [10^20 /m^3]
!     TIL  : ion temperature [kev]
!     ZL   : ion charge number
!     PAL  : ion atomic number

      USE TRCOMM, ONLY : AEE, AMP, EPS0, PI, RKEV
      IMPLICIT NONE
      REAL(8):: ANEL, ANIL, PAL, TIL, ZL
      REAL(8):: COEF

      COEF = 12.D0*PI*SQRT(PI)*EPS0**2*SQRT(PAL*AMP)/(AEE**4*1.D20)
      FTAUI = COEF*(TIL*RKEV)**1.5D0/(ANIL*ZL**4*COULOG(2,2,ANEL,TIL))

      RETURN
      END FUNCTION FTAUI


END MODULE trcalv

