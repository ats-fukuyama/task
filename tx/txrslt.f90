!     $Id$
module results
  use commons
  implicit none
  public

contains
!***********************************************************
!
!           CALCULATE GLOBAL QUANTITIES
!
!***********************************************************

  SUBROUTINE TXGLOB

    use physical_constants, only : AEE, PI, rMU0, rKeV, AME, EION
    use libraries, only : INTG_F, INTG_P, DERIVS

    INTEGER :: I, NS, NF, NR
    REAL(8) :: RKAP, FKAP, RNINT, RPINT, RPEINT, RPIINT, ANFINT, RWINT, POHINT, &
         &     PNBINT, PNFINT, PRFeINT, PRFiINT, PRFeTOT, PRFiTOT, &
         &     AJTINT, AOHINT, ANBINT, SNBINT, FACT, &
         &     BBL, SUMML, SUMPL, PNES, PAI
    REAL(8) :: PIEINT, SIEINT, PCXINT, SUMM, SUMP, SUML
    REAL(8) :: EpsL, FTL, DDX, RL31, RL32, DDD, dPTeV, dPTiV, dPPe, dPPi, &
         &     dPPV, ALFA
    REAL(8), DIMENSION(1:NRMAX) :: BP, BETA, BETAP, BETAL, BETAPL, BETAQ
    real(8), dimension(0:NRMAX) :: Betadef, dBetadr, PP, BthV2
    real(8) :: dBetaSUM, BPINT
    real(8) :: DERIV3

    !     Line Averaged Density and Temperature
    !     Core Density and Temperature
    !     Electoron and ion Work Quantities

    RKAP = 1.D0
    FKAP = 1.D0

    VOLAVN =  2.D0 / RB**2 * INTG_F( PZ*PNiV(0:NRMAX)+PZ*PNbV(0:NRMAX) &
         &                          +PZ*RATIO(0:NRMAX)*PNbrpV(0:NRMAX)-PNeV(0:NRMAX))

    RNINT = INTG_F(PNeV)
    RPEINT = INTG_F(PeV)
    ANSAV(1) = RNINT*2.D0*PI/(PI*RA*RA)
    ANS0(1)  = PNeV(0)
    IF(RNINT > 0.D0) THEN
       TSAV(1) = RPEINT/RNINT
    ELSE
       TSAV(1)=0.D0
    END IF
    TS0(1) = PTeV(0)
    WST(1) = RPEINT*1.5D0*2.D0*PI*RR*2.D0*PI*rKeV*1.D14

    RNINT = INTG_F(PNiV)
    RPIINT = INTG_F(PiV)
    ANSAV(2) = RNINT*2.D0*PI/(PI*RA*RA)
    ANS0(2)  = PNiV(0)
    IF(RNINT > 0.D0) THEN
       TSAV(2) = RPIINT/RNINT
    ELSE
       TSAV(2)=0.D0
    END IF
    TS0(2) = PTiV(0)
    WST(2) = RPIINT*1.5D0*2.D0*PI*RR*2.D0*PI*rKeV*1.D14

    ANFINT = INTG_F(PNbV)
    RWINT  = INTG_F(SNB)
    WFT(1) = 0.5D0*AMI*RWINT**2.D0*1.5D0*2.D0*PI*RR*2.D0*PI*RKAP*rKeV*1.D14
    ANFAV(1) = ANFINT*2.D0*PI/(PI*RA*RA)
    ANF0(1)  = PNbV(0)
    IF(ANFINT > 0.D0) THEN
       TFAV(1)  = RWINT/ANFINT
    ELSE
       TFAV(1)  = 0.D0
    END IF
    IF(ANF0(1) > 0.D0) THEN
       TF0(1)  = (9.D0*SNB(1)-SNB(2))/8.D0/ANF0(1) ! undertaking
    ELSE
       TF0(1)  = 0.D0
    END IF

    !     Input powers

    POHINT = INTG_F(POH)
    PNBINT = INTG_F(PNB)
    !!  RNFINT = INTG_F(PNF)
    PNFINT = 0.D0
    POHT = POHINT*2.D0*PI*RR*2.D0*PI*RKAP/1.D6
    PNBT = PNBINT*2.D0*PI*RR*2.D0*PI*RKAP/1.D6
    PNFT = PNFINT*2.D0*PI*RR*2.D0*PI*RKAP/1.D6

    PRFeINT = INTG_F(PRFe)
    PRFiINT = INTG_F(PRFi)
    PRFeTOT = PRFeINT*2.D0*PI*RR*2.D0*PI*RKAP/1.D6
    PRFiTOT = PRFiINT*2.D0*PI*RR*2.D0*PI*RKAP/1.D6
    PRFT    = PRFeTOT+PRFiTOT

    !      PFINT = INTG_F(PFIN)
    !      PFINT = PFINT*2.D0*PI*RR*2.D0*PI*RKAP/1.D6
    !      DO NS=1,NSM
    !        PFINT = INTG_F(PFCL(0:NRMAX,NS))
    !        PFCLT(NS) = PFINT*2.D0*PI*RR*2.D0*PI*RKAP/1.D6
    !      END DO

    !     Output powers

    SIE(0:NRMAX) = PNeV(0:NRMAX)*rNuION(0:NRMAX)*1.D20
    PIE(0:NRMAX) = SIE(0:NRMAX)*EION*AEE
    PCX(0:NRMAX) = 1.5D0*PNiV(0:NRMAX)*rNuiCX(0:NRMAX)*1.D20*PTiV(0:NRMAX)*rKeV

!    PRLINT=INTG_F(PRL)
    PIEINT=INTG_F(PIE)
    SIEINT=INTG_F(SIE)
    PCXINT=INTG_F(PCX)

!    PRLT = PRLINT*2.D0*PI*RR*2.D0*PI*RKAP/1.D6
    PCXT = PCXINT*2.D0*PI*RR*2.D0*PI*RKAP/1.D6
    PIET = PIEINT*2.D0*PI*RR*2.D0*PI*RKAP/1.D6

    !     Currents

    AJTINT=INTG_F(AJ  )
    AOHINT=INTG_F(AJOH)
    ANBINT=INTG_F(AJNB)

    AJT   = AJTINT*          2.D0*PI*RKAP/1.D6
    AJOHT = AOHINT*          2.D0*PI*RKAP/1.D6
    AJNBT = ANBINT*          2.D0*PI*RKAP/1.D6

    !     Bootstrap currents

    PP(0:NRMAX) = PeV(0:NRMAX) + PiV(0:NRMAX)
    DO NR = 0, NRMAX
       EpsL  = R(NR) / RR
       BBL = SQRT(BphV(NR)**2 + BthV(NR)**2)
       ! +++ Original model +++
       dPPV = DERIV3(NR,R,PP,NRMAX,NRM,0) * 1.D20 * rKeV
       ALFA = (rNuei1(NR)+rNueNC(NR))/rNuei3(NR)*(BthV(NR)/BphV(NR))**2 &
            & + 2.D0*rNuei2(NR)/rNuei3(NR)*BthV(NR)/BphV(NR)
       AJBS1(NR) = -1.D0 / (1.D0 + ALFA) * BthV(NR) / (BBL * BphV(NR)) * rNueNC(NR) / rNuei3(NR) * dPPV
       ! +++ Hirshman model +++
       dPTeV = DERIV3(NR,R,PTeV,NRMAX,NRM,0) * RA
       dPTiV = DERIV3(NR,R,PTiV,NRMAX,NRM,0) * RA
       dPPe  = DERIV3(NR,R,PeV,NRMAX,NRM,0) * RA
       dPPi  = DERIV3(NR,R,PiV,NRMAX,NRM,0) * RA
       FTL   =(1.46D0 * SQRT(EpsL) + 2.4D0 * EpsL) / (1.D0 - EpsL)**1.5D0
       DDX   = 1.414D0 * PZ + PZ**2 + FTL * (0.754D0 + 2.657D0 * PZ &
            &        + 2.D0 * PZ**2) + FTL**2 * (0.348D0 + 1.243D0 * PZ + PZ**2)
       RL31  = FTL * ( 0.754D0 + 2.21D0 * PZ + PZ**2 + FTL * (0.348D0 + 1.243D0 &
            &            * PZ + PZ**2)) / DDX
       RL32  =-FTL * (0.884D0 + 2.074D0 * PZ) / DDX
       DDD   =-1.172D0 / (1.D0 + 0.462D0 * FTL)
       IF(NR == 0) THEN
          AJBS4(NR) = 0.D0
       ELSE
          AJBS4(NR) =- (PeV(NR) * 1.D20 * rKeV) &
               &    * (RL31 * ((dPPe / PeV(NR)) + (PTiV(NR) / (PZ * PTeV(NR))) &
               &    * ((dPPi / PiV(NR)) + DDD * (dPTiV / PTiV(NR)))) &
               &    + RL32 * (dPTeV / PTeV(NR))) / (RA * BthV(NR))
       END IF
    END DO

    AJBST = INTG_F(AJBS1) * 2.D0 * PI * RKAP * 1.D-6

    !      DRH=0.5D0*DR
    !      DO NS=1,NSM
    !         VNP=AV(NRMAX,NS)
    !         DNP=AD(NRMAX,NS)
    !
    !         VTP=(AVK(NRMAX,NS)/1.5D0+AV(NRMAX,NS))
    !         DTP= AK(NRMAX,NS)/(1.5D0)
    !
    !         VXP=0.D0
    !         DXP=(1.5D0*AD(NRMAX,NS)-AK(NRMAX,NS))*PTS(NS)
    !
    !         SLT(NS) =((     DNP/DRH)*RN(NRMAX,NS)
    !     &            +( VNP-DNP/DRH)*PNSS(NS))
    !     &            *2.D0*PI*RR*2.D0*PI*RA*FKAP
    !
    !         PLT(NS) =((     DXP/DRH)*RN(NRMAX,NS)
    !     &            +(     DTP/DRH)*RN(NRMAX,NS)*RT(NRMAX,NS)*1.5D0
    !     &            +( VXP-DXP/DRH)*PNSS(NS)
    !     &            +( VTP-DTP/DRH)*PNSS(NS)*PTS(NS)*1.5D0)
    !     &            *2.D0*PI*RR*2.D0*PI*RA*FKAP*rKeV*1.D14
    !      END DO
    !
    !!      CALL TXSUMD(SIE,R,NRMAX,SIEINT)
    !!      CALL TXSUMD(SNF,R,NRMAX,SNFINT)
    SNBINT = INTG_F(SNB)

    !      SIET = SIEINT*2.D0*PI*RR*2.D0*PI*DR*RKAP
    !      SNFT = SNFINT*2.D0*PI*RR*2.D0*PI*DR*RKAP
    SNBT = SNBINT*2.D0*PI*RR*2.D0*PI*RKAP

    !      DO NS=1,NSM
    !         SPEINT = INTG_F(SPE(0:NRMAX,NS))
    !         SPET(NS) = SPEINT*2.D0*PI*RR*2.D0*PI*RKAP
    !      END DO

    WBULKT = SUM(WST(1:2))
    WTAILT = SUM(WFT(1:1))

    WPT  = WBULKT+WTAILT
    PINT = PNBT+PRFT+POHT
    POUT = PCXT+PIET
    SINT = SIET+SNBT
!    SOUT = SLST

    IF(ABS(T_TX - TPRE) <= 1.D-70) THEN
       WPDOT = 0.D0
    ELSE
       WPDOT = (WPT - WPPRE)/(T_TX - TPRE)
    END IF
    WPPRE = WPT
    TPRE  = T_TX

    IF(PINT <= 0.D0) THEN
       TAUE1 = 0.D0
       TAUE2 = 0.D0
    ELSE
       TAUE1 = WPT/PINT
       TAUE2 = WPT/(PINT - WPDOT)
    END IF

!    FACT = 1.D0
    FACT = RKAP/FKAP**2

    Betadef(0:NRMAX) = (PNeV(0:NRMAX) * PTeV(0:NRMAX) + PNiV(0:NRMAX) * PTiV(0:NRMAX)) &
         &        * 1.D20 * rKeV /((BphV(0:NRMAX)**2 + BthV(0:NRMAX)**2) / (2.D0 * rMU0))
    CALL DERIVS(R,Betadef,NRMAX,dBetadr)
    RPEINT = 0.D0 ; RPIINT = 0.D0 ; dBetaSUM = 0.D0
    DO NR = 1, NRMAX
       RPEINT = RPEINT + INTG_P(PeV,NR,0)
       RPIINT = RPIINT + INTG_P(PiV,NR,0)
       dBetaSUM = dBetaSUM + INTG_P(dBetadr,NR,0)
       RPINT =(RPEINT + RPIINT)*rKeV*1.D20
       SUMM  = 2.D0*PI*RPINT
       SUMP  = PI*R(NR)**2*dBetaSUM
       BBL   = BthV(NR)**2 + BphV(NR)**2
       BETA  (NR) = 2.D0*SUMM *rMU0/(PI*R(NR)**2*BBL)
       BETAP (NR) = 2.D0*SUMM *rMU0/(PI*R(NR)**2*BthV(NRMAX)**2)*FACT
       BETAL (NR) = 2.D0*RPINT*rMU0/(            BBL)
       BETAPL(NR) = 2.D0*RPINT*rMU0/(            BthV(NRMAX)**2)*FACT
       BETAQ (NR) =-2.D0*SUMP *rMU0/(PI*R(NR)**2*BthV(NR)   **2)*FACT
    END DO

!!$    SUMM = 0.D0
!!$    SUMP = 0.D0
!!$    SUML = 0.D0
!!$    BP(1:NRMAX) = BthV(1:NRMAX)
!!$    DO I=1,NRMAX-1
!!$       BBL   = SQRT(BthV(I)**2+BphV(I)**2)
!!$       SUMML = PNeV(I)*PTeV(I)*rKeV*1.D20+PNiV(I)*PTiV(I)*rKeV*1.D20
!!$       SUMPL = (PNeV(I+1)*PTeV(I+1)-PNeV(I  )*PTeV(I  ) &
!!$            & +PNiV(I+1)*PTiV(I+1)-PNiV(I  )*PTiV(I  )) &
!!$            & *rKeV*1.D20/DR
!!$
!!$       SUMM = SUMM + SUMML*2.D0*PI*R(I)*DR
!!$       SUMP = SUMP + 0.5D0*SUMPL*PI*R(I)*R(I)*DR
!!$       SUML = SUML + BP(I)**2*R(I)
!!$       BETA(I)   = 2.D0*SUMM *rMU0/(PI*(R(I)*BBL)**2)
!!$       BETAL(I)  = 2.D0*SUMML*rMU0/(           BBL **2)
!!$       BETAP(I)  = 2.D0*SUMM *rMU0/(PI*(R(I)*BP(NRMAX))**2)*FACT
!!$       BETAPL(I) = 2.D0*SUMML*rMU0/(           BP(NRMAX)**2)*FACT
!!$       BETAQ(I)  =-2.D0*SUMP *rMU0/(PI*(R(I)*BP(I))**2)*FACT
!!$       SUMP = SUMP + 0.5D0*SUMPL*PI*R(I+1)*R(I+1)*DR
!!$    END DO
!!$
!!$    I=NRMAX
!!$    BBL  = SQRT(BthV(I)**2+BphV(I)**2)
!!$    SUMML =(PNeV(I)*PTeV(I)+PNiV(I)*PTiV(I))*R(I)*rKeV*1.D20
!!$    PNES = PNa * EXP(-(RB-RA) / rLn)
!!$    SUMPL = (PNES*PTea-PNeV(I)*PTeV(I)+PNES/PZ*PTia-PNiV(I)*PTiV(I)) &
!!$         & *rKeV*1.D20/DR*2.D0
!!$
!!$    !  DO NF=1,NFM
!!$    !     SUMML = SUMML +SNB(I)*R(I)*rKeV*1.D20
!!$    !     SUMPL = SUMPL +(0.D0-SNB(I))*rKeV*1.D20/DR*2.D0
!!$    !  END DO
!!$
!!$    SUMM = SUMM + SUMML*2.D0*PI*R(I)*DR
!!$    SUMP = SUMP + 0.5D0*SUMPL*PI*R(I)*R(I)*DR
!!$    SUML = SUML + 0.5D0*BP(I)**2*R(I)
!!$    BETA(I)   = 2.D0*SUMM *rMU0/(PI*(R(I)*BBL)**2)
!!$    BETAL(I)  = 2.D0*SUMML*rMU0/(           BBL **2)
!!$    BETAP(I)  = 2.D0*SUMM *rMU0/(PI*(R(I)*BP(NRMAX))**2)*FACT
!!$    BETAPL(I) = 2.D0*SUMML*rMU0/(           BP(NRMAX)**2)*FACT
!!$    BETAQ(I)  =-2.D0*SUMP *rMU0/(PI*(R(I)*BP(I))**2)*FACT

    BETA0  = (4.D0*BETA(1)  -BETA(2)  )/3.D0
    BETAP0 = (4.D0*BETAPL(1)-BETAPL(2))/3.D0
    BETAQ0 = 0.D0

    BETAPA = BETAP(NRMAX)
    BETAA  = BETA(NRMAX)

    BETAN  = BETAA / (rIP / (RA * BB)) * 100.D0

    BthV2(0:NRMAX) = BthV(0:NRMAX)**2
    BPINT = 0.5D0 * INTG_F(BthV2)
    ALI   = 8.D0*PI**2*BPINT*FKAP**2/((rMU0*rIp*1.D6)**2)
    VLOOP = EphV(NRMAX)*2.D0*PI*RR

    !  PAI=(PA(2)*PN(2)+PA(3)*PN(3)+PA(4)*PN(4))/(PN(2)+PN(3)+PN(4))
    PAI=PA

    IF(PINT <= 0.D0) THEN
       TAUEP=0.D0
       TAUEH=0.D0
       QF=0.D0
    ELSE
       TAUEP=4.8D-2*(rIP**0.85D0) &
            &      *(RR**1.2D0) &
            &      *(RA**0.3D0) &
            &      *(RKAP**0.5D0) &
            &      *(ANSAV(1)**0.1D0) &
            &      *(BB**0.2D0) &
            &      *(PAI**0.5D0) &
            &      *(PINT**(-0.5D0))
       TAUEH=0.0562D0*(rIP**0.93D0) &
            &        *(BB**0.15D0) &
            &        *(PINT**(-0.69D0)) &
            &        *(ANSAV(1)**0.41D0) &
            &        *(PAI**0.19D0) &
            &        *(RR**1.97D0) &
            &        *((RA/RR)**0.58D0) &
            &        *(RKAP**0.78D0)
       QF=5.D0*PNFT/PINT
    END IF

    IF(Q(0) >= 1.D0) THEN
       RQ1=0.D0
    ELSE
       OUTER:DO 
          IF(Q(1) > 1.D0) THEN
             RQ1=SQRT( (1.D0-Q(0) )/(Q(1)-Q(0)) )*H(1)
             EXIT OUTER
          END IF
          DO I=2,NRMAX
             IF(Q(I) > 1.D0) THEN
                RQ1=(R(I)-R(I-1))*(1.D0-Q(I-1))/(Q(I)-Q(I-1))+R(I-1)
                EXIT OUTER
             END IF
          END DO
          RQ1=RA
          EXIT OUTER
       END DO OUTER
    END IF

    !  ZEFF0=(4.D0*ZEFF(1)-ZEFF(2))/3.D0
    ZEFF0=Zeff

    !  Evaluate effective particle diffusivity from particle flux minus Ware pinch
    !  Gamma = n v - Deff dn/dr  --> Deff = (n v - ft n Eph / Bth) / (- dn/dr)

    Deff(0) = 0.d0
    do nr = 1, nrmax
       Deff(nr) = (PNiV(NR)*UirV(NR)-ft(NR)*PNiV(NR)*EphV(NR)/BthV(NR)) &
               & /(-DERIV3(NR,R,PNiV,NRMAX,NRM,0))
       if(Deff(nr) < 0.d0) Deff(nr) = 0.d0
    end do

    RETURN
  END SUBROUTINE TXGLOB

  subroutine TXSTAT

    write(6,'(1X,2(A26,1PD10.3,3X))') "Vol. ave. of neutrality = ", VOLAVN
    write(6,'(1X,2(A26,1PD10.3,3X))') "Inductance              = ", ALI, &
         &                            "Loop voltage            = ", VLOOP
    write(6,'(1X,2(A26,1PD10.3,3X))') "Confinement time 1      = ", TAUE1, &
         &                            "Confinement time 2      = ", TAUE2
    write(6,'(1X,2(A26,1PD10.3,3X))') "L-mode scaling time     = ", TAUEP, &
         &                            "IPB98(y,2) scaling time = ", TAUEH
    write(6,'(1X,2(A26,1PD10.3,3X))') "Beta                    = ", BETAA, &
         &                            "Poloidal beta           = ", BETAPA
    write(6,'(1X,2(A26,1PD10.3,3X))') "Normalized beta         = ", BETAN

  end subroutine TXSTAT

end module results
