!     $Id$

!***********************************************************
!
!           CALCULATE GLOBAL QUANTITIES
!
!***********************************************************

SUBROUTINE TXGLOB

  USE physical_constants, only : AEE, PI, rMU0, rKEV
  INCLUDE 'txcomm.inc'

  INTEGER :: I, NS, NF
  REAL(8) :: RKAP, FKAP, RNSUM, RTSUM, ANFSUM, RWSUM, POHSUM, &
       &     PNBSUM, PNFSUM, PRFeSUM, PRFiSUM, PRFeTOT, PRFiTOT,  &
       &     EION, PIE, SIE, TNU, PCX, &
       &     AJTSUM, AOHSUM, ANBSUM, SNBSUM,  FACT, &
       &     BBL, SUMML, SUMPL, PNES, PAI
  REAL(8) :: PIESUM, SIESUM, PCXSUM, SUMM, SUMP, SUML
  REAL(8), DIMENSION(NRM) :: BP, BETA, BETAP
  REAL(8), DIMENSION(NRM) :: BETAL, BETAPL, BETAQ

!     Line Averaged Density and Temperature
!     Core Density and Temperature
!     Electoron and ion Work Quantities

  RKAP=1.D0
  FKAP=1.D0

  RNSUM=SUM(PNeHI(0:NRMAX)*RHI(0:NRMAX))
  RTSUM=SUM(PTeHI(0:NRMAX)*PNeHI(0:NRMAX)*RHI(0:NRMAX))
  ANSAV(1) = RNSUM*2.D0*PI*DR/(PI*RA*RA)
  ANS0(1)  = PNeI(0)
  IF(RNSUM > 0.D0) THEN
     TSAV(1) = RTSUM/RNSUM
  ELSE
     TSAV(1)=0.D0
  END IF
  TS0(1) = PTeI(0)
  WST(1) = RTSUM*1.5D0*2.D0*PI*RR*2.D0*PI*DR*RKEV*1.D14

  RNSUM=SUM(PNiHI(0:NRMAX)*RHI(0:NRMAX))
  RTSUM=SUM(PTiHI(0:NRMAX)*PNeHI(0:NRMAX)*RHI(0:NRMAX))
  ANSAV(2) = RNSUM*2.D0*PI*DR/(PI*RA*RA)
  ANS0(2)  = PNiI(0)
  IF(RNSUM > 0.D0) THEN
     TSAV(2) = RTSUM/RNSUM
  ELSE
     TSAV(2)=0.D0
  END IF
  TS0(2) = PTiI(0)
  WST(2) = RTSUM*1.5D0*2.D0*PI*RR*2.D0*PI*DR*RKEV*1.D14

  ANFSUM=SUM(PNbHI(0:NRMAX)*RHI(0:NRMAX))
  RWSUM =SUM(SNB(0:NRMAX)*RHI(0:NRMAX))
  WFT(1) = 0.5D0*AMi*RWSUM**2.D0*1.5D0 &
       &           *2.D0*PI*RR*2.D0*PI*DR*RKAP*RKEV*1.D14
  ANFAV(1) = ANFSUM*2.D0*PI*DR/(PI*RA*RA)
  ANF0(1)  = PNbI(0)
  IF(ANFSUM > 0.D0) THEN
     TFAV(1)  = RWSUM/ANFSUM
  ELSE
     TFAV(1)  = 0.D0
  END IF
  IF(ANF0(1) > 0.D0) THEN
     TF0(1)  = (9.D0*SNB(1)-SNB(2))/8.D0/ANF0(1)
  ELSE
     TF0(1)  = 0.D0
  END IF

!     Input powers

  POHSUM=SUM(POH(0:NRMAX)*RHI(0:NRMAX))
  PNBSUM=SUM(PNB(0:NRMAX)*RHI(0:NRMAX))
!!  RNFSUM=SUM(PNF(0:NRMAX)*RHI(0:NRMAX))
  PNFSUM=0.D0
  POHT = POHSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
  PNBT = PNBSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
  PNFT = PNFSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6

  PRFeSUM=SUM(PRFe(0:NRMAX)*RHI(0:NRMAX))
  PRFiSUM=SUM(PRFi(0:NRMAX)*RHI(0:NRMAX))
  PRFeTOT=PRFeSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
  PRFiTOT=PRFiSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
  PRFT=PRFeTOT+PRFiTOT

!      DO NS=1,NSM
!        PRFSUM=SUM(PRF(0:NRMAX)*RHI(0:NRMAX))
!        PRFT = PRFSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
!      END DO

!      PBSUM=SUM(PBIN(0:NRMAX)*RHI(0:NRMAX))
!      PBINT = PBSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
!      DO NS=1,NSM
!        CALL TXSUMD(PBCL,RHI,NRMAX,PBSUM)
!        PBCLT = PBSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
!      END DO

!      PFSUM=SUM(PFIN(0:NRMAX)*RHI(0:NRMAX))
!      PFINT = PFSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
!      DO NS=1,NSM
!        CALL TXSUMD(PFCL(1,NS),RHI,NRMAX,PFSUM)
!        PFCLT(NS) = PFSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
!      END DO

!      PRLSUM=SUM(PRL(0:NRMAX)*RHI(0:NRMAX))

!     Output powers

  PIESUM=0.D0
  SIESUM=0.D0
  PCXSUM=0.D0

  DO I=1,NRMAX
     EION  = 13.64D0
     PIE=PNeHI(I)*rNuION(I)*1.D20*EION*AEE
     SIE=PNeHI(I)*rNuION(I)*1.D20
     PIESUM=PIESUM+PIE*RHI(I)
     SIESUM=SIESUM+SIE*RHI(I)

     TNU=0.D0
     PCX=(-1.5D0*PNeHI(I)*rNuION(I)*TNU/1.D20 &
         &+1.5D0*PNiHI(I)*rNuiCX(I)*(PTiHI(I)-TNU))*RKEV*1.D20
     PCXSUM=PCXSUM+PCX*RHI(I)
  END DO

!      PRLT = PRLSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
  PCXT = PCXSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
  PIET = PIESUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6

!     Currents

  AJTSUM=SUM(AJ  (0:NRMAX)*RHI(0:NRMAX))
  AOHSUM=SUM(AJOH(0:NRMAX)*RHI(0:NRMAX))
  ANBSUM=SUM(AJNB(0:NRMAX)*RHI(0:NRMAX))

  AJT   = AJTSUM*            2.D0*PI*DR*RKAP/1.D6
  AJOHT = AOHSUM*            2.D0*PI*DR*RKAP/1.D6
  AJNBT = ANBSUM*            2.D0*PI*DR*RKAP/1.D6

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
!     &            *2.D0*PI*RR*2.D0*PI*RA*FKAP*RKEV*1.D14
!      END DO
!
!!      CALL TXSUMD(SIE,RHI,NRMAX,SIESUM)
!!      CALL TXSUMD(SNF,RHI,NRMAX,SNFSUM)
  SNBSUM=SUM(SNB(0:NRMAX)*RHI(0:NRMAX))

!      SIET = SIESUM*2.D0*PI*RR*2.D0*PI*DR*RKAP
!      SNFT = SNFSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP
  SNBT = SNBSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP

!      DO NS=1,NSM
!         CALL TXSUMD(SPE(1,NS),RHI,NRMAX,SPESUM)
!         SPET(NS) = SPESUM*2.D0*PI*RR*2.D0*PI*DR
!      END DO

  WBULKT=SUM(WST(1:2))
  WTAILT=SUM(WFT(1:1))

  WPT =WBULKT+WTAILT
  PINT=PNBT+PRFT+POHT
  POUT=PCXT+PIET
  SINT=SIET+SNBT
!      SOUT=SLST

  IF(ABS(T_TX-TPRE) <= 1.D-70) THEN
     WPDOT=0.D0
  ELSE
     WPDOT=(WPT-WPPRE)/(T_TX-TPRE)
  END IF
  WPPRE=WPT
  TPRE=T_TX

  IF(PINT <= 0.D0) THEN
     TAUE1=0.D0
     TAUE2=0.D0
  ELSE
     TAUE1=WPT/PINT
     TAUE2=WPT/(PINT-WPDOT)
  END IF

  SUMM=0.D0
  SUMP=0.D0
  SUML=0.D0
  FACT=1.D0
  FACT=RKAP/FKAP**2

  BP(1:NRMAX)=BthI(1:NRMAX)

  DO I=1,NRMAX-1
     BBL   = SQRT(BthI(I)**2+BphI(I)**2)
     SUMML = PNeHI(I)*PTeHI(I)*RKEV*1.D20+PNiHI(I)*PTiHI(I)*RKEV*1.D20
     SUMPL = (PNeHI(I+1)*PTeHI(I+1)-PNeHI(I  )*PTeHI(I  ) &
          & +PNiHI(I+1)*PTiHI(I+1)-PNiHI(I  )*PTiHI(I  )) &
          & *RKEV*1.D20/DR

     SUMM = SUMM + SUMML*2.D0*PI*RHI(I)*DR
     SUMP = SUMP + 0.5D0*SUMPL*PI*RHI(I)*RHI(I)*DR
     SUML = SUML + BP(I)**2*R(I)
     BETA(I)   = 2.D0*SUMM *rMU0/(PI*(R(I)*BBL)**2)
     BETAL(I)  = 2.D0*SUMML*rMU0/(           BBL **2)
     BETAP(I)  = 2.D0*SUMM *rMU0/(PI*(R(I)*BP(NRMAX))**2)*FACT
     BETAPL(I) = 2.D0*SUMML*rMU0/(           BP(NRMAX)**2)*FACT
     BETAQ(I)  =-2.D0*SUMP *rMU0/(PI*(R(I)*BP(I))**2)*FACT
     SUMP = SUMP + 0.5D0*SUMPL*PI*RHI(I+1)*RHI(I+1)*DR
  END DO

  I=NRMAX
  BBL  = SQRT(BthI(I)**2+BphI(I)**2)
  SUMML =(PNeHI(I)*PTeHI(I)+PNiHI(I)*PTiHI(I))*RHI(I)*RKEV*1.D20
  PNES = PNa * EXP(-(RB-RA) / rLn)
  SUMPL = (PNES*PTea-PNeHI(I)*PTeHI(I)+PNES/PZ*PTia-PNiHI(I)*PTiHI(I)) &
       & *RKEV*1.D20/DR*2.D0

!  DO NF=1,NFM
!     SUMML = SUMML +SNB(I)*RHI(I)*RKEV*1.D20
!     SUMPL = SUMPL +(0.D0-SNB(I))*RKEV*1.D20/DR*2.D0
!  END DO

  SUMM = SUMM + SUMML*2.D0*PI*RHI(I)*DR
  SUMP = SUMP + 0.5D0*SUMPL*PI*RHI(I)*RHI(I)*DR
  SUML = SUML + 0.5D0*BP(I)**2*R(I)
  BETA(I)   = 2.D0*SUMM *rMU0/(PI*(R(I)*BBL)**2)
  BETAL(I)  = 2.D0*SUMML*rMU0/(           BBL **2)
  BETAP(I)  = 2.D0*SUMM *rMU0/(PI*(R(I)*BP(NRMAX))**2)*FACT
  BETAPL(I) = 2.D0*SUMML*rMU0/(           BP(NRMAX)**2)*FACT
  BETAQ(I)  =-2.D0*SUMP *rMU0/(PI*(R(I)*BP(I))**2)*FACT

  BETA0 =(4.D0*BETA(1)  -BETA(2)  )/3.D0
  BETAP0=(4.D0*BETAPL(1)-BETAPL(2))/3.D0
  BETAQ0=0.D0

  BETAPA=BETAP(NRMAX)
  BETAA =BETA(NRMAX)

  ALI=8.D0*PI**2*DR*SUML*FKAP**2/((rMU0*rIp*1.D6)**2)
  VLOOP = EphHI(NRMAX-1)*2*PI*RR

!  PAI=(PA(2)*PN(2)+PA(3)*PN(3)+PA(4)*PN(4))/(PN(2)+PN(3)+PN(4))
  PAI=PA

  IF(PINT <= 0.D0) THEN
     TAUEP=0.D0
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
     QF=5.D0*PNFT/PINT
  END IF

  IF(Q(0) >= 1.D0) THEN
     RQ1=0.D0
  ELSE
     IF(Q(1) > 1.D0) THEN
        RQ1=SQRT( (1.D0-Q(0) )/(Q(1)-Q(0)) )*DR
        GOTO 310
     END IF
     DO I=2,NRMAX
        IF(Q(I) > 1.D0) THEN
           RQ1=(R(I)-R(I-1))*(1.D0-Q(I-1))/(Q(I)-Q(I-1))+R(I-1)
           GOTO 310
        END IF
     END DO
     RQ1=RA
310  CONTINUE
  END IF

!  ZEFF0=(4.D0*ZEFF(1)-ZEFF(2))/3.D0
  ZEFF0=Zeff

  RETURN
END SUBROUTINE TXGLOB

!********************************************************
!
!      RADIAL INTEGRATION  (DOUBLE VARIABLES)
!
!********************************************************

!!$      SUBROUTINE TXSUMD(A,B,NMAX,SUM)
!!$
!!$      IMPLICIT NONE
!!$      INTEGER, INTENT(IN) :: NMAX
!!$      REAL(8), DIMENSION(NMAX), INTENT(IN) :: A, B
!!$      REAL(8), INTENT(OUT) :: SUM
!!$      INTEGER :: N
!!$
!!$      SUM = 0.D0
!!$      DO N=1,NMAX
!!$         SUM=SUM+A(N)*B(N)
!!$      END DO
!!$
!!$      RETURN
!!$      END
