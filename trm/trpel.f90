!     ***********************************************************

!           PELLET INJECTION

!     ***********************************************************

      SUBROUTINE TRPELT

      USE TRCOMM, ONLY : MDLPEL, PELTOT
      IMPLICIT NONE


      IF(PELTOT.LE.0.D0) RETURN

      IF(MDLPEL.EQ.0) RETURN

      WRITE(6,*)'# PELLET INJECTED.'

      IF(MDLPEL.EQ.1) THEN
         CALL TRPELA
      ELSEIF(MDLPEL.EQ.2) THEN
         CALL TRPELB
      ELSEIF(MDLPEL.EQ.3) THEN
         CALL TRPELC
      ENDIF

      RETURN
      END  SUBROUTINE TRPELT

!     ***********************************************************

!           PELLET INJECTION (GAUSSIAN PROFILE)

!     ***********************************************************

      SUBROUTINE TRPELA

      USE TRCOMM, ONLY : DR, DVRHO, NRMAX, NSM, PELPAT, PELR0, PELRW,&
           PELTOT, RA, RM, SPE, rkind
      IMPLICIT NONE
      REAL(rkind)    :: FSUM, S0, SPEL
      INTEGER :: NR, NS


      FSUM = 0.D0
      DO NR=1,NRMAX
         FSUM=FSUM+DEXP(-((RA*RM(NR)-PELR0)/PELRW)**2)*DVRHO(NR)*DR
      ENDDO

      S0  =PELTOT/FSUM
      WRITE(6,*) 'S0=',S0

      DO NR=1,NRMAX
         SPEL=S0*DEXP(-((RA*RM(NR)-PELR0)/PELRW)**2)
      DO NS=1,NSM
         SPE(NR,NS)=PELPAT(NS)*SPEL
      ENDDO
      ENDDO

      WRITE(6,*) 'TRPELA:'
      WRITE(6,'(5ES12.4)') SPE(1:NRMAX,1)

      RETURN
      END  SUBROUTINE TRPELA

!     ***********************************************************

!           PELLET ABLATION MODEL (WAKATANI NAKAMURA MODEL)

!     ***********************************************************

      SUBROUTINE TRPELB

      USE TRCOMM, ONLY : AEE, AME, AMP, ANC, ANFE, DR, DVRHO, NRMAX,&
           NSM, PA, PELPAT, PELRAD, PELVEL, PI, PNBENG, PZ, PZC,&
           PZFE, RA, RKEV, RN, RPE, RT, SNB, SNF, SPE, rkind, &
           NNBMAX,SNB_NNB
      IMPLICIT NONE
      REAL(rkind)    :: A1, A2, A3, AMA, AMB, AMD, AMPEL, AMT, ANE,&
           & ANFAST, ANS, B1, B2, B3, EFAST, P1, P2, PAP, QFAST, ROS,&
           & RP, RPDOT, RPPRE, SPEL, TAUS, TE, VCR, VF, SUM
      REAL(rkind),ALLOCATABLE:: VB_NNB(:)
      INTEGER :: NR, NS, NNB
      REAL(rkind)    :: COULOG   ! FUNCTION

      IF(NNBMAX.GT.0) ALLOCATE(VB_NNB(NNBMAX))

      NR = NRMAX
      RP = PELRAD
      PAP= PA(2)*PELPAT(2)+PA(3)*PELPAT(3)+PA(4)*PELPAT(4)
      ROS= 89.D0*PAP
      AMPEL= AMP*PAP

      AMD = PA(2)*AMP
      AMT = PA(3)*AMP
      AMA = PA(4)*AMP
      AMB = PA(2)*AMP
      ANS = ROS/AMPEL*1.D-20
      VF  = SQRT(2.D0*3.5D3*RKEV/AMA)
      DO NNB=1,NNBMAX
         VB_NNB(NNB)  = SQRT(2.D0*PNBENG(NNB)*RKEV/AMB)
      END DO

 1000 ANE=RN(NR,1)
      TE=RT(NR,1)
      P1   = 3.D0*SQRT(0.5D0*PI)*AME/ANE *(ABS(TE)*RKEV/AME)**1.5D0

      VCR   = (P1*(RN(NR,2)*PZ(2)   **2/AMD + RN(NR,3)*PZ(3)   **2/AMT &
     &           + RN(NR,4)*PZ(4)   **2/AMA + ANFE(NR)*PZFE(NR)**2/52.D0 &
     &           + ANC (NR)*PZC(NR) **2/12.D0))**(1.D0/3.D0)

      TAUS = 0.2D0*PA(2)*ABS(TE)**1.5D0 /(PZ(2)**2*ANE*COULOG(1,2,ANE,TE))

      ANFAST=SNF(NR)*LOG(1.D0+(VF/VCR)**3)
      DO NNB=1,NNBMAX
         ANFAST=ANFAST+SNB_NNB(NNB,NR)*LOG(1.D0+(VB_NNB(NNB)/VCR)**3)
      END DO
      ANFAST=ANFAST*TAUS/3.D0
      
      B1   = SNF(NR)*(0.5D0*AMA*VCR*VCR/AEE)**3
      B2   = 0.5D0*(VF/VCR)**6-(VF/VCR)**3+LOG(1.D0+(VF/VCR)**3)
      SUM=0.D0
      DO NNB=1,NNBMAX
         A1   = SNB_NNB(NNB,NR)*(0.5D0*AMB*VCR*VCR/AEE)**3
         A2   = 0.5D0*(VB_NNB(NNB)/VCR)**6 &
              -(VB_NNB(NNB)/VCR)**3+LOG(1.D0+(VB_NNB(NNB)/VCR)**3)
         SUM=SUM+A1*A2
      END DO
      P2   = (B1*B2+SUM)*TAUS/3.D0

      EFAST= (P2/ANFAST)**(1.D0/3.D0)

      A3=0.D0
      DO NNB=1,NNBMAX
         A3=A3+SNB_NNB(NNB,NR)*1.D20*AMB &
              *((VB_NNB(NNB)/VCR)**3-LOG(1.D0+(VB_NNB(NNB)/VCR)**3))
      END DO
      B3   = SNF(NR)*1.D20*AMA*((VF/VCR)**3-LOG(1.D0+(VF/VCR)**3))

      QFAST= (A3+B3)*TAUS*VCR**3/(AEE*24.D0)

      RPDOT= 4.55D6*(2.D0*AMPEL)**(2.D0/3.D0)*QFAST**(1.D0/3.D0) &
     &       *EFAST**0.4D0*RP**(-2.D0/3.D0)/ROS

      RPPRE=RP
      RP=RP-RPDOT*DR*RA/PELVEL

      IF(RP.LE.0.D0) THEN
         RPE=RA-RA*DR*(NR-1)+RA*DR*RPPRE/(RPPRE-RP)
         GOTO 9000
      ENDIF

      IF(NR.EQ.1) THEN
         RPE=RA
         GOTO 9000
      ENDIF

      SPEL=ANS*RP*0.5D0*(RP+RPPRE)*RPDOT*4.D0*PI*RA/(DVRHO(NR)*PELVEL)
      DO NS=1,NSM
         SPE(NR,NS)=PELPAT(NS)*SPEL
      ENDDO

      NR = NR-1
      GOTO 1000

 9000 WRITE(6,600) RPE
      RETURN

  600 FORMAT(' ','# PENETRATION LENGTH =',1F6.3,'(M)')
      END SUBROUTINE TRPELB

!     ***********************************************************

!           ABLATION MODEL : BY S.K.HO AND L.JOHN PERKINS

!     ***********************************************************

      SUBROUTINE TRPELC

      USE TRCOMM, ONLY : AMP, DR, DVRHO, NRMAX, NSMAX, PA, PELPAT,&
           & PELRAD, PELVEL, PI, RA, RN, RT, SPE, rkind
      IMPLICIT NONE
      REAL(rkind)    :: ADIVE, AMPEL, ANA, ANE, ANP, PAP, PP, ROS, RP,&
           & RPDOT, RPE, RPPRE, SPEL, TE
      INTEGER :: NR

      NR = NRMAX
      RP = PELRAD
      PAP= PA(2)*PELPAT(2)+PA(3)*PELPAT(3)+PA(4)*PELPAT(4)
      ROS= 89.D0*PAP
      AMPEL= AMP*PAP

      ANP = ROS/AMPEL*1.D-20

  100 ANE=RN(NR,1)
      ANA=RN(NR,4)
      TE =RT(NR,1)
      ADIVE = ANA/ANE

      IF(ADIVE.GT.3.D-2) THEN
         PP=2.82D0*(ADIVE-3.D-2)+5.00D0*3.D-2
      ELSE
         PP=5.00D0*ADIVE
      ENDIF

      RPDOT=1.72D-8*(1.D0+PP)*(RP*1.D2)**(-2.D0/3.D0)*(ANE*1.D20)**(1.D0/3.D0)*ABS(TE)**1.64D0

      RPPRE=RP
      RP=RP-RPDOT*DR*RA/PELVEL

      IF(RP.LE.0.D0) THEN
         RPE=RA-RA*DR*(NR-1)+RA*DR*RPPRE/(RPPRE-RP)
         GOTO 9000
      ENDIF

      IF(NR.EQ.1) THEN
         RPE=RA
         GOTO 9000
      ENDIF

      SPEL=ANP*RP*0.5D0*(RP+RPPRE)*RPDOT*4.D0*PI*RA/(DVRHO(NR)*PELVEL)
      SPE(NR,1:NSMAX)=PELPAT(1:NSMAX)*SPEL

      WRITE(6,*) 'TRPELC:'
      WRITE(6,'(5ES12.4)') SPE(1:NRMAX,1)

      NR = NR-1
      GOTO 100

 9000 WRITE(6,600) RPE
      RETURN

  600 FORMAT(' ','# PENETRATION LENGTH =',1F6.3,'(M)')
      END SUBROUTINE TRPELC
