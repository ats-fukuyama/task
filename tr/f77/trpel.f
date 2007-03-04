C     $Id$
C
C     ***********************************************************
C
C           PELLET INJECTION
C
C     ***********************************************************
C
      SUBROUTINE TRPELT
C
      INCLUDE 'trcomm.inc'
C
      IF(PELTOT.LE.0.D0) RETURN
C
      IF(MDLPEL.EQ.0) RETURN
C
      WRITE(6,*)'# PELLET INJECTED.'
C
      IF(MDLPEL.EQ.1) THEN
         CALL TRPELA
      ELSEIF(MDLPEL.EQ.2) THEN
         CALL TRPELB
      ELSEIF(MDLPEL.EQ.3) THEN
         CALL TRPELC
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
C           PELLET INJECTION (GAUSSIAN PROFILE)
C
C     ***********************************************************
C
      SUBROUTINE TRPELA
C
      INCLUDE 'trcomm.inc'
C
      FSUM = 0.D0
      DO NR=1,NRMAX
         FSUM=FSUM+DEXP(-((RA*RM(NR)-PELR0)/PELRW)**2)*DVRHO(NR)*DR
      ENDDO
C
      S0  =PELTOT/FSUM
C
      DO NR=1,NRMAX
         SPEL=S0*DEXP(-((RA*RM(NR)-PELR0)/PELRW)**2)
      DO NS=1,NSM
         SPE(NR,NS)=PELPAT(NS)*SPEL
      ENDDO
      ENDDO
C
      RETURN
      END
C
C     ***********************************************************
C
C           PELLET ABLATION MODEL (WAKATANI NAKAMURA MODEL)
C
C     ***********************************************************
C
      SUBROUTINE TRPELB
C
      INCLUDE 'trcomm.inc'
C
      NR = NRMAX
      RP = PELRAD
      PAP= PA(2)*PELPAT(2)+PA(3)*PELPAT(3)+PA(4)*PELPAT(4)
      ROS= 89.D0*PAP
      AMP= AMM*PAP
C
      AMD = PA(2)*AMM
      AMT = PA(3)*AMM
      AMA = PA(4)*AMM
      AMB = PA(2)*AMM
      ANS = ROS/AMP*1.D-20
      VF  = SQRT(2.D0*3.5D3*RKEV/AMA)
      VB  = SQRT(2.D0*PNBENG*RKEV/AMB)
C
 1000 ANE=RN(NR,1)
      TE=RT(NR,1)
      P1   = 3.D0*SQRT(0.5D0*PI)*AME/ANE
     &      *(ABS(TE)*RKEV/AME)**1.5D0
C
      VCR   = (P1*(RN(NR,2)*PZ(2)   **2/AMD
     &           +RN(NR,3)*PZ(3)   **2/AMT
     &           +RN(NR,4)*PZ(4)   **2/AMA
     &           +ANFE(NR)*PZFE(NR)**2/52.D0
     &           +ANC (NR)*PZC(NR) **2/12.D0))**(1.D0/3.D0)
C
      TAUS = 0.2D0*PA(2)*ABS(TE)**1.5D0
     &      /(PZ(2)**2*ANE*COULOG(1,2,ANE,TE))
C
      ANFAST=(SNB(NR)*LOG(1.D0+(VB/VCR)**3)
     &       +SNF(NR)*LOG(1.D0+(VF/VCR)**3))*TAUS/3.D0
C
      A1   = SNB(NR)*(0.5D0*AMB*VCR*VCR/AEE)**3
      B1   = SNF(NR)*(0.5D0*AMA*VCR*VCR/AEE)**3
      A2   = 0.5D0*(VB/VCR)**6-(VB/VCR)**3+LOG(1.D0+(VB/VCR)**3)
      B2   = 0.5D0*(VF/VCR)**6-(VF/VCR)**3+LOG(1.D0+(VF/VCR)**3)
      P2   = (A1*A2+B1*B2)*TAUS/3.D0
C
      EFAST= (P2/ANFAST)**(1.D0/3.D0)
C
      A3   = SNB(NR)*1.D20*AMB*((VB/VCR)**3-LOG(1.D0+(VB/VCR)**3))
      B3   = SNF(NR)*1.D20*AMA*((VF/VCR)**3-LOG(1.D0+(VF/VCR)**3))
C
      QFAST= (A3+B3)*TAUS*VCR**3/(AEE*24.D0)
C
      RPDOT= 4.55D6*(2.D0*AMP)**(2.D0/3.D0)*QFAST**(1.D0/3.D0)
     &       *EFAST**0.4D0*RP**(-2.D0/3.D0)/ROS
C
      RPPRE=RP
      RP=RP-RPDOT*DR*RA/PELVEL
C
      IF(RP.LE.0.D0) THEN
         RPE=RA-RA*DR*(NR-1)+RA*DR*RPPRE/(RPPRE-RP)
         GOTO 9000
      ENDIF
C
      IF(NR.EQ.1) THEN
         RPE=RA
         GOTO 9000
      ENDIF
C
      SPEL=ANS*RP*0.5D0*(RP+RPPRE)*RPDOT*4.D0*PI*RA/(DVRHO(NR)*PELVEL)
      DO NS=1,NSM
         SPE(NR,NS)=PELPAT(NS)*SPEL
      ENDDO
C
      NR = NR-1
      GOTO 1000
C
 9000 WRITE(6,600) RPE
      RETURN
C
  600 FORMAT(' ','# PENETRATION LENGTH =',1F6.3,'(M)')
      END
C
C     ***********************************************************
C
C           ABLATION MODEL : BY S.K.HO AND L.JOHN PERKINS
C
C     ***********************************************************
C
      SUBROUTINE TRPELC
C
      INCLUDE 'trcomm.inc'
C
      NR = NRMAX
      RP = PELRAD
      PAP= PA(2)*PELPAT(2)+PA(3)*PELPAT(3)+PA(4)*PELPAT(4)
      ROS= 89.D0*PAP
      AMP= AMM*PAP
C
      ANP = ROS/AMP*1.D-20
C
  100 ANE=RN(NR,1)
      ANA=RN(NR,4)
      TE =RT(NR,1)
      ADIVE = ANA/ANE
C
      IF(ADIVE.GT.3.D-2) THEN
         PP=2.82D0*(ADIVE-3.D-2)+5.00D0*3.D-2
      ELSE
         PP=5.00D0*ADIVE
      ENDIF
C
      RPDOT=1.72D-8*(1.D0+PP)*(RP*1.D2)**(-2.D0/3.D0)
     &     *(ANE*1.D20)**(1.D0/3.D0)*ABS(TE)**1.64D0
C
      RPPRE=RP
      RP=RP-RPDOT*DR*RA/PELVEL
C
      IF(RP.LE.0.D0) THEN
         RPE=RA-RA*DR*(NR-1)+RA*DR*RPPRE/(RPPRE-RP)
         GOTO 9000
      ENDIF
C
      IF(NR.EQ.1) THEN
         RPE=RA
         GOTO 9000
      ENDIF
C
      SPEL=ANP*RP*0.5D0*(RP+RPPRE)*RPDOT*4.D0*PI*RA/(DVRHO(NR)*PELVEL)
      DO NS=1,NSM
         SPE(NR,NS)=PELPAT(NS)*SPEL
      ENDDO
C
      NR = NR-1
      GOTO 100
C
 9000 WRITE(6,600) RPE
      RETURN
C
  600 FORMAT(' ','# PENETRATION LENGTH =',1F6.3,'(M)')
      END
