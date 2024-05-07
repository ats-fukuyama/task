!     ***********************************************************

!           PELLET INJECTION

!     ***********************************************************

      SUBROUTINE TRPELT

      USE TRCOMM
      IMPLICIT NONE
      INTEGER:: npel,ns,nr

      IF(npelmax.EQ.0) RETURN

      SPE(1:nsmax,1:nrmax)=0.D0
      DO npel=1,npelmax
         SPE_nsnpelnr(1:nsmax,npel,1:nrmax)=0.D0
      END DO
      
      DO npel=1,npelmax
         IF(PELIN(npel).LE.0.D0) exit
         IF(MDLPEL(npel).EQ.0) exit

         WRITE(6,*)'# PELLET INJECTED.'

         IF(MDLPEL(npel).EQ.1) THEN
            CALL TRPELA(npel)
         ELSEIF(MDLPEL(npel).EQ.2) THEN
            CALL TRPELB(npel)
         ELSEIF(MDLPEL(npel).EQ.3) THEN
            CALL TRPELC(npel)
         ENDIF
      END DO

      DO NR=1,NRMAX
         DO NS=1,NSMAX
            SPE(ns,nr)=SPE(ns,nr)+SUM(SPE_NSNPELNR(ns,1:npelmax,nr))
         END DO
      END DO
      RETURN
      END  SUBROUTINE TRPELT

!     ***********************************************************

!           PELLET INJECTION (GAUSSIAN PROFILE)

!     ***********************************************************

      SUBROUTINE TRPELA(npel)

      USE TRCOMM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: npel
      REAL(rkind):: FSUM,S0,SPEL
      INTEGER:: NR,NS


      FSUM = 0.D0
      DO NR=1,NRMAX
         FSUM=FSUM &
              +DEXP(-((RA*RM(NR)-PELR0(npel))/PELRW(npel))**2)*DVRHO(NR)*DR
      ENDDO

      S0=PELIN(npel)/FSUM

      DO NR=1,NRMAX
         SPEL=S0*DEXP(-((RA*RM(NR)-PELR0(npel))/PELRW(npel))**2)
      DO NS=1,NSMAX
         SPE_NSNPELNR(NS,NPEL,NR)=PELPAT(NS,npel)*SPEL
      ENDDO
      ENDDO

      WRITE(6,*) 'TRPELA:'
      WRITE(6,'(5ES12.4)') SPE(1:NRMAX,1)

      RETURN
      END  SUBROUTINE TRPELA

!     ***********************************************************

!           PELLET ABLATION MODEL (WAKATANI NAKAMURA MODEL)

!     ***********************************************************

      SUBROUTINE TRPELB(npel)

      USE TRCOMM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: npel
      REAL(rkind)    :: A1, A2, A3, AMB, AMF, AMPEL, ANE,&
            ANFAST, ANS, B1, B2, B3, EFAST, P1, P2, PAP, QFAST, ROS,&
            RPEL, RPELDOT, RPELPRE, SPEL, TAUS, TE, VCR, VB, VF, SUM, RNZZM
      REAL(rkind),ALLOCATABLE:: VB_NNB(:)
      INTEGER :: NR,NS,NNB,NNF,NSB,NSF
      REAL(rkind)    :: COULOG   ! FUNCTION

      IF(NNBMAX.GT.0) ALLOCATE(VB_NNB(NNBMAX))

      NR = NRMAX
      RPEL = PELRAD(npel)
      PAP=0.D0
      DO NS=1,NSMAX
         IF(NS.NE.NS_e) THEN
            PAP=PAP+PM(NS)*PELPAT(NS,npel)
         END IF
      END DO
      ROS= 89.D0*PAP
      AMPEL= AMP*PAP
      ANS = ROS/AMPEL*1.D-20

 1000 ANE=RN(NR,NS_e)
      TE=RT(NR,NS_e)
      P1   = 3.D0*SQRT(0.5D0*PI)*AME/ANE *(ABS(TE)*RKEV/AME)**1.5D0
      RNZZM=0.D0
      DO NS=1,NSMAX
         RNZZM=RNZZM+RN(NR,NS)*PZ(NS)**2/(PM(NS)*AMP)
      END DO
      VCR=(P1*RNZZM)**(1.D0/3.D0)

      ANFAST=0.D0
      SUM=0.D0
      A3=0.D0
      B3=0.D0
      DO NNB=1,NNBMAX
         NSB=NS_NNB(NNB)
         AMB=PM(NSB)*AMP
         VB=SQRT(2.D0*PNBENG(NNB)*RKEV/AMB)
         TAUS = 0.2D0*PM(NSB)*ABS(TE)**1.5D0 &
              /(PZ(NNB)**2*ANE*COULOG(NS_e,NSB,ANE,TE))
         ANFAST=ANFAST+SNB_NNBNR(NNB,NR)*LOG(1.D0+(VB/VCR)**3)*TAUS/3.D0
         A1=SNB_NNBNR(NNB,NR)*(0.5D0*AMB*VCR*VCR/AEE)**3
         A2=0.5D0*(VB/VCR)**6-(VB/VCR)**3+LOG(1.D0+(VB/VCR)**3)
         SUM=SUM+A1*A2*TAUS/3.D0
         A3=A3+SNB_NNBNR(NNB,NR)*1.D20*AMB*((VB/VCR)**3-LOG(1.D0+(VB/VCR)**3))
      END DO
      DO NNF=1,NNFMAX
         NSF=NS_NNF(NNF)
         AMF=PM(NSF)*AMP
         TAUS = 0.2D0*PM(NSF)*ABS(TE)**1.5D0 &
              /(PZ(NSF)**2*ANE*COULOG(NS_e,NSF,ANE,TE))
         ANFAST=ANFAST+SNF_NNFNR(NNF,NR)*LOG(1.D0+(VF/VCR)**3)*TAUS/3.D0
         VF=SQRT(2.D0*ENF_NNF(NNF)*RKEV/AMF)
         B1=SNF_NNFNR(NNF,NR)*(0.5D0*AMF*VCR*VCR/AEE)**3
         B2=0.5D0*(VF/VCR)**6-(VF/VCR)**3+LOG(1.D0+(VF/VCR)**3)
         SUM=SUM+B1*B2*TAUS/3.D0
         B3=B3+SNF_NNFNR(NNF,NR)*1.D20*AMF*((VF/VCR)**3-LOG(1.D0+(VF/VCR)**3))
      END DO
      P2=SUM
      EFAST= (P2/ANFAST)**(1.D0/3.D0)
      QFAST= (A3+B3)*TAUS*VCR**3/(AEE*24.D0)

      RPELDOT= 4.55D6*(2.D0*AMPEL)**(2.D0/3.D0)*QFAST**(1.D0/3.D0) &
     &       *EFAST**0.4D0*RPEL**(-2.D0/3.D0)/ROS

      RPELPRE=RPEL
      RPEL=RPEL-RPELDOT*DR*RA/PELVEL(npel)

      IF(RPEL.LE.0.D0) THEN
         RPEL=RA-RA*DR*(NR-1)+RA*DR*RPELPRE/(RPELPRE-RPEL)
         GOTO 9000
      ENDIF

      IF(NR.EQ.1) THEN
         RPEL=RA
         GOTO 9000
      ENDIF

      SPEL=ANS*RPEL*0.5D0*(RPEL+RPELPRE)*RPELDOT*4.D0*PI*RA &
           /(DVRHO(NR)*PELVEL(npel))
      DO NS=1,NSMAX
         SPE(NR,NS)=PELPAT(NS,npel)*SPEL
      ENDDO

      NR = NR-1
      GOTO 1000

 9000 WRITE(6,600) RPEL
      RETURN

  600 FORMAT(' ','# PENETRATION LENGTH =',1F6.3,'(M)')
      END SUBROUTINE TRPELB

!     ***********************************************************

!           ABLATION MODEL : BY S.K.HO AND L.JOHN PERKINS

!     ***********************************************************

      SUBROUTINE TRPELC(npel)

      USE TRCOMM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: npel
      REAL(rkind)    :: ADIVE, AMPEL, ANA, ANE, ANP, PAP, PP, ROS, &
           & RPELDOT, RPEL, RPELPRE, SPEL, TE
      INTEGER :: NR,NS

      NR = NRMAX
      RPEL=PELRAD(npel)
      PAP=0.D0
      DO NS=1,NSMAX
         IF(NS.NE.NS_e) THEN
            PAP=PAP+PM(NS)*PELPAT(NS,npel)
         END IF
      END DO
      ROS= 89.D0*PAP
      AMPEL= AMP*PAP
      ANP = ROS/AMPEL*1.D-20

  100 ANE=RN(NR,NS_e)
      TE =RT(NR,NS_e)
      ANA=RN(NR,NS_A)
      ADIVE = ANA/ANE

      IF(ADIVE.GT.3.D-2) THEN
         PP=2.82D0*(ADIVE-3.D-2)+5.00D0*3.D-2
      ELSE
         PP=5.00D0*ADIVE
      ENDIF

      RPELDOT=1.72D-8*(1.D0+PP)*(RPEL*1.D2)**(-2.D0/3.D0) &
                     *(ANE*1.D20)**(1.D0/3.D0)*ABS(TE)**1.64D0

      RPELPRE=RPEL
      RPEL=RPEL-RPELDOT*DR*RA/PELVEL(npel)

      IF(RPEL.LE.0.D0) THEN
         RPEL=RA-RA*DR*(NR-1)+RA*DR*RPELPRE/(RPELPRE-RPEL)
         GOTO 9000
      ENDIF

      IF(NR.EQ.1) THEN
         RPEL=RA
         GOTO 9000
      ENDIF

      SPEL=ANP*RPEL*0.5D0*(RPEL+RPELPRE)*RPELDOT*4.D0*PI*RA &
           /(DVRHO(NR)*PELVEL(npel))
      SPE(NR,1:NSMAX)=PELPAT(1:NSMAX,npel)*SPEL

      WRITE(6,*) 'TRPELC:'
      WRITE(6,'(5ES12.4)') SPE(1:NRMAX,1)

      NR = NR-1
      GOTO 100

 9000 WRITE(6,600) RPEL
      RETURN

  600 FORMAT(' ','# PENETRATION LENGTH =',1F6.3,'(M)')
      END SUBROUTINE TRPELC
