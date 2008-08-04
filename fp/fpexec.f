C     $Id$
C
C     **************************
C        EXECUTE TIME ADVANCE
C     **************************
C
      SUBROUTINE FPEXEC(NSA,IERR)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION AAM(4*NTHM-1,NTHM*NPM)
C
      DO NR=1,NRMAX
C
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            NMA(NTH,NP,1)=NTH+NTHMAX*(NP-1)
         ENDDO
         ENDDO
         NMMAX=NTHMAX*NPMAX
         NWMAX=4*NTHMAX-1
         NWCEN=2*NTHMAX
         DO NM=1,NMMAX
            DO NL=1,NWMAX
               AAM(NL,NM)=0.D0
            ENDDO
            DO NL=1,NLM
               LL(NM,NL)=0
               AL(NM,NL)=0.D0
            ENDDO
         ENDDO
C
         NLMAX=0
         IF(MODELA.EQ.0) THEN
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               CALL FPSETM(NTH,NP,NR,NSA,NL)
               NLMAX=MAX(NLMAX,NL)
            ENDDO
            ENDDO
         ELSE
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX/2
                  CALL FPSETM(NTH,NP,NR,NSA,NL)
                  NLMAX=MAX(NLMAX,NL)
               ENDDO
               DO NTH=ITU(NR)+1,NTHMAX
                  CALL FPSETM(NTH,NP,NR,NSA,NL)
                  NLMAX=MAX(NLMAX,NL)
               ENDDO
            ENDDO
         ENDIF
C
         IF(MODELA.EQ.0) THEN
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,1)
               AAM(NWCEN,NM)=1.D0-RIMPL*DELT*DL(NM)
               FM(NM)=F(NTH,NP,NR)
               BM(NM)=(1.D0+(1.D0-RIMPL)*DELT*DL(NM))*FM(NM)
     &               +DELT*SP(NTH,NP,NR,NSA)
            ENDDO
            ENDDO
         ELSE
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX/2
                  NM=NMA(NTH,NP,1)
                  AAM(NWCEN,NM)= RLAMDA(NTH,NR)-RIMPL       *DELT*DL(NM)
                  FM(NM)=F(NTH,NP,NR)
                  BM(NM)=(RLAMDA(NTH,NR)+(1.D0-RIMPL)*DELT*DL(NM))
     &                   *FM(NM)
     &                  +DELT*SP(NTH,NP,NR,NSA)
               ENDDO
               DO NTH=NTHMAX/2+1,ITU(NR)
                  NTHS=NTHMAX+1-NTH
                  NM=NMA(NTH,NP,1)
                  AAM(NWCEN,NM)=1.D0
                  AAM(NWCEN+NTHS-NTH,NM)=-1.D0
                  FM(NM)=F(NTH,NP,NR)
                  BM(NM)=0.D0
               ENDDO
               DO NTH=ITU(NR)+1,NTHMAX
                  NM=NMA(NTH,NP,1)
                  AAM(NWCEN,NM)= RLAMDA(NTH,NR)-RIMPL       *DELT*DL(NM)
                  FM(NM)=F(NTH,NP,NR)
                  BM(NM)=(RLAMDA(NTH,NR)+(1.D0-RIMPL)*DELT*DL(NM))
     &                   *FM(NM)
     &                  +DELT*SP(NTH,NP,NR,NSA)
               ENDDO
            ENDDO
         ENDIF     
C
         DO NM=1,NMMAX
         DO NL=1,NLMAX
            IF(LL(NM,NL).NE.0) THEN
               NP=(LL(NM,NL)-1)/NTHMAX+1
               NTH=MOD(LL(NM,NL)-1,NTHMAX)+1
               NLL=NTH+NTHMAX*(NP-1)-NM+NWCEN
               AAM(NLL,NM)=-RIMPL*DELT*AL(NM,NL)
               BM(NM)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,NL)*FM(LL(NM,NL))
            ENDIF
         ENDDO
         ENDDO
C
         CALL BANDRD(AAM,BM,NMMAX,NWMAX,4*NTHM-1,IERR)
C
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            NM=NMA(NTH,NP,1)
            F1(NTH,NP,NR)=BM(NM)
         ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C
C ***********************
C     SET ONE ELEMENT
C ***********************
C
      SUBROUTINE FPSETM(NTH,NP,NR,NSA,NL)
C
      INCLUDE 'fpcomm.inc'
C
      IF(NTH.EQ.ITL(NR)) THEN
         NTB=ITU(NR)
      ELSE
         NTB=0
      ENDIF
C
      NL=0
      NM=NMA(NTH,NP,1)
      PL=PM(NP)
      DPPM=PG(NP  )**2
      DPPP=PG(NP+1)**2
      DPTM=PG(NP  )
      DPTP=PG(NP+1)
      SL=SINM(NTH)
      DTPM=SING(NTH  )
      DTPP=SING(NTH+1)
      DTTM=SING(NTH  )/PL
      DTTP=SING(NTH+1)/PL
      IF(NTB.NE.0) THEN
         DTPM=0.5D0*DTPM
         DTTM=0.5D0*DTTM
      ENDIF
C      RL=RM(NR)
C      DRRM=RG(NR  )
C      DRRP=RG(NR+1)
      WPM=WEIGHP(NTH  ,NP  ,NR  ,NSA)
      WPP=WEIGHP(NTH  ,NP+1,NR  ,NSA)
      WTM=WEIGHT(NTH  ,NP  ,NR  ,NSA)
      WTP=WEIGHT(NTH+1,NP  ,NR  ,NSA)
C      WRM=WEIGHR(NTH  ,NP  ,NR  ,NSA)
C      WRP=WEIGHR(NTH  ,NP  ,NR+1,NSA)
      VPM=1.D0-WPM
      VPP=1.D0-WPP
      VTM=1.D0-WTM
      VTP=1.D0-WTP
C      VRM=1.D0-WPM
C      VRP=1.D0-WPP
      IF(NTB.NE.0) THEN
         WTB=WEIGHT(NTB+1,NP  ,NR  ,NSA)
         VTB=1.D0-WTB
      ENDIF
      DIVDPP=1.D0/(     PL*PL*DELP *DELP)
      DIVDPT=1.D0/(2.D0*PL*PL*DELP *DELTH)
      DIVDTP=1.D0/(2.D0*PL*SL*DELTH*DELP)
      DIVDTT=1.D0/(     PL*SL*DELTH*DELTH)
      DIVFPP=1.D0/(     PL*PL*DELP)
      DIVFTH=1.D0/(     PL*SL*DELTH)
C      DIVDRR=1.D0/(     RL   *DELR *DELR)
C      DIVFRR=1.D0/(     RL   *DELR)
      IF(NP.EQ.NPMAX) THEN
         DIVDTP=2.D0*DIVDTP
      ENDIF
C
C      IF(NR.GT.1) THEN
C         NL=NL+1
C         LL(NM,NL)=NMA(NTH,NP,NR-1)
C         AL(NM,NL)=+DRR(NTH  ,NP  ,NR,NSA)    *DIVDRR*DRRM
C     &             +FRR(NTH  ,NP  ,NR,NSA)*WRM*DIVFRR*DRRM
C         IF(ABS(AL(NM,NL)).LT.1.D-70) NL=NL-1
C      ENDIF
C
      IF(NP.NE.1.AND.NTH.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH-1,NP-1,1)
         AL(NM,NL)=+DPT(NTH  ,NP  ,NR,NSA)*WPM*DIVDPT*DPTM
     &             +DTP(NTH  ,NP  ,NR,NSA)*WTM*DIVDTP*DTPM
         IF(ABS(AL(NM,NL)).LT.1.D-70) NL=NL-1
      ENDIF
C
      IF(NP.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH,NP-1,1)
         AL(NM,NL)=+DPP(NTH  ,NP  ,NR,NSA)    *DIVDPP*DPPM
     &             +FPP(NTH  ,NP  ,NR,NSA)*WPM*DIVFPP*DPPM
     &             -DTP(NTH+1,NP  ,NR,NSA)*WTP*DIVDTP*DTPP
     &             +DTP(NTH  ,NP  ,NR,NSA)*VTM*DIVDTP*DTPM
         IF(NTB.NE.0) THEN
            AL(NM,NL)=AL(NM,NL)
     &             -DTP(NTB+1,NP  ,NR,NSA)*VTB*DIVDTP*DTPM
         ENDIF
      ENDIF
C
      IF(NP.NE.1.AND.NTH.NE.NTHMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP-1,1)
         AL(NM,NL)=-DPT(NTH  ,NP  ,NR,NSA)*WPM*DIVDPT*DPTM
     &             -DTP(NTH+1,NP  ,NR,NSA)*VTP*DIVDTP*DTPP
         IF(ABS(AL(NM,NL)).LT.1.D-70) NL=NL-1
      ENDIF
C
      IF(NP.NE.1.AND.NTB.NE.0) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTB+1,NP-1,1)
         AL(NM,NL)=-DTP(NTB+1,NP  ,NR,NSA)*WTB*DIVDTP*DTPM
         IF(ABS(AL(NM,NL)).LT.1.D-70) NL=NL-1
      ENDIF
C
      IF(NTH.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH-1,NP,1)
         AL(NM,NL)=+DTT(NTH  ,NP  ,NR,NSA)    *DIVDTT*DTTM
     &             +FTH(NTH  ,NP  ,NR,NSA)*WTM*DIVFTH*DTPM
     &             -DPT(NTH  ,NP+1,NR,NSA)*WPP*DIVDPT*DPTP
     &             +DPT(NTH  ,NP  ,NR,NSA)*VPM*DIVDPT*DPTM
         IF(NP.EQ.NPMAX) THEN
            AL(NM,NL)=AL(NM,NL)
     &             -DTP(NTH  ,NP  ,NR,NSA)*WTM*DIVDTP*DTPM
         ENDIF
      ENDIF
C
      IF(NTH.NE.NTHMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP,1)
         AL(NM,NL)=+DTT(NTH+1,NP  ,NR,NSA)    *DIVDTT*DTTP
     &             -FTH(NTH+1,NP  ,NR,NSA)*VTP*DIVFTH*DTPP
     &             +DPT(NTH  ,NP+1,NR,NSA)*WPP*DIVDPT*DPTP
     &             -DPT(NTH  ,NP  ,NR,NSA)*VPM*DIVDPT*DPTM
         IF(NP.EQ.NPMAX) THEN
            AL(NM,NL)=AL(NM,NL)
     &             +DTP(NTH+1,NP  ,NR,NSA)*VTP*DIVDTP*DTPP
         ENDIF
      ENDIF
C
      IF(NTB.NE.0) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTB+1,NP,1)
         AL(NM,NL)=+DTT(NTB+1,NP  ,NR,NSA)    *DIVDTT*DTTM
     &             -FTH(NTB+1,NP  ,NR,NSA)*WTB*DIVFTH*DTPM
         IF(NP.EQ.NPMAX) THEN
            AL(NM,NL)=AL(NM,NL)
     &             +DTP(NTB+1,NP  ,NR,NSA)*WTB*DIVDTP*DTPM
         ENDIF
      ENDIF
C
      IF(NP.NE.NPMAX.AND.NTH.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH-1,NP+1,1)
         AL(NM,NL)=-DPT(NTH  ,NP+1,NR,NSA)*VPP*DIVDPT*DPTP
     &             -DTP(NTH  ,NP  ,NR,NSA)*WTM*DIVDTP*DTPM
         IF(ABS(AL(NM,NL)).LT.1.D-70) NL=NL-1
      ENDIF
C
      IF(NP.NE.NPMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH,NP+1,1)
         AL(NM,NL)=+DPP(NTH  ,NP+1,NR,NSA)    *DIVDPP*DPPP
     &             -FPP(NTH  ,NP+1,NR,NSA)*VPP*DIVFPP*DPPP
     &             +DTP(NTH+1,NP  ,NR,NSA)*WTP*DIVDTP*DTPP
     &             -DTP(NTH  ,NP  ,NR,NSA)*VTM*DIVDTP*DTPM
         IF(NTB.NE.0) THEN
            AL(NM,NL)=AL(NM,NL)
     &             +DTP(NTB+1,NP  ,NR,NSA)*VTB*DIVDTP*DTPM
         ENDIF
      ENDIF
C
      IF(NP.NE.NPMAX.AND.NTH.NE.NTHMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP+1,1)
         AL(NM,NL)=+DPT(NTH  ,NP+1,NR,NSA)*VPP*DIVDPT*DPTP
     &             +DTP(NTH+1,NP  ,NR,NSA)*VTP*DIVDTP*DTPP
         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
            NL=NL-1
         ENDIF
      ENDIF
C
      IF(NP.NE.NPMAX.AND.NTB.NE.0) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTB+1,NP+1,1)
         AL(NM,NL)=+DTP(NTB+1,NP  ,NR,NSA)*WTB*DIVDTP*DTPM
         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
            NL=NL-1
         ENDIF
      ENDIF
C
C      IF(NR.LT.NRMAX) THEN
C         NL=NL+1
C         LL(NM,NL)=NMA(NTH,NP,NR+1)
C         AL(NM,NL)=+DRR(NTH  ,NP  ,NR+1,NSA)    *DIVDRR*DRRP
C     &             -FRR(NTH  ,NP  ,NR+1,NSA)*VRP*DIVFRR*DRRP
C         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
C            LL(NM,NL)=0
C            AL(NM,NL)=0.D0
C            NL=NL-1
C         ENDIF
C      ENDIF
C
      DL(NM)=-DPP(NTH  ,NP+1,NR  ,NSA)    *DIVDPP*DPPP
     &       -FPP(NTH  ,NP+1,NR  ,NSA)*WPP*DIVFPP*DPPP
     &       -DPP(NTH  ,NP  ,NR  ,NSA)    *DIVDPP*DPPM
     &       +FPP(NTH  ,NP  ,NR  ,NSA)*VPM*DIVFPP*DPPM
     &       -DTT(NTH+1,NP  ,NR  ,NSA)    *DIVDTT*DTTP
     &       -FTH(NTH+1,NP  ,NR  ,NSA)*WTP*DIVFTH*DTPP
     &       -DTT(NTH  ,NP  ,NR  ,NSA)    *DIVDTT*DTTM
     &       +FTH(NTH  ,NP  ,NR  ,NSA)*VTM*DIVFTH*DTPM
C     &       -DRR(NTH  ,NP  ,NR+1,NSA)    *DIVDRR*DRRP
C     &       -FRR(NTH  ,NP  ,NR+1,NSA)*WRP*DIVFRR*DRRP
C     &       -DRR(NTH  ,NP  ,NR  ,NSA)    *DIVDRR*DRRM
C     &       +FRR(NTH  ,NP  ,NR  ,NSA)*VRM*DIVFRR*DRRM
      IF(NTB.NE.0) THEN
         DL(NM)=DL(NM)
     &       -DTT(NTB+1,NP  ,NR  ,NSA)    *DIVDTT*DTTM
     &       -FTH(NTB+1,NP  ,NR  ,NSA)*VTB*DIVFTH*DTPM
      ENDIF
      IF(NP.EQ.NPMAX) THEN
         DL(NM)=DL(NM)
     &       +DTP(NTH+1,NP  ,NR  ,NSA)*WTP*DIVDTP*DTPP
     &       -DTP(NTH  ,NP  ,NR  ,NSA)*VTM*DIVDTP*DTPM
         IF(NTB.NE.0) THEN
            DL(NM)=DL(NM)
     &       +DTP(NTB+1,NP  ,NR  ,NSA)*VTB*DIVDTP*DTPM
         ENDIF
      ENDIF
C      IF(NR.EQ.1) THEN
C         SP(NTH,NP,NR,NSA)=SP(NTH,NP,NR,NSA)+FS1(NTH,NP)
C     &         *(DRR(NTH  ,NP  ,NR  ,NSA)    *DIVDRR*DRRM
C     &          +FRR(NTH  ,NP  ,NR  ,NSA)*WRM*DIVFRR*DRRM)
C      ENDIF
C      IF(NR.EQ.NRMAX) THEN
C         SP(NTH,NP,NR,NSA)=SP(NTH,NP,NR,NSA)+FS2(NTH,NP,NSA)
C     &         *(DRR(NTH  ,NP  ,NR+1,NSA)    *DIVDRR*DRRP
C     &          -FRR(NTH  ,NP  ,NR+1,NSA)*VRP*DIVFRR*DRRP)
C      ENDIF
      RETURN
      END
