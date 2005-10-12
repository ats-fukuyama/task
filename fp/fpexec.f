C     $Id$
C
C     **************************
C        EXECUTE TIME ADVANCE
C     **************************
C
      SUBROUTINE FPEXECX(IERR)
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
         ENDDO
C
         NLMAX=0
         IF(MODELA.EQ.0) THEN
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               CALL FPSETM(NTH,NP,NR,NL)
               NLMAX=MAX(NLMAX,NL)
            ENDDO
            ENDDO
         ELSE
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX/2
                  CALL FPSETM(NTH,NP,NR,NL)
                  NLMAX=MAX(NLMAX,NL)
               ENDDO
               DO NTH=ITU(NR)+1,NTHMAX
                  CALL FPSETM(NTH,NP,NR,NL)
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
     &               +DELT*SP(NTH,NP,NR)
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
     &                  +DELT*SP(NTH,NP,NR)
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
     &                  +DELT*SP(NTH,NP,NR)
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
      SUBROUTINE FPEXEC(NOCONV)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION EPS(2)
C
      DO NR=1,NRMAX
C
         IF(MODELA.EQ.0) THEN
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               NMA(NTH,NP,1)=NTH+NTHMAX*(NP-1)
            ENDDO
            ENDDO
            NMMAX=NTHMAX*NPMAX
         ELSE
            NM=0
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX/2
                  NM=NM+1
                  NMA(NTH,NP,1)=NM
               ENDDO
               DO NTH=ITU(NR)+1,NTHMAX
                  NM=NM+1
                  NMA(NTH,NP,1)=NM
               ENDDO
               DO NTH=NTHMAX/2+1,ITU(NR)
                  NMA(NTH,NP,1)=NMA(NTHMAX-NTH+1,NP,1)
               ENDDO
            ENDDO
            NMMAX=NM
         ENDIF
C
         IF(MODELA.EQ.0) THEN
            NLMAX=0
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               CALL FPSETM(NTH,NP,NR,NL)
               NLMAX=MAX(NLMAX,NL)
            ENDDO
            ENDDO
         ELSE
            NLMAX=0
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX/2
                  CALL FPSETM(NTH,NP,NR,NL)
                  NLMAX=MAX(NLMAX,NL)
               ENDDO
               DO NTH=ITU(NR)+1,NTHMAX
                  CALL FPSETM(NTH,NP,NR,NL)
                  NLMAX=MAX(NLMAX,NL)
               ENDDO
            ENDDO
         ENDIF
C
         NOCONV=0
         IF(MODELA.EQ.0) THEN
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,1)
               DM(NM)=1.D0-RIMPL*DELT*DL(NM)
               FM(NM)=F(NTH,NP,NR)
               BM(NM)=(1.D0+(1.D0-RIMPL)*DELT*DL(NM))*FM(NM)
     &               +DELT*SP(NTH,NP,NR)
            ENDDO
            ENDDO
         ELSE
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX/2
                  NM=NMA(NTH,NP,1)
                  DM(NM)= RLAMDA(NTH,NR)-RIMPL       *DELT*DL(NM)
                  FM(NM)=F(NTH,NP,NR)
                  BM(NM)=(RLAMDA(NTH,NR)+(1.D0-RIMPL)*DELT*DL(NM))
     &                   *FM(NM)
     &                  +DELT*SP(NTH,NP,NR)
               ENDDO
               DO NTH=ITU(NR)+1,NTHMAX
                  NM=NMA(NTH,NP,1)
                  DM(NM)= RLAMDA(NTH,NR)-RIMPL       *DELT*DL(NM)
                  FM(NM)=F(NTH,NP,NR)
                  BM(NM)=(RLAMDA(NTH,NR)+(1.D0-RIMPL)*DELT*DL(NM))
     &                   *FM(NM)
     &                  +DELT*SP(NTH,NP,NR)
               ENDDO
            ENDDO
         ENDIF     
C
         DO NM=1,NMMAX
         DO NL=1,NLMAX
            LM(NM,NL)=LL(NM,NL)
            IF(LL(NM,NL).EQ.0) THEN
               AM(NM,NL)=0.D0
            ELSE
               AM(NM,NL)=-RIMPL*DELT*AL(NM,NL)
               BM(NM)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,NL)*FM(LL(NM,NL))
            ENDIF
         ENDDO
         ENDDO
C
         EPS(1)=EPSM
         PARM=0.D0
         ITR=0
C
C         WRITE(6,*) 'PCGPME : NMMAX = ',NMMAX,'   NLMAX = ',NLMAX
         CALL PCGPME(AM,NLMAX,NMM,LM,DM,NMMAX,NMM2,BM,FM,EPS,PARM,
     &               ITR,WK1,WK2,IER)
C
         IF(IER.NE.0) THEN
            WRITE(6,*) 'XX PCGPME ERROR: IER = ',IER
            IF(IER.EQ.3000) THEN
               NOCONV=1
            ENDIF
         ENDIF
         IF(EPS(1).GT.EPSM) THEN
            WRITE(6,*) 'PCGPME : ITR = ',ITR,'   ERROR = ',EPS(1)
         ENDIF
C
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            NM=NMA(NTH,NP,1)
            F1(NTH,NP,NR)=FM(NM)
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
      SUBROUTINE FPSETM(NTH,NP,NR,NL)
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
      WPM=WEIGHP(NTH  ,NP  ,NR  )
      WPP=WEIGHP(NTH  ,NP+1,NR  )
      WTM=WEIGHT(NTH  ,NP  ,NR  )
      WTP=WEIGHT(NTH+1,NP  ,NR  )
C      WRM=WEIGHR(NTH  ,NP  ,NR  )
C      WRP=WEIGHR(NTH  ,NP  ,NR+1)
      VPM=1.D0-WPM
      VPP=1.D0-WPP
      VTM=1.D0-WTM
      VTP=1.D0-WTP
C      VRM=1.D0-WPM
C      VRP=1.D0-WPP
      IF(NTB.NE.0) THEN
         WTB=WEIGHT(NTB+1,NP  ,NR  )
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
      DO 1000 I=1,NLM
         LL(NM,I)=0
         AL(NM,I)=0.D0
 1000 CONTINUE
C
C      IF(NR.GT.1) THEN
C         NL=NL+1
C         LL(NM,NL)=NMA(NTH,NP,NR-1)
C         AL(NM,NL)=+DRR(NTH  ,NP  ,NR)    *DIVDRR*DRRM
C     &             +FRR(NTH  ,NP  ,NR)*WRM*DIVFRR*DRRM
C         IF(ABS(AL(NM,NL)).LT.1.D-70) NL=NL-1
C      ENDIF
C
      IF(NP.NE.1.AND.NTH.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH-1,NP-1,1)
         AL(NM,NL)=+DPT(NTH  ,NP  ,NR)*WPM*DIVDPT*DPTM
     &             +DTP(NTH  ,NP  ,NR)*WTM*DIVDTP*DTPM
         IF(ABS(AL(NM,NL)).LT.1.D-70) NL=NL-1
      ENDIF
C
      IF(NP.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH,NP-1,1)
         AL(NM,NL)=+DPP(NTH  ,NP  ,NR)    *DIVDPP*DPPM
     &             +FPP(NTH  ,NP  ,NR)*WPM*DIVFPP*DPPM
     &             -DTP(NTH+1,NP  ,NR)*WTP*DIVDTP*DTPP
     &             +DTP(NTH  ,NP  ,NR)*VTM*DIVDTP*DTPM
         IF(NTB.NE.0) THEN
            AL(NM,NL)=AL(NM,NL)
     &             -DTP(NTB+1,NP  ,NR)*VTB*DIVDTP*DTPM
         ENDIF
      ENDIF
C
      IF(NP.NE.1.AND.NTH.NE.NTHMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP-1,1)
         AL(NM,NL)=-DPT(NTH  ,NP  ,NR)*WPM*DIVDPT*DPTM
     &             -DTP(NTH+1,NP  ,NR)*VTP*DIVDTP*DTPP
         IF(ABS(AL(NM,NL)).LT.1.D-70) NL=NL-1
      ENDIF
C
      IF(NP.NE.1.AND.NTB.NE.0) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTB+1,NP-1,1)
         AL(NM,NL)=-DTP(NTB+1,NP  ,NR)*WTB*DIVDTP*DTPM
         IF(ABS(AL(NM,NL)).LT.1.D-70) NL=NL-1
      ENDIF
C
      IF(NTH.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH-1,NP,1)
         AL(NM,NL)=+DTT(NTH  ,NP  ,NR)    *DIVDTT*DTTM
     &             +FTH(NTH  ,NP  ,NR)*WTM*DIVFTH*DTPM
     &             -DPT(NTH  ,NP+1,NR)*WPP*DIVDPT*DPTP
     &             +DPT(NTH  ,NP  ,NR)*VPM*DIVDPT*DPTM
         IF(NP.EQ.NPMAX) THEN
            AL(NM,NL)=AL(NM,NL)
     &             -DTP(NTH  ,NP  ,NR)*WTM*DIVDTP*DTPM
         ENDIF
      ENDIF
C
      IF(NTH.NE.NTHMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP,1)
         AL(NM,NL)=+DTT(NTH+1,NP  ,NR)    *DIVDTT*DTTP
     &             -FTH(NTH+1,NP  ,NR)*VTP*DIVFTH*DTPP
     &             +DPT(NTH  ,NP+1,NR)*WPP*DIVDPT*DPTP
     &             -DPT(NTH  ,NP  ,NR)*VPM*DIVDPT*DPTM
         IF(NP.EQ.NPMAX) THEN
            AL(NM,NL)=AL(NM,NL)
     &             +DTP(NTH+1,NP  ,NR)*VTP*DIVDTP*DTPP
         ENDIF
      ENDIF
C
      IF(NTB.NE.0) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTB+1,NP,1)
         AL(NM,NL)=+DTT(NTB+1,NP  ,NR)    *DIVDTT*DTTM
     &             -FTH(NTB+1,NP  ,NR)*WTB*DIVFTH*DTPM
         IF(NP.EQ.NPMAX) THEN
            AL(NM,NL)=AL(NM,NL)
     &             +DTP(NTB+1,NP  ,NR)*WTB*DIVDTP*DTPM
         ENDIF
      ENDIF
C
      IF(NP.NE.NPMAX.AND.NTH.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH-1,NP+1,1)
         AL(NM,NL)=-DPT(NTH  ,NP+1,NR)*VPP*DIVDPT*DPTP
     &             -DTP(NTH  ,NP  ,NR)*WTM*DIVDTP*DTPM
         IF(ABS(AL(NM,NL)).LT.1.D-70) NL=NL-1
      ENDIF
C
      IF(NP.NE.NPMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH,NP+1,1)
         AL(NM,NL)=+DPP(NTH  ,NP+1,NR)    *DIVDPP*DPPP
     &             -FPP(NTH  ,NP+1,NR)*VPP*DIVFPP*DPPP
     &             +DTP(NTH+1,NP  ,NR)*WTP*DIVDTP*DTPP
     &             -DTP(NTH  ,NP  ,NR)*VTM*DIVDTP*DTPM
         IF(NTB.NE.0) THEN
            AL(NM,NL)=AL(NM,NL)
     &             +DTP(NTB+1,NP  ,NR)*VTB*DIVDTP*DTPM
         ENDIF
      ENDIF
C
      IF(NP.NE.NPMAX.AND.NTH.NE.NTHMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP+1,1)
         AL(NM,NL)=+DPT(NTH  ,NP+1,NR)*VPP*DIVDPT*DPTP
     &             +DTP(NTH+1,NP  ,NR)*VTP*DIVDTP*DTPP
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
         AL(NM,NL)=+DTP(NTB+1,NP  ,NR)*WTB*DIVDTP*DTPM
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
C         AL(NM,NL)=+DRR(NTH  ,NP  ,NR+1)    *DIVDRR*DRRP
C     &             -FRR(NTH  ,NP  ,NR+1)*VRP*DIVFRR*DRRP
C         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
C            LL(NM,NL)=0
C            AL(NM,NL)=0.D0
C            NL=NL-1
C         ENDIF
C      ENDIF
C
      DL(NM)=-DPP(NTH  ,NP+1,NR  )    *DIVDPP*DPPP
     &       -FPP(NTH  ,NP+1,NR  )*WPP*DIVFPP*DPPP
     &       -DPP(NTH  ,NP  ,NR  )    *DIVDPP*DPPM
     &       +FPP(NTH  ,NP  ,NR  )*VPM*DIVFPP*DPPM
     &       -DTT(NTH+1,NP  ,NR  )    *DIVDTT*DTTP
     &       -FTH(NTH+1,NP  ,NR  )*WTP*DIVFTH*DTPP
     &       -DTT(NTH  ,NP  ,NR  )    *DIVDTT*DTTM
     &       +FTH(NTH  ,NP  ,NR  )*VTM*DIVFTH*DTPM
C     &       -DRR(NTH  ,NP  ,NR+1)    *DIVDRR*DRRP
C     &       -FRR(NTH  ,NP  ,NR+1)*WRP*DIVFRR*DRRP
C     &       -DRR(NTH  ,NP  ,NR  )    *DIVDRR*DRRM
C     &       +FRR(NTH  ,NP  ,NR  )*VRM*DIVFRR*DRRM
      IF(NTB.NE.0) THEN
         DL(NM)=DL(NM)
     &       -DTT(NTB+1,NP  ,NR  )    *DIVDTT*DTTM
     &       -FTH(NTB+1,NP  ,NR  )*VTB*DIVFTH*DTPM
      ENDIF
      IF(NP.EQ.NPMAX) THEN
         DL(NM)=DL(NM)
     &       +DTP(NTH+1,NP  ,NR  )*WTP*DIVDTP*DTPP
     &       -DTP(NTH  ,NP  ,NR  )*VTM*DIVDTP*DTPM
         IF(NTB.NE.0) THEN
            DL(NM)=DL(NM)
     &       +DTP(NTB+1,NP  ,NR  )*VTB*DIVDTP*DTPM
         ENDIF
      ENDIF
C      IF(NR.EQ.1) THEN
C         SP(NTH,NP,NR)=SP(NTH,NP,NR)+FS1(NTH,NP)
C     &         *(DRR(NTH  ,NP  ,NR  )    *DIVDRR*DRRM
C     &          +FRR(NTH  ,NP  ,NR  )*WRM*DIVFRR*DRRM)
C      ENDIF
C      IF(NR.EQ.NRMAX) THEN
C         SP(NTH,NP,NR)=SP(NTH,NP,NR)+FS2(NTH,NP)
C     &         *(DRR(NTH  ,NP  ,NR+1)    *DIVDRR*DRRP
C     &          -FRR(NTH  ,NP  ,NR+1)*VRP*DIVFRR*DRRP)
C      ENDIF
      RETURN
      END
