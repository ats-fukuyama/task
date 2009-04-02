C     $Id$
C
C ***********************************************************
C
C                    SELECTION OF GRAPH
C
C ***********************************************************
C
      SUBROUTINE FPGRAF
C
      INCLUDE 'fpcomm.inc'
      SAVE NSA,NSB
      DATA NSA,NSB/1,0/
C
      DIMENSION TEMP(NTHM,NPM,NRM)
      CHARACTER KID*5,KID1*1,KID2*3
      CHARACTER(LEN=80):: STRING1
C
    1 WRITE(6,*)'INPUT GRAPH TYPE :'
      WRITE(6,*)'    : F/FX/FR/FS1/FS2 1/2, Nn:NSA=n, Nnm:NSA=n,NSB=m,'
      WRITE(6,*)'    : D/DC/DW PP/PT/TP/TT/RR, F/FC/FE P/T/R'
      WRITE(6,*)'    : PD/PDC/PDW PP/PT/TP/TT/RR, PF/PFC/FE P/T/R'
      WRITE(6,*)'    : R/T N/I/W/PC/PW/PE/T/Q/E,  Gn,  X:exit'
      READ(5,'(A5)',ERR=1,END=9000) KID
      KID1=KID(1:1)
      CALL GUCPTL(KID1)

      IF(KID1.NE.'P') THEN
         KID2=KID(2:4)
         CALL GUCPTL(KID2(1:1))
         CALL GUCPTL(KID2(2:2))
         CALL GUCPTL(KID2(3:3))
C
      IF (KID1.EQ.'F') THEN
         IF(KID2.EQ.'1  ') THEN
            CALL FPGRAPA('F1',FNS,NS_NSA(NSA))
         ELSE IF(KID2.EQ.'R1 ') THEN
            CALL FPGRAPRA('FR1',FNS,NS_NSA(NSA))
         ELSE IF(KID2.EQ.'2  ') THEN
            CALL FPGRACA('F2',FNS,4,NS_NSA(NSA))
         ELSE IF(KID2.EQ.'X2  ') THEN
            DO NR=1,NRMAX
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               NS=NS_NSA(NSA)
               TEMP(NTH,NP,NR)=FNS(NTH,NP,NR,NS)
     &                        -FNS(NTHMAX+1-NTH,NP,NR,NS)
            ENDDO
            ENDDO
            ENDDO
            WRITE(STRING1,'(A,A1,I2,A1)') 'FX2','(',NSA,')'
            CALL FPGRAC(TRIM(STRING1),TEMP,4)
         ELSE IF(KID2.EQ.'Y2  ') THEN
            DO NR=1,NRMAX
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               NS=NS_NSA(NSA)
               TEMP(NTH,NP,NR)=FNS(NTH,NP,NR,NS)
     &                        -FNS(NTHMAX+1-NTH,NP,NR,NS)
            ENDDO
            ENDDO
            ENDDO
            WRITE(STRING1,'(A,A1,I2,A1)') 'FY2','(',NSA,')'
            CALL FPGRAC(TRIM(STRING1),TEMP,0)
         ELSE IF(KID2.EQ.'S11') THEN
            CALL FPGRAPA('FS11',FS1,NSA)
         ELSE IF(KID2.EQ.'S12') THEN
            CALL FPGRACA('FS12',FS1,4,NSA)
         ELSE IF(KID2.EQ.'S21') THEN
            CALL FPGRAPA('FS21',FS2,NSA)
         ELSE IF(KID2.EQ.'S22') THEN
            CALL FPGRACA('FS22',FS2,4,NSA)
         ELSE IF(KID2.EQ.'P  ') THEN
            CALL FPGRACA('FP  ',FPP,1,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRACA('FT  ',FTH,2,NSA)
         ELSE IF(KID2.EQ.'R  ') THEN
            CALL FPGRACA('FR  ',FRR,0,NSA)
         ELSE IF(KID2.EQ.'PP ') THEN
            CALL FPGRACA('FPP ',FPP,1,NSA)
         ELSE IF(KID2.EQ.'TH ') THEN
            CALL FPGRACA('FTH ',FTH,2,NSA)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL FPGRACA('FRR ',FRR,0,NSA)
         ELSE IF(KID2.EQ.'CP ') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACA('FCP ',FCPP,1,NSA)
            ELSE
               CALL FPGRACAB('FCP ',FCPP2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CT ') THEN
            IF(NSB.eq.0) THEN
               CALL FPGRACA('FCT ',FCTH,2,NSA)
            ELSE
               CALL FPGRACAB('FCT ',FCTH2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'EP ') THEN
            CALL FPGRACA('FEP ',FEPP,1,NSA)
         ELSE IF(KID2.EQ.'ET ') THEN
            CALL FPGRACA('FET ',FETH,2,NSA)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'R') THEN
         IF(KID2.EQ.'N  ') THEN
            CALL FPGRARA('RN  ',RNT,NSA)
         ELSE IF(KID2.EQ.'I  ') THEN
            CALL FPGRARA('RI  ',RJT,NSA)
         ELSE IF(KID2.EQ.'W  ') THEN
            CALL FPGRARA('RW  ',RWT,NSA)
         ELSE IF(KID2.EQ.'PC ') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRARA('RPC ',RPCT,NSA)
            ELSE
               CALL FPGRARAB('RPC ',RPCT2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'PW ') THEN
            CALL FPGRARA('RPW ',RPWT,NSA)
         ELSE IF(KID2.EQ.'PE ') THEN
            CALL FPGRARA('RPE ',RPET,NSA)
         ELSE IF(KID2.EQ.'LH ') THEN
            CALL FPGRARA('RLH ',RLHT,NSA)
         ELSE IF(KID2.EQ.'FW ') THEN
            CALL FPGRARA('RFW ',RFWT,NSA)
         ELSE IF(KID2.EQ.'EC ') THEN
            CALL FPGRARA('REC ',RECT,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRARA('RT  ',RTT,NSA)
         ELSE IF(KID2.EQ.'Q  ') THEN
            CALL FPGRAR('RQ  ',RQT)
         ELSE IF(KID2.EQ.'E  ') THEN
            CALL FPGRAR('RE  ',RET)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'T') THEN
         IF(KID2.EQ.'N  ') THEN
            CALL FPGRATA('TN  ',PNT,NSA)
         ELSE IF(KID2.EQ.'I  ') THEN
            CALL FPGRATA('TI  ',PIT,NSA)
         ELSE IF(KID2.EQ.'W  ') THEN
            CALL FPGRATA('TW  ',PWT,NSA)
         ELSE IF(KID2.EQ.'PC ') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRATA('TPC ',PPCT,NSA)
            ELSE
               CALL FPGRATAB('TPC ',PPCT2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'PW ') THEN
            CALL FPGRATA('TPW ',PPWT,NSA)
         ELSE IF(KID2.EQ.'PE ') THEN
            CALL FPGRATA('TPE ',PPET,NSA)
         ELSE IF(KID2.EQ.'LH ') THEN
            CALL FPGRATA('TLH ',PLHT,NSA)
         ELSE IF(KID2.EQ.'FW ') THEN
            CALL FPGRATA('TFW ',PFWT,NSA)
         ELSE IF(KID2.EQ.'EC ') THEN
            CALL FPGRATA('TEC ',PECT,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRATA('TT  ',PTT,NSA)
         ELSE IF(KID2.EQ.'Q  ') THEN
            CALL FPGRAT('TQ  ',PQT)
         ELSE IF(KID2.EQ.'E  ') THEN
            CALL FPGRAT('TE  ',PET)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'D') THEN
         IF(KID2.EQ.'PP ') THEN
            CALL FPGRACA('DPP ',DPP ,1,NSA)
         ELSE IF(KID2.EQ.'PT ') THEN
            CALL FPGRACA('DPT ',DPT ,1,NSA)
         ELSE IF(KID2.EQ.'TP ') THEN
            CALL FPGRACA('DTP ',DTP ,2,NSA)
         ELSE IF(KID2.EQ.'TT ') THEN
            CALL FPGRACA('DTT ',DTT ,2,NSA)
         ELSE IF(KID2.EQ.'CPP') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACA('DCPP',DCPP,1,NSA)
            ELSE
               CALL FPGRACAB('DCPP',DCPP2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CPT') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACA('DCPT',DCPT,1,NSA)
            ELSE
               CALL FPGRACAB('DCPT',DCPT2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CTP') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACA('DCTP',DCTP,2,NSA)
            ELSE
               CALL FPGRACAB('DCTP',DCTP2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CTT') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACA('DCTT',DCTT,2,NSA)
            ELSE
               CALL FPGRACAB('DCTT',DCTT2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'WPP') THEN
            CALL FPGRACA('DWPP',DWPP,1,NSA)
         ELSE IF(KID2.EQ.'WPT') THEN
            CALL FPGRACA('DWPT',DWPT,1,NSA)
         ELSE IF(KID2.EQ.'WTP') THEN
            CALL FPGRACA('DWTP',DWTP,2,NSA)
         ELSE IF(KID2.EQ.'WTT') THEN
            CALL FPGRACA('DWTT',DWTT,2,NSA)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL FPGRACA('DRR ',DRR ,0,NSA)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'N') THEN
         READ(KID2,*) NSA1
         IF(NSA1.GE.10) THEN
            NSA=NSA1/10
            NSB=MOD(NSA1,10)
         ELSE
            NSA=NSA1
            NSB=0
         ENDIF
         WRITE(6,'(A,2I3)') '# NSA and NSB are changed to',NSA,NSB
      ELSE IF (KID1.EQ.'W') THEN
         IF(KID2.EQ.'P  ') THEN
            CALL FPGRACA('WP  ',WEIGHP,1,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRACA('WT  ',WEIGHT,2,NSA)
         ELSE IF(KID2.EQ.'R  ') THEN
            CALL FPGRACA('WR  ',WEIGHR,2,NSA)
         ENDIF
      ELSE IF (KID1.EQ.'G') THEN
         READ(KID2,*,ERR=1,END=1) NGRAPH
         WRITE(6,'(A,I3)') '# NGRAPH is changed to',NGRAPH
      ELSE IF (KID1.EQ.'X') THEN
         GO TO 9000
      ELSE IF (KID1.EQ.'Q') THEN
         GO TO 9000
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID1'
      END IF
C
      ELSE
         KID1=KID(2:2)
         CALL GUCPTL(KID1)
         KID2=KID(3:5)
         CALL GUCPTL(KID2(1:1))
         CALL GUCPTL(KID2(2:2))
         CALL GUCPTL(KID2(3:3))
C
      IF (KID1.EQ.'F') THEN
         IF(KID2.EQ.'P  ') THEN
            CALL FPGRACPA('FP  ',FPP,1,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRACPA('FT  ',FTH,2,NSA)
         ELSE IF(KID2.EQ.'R  ') THEN
            CALL FPGRACPA('FR  ',FRR,0,NSA)
         ELSE IF(KID2.EQ.'PP ') THEN
            CALL FPGRACPA('FPP ',FPP,1,NSA)
         ELSE IF(KID2.EQ.'TH ') THEN
            CALL FPGRACPA('FTH ',FTH,2,NSA)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL FPGRACPA('FRR ',FRR,0,NSA)
         ELSE IF(KID2.EQ.'CP ') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACPA('FCP ',FCPP,1,NSA)
            ELSE
               CALL FPGRACPAB('FCP ',FCPP2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CT ') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACPA('FCT ',FCTH,2,NSA)
            ELSE
               CALL FPGRACPAB('FCT ',FCTH2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'EP ') THEN
            CALL FPGRACPA('FEP ',FEPP,1,NSA)
         ELSE IF(KID2.EQ.'ET ') THEN
            CALL FPGRACPA('FET ',FETH,2,NSA)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID3'
         ENDIF
      ELSE IF (KID1.EQ.'D') THEN
         IF(KID2.EQ.'PP ') THEN
            CALL FPGRACPA('DPP ',DPP ,1,NSA)
         ELSE IF(KID2.EQ.'PT ') THEN
            CALL FPGRACPA('DPT ',DPT ,1,NSA)
         ELSE IF(KID2.EQ.'TP ') THEN
            CALL FPGRACPA('DTP ',DTP ,2,NSA)
         ELSE IF(KID2.EQ.'TT ') THEN
            CALL FPGRACPA('DTT ',DTT ,2,NSA)
         ELSE IF(KID2.EQ.'CPP') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACPA('DCPP',DCPP,1,NSA)
            ELSE
               CALL FPGRACPAB('DCPP',DCPP2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CPT') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACPA('DCPT',DCPT,1,NSA)
            ELSE
               CALL FPGRACPAB('DCPT',DCPT2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CTP') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACPA('DCTP',DCTP,2,NSA)
            ELSE
               CALL FPGRACPAB('DCTP',DCTP2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CTT') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACPA('DCTT',DCTT,2,NSA)
            ELSE
               CALL FPGRACPAB('DCTT',DCTT2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'WPP') THEN
            CALL FPGRACPA('DWPP',DWPP,1,NSA)
         ELSE IF(KID2.EQ.'WPT') THEN
            CALL FPGRACPA('DWPT',DWPT,1,NSA)
         ELSE IF(KID2.EQ.'WTP') THEN
            CALL FPGRACPA('DWTP',DWTP,2,NSA)
         ELSE IF(KID2.EQ.'WTT') THEN
            CALL FPGRACPA('DWTT',DWTT,2,NSA)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL FPGRACPA('DRR ',DRR ,0,NSA)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID3'
         ENDIF
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID2'
      END IF
      END IF
C
      GO TO 1
C
 9000 RETURN
      END
C
C ***********************************************************
C
C                        R-GRAPHIC
C
C ***********************************************************
C
      SUBROUTINE FPGRARA(STRING,FRA,NSA)
      INCLUDE 'fpcomm.inc'
      DIMENSION FRA(NRM,NSAM,NTG1M),TEMP(NRM,NTG1M)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      DO NT1=1,NTG1
         DO NR=1,NRMAX
            TEMP(NR,NT1)=FRA(NR,NSA,NT1)
         ENDDO
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRAR(TRIM(STRING1),TEMP)
      RETURN
      END
C
      SUBROUTINE FPGRARAB(STRING,FRA,NSA,NSB)
      INCLUDE 'fpcomm.inc'
      DIMENSION FRAB(NRM,NSBM,NSAM,NTG1M),TEMP(NRM,NTG1M)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      DO NT1=1,NTG1
         DO NR=1,NRMAX
            TEMP(NR,NT1)=FRAB(NR,NSB,NSA,NT1)
         ENDDO
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1,I2,A1)') STRING,'(',NSA,',',NSB,')'
      CALL FPGRAR(TRIM(STRING1),TEMP)
      RETURN
      END
C
      SUBROUTINE FPGRAR(STRING,FR)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION FR(NRM,NTG1M),GX(NRM),GY(NRM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
C
      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTR(STRING,FR)
         RETURN
      ENDIF

      CALL PAGES
      CALL SETLNW(0.03)
      CALL SETCHS(0.3,0.)
      CALL SETFNT(32)

      DO NR=1,NRMAX
          GX(NR)=GUCLIP(RM(NR))
      ENDDO
      CALL GMNMX1(GX,1,NRMAX,1,GXMIN,GXMAX)
      IF(RGMIN.NE.0.0) GXMIN=GUCLIP(RGMIN)
      IF(RGMAX.NE.1.0) GXMIN=GUCLIP(RGMAX)

      GYMAX0=-1.E30
      GYMIN0= 1.E30
      DO 230 NT=1,NTG1
        DO 240 NR=1,NRMAX
          GX(NR)=GUCLIP(RM(NR))
          GY(NR)=GUCLIP(FR(NR,NT))
  240   CONTINUE
        CALL GMNMX1(GX,1,NRMAX,1,GXMIN,GXMAX)
        CALL GMNMX1(GY,1,NRMAX,1,GYMIN,GYMAX)
        GYMAX0=MAX(GYMAX0,GYMAX)
        GYMIN0=MIN(GYMIN0,GYMIN)
  230 CONTINUE
      IF(GYMIN0 .GE. 0.) GYMIN0=0.
      IF(GYMAX0 .LE. 0.) GYMAX0=0.

      CALL GQSCAL(GYMIN0,GYMAX0,GYMIN1,GYMAX1,GYSTEP)
      CALL GQSCAL(GXMIN,GXMAX,GXMIN1,GXMAX1,GXSTEP)
      CALL GDEFIN(3.0,23.0,2.0,17.0,GXMIN1,GXMAX1,GYMIN1,GYMAX1)
      CALL GFRAME
      GXORG=(INT(GXMIN1/(2*GXSTEP)+1))*2*GXSTEP
      CALL GSCALE(GXORG,GXSTEP,0.0,GYSTEP,0.1,9)    
      CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
      CALL GVALUE(0.,0.0,0.0,2*GYSTEP,NGULEN(2*GYSTEP))
      DO 250 NT=1,NTG1
        DO 260 NR=1,NRMAX
          GY(NR)=GUCLIP(FR(NR,NT))
          GX(NR)=GUCLIP(RM(NR))
  260   CONTINUE
        CALL SETLIN(0,2,7-MOD(NT-1,5))
        CALL GPLOTP(GX,GY,1,NRMAX,1,0,0,0)
  250 CONTINUE
      CALL SETLIN(0,2,7)
      CALL MOVE(1.0,17.5)
      CALL TEXT(STRING,LEN(STRING))
      CALL MOVE(24.0,1.0)
      CALL TEXT('r',1)
      CALL PAGEE
C
      RETURN
      END
C
C ***********************************************************
C
C                        T-GRAPHIC
C
C ***********************************************************
C
      SUBROUTINE FPGRATA(STRING,FTA,NSA)
      INCLUDE 'fpcomm.inc'
      DIMENSION FTA(NSAM,NTG2M),TEMP(NTG2M)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      DO NT2=1,NTG2
         TEMP(NT2)=FTA(NSA,NT2)
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRAT(TRIM(STRING1),TEMP)
      RETURN
      END
C
      SUBROUTINE FPGRATAB(STRING,FTA,NSA,NSB)
      INCLUDE 'fpcomm.inc'
      DIMENSION FTAB(NSBM,NSAM,NTG2M),TEMP(NTG2M)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      DO NT2=1,NTG2
         TEMP(NT2)=FTAB(NSB,NSA,NT2)
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1,I2,A1)') STRING,'(',NSA,',',NSB,')'
      CALL FPGRAT(TRIM(STRING1),TEMP)
      RETURN
      END
C
      SUBROUTINE FPGRAT(STRING,FT)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION  FT(NTG2M),GX(NTG2M),GY(NTG2M)
      CHARACTER(LEN=*),INTENT(IN):: STRING

      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTT(STRING,FT)
         RETURN
      ENDIF

      CALL PAGES
      CALL SETLNW(0.03)
      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
C
      DO 10 N=1,NTG2
        GX(N)=GUCLIP(PTG(N))
        GY(N)=GUCLIP(FT(N))
   10 CONTINUE
      CALL GMNMX1(GY,1,NTG2,1,GYMIN,GYMAX)
      IF(GYMIN.GT.0.) GYMIN=0.0
      IF(GYMAX.LT.0.) GYMAX=0.0
      CALL GQSCAL(GYMIN,GYMAX,GYMIN1,GYMAX1,GYSTEP)
      CALL GQSCAL(0.,GX(NTG2),GXMIN,GXMAX,GXSTEP)
      CALL GDEFIN(3.,23.,2.,17.,0.,GX(NTG2),GYMIN1,GYMAX1)
      CALL SETLNW(0.01)
      CALL GSCALE(0.,GXSTEP,0.,GYSTEP,1.0,0)
      CALL SETLNW(0.03)
      CALL GFRAME
      CALL GVALUE(0.,GXSTEP*2,0.,0.,NGULEN(2*GXSTEP))
      CALL GVALUE(0.,0.,0.,GYSTEP*2,NGULEN(2*GYSTEP))
      CALL GPLOTP(GX,GY,1,NTG2,1,0,0,0)
      CALL MOVE(1.0,17.5)
      CALL TEXT(STRING,LEN(STRING))
      CALL MOVE(24.0,1.0)
      CALL TEXT('t',1)
      CALL PAGEE
      RETURN
      END
C
C ***********************************************************
C
C                        P-GRAPHIC
C
C ***********************************************************
C
      SUBROUTINE FPGRAPA(STRING,FGA,NS)
      INCLUDE 'fpcomm.inc'
      DIMENSION FGA(NTHM,NPM,NRM,*),TEMP(NTHM,NPM,NRM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      DO NR=1,NRMAX
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=FGA(NTH,NP,NR,NS)
            ENDDO
         ENDDO
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NS,')'
      CALL FPGRAP(TRIM(STRING1),TEMP)
      RETURN
      END
C
      SUBROUTINE FPGRAPAB(STRING,FGAB,NSA,NSB)
      INCLUDE 'fpcomm.inc'
      DIMENSION FGAB(NTHM,NPM,NRM,NSBM,NSAM),TEMP(NTHM,NPM,NRM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      DO NR=1,NRMAX
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=FGAB(NTH,NP,NR,NSB,NSA)
            ENDDO
         ENDDO
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1,I2,A1)') STRING,'(',NSA,',',NSB,')'
      CALL FPGRAP(TRIM(STRING1),TEMP)
      RETURN
      END
C
      SUBROUTINE FPGRAP(STRING,FG)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION FG(NTHM,NPM,NRM),GX(NPM),GY(NPM,NRM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1

      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTP(STRING,FG)
         RETURN
      ENDIF

    1 WRITE(6,'(A,I4,A)') '# INPUT NTH (1..',NTHMAX,' or 0) :'
      READ(5,*,ERR=1,END=9000) NTH
      IF(NTH.LT.1) GOTO 9000
      IF(NTH.GT.NTHMAX) GOTO 1
      WRITE(STRING1,'(A,A,I4)') STRING,' : NTH=',NTH

      CALL PAGES
      CALL SETLNW(0.03)
      CALL SETCHS(0.3,0.)
      CALL SETFNT(32)
      DO 100 NR=1,NRMAX
      DO 100 NP=1,NPMAX
         IF(FG(NTH,NP,NR).LT.1.D-14) THEN
            GY(NP,NR)=-14.0
         ELSE
            GY(NP,NR)=GUCLIP(LOG10(FG(NTH,NP,NR)))
         ENDIF
  100 CONTINUE
      DO NP=1,NPMAX
         GX(NP)=GUCLIP(PM(NP)**2)
      ENDDO
      GXMAX=GUCLIP(PMAX**2)
      IF(PGMAX.NE.0.0) GXMAX=GUCLIP(PGMAX**2)

      CALL GMNMX2(GY,NPM,1,NPMAX,1,1,NRMAX,1,GYMIN,GYMAX)
      CALL GQSCAL(GYMIN,GYMAX,GYMINP,GYMAXP,GYSTEP)
      CALL GQSCAL(0.0,GXMAX,GXMINP,GXMAXP,GXSTEP)
C
      CALL GDEFIN(3.0,23.0,2.0,17.0,0.0,GXMAX,-8.0,0.0)
      CALL GFRAME
      CALL SETLNW(0.01)
      CALL GSCALE(0.,GXSTEP*2,0.0,0.0,1.0,0)
      CALL SETLNW(0.03)
      CALL GSCALE(0.,GXSTEP,0.0,0.0,0.2,9)
      CALL GSCALL(0.,0,0.0,1,1.0,0)
      CALL GSCALL(0.,0,0.0,2,0.2,9)
      CALL GFRAME
      CALL GVALUE(0.,GXSTEP*2,0.0,0.0,NGULEN(2*GXSTEP))
      CALL GVALUL(0.,0.0,0.0,1,0)
      DO 250 NR=1,NRMAX
        CALL SETLIN(0,0,7-MOD(NR-1,5))
        CALL GPLOTP(GX,GY(1,NR),1,NPMAX,1,0,0,0)
  250 CONTINUE
      CALL SETLIN(0,0,7)
      CALL MOVE(1.0,17.5)
      CALL TEXT(STRING1,LEN(STRING1))
      CALL MOVE(24.0,1.0)
      CALL TEXT('p**2',4)
      CALL PAGEE
C
 9000 RETURN
      END
C
C ***********************************************************
C
C                        PR-GRAPHIC
C
C ***********************************************************
C
      SUBROUTINE FPGRAPRA(STRING,FGA,NS)
      INCLUDE 'fpcomm.inc'
      DIMENSION FGA(NTHM,NPM,NRM,*),TEMP(NTHM,NPM,NRM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      DO NR=1,NRMAX
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=FGA(NTH,NP,NR,NS)
            ENDDO
         ENDDO
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NS,')'
      CALL FPGRAPR(TRIM(STRING1),TEMP)
      RETURN
      END
C
      SUBROUTINE FPGRAPRAB(STRING,FGAB,NSA,NSB)
      INCLUDE 'fpcomm.inc'
      DIMENSION FGAB(NTHM,NPM,NRM,NSBM,NSAM),TEMP(NTHM,NPM,NRM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      DO NR=1,NRMAX
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=FGAB(NTH,NP,NR,NSB,NSA)
            ENDDO
         ENDDO
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1,I2,A1)') STRING,'(',NSA,',',NSB,')'
      CALL FPGRAPR(TRIM(STRING1),TEMP)
      RETURN
      END
C
      SUBROUTINE FPGRAPR(STRING,FG)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION FG(NTHM,NPM,NRM),GX(NPM),GY(NRM),GZ(NPM,NRM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1

      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTP(STRING,FG)
         RETURN
      ENDIF

    1 WRITE(6,'(A,I4,A)') '# INPUT NTH (1..',NTHMAX,' or 0) :'
      READ(5,*,ERR=1,END=9000) NTH
      IF(NTH.LT.1) GOTO 9000
      IF(NTH.GT.NTHMAX) GOTO 1
      WRITE(STRING1,'(A,A,I4)') STRING,' : NTH=',NTH

      NPGMAX=1
      DO NP=1,NPMAX
         IF(PGMAX.NE.0.0) THEN
            IF(PM(NP).GT.PGMAX) GOTO 10
         ENDIF
         NPGMAX=NP
         GX(NP)=GUCLIP(PM(NP)**2)
      ENDDO
   10 CONTINUE
      DO NR=1,NRMAX
         GY(NR)=GUCLIP(RM(NR))
      ENDDO
      DO NR=1,NRMAX
      DO NP=1,NPMAX
         IF(FG(NTH,NP,NR).LT.1.D-14) THEN
            GZ(NP,NR)=-14.0
         ELSE
            GZ(NP,NR)=GUCLIP(LOG10(FG(NTH,NP,NR)))
         ENDIF
      ENDDO
      ENDDO

      GX1=3.0
      GX2=20.0
      GY1=2.0
      GY2=17.0

      CALL PAGES
      CALL FPGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,NPM,NPGMAX,NRMAX,
     &     GUCLIP(PGMAX**2),GUCLIP(RGMIN),GUCLIP(RGMAX),TRIM(STRING1))
      CALL PAGEE

 9000 RETURN
      END

!     ***********************************************************

!           SUBPROGRAM FOR 2D PROFILE

!     ***********************************************************

      SUBROUTINE FPGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,NXM,NXMAX,NYMAX,
     &                  GXMAX1,GYMIN1,GYMAX1,STR)

      IMPLICIT NONE
      REAL(4),        INTENT(IN):: GX1,GX2,GY1,GY2,GXMAX1,GYMIN1,GYMAX1
      INTEGER(4),     INTENT(IN):: NXM,NXMAX,NYMAX
      REAL(4),DIMENSION(NXMAX),      INTENT(IN):: GX
      REAL(4),DIMENSION(NYMAX),      INTENT(IN):: GY
      REAL(4),DIMENSION(NXM,NYMAX),  INTENT(IN):: GZ
      CHARACTER(LEN=*),             INTENT(IN):: STR
      INTEGER(4) :: I, NGULEN
      REAL(4)    :: GOX, GOY, GOZ, GPHI, GRADIUS, GSTEPX, GSTEPY, 
     &              GSTEPZ, GSXMAX, GSXMIN, GSYMAX, GSYMIN, GSZMAX,
     &              GSZMIN, GTHETA, GXL, GXMAX, GXMIN, GXORG, GYL, 
     &              GYMAX, GYMIN, GYORG, GZL, GZMAX, GZMIN, GZVAL
      EXTERNAL R2W2B,R2Y2W,W2G2B,WHITE

      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.035)
      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLIN(-1,-1,7)

      CALL MOVE(GX1,GY2+0.2)
      CALL TEXT(STR,LEN(STR))

      CALL GMNMX2(GZ,NXM,1,NXMAX,1,1,NYMAX,1,GZMIN,GZMAX)
      GZVAL=0.5*(ABS(GZMIN)+ABS(GZMAX))
      IF((GZVAL.LE.1.D-6).OR.(ABS(GZMAX-GZMIN)/GZVAL.LT.1.D-6)) THEN
         WRITE(6,*) GZMIN,GZMAX
         CALL GUFLSH
         RETURN
      ENDIF

      IF(GZMIN.GE.0.0) THEN
         GZMIN=0.0
      ELSEIF(GZMAX.LE.0.0) THEN
         GZMAX=0.0
      ENDIF

      CALL GMNMX1(GX,1,NXMAX,1,GXMIN,GXMAX)
      CALL GMNMX1(GY,1,NYMAX,1,GYMIN,GYMAX)

      IF(GXMAX1.NE.0.0) GXMAX=GXMAX1
      IF(GYMIN1.NE.0.0) GYMIN=GYMIN1
      IF(GYMAX1.NE.1.0) GYMAX=GYMAX1

      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSTEPY)
      CALL GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GSTEPZ)
!      WRITE(6,*) GSZMIN,GSZMAX,GSTEPZ
!      CALL GUFLSH

      IF(GXMIN*GXMAX.LT.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GSXMIN
      ENDIF
      IF(GYMIN*GYMAX.LT.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=GSYMIN
      ENDIF
      IF(GZMIN*GZMAX.LT.0.0) THEN
         GSZMAX=MAX(ABS(GSZMIN),ABS(GSZMAX))
         GSZMIN=-GSZMAX
      ENDIF

      GXL=10.0*1.5
      GYL=20.0*1.5
      GZL=10.0*1.5
      CALL GDEFIN3D(GX1,GX2,GY1,GY2,GXL,GYL,GZL)
      GPHI=-60.0
      GTHETA=65.0
      GRADIUS=100.0
      GOX=0.5*(GSXMIN+GSXMAX)
      GOY=0.5*(GYMIN+GYMAX)
      GOZ=0.5*(GSZMIN+GSZMAX)
      CALL GVIEW3D(GPHI,GTHETA,GRADIUS,1.0,1,GOX,GOY,GOZ)
      CALL GDATA3D3(GZ,GX,GY,NXM,NXMAX,NYMAX,
     &              GSXMIN,GSXMAX,GSYMIN,GSYMAX,GSZMIN,GSZMAX)
      CALL SETCHS(0.2,0.0)
      CALL SETLIN(0,0,7)

      CALL GSCALE3DX(GSXMIN,GSTEPX,0.3,0)
      CALL GSCALE3DY(GSYMIN,GSTEPY,0.3,0)
      CALL GSCALE3DZ(GSZMIN,GSTEPZ,0.3,0)
      CALL GVALUE3DX(GSXMIN,GSTEPX,1,NGULEN(GSTEPX))
      CALL GVALUE3DY(GSYMIN,GSTEPY,1,NGULEN(GSTEPY))
      CALL GVALUE3DZ(GSZMIN,GSTEPZ,2,NGULEN(GSTEPZ))

      IF(GZMIN*GZMAX.LT.0.0) THEN
         CALL CPLOT3D1(7,R2W2B)
      ELSEIF(GZMIN.LT.0.0) THEN
         CALL CPLOT3D1(7,W2G2B)
C         CALL CPLOT3D1(7,WHITE)
      ELSE
         CALL CPLOT3D1(7,R2Y2W)
      ENDIF

      CALL GAXIS3D(0)

      CALL SETLIN(0,0,7)
      RETURN
      END SUBROUTINE FPGR3D
C
C ***********************************************************
C
C                        C-GRAPHIC
C
C ***********************************************************
C
      SUBROUTINE FPGRACA(STRING,FGA,MODE,NSA)
      INCLUDE 'fpcomm.inc'
      DIMENSION FGA(*),TEMP(NMM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
C
      IF(MODE.EQ.0.OR.MODE.EQ.4) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*(NSA-1)+NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGA(NM1)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX+1
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*(NSA-1)+NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGA(NM1)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX+1
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*(NSA-1)+NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGA(NM1)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRAC(TRIM(STRING1),TEMP,MODE)
      RETURN
      END
C
      SUBROUTINE FPGRACAB(STRING,FGAB,MODE,NSA,NSB)
      INCLUDE 'fpcomm.inc'
      DIMENSION FGAB(*),TEMP(NMM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
C
      IF(MODE.EQ.0.OR.MODE.EQ.4) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*NSBM*(NSA-1)+NPM*NTHM*NRM*(NSB-1)
     &         +NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NM1)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX+1
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*NSBM*(NSA-1)+NPM*NTHM*NRM*(NSB-1)
     &         +NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NM1)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX+1
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*NSBM*(NSA-1)+NPM*NTHM*NRM*(NSB-1)
     &         +NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NM1)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      WRITE(STRING1,'(A,A1,I2,A1,I2,A1)') STRING,'(',NSA,',',NSB,')'
      CALL FPGRAC(TRIM(STRING1),TEMP,MODE)
      RETURN
      END
C
      SUBROUTINE FPGRAC(STRING,FG,MODE)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION FG(*)
      DIMENSION GF(NPM,NTHM),GP(NPM),GTH(NTHM),KA(8,NPM,NTHM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      INTEGER,PARAMETER:: NGLM=30
      REAL(4):: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER:: ILN(NGLM)
      REAL(4):: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX

      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTC(STRING,FG,MODE)
         RETURN
      ENDIF

    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.0) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,'(A,I4,A)') '# INPUT NR (1..',NRG,' or 0) :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF
      WRITE(STRING1,'(A,A,I4)') STRING,' : NR=',NR
C
      LMODE=MODE/4
C
      IF(MOD(MODE,2).EQ.0) THEN
         DO 10 NP=1,NPMAX
            GP(NP)=GUCLIP(PM(NP))
   10    CONTINUE
         NPG=NPMAX
      ELSE
         DO 20 NP=1,NPMAX+1
            GP(NP)=GUCLIP(PG(NP))
   20    CONTINUE
         NPG=NPMAX+1
      ENDIF
C
      IF(MOD(MODE/2,2).EQ.0) THEN
         DO 30 NTH=1,NTHMAX
            GTH(NTH)=GUCLIP(THM(NTH))
   30    CONTINUE
         NTHG=NTHMAX
      ELSE
         DO 40 NTH=1,NTHMAX+1
            GTH(NTH)=GUCLIP(THG(NTH))
   40    CONTINUE
         NTHG=NTHMAX+1
      ENDIF
C
      IF(MODE.EQ.0) THEN
         DO 50 NTH=1,NTHMAX
         DO 50 NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            GF(NP,NTH)=GUCLIP(FG(NM))
   50    CONTINUE
      ELSEIF(MODE.EQ.1) THEN
         DO 60 NTH=1,NTHMAX
         DO 60 NP=1,NPMAX+1
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            GF(NP,NTH)=GUCLIP(FG(NM))
   60    CONTINUE
      ELSEIF(MODE.EQ.2) THEN
         DO 70 NTH=1,NTHMAX+1
         DO 70 NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            GF(NP,NTH)=GUCLIP(FG(NM))
   70    CONTINUE
      ELSEIF(MODE.EQ.4) THEN
         DO 80 NTH=1,NTHMAX
         DO 80 NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            IF(FG(NM).LT.1.D-70) THEN
               GF(NP,NTH)=-70.0
            ELSE
               GF(NP,NTH)=GUCLIP(LOG10(ABS(FG(NM))))
            ENDIF
   80    CONTINUE
      ENDIF
      GPMAX=GUCLIP(PMAX)
C
      CALL PAGES
      CALL SETLNW(0.03)
      CALL SETCHS(.3,0.)
      CALL SETFNT(32)
C
      CALL GMNMX2(GF,NPM,1,NPG,1,1,NTHG,1,GFMIN,GFMAX)
      CALL GQSCAL(GFMIN,GFMAX,GFMIN1,GFMAX1,GFSTEP)
      CALL GQSCAL(0.0,GPMAX,GPMIN1,GPMAX1,GPSTEP)
C
      CALL GDEFIN(3.,23.,2.,12.,-GPMAX,GPMAX,0.,GPMAX)
      CALL GFRAME
      CALL SETLNW(0.01)
      CALL GSCALE(0.,GPSTEP,0.,GPSTEP,1.0,0)
      CALL SETLNW(0.03)
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,NGULEN(2*GPSTEP))
C
      IF(LMODE.EQ.0) THEN
         IF(GFMIN*GFMAX.GE.0.0) THEN
            IF(GFMIN.GE.0.0) THEN
               NGLMAX=INT((GFMAX1-GFMIN1)/(0.5*GFSTEP))
               DO NGL=1,NGLMAX
                  ZL(NGL)=GFMIN1+0.5*GFSTEP*(NGL-1)
                  RGB(1,NGL)=1.D0
                  RGB(2,NGL)=FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  ILN(NGL)=0
                  IF(MOD(NGL,5).EQ.0) THEN
                     WLN(NGL)=0.06
                  ELSE
                     WLN(NGL)=0.03
                  ENDIF
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,
     &                     ZL,RGB,ILN,WLN,NGLMAX,0,0)
c$$$               CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
c$$$     &                     GFMIN1,0.5*GFSTEP,30,0,KA)
            ELSE
               NGLMAX=INT((GFMAX1-GFMIN1)/(0.5*GFSTEP))
               DO NGL=1,NGLMAX
                  ZL(NGL)=GFMAX1-0.5*GFSTEP*(NGL-1)
                  RGB(1,NGL)=FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(2,NGL)=FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=1.D0
                  ILN(NGL)=0
                  IF(MOD(NGL,5).EQ.0) THEN
                     WLN(NGL)=0.06
                  ELSE
                     WLN(NGL)=0.03
                  ENDIF
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,
     &                     ZL,RGB,ILN,WLN,NGLMAX,0,0)
c$$$               CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
c$$$     &                     GFMIN1,0.5*GFSTEP,30,2,KA)
            ENDIF
         ELSE
               GFFMAX=MAX(ABS(GFMAX1),ABS(GFMIN1))
               NGLMAX=INT(GFFMAX/(0.5*GFSTEP))
               DO NGL=1,NGLMAX
                  ZL(NGL)=-0.25*GFSTEP-0.5*GFSTEP*(NGL-1)
                  RGB(1,NGL)=FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(2,NGL)=FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=1.D0
                  ILN(NGL)=0
                  IF(MOD(NGL,5).EQ.0) THEN
                     WLN(NGL)=0.06
                  ELSE
                     WLN(NGL)=0.03
                  ENDIF
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,
     &                     ZL,RGB,ILN,WLN,NGLMAX,0,0)
               DO NGL=1,NGLMAX
                  ZL(NGL)= 0.25*GFSTEP+0.5*GFSTEP*(NGL-1)
                  RGB(1,NGL)=1.D0
                  RGB(2,NGL)=FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  ILN(NGL)=0
                  IF(MOD(NGL,5).EQ.0) THEN
                     WLN(NGL)=0.06
                  ELSE
                     WLN(NGL)=0.03
                  ENDIF
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,
     &                     ZL,RGB,ILN,WLN,NGLMAX,0,0)
c$$$            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
c$$$     &                   0.25*GFSTEP, 0.5*GFSTEP,15,0,KA)
c$$$            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
c$$$     &                  -0.25*GFSTEP,-0.5*GFSTEP,15,2,KA)
         ENDIF
      ELSE
         DO 100 I=1,NGLINE
           GLIN=GFMAX-0.020*(I-1)**2
           CALL SETLIN(0,0,7-MOD(I-1,5))
           CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
     &                 GLIN,GFSTEP*100,1,0,KA)
  100    CONTINUE
      ENDIF
C
      CALL SETLIN(0,2,7)
      CALL MOVE(24.0,1.0)
      CALL TEXT('PPARA',5)
      CALL MOVE(1.0,13.5)
      CALL TEXT('PPERP',5)
C
      CALL MOVE(3.0,12.5)
      CALL TEXT(STRING1,LEN(STRING1))
      CALL MOVE(8.0,12.5)
      CALL TEXT('FMIN =',6)
      CALL NUMBR(GFMIN,'(1PE12.4)',12)
      CALL MOVE(13.0,12.5)
      CALL TEXT('FMAX =',6)
      CALL NUMBR(GFMAX,'(1PE12.4)',12)
      IF(LMODE.EQ.0) THEN
         CALL MOVE(18.0,12.5)
         CALL TEXT('STEP =',6)
         CALL NUMBR(0.5*GFSTEP,'(1PE12.4)',12)
      ENDIF
      CALL PAGEE
      IF(NRMAX.GT.1) GOTO 1
C
 9000 RETURN
      END
C
C ***********************************************************
C
C                        C-GRAPHIC muliplied p
C
C ***********************************************************
C
      SUBROUTINE FPGRACPA(STRING,FGA,MODE,NSA)
      INCLUDE 'fpcomm.inc'
      DIMENSION FGA(*),TEMP(NMM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
C
      IF(MODE.EQ.0.OR.MODE.EQ.4) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*(NSA-1)+NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGA(NM1)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX+1
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*(NSA-1)+NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGA(NM1)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX+1
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*(NSA-1)+NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGA(NM1)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRACP(TRIM(STRING1),TEMP,MODE)
      RETURN
      END
C
      SUBROUTINE FPGRACPAB(STRING,FGAB,MODE,NSA,NSB)
      INCLUDE 'fpcomm.inc'
      DIMENSION FGAB(*),TEMP(NMM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
C
      IF(MODE.EQ.0.OR.MODE.EQ.4) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*NSBM*(NSA-1)+NPM*NTHM*NRM*(NSB-1)
     &         +NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NM1)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX+1
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*NSBM*(NSA-1)+NPM*NTHM*NRM*(NSB-1)
     &         +NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NM1)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX+1
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*NSBM*(NSA-1)+NPM*NTHM*NRM*(NSB-1)
     &         +NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NM1)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRACP(TRIM(STRING1),TEMP,MODE)
      RETURN
      END
C
      SUBROUTINE FPGRACP(STRING,FG,MODE)
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION FG(*)
      DIMENSION GF(NPM,NTHM),GP(NPM),GTH(NTHM),KA(8,NPM,NTHM)
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      INTEGER,PARAMETER:: NGLM=30
      REAL(4):: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER:: ILN(NGLM)
      REAL(4):: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
C
      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTC(STRING,FG,MODE)
         RETURN
      ENDIF

    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.0) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,'(A,I4,A)') '# INPUT NR (1..',NRG,' or 0) :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF
      WRITE(STRING1,'(A,A,I4)') STRING,' : NR=',NR
C
      LMODE=MODE/4
C
      IF(MOD(MODE,2).EQ.0) THEN
         DO 10 NP=1,NPMAX
            GP(NP)=GUCLIP(PM(NP))
   10    CONTINUE
         NPG=NPMAX
      ELSE
         DO 20 NP=1,NPMAX+1
            GP(NP)=GUCLIP(PG(NP))
   20    CONTINUE
         NPG=NPMAX+1
      ENDIF
C
      IF(MOD(MODE/2,2).EQ.0) THEN
         DO 30 NTH=1,NTHMAX
            GTH(NTH)=GUCLIP(THM(NTH))
   30    CONTINUE
         NTHG=NTHMAX
      ELSE
         DO 40 NTH=1,NTHMAX+1
            GTH(NTH)=GUCLIP(THG(NTH))
   40    CONTINUE
         NTHG=NTHMAX+1
      ENDIF
C
      IF(MODE.EQ.0) THEN
         DO 50 NTH=1,NTHMAX
         DO 50 NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            GF(NP,NTH)=GUCLIP(FG(NM)*PM(NP))
   50    CONTINUE
      ELSEIF(MODE.EQ.1) THEN
         DO 60 NTH=1,NTHMAX
         DO 60 NP=1,NPMAX+1
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            GF(NP,NTH)=GUCLIP(FG(NM)*PG(NP))
   60    CONTINUE
      ELSEIF(MODE.EQ.2) THEN
         DO 70 NTH=1,NTHMAX+1
         DO 70 NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            GF(NP,NTH)=GUCLIP(FG(NM)*PM(NP))
   70    CONTINUE
      ELSEIF(MODE.EQ.4) THEN
         DO 80 NTH=1,NTHMAX
         DO 80 NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            IF(FG(NM)*PM(NP).LT.1.D-70) THEN
               GF(NP,NTH)=-70.0
            ELSE
               GF(NP,NTH)=GUCLIP(LOG10(ABS(FG(NM)*PM(NP))))
            ENDIF
   80    CONTINUE
      ENDIF
      GPMAX=GUCLIP(PMAX)
C
      CALL PAGES
      CALL SETLNW(0.03)
      CALL SETCHS(.3,0.)
      CALL SETFNT(32)
C
      CALL GMNMX2(GF,NPM,1,NPG,1,1,NTHG,1,GFMIN,GFMAX)
      CALL GQSCAL(GFMIN,GFMAX,GFMIN1,GFMAX1,GFSTEP)
      CALL GQSCAL(0.0,GPMAX,GPMIN1,GPMAX1,GPSTEP)
C
      CALL GDEFIN(3.,23.,2.,12.,-GPMAX,GPMAX,0.,GPMAX)
      CALL GFRAME
      CALL SETLNW(0.01)
      CALL GSCALE(0.,GPSTEP,0.,GPSTEP,1.0,0)
      CALL SETLNW(0.03)
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,NGULEN(2*GPSTEP))
C
      IF(LMODE.EQ.0) THEN
         IF(GFMIN*GFMAX.GE.0.0) THEN
            IF(GFMIN.GE.0.0) THEN
               NGLMAX=INT((GFMAX-GFMIN1)/(0.5*GFSTEP))
               DO NGL=1,NGLMAX
                  RGB(1,NGL)=1.D0
                  RGB(2,NGL)=FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  ILN(NGL)=0
                  IF(MOD(NGL,5).EQ.0) THEN
                     WLN(NGL)=0.06
                  ELSE
                     WLN(NGL)=0.03
                  ENDIF
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,
     &                     ZL,RGB,ILN,WLN,NGLMAX,0,0)
            ELSE
               CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
     &                     GFMIN1,0.5*GFSTEP,30,2,KA)
            ENDIF
         ELSE
            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
     &                   0.25*GFSTEP, 0.5*GFSTEP,15,0,KA)
            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
     &                  -0.25*GFSTEP,-0.5*GFSTEP,15,2,KA)
         ENDIF
      ELSE
         DO 100 I=1,NGLINE
           GLIN=GFMAX-0.020*(I-1)**2
           CALL SETLIN(0,0,7-MOD(I-1,5))
           CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
     &                 GLIN,GFSTEP*100,1,0,KA)
  100    CONTINUE
      ENDIF
C
      CALL SETLIN(0,2,7)
      CALL MOVE(24.0,1.0)
      CALL TEXT('PPARA',5)
      CALL MOVE(1.0,13.5)
      CALL TEXT('PPERP',5)
C
      CALL MOVE(3.0,12.5)
      CALL TEXT(STRING1,LEN(STRING1))
      CALL MOVE(8.0,12.5)
      CALL TEXT('FMIN =',6)
      CALL NUMBR(GFMIN,'(1PE12.4)',12)
      CALL MOVE(13.0,12.5)
      CALL TEXT('FMAX =',6)
      CALL NUMBR(GFMAX,'(1PE12.4)',12)
      IF(LMODE.EQ.0) THEN
         CALL MOVE(18.0,12.5)
         CALL TEXT('STEP =',6)
         CALL NUMBR(0.5*GFSTEP,'(1PE12.4)',12)
      ENDIF
      CALL PAGEE
      IF(NRMAX.GT.1) GOTO 1
C
 9000 RETURN
      END
C>>>>>>
C
C     ****** CONTOUR PLOT : R-THETA, VARIABLE STEP, PATTERN ******
C
      SUBROUTINE CONTG4X(Z,R,T,NXA,NXMAX,NYMAX,
     &                  ZL,RGB,ILN,WLN,NSTEP,ISPL)
C
      IMPLICIT LOGICAL(L)
      COMMON /GSCTR4/ RMAX,RT,TT,XT,YT
C
      EXTERNAL CONTV2X
      DIMENSION Z(NXA,NYMAX),R(NXMAX),T(NYMAX)
      DIMENSION ZL(NSTEP),RGB(3,NSTEP),ILN(NSTEP),WLN(NSTEP)
C
      RMAX=R(NXMAX)
C
      CALL CONTG0X(Z,R,T,NXA,NXMAX,NYMAX,
     &            ZL,RGB,ILN,WLN,NSTEP,ISPL,0,3,CONTV2X)
C
      RETURN
      END
C
C     ****** CONTOUR PLOT : COMMON SUB ******
C
      SUBROUTINE CONTG0X(Z,X,Y,NXA,NXMAX,NYMAX,
     &                  ZL,RGB,ILN,WLN,NSTEP,
     &                  ISPL,IPRD,INDX,SUBV)
C
      IMPLICIT LOGICAL(L)
      EXTERNAL SUBV
      DIMENSION Z(NXA,NYMAX),X(NXMAX),Y(NYMAX),KA(2,NXMAX*NYMAX)
      DIMENSION ZL(NSTEP),RGB(3,NSTEP),ILN(NSTEP),WLN(NSTEP)
      PARAMETER(NH=101)
      DIMENSION ZLS(NH),RGBS(3,NH),ILNS(NH),WLNS(NH)
C      PARAMETER (NFMAX=2000,NGMAX=4000)
      PARAMETER (NFMAX=200,NGMAX=400)
      DIMENSION XF(NFMAX),YF(NFMAX)
      DIMENSION XP(NFMAX),YP(NFMAX)
      DIMENSION XG(NGMAX),YG(NGMAX)
C
      IF(ISPL.GE.0) THEN
         CALL INQRGB(RS,GS,BS)
         CALL INQLNW(WS)
      ENDIF
C
      KMAX=NSTEP
      IF(KMAX.GT.NH) KMAX=NH
C
      IF(ZL(1).LE.ZL(KMAX)) THEN
         ZORG=ZL(1)
         DO 10 I=1,KMAX
            ZLS(I)=ZL(I)-ZORG
   10    CONTINUE
      ELSE
         ZORG=ZL(KMAX)
         DO 20 I=1,KMAX
            ZLS(I)=ZL(KMAX-I+1)-ZORG
   20    CONTINUE
      ENDIF
C
      IF(ISPL.GE.0) THEN
         IF(ZL(1).LE.ZL(KMAX)) THEN
            DO 30 I=1,KMAX
               RGBS(1,I)=RGB(1,I)
               RGBS(2,I)=RGB(2,I)
               RGBS(3,I)=RGB(3,I)
               ILNS(I)=ILN(I)
               WLNS(I)=WLN(I)
   30       CONTINUE
         ELSE
            DO 40 I=1,KMAX
               RGBS(1,I)=RGB(1,KMAX-I+1)
               RGBS(2,I)=RGB(2,KMAX-I+1)
               RGBS(3,I)=RGB(3,KMAX-I+1)
               ILNS(I)=ILN(KMAX-I+1)
               WLNS(I)=WLN(KMAX-I+1)
   40       CONTINUE
         ENDIF
      ELSE
         DO 50 I=1,KMAX
            ILNS(I)=ILN(1)
            WLNS(I)=WLN(1)
   50    CONTINUE
      ENDIF
C
      NXM=NXMAX-1
      NYM=NYMAX-1
      IF(IPRD.EQ.1.OR.IPRD.EQ.3) NXM=NXMAX
      IF(IPRD.EQ.2.OR.IPRD.EQ.3) NYM=NYMAX
C
      DO 70 NY=1,NYM
      DO 70 NX=1,NXM
         IE=NXM*(NY-1)+NX
         NXP=NX+1
         NYP=NY+1
         IF(NXP.GT.NXMAX) NXP=1
         IF(NYP.GT.NYMAX) NYP=1
         U1=Z(NX ,NY )-ZORG
         U2=Z(NXP,NY )-ZORG
         U3=Z(NXP,NYP)-ZORG
         U4=Z(NX ,NYP)-ZORG
         UMAX=MAX(U1,U2,U3,U4)
         UMIN=MIN(U1,U2,U3,U4)
         KA(1,IE)=1
         DO 60 K=1,KMAX
            IF(ZLS(K).LT.UMIN) KA(1,IE)=K+1
            IF(ZLS(K).LE.UMAX) KA(2,IE)=K
   60    CONTINUE
   70 CONTINUE
C
      DO 1000 K=1,KMAX
         U0=ZLS(K)
         IF(ISPL.GE.0) THEN
            CALL SETRGB(RGBS(1,K),RGBS(2,K),RGBS(3,K))
            CALL SETLNW(WLNS(K))
         ENDIF
      DO 1000 NY=1,NYM
      DO 1000 NX=1,NXM
         IE=NXM*(NY-1)+NX
         IF(KA(1,IE).LE.K.AND.KA(2,IE).GE.K) THEN
            N1X=NX
            N1Y=NY
            N2X=NX+1
            N2Y=NY
            N3X=NX+1
            N3Y=NY+1
            N4X=NX
            N4Y=NY+1
C     
            U1=Z(N1X,N1Y)-ZORG
            I2X=N2X
            I2Y=N2Y
            IF(I2X.GT.NXMAX) I2X=1
            IF(I2Y.GT.NYMAX) I2Y=1
            U2=Z(I2X,I2Y)-ZORG
            I3X=N3X
            I3Y=N3Y
            IF(I3X.GT.NXMAX) I3X=1
            IF(I3Y.GT.NYMAX) I3Y=1
            U3=Z(I3X,I3Y)-ZORG
            I4X=N4X
            I4Y=N4Y
            IF(I4X.GT.NXMAX) I4X=1
            IF(I4Y.GT.NYMAX) I4Y=1
            U4=Z(I4X,I4Y)-ZORG
C
            IF((U1.GT.U0.AND.U2.LE.U0).OR.
     &         (U1.LE.U0.AND.U2.GT.U0)) THEN
               NAX=N1X
               NAY=N1Y
               UA=U1
               NBX=N2X
               NBY=N2Y
               UB=U2
               MODE=1
               IF((U2.GT.U0.AND.U3.LE.U0).OR.
     &            (U2.LE.U0.AND.U3.GT.U0)) THEN
                  NSAX=N2X
                  NSAY=N2Y
                  USA=U2
                  NSBX=N3X
                  NSBY=N3Y
                  USB=U3
                  MODES=2
               ELSEIF((U3.GT.U0.AND.U4.LE.U0).OR.
     &                (U3.LE.U0.AND.U4.GT.U0)) THEN
                  NSAX=N3X
                  NSAY=N3Y
                  USA=U3
                  NSBX=N4X
                  NSBY=N4Y
                  USB=U4
                  MODES=3
               ELSEIF((U4.GT.U0.AND.U1.LE.U0).OR.
     &                (U4.LE.U0.AND.U1.GT.U0)) THEN
                  NSAX=N4X
                  NSAY=N4Y
                  USA=U4
                  NSBX=N1X
                  NSBY=N1Y
                  USB=U1
                  MODES=4
               ELSE
                  GOTO 1000
               ENDIF
            ELSEIF((U2.GT.U0.AND.U3.LE.U0).OR.
     &             (U2.LE.U0.AND.U3.GT.U0)) THEN
               NAX=N2X
               NAY=N2Y
               UA=U2
               NBX=N3X
               NBY=N3Y
               UB=U3
               MODE=2
               IF((U3.GT.U0.AND.U4.LE.U0).OR.
     &            (U3.LE.U0.AND.U4.GT.U0)) THEN
                  NSAX=N3X
                  NSAY=N3Y
                  USA=U3
                  NSBX=N4X
                  NSBY=N4Y
                  USB=U4
                  MODES=3
               ELSEIF((U4.GT.U0.AND.U1.LE.U0).OR.
     &                (U4.LE.U0.AND.U1.GT.U0)) THEN
                  NSAX=N4X
                  NSAY=N4Y
                  USA=U4
                  NSBX=N1X
                  NSBY=N1Y
                  USB=U1
                  MODES=4
               ELSE
                  GOTO 1000
               ENDIF
            ELSEIF((U3.GT.U0.AND.U4.LE.U0).OR.
     &             (U3.LE.U0.AND.U4.GT.U0)) THEN
               NAX=N3X
               NAY=N3Y
               UA=U3
               NBX=N4X
               NBY=N4Y
               UB=U4
               MODE=3
               IF((U4.GT.U0.AND.U1.LE.U0).OR.
     &            (U4.LE.U0.AND.U1.GT.U0)) THEN
                  NSAX=N4X
                  NSAY=N4Y
                  USA=U4
                  NSBX=N1X
                  NSBY=N1Y
                  USB=U1
                  MODES=4
               ELSE
                  GOTO 1000
               ENDIF
            ELSE
               GOTO 1000
            ENDIF
C
            IF(INDX.EQ.0.OR.INDX.EQ.2) THEN
               XA=X(1)*(NAX-1)+X(2)
               XB=X(1)*(NBX-1)+X(2)
            ELSE
               XA=X(NAX)
               XB=X(NBX)
            ENDIF
            IF(INDX.EQ.0.OR.INDX.EQ.1) THEN
               YA=Y(1)*(NAY-1)+Y(2)
               YB=Y(1)*(NBY-1)+Y(2)
            ELSE
               YA=Y(NAY)
               YB=Y(NBY)
            ENDIF
            J=1
            RT=(U0-UA)/(UB-UA)
            XF(J)=(XB-XA)*RT+XA
            YF(J)=(YB-YA)*RT+YA
C
            IEL=IE
            NXL=NX
            NYL=NY
            MODEL=MODE
            LINV=.FALSE.
            LEND=.FALSE.
C
            NXL1=NX
            NYL1=NY
            IF(MODEL.EQ.1) NYL1=NYL-1
            IF(MODEL.EQ.2) NXL1=NXL+1
            IF(MODEL.EQ.3) NYL1=NYL+1
            IF(MODEL.EQ.4) NXL1=NXL-1
C
            IF(NXL1.GE.1.AND.NXL1.LE.NXM.AND.
     &         NYL1.GE.1.AND.NYL1.LE.NYM) THEN
               IE1=NXM*(NYL1-1)+NXL1
            ELSE
               IE1=IE
            ENDIF
C
  200       IF(MODEL.EQ.1) NYL=NYL-1
            IF(MODEL.EQ.2) NXL=NXL+1
            IF(MODEL.EQ.3) NYL=NYL+1
            IF(MODEL.EQ.4) NXL=NXL-1
C
            IF(NXL.GE.1.AND.NXL.LE.NXM.AND.
     &         NYL.GE.1.AND.NYL.LE.NYM) THEN
               IEL=NXM*(NYL-1)+NXL
               N1X=NXL
               N1Y=NYL
               N2X=NXL+1
               N2Y=NYL
               N3X=NXL+1
               N3Y=NYL+1
               N4X=NXL
               N4Y=NYL+1
               U1=Z(N1X,N1Y)-ZORG
               I2X=N2X
               I2Y=N2Y
               IF(I2X.GT.NXMAX) I2X=1
               IF(I2Y.GT.NYMAX) I2Y=1
               U2=Z(I2X,I2Y)-ZORG
               I3X=N3X
               I3Y=N3Y
               IF(I3X.GT.NXMAX) I3X=1
               IF(I3Y.GT.NYMAX) I3Y=1
               U3=Z(I3X,I3Y)-ZORG
               I4X=N4X
               I4Y=N4Y
               IF(I4X.GT.NXMAX) I4X=1
               IF(I4Y.GT.NYMAX) I4Y=1
               U4=Z(I4X,I4Y)-ZORG
C
               IF(MODEL.EQ.1) THEN
                  IF((U4.GT.U0.AND.U1.LE.U0).OR.
     &               (U4.LE.U0.AND.U1.GT.U0)) THEN
                     NAX=N4X
                     NAY=N4Y
                     UA=U4
                     NBX=N1X
                     NBY=N1Y
                     UB=U1
                     MODEL=4
                  ELSEIF((U1.GT.U0.AND.U2.LE.U0).OR.
     &                   (U1.LE.U0.AND.U2.GT.U0)) THEN
                     NAX=N1X
                     NAY=N1Y
                     UA=U1
                     NBX=N2X
                     NBY=N2Y
                     UB=U2
                     MODEL=1
                  ELSEIF((U2.GT.U0.AND.U3.LE.U0).OR.
     &                   (U2.LE.U0.AND.U3.GT.U0)) THEN
                     NAX=N2X
                     NAY=N2Y
                     UA=U2
                     NBX=N3X
                     NBY=N3Y
                     UB=U3
                     MODEL=2
                  ELSE
                     LEND=.TRUE.
                  ENDIF
               ELSE IF(MODEL.EQ.2) THEN
                  IF((U1.GT.U0.AND.U2.LE.U0).OR.
     &               (U1.LE.U0.AND.U2.GT.U0)) THEN
                     NAX=N1X
                     NAY=N1Y
                     UA=U1
                     NBX=N2X
                     NBY=N2Y
                     UB=U2
                     MODEL=1
                  ELSEIF((U2.GT.U0.AND.U3.LE.U0).OR.
     &                   (U2.LE.U0.AND.U3.GT.U0)) THEN
                     NAX=N2X
                     NAY=N2Y
                     UA=U2
                     NBX=N3X
                     NBY=N3Y
                     UB=U3
                     MODEL=2
                  ELSEIF((U3.GT.U0.AND.U4.LE.U0).OR.
     &                   (U3.LE.U0.AND.U4.GT.U0)) THEN
                     NAX=N3X
                     NAY=N3Y
                     UA=U3
                     NBX=N4X
                     NBY=N4Y
                     UB=U4
                     MODEL=3
                  ELSE
                     LEND=.TRUE.
                  ENDIF
               ELSEIF(MODEL.EQ.3) THEN
                  IF((U2.GT.U0.AND.U3.LE.U0).OR.
     &               (U2.LE.U0.AND.U3.GT.U0)) THEN
                     NAX=N2X
                     NAY=N2Y
                     UA=U2
                     NBX=N3X
                     NBY=N3Y
                     UB=U3
                     MODEL=2
                  ELSEIF((U3.GT.U0.AND.U4.LE.U0).OR.
     &                   (U3.LE.U0.AND.U4.GT.U0)) THEN
                     NAX=N3X
                     NAY=N3Y
                     UA=U3
                     NBX=N4X
                     NBY=N4Y
                     UB=U4
                     MODEL=3
                  ELSEIF((U4.GT.U0.AND.U1.LE.U0).OR.
     &                   (U4.LE.U0.AND.U1.GT.U0)) THEN
                     NAX=N4X
                     NAY=N4Y
                     UA=U4
                     NBX=N1X
                     NBY=N1Y
                     UB=U1
                     MODEL=4
                  ELSE
                     LEND=.TRUE.
                  ENDIF
               ELSEIF(MODEL.EQ.4) THEN
                  IF((U3.GT.U0.AND.U4.LE.U0).OR.
     &                   (U3.LE.U0.AND.U4.GT.U0)) THEN
                     NAX=N3X
                     NAY=N3Y
                     UA=U3
                     NBX=N4X
                     NBY=N4Y
                     UB=U4
                     MODEL=3
                  ELSEIF((U4.GT.U0.AND.U1.LE.U0).OR.
     &                   (U4.LE.U0.AND.U1.GT.U0)) THEN
                     NAX=N4X
                     NAY=N4Y
                     UA=U4
                     NBX=N1X
                     NBY=N1Y
                     UB=U1
                     MODEL=4
                  ELSEIF((U1.GT.U0.AND.U2.LE.U0).OR.
     &               (U1.LE.U0.AND.U2.GT.U0)) THEN
                     NAX=N1X
                     NAY=N1Y
                     UA=U1
                     NBX=N2X
                     NBY=N2Y
                     UB=U2
                     MODEL=1
                  ELSE
                     LEND=.TRUE.
                  ENDIF
               ENDIF
            ELSE
               IF(LINV) THEN
                  LEND=.TRUE.
               ELSE
C
                  DO 300 NF=NFMAX,NFMAX-J+1,-1
                     XF(NF)=XF(NF-NFMAX+J)
                     YF(NF)=YF(NF-NFMAX+J)
  300             CONTINUE
                  J=NFMAX-J+1
C
                  IEL=IE
                  NAX=NSAX
                  NAY=NSAY
                  UA=USA
                  NBX=NSBX
                  NBY=NSBY
                  UB=USB
                  NXL=NX
                  NYL=NY
                  MODEL=MODES
                  LINV=.TRUE.
               ENDIF
            ENDIF
C
            IF(.NOT.LEND) THEN
C
               IF(INDX.EQ.0.OR.INDX.EQ.2) THEN
                  XA=X(1)*(NAX-1)+X(2)
                  XB=X(1)*(NBX-1)+X(2)
               ELSE
                  XA=X(NAX)
                  XB=X(NBX)
               ENDIF
               IF(INDX.EQ.0.OR.INDX.EQ.1) THEN
                  YA=Y(1)*(NAY-1)+Y(2)
                  YB=Y(1)*(NBY-1)+Y(2)
               ELSE
                  YA=Y(NAY)
                  YB=Y(NBY)
               ENDIF
               IF(.NOT.LINV) THEN
                  IF(J.EQ.NFMAX) THEN
c$$$                     WRITE(6,*) 'TYPE A'
c$$$                     DO I=1,J
c$$$                        WRITE(6,'(I5,1P2E12.4,5X,I5,1P2E12.4)')
c$$$     &                       I,XF(I),YF(I)
c$$$                     ENDDO
c$$$                     CALL GUFLSH
c$$$                     PAUSE
                     CALL GUSP2DX(XF(1),YF(1),J,XP,YP,NFMAX,NP,ISPL)
                     CALL SUBV(XP,YP,NP,XG,YG,NGMAX,NN)
                     CALL LINEPTX(XG,YG,NN,ILNS(K))            
                     J=1
                     XF(1)=XF(NFMAX)
                     YF(1)=YF(NFMAX)
c$$$                     WRITE(6,'(A,7I5,1P2E12.4)') 
c$$$     &                    '-3-',J,IEL,IE,NAX,NAY,NBX,NBY,XF(J),YF(J)
                  ENDIF
                  J=J+1
               ELSE
                  IF(J.EQ.1) THEN
c$$$                     WRITE(6,*) 'TYPE B'
c$$$                     DO I=J,NFMAX
c$$$                        WRITE(6,'(I5,1P2E12.4)')
c$$$     &                       I,XF(I),YF(I)
c$$$                     ENDDO
c$$$                     CALL GUFLSH
c$$$                     PAUSE
                     CALL GUSP2DX(XF(J),YF(J),NFMAX-J+1,
     &                           XP,YP,NFMAX,NP,ISPL)
                     CALL SUBV(XP,YP,NP,XG,YG,NGMAX,NN)
                     CALL LINEPTX(XG,YG,NN,ILNS(K))
                     J=NFMAX
                     XF(NFMAX)=XF(1)
                     YF(NFMAX)=YF(1)
c$$$                     WRITE(6,'(A,7I5,1P2E12.4)') 
c$$$     &                    '-2-',J,IEL,IE,NAX,NAY,NBX,NBY,XF(J),YF(J)
                  ENDIF
                  J=J-1
               ENDIF
               RT=(U0-UA)/(UB-UA)
               XF(J)=(XB-XA)*RT+XA
               YF(J)=(YB-YA)*RT+YA
c$$$                     WRITE(6,'(A,7I5,1P2E12.4)') 
c$$$     &                    '-1-',J,IEL,IE,NAX,NAY,NBX,NBY,XF(J),YF(J)
C
               KA(1,IEL)=KA(1,IEL)+1
               IF(IEL.EQ.IE) LEND=.TRUE.
C               IF(IEL.EQ.IE.AND..NOT.LINV) LEND=.TRUE.
C               IF(IEL.EQ.IE1.AND..NOT.LINV) LEND=.TRUE.
            ENDIF
C
            IF(.NOT.LEND) GOTO 200
C
            IF(.NOT.LINV) THEN
c$$$                     WRITE(6,*) 'TYPE C'
c$$$                     DO I=1,J
c$$$                        WRITE(6,'(I5,1P2E12.4,5X,I5,1P2E12.4)')
c$$$     &                       I,XF(I),YF(I)
c$$$                     ENDDO
c$$$                     CALL GUFLSH
               CALL GUSP2DX(XF(1),YF(1),J,XP,YP,NFMAX,NP,ISPL)
            ELSE
c$$$                     WRITE(6,*) 'TYPE D'
c$$$                     DO I=J,NFMAX
c$$$                        WRITE(6,'(I5,1P2E12.4,5X,I5,1P2E12.4)')
c$$$     &                       I,XF(I),YF(I)
c$$$                     ENDDO
c$$$                     CALL GUFLSH
               CALL GUSP2DX(XF(J),YF(J),NFMAX-J+1,XP,YP,NFMAX,NP,ISPL)
            ENDIF
            CALL SUBV(XP,YP,NP,XG,YG,NGMAX,NN)
            CALL LINEPTX(XG,YG,NN,ILNS(K))            
C
C            IF(.NOT.LINV) THEN
C               WRITE(6,*) J,NP,NN
C            ELSE
C               WRITE(6,*) NFMAX-J+1,NP,NN
C            ENDIF
C            CALL GUFLSH
C
         ENDIF
 1000 CONTINUE
      IF(ISPL.GE.0) THEN
         CALL SETLNW(WS)
         CALL SETRGB(RS,GS,BS)
      ENDIF
      RETURN
      END
C
C     ****** CONTOUR PLOT SLAVE ROUTINE : CONV RT ******
C
      SUBROUTINE CONTV2X(RA,TA,N,XB,YB,M,NN)
C
      IMPLICIT LOGICAL(L)
      COMMON /GSGFXY/ DX,DY,PXS,PYS,PXE,PYE,GXS,GYS,GXE,GYE,LGF
      COMMON /GSCTR4/ RMAX,RT,TT,XT,YT
      DIMENSION RA(N),TA(N),XB(M),YB(M)
C
      RT=RA(1)
      TT=TA(1)
      XT=DX*(RT*COS(TT)-GXS)+PXS
      YT=DY*(RT*SIN(TT)-GYS)+PYS
      J=1
      XB(J)=XT
      YB(J)=YT
C
      DO 100 I=2,N
         RS=RA(I)
         TS=TA(I)
         XS=DX*(RS*COS(TS)-GXS)+PXS
         YS=DY*(RS*SIN(TS)-GYS)+PYS
C
         IMAX=INT(SQRT((XS-XT)**2+(YS-YT)**2)*8/RMAX)+1
         DELR=(RS-RT)/IMAX
         DELT= TS-TT
         IF(DELT.GT. 4.0) DELT=DELT-2*3.1415926
         IF(DELT.LT.-4.0) DELT=DELT+2*3.1415926
         DELT=DELT/IMAX
         DO 80 K=1,IMAX
            RI=DELR*K+RT
            TI=DELT*K+TT
            XI=DX*(RI*COS(TI)-GXS)+PXS
            YI=DY*(RI*SIN(TI)-GYS)+PYS
            IF(J.GE.M) THEN
               WRITE(6,*) 'XX GSAF CONTV2 ERROR: BUFFER OVER'
               NN=J
               RETURN
            ENDIF
            J=J+1
            XB(J)=XI
            YB(J)=YI
   80    CONTINUE
         RT=RS
         TT=TS
         XT=XS
         YT=YS
  100 CONTINUE
      NN=J
C
C     ***** slightly expand the region *****
C
      X0=0.D0
      Y0=0.D0
      DO J=1,NN
         X0=X0+XB(J)
         Y0=Y0+YB(J)
      ENDDO
      X0=X0/NN
      Y0=Y0/NN
      DO J=1,NN
         S=SQRT((XB(J)-X0)**2+(YB(J)-Y0)**2)
         IF(ABS(S).LE.1.E-32) S=1.D0
         XB(J)=XB(J)+0.02*(XB(J)-X0)/S
         YB(J)=YB(J)+0.02*(YB(J)-Y0)/S
      ENDDO
      RETURN
      END
C
C     ****** DRAW LINES WITH PATTERN ******
C
      SUBROUTINE LINEPTX(XG,YG,N,IPAT)
C
      DIMENSION XG(N),YG(N)
C
      CALL MOVEPT(XG(1),YG(1),IPAT)
      DO 100 I=2,N
         CALL DRAWPT(XG(I),YG(I))
  100 CONTINUE
      RETURN
      END
C
C     ****** SPLINE INTERPOLATION OF 2D LINES ******
C
      SUBROUTINE GUSP2DX(XH,YH,N,XP,YP,NPM,NP,ISPL)
C
      PARAMETER(NPA=2001)
      PARAMETER(M=3)
      DIMENSION XH(N),YH(N),XP(NPM),YP(NPM)
      DIMENSION IKN(0-M:NPA+M)
C
      IF(ISPL.EQ.0) THEN
         DO I=1,N
            XP(I)=XH(I)
            YP(I)=YH(I)
         ENDDO
         NP=N
         RETURN
      ENDIF
C
      NQ=NPM
      IF(NQ.GT.NPA) NQ=NPA
C
      IF((ABS(XH(N)-XH(1)).LT.1.E-30).AND.
     &   (ABS(YH(N)-YH(1)).LT.1.E-30)) THEN
         IOC=1
      ELSE
         IOC=0
      ENDIF
C
      IF((ISPL.LE.0).OR.((N.EQ.2).AND.(IOC.EQ.0))) THEN
         DO 1 I=1,N
            XP(I)=XH(I)
            YP(I)=YH(I)
    1    CONTINUE
         NP=N
      ELSEIF(N.GT.2) THEN
         NP=ISPL*N
         IF(NP.GT.NQ) NP=NQ
C
         IF(IOC.EQ.0) THEN
            DO 10 I=0-M,0
               IKN(I)=0
   10       CONTINUE
C
            DO 20 I=1,N-1
               IKN(I)=I
   20       CONTINUE
C
            DO 30 I=N,N+M-1
               IKN(I)=N-1
   30       CONTINUE
         ELSEIF(IOC.EQ.1) THEN
            DO 110 I=0-M,N+M-1
               IKN(I)=I
  110       CONTINUE
         ENDIF
         CALL GUCSPLX(N-1,XH,YH,IKN,IOC,NP,XP,YP)
      ENDIF
      RETURN
      END
C
C     ****** SPLINE ******
C
      SUBROUTINE GUCSPLX(N,X,Y,IKN,IOC,NP,XP,YP)
C
      PARAMETER (M=3)
      DIMENSION IKN(0-M:N+M),B(0-M:M,0:M)
      DIMENSION X(0:N),Y(0:N),XP(NP),YP(NP)
C
      H=(IKN(N)-IKN(0))/REAL(NP-1)
C
      DO 100 J=1,NP
         TP=IKN(0)+H*(J-1)
         CALL GUBSPLX(TP,ITM,N,IKN,M,B)
         XV=0
         YV=0
         DO 10 I=0,0-M,-1
            IV=ITM+I+2
            IF(IOC.EQ.1) THEN
               IF(IV.LT.0) THEN
                  IV=IV+N
               ELSEIF(IV.GT.N) THEN
                  IV=IV-N
               ENDIF
            ENDIF
            XV=XV+X(IKN(IV))*B(I,M)
            YV=YV+Y(IKN(IV))*B(I,M)
   10    CONTINUE
         XP(J)=XV
         YP(J)=YV
  100 CONTINUE
      RETURN
      END
C
C     ****** B SPLINE ******
C
      SUBROUTINE GUBSPLX(TP,ITM,N,IKN,M,B)
C
      DIMENSION IKN(0-M:N+M),B(0-M:M,0:M)
C
      DO 10 JT=N-1,0,-1
         IF(TP.GE.REAL(IKN(JT))) THEN
            ITM=JT
            GOTO 11
         ENDIF
   10 CONTINUE
   11 CONTINUE
C
      DO 20 K=0,M-1
      DO 20 I=-1-K,M-K
         B(I,K)=0.0
   20 CONTINUE
C
      B(0,0)=1.0
      DO 30 K=1,M
      DO 30 I=0-K,0
         IV=ITM+I
         IF(IKN(IV+K).GT.IKN(IV)) THEN
            B0=(TP-IKN(IV))/(IKN(IV+K)-IKN(IV)+1.E-30) 
     &           * B(I,K-1)
         ELSE
            B0=0
         ENDIF
C
         IF(IKN(IV+K+1).GT.IKN(IV+1)) THEN
            B1=(IKN(IV+K+1)-TP)/(IKN(IV+K+1)-IKN(IV+1)+1.E-30)
     &           *B(I+1,K-1)
         ELSE
            B1=0
         ENDIF
C
         B(I,K)=B0+B1
   30 CONTINUE
      RETURN
      END
