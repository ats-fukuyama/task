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
C
      DIMENSION TEMP(NTHM,NPM,NRM)
      CHARACTER KID*4,KID1*1,KID2*3
C
    1 WRITE(6,*)'INPUT GRAPH TYPE : F/FX/FS1/FS2 1/2, X:exit,'
      WRITE(6,*)'             : D/DC/DW PP/PT/TP/TT/RR, F/FC/FE P/T/R'
      WRITE(6,*)'             : R/T N/I/W/PC/PW/PE/T/Q/E'
      READ(5,'(A4)',ERR=1,END=9000) KID
      KID1=KID(1:1)
      KID2=KID(2:4)
      CALL GUCPTL(KID1)
      CALL GUCPTL(KID2(1:1))
      CALL GUCPTL(KID2(2:2))
      CALL GUCPTL(KID2(3:3))
C
      IF (KID1.EQ.'F') THEN
         IF(KID2.EQ.'1  ') THEN
            CALL FPGRAP('F1  ',F)
         ELSE IF(KID2.EQ.'2  ') THEN
            CALL FPGRAC('F2  ',F,4)
         ELSE IF(KID2.EQ.'X2  ') THEN
            DO NR=1,NRMAX
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=F(NTH,NP,NR)-F(NTHMAX+1-NTH,NP,NR)
            ENDDO
            ENDDO
            ENDDO
            CALL FPGRAC('FX2 ',TEMP,4)
         ELSE IF(KID2.EQ.'Y2  ') THEN
            DO NR=1,NRMAX
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=F(NTH,NP,NR)-F(NTHMAX+1-NTH,NP,NR)
            ENDDO
            ENDDO
            ENDDO
            CALL FPGRAC('FY2 ',TEMP,0)
         ELSE IF(KID2.EQ.'S11') THEN
            CALL FPGRAP('FS11',FS1)
         ELSE IF(KID2.EQ.'S12') THEN
            CALL FPGRAC('FS12',FS1,4)
         ELSE IF(KID2.EQ.'S21') THEN
            CALL FPGRAP('FS21',FS2)
         ELSE IF(KID2.EQ.'S22') THEN
            CALL FPGRAC('FS22',FS2,4)
         ELSE IF(KID2.EQ.'P  ') THEN
            CALL FPGRAC('FP  ',FPP,1)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRAC('FT  ',FTH,2)
         ELSE IF(KID2.EQ.'R  ') THEN
            CALL FPGRAC('FR  ',FRR,0)
         ELSE IF(KID2.EQ.'PP ') THEN
            CALL FPGRAC('FPP ',FPP,1)
         ELSE IF(KID2.EQ.'TH ') THEN
            CALL FPGRAC('FTH ',FTH,2)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL FPGRAC('FRR ',FRR,0)
         ELSE IF(KID2.EQ.'CP ') THEN
            CALL FPGRAC('FCP ',FCPP,1)
         ELSE IF(KID2.EQ.'CT ') THEN
            CALL FPGRAC('FCT ',FCTH,2)
         ELSE IF(KID2.EQ.'EP ') THEN
            CALL FPGRAC('FEP ',FEPP,1)
         ELSE IF(KID2.EQ.'ET ') THEN
            CALL FPGRAC('FET ',FETH,2)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'R') THEN
         IF(KID2.EQ.'N  ') THEN
            CALL FPGRAR('RN  ',RNT)
         ELSE IF(KID2.EQ.'I  ') THEN
            CALL FPGRAR('RI  ',RJT)
         ELSE IF(KID2.EQ.'W  ') THEN
            CALL FPGRAR('RW  ',RWT)
         ELSE IF(KID2.EQ.'PC ') THEN
            CALL FPGRAR('RPC ',RPCT)
         ELSE IF(KID2.EQ.'PW ') THEN
            CALL FPGRAR('RPW ',RPWT)
         ELSE IF(KID2.EQ.'PE ') THEN
            CALL FPGRAR('RPE ',RPET)
         ELSE IF(KID2.EQ.'LH ') THEN
            CALL FPGRAR('RLH ',RLHT)
         ELSE IF(KID2.EQ.'FW ') THEN
            CALL FPGRAR('RFW ',RFWT)
         ELSE IF(KID2.EQ.'EC ') THEN
            CALL FPGRAR('REC ',RECT)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRAR('RT  ',RTT)
         ELSE IF(KID2.EQ.'Q  ') THEN
            CALL FPGRAR('RQ  ',RQT)
         ELSE IF(KID2.EQ.'E  ') THEN
            CALL FPGRAR('RE  ',RET)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'T') THEN
         IF(KID2.EQ.'N  ') THEN
            CALL FPGRAT('TN  ',PNT)
         ELSE IF(KID2.EQ.'I  ') THEN
            CALL FPGRAT('TI  ',PIT)
         ELSE IF(KID2.EQ.'W  ') THEN
            CALL FPGRAT('TW  ',PWT)
         ELSE IF(KID2.EQ.'PC ') THEN
            CALL FPGRAT('TPC ',PPCT)
         ELSE IF(KID2.EQ.'PW ') THEN
            CALL FPGRAT('TPW ',PPWT)
         ELSE IF(KID2.EQ.'PE ') THEN
            CALL FPGRAT('TPE ',PPET)
         ELSE IF(KID2.EQ.'LH ') THEN
            CALL FPGRAT('TLH ',PLHT)
         ELSE IF(KID2.EQ.'FW ') THEN
            CALL FPGRAT('TFW ',PFWT)
         ELSE IF(KID2.EQ.'EC ') THEN
            CALL FPGRAT('TEC ',PECT)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRAT('TT  ',PTT)
         ELSE IF(KID2.EQ.'Q  ') THEN
            CALL FPGRAT('TQ  ',PQT)
         ELSE IF(KID2.EQ.'E  ') THEN
            CALL FPGRAT('TE  ',PET)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'D') THEN
         IF(KID2.EQ.'PP ') THEN
            CALL FPGRAC('DPP ',DPP ,1)
         ELSE IF(KID2.EQ.'PT ') THEN
            CALL FPGRAC('DPT ',DPT ,1)
         ELSE IF(KID2.EQ.'TP ') THEN
            CALL FPGRAC('DTP ',DTP ,2)
         ELSE IF(KID2.EQ.'TT ') THEN
            CALL FPGRAC('DTT ',DTT ,2)
         ELSE IF(KID2.EQ.'CPP') THEN
            CALL FPGRAC('DCPP',DCPP,1)
         ELSE IF(KID2.EQ.'CPT') THEN
            CALL FPGRAC('DCPT',DCPT,1)
         ELSE IF(KID2.EQ.'CTP') THEN
            CALL FPGRAC('DCTP',DCTP,2)
         ELSE IF(KID2.EQ.'CTT') THEN
            CALL FPGRAC('DCTT',DCTT,2)
         ELSE IF(KID2.EQ.'W  ') THEN
            CALL FPGRAC2('DW  ',DWPP,DWTT,0)
         ELSE IF(KID2.EQ.'WPP') THEN
            CALL FPGRAC('DWPP',DWPP,1)
         ELSE IF(KID2.EQ.'WPT') THEN
            CALL FPGRAC('DWPT',DWPT,1)
         ELSE IF(KID2.EQ.'WTP') THEN
            CALL FPGRAC('DWTP',DWTP,2)
         ELSE IF(KID2.EQ.'WTT') THEN
            CALL FPGRAC('DWTT',DWTT,2)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL FPGRAC('DRR ',DRR ,0)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'X') THEN
         GO TO 9000
      ELSE IF (KID1.EQ.'Q') THEN
         GO TO 9000
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID1'
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
      SUBROUTINE FPGRAR(STRING,FR)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION FR(NRM,NTG1M),GX(NRM),GY(NRM)
      CHARACTER STRING*4
C
      CALL PAGES
      CALL SETCHS(0.3,0.)
      CALL SETFNT(32)
      GYMAX0=-1.E30
      GYMIN0= 1.E30
      DO 230 NT=1,NTG1
        DO 240 NR=1,NRMAX
          GY(NR)=GUCLIP(FR(NR,NT))
          GX(NR)=GUCLIP(RM(NR))
  240   CONTINUE
        CALL GMNMX1(GY,1,NRMAX,1,GYMIN,GYMAX)
        GYMAX0=MAX(GYMAX0,GYMAX)
        GYMIN0=MIN(GYMIN0,GYMIN)
  230 CONTINUE
C
      IF(GYMIN0 .GE. 0.) GYMIN0=0.
      IF(GYMAX0 .LE. 0.) GYMAX0=0.
      CALL GQSCAL(GYMIN0,GYMAX0,GYMIN1,GYMAX1,GYSTEP)
      CALL GQSCAL(GUCLIP(RHOGMN),GUCLIP(RHOGMX),GXMIN1,GXMAX1,GXSTEP)
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
      CALL TEXT(STRING,4)
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
      SUBROUTINE FPGRAT(STRING,FT)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION  FT(NTG2M),GX(NTG2M),GY(NTG2M)
      CHARACTER STRING*4
C
      CALL PAGES
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
      CALL GSCALE(0.,GXSTEP,0.,GYSTEP,1.0,0)
      CALL GFRAME
      CALL GVALUE(0.,GXSTEP*2,0.,0.,NGULEN(2*GXSTEP))
      CALL GVALUE(0.,0.,0.,GYSTEP*2,NGULEN(2*GYSTEP))
      CALL GPLOTP(GX,GY,1,NTG2,1,0,0,0)
      CALL MOVE(1.0,17.5)
      CALL TEXT(STRING,4)
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
      SUBROUTINE FPGRAP(STRING,FG)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION FG(NTHM,NPM,NRM),GX(NPM),GY(NPM,NRM)
      CHARACTER STRING*4
C
    1 WRITE(6,*) '# INPUT NTH (1..',NTHMAX,') :'
      READ(5,*,ERR=1,END=9000) NTH
      IF(NTH.LT.1) GOTO 9000
      IF(NTH.GT.NTHMAX) GOTO 1
C
      CALL PAGES
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
      DO 200 NP=1,NPMAX
         GX(NP)=GUCLIP(PM(NP)**2)
  200 CONTINUE
      GXMAX=GUCLIP(PMAX**2)
CXX
      GXMAX=100.0
CXX
      CALL GMNMX2(GY,NPM,1,NPMAX,1,1,NRMAX,1,GYMIN,GYMAX)
      CALL GQSCAL(GYMIN,GYMAX,GYMINP,GYMAXP,GYSTEP)
      CALL GQSCAL(0.0,GXMAX,GXMINP,GXMAXP,GXSTEP)
C
      CALL GDEFIN(3.0,23.0,2.0,17.0,0.0,GXMAX,-8.0,0.0)
      CALL GFRAME
      CALL GSCALE(0.,GXSTEP*2,0.0,0.0,1.0,0)
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
      CALL TEXT(STRING,4)
      CALL MOVE(24.0,1.0)
      CALL TEXT('p**2',4)
      CALL PAGEE
C
 9000 RETURN
      END
C
C ***********************************************************
C
C                        C-GRAPHIC
C
C ***********************************************************
C
      SUBROUTINE FPGRAC(STRING,FG,MODE)
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION FG(*)
      DIMENSION GF(NPMP,NTHMP),GP(NPMP),GTH(NTHMP),KA(8,NPMP,NTHMP)
      CHARACTER STRING*4
C
    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.0) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,*) '# INPUT NR (1..',NRG,') :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF
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
            NM=NPMP*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            GF(NP,NTH)=GUCLIP(FG(NM))
   60    CONTINUE
      ELSEIF(MODE.EQ.2) THEN
         DO 70 NTH=1,NTHMAX+1
         DO 70 NP=1,NPMAX
            NM=NPM*NTHMP*(NR-1)+NTHMP*(NP-1)+NTH
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
C
      GPMAX=GUCLIP(PMAX)
C
      CALL PAGES
      CALL SETCHS(.3,0.)
      CALL SETFNT(32)
C
      CALL GMNMX2(GF,NPMP,1,NPG,1,1,NTHG,1,GFMIN,GFMAX)
      CALL GQSCAL(GFMIN,GFMAX,GFMIN1,GFMAX1,GFSTEP)
      CALL GQSCAL(0.0,GPMAX,GPMIN1,GPMAX1,GPSTEP)
C
      CALL GDEFIN(3.,23.,2.,12.,-GPMAX,GPMAX,0.,GPMAX)
      CALL GFRAME
      CALL GSCALE(0.,GPSTEP,0.,GPSTEP,1.0,0)
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,NGULEN(2*GPSTEP))
C
      IF(LMODE.EQ.0) THEN
         IF(GFMIN*GFMAX.GE.0.0) THEN
            IF(GFMIN.GE.0.0) THEN
               CALL CONTQ4(GF,GP,GTH,NPMP,NPG,NTHG,
     &                     GFMIN1,0.5*GFSTEP,30,0,KA)
            ELSE
               CALL CONTQ4(GF,GP,GTH,NPMP,NPG,NTHG,
     &                     GFMIN1,0.5*GFSTEP,30,2,KA)
            ENDIF
         ELSE
            CALL CONTQ4(GF,GP,GTH,NPMP,NPG,NTHG,
     &                   0.25*GFSTEP, 0.5*GFSTEP,15,0,KA)
            CALL CONTQ4(GF,GP,GTH,NPMP,NPG,NTHG,
     &                  -0.25*GFSTEP,-0.5*GFSTEP,15,2,KA)
         ENDIF
      ELSE
         DO 100 I=1,NGLINE
           GLIN=GFMAX-0.020*(I-1)**2
           CALL SETLIN(0,0,7-MOD(I-1,5))
           CALL CONTQ4(GF,GP,GTH,NPMP,NPG,NTHG,
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
      CALL TEXT(STRING,4)
      CALL MOVE(5.0,12.5)
      CALL TEXT('FMIN =',6)
      CALL NUMBR(GFMIN,'(1PE12.4)',12)
      CALL MOVE(11.0,12.5)
      CALL TEXT('FMAX =',6)
      CALL NUMBR(GFMAX,'(1PE12.4)',12)
      IF(LMODE.EQ.0) THEN
         CALL MOVE(17.0,12.5)
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
C                        C-GRAPHIC
C
C ***********************************************************
C
      SUBROUTINE FPGRAC2(STRING,FG,FH,MODE)
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION FG(*),FH(*)
      DIMENSION GF(NPMP,NTHMP),GP(NPMP),GTH(NTHMP),KA(8,NPMP,NTHMP)
      CHARACTER STRING*4
C
    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.0) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,*) '# INPUT NR (1..',NRG,') :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF
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
            NMP1=NPMP*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NMP2=NPMP*NTHM*(NR-1)+NTHM* NP   +NTH
            FGA=0.5D0*(FG(NMP1)+FG(NMP2))
            NMT1=NPM*NTHMP*(NR-1)+NTHMP*(NP-1)+NTH
            NMT2=NPM*NTHMP*(NR-1)+NTHMP*(NP-1)+NTH+1
            FHA=0.5D0*(FH(NMT1)+FH(NMT2))
            GF(NP,NTH)=GUCLIP(SQRT(FGA**2+FHA**2))
   50    CONTINUE
      ELSEIF(MODE.EQ.1) THEN
         DO 60 NTH=1,NTHMAX
         DO 60 NP=1,NPMAX+1
            NM=NPMP*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            GF(NP,NTH)=GUCLIP(SQRT(FG(NM)**2+FH(NM)**2))
   60    CONTINUE
      ELSEIF(MODE.EQ.2) THEN
         DO 70 NTH=1,NTHMAX+1
         DO 70 NP=1,NPMAX
            NM=NPM*NTHMP*(NR-1)+NTHMP*(NP-1)+NTH
            GF(NP,NTH)=GUCLIP(SQRT(FG(NM)**2+FH(NM)**2))
   70    CONTINUE
      ELSEIF(MODE.EQ.4) THEN
         DO 80 NTH=1,NTHMAX
         DO 80 NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            IF(FG(NM).LT.1.D-70) THEN
               GF(NP,NTH)=-70.0
            ELSE
               GF(NP,NTH)=GUCLIP(LOG10(ABS(SQRT(FG(NM)**2+FH(NM)**2))))
            ENDIF
   80    CONTINUE
      ENDIF         
C
      GPMAX=GUCLIP(PMAX)
C
      CALL PAGES
      CALL SETCHS(.3,0.)
      CALL SETFNT(32)
C
      CALL GMNMX2(GF,NPMP,1,NPG,1,1,NTHG,1,GFMIN,GFMAX)
      CALL GQSCAL(GFMIN,GFMAX,GFMIN1,GFMAX1,GFSTEP)
      CALL GQSCAL(0.0,GPMAX,GPMIN1,GPMAX1,GPSTEP)
C
      CALL GDEFIN(3.,23.,2.,12.,-GPMAX,GPMAX,0.,GPMAX)
      CALL GFRAME
      CALL GSCALE(0.,GPSTEP,0.,GPSTEP,1.0,0)
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,NGULEN(2*GPSTEP))
CXX
      GFSTEP=2000.0
CXX
      IF(LMODE.EQ.0) THEN
         IF(GFMIN*GFMAX.GE.0.0) THEN
            IF(GFMIN.GE.0.0) THEN
               CALL CONTQ4(GF,GP,GTH,NPMP,NPG,NTHG,
     &                     GFMIN1,0.5*GFSTEP,30,0,KA)
            ELSE
               CALL CONTQ4(GF,GP,GTH,NPMP,NPG,NTHG,
     &                     GFMIN1,0.5*GFSTEP,30,2,KA)
            ENDIF
         ELSE
            CALL CONTQ4(GF,GP,GTH,NPMP,NPG,NTHG,
     &                   0.25*GFSTEP, 0.5*GFSTEP,20,0,KA)
            CALL CONTQ4(GF,GP,GTH,NPMP,NPG,NTHG,
     &                  -0.25*GFSTEP,-0.5*GFSTEP,20,2,KA)
         ENDIF
      ELSE
         DO 100 I=1,NGLINE
           GLIN=GFMAX-0.020*(I-1)**2
           CALL SETLIN(0,2,7-MOD(I-1,5))
           CALL CONTQ4(GF,GP,GTH,NPMP,NPG,NTHG,
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
      CALL TEXT(STRING,4)
      CALL MOVE(5.0,12.5)
      CALL TEXT('FMIN =',6)
      CALL NUMBR(GFMIN,'(1PE12.4)',12)
      CALL MOVE(11.0,12.5)
      CALL TEXT('FMAX =',6)
      CALL NUMBR(GFMAX,'(1PE12.4)',12)
      IF(LMODE.EQ.0) THEN
         CALL MOVE(17.0,12.5)
         CALL TEXT('STEP =',6)
         CALL NUMBR(0.5*GFSTEP,'(1PE12.4)',12)
      ENDIF
      CALL PAGEE
      IF(NRMAX.GT.1) GOTO 1
C
 9000 RETURN
      END
