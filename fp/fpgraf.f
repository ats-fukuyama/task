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
      INCLUDE 'fpcomm.h'
C
      DIMENSION TEMP(NTHM,NPM,NRM)
      CHARACTER KID*4,KID1*1,KID2*3
C
    1 WRITE(6,*)'INPUT GRAPH TYPE : F/FX/FS1/FS2 1/2, X,'
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
         IF(KID2.EQ.'1  ') CALL FPGRAP('F   ',F)
         IF(KID2.EQ.'2  ') CALL FPGRAC('F   ',F,4)
         IF(KID2.EQ.'X2  ') THEN
            DO NR=1,NRMAX
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=F(NTH,NP,NR)-F(NTHMAX+1-NTH,NP,NR)
            ENDDO
            ENDDO
            ENDDO
            CALL FPGRAC('FX  ',TEMP,4)
         ENDIF
         IF(KID2.EQ.'Y2  ') THEN
            DO NR=1,NRMAX
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=F(NTH,NP,NR)-F(NTHMAX+1-NTH,NP,NR)
            ENDDO
            ENDDO
            ENDDO
            CALL FPGRAC('FX  ',TEMP,0)
         ENDIF
         IF(KID2.EQ.'S11') CALL FPGRAP('F   ',FS1)
         IF(KID2.EQ.'S12') CALL FPGRAC('F   ',FS1,4)
         IF(KID2.EQ.'S21') CALL FPGRAP('F   ',FS2)
         IF(KID2.EQ.'S22') CALL FPGRAC('F   ',FS2,4)
         IF(KID2.EQ.'P  ') CALL FPGRAC('FP  ',FPP,1)
         IF(KID2.EQ.'T  ') CALL FPGRAC('FTH ',FTH,2)
         IF(KID2.EQ.'R  ') CALL FPGRAC('FR  ',FRR,0)
         IF(KID2.EQ.'PP ') CALL FPGRAC('FP  ',FPP,1)
         IF(KID2.EQ.'TH ') CALL FPGRAC('FTH ',FTH,2)
         IF(KID2.EQ.'RR ') CALL FPGRAC('FR  ',FRR,0)
         IF(KID2.EQ.'CP ') CALL FPGRAC('FCP ',FCPP,1)
         IF(KID2.EQ.'CT ') CALL FPGRAC('FCTH',FCTH,2)
         IF(KID2.EQ.'EP ') CALL FPGRAC('FEP ',FEPP,1)
         IF(KID2.EQ.'ET ') CALL FPGRAC('FETH',FETH,2)
      ELSE IF (KID1.EQ.'R') THEN
         IF(KID2.EQ.'N  ') CALL FPGRAR('N   ',RNT)
         IF(KID2.EQ.'I  ') CALL FPGRAR('I   ',RJT)
         IF(KID2.EQ.'W  ') CALL FPGRAR('W   ',RWT)
         IF(KID2.EQ.'PC ') CALL FPGRAR('PC  ',RPCT)
         IF(KID2.EQ.'PW ') CALL FPGRAR('PW  ',RPWT)
         IF(KID2.EQ.'PE ') CALL FPGRAR('PE  ',RPET)
         IF(KID2.EQ.'LH ') CALL FPGRAR('LH  ',RLHT)
         IF(KID2.EQ.'FW ') CALL FPGRAR('FW  ',RFWT)
         IF(KID2.EQ.'EC ') CALL FPGRAR('EC  ',RECT)
         IF(KID2.EQ.'T  ') CALL FPGRAR('T   ',RTT)
         IF(KID2.EQ.'Q  ') CALL FPGRAR('Q   ',RQT)
         IF(KID2.EQ.'E  ') CALL FPGRAR('E   ',RET)
      ELSE IF (KID1.EQ.'T') THEN
         IF(KID2.EQ.'N  ') CALL FPGRAT('N   ',PNT)
         IF(KID2.EQ.'I  ') CALL FPGRAT('I   ',PIT)
         IF(KID2.EQ.'W  ') CALL FPGRAT('W   ',PWT)
         IF(KID2.EQ.'PC ') CALL FPGRAT('PC  ',PPCT)
         IF(KID2.EQ.'PW ') CALL FPGRAT('PW  ',PPWT)
         IF(KID2.EQ.'PE ') CALL FPGRAT('PE  ',PPET)
         IF(KID2.EQ.'LH ') CALL FPGRAT('LH  ',PLHT)
         IF(KID2.EQ.'FW ') CALL FPGRAT('FW  ',PFWT)
         IF(KID2.EQ.'EC ') CALL FPGRAT('EC  ',PECT)
         IF(KID2.EQ.'T  ') CALL FPGRAT('T   ',PTT)
         IF(KID2.EQ.'Q  ') CALL FPGRAT('Q   ',PQT)
         IF(KID2.EQ.'E  ') CALL FPGRAT('E   ',PET)
      ELSE IF (KID1.EQ.'D') THEN
         IF(KID2.EQ.'PP ') CALL FPGRAC('DPP ',DPP ,1)
         IF(KID2.EQ.'PT ') CALL FPGRAC('DPT ',DPT ,1)
         IF(KID2.EQ.'TP ') CALL FPGRAC('DTP ',DTP ,2)
         IF(KID2.EQ.'TT ') CALL FPGRAC('DTT ',DTT ,2)
         IF(KID2.EQ.'CPP') CALL FPGRAC('DCPP',DCPP,1)
         IF(KID2.EQ.'CPT') CALL FPGRAC('DCPT',DCPT,1)
         IF(KID2.EQ.'CTP') CALL FPGRAC('DCTP',DCTP,2)
         IF(KID2.EQ.'CTT') CALL FPGRAC('DCTT',DCTT,2)
         IF(KID2.EQ.'W  ') CALL FPGRAC2('DW  ',DWPP,DWTT,0)
         IF(KID2.EQ.'WPP') CALL FPGRAC('DWPP',DWPP,1)
         IF(KID2.EQ.'WPT') CALL FPGRAC('DWPT',DWPT,1)
         IF(KID2.EQ.'WTP') CALL FPGRAC('DWTP',DWTP,2)
         IF(KID2.EQ.'WTT') CALL FPGRAC('DWTT',DWTT,2)
         IF(KID2.EQ.'RR ') CALL FPGRAC('DRR ',DRR ,0)
      ELSE IF (KID1.EQ.'X') THEN
         GO TO 9000
      ELSE IF (KID1.EQ.'Q') THEN
         GO TO 9000
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
      INCLUDE 'fpcomm.h'
      DIMENSION FR(NRM,NTG1M),GX(NRM),GY(NRM)
      CHARACTER STRING*4
C
      CALL PAGES
      CALL SETCHS(0.3,0.)
      GYMAX0=-1.E30
      GYMIN0= 1.E30
      DO 230 NT=1,NTG1
        DO 240 NR=1,NRMAX
          GY(NR)=GCLIP(FR(NR,NT))
          GX(NR)=GCLIP(RM(NR))
  240   CONTINUE
        CALL GMNMX1(GY,1,NRMAX,1,GYMIN,GYMAX)
        GYMAX0=MAX(GYMAX0,GYMAX)
        GYMIN0=MIN(GYMIN0,GYMIN)
  230 CONTINUE
C
      IF(GYMIN0 .GE. 0.) GYMIN0=0.
      IF(GYMAX0 .LE. 0.) GYMAX0=0.
      CALL GQSCAL(GYMIN0,GYMAX0,GYMIN1,GYMAX1,GYSTEP)
      CALL GQSCAL(0.0,1.0,GXMIN1,GXMAX1,GXSTEP)
      CALL GDEFIN(3.0,23.0,2.0,17.0,0.0,1.0,GYMIN1,GYMAX1)
      CALL GSCALE(0.,GXSTEP,0.0,GYSTEP,1.0,0)
      CALL GFRAME
      CALL GVALUE(0.,2*GXSTEP,0.0,0.0,NGVLEN(2*GXSTEP))
      CALL GVALUE(0.,0.0,0.0,2*GYSTEP,NGVLEN(2*GYSTEP))
      DO 250 NT=1,NTG1
        DO 260 NR=1,NRMAX
          GY(NR)=GCLIP(FR(NR,NT))
          GX(NR)=GCLIP(RM(NR))
  260   CONTINUE
        CALL SETLIN(0,0,7-MOD(NT-1,5))
        CALL GPLOTP(GX,GY,1,NRMAX,1,0,0,0)
  250 CONTINUE
      CALL SETLIN(0,0,7)
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
      INCLUDE 'fpcomm.h'
      DIMENSION  FT(NTG2M),GX(NTG2M),GY(NTG2M)
      CHARACTER STRING*4
C
      CALL PAGES
      CALL SETCHS(0.3,0.0)
C
      DO 10 N=1,NTG2
        GX(N)=GCLIP(PTG(N))
        GY(N)=GCLIP(FT(N))
   10 CONTINUE
      CALL GMNMX1(GY,1,NTG2,1,GYMIN,GYMAX)
      IF(GYMIN.GT.0.) GYMIN=0.0
      IF(GYMAX.LT.0.) GYMAX=0.0
      CALL GQSCAL(GYMIN,GYMAX,GYMIN1,GYMAX1,GYSTEP)
      CALL GQSCAL(0.,GX(NTG2),GXMIN,GXMAX,GXSTEP)
      CALL GDEFIN(3.,23.,2.,17.,0.,GX(NTG2),GYMIN1,GYMAX1)
      CALL GSCALE(0.,GXSTEP,0.,GYSTEP,1.0,0)
      CALL GFRAME
      CALL GVALUE(0.,GXSTEP*2,0.,0.,NGVLEN(2*GXSTEP))
      CALL GVALUE(0.,0.,0.,GYSTEP*2,NGVLEN(2*GYSTEP))
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
      INCLUDE 'fpcomm.h'
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
      DO 100 NR=1,NRMAX
      DO 100 NP=1,NPMAX
         IF(FG(NTH,NP,NR).LT.1.D-14) THEN
            GY(NP,NR)=-14.0
         ELSE
            GY(NP,NR)=GCLIP(LOG10(FG(NTH,NP,NR)))
         ENDIF
  100 CONTINUE
      DO 200 NP=1,NPMAX
         GX(NP)=GCLIP(PM(NP)**2)
  200 CONTINUE
      GXMAX=GCLIP(PMAX**2)
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
      CALL GVALUE(0.,GXSTEP*2,0.0,0.0,NGVLEN(2*GXSTEP))
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
      INCLUDE 'fpcomm.h'
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
            GP(NP)=GCLIP(PM(NP))
   10    CONTINUE
         NPG=NPMAX
      ELSE
         DO 20 NP=1,NPMAX+1
            GP(NP)=GCLIP(PG(NP))
   20    CONTINUE
         NPG=NPMAX+1
      ENDIF
C
      IF(MOD(MODE/2,2).EQ.0) THEN
         DO 30 NTH=1,NTHMAX
            GTH(NTH)=GCLIP(THM(NTH))
   30    CONTINUE
         NTHG=NTHMAX
      ELSE
         DO 40 NTH=1,NTHMAX+1
            GTH(NTH)=GCLIP(THG(NTH))
   40    CONTINUE
         NTHG=NTHMAX+1
      ENDIF
C
      IF(MODE.EQ.0) THEN
         DO 50 NTH=1,NTHMAX
         DO 50 NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            GF(NP,NTH)=GCLIP(FG(NM))
   50    CONTINUE
      ELSEIF(MODE.EQ.1) THEN
         DO 60 NTH=1,NTHMAX
         DO 60 NP=1,NPMAX+1
            NM=NPMP*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            GF(NP,NTH)=GCLIP(FG(NM))
   60    CONTINUE
      ELSEIF(MODE.EQ.2) THEN
         DO 70 NTH=1,NTHMAX+1
         DO 70 NP=1,NPMAX
            NM=NPM*NTHMP*(NR-1)+NTHMP*(NP-1)+NTH
            GF(NP,NTH)=GCLIP(FG(NM))
   70    CONTINUE
      ELSEIF(MODE.EQ.4) THEN
         DO 80 NTH=1,NTHMAX
         DO 80 NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            IF(FG(NM).LT.1.D-70) THEN
               GF(NP,NTH)=-70.0
            ELSE
               GF(NP,NTH)=GCLIP(LOG10(ABS(FG(NM))))
            ENDIF
   80    CONTINUE
      ENDIF         
C
      GPMAX=GCLIP(PMAX)
C
      CALL PAGES
      CALL SETCHS(.3,0.)
C
      CALL GMNMX2(GF,NPMP,1,NPG,1,1,NTHG,1,GFMIN,GFMAX)
      CALL GQSCAL(GFMIN,GFMAX,GFMIN1,GFMAX1,GFSTEP)
      CALL GQSCAL(0.0,GPMAX,GPMIN1,GPMAX1,GPSTEP)
C
      CALL GDEFIN(3.,23.,2.,12.,-GPMAX,GPMAX,0.,GPMAX)
      CALL GFRAME
      CALL GSCALE(0.,GPSTEP,0.,GPSTEP,1.0,0)
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,NGVLEN(2*GPSTEP))
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
      CALL SETLIN(0,0,7)
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
      INCLUDE 'fpcomm.h'
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
            GP(NP)=GCLIP(PM(NP))
   10    CONTINUE
         NPG=NPMAX
      ELSE
         DO 20 NP=1,NPMAX+1
            GP(NP)=GCLIP(PG(NP))
   20    CONTINUE
         NPG=NPMAX+1
      ENDIF
C
      IF(MOD(MODE/2,2).EQ.0) THEN
         DO 30 NTH=1,NTHMAX
            GTH(NTH)=GCLIP(THM(NTH))
   30    CONTINUE
         NTHG=NTHMAX
      ELSE
         DO 40 NTH=1,NTHMAX+1
            GTH(NTH)=GCLIP(THG(NTH))
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
            GF(NP,NTH)=GCLIP(SQRT(FGA**2+FHA**2))
   50    CONTINUE
      ELSEIF(MODE.EQ.1) THEN
         DO 60 NTH=1,NTHMAX
         DO 60 NP=1,NPMAX+1
            NM=NPMP*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            GF(NP,NTH)=GCLIP(SQRT(FG(NM)**2+FH(NM)**2))
   60    CONTINUE
      ELSEIF(MODE.EQ.2) THEN
         DO 70 NTH=1,NTHMAX+1
         DO 70 NP=1,NPMAX
            NM=NPM*NTHMP*(NR-1)+NTHMP*(NP-1)+NTH
            GF(NP,NTH)=GCLIP(SQRT(FG(NM)**2+FH(NM)**2))
   70    CONTINUE
      ELSEIF(MODE.EQ.4) THEN
         DO 80 NTH=1,NTHMAX
         DO 80 NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            IF(FG(NM).LT.1.D-70) THEN
               GF(NP,NTH)=-70.0
            ELSE
               GF(NP,NTH)=GCLIP(LOG10(ABS(SQRT(FG(NM)**2+FH(NM)**2))))
            ENDIF
   80    CONTINUE
      ENDIF         
C
      GPMAX=GCLIP(PMAX)
C
      CALL PAGES
      CALL SETCHS(.3,0.)
C
      CALL GMNMX2(GF,NPMP,1,NPG,1,1,NTHG,1,GFMIN,GFMAX)
      CALL GQSCAL(GFMIN,GFMAX,GFMIN1,GFMAX1,GFSTEP)
      CALL GQSCAL(0.0,GPMAX,GPMIN1,GPMAX1,GPSTEP)
C
      CALL GDEFIN(3.,23.,2.,12.,-GPMAX,GPMAX,0.,GPMAX)
      CALL GFRAME
      CALL GSCALE(0.,GPSTEP,0.,GPSTEP,1.0,0)
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,NGVLEN(2*GPSTEP))
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
           CALL SETLIN(0,0,7-MOD(I-1,5))
           CALL CONTQ4(GF,GP,GTH,NPMP,NPG,NTHG,
     &                 GLIN,GFSTEP*100,1,0,KA)
  100    CONTINUE
      ENDIF
C
      CALL SETLIN(0,0,7)
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
