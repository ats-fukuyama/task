C     $Id$
C
C *****************************
C     PREPARATION OF FPLOOP
C *****************************
C
      SUBROUTINE FPPREP(IERR)
C
      INCLUDE 'fpcomm.inc'
C
      EXTERNAL FPFN0U,FPFN0T,FPFN1A,FPFN2A
C
      IF(MODELG.EQ.3) THEN
         CALL EQLOAD(3,KNAMEQ,IERR)
         IF(IERR.EQ.0) THEN
            CALL EQSETP
            CALL EQPSIC(51,32,64,IERR)
            CALL EQGETB(BB,RR,RIP,RA,RKAP,RDEL,RB)
         ENDIF
      ENDIF
C
      IF(NRMAX.EQ.1) THEN
         DELR=DELR1
      ELSE
         DELR=(RMAX-RMIN)/NRMAX
      ENDIF
C
C      BETAN=LOG(RNE0/RNES)/RA**ALPHN
C      BETAT=LOG(TE0 /TES )/RA**ALPHT
C
      IF(NRMAX.EQ.1) THEN
         RM(1)=R1
         RG(1)=R1-0.5D0*DELR
         RG(2)=R1+0.5D0*DELR
      ELSE
         DO 10 NR=1,NRMAX
            RM(NR)=RMIN+DELR*(NR-1)+0.5D0*DELR
            RG(NR)=RMIN+DELR*(NR-1)
   10    CONTINUE
         RG(NRMAX+1)=RMAX
      ENDIF
C
      IF(MODELW.EQ.2.OR.MODELW.EQ.3) THEN
C      IF(MODELW.EQ.2) THEN
         CALL FPLDWR(IERR)
         IF(IERR.NE.0) RETURN
      ELSE
         RNE0=PN(1)
         RNES=PNS(1)
         TE0=(PTPR(1)+2.D0*PTPP(1))/3.D0
         TES=PTS(1)
      ENDIF
C
      DO 12 NR=1,NRMAX
C
C     ***** EQUI-SPACE PSI *****
C
C         PSIL=RM(NR)
C
C     ***** EQUI-SPACE RADIUS *****
C
         PSIL=RM(NR)**2
C
         CALL PLPROF(PSIL)
         RNE(NR)=RN(1)/RNE0
         TE(NR) =RTPR(1)/TE0
C
C         RNE(NR)=EXP(-BETAN*RM(NR)**ALPHN)
C         TE(NR) =EXP(-BETAT*RM(NR)**ALPHT)
C
         VTE=SQRT(TE0*TE(NR)*1.D3*AEE/AME)
         WPE(NR)=SQRT(RNE0*RNE(NR)*1.D20*AEE**2/(AME*EPS0))
         RNU(NR)=15.D0*WPE(NR)**4/(4.D0*PI*RNE0*RNE(NR)*1.D20*VTE**3)
         PTH(NR)=SQRT(TE0*TE(NR)*1.D3*AEE*AME)
         EPSR(NR)=RSPSIN(PSIL)/RR
   12 CONTINUE
C
      DO 13 NR=1,NRMAX+1
         CALL FPSETB(RG(NR)*RG(NR),0.5D0*PI,BT,BP(NR))
   13 CONTINUE
C
      DO 20 NR=1,NRMAX
         RJ1(NR)=(RG(NR+1)*BP(NR+1)-RM(NR)*BP(NR))
     &          /(RMU0*RM(NR)*DELR)
   20 CONTINUE
C
      DO 30 NR=1,NRMAX
         E1(NR)=E0
         E2(NR)=E0
         RJ2(NR)=RJ1(NR)
   30 CONTINUE
C
      VTE0=SQRT(TE0*1.D3*AEE/AME)
      WPE0=SQRT(RNE0*1.D20*AEE**2/(AME*EPS0))
      RNU0=15.D0*WPE0**4/(4.D0*PI*RNE0*1.D20*VTE0**3)
      PTH0=SQRT(TE0*1.D3*AEE*AME)
      DELP =PMAX/NPMAX
      DELTH=PI/NTHMAX
C
      DO 40 NP=1,NPMAX
        PG(NP)=DELP*(NP-1)
        PM(NP)=DELP*(NP-0.5D0)
   40 CONTINUE
      PG(NPMAX+1)=PMAX
C
      DO 50 NTH=1,NTHMAX
         THG(NTH)=DELTH*(NTH-1)
         THM(NTH)=DELTH*(NTH-0.5D0)
C
         SINM(NTH)=SIN(THM(NTH))
         COSM(NTH)=COS(THM(NTH))
         SING(NTH)=SIN(THG(NTH))
         COSG(NTH)=COS(THG(NTH))
   50 CONTINUE
         THG(NTHMAX+1)=PI
         SING(NTHMAX+1)=0.D0
         COSG(NTHMAX+1)=-1.D0
C
      DO 60 NP=1,NPMAX
      DO 60 NTH=1,NTHMAX
         VOL(NTH,NP)=2.D0*PI*SINM(NTH)*PM(NP)**2*DELP*DELTH
   60 CONTINUE
C
C =============   BOUNDARY DISTRIBUTION  =============
C
      IF (MODELR.EQ.0) THEN
C
         THETA0=0.D0
         DO 100 NR=1,NRMAX
            THETA(NR)=0.D0
            DKBSR(NR)=0.D0
  100    CONTINUE
C
      ELSE
C
         THETA0=TE0/511.D0
         DO 110 NR=1,NRMAX
            THETA(NR)=TE0*TE(NR)/511.D0
            Z=1.D0/THETA(NR)
C            IF(Z.LE.100.D0) THEN
               DKBSR(NR)=DKBES(2,Z)
C            ELSE
C               DKBSR(NR)=SQRT(0.5D0*PI/Z)*(1.D0+15.D0/(8.D0*Z))
C            ENDIF
  110    CONTINUE
      ENDIF
C
      DO 120 NP=1,NPMAX
         FL=FPMXWL(PM(NP),0)
      DO 120 NTH=1,NTHMAX
         FS1(NTH,NP)=FL
  120 CONTINUE
C
      DO 130 NP=1,NPMAX
         FL=FPMXWL(PM(NP),NRMAX+1)
      DO 130 NTH=1,NTHMAX
         FS2(NTH,NP)=FL
  130 CONTINUE
C
C ================  BOUNCE AVERAGE  ==========================
C
      IF (MODELA.NE.0) THEN
         DO 290 NR=1,NRMAX
            A1=ACOS(SQRT(2.D0*EPSR(NR)/(1.D0+EPSR(NR))))
            DO 200 NTH=1,NTHMAX/2
               IF (THG(NTH).LE.A1.AND.THG(NTH+1).GE.A1) GOTO 201
  200       CONTINUE
C
  201       CONTINUE
C
            WRITE(6,*) 'NR,NTHC=',NR,NTH
            IF(NTH.EQ.NTHMAX/2) NTH=NTHMAX/2-1
C
            ITL(NR)=NTH
            ITU(NR)=NTHMAX-NTH+1
            EPSL=COSM(ITL(NR))**2/(2.D0-COSM(ITL(NR))**2)
            EPSR(NR)=EPSL
            FACT=(1.D0+EPSL)/(2.D0*EPSL)
C
            DO 210 NTH=1,ITL(NR)
               ETAM(NTH,NR)=PI/2.D0
  210       CONTINUE
C
            DO 220 NTH=ITL(NR)+1,ITU(NR)-1
               A1=FACT*COSM(NTH)**2
               ETAM(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
  220       CONTINUE
C
            DO 230 NTH=ITU(NR),NTHMAX
               ETAM(NTH,NR)=PI/2.D0
  230       CONTINUE
C
            DO 240 NTH=1,ITL(NR)
               ETAG(NTH,NR)=PI/2.D0
  240       CONTINUE
C
            DO 250 NTH=ITL(NR)+1,ITU(NR)
               A1=FACT*COSG(NTH)**2
               ETAG(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
  250       CONTINUE
C
            DO 260 NTH=ITU(NR)+1,NTHMAX+1
               ETAG(NTH,NR)=PI/2.D0
  260       CONTINUE
C
C -----------------  CALCULATION OF RLAMDA  ------------------
C
            NRX=NR
            DO 270 NTH=1,NTHMAX/2
               NTHX=NTH
               IF(NTH.LT.ITL(NR)) THEN
                  CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPFN0U)
               ELSEIF(NTH.GT.ITL(NR)) THEN
                  CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPFN0T)
               ELSE
                  RINT0=0.D0
               ENDIF
               CALL DEFT(RINT2,ES,H0DE,EPSDE,0,FPFN2A)
               RLAMDA(NTH,NR)=RINT0*ABS(COSM(NTH))/PI
               RLAMDC(NTH,NR)=RINT2/(PI*(1.D0+EPSR(NR))*ABS(COSG(NTH)))
  270       CONTINUE
            RLAMDA(ITL(NR),NR)=0.5D0*(RLAMDA(ITL(NR)-1,NR)
     &                               +RLAMDA(ITL(NR)+1,NR))
            DO 280 NTH=1,NTHMAX/2
               RLAMDA(NTHMAX-NTH+1,NR)=RLAMDA(NTH,NR)
               RLAMDC(NTHMAX-NTH+2,NR)=RLAMDC(NTH,NR)
  280       CONTINUE
               RLAMDC(NTHMAX/2+1,NR)=0.D0
C
  290    CONTINUE
      ELSE
         DO 320 NR=1,NRMAX
            ITL(NR)=0
            ITU(NR)=0
            DO 300 NTH=1,NTHMAX
               ETAM(NTH,NR)=PI/2.D0
               RLAMDA(NTH,NR)=1.D0
  300       CONTINUE
            DO 310 NTH=1,NTHMAX+1
               ETAG(NTH,NR)=PI/2.D0
  310       CONTINUE
  320    CONTINUE
      END IF
C
      IERR=0
      RETURN
      END
C
C ****************************************
C     INITIALIZE VELOCITY DISTRIBUTION    
C ****************************************
C
      SUBROUTINE FPFINI
C
      INCLUDE 'fpcomm.inc'
C
      DO 100 NR=1,NRMAX
      DO 100 NP=1,NPMAX
         FL=FPMXWL(PM(NP),NR)
      DO 100 NTH=1,NTHMAX
         F(NTH,NP,NR)=FL
  100 CONTINUE
      RETURN
      END
C
C ****************************************
C     MAXWELLIAN VELOCITY DISTRIBUTION
C ****************************************
C
      FUNCTION FPMXWL(PML,NR)
C
      INCLUDE 'fpcomm.inc'
C
      IF(NR.EQ.0) THEN
         RL=RM(1)-DELR
         PSIN=RL**2
         CALL PLPROF(PSIN)
         RNEL=RN(1)/RNE0
         TEL =RTPR(1)/TE0
      ELSEIF(NR.EQ.NRMAX+1) THEN
         RL=RM(NRMAX)+DELR
         PSIN=RL**2
         CALL PLPROF(PSIN)
         RNEL=RN(1)/RNE0
         TEL =RTPR(1)/TE0
      ELSE
         RNEL=RNE(NR)
         TEL =TE (NR)
      ENDIF
C
      IF (MODELR.EQ.0) THEN
C
         FACT=RNEL/SQRT(2.D0*PI*TEL)**3
         EX=PML**2/(2.D0*TEL)            
         IF(EX.GT.100.D0) THEN
            FPMXWL=0.D0
         ELSE
            FPMXWL=FACT*EXP(-EX)
         ENDIF
C
      ELSE
C
         IF(NR.EQ.0.OR.NR.EQ.NRMAX+1) THEN
            THETAL=TE0*TEL/511.D0
            Z=1.D0/THETAL
C            IF(Z.LE.100.D0) THEN
               DKBSL=DKBES(2,Z)
C            ELSE
C               DKBSL=SQRT(0.5D0*PI/Z)*(1.D0+15.D0/(8.D0*Z))
C            ENDIF
         ELSE
            THETAL=THETA(NR)
            DKBSL=DKBSR(NR)
         ENDIF
         FACT=RNEL*SQRT(THETA0)/(4.D0*PI*TEL*DKBSL)
         EX=(1.D0-SQRT(1.D0+PML**2*THETA0))/THETAL
         IF(EX.LT.-100.D0) THEN
            FPMXWL=0.D0
         ELSE
            FPMXWL=FACT*EXP(EX)
         ENDIF
      ENDIF
      RETURN
      END
C
C *************************
C     INITIAL DATA SAVE
C *************************
C
      SUBROUTINE FPSAVI
C
      INCLUDE 'fpcomm.inc'
C
      ISAVE=0
      DO NR=1,NRMAX
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX
         DCPP(NTH,NP,NR)=0.D0
         DCPT(NTH,NP,NR)=0.D0
         FCPP(NTH,NP,NR)=0.D0
         DWPP(NTH,NP,NR)=0.D0
         DWPT(NTH,NP,NR)=0.D0
         FEPP(NTH,NP,NR)=0.D0
         DWLHPP(NTH,NP,NR)=0.D0
         DWLHPT(NTH,NP,NR)=0.D0
         DWFWPP(NTH,NP,NR)=0.D0
         DWFWPT(NTH,NP,NR)=0.D0
         DWECPP(NTH,NP,NR)=0.D0
         DWECPT(NTH,NP,NR)=0.D0
      ENDDO
      ENDDO
      ENDDO
      CALL FPSPRF
      CALL FPSGLB
      CALL FPWRIT
C
      RETURN
      END
C
C ***************************************************************
C
C                       SET OF INTEGRAND
C
C ***************************************************************
C
C ============================================================
C
      REAL*8 FUNCTION  FPFN0U(X,XM,XP)
C
      INCLUDE 'fpcomm.inc'
C
      XX=X
      XX=XM
      A0=ETAM(NTHX,NRX)
      A1=1.D0+EPSR(NRX)*COS(A0*XP)
      A2=(1.D0+EPSR(NRX))*SINM(NTHX)**2
      FPFN0U=A0*SQRT(A1/(A1-A2))
C
      RETURN
      END
C
C ============================================================
C
      REAL*8 FUNCTION FPFN0T(X,XM,XP)
C
      INCLUDE 'fpcomm.inc'
C
      XX=X
      A0=ETAM(NTHX,NRX)
      A1=1.D0+EPSR(NRX)*COS(A0*XP)
      A2=2.D0*EPSR(NRX)*SIN(0.5D0*A0*XM)*SIN(0.5D0*A0*(XP+2.D0))
      FPFN0T=A0*SQRT(A1/A2)
C
      RETURN
      END
C
C ============================================================
C
      REAL*8 FUNCTION FPFN1A(X,XM,XP)
C
      INCLUDE 'fpcomm.inc'
C
      XX=X
      XX=XM
      A0=ETAM(NTHX,NRX)
      A1=1.D0+EPSR(NRX)*COS(A0*XP)
      FPFN1A=A0*SQRT(A1)
C
      RETURN
      END
C
C ============================================================
C
      REAL*8 FUNCTION FPFN2A(X,XM,XP)
C
      INCLUDE 'fpcomm.inc'
C
      XX=X
      XX=XM
      A0=ETAG(NTHX,NRX)
      A1=1.D0+EPSR(NRX)*COS(A0*XP)
      A2=A1-(1.D0+EPSR(NRX))*SING(NTHX)**2
      FPFN2A=A0*SQRT(A1*A2)
C
      RETURN
      END
C
C *****************
C     MAIN LOOP
C *****************
C
      SUBROUTINE FPLOOP
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION RJN(NRM),RJ3(NRM),E3(NRM),DELE(NRM)
      
C
      IF(MODELE.NE.0) CALL FPNEWE
C
      DO 1000 NT=1,NTMAX
C
         L=0
C
         IF(MODELE.NE.0) THEN
            DO 100 NR=1,NRMAX
               E3(NR)=0.D0
               RJ3(NR)=0.D0
  100       CONTINUE
         ENDIF
C
    1    L=L+1
C
         IF (MOD(NT-1,NTSTPC).EQ.0) CALL FPCOEF
C
         CALL FPEXEC(NOCONV)
         IF(NOCONV.NE.0) GOTO 250
C
         IF(MODELE.NE.0) THEN
            DO 210 NR=2,NRMAX
               RSUM=0.D0
               DO 200 NP=1,NPMAX
               DO 200 NTH=1,NTHMAX
                  RSUM=RSUM+VOL(NTH,NP)*F1(NTH,NP,NR)*PM(NP)
  200          CONTINUE
               RJN(NR)=-AEE*RNE0*1.D20*PTH0*DELP*RSUM/(AME*RM(NR)*RA)
  210       CONTINUE
            RJN(1)=(4.D0*RJN(2)-RJN(3))/3.D0
C
            DELEM=0.D0
            DO 220 NR=1,NRMAX
               IF(ABS(RJN(NR)-RJ3(NR)).GT.1.D-20) THEN
                  DELE(NR)=(RJN(NR)-RJ2(NR))*(E2(NR)-E3(NR))
     &                    /(RJN(NR)-RJ3(NR))
                  E3(NR)=E2(NR)
                  RJ3(NR)=RJN(NR)
                  E2(NR)=E2(NR)-DELE(NR)
                  DELEM=MAX(ABS(DELE(NR))/MAX(ABS(E1(NR)),1.D-6),DELEM)
               ENDIF
  220       CONTINUE
C
            IF (L.LT.LMAXE.AND.DELEM.GT.EPSE) GO TO 1
            IF (L.GE.LMAXE) WRITE(6,*) 'L IS LARGER THAN LMAXE'
C
            DO 230 NR=1,NRMAX
               E1(NR)=E2(NR)
               RJ1(NR)=RJN(NR)
  230       CONTINUE
            CALL FPNEWE
         ENDIF
C
  250    CONTINUE
         DO 300 NR=1,NRMAX
         DO 300 NP=1,NPMAX
         DO 300 NTH=1,NTHMAX
            F(NTH,NP,NR)=F1(NTH,NP,NR)
  300    CONTINUE
C
         TIMEFP=TIMEFP+DELT
C
         ISAVE=0
         IF (MOD(NT,NTSTP1).EQ.0) THEN
            CALL FPSPRF
         ENDIF
         IF (MOD(NT,NTSTP2).EQ.0) THEN
            CALL FPSGLB
            CALL FPWRIT
         ENDIF
C
         IF(NOCONV.NE.0) GOTO 1100
 1000 CONTINUE
 1100 CONTINUE
C
      RETURN
      END
C
C ************************************
C     PREDICTION OF ELECTRIC FIELD
C ************************************
C

      SUBROUTINE FPNEWE
C
      INCLUDE 'fpcomm.inc'
C
      DO 10 NR=2,NRMAX
         BP(NR)=BP(NR)+(E1(NR)-E1(NR-1))*DELT/(RA*DELR)
   10 CONTINUE
C
      DO 20 NR=1,NRMAX
         RJ2(NR)=(RG(NR+1)*BP(NR+1)-RG(NR)*BP(NR))
     &           /(RMU0*RM(NR)*DELR*RA)
   20 CONTINUE
C
      DO 30 NR=1,NRMAX
         E2(NR)=RJ2(NR)*E1(NR)/RJ1(NR)
   30 CONTINUE
C
      RETURN
      END
