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
C     ----- exec EQ -----
C
      IF(MODELG.EQ.3) THEN
         CALL EQLOAD(3,KNAMEQ,IERR)
         IF(IERR.EQ.0) THEN
            CALL EQSETP
            CALL EQPSIC(51,32,64,IERR)
            CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
         ENDIF
      ENDIF
C
C     ----- set radial mesh -----
C
      IF(NRMAX.EQ.1) THEN
         DELR=DELR1
      ELSE
         DELR=(RMAX-RMIN)/NRMAX
      ENDIF
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
C     ----- load WR resluts -----
C
      IF(MODELW.EQ.2.OR.MODELW.EQ.3) THEN
         CALL FPLDWR(IERR)
         IF(IERR.NE.0) RETURN
      ENDIF
C
C     ----- set parameters for target species -----
C
      AEFP=PZ(NSFP)*AEE
      AMFP=PA(NSFP)*AMP
      RNFP0=PN(NSFP)
      RNFPS=PNS(NSFP)
      RTFP0=(PTPR(NSFP)+2.D0*PTPP(NSFP))/3.D0
      RTFPS=PTS(NSFP)
C
      IF(MODELR.EQ.0) THEN
         PTFP0=SQRT(RTFP0*1.D3*AEE*AMFP)
         VTFP0=SQRT(RTFP0*1.D3*AEE/AMFP)
      ELSE
         RKE=RTFP0*1.D3*AEE
         PTFP0=SQRT(RKE*RKE+2.D0*RKE*AMFP*VC*VC)/VC
         VTFP0=PTFP0/SQRT(AMFP**2+PTFP0**2/VC**2)
      ENDIF
C
C     ----- set profile data -----
C
      DO NR=1,NRMAX
C
         PSIL=RM(NR)**2
         CALL PLPROF(PSIL)
C
         RNFP(NR)=RN(NSFP)
         RTFP(NR)=(RTPR(NSFP)+2.D0*RTPP(NSFP))/3.D0
         IF(MODELR.EQ.0) THEN
            PTFP(NR)=SQRT(RTFP(NR)*1.D3*AEE*AMFP)
            VTFP(NR)=SQRT(RTFP(NR)*1.D3*AEE/AMFP)
         ELSE
            RKE=RTFP(NR)*1.D3*AEE
            PTFP(NR)=SQRT(RKE*RKE+2.D0*RKE*AMFP*VC*VC)/VC
            VTFP(NR)=PTFP(NR)/SQRT(AMFP**2+PTFP(NR)**2/VC**2)
         ENDIF
         RNE=RN(1)
         RTE=(RTPR(1)+2.D0*RTPP(1))/3.D0
         DO NS=1,NSMAX
            AEFD=PZ(NS)*AEE
            AMFD=PA(NS)*AMP
            RNFD(NR,NS)=RN(NS)
            RTFD(NR,NS)=(RTPR(NS)+2.D0*RTPP(NS))/3.D0
            IF(MODELR.EQ.0) THEN
               PTFD(NR,NS)=SQRT(RTFD(NR,NS)*1.D3*AEE*AMFD)
               VTFD(NR,NS)=SQRT(RTFD(NR,NS)*1.D3*AEE/AMFD)
            ELSE
               RKE=RTFD(NR,NS)*1.D3*AEE
               PTFD(NR,NS)=SQRT(RKE*RKE+2.D0*RKE*AMFD*VC*VC)/VC
               VTFD(NR,NS)=PTFD(NR,NS)
     &                    /SQRT(AMFD**2+PTFD(NR,NS)**2/VC**2)
            ENDIF
            IF(NSFP.EQ.1.AND.NS.EQ.1) THEN
               RLNRL=14.9D0-0.5D0*LOG(RNE)+LOG(RTE)
            ELSEIF(NSFP.EQ.1.OR.NS.EQ.1) THEN
               RLNRL=15.2D0-0.5D0*LOG(RNE)+LOG(RTE)
            ELSE
               RLNRL=17.3D0-0.5D0*LOG(RNE)+1.5D0*LOG(RTFD(NR,NS))
            ENDIF
            FACT=AEFP**2*AEFD**2*RLNRL/(4.D0*PI*EPS0**2)
            RNUF(NR,NS)=FACT*RNFD(1,NS)*1.D20
     &                 /(2.D0*AMFD*VTFD(1,NS)**2*PTFP(1))
            RNUD(NR,NS)=FACT*RNFD(1,NS)*1.D20
     &                 /(SQRT(2.D0)*VTFD(1,NS)*PTFP(1)**2)
         ENDDO
      ENDDO
C
C     ----- set poloidal magneticl field -----
C
      DO NR=1,NRMAX+1
         PSIL=RG(NR)**2
         CALL FPSETB(PSIL,0.5D0*PI,BT,BP(NR))
         EPSR(NR)=RSPSIN(PSIL)/RR
      ENDDO
C
C     ----- set parallel current density -----
C
      DO NR=1,NRMAX
         RJ1(NR)=(RG(NR+1)*BP(NR+1)-RM(NR)*BP(NR))
     &          /(RMU0*RM(NR)*DELR)
      ENDDO
C
C     ----- set parallel electric field -----
C
      DO NR=1,NRMAX
         E1(NR)=E0
         E2(NR)=E0
         RJ2(NR)=RJ1(NR)
      ENDDO
C
C     ----- set momentum space mesh -----
C
      DELP =PMAX/NPMAX
      DELTH=PI/NTHMAX
C
      DO NP=1,NPMAX
        PG(NP)=DELP*(NP-1)
        PM(NP)=DELP*(NP-0.5D0)
      ENDDO
      PG(NPMAX+1)=PMAX
C
      DO NTH=1,NTHMAX
         THG(NTH)=DELTH*(NTH-1)
         THM(NTH)=DELTH*(NTH-0.5D0)
C
         SINM(NTH)=SIN(THM(NTH))
         COSM(NTH)=COS(THM(NTH))
         SING(NTH)=SIN(THG(NTH))
         COSG(NTH)=COS(THG(NTH))
      ENDDO
      THG(NTHMAX+1)=PI
      SING(NTHMAX+1)=0.D0
      COSG(NTHMAX+1)=-1.D0
C
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX
         VOL(NTH,NP)=2.D0*PI*SINM(NTH)*PM(NP)**2*DELP*DELTH
      ENDDO
      ENDDO
C
C     ----- set relativistic parameters -----
C
      IF (MODELR.EQ.0) THEN
C
         THETA0=0.D0
         DO NR=1,NRMAX
            THETA(NR)=0.D0
            DKBSR(NR)=0.D0
         ENDDO
C
      ELSE
C
         THETA0=RTFP0*1.D3*AEE/(AMFP*VC*VC)
         DO NR=1,NRMAX
            THETA(NR)=THETA0*RTFP(NR)
            Z=1.D0/THETA(NR)
C            IF(Z.LE.100.D0) THEN
               DKBSR(NR)=DKBES(2,Z)
C            ELSE
C               DKBSR(NR)=SQRT(0.5D0*PI/Z)*(1.D0+15.D0/(8.D0*Z))
C            ENDIF
         ENDDO
      ENDIF
C
C     ----- set boundary distribution functions -----
C
      DO NP=1,NPMAX
         FL=FPMXWL(PM(NP),0)
         DO NTH=1,NTHMAX
            FS1(NTH,NP)=FL
         ENDDO
      ENDDO
C
      DO NP=1,NPMAX
         FL=FPMXWL(PM(NP),NRMAX+1)
         DO NTH=1,NTHMAX
            FS2(NTH,NP)=FL
         ENDDO
      ENDDO
C
C     ----- set bounce-average parameters -----
C
      IF (MODELA.NE.0) THEN
         DO NR=1,NRMAX
            A1=ACOS(SQRT(2.D0*EPSR(NR)/(1.D0+EPSR(NR))))
            DO NTH=1,NTHMAX/2
               IF (THG(NTH).LE.A1.AND.THG(NTH+1).GE.A1) GOTO 201
            ENDDO
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
            DO NTH=1,ITL(NR)
               ETAM(NTH,NR)=PI/2.D0
            ENDDO
C
            DO NTH=ITL(NR)+1,ITU(NR)-1
               A1=FACT*COSM(NTH)**2
               ETAM(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
            ENDDO
C
            DO NTH=ITU(NR),NTHMAX
               ETAM(NTH,NR)=PI/2.D0
            ENDDO
C
            DO NTH=1,ITL(NR)
               ETAG(NTH,NR)=PI/2.D0
            ENDDO
C
            DO NTH=ITL(NR)+1,ITU(NR)
               A1=FACT*COSG(NTH)**2
               ETAG(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
            ENDDO
C
            DO NTH=ITU(NR)+1,NTHMAX+1
               ETAG(NTH,NR)=PI/2.D0
            ENDDO
C
            NRX=NR
            DO NTH=1,NTHMAX/2
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
            ENDDO
            RLAMDA(ITL(NR),NR)=0.5D0*(RLAMDA(ITL(NR)-1,NR)
     &                               +RLAMDA(ITL(NR)+1,NR))
            DO NTH=1,NTHMAX/2
               RLAMDA(NTHMAX-NTH+1,NR)=RLAMDA(NTH,NR)
               RLAMDC(NTHMAX-NTH+2,NR)=RLAMDC(NTH,NR)
            ENDDO
            RLAMDC(NTHMAX/2+1,NR)=0.D0
         ENDDO
C
      ELSE
         DO NR=1,NRMAX
            ITL(NR)=0
            ITU(NR)=0
            DO NTH=1,NTHMAX
               ETAM(NTH,NR)=PI/2.D0
               RLAMDA(NTH,NR)=1.D0
            ENDDO
            DO NTH=1,NTHMAX+1
               ETAG(NTH,NR)=PI/2.D0
            ENDDO
         ENDDO
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
      DO NR=1,NRMAX
      DO NP=1,NPMAX
         FL=FPMXWL(PM(NP),NR)
         DO NTH=1,NTHMAX
            F(NTH,NP,NR)=FL
         ENDDO
      ENDDO
      ENDDO
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
         RNFPL=RN(NSFP)/RNFP0
         RTFPL=RTPR(NSFP)/RTFP0
      ELSEIF(NR.EQ.NRMAX+1) THEN
         RL=RM(NRMAX)+DELR
         PSIN=RL**2
         CALL PLPROF(PSIN)
         RNFPL=RN(NSFP)/RNFP0
         RTFPL=RTPR(NSFP)/RTFP0
      ELSE
         RNFPL=RNFP(NR)
         RTFPL=RTFP(NR)
      ENDIF
C
      IF (MODELR.EQ.0) THEN
C
         FACT=RNFPL/SQRT(2.D0*PI*RTFPL)**3
         EX=PML**2/(2.D0*RTFPL)            
         IF(EX.GT.100.D0) THEN
            FPMXWL=0.D0
         ELSE
            FPMXWL=FACT*EXP(-EX)
         ENDIF
C
      ELSE
C
         IF(NR.EQ.0.OR.NR.EQ.NRMAX+1) THEN
            THETAL=THETA0*RTFPL
            Z=1.D0/THETAL
            DKBSL=DKBES(2,Z)
         ELSE
            THETAL=THETA(NR)
            DKBSL=DKBSR(NR)
         ENDIF
         FACT=RNFPL*SQRT(THETA0)/(4.D0*PI*RTFPL*DKBSL)
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
      DO NT=1,NTMAX
C
         L=0
C
         IF(MODELE.NE.0) THEN
            DO NR=1,NRMAX
               E3(NR)=0.D0
               RJ3(NR)=0.D0
            ENDDO
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
            DO NR=2,NRMAX
               RSUM=0.D0
               DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  RSUM=RSUM+VOL(NTH,NP)*F1(NTH,NP,NR)*PM(NP)
               ENDDO
               ENDDO
               RJN(NR)=AEFP*RNFP0*1.D20*PTH0*DELP*RSUM/(AMFP*RM(NR)*RA)
            ENDDO
            RJN(1)=(4.D0*RJN(2)-RJN(3))/3.D0
C
            DELEM=0.D0
            DO NR=1,NRMAX
               IF(ABS(RJN(NR)-RJ3(NR)).GT.1.D-20) THEN
                  DELE(NR)=(RJN(NR)-RJ2(NR))*(E2(NR)-E3(NR))
     &                    /(RJN(NR)-RJ3(NR))
                  E3(NR)=E2(NR)
                  RJ3(NR)=RJN(NR)
                  E2(NR)=E2(NR)-DELE(NR)
                  DELEM=MAX(ABS(DELE(NR))/MAX(ABS(E1(NR)),1.D-6),DELEM)
               ENDIF
            ENDDO
C
            IF (L.LT.LMAXE.AND.DELEM.GT.EPSE) GO TO 1
            IF (L.GE.LMAXE) WRITE(6,*) 'L IS LARGER THAN LMAXE'
C
            DO NR=1,NRMAX
               E1(NR)=E2(NR)
               RJ1(NR)=RJN(NR)
            ENDDO
            CALL FPNEWE
         ENDIF
C
  250    CONTINUE
         DO NR=1,NRMAX
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            F(NTH,NP,NR)=F1(NTH,NP,NR)
         ENDDO
         ENDDO
         ENDDO
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
      ENDDO
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
      DO NR=2,NRMAX
         BP(NR)=BP(NR)+(E1(NR)-E1(NR-1))*DELT/(RA*DELR)
      ENDDO
C
      DO NR=1,NRMAX
         RJ2(NR)=(RG(NR+1)*BP(NR+1)-RG(NR)*BP(NR))
     &           /(RMU0*RM(NR)*DELR*RA)
      ENDDO
C
      DO NR=1,NRMAX
         E2(NR)=RJ2(NR)*E1(NR)/RJ1(NR)
      ENDDO
C
      RETURN
      END
