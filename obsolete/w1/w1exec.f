C
C   ***** TASK/W1 W1EXEC: Calculate waves *****
C
      SUBROUTINE W1EXEC
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      PARAMETER (NHM=2*NHARMM+1,NHM1=NHARMM+1)
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM3/ PA(ISM),PZ(ISM),PN(ISM),PTPP(ISM),PTPR(ISM),
     &                PNS(ISM),PTS(ISM),PU(ISM),PZCL(ISM)
      COMMON /W1PRM4/ RZ,DZ,NZP
      COMMON /W1PRM6/ IELEC(ISM)
      COMMON /W1ANT1/ AJYL(IAM),AJYH(IAM),AJZL(IAM),AJZH(IAM),
     &                ALYL(IAM),ALYH(IAM),APYL(IAM),APYH(IAM)
      COMMON /W1BND1/ CGIN(3,5),CGOT(3,5),CFJY1,CFJY2,CFJZ1,CFJZ2
      COMMON /W1EF2D/ CE2DA(NZPM,NXTM,3)
      COMMON /W1ZDAT/ CJ1(NZPM),CJ2(NZPM),CJ3(NZPM),CJ4(NZPM),
     &                ZA(NZPM),AKZ(NZPM),NZANT1(IAM),NZANT2(IAM),NANT
      COMMON /W1CTRL/ DRF,DRKZ,DXFACT,DXWDTH,APRFPN,APRFTR,APRFTP,
     &                NPRINT,NFILE,NGRAPH,NLOOP,NSYM,NFLR,NALPHA,
     &                NSYS,NDISP,NXABS
      COMMON /W1QCTL/ XDMAX,DXD,NDMAX,NHARM,IHARM(ISM),MATL,NMODEL
      COMMON /W1PRM5/ EPSH
      COMMON /W1XDAT/ XA(NXQ),XAM(NXTM)
      COMMON /W1CDDT/ AJCDX(NXPM),AJCDK(NZPM),AJCDT,ZEFF,WVYSIZ,NCDTYP
C
C     ******* 2-DIMENSIONAL ANALYSIS *******
C
      RFSAVE=RF
      RKSAVE=RKZ
      DO 2000 NL=1 , NLOOP
         IF(NLOOP.NE.1) WRITE(6,630) RF,RKZ
         CALL W1SETZ(IERR)
            IF(IERR.NE.0) GOTO 2000
         CALL W1ANTS
         CALL W1SETX(IERR)
            IF(IERR.NE.0) GOTO 2000
         CALL W1PROF
         CALL W1PWRI
C
C     ******* FOURIER TRANSFORM OF ANTENNA CURRENT *******
C
         CALL W1FFTL(CJ1,NZP,0)
         CALL W1FFTL(CJ2,NZP,0)
         CALL W1FFTL(CJ3,NZP,0)
         CALL W1FFTL(CJ4,NZP,0)
C
C     ******* CALCULATION FOR EACH KZ *******
C
         DO 1000 NZDO = 1 , NZP
            NZ=NZDO
            IF(NZ.LE.(NZP/2+1).OR.NSYM.EQ.0) THEN
               RKZ   = AKZ(NZ)
               CFJY1 = CJ1(NZ)
               CFJY2 = CJ2(NZ)
               CFJZ1 = CJ3(NZ)
               CFJZ2 = CJ4(NZ)
               CALL W1BCND
               IF(NMODEL.LE.5) THEN
                  CALL W1DSPA(NALPHA)
               ELSE
                  CALL W1DSPQ(ICL)
               ENDIF
               IF(NMODEL.EQ.0.OR.NMODEL.EQ.2) THEN
                  CALL W1BNDA(IERR)
                     IF(IERR.NE.0) GOTO 2000
                  CALL W1EPWA(NZ)
               ELSEIF(NMODEL.EQ.1.OR.NMODEL.EQ.3) THEN
                  CALL W1WKXB
                  CALL W1BNDB(IERR,NXABS)
                     IF(IERR.NE.0) GOTO 2000
                  CALL W1EPWB(NZ)
                  CALL W1HELD(4)
               ELSEIF(NMODEL.EQ.4) THEN
                  CALL W1BNDC(IERR)
                     IF(IERR.NE.0) GOTO 2000
                  CALL W1EPWC(NZ)
               ELSEIF(NMODEL.EQ.5) THEN
                  CALL W1WKXD
                  CALL W1BNDD(IERR,NXABS)
                     IF(IERR.NE.0) GOTO 2000
                  CALL W1EPWD(NZ)
                  CALL W1HELD(6)
               ELSEIF(NMODEL.EQ.6) THEN
                  CALL W1BNDQ(IERR)
                     IF(IERR.NE.0) GOTO 2000
                  CALL W1EPWQ(NZ)
               ENDIF
               CALL W1EVAC(NZ,NSYM)
               CALL W1CLCD(NZ)
               CALL W1CLPW(NZ,NSYM)
            ELSE
               CALL W1SYMS(NZ,NSYM)
            ENDIF
 1000    CONTINUE
         IF(NMODEL.EQ.6) WRITE(6,616) MATL,MATLM,ICL,NCLM,NDMAX,XDMAX
C
C     ******* INVERSE FOURIER TRANSFORM *******
C
         CALL W1FFTL(CJ1,NZP,1)
         CALL W1FFTL(CJ2,NZP,1)
         CALL W1FFTL(CJ3,NZP,1)
         CALL W1FFTL(CJ4,NZP,1)
C
         DO 1100 NX=1,NXT
         DO 1100 IC=1,3
            CALL W1FFTL(CE2DA(1,NX,IC),NZP,1)
 1100    CONTINUE
C
C     ******* POWER ABSORPTION AND OUTPUT *******
C
         CALL W1PWRS
         CALL W1PRNT(NPRINT)
         CALL W1FILE(NFILE)
         IF(NGRAPH.GT.0) THEN
            IF(NZP.EQ.1) THEN
               CALL W1GR1D(MOD(NGRAPH,   4))
               CALL W1GR1B(MOD(NGRAPH/4, 2))
               CALL W1GR1F(MOD(NGRAPH/8, 2))
               CALL W1GR1H(MOD(NGRAPH/16,2))
               CALL W1GRUD(MOD(NGRAPH/32,2))
               CALL W1GRUF(MOD(NGRAPH/64,2))
            ELSE
               CALL W1GR2D(MOD(NGRAPH,   4))
               CALL W1GR1D(MOD(NGRAPH,   4))
            ENDIF
         ENDIF
         IF(NLOOP.NE.1) THEN
            RF   = RF   + DRF
            RKZ  = RKZ  + DRKZ
         ENDIF
 2000 CONTINUE
C
      RF =RFSAVE
      RKZ=RKSAVE
      CALL FCLOCK(TT2)
      WRITE(6,699) (TT2-TT1)
  699 FORMAT(1H ,'## CPU TIME = ',F15.3,' SEC')
      RETURN
  616 FORMAT('## INTEGRO-DIFF EQ. : ',
     &           'MATL    MATLM    ICL     NCLM    NDMAX   XDMAX'/
     &       22X,5I8,1P1D12.4)
  630 FORMAT(1H /
     &       1H ,'** RF = ',1PD12.4,' (MHZ) , RKZ = ',1PD12.4,
     &           ' (/M) **')
      END
