C     $Id$
C
C     ***** INITIALIZE AK-EQ INTERFACE *****
C
      SUBROUTINE AKEQIN(KNAMEQ1,NRMAX,NTHMAX,NSUMAX,
     &                  ALPMAX,RAXIS,ZAXIS,IERR)
C
      IMPLICIT REAL*8(A-F,H,O-Z)
      CHARACTER KNAMEQ1*32
C
      IERR=0
      CALL EQLOAD(1,KNAMEQ1,IERR)
      IF(IERR.NE.0) RETURN
C
      CALL EQSETP
      CALL EQPSIC(NRMAX,NTHMAX,NSUMAX,IERR)
      IF(IERR.NE.0) RETURN
      ALPMAX=FNFTS(FNPSIN(1.D0))
      CALL GETAXS(RAXIS,ZAXIS)
C
      CALL EQCALA(IERR)
      RETURN
      END
C
C     ***** INITIALIZE RZ-AB CONVERSION *****
C
      SUBROUTINE EQCALA(IERR)
C
      INCLUDE '../eq/eqcomq.h'
      COMMON /EQAKV1/ ALPHG(NRM),BETAG(NRM)
      COMMON /EQAKV2/ UABR(4,NTHM,NRM),UABZ(4,NTHM,NRM)
      DIMENSION RPSA(NTHM,NRM),RPSB(NTHM,NRM),RPSAB(NTHM,NRM)
      DIMENSION ZPSA(NTHM,NRM),ZPSB(NTHM,NRM),ZPSAB(NTHM,NRM)
C
      DTH=2.D0*PI/NTHMAX
      DO NR=1,NRMAX
         ALPHG(NR)=SQRT(FTS(NR))
      ENDDO
      DO NTH=1,NTHMAX
         BETAG(NTH)=DTH*(NTH-1)
      ENDDO
C
      DO NR=1,NRMAX
         RPSB(     1,NR)=(RPS(2,NR)-RPS(NTHMAX  ,NR))/(2.D0*DTH)
         RPSB(NTHMAX,NR)=(RPS(1,NR)-RPS(NTHMAX-1,NR))/(2.D0*DTH)
      ENDDO
      RPSAB(     1,    1)=(RPSB(     1,    2)-RPSB(1     ,      1))
     &                   /(ALPHG(    2)-ALPHG(      1))
      RPSAB(NTHMAX,    1)=(RPSB(NTHMAX,    2)-RPSB(NTHMAX,      1))
     &                   /(ALPHG(    2)-ALPHG(      1))
      RPSAB(     1,NRMAX)=(RPSB(     1,NRMAX)-RPSB(1     ,NRMAX-1))
     &                   /(ALPHG(NRMAX)-ALPHG(NRMAX-1))
      RPSAB(NTHMAX,NRMAX)=(RPSB(NTHMAX,NRMAX)-RPSB(NTHMAX,NRMAX-1))
     &                   /(ALPHG(NRMAX)-ALPHG(NRMAX-1))
C
      CALL SPL2D(BETAG,ALPHG,RPS,RPSB,RPSA,RPSAB,UABR,
     &           NTHM,NTHMAX,NRMAX,0,0,IERR)
      IF(IERR.NE.0) THEN
         IERR=IERR+10000
         RETURN
      ENDIF
C
      DO NR=1,NRMAX
         ZPSB(     1,NR)=(ZPS(2,NR)-ZPS(NTHMAX  ,NR))/(2.D0*DTH)
         ZPSB(NTHMAX,NR)=(ZPS(1,NR)-ZPS(NTHMAX-1,NR))/(2.D0*DTH)
      ENDDO
      ZPSAB(     1,    1)=(ZPSB(     1,    2)-ZPSB(1     ,      1))
     &                   /(ALPHG(    2)-ALPHG(      1))
      ZPSAB(NTHMAX,    1)=(ZPSB(NTHMAX,    2)-ZPSB(NTHMAX,      1))
     &                   /(ALPHG(    2)-ALPHG(      1))
      ZPSAB(     1,NRMAX)=(ZPSB(     1,NRMAX)-ZPSB(1     ,NRMAX-1))
     &                   /(ALPHG(NRMAX)-ALPHG(NRMAX-1))
      ZPSAB(NTHMAX,NRMAX)=(ZPSB(NTHMAX,NRMAX)-ZPSB(NTHMAX,NRMAX-1))
     &                   /(ALPHG(NRMAX)-ALPHG(NRMAX-1))
C
      CALL SPL2D(BETAG,ALPHG,ZPS,ZPSB,ZPSA,ZPSAB,UABZ,
     &           NTHM,NTHMAX,NRMAX,0,0,IERR)
      IF(IERR.NE.0) THEN
         IERR=IERR+20000
         RETURN
      ENDIF
C
      RETURN
      END
C
C     ***** CONVERSION (R,Z) to (Alpha,Beta) *****
C
      SUBROUTINE RZTOAB(R,Z,ALPHA,BETA)
C
      INCLUDE '../eq/eqcomq.h'
      COMMON /EQAKF1/ ALPHAF1,RF1,ZF1
      EXTERNAL FNBETA
C
      PSIL=PSIG(R,Z)
C      WRITE(6,'(A,1PE12.4)') 'PSIL=',PSIL
C      CALL GUFLSH
      CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
C      WRITE(6,'(A,1PE12.4)') 'FTL=',FTL
C      CALL GUFLSH
      IF(FTL.LT.0.D0) FTL=0.D0
      ALPHA=FTL
      ALPHAF1=ALPHA
      RF1=R
      ZF1=Z
C
      EPSZ=1.D-8
      BETAIN=ATAN2((Z-ZAXIS)/RKAP,R-RAXIS)
C      WRITE(6,'(A,1PE12.4)') 'BETAIN=',BETAIN
C      CALL GUFLSH
      BETA=ZBRENTX(FNBETA,BETAIN-0.5D0,BETAIN+0.5D0,EPSZ)
C      FX=FNBETA(BETA)
C      WRITE(6,'(A,1P2E12.4)') 'BETA,FNBETA=',BETA,FX
C      CALL GUFLSH
      RETURN
      END
C
C     ----- FNBETA=0 FOR BETA CLOSEST TO (RF1,ZF1) FOR FIXED ALPHA
C
      FUNCTION FNBETA(BETA)
C
      INCLUDE '../eq/eqcomq.h'
      COMMON /EQAKF1/ ALPHAF1,RF1,ZF1
C
      DBETA=1.D-6
      ALPHA1=ALPHAF1
      BETA1=BETA+0.5D0*DBETA
      CALL ABTORZ(ALPHA1,BETA1,R1,Z1,IERR)
      FN1=(R1-RF1)**2+(Z1-ZF1)**2
      BETA2=BETA-0.5D0*DBETA
      CALL ABTORZ(ALPHA1,BETA2,R2,Z2,IERR)
      FN2=(R2-RF1)**2+(Z2-ZF1)**2
      FNBETA=(FN2-FN1)/DBETA
C      WRITE(6,'(A,1P5E11.3)') 
C     &     'BETA,R1,Z1,FN1,FNBETA=',BETA,R1,Z1,FN1,FNBETA
C      CALL GUFLSH
      RETURN
      END
C
C     ***** CONVERSION (Alpha,Beta) TO (R,Z) *****
C
      SUBROUTINE ABTORZ(ALPHA,BETA,R,Z,IERR)
C
      INCLUDE '../eq/eqcomq.h'
      COMMON /EQAKV1/ ALPHG(NRM),BETAG(NRM)
      COMMON /EQAKV2/ UABR(4,NTHM,NRM),UABZ(4,NTHM,NRM)
C
      IERR=0
      IF(ALPHA.LE.0.D0) THEN
         ALPH=0.D0
      ELSE
         ALPH=SQRT(ALPHA)
      ENDIF
      CALL SPL2DF(BETA,ALPH,R,
     &            BETAG,ALPHG,UABR,NTHM,NTHMAX,NRMAX,IERR)
      IF(IERR.NE.0) THEN
         IERR=IERR+10000
         RETURN
      ENDIF
      CALL SPL2DF(BETA,ALPH,Z,
     &            BETAG,ALPHG,UABZ,NTHM,NTHMAX,NRMAX,IERR)
      IF(IERR.NE.0) THEN
         IERR=IERR+20000
         RETURN
      ENDIF
      RETURN
      END
C
C     ***** CALCULATE q AND dq/dalpha *****
C
      SUBROUTINE SUBQPA(ALPHA,Q,DQDA)
C
      INCLUDE '../eq/eqcomq.h'
C
      FTL=ALPHA
      CALL SPL1DD(FTL,PSIL,DPSIL,FTS,UFTT,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SUBQPA: SPL1DD ERROR1 : IERR=',IERR
      CALL SPL1DD(PSIL,QPL,DQPL,PSS,UQPS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SUBQPA: SPL1DD ERROR2 : IERR=',IERR
      Q=QPL
      IF(DPSIL.EQ.0.D0) THEN
         DQDA=0.D0
      ELSE
         DQDA=DQPL/DPSIL
      ENDIF
      RETURN
      END
C
C     ***** CALCULATE Bmax AND dBmax/dalpha *****
C
      SUBROUTINE SUBBMX(ALPHA,BMX,DBMXDA)
C
      INCLUDE '../eq/eqcomq.h'
C
      FTL=ALPHA
      CALL SPL1DD(FTL,PSIL,DPSIL,FTS,UFTT,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SUBBMX: SPL1DD ERROR1 : IERR=',IERR
      CALL SPL1DD(PSIL,BBMAXL,DBBMAXL,PSS,UBBMAX,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SUBBMX: SPL1DD ERROR2 : IERR=',IERR
      BMX=BBMAXL
      IF(DPSIL.EQ.0.D0) THEN
         DBMXDA=0.D0
      ELSE
         DBMXDA=DBBMAXL/DPSIL
      ENDIF
      RETURN
      END
C
C     ***** CALCULATE Bmin AND dBmin/dalpha *****
C
      SUBROUTINE SUBBMN(ALPHA,BMN,DBMNDA)
C
      INCLUDE '../eq/eqcomq.h'
C
      FTL=ALPHA
      CALL SPL1DD(FTL,PSIL,DPSIL,FTS,UFTT,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SUBBMN: SPL1DD ERROR1 : IERR=',IERR
      CALL SPL1DD(PSIL,BBMINL,DBBMINL,PSS,UBBMIN,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SUBBMN: SPL1DD ERROR2 : IERR=',IERR
      BMN=BBMINL
      IF(DPSIL.EQ.0.D0) THEN
         DBMNDA=0.D0
      ELSE
         DBMNDA=DBBMINL/DPSIL
      ENDIF
      RETURN
      END
C
C     ***** CALCULATE MAGNETIC FIELD *****
C
      SUBROUTINE SUBMAG(ALPHA,BETA,BR,BZ,BT,BTOT)
C
      INCLUDE '../eq/eqcomq.h'
C
      CALL ABTORZ(ALPHA,BETA,RP,ZP,IERR)
      IF(IERR.NE.0) WRITE(6,*) 
     &        'XX BTOTAB: ABTORZ ERROR : IERR=',IERR
C
      CALL SPL2DD(RP,ZP,PSIL,PSIR,PSIZ,
     &            RG,ZG,URZ,NRGM,NRGMAX,NZGMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 
     &        'XX BTOTAB: SPL1DF ERROR : IERR=',IERR
C
      CALL SPL1DF(PSIL,TTL,PSS,UTTS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 
     &        'XX BTOTAB: SPL1DF ERROR : IERR=',IERR
C
      BT= TTL/RP
      BR=-PSIZ/RP
      BZ= PSIR/RP
      BTOT=SQRT(BT**2+BR**2+BZ**2)
      RETURN
      END
C
C     ***** GET MAGNETIC AXIS *****
C
      SUBROUTINE GETAXS(RAXIS1,ZAXIS1)
C
      INCLUDE '../eq/eqcomq.h'
C
      RAXIS1=RAXIS
      ZAXIS1=ZAXIS
      RETURN
      END
C
C     ***** GET PLASMA BOUNDARY POSITION *****
C
      SUBROUTINE GETRSU(RSU1,ZSU1,N,NSUMAX1)
C
      INCLUDE '../eq/eqcomq.h'
      DIMENSION RSU1(N),ZSU1(N)
C
      NSUMAX1=NSUMAX
      DO NSU=1,MIN(N,NSUMAX)
         RSU1(NSU)=RSU(NSU)
         ZSU1(NSU)=ZSU(NSU)
      ENDDO
      RETURN
      END
C
      FUNCTION ZBRENTX(FUNC,X1,X2,TOL)
C
      INTEGER ITMAX
      REAL*8 ZBRENTX,TOL,X1,X2,FUNC,EPS
      EXTERNAL FUNC
      PARAMETER (ITMAX=100,EPS=1.D-15)
      INTEGER ITER
      REAL*8 A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF((FA.GT.0..AND.FB.GT.0.).OR.(FA.LT.0..AND.FB.LT.0.)) THEN
         WRITE(6,'(A,1P3E12.4)') 'XX ZBRENT: ROOT MUST BE BETWEEN',X1,X2
         ZBRENTX=A
         RETURN
      ENDIF
      C=B
      FC=FB
      DO 11 ITER=1,ITMAX
        IF((FB.GT.0..AND.FC.GT.0.).OR.(FB.LT.0..AND.FC.LT.0.))THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.D0*EPS*ABS(B)+0.5D0*TOL
        XM=0.5D0*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.D0)THEN
          ZBRENTX=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.D0*XM*S
            Q=1.D0-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.D0*XM*Q*(Q-R)-(B-A)*(R-1.D0))
            Q=(Q-1.D0)*(R-1.D0)*(S-1.D0)
          ENDIF
          IF(P.GT.0.D0) Q=-Q
          P=ABS(P)
          IF(2.D0*P .LT. MIN(3.D0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      WRITE(6,*) 'ZBRENT EXCEEDING MAXIMUM ITERATIONS'
      ZBRENTX=B
      RETURN
      END
