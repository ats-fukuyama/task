C     $Id$
C
C     **********************************
C
C     CROSS SECTION OF ATOMIC PROCESSES
C
C     **********************************
C
C     ***** INITIALIZE ROUTINE *****
C
      SUBROUTINE ATINIT(KID,NN)
C
      USE libspl1d
      IMPLICIT REAL*8(A-B,D-F,H,O-Z)
      PARAMETER (NNM=3)
      PARAMETER (NENM=45)
      PARAMETER (NTEM=101)
      COMMON /ATCOM1/ ENLA(NENM,NNM),SIGM(NENM,NNM),NENMAX(NNM)
      COMMON /ATCOM2/ TELA(NTEM,NNM),SIGV(NTEM,NNM),NTEMAX
      COMMON /ATCOM3/ UEN(4,NENM,NNM)
      COMMON /ATCOM4/ UTE(4,NTEM,NNM)
      COMMON /ATCOM5/ TELDE,NNDE
      DIMENSION FX1(NENM),FX2(NTEM)
      CHARACTER KID*4
C
      EXTERNAL ATFUNC
C
      IF(NN.GT.NNM) GOTO 9999
C
      PI  = ASIN(1.D0)*2.D0
      AME = 9.1093897D-31
      AEE = 1.60217733D-19
C
      IF (KID.EQ.'AR  ') THEN
         CALL DATAAR(ENLA(1,NN),SIGM(1,NN),NENM,NENMAX(NN))
      ELSEIF (KID.EQ.'H2  ') THEN
         CALL DATAH2(ENLA(1,NN),SIGM(1,NN),NENM,NENMAX(NN))
      ELSEIF (KID.EQ.'CF4 ') THEN
         CALL DATACF4(ENLA(1,NN),SIGM(1,NN),NENM,NENMAX(NN))
      ELSE
         WRITE(6,*) 'XX ATINIT: UNKNOWN KID: KID=',KID
         GO TO 9000
      ENDIF
C
      DO NEN=1,NENMAX(NN)
         FX1(NEN)=0.D0
      ENDDO
C
      CAll SPL1D(ENLA(:,NN),SIGM(:,NN),FX1,
     &           UEN(:,:,NN),NENMAX(NN),0,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX ATINIT: SPLINE 1 FAILED: IERR=',IERR
         GOTO 9000
      ENDIF
C
      TELMIN=-4.D0
      TELMAX= 4.D0
      NTEMAX= 101
C
      DTEL=(TELMAX-TELMIN)/(NTEMAX-1)
      DO NTE=1,NTEMAX
         TELA(NTE,NN)=TELMIN+DTEL*(NTE-1)
      ENDDO
C
      H0=0.5D0
      EPS=1.D-8
      ILST=0
C
      DO NTE=1,NTEMAX
         TELDE=TELA(NTE,NN)
         NNDE=NN
         CALL DEHIFE(ATFUNC,FS,ES,H0,EPS,ILST)
C
C         IF(ILST.EQ.0) WRITE(6,602) CS,ES
C  602    FORMAT(1H ,13X,16X,1PD24.15,1PD14.5)
C
         TE=10.D0**TELA(NTE,NN)
         FACT=DSQRT(8.D0*TE*AEE/(PI*AME))
C
         SIGV(NTE,NN)=LOG10(FS*FACT)
      ENDDO
C
      DO NTE=1,NTEMAX
         FX2(NTE)=0.D0
      ENDDO
C
      CAll SPL1D(TELA(1,NN),SIGV(1,NN),FX2,
     &           UTE(:,:,NN),NTEMAX,0,IERR)
C
 9000 RETURN
C
 9999 WRITE(6,*) 'XX ATINIT: NN.GT.NNM: NN,NNM=',NN,NNM
      RETURN
      END
C
C     ***** FUNCTION SIGMA V *****
C
      REAL*8 FUNCTION ATFUNC(Z)
C
      IMPLICIT REAL*8(A-B,D-F,H,O-Z)
      COMMON /ATCOM5/ TELDE,NNDE
C
      TE=10.D0**TELDE
      EN=TE*Z
      NN=NNDE
      CALL ATSIGM(EN,SIGMA,NN)
      ATFUNC=SIGMA*Z*EXP(-Z)
      RETURN
      END
C
C     ***** INETRPOLATE <SIGMA V> *****
C
      SUBROUTINE ATSIGV(TE,SIGVA,NN)
C
      USE libspl1d
      IMPLICIT REAL*8(A-B,D-F,H,O-Z)
      PARAMETER (NNM=3)
      PARAMETER (NTEM=101)
      COMMON /ATCOM2/ TELA(NTEM,NNM),SIGV(NTEM,NNM),NTEMAX
      COMMON /ATCOM4/ UTE(4,NTEM,NNM)
C
      TEL=DLOG10(TE)
      IF(TEL.LT.TELA(1,NN)) THEN
         SIGVAL=SIGV(1,NN)
      ELSEIF(TEL.GT.TELA(NTEMAX,NN)) THEN
         SIGVAL=SIGV(NTEMAX,NN)
      ELSE
         CALL SPL1DX(TEL,SIGVAL,TELA(:,NN),UTE(:,:,NN),NTEMAX,IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX ATSIGM: IERR=',IERR
         ENDIF
      ENDIF
C
      SIGVA=10.D0**SIGVAL
      RETURN
      END
C
C     ***** INTERPOLATE SIGMA *****
C
      SUBROUTINE ATSIGM(EN,SIGMA,NN)
C
      USE libspl1d
      IMPLICIT REAL*8(A-B,D-F,H,O-Z)
      PARAMETER (NNM=3)
      PARAMETER (NENM=45)
      COMMON /ATCOM1/ ENLA(NENM,NNM),SIGM(NENM,NNM),NENMAX(NNM)
      COMMON /ATCOM3/ UEN(4,NENM,NNM)
C
      ENL=DLOG10(EN)
      IF(ENL.LT.ENLA(1,NN)) THEN
         SIGMAL=SIGM(1,NN)
      ELSEIF(ENL.GT.ENLA(NENMAX(NN),NN)) THEN
         SIGMAL=SIGM(NENMAX(NN),NN)
      ELSE
         CALL SPL1DX(ENL,SIGMAL,ENLA(:,NN),
     &               UEN(:,:,NN),NENMAX(NN),IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX ATSIGM: IERR=',IERR
         ENDIF
      ENDIF
C
      SIGMA=10.D0**SIGMAL
      RETURN
      END
C
C     ***** INQUIRE ENERGY RANGE OF DATA *****
C
      SUBROUTINE ATINQR(NN,ENMIN,ENMAX,TEMIN,TEMAX)
C
      IMPLICIT REAL*8(A-B,D-F,H,O-Z)
      PARAMETER (NNM=3)
      PARAMETER (NENM=45)
      PARAMETER (NTEM=101)
      COMMON /ATCOM1/ ENLA(NENM,NNM),SIGM(NENM,NNM),NENMAX(NNM)
      COMMON /ATCOM2/ TELA(NTEM,NNM),SIGV(NTEM,NNM),NTEMAX
C
      NEN=NENMAX(NN)
      ENMIN=10.D0**ENLA(1,NN)
      ENMAX=10.D0**ENLA(NEN,NN)
      TEMIN=10.D0**TELA(1,NN)
      TEMAX=10.D0**TELA(NTEMAX,NN)
      RETURN
      END
C
C     ***** CROSS SECTION DATA FOR Ar *****
C
      SUBROUTINE DATAAR(X,Y,NENM1,NENMAX)
C
      IMPLICIT REAL*8(A-B,D-F,H,O-Z)
      PARAMETER (NENM=45)
      DIMENSION ENT(NENM)
      DIMENSION SIGMAT(NENM)
      DIMENSION X(NENM1),Y(NENM1)
C
      DATA ENT/
     &     1.00D-2, 3.00D-2, 1.00D-1, 3.00D-1, 1.00D+0,
     &     1.20D+0, 1.50D+0, 2.00D+0, 2.50D+0, 3.00D+0,
     &     4.00D+0, 5.00D+0, 6.00D+0, 8.00D+0, 1.00D+1,
     &     1.20D+1, 1.50D+1, 2.00D+1, 2.50D+1, 3.00D+1,
     &     4.00D+1, 5.00D+1, 6.00D+1, 8.00D+1, 1.00D+2,
     &     1.20D+2, 1.50D+2, 2.00D+2, 2.50D+2, 3.00D+2,
     &     4.00D+2, 5.00D+2, 6.00D+2, 8.00D+2, 1.00D+3,
     &     1.20D+3, 1.50D+3, 2.00D+3, 2.50D+3, 3.00D+3,
     &     4.00D+3, 5.00D+3, 6.00D+3, 8.00D+3, 1.00D+4 /
C
      DATA SIGMAT/
     &     4.50D-16, 2.57D-16, 6.20D-17, 1.51D-17, 1.37D-16,
     &     1.66D-16, 2.05D-16, 2.70D-16, 3.43D-16, 4.20D-16,
     &     5.70D-16, 7.50D-16, 9.00D-16, 1.20D-15, 1.40D-15,
     &     1.39D-15, 1.24D-15, 9.20D-16, 6.90D-16, 5.60D-16, 
     &     4.10D-16, 3.15D-16, 2.65D-16, 2.00D-16, 1.62D-16,
     &     1.36D-16, 1.10D-16, 8.70D-17, 7.10D-17, 6.00D-17,
     &     4.60D-17, 3.70D-17, 3.10D-17, 2.25D-17, 1.70D-17,
     &     1.35D-17, 9.80D-18, 6.50D-18, 4.65D-18, 3.50D-18,
     &     2.20D-18, 1.52D-18, 1.13D-18, 7.00D-19, 4.90D-19 /
C
      NENMAX=MIN(NENM,NENM1)
      DO NEN=1,NENMAX
         X(NEN)=DLOG10(ENT(NEN))
         Y(NEN)=DLOG10(SIGMAT(NEN)*1.0D-4)
      ENDDO
C
      RETURN
      END
C
C     ***** CROSS SECTION DATA FOR H2 *****
C
      SUBROUTINE DATAH2(X,Y,NENM1,NENMAX)
C
      IMPLICIT REAL*8(A-B,D-F,H,O-Z)
      PARAMETER (NENM=37)
      DIMENSION ENT(NENM)
      DIMENSION SIGMAT(NENM)
      DIMENSION X(NENM1),Y(NENM1)
C
      DATA ENT/
     &     1.00D-3, 3.00D-3, 1.00D-2, 3.00D-2, 1.00D-1,
     &     3.00D-1, 1.00D+0, 1.20D+0, 1.50D+0, 2.00D+0,
     &     2.50D+0, 3.00D+0, 4.00D+0, 5.00D+0, 6.00D+0,
     &     8.00D+0, 1.00D+1, 1.20D+1, 1.50D+1, 2.00D+1,
     &     2.50D+1, 3.00D+1, 4.00D+1, 5.00D+1, 6.00D+1,
     &     8.00D+1, 1.00D+2, 1.20D+2, 1.50D+2, 2.00D+2,
     &     2.50D+2, 3.00D+2, 4.00D+2, 5.00D+2, 6.00D+2,
     &     8.00D+2, 1.00D+3 /
C
      DATA SIGMAT/
     &     6.55D-16, 6.75D-16, 7.30D-16, 8.50D-16, 1.05D-15,
     &     1.30D-15, 1.74D-15, 1.80D-15, 1.82D-15, 1.80D-15,
     &     1.69D-15, 1.54D-15, 1.30D-15, 1.06D-15, 8.90D-16,
     &     6.50D-16, 5.00D-16, 3.90D-16, 2.90D-16, 1.88D-16,
     &     1.35D-16, 1.00D-16, 6.20D-17, 4.40D-17, 3.10D-17,
     &     1.95D-17, 1.32D-17, 1.04D-17, 7.30D-18, 4.50D-18,
     &     3.05D-18, 2.20D-18, 1.35D-18, 9.10D-19, 6.60D-19,
     &     4.00D-19, 2.70D-19 /
C
      NENMAX=MIN(NENM,NENM1)
      DO NEN=1,NENMAX
         X(NEN)=DLOG10(ENT(NEN))
         Y(NEN)=DLOG10(SIGMAT(NEN)*1.0D-4)
      ENDDO
C
      RETURN
      END
C
C     ***** CROSS SECTION DATA FOR CF4 *****
C
      SUBROUTINE DATACF4(X,Y,NENM1,NENMAX)
C
      IMPLICIT REAL*8(A-B,D-F,H,O-Z)
      PARAMETER (NENM=11)
      DIMENSION ENT(NENM)
      DIMENSION SIGMAT(NENM)
      DIMENSION X(NENM1),Y(NENM1)
C
      DATA ENT/
     &     2.00D-2, 5.00D-2, 1.00D-1, 2.00D-1, 5.00D-1,
     &     1.00D+0, 5.00D+0, 1.00D+1, 3.00D+1, 1.00D+2,
     &     2.00D+2 /

      DATA SIGMAT/
     &     3.31D-16, 1.25D-16, 8.2D-17, 5.5D-17, 4.9D-17,
     &     9.2D-17, 5.6D-16, 9.9D-16, 9.2D-16, 5.4D-16,
     &     4.6D-16 /
C
      NENMAX=MIN(NENM,NENM1)
      DO NEN=1,NENMAX
         X(NEN)=DLOG10(ENT(NEN))
         Y(NEN)=DLOG10(SIGMAT(NEN)*1.0D-4)
      ENDDO
C
      RETURN
      END
C
C     ****** One-Dimensional Spline Interpolation ******
C             **** arbitrary data interval  ****
C
      SUBROUTINE SPL1DX(X0,F0,X,U,NXMAX,IERR)
C
      IMPLICIT REAL*8(A-B,D-F,H,O-Z)
C
      DIMENSION X(NXMAX)
      DIMENSION U(4,NXMAX)
C
      IF(X0.LT.X(1)) GOTO 9001
      IF(X0.GT.X(NXMAX)) GOTO 9002
C
      NX=INT((X0-X(1))/(X(NXMAX)-X(1))*(NXMAX-1)-1.D-14)+2
c      write(6,*) x0, x(nx-1), x(nx), nx
      IF(NX.LT.2) NX=2
      IF(NX.GT.NXMAX) NX=NXMAX
C
 10   IF(X0.LT.X(NX-1)) THEN
         NX=NX-1
         GOTO 10
      ENDIF
      IF(X0.GT.X(NX)) THEN
         NX=NX+1
         GOTO 10
      ENDIF
      CONTINUE
C
      DX=X0-X(NX-1)
C
      F0= U(1,NX)
     &  + U(2,NX)*DX
     &  + U(3,NX)*DX*DX
     &  + U(4,NX)*DX*DX*DX
      IERR=0
      RETURN
C
 9001 IERR=1
      RETURN
 9002 IERR=2
      RETURN
      END
C
C****************************************************************
C
C      HALF INFINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA 
C               FOR INTEGRAND WITH FACTOR EXP(-X)
C                    (0, +INFINITE)
C         INTEGRAND SHOULD BE DEFINED BY FUNC(X)
C
C            INPUT:  H0:   INITIAL STEP SIZE
C                    EPS:  CONVERGENCE CRITERION
C                    ILST: 0 : NO PRINT OUT
C                          1 : PRINT FINAL RESULT
C                          2 : PRINT INTERMEDIATE RESULT
C                    FUNC: INTEGRAND FUNCTION
C            OUTPUT: CS:   INTEGRAL
C                    ES:   ESTIMATED ERROR
C
C****************************************************************
C
      SUBROUTINE DEHIFE(FUNC,CS,ES,H0,EPS,ILST)
      IMPLICIT REAL*8(A-F,H,O-Z)
      CHARACTER KL*1
      EXTERNAL FUNC
C
      EPS1=EPS**0.75D0
      H=H0
      X=EXP(-1.D0)
      CSI=2.D0*X*FUNC(X)
      CS=H*CSI
      CSP=0.D0
      NP=0
      NM=0
      NPMIN=1
      NMMIN=1
C
    5 IND=0
      ATP=ABS(CSI)
      ATM=ATP 
      NPD=2
      IF(NP.EQ.0) NPD=1
      NMD=2
      IF(NM.EQ.0) NMD=1
C
   10 IF(IND.NE.1) THEN
         IF(NP.EQ.NPMIN+2) NPD=1
         NP=NP+NPD
         HN=DBLE(NP)*H
         HS=EXP(-HN)
         X=EXP( HN-HS)
         CT=H*(1.D0+HS)*X*FUNC(X)
         CS=CS+CT
         AT=ATP
         ATP=ABS(CT)/H
         IF(NP.GE.NPMIN) THEN
            IF(AT+ATP.LE.EPS1*ABS(CS)) THEN
               IF(IND.EQ.-1) GO TO 100
               IND=1
            ENDIF
         ENDIF
      ENDIF
C     
      IF(IND.NE.-1) THEN
         IF(NM.EQ.NMMIN+2) NMD=1
         NM=NM+NMD
         HN=DBLE(NM)*H
         HS=EXP( HN)
         X=EXP(-HN-HS)
         CT=H*(1.D0+HS)*X*FUNC(X)
         CS=CS+CT
         AT=ATM
         ATM=ABS(CT)/H
         IF(NM.GE.NMMIN) THEN
            IF(AT+ATM.LE.EPS1*ABS(CS)) THEN
               IF(IND.EQ.1) GO TO 100
               IND=-1
            ENDIF
         ENDIF
      ENDIF
      GO TO 10
C
  100 ES=ABS(CS-CSP)
      IF(ILST.EQ.2) THEN
         IF(H.GE.H0) WRITE(6,601) H,NP,NM,CS
         IF(H.LT.H0) WRITE(6,602) H,NP,NM,CS,ES
      ENDIF
      CSP=CS
      IF(ES.LT.EPS1*ABS(CS)) GO TO 200
      NMAX=MAX0(NP,NM)
      IF(NMAX.GT.1000) THEN
  110    WRITE(6,603)
         READ(5,501,END=120) KL
         IF(KL.EQ.'Q') GOTO 200
         IF(KL.NE.'C') GOTO 110
      ENDIF
  120 CONTINUE
      H=0.5D0*H
      CS=0.5D0*CS
      NPMIN=NP*2-1
      NMMIN=NM*2-1
      NP=-1
      NM=-1
      GO TO 5
C
  200 IF(ILST.EQ.1) THEN
         WRITE(6,602) H,NP,NM,CS,ES
      ENDIF
      RETURN
C
  501 FORMAT(A1)
  601 FORMAT(1H ,1PD13.5,2I8,1PD24.15)
  602 FORMAT(1H ,1PD13.5,2I8,1PD24.15,1PD14.5)
  603 FORMAT(1H ,'# SUB DEHIFE # C or CR : CONTINUE / Q : QUIT')
      END
