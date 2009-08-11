C     $Id$
C
C     ************************************************************
C        FINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
C                    (-1.D0, +1.D0)
C         INTEGRAND SHOULD BE DEFINED BY FUNC(X,1-X,1+X)
C     ************************************************************
C
C        For integral with respect to y from a to b,
C            variable transformation should be
C                 x = (2*y-a-b)/(b-a)
C
C                 near y ~ a   y-a=(b-a)*(1+x)/2
C                 near y ~ b   b-y=(b-a)*(1-x)/2
C                 otherwize      y=(b-a)*x/2+(a+b)/2
C
      SUBROUTINE DEFT(CS,ES,H0,EPS,ILST,FUNC,KID)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUNC
      CHARACTER*(*) KID
      DATA HP/1.5707 96326 79489 66192D0/
C
      EPS1=EPS**0.75D0
      H=H0
      X=0.D0
      CSI=HP*FUNC(X,1.D0-X,1.D0+X)
      CS=H*CSI
      CSP=0.D0
      NP=0
      NM=0
      NPMIN=1
      NMMIN=1
C
    1 IND=0
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
         HC=HP*H*COSH(HN)
         HS=HP*SINH(-HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=CC*EXP(-HS)
         XP=CC*EXP( HS)
         CT=HC*FUNC(X,XM,XP)*CC*CC
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
         HC=HP*H*COSH(HN)
         HS=HP*SINH( HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=CC*EXP(-HS)
         XP=CC*EXP( HS)
         CT=HC*FUNC(X,XM,XP)*CC*CC
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
      IF(ILST.NE.0) THEN
         IF(H.GE.H0) WRITE(6,601) H,NP,NM,CS
         IF(H.LT.H0) WRITE(6,602) H,NP,NM,CS,ES
      ENDIF
      CSP=CS
      IF(ES.LE.EPS1*ABS(CS)) GO TO 200
      NMAX=MAX0(NP,NM)
      IF(NMAX.GT.1000) THEN
         WRITE(6,603) TRIM(KID)
         GO TO 9999
      ENDIF
      H=0.5D0*H
      CS=0.5D0*CS
      NPMIN=NP*2-1
      NMMIN=NM*2-1
      NP=-1
      NM=-1
      GO TO 1
C
  200 RETURN
C
  501 FORMAT(A1)
  601 FORMAT(1H ,1PD13.5,2I8,1PD24.15)
  602 FORMAT(1H ,1PD13.5,2I8,1PD24.15,1PD14.5)
  603 FORMAT(1H ,'XX DEFT: NMAX EXCEEDS 1000: FUNC=',A)
 9999 STOP
      END
C
C     *************************************************************
C      HALF INFINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
C                    (0, +INFINITE)
C         INTEGRAND SHOULD BE DEFINED BY FUNC(X)
C     *************************************************************
C
      SUBROUTINE DEHIFT(CS,ES,H0,EPS,ILST,FUNC,KID)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUNC
      CHARACTER*(*) KID
      DATA HP/1.5707 96326 79489 66192D0/
C
      EPS1=EPS**0.75D0
      H=H0
      X=1.D0
      CSI=HP*FUNC(X)
      CS=H*CSI
      CSP=0.D0
      NP=0
      NM=0
      NPMIN=1
      NMMIN=1
C
    1 IND=0
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
         HC=HP*H*COSH(HN)
         HS=HP*SINH(-HN)
         X=EXP(HS)
         CT=HC*X*FUNC(X)
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
         HC=HP*H*COSH(HN)
         HS=HP*SINH( HN)
         X=EXP(HS)
         CT=HC*X*FUNC(X)
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
      IF(ILST.NE.0) THEN
         IF(H.GE.H0) WRITE(6,601) H,NP,NM,CS
         IF(H.LT.H0) WRITE(6,602) H,NP,NM,CS,ES
      ENDIF
      CSP=CS
      IF(ES.LE.EPS1*ABS(CS)) GO TO 200
      NMAX=MAX0(NP,NM)
      IF(NMAX.GT.1000) THEN
         WRITE(6,603) TRIM(KID)
         GO TO 9999
      ENDIF
      H=0.5D0*H
      CS=0.5D0*CS
      NPMIN=NP*2-1
      NMMIN=NM*2-1
      NP=-1
      NM=-1
      GO TO 1
C
  200 RETURN
C
  501 FORMAT(A1)
  601 FORMAT(1H ,1PD13.5,2I8,1PD24.15)
  602 FORMAT(1H ,1PD13.5,2I8,1PD24.15,1PD14.5)
  603 FORMAT(1H ,'XX DEHIFT: NMAX EXCEEDS 1000: FUNC=',A)
 9999 STOP
      END
C
C     *************************************************************
C      HALF INFINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA 
C               FOR INTEGRAND WITH FACTOR EXP(-X)
C                    (0, +INFINITE)
C         INTEGRAND SHOULD BE DEFINED BY FUNC(X)
C     *************************************************************
C
      SUBROUTINE DEHIFE(CS,ES,H0,EPS,ILST,FUNC,KID)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*(*) KID
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
    1 IND=0
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
         WRITE(6,603) TRIM(KID)
         GOTO 9999
      ENDIF
      H=0.5D0*H
      CS=0.5D0*CS
      NPMIN=NP*2-1
      NMMIN=NM*2-1
      NP=-1
      NM=-1
      GO TO 1
C
  200 IF(ILST.EQ.1) THEN
         WRITE(6,602) H,NP,NM,CS,ES
      ENDIF
      RETURN
C
  501 FORMAT(A1)
  601 FORMAT(1H ,1PD13.5,2I8,1PD24.15)
  602 FORMAT(1H ,1PD13.5,2I8,1PD24.15,1PD14.5)
  603 FORMAT(1H ,'XX DEHIFE: NMAX EXCEEDS 1000: FUNC=',A)
 9999 STOP
      END
C
C     *************************************************************
C      HALF INFINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA 
C               FOR INTEGRAND WITH FACTOR EXP(-X)
C                    (0, +INFINITE)
C         INTEGRAND SHOULD BE DEFINED BY FUNC(X)
C     *************************************************************
C
      SUBROUTINE DEHIFEC(CS,ES,H0,EPS,ILST,CFUNC,KID)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 CFUNC,CS,CSI,CSP,CT
      CHARACTER*(*) KID
      EXTERNAL CFUNC
C
      EPS1=EPS**0.75D0
      H=H0
      X=EXP(-1.D0)
      CSI=2.D0*X*CFUNC(X)
      CS=H*CSI
      CSP=0.D0
      NP=0
      NM=0
      NPMIN=1
      NMMIN=1
C
    1 IND=0
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
         CT=H*(1.D0+HS)*X*CFUNC(X)
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
         CT=H*(1.D0+HS)*X*CFUNC(X)
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
         WRITE(6,603) TRIM(KID)
         GOTO 9999
      ENDIF
      H=0.5D0*H
      CS=0.5D0*CS
      NPMIN=NP*2-1
      NMMIN=NM*2-1
      NP=-1
      NM=-1
      GO TO 1
C
  200 IF(ILST.EQ.1) THEN
         WRITE(6,602) H,NP,NM,CS,ES
      ENDIF
      RETURN
C
  501 FORMAT(A1)
  601 FORMAT(1H ,1PD12.4,2I8,1P2D12.4)
  602 FORMAT(1H ,1PD12.4,2I8,1P2D12.4,1PD12.4)
  603 FORMAT(1H ,'XX DEHIFEC: NMAX EXCEEDS 1000: FUNC=',A)
 9999 STOP
      END
