C     $Id$
C
C     ************************************************************
C        FINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
C                    (-1.D0, +1.D0)
C         INTEGRAND SHOULD BE DEFINED BY FUNC(X,1-X,1+X)
C     ************************************************************
C
      SUBROUTINE DEFT(CS,ES,H0,EPS,ILST,FUNC)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUNC
C     CHARACTER KL*1
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
  110    WRITE(6,603)
         GO TO 9
C        READ(5,501,END=120) KL
C        IF(KL.EQ.'Q') GOTO 200
C        IF(KL.NE.'C') GOTO 110
      ENDIF
  130 H=0.5D0*H
      CS=0.5D0*CS
      NPMIN=NP*2-1
      NMMIN=NM*2-1
      NP=-1
      NM=-1
      GO TO 5
C 120 REWIND 5
C     GO TO 130
C
  200 RETURN
C
  501 FORMAT(A1)
  601 FORMAT(1H ,1PD13.5,2I8,1PD24.15)
  602 FORMAT(1H ,1PD13.5,2I8,1PD24.15,1PD14.5)
  603 FORMAT(1H ,'# SUB DEFT # C or CR : CONTINUE / Q : QUIT')
    9 STOP
      END
C
C     *************************************************************
C      HALF INFINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
C                    (0, +INFINITE)
C         INTEGRAND SHOULD BE DEFINED BY FUNC(X)
C     *************************************************************
C
      SUBROUTINE DEHIFT(CS,ES,H0,EPS,ILST,FUNC)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUNC
C     CHARACTER KL*1
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
  110    WRITE(6,603)
         GO TO 9
C        READ(5,501,END=120) KL
C        IF(KL.EQ.'Q') GOTO 200
C        IF(KL.NE.'C') GOTO 110
      ENDIF
  130 H=0.5D0*H
      CS=0.5D0*CS
      NPMIN=NP*2-1
      NMMIN=NM*2-1
      NP=-1
      NM=-1
      GO TO 5
C 120 REWIND 5
C     GO TO 130
C
  200 RETURN
C
  501 FORMAT(A1)
  601 FORMAT(1H ,1PD13.5,2I8,1PD24.15)
  602 FORMAT(1H ,1PD13.5,2I8,1PD24.15,1PD14.5)
  603 FORMAT(1H ,'# SUB DEHIFT # C or CR : CONTINUE / Q : QUIT')
    9 STOP
      END
