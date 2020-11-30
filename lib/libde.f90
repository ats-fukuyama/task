!
MODULE libde

  USE task_kinds,ONLY: dp
  PRIVATE
  PUBLIC DEFT,DEFTC,DEHIFT,DEHIFTC,DEHIFE,DEHIFEC

CONTAINS
!     ************************************************************
!        FINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
!                    (-1.D0, +1.D0)
!         INTEGRAND SHOULD BE DEFINED BY FUNC(X,1-X,1+X)
!     ************************************************************
!
!        For integral with respect to y from a to b,
!            variable transformation should be
!                 x = (2*y-a-b)/(b-a)
!
!                 near y ~ a   y-a=(b-a)*(1+x)/2
!                 near y ~ b   b-y=(b-a)*(1-x)/2
!                 otherwize      y=(b-a)*x/2+(a+b)/2

  SUBROUTINE DEFT(CS,ES,H0,EPS,ILST,FUNC,KID)
    IMPLICIT NONE
    REAL(dp),INTENT(OUT):: CS   ! Integral
    REAL(dp),INTENT(OUT):: ES   ! Estimated error 
    REAL(dp),INTENT(IN)::  H0   ! Initial step size
    REAL(dp),INTENT(IN)::  EPS  ! Convergence thrshold
    INTEGER,INTENT(IN)::  ILST ! print out control: 0 for no print out
    INTERFACE
       FUNCTION FUNC(X,XM,XP)
         USE task_kinds,ONLY: dp
         REAL(dp):: FUNC
         REAL(dp),INTENT(IN):: X,XM,XP
       END FUNCTION FUNC
    END INTERFACE
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: KID   ! function identifier string
    REAL(dp),PARAMETER:: HP=1.5707963267948966192D0

    REAL(dp):: EPS1,H,X,CSI,CSP,ATP,ATM,HN,HC,HS,CC,XM,XP,CT,AT
    INTEGER:: NP,NM,NPMIN,NMMIN,IND,NPD,NMD,NMAX

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

    1 IND=0
      ATP=ABS(CSI)
      ATM=ATP
      NPD=2
      IF(NP.EQ.0) NPD=1
      NMD=2
      IF(NM.EQ.0) NMD=1

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

  200 RETURN

  501 FORMAT(A1)
  601 FORMAT(1H ,ES13.5,2I8,ES24.15)
  602 FORMAT(1H ,ES13.5,2I8,ES24.15,ES14.5)
  603 FORMAT(1H ,'XX DEFT: NMAX EXCEEDS 1000: FUNC=',A)
 9999 STOP
 END SUBROUTINE DEFT

!     ************************************************************
!        FINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
!                    (-1.D0, +1.D0)
!         INTEGRAND SHOULD BE DEFINED BY FUNC(X,1-X,1+X)
!     ************************************************************
!
!        For integral with respect to y from a to b,
!            variable transformation should be
!                 x = (2*y-a-b)/(b-a)
!
!                 near y ~ a   y-a=(b-a)*(1+x)/2
!                 near y ~ b   b-y=(b-a)*(1-x)/2
!                 otherwize      y=(b-a)*x/2+(a+b)/2

  SUBROUTINE DEFTC(CS,ES,H0,EPS,ILST,CFUNC,KID)
    IMPLICIT NONE
    COMPLEX(dp),INTENT(OUT):: CS   ! Integral
    REAL(dp),INTENT(OUT):: ES   ! Estimated error 
    REAL(dp),INTENT(IN)::  H0   ! Initial step size
    REAL(dp),INTENT(IN)::  EPS  ! Convergence thrshold
    INTEGER,INTENT(IN)::  ILST ! print out control: 0 for no print out
    INTERFACE
       FUNCTION CFUNC(X,XM,XP)
         USE task_kinds,ONLY: dp
         COMPLEX(dp):: CFUNC
         REAL(dp),INTENT(IN):: X,XM,XP
       END FUNCTION CFUNC
    END INTERFACE
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: KID   ! function identifier string
    REAL(dp),PARAMETER:: HP=1.5707963267948966192D0

    COMPLEX(dp):: CSI,CSP,CT
    REAL(dp):: EPS1,H,X,ATP,ATM,HN,HC,HS,CC,XM,XP,AT
    INTEGER:: NP,NM,NPMIN,NMMIN,IND,NPD,NMD,NMAX

      EPS1=EPS**0.75D0
      H=H0
      X=0.D0
      CSI=HP*CFUNC(X,1.D0-X,1.D0+X)
      CS=H*CSI
      CSP=0.D0
      NP=0
      NM=0
      NPMIN=1
      NMMIN=1

    1 IND=0
      ATP=ABS(CSI)
      ATM=ATP
      NPD=2
      IF(NP.EQ.0) NPD=1
      NMD=2
      IF(NM.EQ.0) NMD=1

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
         CT=HC*CFUNC(X,XM,XP)*CC*CC
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
         CT=HC*CFUNC(X,XM,XP)*CC*CC
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

  200 RETURN

  501 FORMAT(A1)
  601 FORMAT(1H ,ES13.5,2I8,ES24.15)
  602 FORMAT(1H ,ES13.5,2I8,ES24.15,ES14.5)
  603 FORMAT(1H ,'XX DEFT: NMAX EXCEEDS 1000: FUNC=',A)
 9999 STOP
 END SUBROUTINE DEFTC

!     *************************************************************
!      HALF INFINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
!                    (0, +INFINITE)
!         INTEGRAND SHOULD BE DEFINED BY FUNC(X)
!     *************************************************************

 SUBROUTINE DEHIFT(CS,ES,H0,EPS,ILST,FUNC,KID)
    IMPLICIT NONE
    REAL(dp),INTENT(OUT):: CS   ! Integral
    REAL(dp),INTENT(OUT):: ES   ! Estimated error 
    REAL(dp),INTENT(IN)::  H0   ! Initial step size
    REAL(dp),INTENT(IN)::  EPS  ! Convergence thrshold
    INTEGER,INTENT(IN)::  ILST ! print out control: 0 for no print out
    INTERFACE
       FUNCTION FUNC(X)
         USE task_kinds,ONLY: dp
         REAL(dp):: FUNC
         REAL(dp),INTENT(IN):: X
       END FUNCTION FUNC
    END INTERFACE
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: KID   ! function identifier string
    REAL(dp),PARAMETER:: HP=1.5707963267948966192D0

    REAL(dp):: EPS1,H,X,CSI,CSP,ATP,ATM,HN,HC,HS,CC,XM,XP,CT,AT
    INTEGER:: NP,NM,NPMIN,NMMIN,IND,NPD,NMD,NMAX

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

    1 IND=0
      ATP=ABS(CSI)
      ATM=ATP
      NPD=2
      IF(NP.EQ.0) NPD=1
      NMD=2
      IF(NM.EQ.0) NMD=1

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

  200 RETURN

  501 FORMAT(A1)
  601 FORMAT(1H ,ES13.5,2I8,ES24.15)
  602 FORMAT(1H ,ES13.5,2I8,ES24.15,ES14.5)
  603 FORMAT(1H ,'XX DEHIFT: NMAX EXCEEDS 1000: FUNC=',A)
 9999 STOP
  END SUBROUTINE DEHIFT

!     *************************************************************
!      HALF INFINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
!                    (0, +INFINITE)
!         INTEGRAND SHOULD BE DEFINED BY FUNC(X)
!     *************************************************************

  SUBROUTINE DEHIFTC(CS,ES,H0,EPS,ILST,CFUNC,KID)
    IMPLICIT NONE
    COMPLEX(dp),INTENT(OUT):: CS   ! Integral
    REAL(dp),INTENT(OUT):: ES   ! Estimated error 
    REAL(dp),INTENT(IN)::  H0   ! Initial step size
    REAL(dp),INTENT(IN)::  EPS  ! Convergence thrshold
    INTEGER,INTENT(IN)::  ILST ! print out control: 0 for no print out
    INTERFACE
       FUNCTION CFUNC(X)
         USE task_kinds,ONLY: dp
         COMPLEX(dp):: CFUNC
         REAL(dp),INTENT(IN):: X
       END FUNCTION CFUNC
    END INTERFACE
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: KID   ! function identifier string
    REAL(dp),PARAMETER:: HP=1.5707963267948966192D0

    COMPLEX(dp):: CSI,CSP,CT
    REAL(dp):: EPS1,H,X,ATP,ATM,HN,HC,HS,XM,XP,AT
    INTEGER:: NP,NM,NPMIN,NMMIN,IND,NPD,NMD,NMAX

      EPS1=EPS**0.75D0
      H=H0
      X=1.D0
      CSI=HP*CFUNC(X)
      CS=H*CSI
      CSP=0.D0
      NP=0
      NM=0
      NPMIN=1
      NMMIN=1

    1 IND=0
      ATP=ABS(CSI)
      ATM=ATP
      NPD=2
      IF(NP.EQ.0) NPD=1
      NMD=2
      IF(NM.EQ.0) NMD=1

   10 IF(IND.NE.1) THEN
         IF(NP.EQ.NPMIN+2) NPD=1
         NP=NP+NPD
         HN=DBLE(NP)*H
         HC=HP*H*COSH(HN)
         HS=HP*SINH(-HN)
         X=EXP(HS)
         CT=HC*X*CFUNC(X)
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

      IF(IND.NE.-1) THEN
         IF(NM.EQ.NMMIN+2) NMD=1
         NM=NM+NMD
         HN=DBLE(NM)*H
         HC=HP*H*COSH(HN)
         HS=HP*SINH( HN)
         X=EXP(HS)
         CT=HC*X*CFUNC(X)
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

  200 RETURN

  501 FORMAT(A1)
  601 FORMAT(1H ,ES13.5,2I8,ES24.15)
  602 FORMAT(1H ,ES13.5,2I8,ES24.15,ES14.5)
  603 FORMAT(1H ,'XX DEHIFT: NMAX EXCEEDS 1000: FUNC=',A)
 9999 STOP
  END SUBROUTINE DEHIFTC

!     *************************************************************
!      HALF INFINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA 
!               FOR INTEGRAND WITH FACTOR EXP(-X)
!                    (0, +INFINITE)
!         INTEGRAND SHOULD BE DEFINED BY FUNC(X)
!     *************************************************************

  SUBROUTINE DEHIFE(CS,ES,H0,EPS,ILST,FUNC,KID)
    IMPLICIT NONE
    REAL(dp),INTENT(OUT):: CS   ! Integral
    REAL(dp),INTENT(OUT):: ES   ! Estimated error 
    REAL(dp),INTENT(IN)::  H0   ! Initial step size
    REAL(dp),INTENT(IN)::  EPS  ! Convergence thrshold
    INTEGER,INTENT(IN)::  ILST ! print out control: 0 for no print out
    INTERFACE
       FUNCTION FUNC(X)
         USE task_kinds,ONLY: dp
         REAL(dp):: FUNC
         REAL(dp),INTENT(IN):: X
       END FUNCTION FUNC
    END INTERFACE
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: KID   ! function identifier string
    REAL(dp),PARAMETER:: HP=1.5707963267948966192D0

    REAL(dp):: EPS1,H,X,CSI,CSP,ATP,ATM,HN,HC,HS,CC,XM,XP,CT,AT
    INTEGER:: NP,NM,NPMIN,NMMIN,IND,NPD,NMD,NMAX

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

    1 IND=0
      ATP=ABS(CSI)
      ATM=ATP 
      NPD=2
      IF(NP.EQ.0) NPD=1
      NMD=2
      IF(NM.EQ.0) NMD=1

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

  200 IF(ILST.EQ.1) THEN
         WRITE(6,602) H,NP,NM,CS,ES
      ENDIF
      RETURN

  501 FORMAT(A1)
  601 FORMAT(1H ,ES13.5,2I8,ES24.15)
  602 FORMAT(1H ,ES13.5,2I8,ES24.15,ES14.5)
  603 FORMAT(1H ,'XX DEHIFE: NMAX EXCEEDS 1000: FUNC=',A)
 9999 STOP
  END SUBROUTINE DEHIFE

!     *************************************************************
!      HALF INFINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA 
!               FOR INTEGRAND WITH FACTOR EXP(-X)
!                    (0, +INFINITE)
!         INTEGRAND SHOULD BE DEFINED BY FUNC(X)
!     *************************************************************

  SUBROUTINE DEHIFEC(CS,ES,H0,EPS,ILST,CFUNC,KID)
    IMPLICIT NONE
    COMPLEX(dp),INTENT(OUT):: CS   ! Integral
    REAL(dp),INTENT(OUT):: ES   ! Estimated error 
    REAL(dp),INTENT(IN)::  H0   ! Initial step size
    REAL(dp),INTENT(IN)::  EPS  ! Convergence thrshold
    INTEGER,INTENT(IN)::  ILST ! print out control: 0 for no print out
    INTERFACE
       FUNCTION CFUNC(X)
         USE task_kinds,ONLY: dp
         COMPLEX(dp):: CFUNC
         REAL(dp),INTENT(IN):: X
       END FUNCTION CFUNC
    END INTERFACE
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: KID   ! function identifier string
    REAL(dp),PARAMETER:: HP=1.5707963267948966192D0

    COMPLEX(dp):: CSI,CSP,CT
    REAL(dp):: EPS1,H,X,ATP,ATM,HN,HC,HS,CC,XM,XP,AT
    INTEGER:: NP,NM,NPMIN,NMMIN,IND,NPD,NMD,NMAX

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

    1 IND=0
      ATP=ABS(CSI)
      ATM=ATP 
      NPD=2
      IF(NP.EQ.0) NPD=1
      NMD=2
      IF(NM.EQ.0) NMD=1

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

  200 IF(ILST.EQ.1) THEN
         WRITE(6,602) H,NP,NM,CS,ES
      ENDIF
      RETURN

  501 FORMAT(A1)
  601 FORMAT(1H ,ES12.4,2I8,2ES12.4)
  602 FORMAT(1H ,ES12.4,2I8,2ES12.4,ES12.4)
  603 FORMAT(1H ,'XX DEHIFEC: NMAX EXCEEDS 1000: FUNC=',A)
 9999 STOP
  END SUBROUTINE DEHIFEC
END MODULE libde
