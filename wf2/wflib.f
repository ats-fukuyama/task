C     $Id$
C
C     ***** REGULA FALSI METHOD *****
C
      SUBROUTINE FRGFLS(XS,XE,DX,XR,FUNC,EPS,ILL)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
C
      EXTERNAL FUNC
      ILMAX=-ILL
      IF(ILMAX.EQ.0) ILMAX=-30
      ILL=0
      IF(ABS(DX).LE.1.E-15) THEN
C        WRITE(6,600)
         ILL=900
         RETURN
      ENDIF
      IF((XE-XS)*DX.LE.0.) THEN
C        WRITE(6,601)
         ILL=901
         RETURN
      ENDIF
      D=DX
C
      X=XS
      Y=FUNC(X)
      IF(ABS(Y).LE.EPS) THEN
         XR=X
         RETURN
      ENDIF
C
   10 IF((XE-X)*D.LE.0.) THEN
C        WRITE(6,602)
         ILL=902
         RETURN
      ENDIF
      XX=X
      YY=Y
      X=XX+D
      Y=FUNC(X)
      IF(ABS(Y).LE.EPS) THEN
         XR=X
         RETURN
      ENDIF
      IF(Y*YY.GT.0.) GOTO 10
C
   20 ILL=ILL-1
      IF(ILL.LE.ILMAX) THEN
         XR=X
C        WRITE(6,603)
         ILL=903
         RETURN
      ENDIF
      IF(Y*YY.LE.0.) THEN
         YYY=Y
         DD=-D
      ELSE
         XX=X
         YY=Y
         DD=D+DD
      ENDIF
      D=DD*YY/(YYY-YY)
      X=XX+D
      Y=FUNC(X)
      IF(ABS(Y).GT.EPS.AND.ABS(D).GT.EPS) GOTO 20
C
      XR=X
      RETURN
C
C 600 FORMAT(1H ,'## FRGFLS : ABS(DX) .LE. 1.E-14 !!')
C 601 FORMAT(1H ,'## FRGFLS : (XE-XS)*DX .LE. 0. !!')
C 602 FORMAT(1H ,'## FRGFLS : NO ROOT BETWEEN XS AND XE !!')
C 603 FORMAT(1H ,'## FRGFLS : NO CONVERGENSE !!')
      END
