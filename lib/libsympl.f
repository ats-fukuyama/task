C     $Id$
C
C     ***** SYMPLECTIC METHOD FOR ODE DERIVED FROM HAMILTONIAN *****
C
      SUBROUTINE SYMPLECTIC(Y,H,SUB,NEQMAX,
     &                      NLPMAX,EPS,NLPX,ERROR,IERR)
C
      IMPLICIT NONE
      INTEGER NEQM
      PARAMETER(NEQM=12)
      INTEGER NEQMAX,NLPMAX,IERR
      REAL*8 Y(NEQM),H,EPS,ERROR
      REAL*8 YA(NEQM),Y1(NEQM),Y2(NEQM),Y1N(NEQM),Y2N(NEQM)
      REAL*8 F1(NEQM),F2(NEQM)
      REAL*8 A11,A12,A21,A22,DY1,DY2
      INTEGER NEQ,NLP,NLPX
      EXTERNAL SUB
C
      IF(NEQMAX.GT.NEQM) THEN
         WRITE(6,*) 'XX SYMPLECTIC: NEQMAX > NEQM'
         IERR=2
         RETURN
      ENDIF
C
      A11=0.25D0
      A12=0.25D0-SQRT(3.D0)/6.D0
      A21=0.25D0+SQRT(3.D0)/6.D0
      A22=0.25D0
C
      DO NEQ=1,NEQMAX
         Y1(NEQ)=Y(NEQ)
         Y2(NEQ)=Y(NEQ)
      ENDDO
      CALL SUB(Y1,F1,NEQMAX)
      DO NEQ=1,NEQMAX
         F2(NEQ)=F1(NEQ)
         YA(NEQ)=MAX(ABS(Y(NEQ)),1.D0)
      ENDDO
C
      DO NLP=1,NLPMAX
         ERROR=0.D0
         DO NEQ=1,NEQMAX
            Y1N(NEQ)=Y(NEQ)+H*(A11*F1(NEQ)+A12*F2(NEQ))
            Y2N(NEQ)=Y(NEQ)+H*(A21*F1(NEQ)+A22*F2(NEQ))
            DY1=(Y1N(NEQ)-Y1(NEQ))/YA(NEQ)
            DY2=(Y2N(NEQ)-Y2(NEQ))/YA(NEQ)
            ERROR=DY1*DY1+DY2*DY2
         ENDDO
         ERROR=ERROR/(2*NEQMAX)
         DO NEQ=1,NEQMAX
            Y1(NEQ)=Y1N(NEQ)
            Y2(NEQ)=Y2N(NEQ)
         ENDDO
         CALL SUB(Y1,F1,NEQMAX)
         CALL SUB(Y2,F2,NEQMAX)
         IF(ERROR.LT.EPS) THEN
            NLPX=NLP
            GOTO 9000
         ENDIF
         
      ENDDO
      NLPX=NLP
      WRITE(6,*) 'XX SYMPLECTIC: NLPMAX is not enough'
      IERR=1
      RETURN
C
 9000 DO NEQ=1,NEQMAX
         Y(NEQ)=Y(NEQ)+0.5D0*H*(F1(NEQ)+F2(NEQ))
      ENDDO
      IERR=0
      RETURN
      END
         
