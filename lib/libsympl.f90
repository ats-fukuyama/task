!     ***** SYMPLECTIC METHOD FOR ODE DERIVED FROM HAMILTONIAN *****

      SUBROUTINE SYMPLECTIC(Y,H,SUB,NEQMAX, NLPMAX,EPS,NLPX,ERROR,IERR)

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER, PARAMETER    :: NEQM = 12
      INTEGER, INTENT(IN)   :: NEQMAX,NLPMAX
      INTEGER, INTENT(OUT)  :: IERR
      REAL(dp), INTENT(INOUT), DIMENSION(NEQM) ::  Y
      REAL(dp), INTENT(IN)      :: H, EPS
      REAL(dp), INTENT(OUT)     :: ERROR
      INTEGER, INTENT(OUT)  :: NLPX
      REAL(dp), DIMENSION(NEQM) :: YA, Y1, Y2, Y1N, Y2N, F1, F2
      REAL(dp)    :: A11, A12, A21, A22, DY1, DY2
      INTEGER :: NEQ, NLP
      EXTERNAL SUB

      IF(NEQMAX.GT.NEQM) THEN
         WRITE(6,*) 'XX SYMPLECTIC: NEQMAX > NEQM'
         IERR=2
         RETURN
      ENDIF

      A11=0.25D0
      A12=0.25D0-SQRT(3.D0)/6.D0
      A21=0.25D0+SQRT(3.D0)/6.D0
      A22=0.25D0

      DO NEQ=1,NEQMAX
         Y1(NEQ)=Y(NEQ)
         Y2(NEQ)=Y(NEQ)
      ENDDO
      CALL SUB(Y1,F1,NEQMAX)
      DO NEQ=1,NEQMAX
         F2(NEQ)=F1(NEQ)
         YA(NEQ)=MAX(ABS(Y(NEQ)),1.D0)
      ENDDO

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

 9000 DO NEQ=1,NEQMAX
         Y(NEQ)=Y(NEQ)+0.5D0*H*(F1(NEQ)+F2(NEQ))
      ENDDO
      IERR=0
      RETURN
      END SUBROUTINE SYMPLECTIC
