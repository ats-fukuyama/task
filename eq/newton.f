C
      SUBROUTINE NEWTON4(N,X0,SUB,STATUS,COUNT,XLAST)
C
C      NEED  SUBROUTINE SUB(X(N),FX(N),N)
C      NEED  M=N
C
      IMPLICIT NONE
      INTEGER M,N,STATUS,COUNT
      INTEGER MAXITS,ITRATE,SOLVED,LIMIT,FLAT,I,J,TOLCT,NWEQM,NWEQM2
      PARAMETER (NWEQM=4,NWEQM2=2*NWEQM)
      REAL*8 TINY,TOL,DFDXMIN
      REAL*8 X0(NWEQM),OLDX(NWEQM),NEWX(NWEQM),XLAST(NWEQM)
      REAL*8 DFDX(NWEQM,NWEQM),INDFDX(NWEQM,NWEQM)
      REAL*8 DX(NWEQM),FX(NWEQM),WORK(NWEQM,NWEQM2)
      REAL*8 WORKS(NWEQM,3)
      PARAMETER (TINY=1E-10)
      PARAMETER (ITRATE=-1,SOLVED=0,LIMIT=1,FLAT=2)
      EXTERNAL SUB
      INTRINSIC ABS
C
C      TOL=5.D-7
      TOL=1.D-4
      MAXITS=100
      DFDXMIN=1.D0
      M=N
C
      DO I=1,N
         NEWX(I)=X0(I)
      END DO 

      STATUS=ITRATE
      DO COUNT=1,MAXITS
C
C        WRITE(*,*) 'QQQQQ',NEWX(1),NEWX(2),NEWX(3),NEWX(4)
C
         CALL NEWTON_SUB(N,NWEQM,NEWX,DFDX,SUB,WORKS)
         DO I=1,4
         DO J=1,4
            WRITE(*,*) 'DDDDDDDD',DFDX(I,J)
         END DO
         END DO 
         CALL INVEMATRIX(N,NWEQM,DFDX,INDFDX,WORK)
         DO I=1,4
         DO J=1,4
            WRITE(*,*) 'BBBBBBBBB',INDFDX(I,J)
         END DO
         END DO
         DO I=1,M
C         DO J=1,N
C            DFDXMIN=MIN(DFDXMIN,ABS(INDFDX(I,J)))
            DFDXMIN=MIN(DFDXMIN,ABS(INDFDX(I,I)))
C         END DO
         END DO

         IF (ABS(DFDXMIN).LE.TINY) THEN
C          {flat spot}
            STATUS=FLAT
         ELSE
C          {perform Newton algorithm}
            DO I=1,M
               OLDX(I)=NEWX(I)
            END DO
            CALL SUB(OLDX,FX,N)
C
            WRITE(*,*) 'BBBBB',FX(1)
            DO I=1,M
               DX(I)=0.D0
            END DO
            DO I=1,M
            DO J=1,N
               DX(I)=DX(I)-INDFDX(I,J)*FX(J)
            END DO
            END DO
            DO I=1,M
               NEWX(I)=OLDX(I)+DX(I)
            END DO
C
            TOLCT=0
            DO I=1,M
               IF (ABS(DX(I)).LE.ABS(OLDX(I))*TOL) TOLCT=TOLCT+1
            END DO
            IF (ABS(TOLCT-M).LE.TOL) STATUS=SOLVED
         END IF
         IF (STATUS.NE.ITRATE) GO TO 11
      END DO
      STATUS=LIMIT

   11 DO I=1,M
         XLAST(I)=NEWX(I)
      END DO
      WRITE(6,*) 'STATUS=',STATUS
      RETURN
      END
c
c************************************************************
      SUBROUTINE NEWTON_SUB(N,NM,X,DFDX,SUB,WORKS)
C
      IMPLICIT NONE
      INTEGER I,J,M,N,NM
      REAL*8 X(N),DFDX(NM,N)
      REAL*8 WORKS(NM,3)
      REAL*8 DD
      EXTERNAL SUB
C
C     WORKS(I,1) = XNEW(I)
C     WORKS(I,2) = RESULT0(I)
C     WORKS(I,3) = RESULT(I)
C
      DD=1.D-2
C
      CALL SUB(X,WORKS(1,2),N)
      DO J=1,N
         DO I=1,N
            WORKS(I,1)=X(I)
         ENDDO
         WORKS(J,1)=X(J)+DD
         CALL SUB(WORKS(1,1),WORKS(1,3),N)
         WORKS(J,1)=X(J)
         DO I=1,N
            DFDX(I,J)=(WORKS(I,3)-WORKS(I,2))/DD
            WRITE(*,*) 'NNNNNNN',DFDX(I,J)
         END DO
      END DO 
      RETURN
      END

C
C******************************************************
      SUBROUTINE FUNC(X,F,N)
C
      IMPLICIT NONE
      INTEGER N
      REAL*8 X,F
      DIMENSION X(4),F(4)
C
      F(1)=X(1)**2+2*X(2)-5*X(1)+1
      F(2)=3*X(2)**2-4*X(1)**2+X(2)-2
      F(3)=X(3)**2+2*X(4)-5*X(3)+1
      F(4)=3*X(4)**2-4*X(3)**2+X(4)-2
C
      WRITE(6,'(A,1P4E12.4)') 'X:',X(1),X(2),X(3),X(4)
      WRITE(6,'(A,1P4E12.4)') 'F:',F(1),F(2),F(3),F(4)
      RETURN
      END
