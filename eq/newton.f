C
      SUBROUTINE NEWTON4(N,X0,XLAST,SUB,STATUS,COUNT,DD,TOL,MAXITS)
C
C      NEED  SUBROUTINE SUB(X(N),FX(N),N)
C      NEED  M=N
C
      IMPLICIT NONE
      INTEGER M,N,STATUS,COUNT,IERR
      INTEGER MAXITS,ITRATE,SOLVED,LIMIT,FLAT,ERROR
      INTEGER I,J,NWEQM,NWEQM2
      PARAMETER (NWEQM=4,NWEQM2=2*NWEQM)
      REAL*8 TINY,TOL,DFDXMIN,SUM,DD,DDTEMP
      REAL*8 X0(NWEQM),OLDX(NWEQM),NEWX(NWEQM),XLAST(NWEQM)
      REAL*8 DFDX(NWEQM,NWEQM),INDFDX(NWEQM,NWEQM)
      REAL*8 DX(NWEQM),FX(NWEQM),WORK(NWEQM,NWEQM2)
      REAL*8 WORKS(NWEQM,3)
      PARAMETER (TINY=1D-10)
      PARAMETER (ITRATE=-1,SOLVED=0,LIMIT=1,FLAT=2,ERROR=3)
      EXTERNAL SUB
      INTRINSIC ABS
C
C      TOL=5.D-7
C      TOL=1.D-4
C      MAXITS=100
C
      SUM=1.D0
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
         DDTEMP=DD
         CALL NEWTON_SUB(N,NWEQM,NEWX,DFDX,SUB,DDTEMP,WORKS,IERR)
         IF(IERR.NE.0) THEN
            STATUS=ERROR
            GOTO 11
         ENDIF
C
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
            CALL SUB(OLDX,FX,N,IERR)
            IF(IERR.NE.0) THEN
               STATUS=ERROR
               GOTO 11
            ENDIF
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
            SUM=0.D0
            DO I=1,M
               SUM=SUM+(DX(I)/MAX(ABS(OLDX(I)),1.D0))**2
            END DO
            SUM=SQRT(SUM/DFLOAT(M))
            WRITE(6,'(A,1PE12.4)') '** NEWTON4: SUM=',SUM
            IF (SUM.LE.TOL) STATUS=SOLVED
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
      SUBROUTINE NEWTON_SUB(N,NM,X,DFDX,SUB,DD,WORKS,IERR)
C
      IMPLICIT NONE
      INTEGER I,J,N,NM,IERR
      REAL*8 X(N),DFDX(NM,N)
      REAL*8 WORKS(NM,3)
      REAL*8 DD
      EXTERNAL SUB
C
C     WORKS(I,1) = XNEW(I)
C     WORKS(I,2) = RESULT0(I)
C     WORKS(I,3) = RESULT(I)
C
C      DD=1.D-2
C
      CALL SUB(X,WORKS(1,2),N,IERR)
      IF(IERR.NE.0) RETURN
      DO J=1,N
         DO I=1,N
            WORKS(I,1)=X(I)
         ENDDO
         WORKS(J,1)=X(J)+DD
         CALL SUB(WORKS(1,1),WORKS(1,3),N,IERR)
         IF(IERR.NE.0) RETURN
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
      SUBROUTINE TESTSUB(X,F,N,IERR)
C
      IMPLICIT NONE
      INTEGER N,IERR
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
      IERR=0
      RETURN
      END
