C     $Id$
C
C     ***************************************************************
C
C        BAND MATRIX SOLVER BY GAUSSIAN ELIMINATION
C
C     ***************************************************************
C
      SUBROUTINE BANDRD(A , X , N , L , LA , IERR)
C
      REAL*8 A(LA,N), X(N), TEMP
      REAL*8 EPS, ABS1, ABS2
      DATA EPS/ 1.D-70 /
C
CCCCC
CC      WRITE(6,*) 'A(L,N)'
CC      WRITE(6,*) ((A(I,J),I=1,L),J=1,N)
CC      WRITE(6,*) 'X(N)'
CC      WRITE(6,*) (X(I),I=1,N)
CCCCC
      IF (MOD(L,2) .EQ. 0) GOTO 110
      LH  = (L+1) / 2
      LHM = LH - 1
      NM  = N - 1
      DO 20 K = 1, LHM
         LHMK = LH - K
         NPMK = N+1 - K
         DO 20 I = 1, LHMK
            LPMI = L+1 - I
            DO 10 J = 2, L
               A(J-1,K) = A(J,K)
   10       CONTINUE
            A(L,K)       = 0.D0
            A(LPMI,NPMK) = 0.D0
   20 CONTINUE
      DO 80 I = 1, NM
         IPIVOT = I
         IP     = I+1
         ABS2   = DABS(A(1,IPIVOT))
         DO 30 K = IP, LH
            ABS1 = DABS(A(1,K))
            IF (ABS1 .GT. ABS2) THEN
               IPIVOT = K
               ABS2 = ABS1
            ENDIF
   30    CONTINUE
         IF (ABS2 .LT. EPS) GOTO 120
         IF (IPIVOT .NE. I) THEN
C            write(6,*) 'PIVOT : ',I,IPIVOT
            TEMP      = X(I)
            X(I)      = X(IPIVOT)
            X(IPIVOT) = TEMP
            DO 40 J = 1, L
               TEMP        = A(J,I)
               A(J,I)      = A(J,IPIVOT)
               A(J,IPIVOT) = TEMP
   40       CONTINUE
         ENDIF
         TEMP = 1.D0 / A(1,I)
         X(I) = X(I)*TEMP
         DO 50 J = 2, L
            A(J,I) = A(J,I)*TEMP
   50    CONTINUE
         DO 70 K = IP, LH
            TEMP = A(1,K)
            X(K) = X(K) - X(I)*TEMP
            DO 60 J = 2, L
               A(J-1,K) = A(J,K) - A(J,I)*TEMP
   60       CONTINUE
            A(L,K) = 0.D0
   70    CONTINUE
         IF (LH .LT. N)  LH = LH+1
   80 CONTINUE
      IF (DABS(A(1,N)) .LT. EPS) GOTO 120
      X(N) = X(N) / A(1,N)
      JJ = 2
      DO 100 I = 1, NM
         K = N-I
         TEMP = 0.D0
         DO 90 J = 2, JJ
            TEMP = A(J,K)*X(K-1+J) + TEMP
   90    CONTINUE
         X(K) = X(K) - TEMP
         IF (JJ .LT. L) JJ = JJ+1
  100 CONTINUE
      IERR = 0
CCCCC
CC      WRITE(6,*) 'X(N)'
CC      WRITE(6,*) (X(I),I=1,N)
CCCCC
      RETURN
  110 IERR = 10000
      RETURN
  120 IERR = 30000
      WRITE(6,*) '## BAND RD ERROR AT I = ', I
      RETURN
      END
