C     $Id$
C
C     ***********************************************************
C
C           BAND MATRIX SOLVER BY GAUSSIAN ELIMINATION
C
C     ***********************************************************
C
      SUBROUTINE BANDRD( A , X , N , L , LA , IERR )
C
      REAL*8  A(LA,N),X(N),TEMP
      REAL*8  EPS,ABS1,ABS2
      DATA EPS/ 1.D-70 /
C
      IF(MOD(L,2).EQ.0) GOTO 9000
      LH  = (L+1)/2
      LHM = LH-1
      NM  = N -1
      DO 30 K=1,LHM
         LHMK = LH-K
         NPMK = N+1-K
         DO 30 I=1,LHMK
            LPMI = L+1-I
            DO 40 J=2,L
               A(J-1,K) = A(J,K)
   40       CONTINUE
            A(L,K)       = 0.D0
            A(LPMI,NPMK) = 0.D0
   30 CONTINUE
      DO 50 I=1,NM
         IPIVOT = I
         IP     = I+1
         ABS2   = DABS( A(1,IPIVOT) )
         DO 60 K=IP,LH
            ABS1 = DABS( A(1,K) )
            IF(ABS1.GT.ABS2) THEN
                IPIVOT = K
                ABS2 = ABS1
            ENDIF
   60    CONTINUE
         IF(ABS2.LT.EPS) GOTO 9002
         IF(IPIVOT.NE.I) THEN
            TEMP      = X(I)
            X(I)      = X(IPIVOT)
            X(IPIVOT) = TEMP
            DO 90 J=1,L
               TEMP        = A(J,I)
               A(J,I)      = A(J,IPIVOT)
               A(J,IPIVOT) = TEMP
   90       CONTINUE
         ENDIF
         TEMP   = 1.D0/A(1,I)
         X(I)=X(I)*TEMP
         DO 120 J=2,L
            A(J,I) = A(J,I)*TEMP
  120    CONTINUE
         DO 130 K=IP,LH
            TEMP = A(1,K)
            X(K) = X(K)-X(I)*TEMP
            DO 140 J=2,L
               A(J-1,K) = A(J,K)-A(J,I)*TEMP
  140       CONTINUE
            A(L,K) = 0.D0
  130    CONTINUE
         IF(LH.LT.N)  LH = LH+1
   50 CONTINUE
      IF(DABS(A(1,N)).LT.EPS) GOTO 9002
      X(N) = X(N)/A(1,N)
      JJ = 2
      DO 160 I=1,NM
         K = N-I
         TEMP = 0.D0
         DO 170 J=2,JJ
            TEMP = A(J,K)*X(K-1+J)+TEMP
  170    CONTINUE
         X(K) = X(K)-TEMP
         IF(JJ.LT.L)  JJ = JJ+1
  160 CONTINUE
      IERR = 0
      RETURN
 9000 IERR = 10000
      RETURN
 9002 IERR = 30000
      RETURN
      END 
