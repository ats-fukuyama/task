C     $Id$
C
C     ****** REAL MATRIX INVERSION ******
C
      SUBROUTINE INVMRD(A,N,NA,ILL)
C
      REAL*8 A(NA,NA),W
      REAL*8 P,EPS
      DIMENSION  NM(64)
      DATA EPS/1.D-14/
C
      IF(N.GT.64) GOTO 999
C
      DO 10 NN=1,N
         NM(NN)=NN
   10 CONTINUE
C
      DO 100 NN=1,N
         P=0.D0
         IP=NN
         DO 20 I=NN,N
            IF(ABS(A(I,1)).GT.P) THEN
               P=ABS(A(I,1))
               IP=I
            ENDIF
   20    CONTINUE
         IF(P.LE.EPS) GO TO 900
         IF(IP.NE.NN) THEN
            NW=NM(IP)
            NM(IP)=NM(NN)
            NM(NN)=NW
            DO 30 J=1,N
               W=A(IP,J)
               A(IP,J)=A(NN,J)
               A(NN,J)=W
   30       CONTINUE
         ENDIF
         W=1.D0/A(NN,1)
         DO 40 J=2,N
            A(NN,J-1)=A(NN,J)*W
   40    CONTINUE
         A(NN,N)=W
         DO 50 I=1,N
            IF(I.NE.NN) THEN
               W=A(I,1)
               DO 60 J=2,N
                  A(I,J-1)=A(I,J)-W*A(NN,J-1)
   60          CONTINUE
               A(I,N)=-W*A(NN,N)
            ENDIF
   50    CONTINUE
  100 CONTINUE
C
      DO 150 NN=1,N
         DO 110 J=NN,N
            IF(NM(J).EQ.NN) GO TO 120
  110    CONTINUE
         GO TO 900
C
  120    IF(J.NE.NN) THEN
            NM(J)=NM(NN)
            DO 130 I=1,N
               W=A(I,J)
               A(I,J)=A(I,NN)
               A(I,NN)=W
  130       CONTINUE
         ENDIF
  150 CONTINUE
      ILL=0
      RETURN
C
  900 ILL=900
      RETURN
C
  999 WRITE(6,*) 'XX INVMRD: N.GT.64: N=',N
      ILL=999
      RETURN
      END
C
C     ***** COMPLEX MATRIX INVERSION *****
C
      SUBROUTINE INVMCD(A,N,NA,NM,ILL)
C
      COMPLEX*16 A(NA,NA),W
      REAL*8     P,EPS
      DIMENSION  NM(N)
      DATA EPS/1.D-14/
C
      DO 10 NN=1,N
         NM(NN)=NN
   10 CONTINUE
C
      DO 100 NN=1,N
         P=0.D0
         IP=NN
         DO 20 I=NN,N
            IF(ABS(A(I,1)).GT.P) THEN
               P=ABS(A(I,1))
               IP=I
            ENDIF
   20    CONTINUE
         IF(P.LE.EPS) GO TO 900
         IF(IP.NE.NN) THEN
            NW=NM(IP)
            NM(IP)=NM(NN)
            NM(NN)=NW
            DO 30 J=1,N
               W=A(IP,J)
               A(IP,J)=A(NN,J)
               A(NN,J)=W
   30       CONTINUE
         ENDIF
         W=(1.D0,0.D0)/A(NN,1)
         DO 40 J=2,N
            A(NN,J-1)=A(NN,J)*W
   40    CONTINUE
         A(NN,N)=W
         DO 50 I=1,N
            IF(I.NE.NN) THEN
               W=A(I,1)
               DO 60 J=2,N
                  A(I,J-1)=A(I,J)-W*A(NN,J-1)
   60          CONTINUE
               A(I,N)=-W*A(NN,N)
            ENDIF
   50    CONTINUE
  100 CONTINUE
C
      DO 150 NN=1,N
         DO 110 J=NN,N
            IF(NM(J).EQ.NN) GO TO 120
  110    CONTINUE
         GO TO 900
C
  120    IF(J.NE.NN) THEN
            NM(J)=NM(NN)
            DO 130 I=1,N
               W=A(I,J)
               A(I,J)=A(I,NN)
               A(I,NN)=W
  130       CONTINUE
         ENDIF
  150 CONTINUE
      ILL=0
      RETURN
C
  900 ILL=900
      RETURN
C
      END
