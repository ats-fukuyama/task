C     $Id$
C
C
C     ====== FFT USING LIST VECTOR ======
C                                         CODING BY K.HAMAMATSU
C
      SUBROUTINE FFT2L ( A , B , T , LIST , N2 , IND , KEY , LP )
C
      IMPLICIT COMPLEX * 16 ( A-C , E-H , O-Z )
      REAL * 8     DIV
      DIMENSION    A ( N2*2 ) , B( N2*2 ) , T( N2*LP , 2 ) ,
     &             LIST( N2*LP )
      SAVE DIV
C
C     IF( IND .EQ. 1 ) THEN SET TALBES AND LIST
      IF( IND .EQ. 1 ) THEN
           CALL SETTBL ( N2 , T , B , LP )
           CALL SETLST ( LIST , N2 , LP )
           DIV = 1.D0 / (N2*2)
           IND = 0
      ENDIF
      K =  1
      L = N2
      DO 10 I = 1 , LP
           M = (I-1)*N2 + 1
           IF( K .EQ. 1 ) THEN
                CALL FFTSUB( A , B , T(M,KEY) , LIST(M) , L , N2 )
           ELSE
                CALL FFTSUB( B , A , T(M,KEY) , LIST(M) , L , N2 )
           ENDIF
           K = K * (-1)
           L = L / 2
 10   CONTINUE
C     IF( KEY .EQ. 2 ) THEN INVERSE TRANSFORMATION
C     ELSE    NORMAL TRANSFORMATION
      IF( KEY .EQ. 1 ) THEN
           IF( K .EQ. 1 ) THEN
                DO 20 I = 1 , N2*2
 20                 B( I ) = A( I ) * DIV
            ELSE
                DO 30 I = 1 , N2*2
 30                 B( I ) = B( I ) * DIV
            ENDIF
      ELSEIF( K .EQ. 1 ) THEN
            DO 40 I = 1 , N2*2
 40             B( I ) = A ( I )
      ENDIF
      RETURN
      END
C
C     ====== TRANSFORMATION FOR FFT ======
C                                          CODING BY K.HAMAMATSU
C
      SUBROUTINE FFTSUB ( A , B , T , LIST , L , N2 )
C
      IMPLICIT COMPLEX * 16 ( A-H , O-Z )
      DIMENSION A( N2*2 ) , B( N2 , 2 ) , T( N2 ) , LIST( N2 )
      DO 10 I = 1 , N2
          B( I , 1 ) = A( LIST( I ) ) + A( LIST( I ) + L ) * T( I )
 10   CONTINUE
      DO 20 I = 1 , N2
          B( I , 2 ) = A( LIST( I ) ) - A( LIST( I ) + L ) * T( I )
 20   CONTINUE
      RETURN
      END
C
C     ====== TABLE SETTING FOR FFT ======
C                                          CODING BY K.HAMAMATSU
C
      SUBROUTINE SETTBL ( N2 , T , B , LP )
C
      COMPLEX * 16     T , B
      REAL    *  8     TR , TI , PAI
      DIMENSION       T( N2 , LP , 2 ) , B( N2 , 2 )
      DATA PAI / 3.1415926535898D0 /
C
      DO 10 I = 1 , N2
          TR = DCOS( 2.D0 * PAI * (I-1) / (N2*2) )
          TI = DSIN( 2.D0 * PAI * (I-1) / (N2*2) )
          B( I , 1 ) = DCMPLX( TR , -TI )
          B( I , 2 ) = DCMPLX( TR ,  TI )
 10   CONTINUE
      K  =  1
      NN = N2
      DO 40 L = 1 , LP
          DO 30 J = 0 , K-1
              DO 20 I = 1 , NN
                  T( I+J*NN , L , 1 ) = B( 1+NN*J , 1 )
                  T( I+J*NN , L , 2 ) = B( 1+NN*J , 2 )
 20           CONTINUE
 30       CONTINUE
          K  = K * 2
          NN = NN / 2
 40   CONTINUE
      RETURN
      END
C
C     ====== LIST SETTING FOR FFT ======
C                                          CODING BY K.HAMAMATSU
C
      SUBROUTINE SETLST ( LIST , N2 , LP )
C
      DIMENSION  LIST ( N2 , LP )
      N1 = N2
      NN =  1
      DO 30 K = 1 , LP
          M = 0
          DO 20 J = 1 , NN
              DO 10 I = 1 , N1
                  M = M + 1
                  LIST ( M , K ) = I + (J-1) * 2 * N1
 10           CONTINUE
 20       CONTINUE
          N1 = N1 / 2
          NN = NN * 2
 30   CONTINUE
      RETURN
      END
