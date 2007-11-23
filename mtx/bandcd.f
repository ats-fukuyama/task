C     $Id$
C
C     ****** SOLUTION OF BAND MATRIX (GAUSSIAN ELIMINATION) ******
C
      SUBROUTINE BANDCDNB( A , X , N , L , LA , IERR )
C
C     COMPLEX * 16    A( LA , N ) , X( N ) , ATMP( LMAX ) , TEMP
      COMPLEX * 16    A( LA , N ) , X( N ) , TEMP
      REAL    *  8    EPS , ABS1 , ABS2
      DATA EPS/ 1.D-70 /
C
      NRP=0
      DO MS=1,N
         X(MS)=(0.D0,0.D0)
         DO MB=1,L
            A(MB,MS)=(0.D0,0.D0)
         ENDDO
         CALL WMSETM(A(1,MS),X(MS),MS,LA,NRP)
      ENDDO
C
      IF( MOD(L,2) .EQ. 0 ) GOTO 9000
      LH  = (L+1)/2
      LHM = LH-1
      NM  = N -1
C
      DO K = 1 , LHM
         LHMK = LH-K
         NPMK = N+1-K
         DO I = 1 , LHMK
            LPMI = L+1-I
            DO J = 2 , L
               A( J-1 , K ) = A( J , K )
            ENDDO
            A( L    , K    ) = ( 0.D0 , 0.D0 )
            A( LPMI , NPMK ) = ( 0.D0 , 0.D0 )
         ENDDO
      ENDDO
C
      DO I = 1 , NM
         IPIVOT = I
         IP     = I+1
         ABS2   = CDABS( A(1,IPIVOT) )
         DO K = IP , LH
            ABS1 = CDABS( A(1,K) )
            IF( ABS1 .GT. ABS2 ) THEN
                IPIVOT = K
                ABS2 = ABS1
            ENDIF
         ENDDO
C
         IF( CDABS(A(1,IPIVOT)) .LT. EPS ) THEN
            write(6,'(A,1P3E12.4)') 'A(1,IPIVOT),EPS=',A(1,IPIVOT),EPS
            GOTO 9001
         END IF
C
         IF( IPIVOT .NE. I ) THEN
            TEMP        = X( I      )
            X( I      ) = X( IPIVOT )
            X( IPIVOT ) = TEMP
            DO J = 1 , L
C              ATMP( J )          = A   ( J , I      )
               TEMP               = A   ( J , I      )
               A   ( J , I      ) = A   ( J , IPIVOT )
C              A   ( J , IPIVOT ) = ATMP( J )
               A   ( J , IPIVOT ) = TEMP
            ENDDO
         END IF
C
         TEMP   = 1.D0   / A( 1 , I )
         X( I ) = X( I ) * TEMP
C
         DO J = 2 , L
            A( J , I ) = A( J , I ) * TEMP
         ENDDO
C
         DO K = IP , LH
            TEMP   = A( 1 , K )
            X( K ) = X( K ) - X( I ) * TEMP
            DO J = 2 , L
               A( J-1 , K ) = A( J , K ) - A( J , I ) * TEMP
            ENDDO
C
            A( L , K ) = ( 0.D0 , 0.D0 )
         ENDDO
         IF( LH .LT. N ) LH = LH + 1
      ENDDO
C
      IF( CDABS(A(1,N)) .LT. EPS ) THEN
         write(6,'(A,1P3E12.4)') 'A(1,N),EPS=',A(1,N),EPS
         GOTO 9002
      END IF
C
      X( N ) = X( N ) / A( 1 , N )
      JJ = 2
      DO I = 1 , NM
         K = N-I
         TEMP = ( 0.D0 , 0.D0 )
         DO J = 2 , JJ
            TEMP = A( J , K ) * X( K-1+J ) + TEMP
         ENDDO
         X( K ) = X( K ) - TEMP
         IF( JJ .LT. L ) JJ = JJ + 1
      ENDDO
C
      IERR = 0
      RETURN
C
 9000 IERR = 10000
      RETURN
 9001 IERR = 20000+I
      RETURN
 9002 IERR = 30000+I
      RETURN
C
      END
