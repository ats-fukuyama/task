!     *****************************************************
!     **     numerical solution of band matrix (SUB)     **
!     **          Gaussian Elimination Method            **
!     *****************************************************

      SUBROUTINE BANDRD( A , X , N , L , LA , IERR )

!          INPUT : A(LA,N) : COEFFICIENT MATRIX
!                  X(N)    : RIGHT-HAND-SIDE VECTOR
!                  N       : MATRIX SIZE
!                  L       : BAND WIDTH (L.LE.LA)
!                  LA      : SIZE OF ARRAY A
!          OUTPUT: X(N)    : SOLUTION VECTOR
!                  IERR    : ERROR CODE : 0 : NORMAL END
!                                         10000 : L IS EVEN
!                                         30000 : SINGULAR MATRIX
!          NOTICE: ARRAY A AND X WILL BE DESTROYED.

      IMPLICIT NONE
      REAL(8), DIMENSION(LA, N),INTENT(INOUT) :: A
      REAL(8), DIMENSION(N),    INTENT(INOUT) :: X
      INTEGER(4),               INTENT(IN)    :: N, L, LA
      INTEGER(4),               INTENT(OUT)   :: IERR
      REAL(8)            :: ABS1, ABS2, TEMP
      REAL(8), PARAMETER :: EPS = 1.D-70
      INTEGER(4)         :: I, J, K, LH, LHM, NM, LHMK, NPMK, LPMI, IPIVOT, IP, JJ

      IF( MOD(L,2) .EQ. 0 ) GO TO 9000
      LH  = (L+1)/2
      LHM = LH-1
      NM  = N -1

      DO K = 1 , LHM
         LHMK = LH-K
         NPMK = N+1-K
         DO I = 1 , LHMK
            LPMI = L+1-I
            A( 1:L-1 , K ) = A( 2:L , K )
            A( L    , K    ) = 0.D0
            A( LPMI , NPMK ) = 0.D0
         ENDDO
      ENDDO

      DO I = 1 , NM
         IPIVOT = I
         IP     = I+1
         ABS2   = ABS( A(1,IPIVOT) )
         DO K = IP , LH
            ABS1 = ABS( A(1,K) )
            IF( ABS1 .GT. ABS2 ) THEN
                IPIVOT = K
                ABS2 = ABS1
            ENDIF
         ENDDO

         IF( ABS2 .LT. EPS ) GO TO 9002
         IF( IPIVOT .NE. I ) THEN
            TEMP        = X( I      )
            X( I      ) = X( IPIVOT )
            X( IPIVOT ) = TEMP
            DO J = 1 , L
               TEMP            = A( J , I      )
               A( J , I      ) = A( J , IPIVOT )
               A( J , IPIVOT ) = TEMP
            ENDDO
         END IF

         TEMP   = 1.D0   / A( 1 , I )
         X( I ) = X( I ) * TEMP

         A( 2:L , I ) = A( 2:L , I ) * TEMP

         DO K = IP , LH
            TEMP   = A( 1 , K )
            X( K ) = X( K ) - X( I ) * TEMP
            A( 1:L-1 , K ) = A( 2:L , K ) - A( 2:L , I ) * TEMP

            A( L , K ) = 0.D0
         ENDDO
         IF( LH .LT. N ) LH = LH + 1
      ENDDO

      IF( ABS(A(1,N)) .LT. EPS ) GO TO 9002
      X( N ) = X( N ) / A( 1 , N )
      JJ = 2
      DO I = 1 , NM
         K = N-I
         TEMP = 0.D0
         DO J = 2 , JJ
            TEMP = A( J , K ) * X( K-1+J ) + TEMP
         ENDDO
         X( K ) = X( K ) - TEMP
         IF( JJ .LT. L ) JJ = JJ + 1
      ENDDO

      IERR = 0
      RETURN

 9000 IERR = 10000
      RETURN
 9002 IERR = 30000+I
      RETURN
      END SUBROUTINE BANDRD

!----- complex version ----

      SUBROUTINE BANDCD( A , X , N , L , LA , IERR )

!          INPUT : A(LA,N) : COEFFICIENT MATRIX
!                  X(N)    : RIGHT-HAND-SIDE VECTOR
!                  N       : MATRIX SIZE
!                  L       : BAND WIDTH (L.LE.LA)
!                  LA      : SIZE OF ARRAY A
!          OUTPUT: X(N)    : SOLUTION VECTOR
!                  IERR    : ERROR CODE : 0 : NORMAL END
!                                         10000 : L IS EVEN
!                                         30000 : SINGULAR MATRIX
!          NOTICE: ARRAY A AND X WILL BE DESTROYED.

      IMPLICIT NONE
      COMPLEX(8), DIMENSION(LA, N),INTENT(INOUT) :: A
      COMPLEX(8), DIMENSION(N),    INTENT(INOUT) :: X
      INTEGER(4),               INTENT(IN)    :: N, L, LA
      INTEGER(4),               INTENT(OUT)   :: IERR
      REAL(8) :: ABS1, ABS2
      COMPLEX(8) :: TEMP
      REAL(8), PARAMETER :: EPS = 1.D-70
      INTEGER(4) :: I, J, K, LH, LHM, NM, LHMK, NPMK, LPMI, IPIVOT, IP, JJ

      IF( MOD(L,2) .EQ. 0 ) GO TO 9000
      LH  = (L+1)/2
      LHM = LH-1
      NM  = N -1

      DO K = 1 , LHM
         LHMK = LH-K
         NPMK = N+1-K
         DO I = 1 , LHMK
            LPMI = L+1-I
            A( 1:L-1 , K ) = A( 2:L , K )
            A( L    , K    ) = 0.D0
            A( LPMI , NPMK ) = 0.D0
         ENDDO
      ENDDO

      DO I = 1 , NM
         IPIVOT = I
         IP     = I+1
         ABS2   = ABS( A(1,IPIVOT) )
         DO K = IP , LH
            ABS1 = ABS( A(1,K) )
            IF( ABS1 .GT. ABS2 ) THEN
                IPIVOT = K
                ABS2 = ABS1
            ENDIF
         ENDDO

!         IF( ABS2 .LT. EPS ) GO TO 9002
      IF( ABS2 .LT. EPS ) THEN
         write(6,'(A,I8,1P3E12.4)') 'A(1,I),EPS=',I,A(1,I),EPS
         GOTO 9002
      ENDIF

         IF( IPIVOT .NE. I ) THEN
            TEMP        = X( I      )
            X( I      ) = X( IPIVOT )
            X( IPIVOT ) = TEMP
            DO J = 1 , L
               TEMP            = A( J , I      )
               A( J , I      ) = A( J , IPIVOT )
               A( J , IPIVOT ) = TEMP
            ENDDO
         END IF

         TEMP   = 1.D0   / A( 1 , I )
         X( I ) = X( I ) * TEMP

         A( 2:L , I ) = A( 2:L , I ) * TEMP

         DO K = IP , LH
            TEMP   = A( 1 , K )
            X( K ) = X( K ) - X( I ) * TEMP
            A( 1:L-1 , K ) = A( 2:L , K ) - A( 2:L , I ) * TEMP

            A( L , K ) = 0.D0
         ENDDO
         IF( LH .LT. N ) LH = LH + 1
      ENDDO

!      IF( ABS(A(1,N)) .LT. EPS ) GO TO 9002
!      IF( CDABS(A(1,N)) .LT. EPS ) GOTO 9002
      IF( CDABS(A(1,N)) .LT. EPS ) THEN
         write(6,'(A,1P3E12.4)') 'A(1,N),EPS=',A(1,N),EPS
         GOTO 9002
      ENDIF

      X( N ) = X( N ) / A( 1 , N )
      JJ = 2
      DO I = 1 , NM
         K = N-I
         TEMP = 0.D0
         DO J = 2 , JJ
            TEMP = A( J , K ) * X( K-1+J ) + TEMP
         ENDDO
         X( K ) = X( K ) - TEMP
         IF( JJ .LT. L ) JJ = JJ + 1
      ENDDO

      IERR = 0
      RETURN

 9000 IERR = 10000
      RETURN
 9002 IERR = 30000+I
      RETURN
      END SUBROUTINE BANDCD
