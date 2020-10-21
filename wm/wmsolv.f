C     $Id$
C
C     ****** SOLVE SIMULTANEOUS EQUATIONS ******
C
      SUBROUTINE WMSOLV
C
      INCLUDE 'wmcomm.inc'
C
      EXTERNAL WMSETM
C
      PARAMETER (NBSIZM=3*MDM*NDM)
C
      NBSIZ=3*MDSIZ*NDSIZ
      IF(MODEWG.EQ.0) THEN
         MSIZ=NRMAX*NBSIZ
      ELSE
         MSIZ=NRMAX*NBSIZ+MWGMAX*NAMAX
      ENDIF
      MBND=2*NBSIZ

      IF(MODEEG.EQ.0) THEN
         CALL WMSOLV_MTXP(CFVG,MSIZ,2*MBND-1,WMSETM,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XXX WMSOLV_MTXP: ERROR: IERR=',IERR
      ELSE
         CALL WMSOLV_BAND(CFVG,MSIZ,2*MBND-1,WMSETM,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XXX WMSOLV_BAND: ERROR: IERR=',IERR
      END IF

      RETURN
      END

      SUBROUTINE WMSOLV_MTXP(XX , N , L , WMSETM, IERR )

      USE libmtx
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(N),INTENT(OUT):: XX
      INTEGER,INTENT(IN):: N,L
      INTEGER,INTENT(OUT)::IERR 
      EXTERNAL WMSETM
      COMPLEX(8),DIMENSION(:),ALLOCATABLE:: A
      COMPLEX(8):: X
      INTEGER:: istart,iend,i,j,NRP
      INTEGER:: itype,its
      REAL(8):: tolerance

      ALLOCATE(A(L))        ! local coefficient matrix

      CALL mtxc_setup(N,istart,iend,jwidth=L)

C   ***** CALCULATE MATRIX COEFFICIENTS *****

      NRP=0
      DO i=istart,iend
         X=(0.D0,0.D0)
         A(1:L)=(0.D0,0.D0)
         CALL WMSETM(A,X,i,L,NRP)
         WRITE(21,'(A,2I6)') 'wmsolv:',NRP,i
         WRITE(21,'(6ES12.4)') (A(j),j=1,L)
         WRITE(21,'(2ES12.4)') X
         DO j=MAX(i-(L+1)/2+1,1),MIN(N,i+(L+1)/2-1)
            IF(ABS(A(j-i+(L+1)/2)).GT.0.D0) THEN
               CALL mtxc_set_matrix(i,j,A(j-i+(L+1)/2))
            END IF
         END DO
         CALL mtxc_set_source(i,X)
      END DO

      itype=0
      tolerance=1.D-12
      CALL mtxc_solve(itype,tolerance,its)

      CALL mtxc_gather_vector(XX)

      CALL mtxc_cleanup
C
      DEALLOCATE(A)
      IERR=0

      RETURN
      END
C
C     ****** SOLUTION OF BAND MATRIX (GAUSSIAN ELIMINATION) ******
C
      SUBROUTINE WMSOLV_BAND(XX , N , L , WMSETM, IERR )
C
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(:,:),ALLOCATABLE:: A
      COMPLEX(8),DIMENSION(:),ALLOCATABLE:: X
      COMPLEX(8),DIMENSION(N):: XX
      COMPLEX(8):: TEMP
      REAL(8):: EPS , ABS1 , ABS2
      INTEGER:: N,L,NRP,NR,MS,MB,LH,LHM,NM,K,LHMK,NPMK,I,LPMI,J
      INTEGER:: IPIVOT,IP,JJ,IERR
      EXTERNAL WMSETM
      DATA EPS/ 1.D-70 /
C
      ALLOCATE(A(L,N))
      ALLOCATE(X(N))
C
      NRP=0
      DO MS=1,N
         X(MS)=(0.D0,0.D0)
         DO MB=1,L
            A(MB,MS)=(0.D0,0.D0)
         ENDDO
         CALL WMSETM(A(1,MS),X(MS),MS,L,NRP)
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
      XX(1:N)=X(1:N)
      DEALLOCATE(A,X)
      IERR = 0
      RETURN
C
 9000 IERR = 10000
      DEALLOCATE(A,X)
      RETURN
 9001 IERR = 20000+I
      DEALLOCATE(A,X)
      RETURN
 9002 IERR = 30000+I
      DEALLOCATE(A,X)
      RETURN
C
      END
