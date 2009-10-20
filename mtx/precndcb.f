C     $Id$
C
C     ######## PRECONDITIONING ########
C                  (JACOBI)
C
      SUBROUTINE CALJACB(A,LA,N,M,X,Y,NBSIZ,IERR)
C
      COMPLEX * 16 A(LA,N),X(N),Y(N)
C
      IF(N.LE.0.OR.NBSIZ.LE.0) THEN
         IERR=1
      ENDIF
C
      DO I=1,N
         X(I)=Y(I)/A(M+1,I)
      ENDDO
C
      RETURN
      END
C
C     ######## PRECONDITIONING ROUTINE ########
C       (INCOMPLETE LDU DICOMPOSITION VERSION)
C     
      SUBROUTINE CALCUDB(A,D,TEMP1,NM,LA,LD,N,M,NBSIZ,IERR)
C
      COMMON /WMBLPR/ ICF,IDF,IEF,ISCNF,IECNF
C
      COMPLEX * 16 A(LA,N),D(LD,N)
      COMPLEX * 16 TEMP1(LD,NBSIZ),NM(NBSIZ)
C
      IF(N.LE.0.OR.NBSIZ.LE.0) THEN
         IERR=1
      ENDIF
C
      DO IBN=1,N/NBSIZ
         IBLC=NBSIZ*(IBN-1)+1
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               TEMP1(I,J)=(0.D0,0.D0)
               DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
                  K=IC+NBSIZ-IBLC+1
                  TEMP1(I,J)=TEMP1(I,J)
     &                      +A(ICF +(K-1)-(J-1),IBLC+(J-1))
     &                      *D(I,IC)
               ENDDO
            ENDDO
         ENDDO
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               D(I,IBLC+(J-1))=(0.D0,0.D0)
               DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
                  K=IC+NBSIZ-IBLC+1
                  D(I,IBLC+(J-1))=D(I,IBLC+(J-1))
     &                           +TEMP1(K,J)
     &                           *A(IEF+(I-1)-(K-1),IC)
               ENDDO
               D(I,IBLC+(J-1))=A(IDF+(I-1)-(J-1),IBLC+(J-1))
     &                        -D(I,IBLC+(J-1))
            ENDDO
         ENDDO
C
         CALL INVMCD(D(1,IBLC),NBSIZ,LD,IERR)
C
      ENDDO
C
      RETURN
      END
C
C
C     ######## CALCULATE (LDU)^{1}*A ########
C
      SUBROUTINE CALLDUB(A,D,TEMP2,X,LA,LD,N,M,NBSIZ,IERR)
C
      COMMON /WMBLPR/ ICF,IDF,IEF,ISCNF,IECNF
C
      COMPLEX * 16 A(LA,N),D(LD,N),X(N),TEMP2(NBSIZ)
C
      IF(N.LE.0.OR.NBSIZ.LE.0) THEN
         IERR=1
      ENDIF
C
      DO IBN=1,N/NBSIZ
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=X(IBLC+(I-1))
            DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)-A(ICF+(K-1)-(I-1),IBLC+(I-1))
     &                          *X(IC)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            X(IBLC+(I-1))=(0.D0,0.D0)
            DO K=1,NBSIZ
               X(IBLC+(I-1))=X(IBLC+(I-1))
     &                       +D(K,IBLC+(I-1))*TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO IBN=N/NBSIZ,1,-1
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=(0.D0,0.D0)
            DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,N)
               K=IC-IBLC-NBSIZ+1
               TEMP2(I)=TEMP2(I)+A (IEF+(K-1)-(I-1),IBLC+(I-1))
     &                          *X(IC)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO K=1,NBSIZ
               X(IBLC+(I-1))=X(IBLC+(I-1))
     &                      -D(K,IBLC+(I-1))
     &                      *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
