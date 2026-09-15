! $Id$

MODULE wihot

  PRIVATE
  PUBLIC wi_hot

CONTAINS

  SUBROUTINE wi_hot(iprint,ratea,ierr)

    USE wicomm
    USE libinv
    USE libbnd
    IMPLICIT NONE
    INTEGER(ikind),INTENT(IN):: iprint
    REAL(rkind),INTENT(OUT):: ratea
    INTEGER(ikind),INTENT(OUT):: ierr

    mlmax=nxmax*2+3
    mwmax=4*nwmax+3

    CALL SUBFW    ! calculate elements of kernel function
    CALL SUBCK2   ! calculate coefficient matrix
    CALL SUBINI   ! calculate right-hand-side vector
    IF(NWMAX.EQ.NXMAX) THEN
       CALL INVMCD(CK,mlmax,MLEN,IERR)   ! full matrix solver
       IF(IERR.NE.0) GOTO 9900
       CALL SUBFY                  ! calculate field vector
    ELSE
!       DO ML=1,4
!        WRITE(6,'(I5,1P6E12.4)') ML,(CK(MW,ML),MW=(MWMAX+1)/2-1,(MWMAX+1)/2+1)
!       END DO
!       DO ML=400,410
!        WRITE(6,'(I5,1P6E12.4)') ML,(CK(MW,ML),MW=(MWMAX+1)/2-1,(MWMAX+1)/2+1)
!       END DO
!       DO ML=MLMAX-3,MLMAX
!        WRITE(6,'(I5,1P6E12.4)') ML,(CK(MW,ML),MW=(MWMAX+1)/2-1,(MWMAX+1)/2+1)
!       END DO
       CALL BANDCD(CK,CSO,mlmax,mwmax,MWID,IERR)   ! band matrix solver
          IF(IERR.NE.0) GOTO 9900
       CALL SUBFYW                               ! calculate field vector
    ENDIF
    CALL SUBPOW    ! calculate sbsorbed power
!!!       RATEA=1.D0-ABS(CFY(NXMAX*2+3))**2
!       WRITE(6,'(A,1PE12.4)') 'ABS(CFY(NXMAX*2+2))**2=',ABS(CFY(NXMAX*2+2))**2
!       WRITE(6,'(A,1PE12.4)') 'ABS(CFY(NXMAX*2+3))**2=',ABS(CFY(NXMAX*2+3))**2
!       WRITE(6,'(A,1P2E12.4)') 'CFY(NXMAX*2+2)=',CFY(NXMAX*2+2)
!       WRITE(6,'(A,1P2E12.4)') 'CFY(NXMAX*2+3)=',CFY(NXMAX*2+3)
       RATEA=1.D0-ABS(CFY(NXMAX*2+3))**2
       IF(iprint > 0) WRITE(6,'(A,ES12.4)') '## Absorption rate: ',RATEA
9900  CONTINUE
      RETURN
    END SUBROUTINE wi_hot

!     *****  SET KERNEL FUNCTION  ***** 

    SUBROUTINE SUBFW

      USE wicomm
      USE wieul
      IMPLICIT NONE
      REAL(rkind):: dx,rky,x
      COMPLEX(rkind):: CS
      INTEGER(ikind):: J,L,NW

      DX=(XMAX-XMIN)/NXMAX
      RKY=ANY
      DO J=1,2
         DO NW=0,NWMAX
            X=NW*DX
            CALL wi_eul(X,ALFA,BETA,RKY,CS,J,5,L)
            CU(J, NW)=CS
            CU(J,-NW)=CS
         END DO
      END DO
      RETURN
    END SUBROUTINE SUBFW

!     *****  CALCULATION OF COEFFICIENT MATRIX  ***** 

    SUBROUTINE SUBCK2

      USE wicomm
      IMPLICIT NONE
      COMPLEX(rkind):: ciky,cbb
      REAL(rkind):: rky,rky2,dx,dx2,dky
      REAL(rkind):: ANB,beta0
      INTEGER(ikind):: NDUB,NBAND,NWDUB,NWDDUB,I,J,MM,ID,JD,NS,NE,NN
      INTEGER(ikind):: KK,KD,KS,IOB,IO,I2

      RKY=ANY
      RKY2=RKY**2
      DKY=ANY*ANY
      CIKY=CI*ANY
      ANB=DEXP(-ALFA*xgrid(nxmax))
      CBB=CI/DSQRT(1.D0-ANB-ANY*ANY)
      BETA0=BETA

      NDUB=2*NXMAX
      IF(NWMAX.EQ.NXMAX) THEN
         NBAND=0
         NWDUB=NDUB
         NWDDUB=NDUB
      ELSE
         NBAND=1
         NWDUB=2*NWMAX
         NWDDUB=4*NWMAX
      ENDIF

      DO I=1,NDUB+3
         DO J=1,NWDDUB+3
            CK(J,I)=(0.,0.)
         END DO
      END DO
      DO MM=0,NXMAX-1
         DX=xgrid(MM+1)-xgrid(MM)
         DX2=DX*DX
         DO I=MM,MM+1
            ID=2*I
            DO J=MM,MM+1
               JD=2*J
               IF(NWMAX.NE.NXMAX) JD=2*NWMAX+1+2*J-2*I
               CK(JD+1      ,ID+1)=CK(JD+1      ,ID+1) &
                                    +(DKY-1.D0)*DX*D0(I-MM,J-MM)
               CK(JD+2      ,ID+1)=CK(JD+2      ,ID+1) &
                                    +CIKY*D1(J-MM,I-MM)
               CK(JD+1-NBAND,ID+2)=CK(JD+1-NBAND,ID+2) &
                                    -CIKY*D1(I-MM,J-MM)
               CK(JD+2-NBAND,ID+2)=CK(JD+2-NBAND,ID+2) &
                                    +D2(I-MM,J-MM)/DX &
                                    -DX*D0(I-MM,J-MM)
            END DO
         END DO
      END DO
      CK(NWDUB+3      ,NDUB+2)=-CBB
      CK(NWDUB+2-NBAND,NDUB+3)=1.D0
      CK(NWDUB+3-NBAND,NDUB+3)=-1.D0

      DO MM=0,NXMAX-1
         NS=MM-NWMAX+1
         NE=MM+NWMAX-1
         IF(NS.LE.0) NS=0
         IF(NE.GE.NXMAX-1) NE=NXMAX-1
         IF(XMAX.GE.500.D0.AND.ALFA*XMAX.LT.10.D0) THEN
            IF(XMAX-XGRID(MM).LT.Bwidth) THEN
               BETA=BETA0*(XMAX-XGRID(MM))/Bwidth
            ELSE
               BETA=BETA0
            END IF
         ELSE
            BETA=BETA0
         END IF
         DO NN=NS,NE
            DO I=MM,MM+1
               ID=2*I
               DO J=NN,NN+1
                  JD=2*J
                  IF(NWMAX.NE.NXMAX) JD=2*NWMAX+1+2*J-2*I
                  DO KK=MM,MM+1
                     DO KD=NN,NN+1
                        CK(JD+1,ID+1)=CK(JD+1,ID+1) &
                                     -CWP(KD)*CWE(KK)*CWE(KD)*BETA &
                                     *(DX2*RKY2*CU(1,KK-KD) &
                                     *D0(I-MM,KK-MM)*D0(J-NN,KD-NN) &
                                     +(CU(1,KK-KD)-CI*CU(2,KK-KD)) &
                                     *D1(I-MM,KK-MM)*D1(J-NN,KD-NN))
                        CK(JD+2,ID+1)=CK(JD+2,ID+1) &
                                     -CWP(KD)*CWE(KK)*CWE(KD)*BETA &
                                     *(-DX*RKY*CU(2,KK-KD) &
                                     *D0(I-MM,KK-MM)*D1(J-NN,KD-NN))
                        CK(JD+1-NBAND,ID+2)=CK(JD+1-NBAND,ID+2) &
                                     -CWP(KD)*CWE(KK)*CWE(KD)*BETA &
                                     *(DX*RKY*CU(2,KK-KD) &
                                     *D1(I-MM,KK-MM)*D0(J-NN,KD-NN))
                        CK(JD+2-NBAND,ID+2)=CK(JD+2-NBAND,ID+2) &
                                     -CWP(KD)*CWE(KK)*CWE(KD)*BETA &
                                *(RKY2*DX2*(CU(1,KK-KD)-CI*CU(2,KK-KD)) &
                                     *D0(I-MM,KK-MM)*D0(J-NN,KD-NN) &
                                     +CU(1,KK-KD) &
                                     *D1(I-MM,KK-MM)*D1(J-NN,KD-NN))
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      DO MM=0,NXMAX-1
         DO I=MM,MM+1
            ID=2*I
            DO J=MM,MM+1
               JD=2*J
               IF(NWMAX.NE.NXMAX) JD=2*NWMAX+1+2*J-2*I
               DO KS=MM,MM+1
                  CK(JD+1      ,ID+1)=CK(JD+1      ,ID+1) &
                                     +CWP(KS)*CWE(KS)*CWE(KS)*DX &
                                     *D3(I-MM,J-MM,KS-MM)
                  CK(JD+2-NBAND,ID+2)=CK(JD+2-NBAND,ID+2) &
                                     +CWP(KS)*CWE(KS)*CWE(KS)*DX &
                                     *D3(I-MM,J-MM,KS-MM)
               END DO
            END DO
         END DO
      END DO

      DO IO=1,NWDDUB+3
         IF(NWMAX.NE.NXMAX) THEN
            IOB=2*NWMAX+4-IO
         ELSE
            IOB=2
         ENDIF
         IF(IOB.GE.1) CK(IOB,IO)=(0.D0,0.D0)
         CK(IO,2)=(0.D0,0.D0)
      END DO
      IF(NWMAX.NE.NXMAX) THEN
         I2=2*NWMAX+2
      ELSE
         I2=2
      ENDIF
      CK(I2,2)=(1.D0,0.D0)
      BETA=BETA0
      RETURN
    END SUBROUTINE SUBCK2

!     *****  CALCULATION OF RHS VECTOR  *****   

    SUBROUTINE SUBINI

      USE wicomm
      IMPLICIT NONE
      COMPLEX(rkind):: CBB
      INTEGER(ikind):: ML
      REAL(rkind):: ANB

      ANB=PN0*DEXP(-ALFA*xgrid(nxmax))
      CBB=CI/DSQRT(1.D0-ANB-ANY*ANY)
      DO ML=1,NXMAX*2+1
         CSO(ML)=(0.D0,0.D0)
      END DO
      CSO(NXMAX*2+2)=-CBB*CFYN
      CSO(NXMAX*2+3)=     CFYN
      RETURN
    END SUBROUTINE SUBINI

!     *****  SET FIELD (FULL MATRIX)  ***** 

    SUBROUTINE SUBFY

      USE wicomm
      IMPLICIT NONE
      INTEGER(ikind):: ML,MW

      DO ML=1,NXMAX*2+3
         CFY(ML)=(0.D0,0.D0)
         DO MW=1,NXMAX*2+3
            CFY(ML)=CFY(ML)+CK(MW,ML)*CSO(MW)
         END DO
      END DO
      RETURN
    END SUBROUTINE SUBFY

!     *****  SET FIELD (BAND MATRIX)  ***** 

    SUBROUTINE SUBFYW

      USE wicomm
      IMPLICIT NONE
      INTEGER(ikind):: ML

      DO ML=1,NXMAX*2+3
         CFY(ML)=CSO(ML)
      END DO
      RETURN
    END SUBROUTINE SUBFYW

!     *****  ABSORBED POWER  *****

    SUBROUTINE SUBPOW

      USE wicomm
      IMPLICIT NONE
      COMPLEX(rkind):: cp1,cp2,cp3,cp4,cpa,cpb
      INTEGER(ikind):: NX,ns,ne,nn,i,j,id,jd,kk,kd
      REAL(rkind):: rky,rky2,dx,dx2,AD,BD,BETA0

      RKY=ANY
      RKY2=RKY**2
      BETA0=BETA

      DO NX=0,NXMAX
         CPOWER(NX)=(0.D0,0.D0)
      END DO
      PTOT=0.D0

      DO NX=0,NXMAX-1
         IF(XMAX.GE.500.D0.AND.ALFA*XMAX.LT.10.D0) THEN
            IF(XMAX-XGRID(NX).LT.Bwidth) THEN
               BETA=BETA0*(XMAX-XGRID(NX))/Bwidth
            ELSE
               BETA=BETA0
            END IF
         ELSE
            BETA=BETA0
         END IF
         DX=xgrid(nx+1)-xgrid(nx)
         DX2=DX*DX
         NS=NX-NWMAX+1
         NE=NX+NWMAX-1
         AD=1.D0/(2.D0*DX) 
         BD=1.D0/(2.D0*DX) 
         IF(NX.EQ.0) AD=1.D0/DX
         IF(NX.EQ.NXMAX-1) BD=1.D0/DX
         IF(NS.LE.0) NS=0
         IF(NE.GE.NXMAX-1) NE=NXMAX-1
         DO NN=NS,NE
            DO I=NX,NX+1
               ID=2*I
               DO J=NN,NN+1
                  JD=2*J
                  DO KK=NX,NX+1
                     DO KD=NN,NN+1
                        CP1=DX2*RKY2*CU(1,KK-KD) &
                            *D0(I-NX,KK-NX)*D0(J-NN,KD-NN) &
                           +(CU(1,KK-KD)-CI*CU(2,KK-KD)) &
                            *D1(I-NX,KK-NX)*D1(J-NN,KD-NN)
                        CP2=-DX*RKY*CU(2,KK-KD) &
                            *D0(I-NX,KK-NX)*D1(J-NN,KD-NN)
                        CP3= DX*RKY*CU(2,KK-KD) &
                            *D1(I-NX,KK-NX)*D0(J-NN,KD-NN)
                        CP4=DX2*RKY2*(CU(1,KK-KD)-CI*CU(2,KK-KD)) &
                            *D0(I-NX,KK-NX)*D0(J-NN,KD-NN) &
                           +CU(1,KK-KD)*D1(I-NX,KK-NX)*D1(J-NN,KD-NN)
                        CP1=-CI*CP1
                        CP2=-CI*CP2
                        CP3=-CI*CP3
                        CP4=-CI*CP4
                        CPA=CWP(KD)*CWE(KK)*CWE(KD)*BETA &
                            *(CONJG(CFY(ID+1))*(CP1*CFY(JD+1) &
                                               +CP2*CFY(JD+2))  &
                             +CONJG(CFY(ID+2))*(CP3*CFY(JD+1) &
                                               +CP4*CFY(JD+2)))
                        CPB=CWP(KD)*CWE(KK)*CWE(KD)*BETA &
                            *(CONJG(CFY(ID+1))*(CONJG(CP1)*CFY(JD+1) &
                                               +CONJG(CP3)*CFY(JD+2))  &
                             +CONJG(CFY(ID+2))*(CONJG(CP2)*CFY(JD+1) &
                                               +CONJG(CP4)*CFY(JD+2)))
                        CPOWER(NX  )=CPOWER(NX  )+AD*0.5D0*(CPA+CPB)
                        CPOWER(NX+1)=CPOWER(NX+1)+BD*0.5D0*(CPA+CPB)
                        PTOT=PTOT+REAL(0.5D0*(CPA+CPB))
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      BETA=BETA0
      RETURN
    END SUBROUTINE SUBPOW

  END MODULE wihot
