! $Id$

MODULE wifluid

  PRIVATE
  PUBLIC wi_fluid

CONTAINS

  SUBROUTINE wi_fluid(iprint,ratea,ierr)

    USE wicomm
    IMPLICIT NONE
    INTEGER(ikind),INTENT(IN):: iprint
    REAL(rkind),INTENT(OUT):: ratea
    INTEGER(ikind),INTENT(OUT):: ierr
    INTEGER(ikind):: mlmax,mwmax,ml,mw
    REAL:: GT1,GT2
    
    mlmax=nxmax*2+3
    mwmax=7
    
    CALL INITDS   ! initiaize FEM coefficients
       CALL GUTIME(GT1)
    CALL SUBFW    ! calculate elements of kernel function
       CALL GUTIME(GT2)
       IF(iprint > 0) WRITE(6,601) 'SUBFW ',GT2-GT1
    CALL SUBCK2   ! calculate coefficient matrix
    CALL SUBINI   ! calculate right-hand-side vector
       CALL GUTIME(GT2)
       IF(iprint > 0) WRITE(6,601) 'SUBCK2',GT2-GT1
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
       CALL GUTIME(GT2)
       IF(iprint > 0) WRITE(6,601) 'SOLVER',GT2-GT1
    CALL SUBPOW    ! calculate sbsorbed power
       CALL GUTIME(GT2)
       IF(iprint > 0) WRITE(6,601) 'POWER ',GT2-GT1
       RATEA=1.D0-ABS(CFY(NXMAX*2+3))**2
       IF(iprint > 0) WRITE(6,'(A,F8.5)') '## Absorption rate: ',RATEA
9900  CONTINUE
      RETURN
  601 FORMAT('## END OF ',A6,' ##  CPU TIME = ',F8.3,' SEC')
    END SUBROUTINE wi_fluid

!     *****  INITIALIZE D0,D1,D2,D3  *****

    SUBROUTINE INITDS

      USE wicomm
      IMPLICIT NONE
      INTEGER(ikind):: L,M,N

      D0(0,0)=1.D0/3.D0
      D0(0,1)=1.D0/6.D0
      D0(1,0)=1.D0/6.D0
      D0(1,1)=1.D0/3.D0
      D1(0,0)=-0.5D0
      D1(0,1)=-0.5D0
      D1(1,0)=0.5D0
      D1(1,1)=0.5D0
      D2(0,0)=1.D0
      D2(0,1)=-1.D0
      D2(1,0)=-1.D0
      D2(1,1)=1.D0
      DO L=0,1
         DO M=0,1
            DO N=0,1
               D3(L,M,N)=1.D0/12.D0
            END DO
         END DO
      END DO
      D3(0,0,0)=1.D0/4.D0
      D3(1,1,1)=1.D0/4.D0
      RETURN
    END SUBROUTINE INITDS

!     *****  SET KERNEL FUNCTION  ***** 

    SUBROUTINE SUBFW

      USE wicomm
      USE wigcom
      IMPLICIT NONE
      REAL(rkind):: dx,rky,x
      COMPLEX(rkind):: CS
      INTEGER(ikind):: J,L,NX,NW

      DX=XMAX/NXMAX
      DO NX=0,NXMAX
         CWE(NX)=DEXP(-ALFA*NX*DX)
         CWP(NX)=PN0
      END DO
      RETURN
    END SUBROUTINE SUBFW

!     *****  CALCULATION OF COEFFICIENT MATRIX  ***** 

    SUBROUTINE SUBCK2

      USE wicomm
      IMPLICIT NONE
      COMPLEX(rkind):: ciky,cbb
      REAL(rkind):: rky,rky2,dx,dx2,beta2,dky
      INTEGER(ikind):: MLMAX,MWMAX,ML,MW,I,J,NX,ID,JD
      INTEGER(ikind):: KK,KD,KS,IOB,IO,I2

      RKY=AKY*BETA
      RKY2=RKY**2
      DX=XMAX/DBLE(NXMAX)
      DX2=DX**2
      BETA2=BETA*BETA
      DKY=AKY*AKY
      CIKY=CI*AKY/BETA
      CBB=CI/(DSQRT(1.D0-AKY*AKY)*BETA)

      MLMAX=2*NXMAX+3
      MWMAX=7

      DO ML=1,MLMAX
         DO MW=1,MWMAX
            CK(MW,ML)=(0.D0,0.D0)
         END DO
      END DO

      DO NX=0,NXMAX-1
         DO I=NX,NX+1
            ID=2*I
            DO J=NX,NX+1
               JD=3+2*J-2*I
               CK(JD+1,ID+1)=CK(JD+1,ID+1) &
                            +(DKY-1.D0)*DX*D0(I-NX,J-NX)
               CK(JD+2,ID+1)=CK(JD+2,ID+1) &
                            +CIKY*D1(J-NX,I-NX)
               CK(JD  ,ID+2)=CK(JD  ,ID+2) &
                            -CIKY*D1(I-NX,J-NX)
               CK(JD+1,ID+2)=CK(JD+1,ID+2) &
                            +D2(I-NX,J-NX)/(DX*BETA2) &
                            -DX*D0(I-NX,J-NX)
            END DO
         END DO
      END DO
      CK(5,MLMAX-1)=-CBB
      CK(3,MLMAX)=1.D0
      CK(4,MLMAX)=-1.D0

      DO NX=0,NXMAX-1
         DO I=NX,NX+1
            ID=2*I
            DO J=NX,NX+1
               JD=3+2*J-2*I
               DO KS=NX,NX+1
                  CK(JD+1,ID+1)=CK(JD+1,ID+1) &
                                     +CWP(KS)*CWE(KS)*DX &
                                      *D3(I-NX,J-NX,KS-NX)/(1.D0+CI*PNU)
                  CK(JD+1,ID+2)=CK(JD+1,ID+2) &
                                     +CWP(KS)*CWE(KS)*DX &
                                      *D3(I-NX,J-NX,KS-NX)/(1.D0+CI*PNU)
               END DO
            END DO
         END DO
      END DO

      DO IO=1,7
         IOB=6-IO
         IF(IOB.GE.1) CK(IOB,IO)=(0.D0,0.D0)
         CK(IO,2)=(0.D0,0.D0)
      END DO
      I2=4
      CK(I2,2)=(1.D0,0.D0)
      RETURN
    END SUBROUTINE SUBCK2

!     *****  CALCULATION OF RHS VECTOR  *****   

    SUBROUTINE SUBINI

      USE wicomm
      IMPLICIT NONE
      COMPLEX(rkind):: CBB
      INTEGER(ikind):: ML

      CBB=CI/(DSQRT(1.D0-AKY*AKY)*BETA)

      DO ML=1,NXMAX*2+1
         CSO(ML)=(0.D0,0.D0)
      END DO
      CSO(NXMAX*2+2)=-CBB*CFYN
      CSO(NXMAX*2+3)=     CFYN
      RETURN
    END SUBROUTINE SUBINI

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
      COMPLEX(ikind):: cp1,cp2,cp3,cp4,cpa
      INTEGER(ikind):: NX,i,j,id,jd,kk,kd,ks
      REAL(rkind):: rky,rky2,dx,dx2,AD,BD

      RKY=AKY*BETA
      RKY2=RKY**2
      DX=XMAX/DBLE(NXMAX)
      DX2=DX**2

      DO NX=0,NXMAX
         CPOWER(NX)=(0.D0,0.D0)
      END DO
      PTOT=0.D0

      DO NX=0,NXMAX-1
         AD=1.D0/(2.D0*DX) 
         BD=1.D0/(2.D0*DX) 
         IF(NX.EQ.0) AD=1.D0/DX
         IF(NX.EQ.NXMAX-1) BD=1.D0/DX
         DO I=NX,NX+1
            ID=2*I
            DO J=NX,NX+1
               JD=2*J
               DO KS=NX,NX+1
                  CP1=-DX*D3(I-NX,J-NX,KS-NX)/(1.D0+CI*PNU)
                  CP4=-DX*D3(I-NX,J-NX,KS-NX)/(1.D0+CI*PNU)
                  CPA=CWP(KS)*CWE(KS)*CWE(KS) &
                       *(CONJG(CFY(ID+1))*CP1*CFY(JD+1) &
                        +CONJG(CFY(ID+2))*CP4*CFY(JD+2))
                  CPOWER(NX  )=CPOWER(NX  )-CI*AD*CPA
                  CPOWER(NX+1)=CPOWER(NX+1)-CI*BD*CPA
                  PTOT=PTOT-REAL(CI*CPA)
               END DO
            END DO
         END DO
      END DO
      RETURN
    END SUBROUTINE SUBPOW

  END MODULE wifluid
