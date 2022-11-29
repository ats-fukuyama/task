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
    REAL(qkind),INTENT(OUT):: ratea
    INTEGER(ikind),INTENT(OUT):: ierr

    mlmax=nxmax*2+3
    mwmax=4*nwmax+3

    CALL SUBFW    ! calculate elements of kernel function
    CALL SUBCK2   ! calculate coefficient matrix
    CALL SUBINI   ! calculate right-hand-side vector
    IF(NWMAX.EQ.NXMAX) THEN
       CALL INVMCQ(CK,mlmax,MLEN,IERR)   ! full matrix solver
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
       CALL BANDCQ(CK,CSO,mlmax,mwmax,MWID,IERR)   ! band matrix solver
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
      USE wigcom
      IMPLICIT NONE
      REAL(qkind):: dx,rky,x
      COMPLEX(qkind):: CS
      INTEGER(ikind):: J,L,NW

      DX=(XMAX-XMIN)/NXMAX
      RKY=ANY
      DO J=1,2
         N1=J
         DO NW=0,NWMAX
            X=NW*DX
            CALL EUL(X,RKY,CS,5,L)
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
      COMPLEX(qkind):: ciky,cbb
      REAL(qkind):: rky,rky2,dx,dx2,dky
      REAL(qkind):: ANB,beta0
      INTEGER(ikind):: NDUB,NBAND,NWDUB,NWDDUB,I,J,MM,ID,JD,NS,NE,NN
      INTEGER(ikind):: KK,KD,KS,IOB,IO,I2

      RKY=ANY
      RKY2=RKY**2
      DKY=ANY*ANY
      CIKY=CI*ANY
      ANB=EXP(-ALFA*xgrid(nxmax))
      CBB=CI/SQRT(1.D0-ANB-ANY*ANY)
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
      COMPLEX(qkind):: CBB
      INTEGER(ikind):: ML
      REAL(qkind):: ANB

      ANB=PN0*EXP(-ALFA*xgrid(nxmax))
      CBB=CI/SQRT(1.D0-ANB-ANY*ANY)
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
      COMPLEX(qkind):: cp1,cp2,cp3,cp4,cpa,cpb
      INTEGER(ikind):: NX,ns,ne,nn,i,j,id,jd,kk,kd
      REAL(qkind):: rky,rky2,dx,dx2,AD,BD,BETA0

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

!     *****  EULER TRANSFOMATION  *****

    SUBROUTINE EUL(X,RKY,CS,M,L)

      USE wicomm
      USE wigcom
      IMPLICIT NONE
      REAL(qkind),PARAMETER:: HP=0.5D0*PI
      INTEGER(ikind),PARAMETER:: LMAX=1000
      REAL(qkind),INTENT(IN):: X,RKY
      INTEGER(ikind),INTENT(IN):: M
      INTEGER(ikind),INTENT(INOUT):: L
      COMPLEX(qkind),INTENT(OUT):: CS
      REAL(qkind),DIMENSION(LMAX)::  A,B
      INTEGER(ikind):: ILST,K
      REAL(qkind):: H0,SR1,SI1,SR,SI,ESR,ESI,SR2,SI2,PARITY,SKR,SKI,BETA0

      BETA0=BETA
      IF(XMAX.GE.500.D0.AND.ALFA*XMAX.LT.10.D0) THEN
         IF(XMAX-X.LT.Bwidth) THEN
            BETA=BETA0*(XMAX-X)/Bwidth
         ELSE
            BETA=BETA0
         END IF
      ELSE
         BETA=BETA0
      END IF


      G2=HP
      G3=X/BETA
      IF(MODELA.EQ.0) THEN
         G4=0.5D0*ALFA*BETA
      ELSE
         G4=0.D0
      END IF
      G5=RKY*BETA
      H0=0.5D0
      ILST=0

      SR1=0.D0
      SI1=0.D0
      IF(M.NE.0) THEN 
         DO K=M-1,0,-1
            G1=DBLE(K)
            CALL DEFTC2(SR,SI,ESR,ESI,H0,EPS_KF,ILST)
            SR1=SR1+SR
            SI1=SI1+SI
         END DO
      ENDIF

      G1=DBLE(M)
      CALL DEFTC2(SR,SI,ESR,ESI,H0,EPS_KF,ILST)
      A(1)=SR
      SR2=0.5D0*SR
      B(1)=SI
      SI2=0.5D0*SI
      PARITY=-1.D0
      L=0

   30 CONTINUE
      L=L+1
      IF(L.GE.LMAX) GOTO 9000
      G1=DBLE(M+L)
      CALL DEFTC2(SR,SI,ESR,ESI,H0,EPS_KF,ILST)
      A(L+1)=SR*PARITY
      B(L+1)=SI*PARITY
      DO K=L,1,-1    
         A(K)=A(K+1)-A(K)
         B(K)=B(K+1)-B(K)
      END DO
      SKR=A(1)*PARITY*0.5D0**(L+1)
      SR2=SR2+SKR
      SKI=B(1)*PARITY*0.5D0**(L+1)
      SI2=SI2+SKI
      PARITY=-PARITY
      IF(ABS(SKR).GT.EPS_KF.OR.ABS(SKI).GT.EPS_KF) GOTO 30

      SR=(SR1+SR2)/SQRT(4.D0*G2)
      SI=(SI1+SI2)/SQRT(4.D0*G2)
      CS=CMPLX(SR,SI,qkind)

      BETA=BETA0

      RETURN

 9000 WRITE(6,*) ' ## DIMENSION OVERFLOW IN EULER TRANSFORMATION.'
      RETURN
    END SUBROUTINE EUL


!     *****  REAL PART  *****

    FUNCTION FUNR(X,XM,XP)

      USE wicomm,ONLY: ikind,qkind
      USE wigcom
      IMPLICIT NONE
      REAL(qkind),INTENT(IN):: X,XM,XP
      REAL(qkind):: FUNR
      REAL(qkind):: Y1,T,T2,YY,AN2

      Y1=XM
      IF(INT(G1).EQ.0) THEN 
         T=0.5D0*G2*XP
         T2=T*T
         YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
         IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
            IF(N1.EQ.1) THEN
               AN2=1.D0
            ELSEIF(N1.EQ.2) THEN
               AN2=T
            ENDIF
               FUNR=AN2*0.5*G2*EXP(YY)*COS(T)
         ELSE
            FUNR=0.D0
         ENDIF 
      ELSE
         T=G2*(X+2.D0*G1)
         T2=T*T
         YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
         IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
            IF(N1.EQ.1) THEN
               AN2=1.D0
            ELSEIF(N1.EQ.2) THEN
               AN2=T
            ENDIF
               FUNR=AN2*G2*EXP(YY)*COS(T)
         ELSE
            FUNR=0.D0
         ENDIF
      ENDIF
      RETURN
    END FUNCTION FUNR

!     *****  IMAG PART  *****

    FUNCTION FUNI(X,XM,XP)

      USE wicomm,ONLY: ikind,qkind
      USE wigcom
      IMPLICIT NONE
      REAL(qkind),INTENT(IN):: X,XM,XP
      REAL(qkind):: FUNI
      REAL(qkind):: Y1,Y2,T,T2,YY,AN2

      Y1=X
      Y2=XM
      T=G2*(XP+2.D0*G1)
      T2=T*T
      YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
      IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
         IF(N1.EQ.1) THEN
            AN2=1.D0
         ELSEIF(N1.EQ.2) THEN
            AN2=T
         ENDIF
         FUNI=AN2*G2*EXP(YY)*SIN(T)
      ELSE
         FUNI=0.D0
      ENDIF
      RETURN
    END FUNCTION FUNI

!     *****  DOUBLE EXPONENTIAL FORMULA  *****

    SUBROUTINE DEFTC2(CSR,CSI,ESR,ESI,H0,EPS,ILST)

!        FINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
!                    (-1.D0, +1.D0)
!         INTEGRAND SHOULD BE DEFINED BY FUNC(X,1-X,1+X)

      USE wicomm,ONLY: qkind,ikind
      IMPLICIT NONE
      REAL(qkind),INTENT(OUT):: CSR,CSI ! Integral
      REAL(qkind),INTENT(OUT):: ESR,ESI ! Estimated error
      REAL(qkind),INTENT(IN)::  H0      ! Initial step size
      REAL(qkind),INTENT(IN)::  EPS     ! Convergence thrshold
      INTEGER,INTENT(IN)::  ILST    ! print out control: 0 for no print out
!      INTERFACE
!         FUNCTION FUNR(X,XM,XP)
!           USE wicomm,ONLY: qkind,ikind
!           REAL(qkind):: FUNR
!           REAL(qkind),INTENT(IN):: X,XM,XP
!         END FUNCTION FUNR
!         FUNCTION FUNI(X,XM,XP)
!           USE wicomm,ONLY: qkind,ikind
!           REAL(qkind):: FUNI
!           REAL(qkind),INTENT(IN):: X,XM,XP
!         END FUNCTION FUNI
!      END INTERFACE
      REAL(qkind),PARAMETER:: HP=1.5707963267948966192D0

      REAL(qkind):: EPS1,H,X,CSRP,CSIP,ATPR,ATPI,ATMR,ATMI,EPSI
      REAL(qkind):: HN,HC,HS,CC,XM,XP,CTR,CTI,ATR,ATI
      INTEGER:: N,NP,NM,NMIN,IND,ND

      EPS1=EPS**0.75
      H=H0
      X=0.D0
      CSR=HP*H*FUNR(X,1.D0-X,1.D0+X)
      CSRP=0.D0
      CSI=HP*H*FUNI(X,1.D0-X,1.D0+X)
      CSIP=0.D0
      N=0
      NP=0
      NM=0
      NMIN=1

    5 IND=0
      ATPR=1.D0
      ATPI=1.D0
      ATMR=1.D0
      ATMI=1.D0
      ND=2
      EPSI=MAX(EPS1*H,2.D-17)
      IF(N.EQ.0) ND=1

   10 N=N+ND
      HN=DBLE(N)*H
      HC=HP*H*COSH(HN)
      IF(IND.NE.1) THEN
         HS=HP*SINH(-HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=EXP(-HS)*CC
         XP=EXP( HS)*CC
         CTR=HC*FUNR(X,XM,XP)*CC*CC
         CSR=CSR+CTR
         CTI=HC*FUNI(X,XM,XP)*CC*CC
         CSI=CSI+CTI
         NP=NP+1
         ATR=ATPR
         ATPR=ABS(CTR)
         ATI=ATPI
         ATPI=ABS(CTI)
         IF(N.GE.NMIN) THEN
            IF(SQRT((ATR+ATPR)**2+(ATI+ATPI)**2) &
                 .LT.EPSI*MAX(SQRT(CSR*CSR+CSI*CSI),1.D0)) THEN 
               IF(IND.EQ.-1) GO TO 100
               IND=1
            ENDIF
         ENDIF
      ENDIF

      IF(IND.NE.-1) THEN
         HS=HP*SINH( HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=EXP(-HS)*CC
         XP=EXP( HS)*CC
         CTR=HC*FUNR(X,XM,XP)*CC*CC
         CSR=CSR+CTR
         CTI=HC*FUNI(X,XM,XP)*CC*CC
         CSI=CSI+CTI
         NM=NM+1
         ATR=ATMR
         ATMR=ABS(CTR)
         ATI=ATMI
         ATMI=ABS(CTI)
         IF(N.GE.NMIN) THEN
            IF(SQRT((ATR+ATMR)**2+(ATI+ATMI)**2) &
                 .LT.EPSI*MAX(SQRT(CSR*CSR+CSI*CSI),1.D0)) THEN 
               IF(IND.EQ.1) GO TO 100
               IND=-1
            ENDIF
         ENDIF
      ENDIF
      GO TO 10

  100 CONTINUE
      ESR=ABS(CSR-CSRP)
      CSRP=CSR
      ESI=ABS(CSI-CSIP)
      CSIP=CSI
      IF(ILST.NE.0) THEN
         IF(H.GE.H0) THEN
            WRITE(6,601) H,NP,NM,CSR
            WRITE(6,604) CSI
         ENDIF
         IF(H.LT.H0) THEN
            WRITE(6,602) H,NP,NM,CSR,ESR
            WRITE(6,605) CSI,ESI
         ENDIF
      ENDIF
      IF(SQRT(ESR*ESR+ESI*ESI) &
        .LT.EPS1*MAX(SQRT(CSR*CSR+CSI*CSI),1.D0)) GOTO 200

      IF(N.GT.1000) THEN
         WRITE(6,*) 'XX DEFTC2: Loop count exceeds 1000'
         STOP
      ENDIF

      H=0.5D0*H
      CSR=0.5D0*CSR
      CSI=0.5D0*CSI
      NMIN=N/2
      N=-1
      GOTO 5

  200 RETURN

  601 FORMAT(1H ,1PD13.5,2I8,1PD24.15)
  604 FORMAT(1H ,13X,16X,1PD24.15)
  602 FORMAT(1H ,1PD13.5,2I8,1PD24.15,1PD14.5)
  605 FORMAT(1H ,13X,16X,1PD24.15,1PD14.5)
    END SUBROUTINE DEFTC2
  END MODULE wihot
