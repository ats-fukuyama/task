MODULE decomm
  USE wimcomm, ONLY: rkind
  REAL(rkind):: G1,G2,G3
  INTEGER:: N1
END MODULE decomm

Module wimexec
  PRIVATE
  PUBLIC wim_exec
 
CONTAINS

  SUBROUTINE wim_exec(ierr)

    USE wimcomm
    
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    REAL(4):: RNT1,RNT2

    IERR=0
    CALL INITDS

    CALL GUTIME(RNT1)
    CALL DIVZ(NZMAX) 

    IF(NTMAX.NE.NTSET.OR.ABS(TMAX-TMAXST).GT.1.D-8) THEN
       CALL SUBFW
       CALL GUTIME(RNT2)
       WRITE(6,601) 'SUBFW ',RNT2-RNT1
       NTSET=NTMAX
       TMAXST=TMAX 
    ENDIF
    CALL INITBC
    SELECT CASE(MODELP)
    CASE(0)
       CALL SUBCK0(NZMAX,NWMAX,MODELW)
    CASE(2)
       CALL SUBCK2(NZMAX,NWMAX,MODELW)
    CASE DEFAULT
       WRITE(6,*) 'XX wimexec: unknown MODELP'
       RETURN
    END SELECT
    CALL SUBINI(NZMAX,CER1,CEL1,CER2,CEL2)
    CALL GUTIME(RNT2)
    WRITE(6,601) 'SUBCK2',RNT2-RNT1
    SELECT CASE(MODELP)
    CASE(0)
       CALL BANDCD(CK,CV,NZMAX*3+7,11,NWID,IERR)
       IF(IERR.NE.0) GOTO 9900
       CALL SUBFYW(NZMAX)
    CASE(2)   
       IF(NWMAX.EQ.NZMAX) THEN
          CALL INVMCD(CK,NZMAX*3+7,NZMAX*3+7,IERR)
          IF(IERR.NE.0) THEN
             WRITE(6,*) 'XX INVMCD ERROR: IERR=',IERR
             GOTO 9900
          END IF
          CALL SUBFY(NZMAX)
       ELSE
          CALL BANDCD(CK,CV,NZMAX*3+7,NWMAX*6+7,NWID,IERR)
          IF(IERR.NE.0) GOTO 9900
          CALL SUBFYW(NZMAX)
       ENDIF
    END SELECT
    CALL GUTIME(RNT2)
       WRITE(6,601) 'SOLVER',RNT2-RNT1

    SELECT CASE(MODELP)
    CASE(0)
       CALL SUBPOW0(NZMAX,NWMAX,CPTOT)
    CASE(2)
       CALL SUBPOW2(NZMAX,NWMAX,CPTOT)
    END SELECT
    CALL SUBKEI(NZMAX)
    CALL GUTIME(RNT2)
    WRITE(6,601) 'SUBPOW',RNT2-RNT1
9900 CONTINUE
    RETURN
601 FORMAT (1H ,'## END OF ',A6,' ## CPUTIME = ',F8.3,' SEC')
  END SUBROUTINE wim_exec


!     *****  INITIALIZE D0,D1,D2,D3  *****

  SUBROUTINE INITDS

    USE wimcomm, ONLY: D0,D1,D2,D3
    IMPLICIT NONE
    iNTEGER:: L,M,N

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

!     *****  DIVIDE Z-MESH  *****

  SUBROUTINE DIVZ(NZMAX)
    
    USE wimcomm, ONLY: rkind,ZA,PN0,DBDZ,ANX,BETA,ZMIN,ZMAX,DZMAX,DZWID
    IMPLICIT NONE
    INTEGER:: NZ,NZMAX
    REAL(rkind):: WT,Z,FACT,Z1,DZ

    DZ=(ZMAX-ZMIN)/NZMAX
    IF(DZMAX.EQ.0.D0) THEN
       DO NZ=0,NZMAX
          ZA(NZ)=ZMIN+DZ*NZ/NZMAX
       END DO
    ELSE
       WT=(ZMAX-ZMIN)+DZMAX*DZWID &
               *(ATAN(ZMAX/DZWID)-ATAN(ZMIN/DZWID))
       Z=ZMIN
       DO NZ=0,NZMAX-1
          ZA(NZ)=Z
          FACT=1.D0+DZMAX/(1.D0+(Z/DZWID)**2)
          Z1=Z+0.5D0*WT/(FACT*NZMAX)
          FACT=1.D0+DZMAX/(1.D0+(Z1/DZWID)**2)
          Z=Z+WT/(FACT*NZMAX)
       END DO
       WRITE(6,601) (ZMAX-Z)*NZMAX/(ZMAX-ZMIN)
601    FORMAT(1H ,'## ZDIV:   ERR = ',1PE12.4) 
       ZA(NZMAX)=ZMAX
    ENDIF
    RETURN
  END SUBROUTINE DIVZ

!     *****  CALCULATE KERNEL FUNCTION  *****

  SUBROUTINE SUBFW

    USE wimcomm,ONLY: rkind,TT,TMAX,NTMAX,CFK1,CFKS1,CFK2,CFKS2,PI
    USE decomm,ONLY: N1
    USE libgrf
    IMPLICIT NONE
    INTEGER:: IOPT(2),NT,L,IER
    REAL(rkind):: DT
    COMPLEX(rkind):: CDF(2)
    REAL(rkind):: F(NTMAX,4)

    DT=TMAX/(NTMAX-1) 
    DO NT=1,NTMAX 
       TT(NT)=(NT-1)*DT 
    END DO

    N1=1
    CFK1(1)=(0.D0,1.D0)/SQRT(2.D0*PI)
    DO NT=2,NTMAX
       CALL EUL(TT(NT),CFK1(NT),5,L)
    END DO
    CDF(1)=(-0.5D0,0.D0) 
    IOPT(1)=2
    IOPT(2)=3
    CALL DSPLC(TT,NTMAX,CFK1,CDF,IOPT,CFKS1,NTMAX,IER)

    N1=2
    CFK2(1)=-1.D0/SQRT(2.D0*PI)
    DO NT=2,NTMAX
       CALL EUL(TT(NT),CFK2(NT),5,L)
    END DO
    CDF(1)=(0.D0,0.D0) 
    IOPT(1)=2
    IOPT(2)=3
    CALL DSPLC(TT,NTMAX,CFK2,CDF,IOPT,CFKS2,NTMAX,IER)

!    CALL PAGES
!    F(1:NTMAX,1)=REAL(CFK1(1:NTMAX))
!    F(1:NTMAX,2)=IMAG(CFK1(1:NTMAX))
!    F(1:NTMAX,3)=REAL(CFK2(1:NTMAX))
!    F(1:NTMAX,4)=IMAG(CFK2(1:NTMAX))
!    CALL GRD1D(0,TT,F,NTMAX,NTMAX,4,'CFK')
!    CALL PAGEE

    RETURN
  END SUBROUTINE SUBFW

!     *****  CFN1( V )  *****

  FUNCTION CFN1(V)

    USE wimcomm,ONLY: rkind,PI,CI,TT,TMAX,NTMAX,CFK1,CFKS1
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: V
    REAL(rkind):: VA(1)
    COMPLEX(rkind):: CFN1,CFN1A(1)
    INTEGER:: IER
    REAL(rkind),PARAMETER:: R13=3.33333333333333D-1
    REAL(rkind),PARAMETER:: R23=6.66666666666667D-1
    REAL(rkind),PARAMETER:: SR3=1.73205080756887D0

    IF(V.LE.TMAX) THEN 
       VA(1)=V
       CALL DSPLF(TT,NTMAX,CFK1,CFKS1,NTMAX,VA,1,CFN1A,IER)
       CFN1=CFN1A(1)
    ELSE 
       CFN1=V**R13*EXP(-0.75D0*V**R23  &
                      +CI*(0.75D0*SR3*V**R23+R13*PI))/SR3 
    ENDIF
    RETURN
  END FUNCTION CFN1

!     *****  CFN2( V )  *****

  FUNCTION CFN2(V)
    USE wimcomm,ONLY: rkind,CI,TT,TMAX,NTMAX,CFK2,CFKS2,PI
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: V
    REAL(rkind):: VA(1)
    COMPLEX(rkind):: CFN2,CFN2A(1)
    INTEGER:: IER
    REAL(rkind),PARAMETER:: R23=6.66666666666667D-1
    REAL(rkind),PARAMETER:: SR3=1.73205080756887D0

    IF(V.LE.TMAX) THEN 
       VA(1)=V
       CALL DSPLF(TT,NTMAX,CFK2,CFKS2,NTMAX,VA,1,CFN2A,IER)
       CFN2=CFN2A(1)
    ELSE 
       CFN2=V*EXP(-0.75D0*V**R23 &
                 +CI*(0.75D0*SR3*V**R23+0.5D0*PI))/SR3 
    ENDIF
    RETURN
  END FUNCTION CFN2

!     *****  INITIALISE BOUNDARY CONDITION  *****

  SUBROUTINE INITBC

    USE wimcomm,ONLY: rkind,CI, &
                      CA1,CB1,CD1,CE1,CF1,CG1, &
                      CA2,CB2,CD2,CE2,CF2,CG2, &
                      PN0,DBDZ,ANX,BETA,ZMIN,ZMAX,DZMAX,DZWID
    IMPLICIT NONE
    REAL(rkind):: BETA2,PN,BB,S,D,P,AK3,AK5,AK4,AK1,AK2,TEMP
    COMPLEX(rkind):: CKK1,CKK2

    BETA2=BETA*BETA 
    PN=PN0*(1.D0+DBDZ*ZMIN)
    BB=1.D0+DBDZ*ZMIN
    IF(ABS(BB*BB-1.D0).LT.1.D-14) BB=1.D0+1.D-14
    S=1.D0-PN/(1.D0-BB**2)
    D=PN*BB/(1.D0-BB**2)
    P=1.D0-PN
    IF(ABS(P).LT.1.D-14) P=1.D-14
    AK3=-((S+P)*ANX*ANX-2*S*P)/(2*P) 
    AK5= ((S-P)*(S-P)*ANX**4-4*D*D*P*ANX**2+4*D*D*P*P) 
    IF(AK5.LT.0.D0) AK5=0.D0
    AK4=SQRT(AK5)/(2*P)
    AK1=BETA2*(AK3+AK4) 
    AK2=BETA2*(AK3-AK4) 
    IF(ABS(AK2).GT.ABS(AK1)) THEN
       TEMP=AK1
       AK1=AK2
       AK2=TEMP
    END IF
    CKK1=SQRT(CMPLX(AK1,0.D0))
    CKK2=SQRT(CMPLX(AK2,0.D0))

    CA1=CI*P*CKK1/(BETA2*(P-ANX*ANX)) 
    CB1=CI*P*CKK2/(BETA2*(P-ANX*ANX)) 
    CD1=-CKK1*D/(BETA2*(ANX*ANX+AK1/BETA2-S)) 
    CE1=-CKK2*D/(BETA2*(ANX*ANX+AK2/BETA2-S)) 
    CF1=CI*D/(ANX*ANX+AK1/BETA2-S)
    CG1=CI*D/(ANX*ANX+AK2/BETA2-S)

    PN=PN0*(1.D0+DBDZ*ZMAX)
    BB=1.D0+DBDZ*ZMAX
    IF(ABS(BB*BB-1.D0).LT.1.D-14) BB=1.D0+1.D-14
    S=1.D0-PN/(1.D0-BB**2)
    D=PN*BB/(1.D0-BB**2)
    P=1.D0-PN
    IF(ABS(P).LT.1.D-14) P=1.D-14
    AK3=-((S+P)*ANX*ANX-2*S*P)/(2*P) 
    AK5= ((S-P)*(S-P)*ANX**4-4*D*D*P*ANX**2+4*D*D*P*P) 
    IF(AK5.LT.0.D0) AK5=0.D0
    AK4=SQRT(AK5)/(2*P)
    AK1=BETA2*(AK3+AK4) 
    AK2=BETA2*(AK3-AK4) 
    IF(ABS(AK2).GT.ABS(AK1)) THEN
       TEMP=AK1
       AK1=AK2
       AK2=TEMP
    END IF
    CKK1=SQRT(CMPLX(AK1,0.D0))
    CKK2=SQRT(CMPLX(AK2,0.D0))

    CA2=CI*P*CKK1/(BETA2*(P-ANX*ANX)) 
    CB2=CI*P*CKK2/(BETA2*(P-ANX*ANX)) 
    CD2=-CKK1*D/(BETA2*(ANX*ANX+AK1/BETA2-S)) 
    CE2=-CKK2*D/(BETA2*(ANX*ANX+AK2/BETA2-S)) 
    CF2=CI*D/(ANX*ANX+AK1/BETA2-S)
    CG2=CI*D/(ANX*ANX+AK2/BETA2-S)
    RETURN
  END SUBROUTINE INITBC

!     *****  CALCULATE COEFFICIENT MATRIX  ***** 

  SUBROUTINE SUBCK2(NZMAX,NWMAX,MODELW) 

    USE wimcomm,ONLY: rkind,CI,PN0,DBDZ,ANX,BETA,ZMIN,ZMAX,DZMAX,DZWID,ZA, &
                      CK,CV,CA1,CB1,CD1,CE1,CF1,CG1,CA2,CB2,CD2,CE2,CF2,CG2, &
                      D0,D1,D2,D3
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZMAX,NWMAX,MODELW
    INTEGER:: NZMAX3,NBAND,NWMAX3R,NWMAX6,NWMAX3L, &
              I,J,MM,ID,JD,NS,NE,NN,KI,KJ,I1,I2
    REAL(rkind):: RKX,RKX2,BETA2,ANX2,DZ,DZI,DZJ,YY,YYI,YYJ,VP,VM,VPI,VPJ, &
                  VMI,VMJ,VVP,VVM,VVZ,X1,X4,X5,DPI,DMI,DPJ,DMJ,DPM,DP,DM,DE, &
                  VY,DS
    COMPLEX(rkind):: CIKX,CF1P,CF1M,CF2Z,CF1Z,CP1,CP2,CP3
 
    RKX=ANX*BETA
    RKX2=RKX**2
    NZMAX3=3*NZMAX
    BETA2=BETA*BETA
    ANX2=ANX*ANX
    CIKX=CI*ANX/BETA
 
    IF(NWMAX.EQ.NZMAX) THEN
       NBAND=0
       NWMAX3R=NZMAX3
       NWMAX6=3*NZMAX
       NWMAX3L=0
    ELSE
       NBAND=1
       NWMAX3R=3*NWMAX
       NWMAX6=6*NWMAX
       NWMAX3L=3*NWMAX
    ENDIF

    DO I=1,NZMAX3+7
       DO J=1,NWMAX6+7
          CK(J,I)=(0.D0,0.D0)
       END DO
    END DO

    DO MM=0,NZMAX-1
       DZ=ZA(MM+1)-ZA(MM)
       DO I=MM,MM+1
          ID=3*I+2
          DO J=MM,MM+1
             JD=3*J+2
             IF(NWMAX.NE.NZMAX) JD=3*NWMAX+4+3*J-3*I
             CK(JD+1,ID+1)=CK(JD+1,ID+1) &
                          +D2(I-MM,J-MM)/(BETA2*DZ)  &
                          -DZ*D0(I-MM,J-MM)
             CK(JD+3,ID+1)=CK(JD+3,ID+1) &
                          -CIKX*D1(I-MM,J-MM)
             CK(JD+2-NBAND,ID+2)=CK(JD+2-NBAND,ID+2) &
                                +D2(I-MM,J-MM)/(BETA2*DZ) &
                                +(ANX2-1.D0)*DZ*D0(I-MM,J-MM)
             CK(JD+1-2*NBAND,ID+3)=CK(JD+1-2*NBAND,ID+3) &
                                  +CIKX*D1(J-MM,I-MM)
             CK(JD+3-2*NBAND,ID+3)=CK(JD+3-2*NBAND,ID+3) &
                                  +(ANX2-1.D0)*DZ*D0(I-MM,J-MM)
          END DO
       END DO
    END DO

    CK(NWMAX3L+1+2*NBAND,1)=-1.D0
    CK(NWMAX3L+2+2*NBAND,1)=-1.D0
    CK(NWMAX3L+3+2*NBAND,1)= 1.D0
    CK(NWMAX3L+1+  NBAND,2)=-CF1
    CK(NWMAX3L+2+  NBAND,2)=-CG1
    CK(NWMAX3L+4+  NBAND,2)= 1.D0
    CK(NWMAX3L+1        ,3)=-CA1
    CK(NWMAX3L+2        ,3)=-CB1
    CK(NWMAX3L+1-  NBAND,4)=-CD1
    CK(NWMAX3L+2-  NBAND,4)=-CE1
    CK(NWMAX3R+6+2*NBAND,NZMAX3+3)=-CA2
    CK(NWMAX3R+7+2*NBAND,NZMAX3+3)=-CB2
    CK(NWMAX3R+6+  NBAND,NZMAX3+4)=-CD2
    CK(NWMAX3R+7+  NBAND,NZMAX3+4)=-CE2
    CK(NWMAX3R+3-  NBAND,NZMAX3+6)= 1.D0
    CK(NWMAX3R+6-  NBAND,NZMAX3+6)=-1.D0
    CK(NWMAX3R+7-  NBAND,NZMAX3+6)=-1.D0
    CK(NWMAX3R+4-2*NBAND,NZMAX3+7)= 1.D0
    CK(NWMAX3R+6-2*NBAND,NZMAX3+7)=-CF2
    CK(NWMAX3R+7-2*NBAND,NZMAX3+7)=-CG2

    DO MM=0,NZMAX-1
       DZI=ZA(MM+1)-ZA(MM)
       NS=MM-NWMAX
       NE=MM+NWMAX-1
       IF(NS.LE.0) NS=0
       IF(NE.GE.NZMAX-1) NE=NZMAX-1
       DO NN=NS,NE
          DZJ=ZA(NN+1)-ZA(NN)
          DO KI=MM,MM+1
             DO KJ=NN,NN+1
                YY =1.D0+0.5D0*DBDZ*(ZA(KJ)+ZA(KI)) 
                YYI=1.D0+DBDZ*ZA(KI) 
                YYJ=1.D0+DBDZ*ZA(KJ) 
                VP =1.D0+YY
                VM =1.D0-YY
                VPI=1.D0+YYI
                VPJ=1.D0+YYJ
                VMI=1.D0-YYI
                VMJ=1.D0-YYJ
                VVP=ABS(VP*(ZA(KI)-ZA(KJ)))
                VVM=ABS(VM*(ZA(KI)-ZA(KJ)))
                VVZ=ABS(   (ZA(KI)-ZA(KJ))) 
                X1=(ABS(YYI*YYJ))**1.5/(YY*YY)
                X4=(YYI*YYJ)/YY
                X5=DBDZ*DBDZ/(2*YY*YY)
                DPI=1.5D0*DBDZ/YYI-DBDZ/YY-DBDZ/VPI
                DMI=1.5D0*DBDZ/YYI-DBDZ/YY+DBDZ/VMI
                DPJ=1.5D0*DBDZ/YYJ-DBDZ/YY-DBDZ/VPJ
                DMJ=1.5D0*DBDZ/YYJ-DBDZ/YY+DBDZ/VMJ
                DPM=0.5D0*DBDZ*DBDZ/(YY*YY)
                IF(VP.GT.0.D0) THEN
                   CF1P=-CI*CFN1(VVP)/(VPI*VPJ) 
                ELSE
                   CF1P= CI*CONJG(CFN1(VVP))/(VPI*VPJ) 
                ENDIF
                IF(VM.GT.0.D0) THEN
                   CF1M=-CI*CFN1(VVM)/(VMI*VMJ) 
                ELSE
                   CF1M= CI*CONJG(CFN1(VVM))/(VMI*VMJ) 
                ENDIF
                CF2Z= CFN2(VVZ)
                CF1Z=-CI*CFN1(VVZ)
                DO I=MM,MM+1
                   ID=3*I+2
                   DO J=NN,NN+1
                      JD=3*J+2
                      IF(NWMAX.NE.NZMAX) JD=3*NWMAX+4+3*J-3*I
                      IF(JD.LE.0) THEN
                         WRITE(6,*) 'JD=',JD
                         WRITE(6,*) 'NN,NS,NE=',NN,NS,NE
                         WRITE(6,*) 'MM=',MM
                         WRITE(6,*) 'I=',I
                         WRITE(6,*) 'J=',J
                         WRITE(6,*) 'NWMAX=',NWMAX
                         WRITE(6,*) 'NZMAX=',NZMAX
                         STOP
                      END IF
                      DP=D1(I-MM,KI-MM)*D1(J-NN,KJ-NN) &
                        +D0(I-MM,KI-MM)*D1(J-NN,KJ-NN)*DZI*DPI &
                        +D1(I-MM,KI-MM)*D0(J-NN,KJ-NN)*DZJ*DPJ &
                        +D0(I-MM,KI-MM)*D0(J-NN,KJ-NN)*DZI*DZJ*(DPI*DPJ+DPM) 
                      DM=D1(I-MM,KI-MM)*D1(J-NN,KJ-NN) &
                        +D0(I-MM,KI-MM)*D1(J-NN,KJ-NN)*DZI*DMI &
                        +D1(I-MM,KI-MM)*D0(J-NN,KJ-NN)*DZJ*DMJ &
                        +D0(I-MM,KI-MM)*D0(J-NN,KJ-NN)*DZI*DZJ*(DMI*DMJ+DPM)
                      DE=D0(I-MM,KI-MM)*D0(J-NN,KJ-NN)*DZI*DZJ
                      CP1=-0.5D0*PN0*CI*X1*(CF1P*DP+CF1M*DM)
                      CP2=-0.5D0*PN0   *X1*(CF1P*DP-CF1M*DM)
                      CP3=-      PN0*CI*X4*(CF1Z*DE-X5*CF2Z*DE) 
                      CK(JD+1,ID+1)        =CK(JD+1,ID+1)        +CP1 
                      CK(JD+2,ID+1)        =CK(JD+2,ID+1)        +CP2 
                      CK(JD+1-NBAND,ID+2)  =CK(JD+1-NBAND,ID+2)  -CP2 
                      CK(JD+2-NBAND,ID+2)  =CK(JD+2-NBAND,ID+2)  +CP1 
                      CK(JD+3-2*NBAND,ID+3)=CK(JD+3-2*NBAND,ID+3)+CP3 
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO

    DO MM=0,NZMAX-1
       DZ=ZA(MM+1)-ZA(MM) 
       DO I=MM,MM+1
          ID=3*I+2
          DO J=MM,MM+1
             JD=3*J+2
             IF(NWMAX.NE.NZMAX) JD=3*NWMAX+2+3*J-3*I
             DO KI=MM,MM+1
                YY=1.D0+DBDZ*ZA(KI) 
                VY=YY 
                DS=DZ*D3(I-MM,J-MM,KI-MM)
                CK(JD+1,ID+1)=CK(JD+1,ID+1) &
                             +DS*PN0*YY/(1.D0-VY*VY) 
                CK(JD+2,ID+1)=CK(JD+2,ID+1) &
                             +DS*CI*PN0*YY*VY/(1.D0-VY*VY) 
                CK(JD+1-NBAND,ID+2)=CK(JD+1-NBAND,ID+2) &
                                   -DS*CI*PN0*YY*VY/(1.D0-VY*VY)
                CK(JD+2-NBAND,ID+2)=CK(JD+2-NBAND,ID+2) &
                                   +DS*PN0*YY/(1.D0-VY*VY) 
             END DO
          END DO
       END DO
    END DO

    IF(MODELW.EQ.1) THEN
       DO I=1,NWMAX6+7
          IF(NWMAX.NE.NZMAX) THEN
             I1=NWMAX*3+4-I
             I2=NWMAX*3+5-I
          ELSE
             I1=3
             I2=4
          ENDIF
          IF(I1.GE.1) CK(I1,I)=(0.D0,0.D0)
          CK(I,3)=(0.D0,0.D0)
          IF(I2.GE.1) CK(I2,I)=(0.D0,0.D0)
          CK(I,4)=(0.D0,0.D0)
       END DO
       IF(NWMAX.NE.NZMAX) THEN
          I1=NWMAX*3+3
          I2=NWMAX*3+3
       ELSE
          I1=3
          I2=4
       ENDIF
       CK(I1,3)=(1.D0,0.D0)
       CK(I2,4)=(1.D0,0.D0)
    ENDIF
    RETURN
  END SUBROUTINE SUBCK2

!     *****  CALCULATE COEFFICIENT MATRIX  ***** 

  SUBROUTINE SUBCK0(NZMAX,NWMAX,MODELW) 

    USE wimcomm,ONLY: rkind,CI,PN0,DBDZ,ANX,BETA,ZMIN,ZMAX,DZMAX,DZWID,ZA, &
                      CK,CV,CA1,CB1,CD1,CE1,CF1,CG1,CA2,CB2,CD2,CE2,CF2,CG2, &
                      D0,D1,D2,D3,PZCL 
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZMAX,NWMAX,MODELW
    INTEGER:: NZMAX3,NBAND,NWMAX3R,NWMAX6,NWMAX3L, &
              I,J,MM,ID,JD,NN,KI,I1,I2
    REAL(rkind):: RKX,RKX2,BETA2,ANX2,DZ,DP,DM,DE,DS
    COMPLEX(rkind):: CIKX,CP1,CP2,CP3,CCOL,YY,VY
 
    RKX=ANX*BETA
    RKX2=RKX**2
    NZMAX3=3*NZMAX
    BETA2=BETA*BETA
    ANX2=ANX*ANX
    CIKX=CI*ANX/BETA
 
    NBAND=1
    NWMAX3R=3*NWMAX
    NWMAX6=6*NWMAX
    NWMAX3L=3*NWMAX

    DO I=1,NZMAX3+7
       DO J=1,NWMAX6+7
          CK(J,I)=(0.D0,0.D0)
       END DO
    END DO

    DO MM=0,NZMAX-1
       DZ=ZA(MM+1)-ZA(MM)
       DO I=MM,MM+1
          ID=3*I+2
          DO J=MM,MM+1
             JD=3*NWMAX+4+3*J-3*I
             CK(JD+1,ID+1)=CK(JD+1,ID+1) &
                          +D2(I-MM,J-MM)/(BETA2*DZ)  &
                          -DZ*D0(I-MM,J-MM)
             CK(JD+3,ID+1)=CK(JD+3,ID+1) &
                          -CIKX*D1(I-MM,J-MM)
             CK(JD+2-NBAND,ID+2)=CK(JD+2-NBAND,ID+2) &
                                +D2(I-MM,J-MM)/(BETA2*DZ) &
                                +(ANX2-1.D0)*DZ*D0(I-MM,J-MM)
             CK(JD+1-2*NBAND,ID+3)=CK(JD+1-2*NBAND,ID+3) &
                                  +CIKX*D1(J-MM,I-MM)
             CK(JD+3-2*NBAND,ID+3)=CK(JD+3-2*NBAND,ID+3) &
                                  +(ANX2-1.D0)*DZ*D0(I-MM,J-MM)
          END DO
       END DO
    END DO

    CK(NWMAX3L+1+2*NBAND,1)=-1.D0
    CK(NWMAX3L+2+2*NBAND,1)=-1.D0
    CK(NWMAX3L+3+2*NBAND,1)= 1.D0
    CK(NWMAX3L+1+  NBAND,2)=-CF1
    CK(NWMAX3L+2+  NBAND,2)=-CG1
    CK(NWMAX3L+4+  NBAND,2)= 1.D0
    CK(NWMAX3L+1        ,3)=-CA1
    CK(NWMAX3L+2        ,3)=-CB1
    CK(NWMAX3L+1-  NBAND,4)=-CD1
    CK(NWMAX3L+2-  NBAND,4)=-CE1
    CK(NWMAX3R+6+2*NBAND,NZMAX3+3)=-CA2
    CK(NWMAX3R+7+2*NBAND,NZMAX3+3)=-CB2
    CK(NWMAX3R+6+  NBAND,NZMAX3+4)=-CD2
    CK(NWMAX3R+7+  NBAND,NZMAX3+4)=-CE2
    CK(NWMAX3R+3-  NBAND,NZMAX3+6)= 1.D0
    CK(NWMAX3R+6-  NBAND,NZMAX3+6)=-1.D0
    CK(NWMAX3R+7-  NBAND,NZMAX3+6)=-1.D0
    CK(NWMAX3R+4-2*NBAND,NZMAX3+7)= 1.D0
    CK(NWMAX3R+6-2*NBAND,NZMAX3+7)=-CF2
    CK(NWMAX3R+7-2*NBAND,NZMAX3+7)=-CG2

    CCOL=1.D0+CI*PZCL
    DO MM=0,NZMAX-1
       DZ=ZA(MM+1)-ZA(MM) 
       DO I=MM,MM+1
          ID=3*I+2
          DO J=MM,MM+1
             JD=3*J+2
             IF(NWMAX.NE.NZMAX) JD=3*NWMAX+2+3*J-3*I
             DO KI=MM,MM+1
                YY=(1.D0+DBDZ*ZA(KI))/CCOL
                VY=YY/CCOL
                DS=DZ*D3(I-MM,J-MM,KI-MM)
                CK(JD+1,ID+1)=CK(JD+1,ID+1) &
                             +DS*PN0*YY/(1.D0-VY*VY) 
                CK(JD+2,ID+1)=CK(JD+2,ID+1) &
                             +DS*CI*PN0*YY*VY/(1.D0-VY*VY) 
                CK(JD+1-NBAND,ID+2)=CK(JD+1-NBAND,ID+2) &
                                   -DS*CI*PN0*YY*VY/(1.D0-VY*VY)
                CK(JD+2-NBAND,ID+2)=CK(JD+2-NBAND,ID+2) &
                                   +DS*PN0*YY/(1.D0-VY*VY) 
                CK(JD+3-2*NBAND,ID+3)=CK(JD+2-2*NBAND,ID+3) &
                                     +DS*PN0*YY
             END DO
          END DO
       END DO
    END DO

    IF(MODELW.EQ.1) THEN
       DO I=1,NWMAX6+7
          IF(NWMAX.NE.NZMAX) THEN
             I1=NWMAX*3+4-I
             I2=NWMAX*3+5-I
          ELSE
             I1=3
             I2=4
          ENDIF
          IF(I1.GE.1) CK(I1,I)=(0.D0,0.D0)
          CK(I,3)=(0.D0,0.D0)
          IF(I2.GE.1) CK(I2,I)=(0.D0,0.D0)
          CK(I,4)=(0.D0,0.D0)
       END DO
       IF(NWMAX.NE.NZMAX) THEN
          I1=NWMAX*3+3
          I2=NWMAX*3+3
       ELSE
          I1=3
          I2=4
       ENDIF
       CK(I1,3)=(1.D0,0.D0)
       CK(I2,4)=(1.D0,0.D0)
    ENDIF
    RETURN
  END SUBROUTINE SUBCK0

!     *****  CALCULATE RHS VECTOR  *****

  SUBROUTINE SUBINI(NZMAX,CER1,CEL1,CER2,CEL2)

    USE wimcomm,ONLY: rkind,CV,CA1,CB1,CD1,CE1,CF1,CG1,CA2,CB2,CD2,CE2,CF2,CG2
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZMAX
    COMPLEX(rkind),INTENT(IN):: CER1,CEL1,CER2,CEL2
    REAL(rkind):: FACTOR
    INTEGER:: NZMAX3,K1,I

    NZMAX3=3*NZMAX 
    DO K1=6,NZMAX3+2
       CV(K1)=(0.D0,0.D0)
    END DO
    FACTOR=1.D0/SQRT(2.D0)
    CV(1)=FACTOR*(     CER1+    CEL1)
    CV(2)=FACTOR*( CF1*CER1+CG1*CEL1)
    CV(3)=FACTOR*(-CA1*CER1-CB1*CEL1)
    CV(4)=FACTOR*(-CD1*CER1-CE1*CEL1)
    CV(5)=0.D0
    CV(NZMAX3+3)=FACTOR*(-CA2*CER2-CB2*CEL2)
    CV(NZMAX3+4)=FACTOR*(-CD2*CER2-CE2*CEL2)
    CV(NZMAX3+5)=0.D0
    CV(NZMAX3+6)=FACTOR*(     CER2+    CEL2)
    CV(NZMAX3+7)=FACTOR*( CF2*CER2+CG2*CEL2)

!    DO I=1,5
!       WRITE(6,'(1P4E12.4)') CV(I),CV(NZMAX3+I+2)
!    END DO
    RETURN
  END SUBROUTINE SUBINI

!     *****  SET ELECTRIC FIELD (FULL MATRIX)  ***** 

  SUBROUTINE SUBFY(NZMAX)

    USE wimcomm,ONLY: CE,CK,CV
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZMAX
    INTEGER:: I,J

    DO I=1,3*NZMAX+7
       CE(I)=(0.D0,0.D0)
       DO J=1,3*NZMAX+7
          CE(I)=CE(I)+CK(J,I)*CV(J)
       END DO
    END DO
!    DO I=1,NZMAX
!       WRITE(6,'(I5,1P6E12.4)') I,CE(3*I+3),CE(3*I+4),CE(3*I+5)
!    END DO
    RETURN
  END SUBROUTINE SUBFY

!     *****  SET ELECTRIC FIELD (BAND MATRIX)  ***** 

  SUBROUTINE SUBFYW(NZMAX)

    USE wimcomm,ONLY: CE,CV
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZMAX
    INTEGER:: I

    DO I=1,3*NZMAX+7
       CE(I)=CV(I)
    END DO
    RETURN
  END SUBROUTINE SUBFYW

!     *****  CALCULATE ABSORBED POWER  ***** 

  SUBROUTINE SUBPOW2(NZMAX,NWMAX,CPTOT)

    USE wimcomm,ONLY: rkind,CI,ZA,CE,PN0,DBDZ,ANX,BETA,ZMIN,ZMAX,DZMAX,DZWID, &
                      CPWR,D0,D1,D2,D3
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZMAX,NWMAX
    COMPLEX(rkind),INTENT(OUT):: CPTOT
    INTEGER:: NZ,MM,NN,NS,NE,I,J,ID,JD,KI,KJ
    REAL(rkind):: DZI,DZJ,AD,BD,YY,YYI,YYJ,VP,VM,VPI,VPJ,VMI,VMJ, &
                  VVP,VVM,VVZ,X1,X4,X5,DPI,DMI,DPJ,DMJ,DPM,DP,DM,DE
    COMPLEX(rkind):: CF1P,CF1M,CF2Z,CF1Z,CP1,CP2,CP3,CPA,CPB
     
    CPTOT=0.D0
    DO NZ=0,NZMAX
       CPWR(NZ)=(0.D0,0.D0)
    END DO

    DO MM=0,NZMAX-1
       DZI=ZA(MM+1)-ZA(MM)
       NS=MM-NWMAX
       NE=MM+NWMAX-1
       IF(NS.LE.0) NS=0
       IF(NE.GE.NZMAX-1) NE=NZMAX-1
       AD=0.5D0/DZI 
       BD=0.5D0/DZI 
       IF(MM.EQ.0)    AD=1.D0/DZI 
       IF(MM.EQ.NZMAX-1) BD=1.D0/DZI 
       DO NN=NS,NE
          DZJ=ZA(NN+1)-ZA(NN)
          DO KI=MM,MM+1
             DO KJ=NN,NN+1
                YY =1.D0+0.5D0*DBDZ*(ZA(KJ)+ZA(KI)) 
                YYI=1.D0+DBDZ*ZA(KI) 
                YYJ=1.D0+DBDZ*ZA(KJ) 
                VP =1.D0+YY
                VM =1.D0-YY
                VPI=1.D0+YYI
                VPJ=1.D0+YYJ
                VMI=1.D0-YYI
                VMJ=1.D0-YYJ
                VVP=ABS(VP*(ZA(KI)-ZA(KJ)))
                VVM=ABS(VM*(ZA(KI)-ZA(KJ)))
                VVZ=ABS(   (ZA(KI)-ZA(KJ))) 
                X1=(YYI*YYJ)**1.5/(YY*YY)
                X4=(YYI*YYJ)/YY
                X5=DBDZ*DBDZ/(2*YY*YY)
                DPI=1.5D0*DBDZ/YYI-DBDZ/YY-DBDZ/VPI
                DMI=1.5D0*DBDZ/YYI-DBDZ/YY+DBDZ/VMI
                DPJ=1.5D0*DBDZ/YYJ-DBDZ/YY-DBDZ/VPJ
                DMJ=1.5D0*DBDZ/YYJ-DBDZ/YY+DBDZ/VMJ
                DPM=0.5D0*DBDZ*DBDZ/(YY*YY)
                IF(VP.GT.0.D0) THEN
                   CF1P=-CI*CFN1(VVP)/(VPI*VPJ) 
                ELSE
                   CF1P= CI*CONJG(CFN1(VVP))/(VPI*VPJ) 
                ENDIF
                IF(VM.GT.0.D0) THEN
                   CF1M=-CI*CFN1(VVM)/(VMI*VMJ) 
                ELSE
                   CF1M= CI*CONJG(CFN1(VVM))/(VMI*VMJ) 
                ENDIF
                CF2Z= CFN2(VVZ)
                CF1Z=-CI*CFN1(VVZ)
                DO I=MM,MM+1
                   ID=3*I+2
                   DO J=NN,NN+1
                      JD=3*J+2
                      IF(NWMAX.NE.NZMAX) JD=3*NWMAX+4+3*J-3*I
                      DP=D1(I-MM,KI-MM)*D1(J-NN,KJ-NN) &
                        +D0(I-MM,KI-MM)*D1(J-NN,KJ-NN)*DZI*DPI &
                        +D1(I-MM,KI-MM)*D0(J-NN,KJ-NN)*DZJ*DPJ &
                        +D0(I-MM,KI-MM)*D0(J-NN,KJ-NN)*DZI*DZJ*(DPI*DPJ+DPM) 
                      DM=D1(I-MM,KI-MM)*D1(J-NN,KJ-NN) &
                        +D0(I-MM,KI-MM)*D1(J-NN,KJ-NN)*DZI*DMI &
                        +D1(I-MM,KI-MM)*D0(J-NN,KJ-NN)*DZJ*DMJ &
                        +D0(I-MM,KI-MM)*D0(J-NN,KJ-NN)*DZI*DZJ*(DMI*DMJ+DPM) 
                      DE=D0(I-MM,KI-MM)*D0(J-NN,KJ-NN)*DZI*DZJ
                      CP1=-0.5D0*PN0*CI*X1*(CF1P*DP+CF1M*DM)
                      CP2=-0.5D0*PN0    *X1*(CF1P*DP-CF1M*DM)
                      CP3=-      PN0*CI*X4*(CF1Z*DE-X5*CF2Z*DE) 
                      CPA=CONJG(CE(ID+1))*( CP1*CE(JD+1)+CP2*CE(JD+2)) &
                         +CONJG(CE(ID+2))*(-CP2*CE(JD+1)+CP1*CE(JD+2)) &
                         +CONJG(CE(ID+3))*  CP3*CE(JD+3)
                      CPB=CE(ID+1)*CONJG( CP1*CE(JD+1)-CP2*CE(JD+2)) &
                         +CE(ID+2)*CONJG( CP2*CE(JD+1)+CP1*CE(JD+2)) &
                         +CE(ID+3)*CONJG( CP3*CE(JD+3))
                      CPWR(MM)  =CPWR(MM)  +CI*AD*0.5D0*(CPA+CPB)
                      CPWR(MM+1)=CPWR(MM+1)+CI*BD*0.5D0*(CPA+CPB)
                      CPTOT     =CPTOT     +CI   *0.5D0*(CPA+CPB)
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE SUBPOW2

!     *****  CALCULATE ABSORBED POWER  ***** 

  SUBROUTINE SUBPOW0(NZMAX,NWMAX,CPTOT)

    USE wimcomm,ONLY: rkind,CI,ZA,CE,PN0,DBDZ,ANX,BETA,ZMIN,ZMAX,DZMAX,DZWID, &
                      CPWR,D0,D1,D2,D3,PZCL
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZMAX,NWMAX
    COMPLEX(rkind),INTENT(OUT):: CPTOT
    INTEGER:: NZ,MM,I,J,ID,JD,KI
    REAL(rkind):: DZ,DS,AD,BD
    COMPLEX(rkind):: CP1,CP2,CP3,CPA,CCOL,YY,VY
     
    CPTOT=0.D0
    DO NZ=0,NZMAX
       CPWR(NZ)=(0.D0,0.D0)
    END DO

    CCOL=1.D0+CI*PZCL
    DO MM=0,NZMAX-1
       DZ=ZA(MM+1)-ZA(MM)
       AD=0.5D0/DZ
       BD=0.5D0/DZ
       DO I=MM,MM+1
          ID=3*I+2
          DO J=MM,MM+1
             JD=3*J+2
             IF(NWMAX.NE.NZMAX) JD=3*NWMAX+2+3*J-3*I
             DO KI=MM,MM+1
                YY=(1.D0+DBDZ*ZA(KI))/CCOL
                VY=YY/CCOL
                DS=DZ*D3(I-MM,J-MM,KI-MM)
                CP1=   DS*PN0*YY   /(1.D0-VY*VY)
                CP2=CI*DS*PN0*YY*VY/(1.D0-VY*VY)
                CP3=   DS*PN0*YY
                CPA=CONJG(CE(ID+1))*( CP1*CE(JD+1)+CP2*CE(JD+2)) &
                   +CONJG(CE(ID+2))*(-CP2*CE(JD+1)+CP1*CE(JD+2)) &
                   +CONJG(CE(ID+3))*  CP3*CE(JD+3)
                CPWR(MM)  =CPWR(MM)  +CI*AD*CPA
                CPWR(MM+1)=CPWR(MM+1)+CI*BD*CPA
                CPTOT     =CPTOT     +CI   *CPA 
             END DO
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE SUBPOW0

!     *****  CALCULATE LOCAL WAVE NUMBER  ***** 

  SUBROUTINE SUBKEI(NZMAX)

    USE wimcomm,ONLY: rkind,ZA,AKD1,AKD2, &
                      PN0,DBDZ,ANX,BETA,ZMIN,ZMAX,DZMAX,DZWID 
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZMAX
    INTEGER:: I
    REAL(rkind):: BETA2,Z,PN,WQQ,S,D,P,AKD3,AKD4,AKD5,AK1,AK2,TEMP

    BETA2=BETA*BETA 

    DO I=0,NZMAX 
       Z=ZA(I)
       PN=PN0*(1.D0+DBDZ*Z) 
       WQQ=1.D0+DBDZ*Z
       IF(ABS(WQQ*WQQ-1.D0).LT.1.D-8) WQQ=1.D0+1.D-8
       S=1.D0-PN/(1.D0-WQQ**2)
       D=PN*WQQ/(1.D0-WQQ**2)
       P=1.D0-PN
       IF(ABS(P).LT.1.D-8) P=1.D-8
       AKD3=-((S+P)*ANX*ANX-2*S*P)/(2*P) 
       AKD5= ((S-P)*(S-P)*ANX**4-4*D*D*P*ANX**2+4*D*D*P*P) 
       AKD4=SQRT(AKD5)/(2*P)
       AK1=BETA2*(AKD3+AKD4) 
       AK2=BETA2*(AKD3-AKD4) 
       IF(ABS(AK2).GT.ABS(AK1)) THEN
          TEMP=AK1
          AK1=AK2
          AK2=TEMP
       END IF
       AKD1(I)=AK1
       AKD2(I)=AK2
       IF(AKD1(I).LT.0.D0) THEN
          AKD1(I)=-SQRT(-AKD1(I))
       ELSE
          AKD1(I)= SQRT( AKD1(I))
       ENDIF
       IF(AKD2(I).LT.0.D0) THEN
          AKD2(I)=-SQRT(-AKD2(I))
       ELSE
          AKD2(I)= SQRT( AKD2(I))
       ENDIF
    END DO
    RETURN
  END SUBROUTINE SUBKEI

!     *****  EULER TRANSFOMATION  *****

  SUBROUTINE EUL(Z,CS,M,L)

    USE wimcomm,ONLY: rkind
    USE libde
    USE decomm,ONLY: G1,G2,G3,N1
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: Z
    COMPLEX(rkind),INTENT(OUT):: CS
    INTEGER,INTENT(IN):: M
    INTEGER,INTENT(OUT):: L
    REAL(rkind):: A(200),B(200)
    REAL(rkind),PARAMETER:: HP=1.5707963267948966192D0
    INTEGER:: ILST,K
    REAL(rkind):: H0,EPS,SR1,SI1,SR,SI,ESR,ESI,SR2,SI2,PARITY,SKR,SKI

    G2=HP
    G3=Z
    H0=0.5D0
    EPS=1.D-5
    ILST=0
     
    SR1=0.D0
    SI1=0.D0
    IF(M.NE.0) THEN 
       DO K=M-1,0,-1
          G1=DBLE(K)
          CALL DEFT(SR,ESR,H0,EPS,ILST,FUNR,'FUNR')
          CALL DEFT(SI,ESI,H0,EPS,ILST,FUNI,'FUNI')
          SR1=SR1+SR
          SI1=SI1+SI
       END DO
    ENDIF

    G1=DBLE(M)
    CALL DEFT(SR,ESR,H0,EPS,ILST,FUNR,'FUNR')
    CALL DEFT(SI,ESI,H0,EPS,ILST,FUNI,'FUNI')
    A(1)=SR
    SR2=0.5D0*SR
    B(1)=SI
    SI2=0.5D0*SI
    PARITY=-1.D0
    L=0

30  CONTINUE
    L=L+1
    IF(L.GE.200) GOTO 9000
    G1=DBLE(M+L)
    CALL DEFT(SR,ESR,H0,EPS,ILST,FUNR,'FUNR')
    CALL DEFT(SI,ESI,H0,EPS,ILST,FUNI,'FUNI')
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
    IF(DABS(SKR).GT.EPS.OR.DABS(SKI).GT.EPS) GOTO 30

    SR=(SR1+SR2)/DSQRT(4.D0*G2)
    SI=(SI1+SI2)/DSQRT(4.D0*G2)
    CS=CMPLX(SR,SI)
    RETURN

9000 WRITE(6,*) '## DIMENSION OVERFLOW IN EULER TRANSFORMATION.'
    RETURN
  END SUBROUTINE EUL

!     *****  REAL PART  *****

  FUNCTION FUNR(X,XM,XP)

    USE wimcomm,ONLY: rkind
    USE decomm,ONLY: G1,G2,G3,N1
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,XM,XP
    REAL(rkind):: FUNR,Y1,T,T2,YY,AN2

    Y1=XM
    IF(IDINT(G1).EQ.0) THEN 
       T=0.5D0*G2*XP
       T2=T*T
       YY=-0.5D0*G3*G3/T2
       IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
          IF(N1.EQ.1) THEN
             AN2=1.D0
          ELSEIF(N1.EQ.2) THEN
             AN2=T
          ELSEIF(N1.EQ.3) THEN
             AN2=T2
          ELSEIF(N1.EQ.4) THEN
             AN2=T2*T
          ELSEIF(N1.EQ.5) THEN
             AN2=G3*G3/(T2*T) 
          ENDIF
          FUNR=AN2*0.5*G2*DEXP(YY)*DCOS(T)
       ELSE
          FUNR=0.D0
       ENDIF
    ELSE
       T=G2*(X+2.D0*G1)
       T2=T*T
       YY=-0.5D0*G3*G3/T2
       IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
          IF(N1.EQ.1) THEN
             AN2=1.D0
          ELSEIF(N1.EQ.2) THEN
             AN2=T
          ELSEIF(N1.EQ.3) THEN
             AN2=T2
          ELSEIF(N1.EQ.4) THEN
             AN2=T2*T
          ELSEIF(N1.EQ.5) THEN
             AN2=G3*G3/(T2*T) 
          ENDIF
          FUNR=AN2*G2*DEXP(YY)*DCOS(T)
       ELSE
          FUNR=0.D0
       ENDIF
    ENDIF
    RETURN
  END FUNCTION FUNR

!     *****  IMAG PART  *****

  FUNCTION FUNI(X,XM,XP)

    USE wimcomm,ONLY: rkind
    USE decomm,ONLY: G1,G2,G3,N1
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,XM,XP
    REAL(rkind):: FUNI,Y1,Y2,T,T2,YY,AN2

    Y1=X
    Y2=XM
    T=G2*(XP+2.D0*G1)
    T2=T*T
    YY=-0.5D0*G3*G3/T2
    IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
       IF(N1.EQ.1) THEN
          AN2=1.D0
       ELSEIF(N1.EQ.2) THEN 
          AN2=T
       ELSEIF(N1.EQ.3) THEN
          AN2=T2
       ELSEIF(N1.EQ.4) THEN
          AN2=T2*T
       ELSEIF(N1.EQ.5) THEN
          AN2=G3*G3/(T2*T) 
       ENDIF
       FUNI=AN2*G2*DEXP(YY)*DSIN(T)
    ELSE
       FUNI=0.D0
    ENDIF
    RETURN
  END FUNCTION FUNI


!     *****  CALCULATE CUBIC SPLINE COEFFICIENT  ***** 

      SUBROUTINE DSPLC(X,N,Y,DF,IOPT,C,NC,IER)
      IMPLICIT REAL*8(A-H,O-Z)
 
!     COMPUTE THE COEFICIENTS OF THE CUBIC SPLINE.

      COMPLEX*16 Y(N),DF(2),C(NC,3),EC(4),D(2)
      COMPLEX*16 HY,Y1,Y2,PIV,DY1,DY2
      DIMENSION X(N),IOPT(2)
      IF(N.LT.2)    GOTO 9000
      IF(NC.LT.N-1) GOTO 9000
      IF(IOPT(1).LT.1 .OR. IOPT(1).GT.3) GOTO 9000
      IF(IOPT(2).LT.1 .OR. IOPT(2).GT.3) GOTO 9000
      NM1=N-1
      DO 5 I=1,NM1
         IF(X(I).GE.X(I+1)) GOTO 9100
    5 CONTINUE
      IER=0

!     SET THE END CONDITIONS.

      II=2
      KS=1
      KE=MIN0(4,N)
      IDER=1
      DO 70 I=1,2
         I1=2*I-1
         I2=2*I
         IB=IOPT(I)
         GOTO (10,20,30),IB
   10    EC(I1)=0.D0
         EC(I2)=2.D0*DF(I)
         GOTO 70
   20    D(I)=DF(I)
   25    IF(I.EQ.2) II=N
         H=X(II)-X(II-1)
         EC(I1)=1.D0
         HY=Y(II)-Y(II-1)
         EC(I2)=6.D0*(HY/H-D(I))/H
         IF(I.EQ.2) EC(I2)=-EC(I2)
         GOTO 70
   30    IF(I.EQ.1) GOTO 40
         KS=MAX0(1,N-3)
         KE=N
         IDER=N
   40    A2=0.D0
         D(I)=0.D0
         DO 60 K=KS,KE
            IF(IDER.EQ.K) GOTO 60
            A1=1.D0
            DO 50 J=KS,KE
               IF(J.EQ.IDER .OR. J.EQ.K) GOTO 50
               X1=X(IDER)-X(J)
               X2=X(K)-X(J)
               A1=A1*X1/X2
   50       CONTINUE
            X3=X(K)-X(IDER)
            D(I)=D(I)+A1*Y(K)/X3
            A2=A2-1.D0/X3
   60    CONTINUE
         D(I)=D(I)+Y(IDER)*A2
         GOTO 25
   70 CONTINUE

!     SET THE ELEMENTS FOR THE SYMMETRIC TRIDIAGONAL EQUATION. 

      IF(N.EQ.2) GOTO 90
      H1=X(2)-X(1)
      Y1=Y(2)-Y(1)
      DO 80 I=2,NM1
         H2=X(I+1)-X(I)
         Y2=Y(I+1)-Y(I)
         HH=H1+H2
         C(I,1)=H2/HH
         C(I,2)=1.D0-C(I,1)
         C(I,3)=6.D0*(Y2/H2-Y1/H1)/HH
         H1=H2
         Y1=Y2
   80 CONTINUE

!     SOLVE THE EQUATION.

   90 C(1,1)=-EC(1)*0.5
      C(1,2)= EC(2)*0.5
      IF(N.EQ.2) GOTO 110
      DO 100 K=2,NM1
         PIV=2.D0+C(K,2)*C(K-1,1)
         C(K,1)=-C(K,1)/PIV
         C(K,2)=(C(K,3)-C(K,2)*C(K-1,2))/PIV
  100 CONTINUE
  110 DY1=(EC(4)-EC(3)*C(NM1,2))/(2.D0+EC(3)*C(NM1,1))
      DO 120 I=1,NM1
         K=N-I
         DY2=C(K,1)*DY1+C(K,2)
         H=X(K+1)-X(K)
         C(K,3)=(DY1-DY2)/(6.D0*H)
         C(K,2)=0.5D0*DY2
         C(K,1)=(Y(K+1)-Y(K))/H-(C(K,2)+C(K,3)*H)*H
         DY1=DY2
  120 CONTINUE
         GOTO 9999

!         ERROR

 9000 IER=2
      WRITE(6,1000) N,NC,IOPT(1),IOPT(2)
 1000 FORMAT(1H0,'(SUBR.DSPLC)  N =',I4,', NC =',I4, &
                 ', IOPT(1) =',I3,', IOPT(2) =',I3 &
            /' ','N, NC, IOPT(1) ANDIOPT(2) SHOULD SATISFY' &
           /' ','THE FOLLOWING INEQUALITIES.' &
           /' ' ,'    2<=N, N-1<=NC, 1<=IOPT(1)<=3, 1<=IOPT(2)<=3 ' &
           /' ' ,'RETURN WITH NO CALCULATION.' )
      GOTO 9999
 9100 IER=1
      IP1=I+1
      WRITE(6,1100) I,X(I),IP1,X(IP1)
 1100 FORMAT(1H0,'(SUBR.DSPLC)  X(',I4,')=',F10.5, &
                             ', X(',I4,')=',F10.5 &
            /' ','X SHOULD SATISFY THE FOLLOWING INEQUALITIES.' &
            /' ','    X(1)<X(2)<...<X(N)' )
 9999 RETURN
   END SUBROUTINE DSPLC

!     *****  INTERPOLATE BY CUBIC SPLINE  ***** 

      SUBROUTINE DSPLF(X,N,Y,C,NC,V,M,F,IER)
      IMPLICIT REAL*8(A-H,O-Z)

!     INTERPOLATION BY THE CUBIC SPLINE.

      COMPLEX*16 Y(N),C(NC,3),F(M)
      DIMENSION X(N),V(M)
      IF(N.LT.2)    GOTO 9000
      IF(M.LT.1)    GOTO 9000
      IF(NC.LT.N-1) GOTO 9000
      IER=0
      I=1
      DO 90 K=1,M
         V1=V(K)-X(I)
         IF(V1) 10,30,40
   10    IF(I.GT.1) GOTO 20
         IER=1
         GOTO 80
   20    I=I-1
         V1=V(K)-X(I)
         IF(V1) 10,30,80
   30    F(K)=Y(I)
         GOTO 90
   40    IF(I.LT.N) GOTO 50
         IER=1
         I=N-1
         GOTO 80
   50    V2=V(K)-X(I+1)
         IF(V2) 80,60,70
   60    I=I+1
         GOTO 30
   70    I=I+1
         V1=V2
         GOTO 40
   80    F(K)=Y(I)+V1*(C(I,1)+V1*(C(I,2)+V1*C(I,3)))
   90 CONTINUE
      GOTO 9999
 9000 IER=2
      WRITE(6,1000) N,NC,M
 1000 FORMAT(' ','(SUBR.DSPLF)  N =',I4,', NC =',I4,', M =',I4/ &
             ' ','N, NC AND M SHOULD SATISFY THE FOLLOWING INEQUALITIES.'/ &
             ' ','    2<=N, N-1<=NC, 1<=M' &
             ' ','RETURN WITH NO CALCULATION.')
 9999 RETURN
  END SUBROUTINE DSPLF
END Module wimexec
