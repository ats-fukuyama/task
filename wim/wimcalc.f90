! wimcalc.f90

MODULE wimcalc
  PRIVATE
  PUBLIC wim_calc_coef0
  PUBLIC wim_calc_coef2
  PUBLIC wim_calc_power0
  PUBLIC wim_calc_power2

CONTAINS

!     *****  CALCULATE COEFFICIENT MATRIX  ***** 

  SUBROUTINE wim_calc_coef0(NZMAX,NWMAX,MODELW) 

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
  END SUBROUTINE wim_calc_coef0
      
!     *****  CFN1( V )  *****

  FUNCTION CFN1(V)

    USE wimcomm,ONLY: rkind,PI,CI,TT,TMAX,NTMAX,CFK1,CFKS1
    USE libpdkf
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
    USE libpdkf
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

  !*****  CALCULATE COEFFICIENT MATRIX  ***** 

  SUBROUTINE wim_calc_coef2(NZMAX,NWMAX,MODELW) 

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
  END SUBROUTINE wim_calc_coef2

!     *****  CALCULATE ABSORBED POWER  ***** 

  SUBROUTINE wim_calc_power0(NZMAX,NWMAX,CPTOT)

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
  END SUBROUTINE wim_calc_power0

!     *****  CALCULATE L!     *****  CALCULATE ABSORBED POWER  ***** 

  SUBROUTINE wim_calc_power2(NZMAX,NWMAX,CPTOT)

    USE wimcomm,ONLY: rkind,CI,ZA,CE,PN0,DBDZ,ANX,BETA,ZMIN,ZMAX,DZMAX,DZWID, &
                      CPWR,D0,D1,D2,D3
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZMAX,NWMAX
    COMPLEX(rkind),INTENT(OUT):: CPTOT
    INTEGER:: NZ,MM,NN,NS,NE,I,J,ID,JD,KI,KJ
    REAL(rkind):: DZI,DZJ,AD,BD,YY,YYI,YYJ,VP,VM,VPI,VPJ,VMI,VMJ, &
                  VVP,VVM,VVZ,X1,X4,X5,DPI,DMI,DPJ,DMJ,DPM,DP,DM,DE
    COMPLEX(rkind):: CF1P,CF1M,CF2Z,CF1Z,CP1,CP2,CP3,CPA
     
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
                      CPWR(MM)  =CPWR(MM)  +CI*AD*CPA
                      CPWR(MM+1)=CPWR(MM+1)+CI*BD*CPA
                      CPTOT     =CPTOT     +CI   *CPA 
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE wim_calc_power2
  
END MODULE wimcalc

