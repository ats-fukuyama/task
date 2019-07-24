MODULE libdsp

  PRIVATE
  PUBLIC DSPFNV,DSPFNVA,DSPFN,DSPFNA

CONTAINS

!     ***** PLASMA DISPERSION FUNCTION *****

  SUBROUTINE DSPFNV(X,Z,DZ,DDZ,DDDZ)
    IMPLICIT NONE
    COMPLEX(8),INTENT(IN):: X
    COMPLEX(8),INTENT(OUT):: Z,DZ,DDZ,DDDZ
    COMPLEX(8):: XA(1),ZA(1),DZA(1),DDZA(1),DDDZA(1)
    XA(1)=X
    CALL DSPFNVA(1,XA,ZA,DZA,DDZA,DDDZA)
    Z=ZA(1)
    DZ=DZA(1)
    DDZ=DDZA(1)
    DDDZ=DDDZA(1)
    RETURN
  END SUBROUTINE DSPFNV

!     ***** PLASMA DISPERSION FUNCTION *****

      SUBROUTINE DSPFNVA(N,X,Z,DZ,DDZ,DDDZ)

!      PROGRAMMED BY T. WATANABE (1991/01/09)
!      CODE IN KAKUYUUGOU-KENKYUU 

      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 X(*),Z(*),DZ(*),DDZ(*),DDDZ(*)
      SAVE INIT,H,CPAI1,CPAI2,COEF0,COEF1,COEF2,COEF3
      SAVE H2,HW,H2PAI,HPAI,HH
      SAVE EN11,EN12,EN13,EN14,EN15,EN16,EN17
      SAVE EN18,EN19,EN1A,EN1B,EN1C,EN1D,EN1E
      SAVE FN11,FN12,FN13,FN14,FN15,FN16,FN17
      SAVE FN18,FN19,FN1A,FN1B,FN1C,FN1D,FN1E
      SAVE EE1,EE2,EE3,EE4,EE5,EE6,EE7
      SAVE EE8,EE9,EEA,EEB,EEC,EED,EEE
      SAVE FE1,FE2,FE3,FE4,FE5,FE6,FE7
      SAVE FE8,FE9,FEA,FEB,FEC,FED,FEE
      SAVE EN21,EN22,EN23,EN24,EN25,EN26,EN27
      SAVE EN28,EN29,EN2A,EN2B,EN2C,EN2D,EN2E
      SAVE FN21,FN22,FN23,FN24,FN25,FN26,FN27
      SAVE FN28,FN29,FN2A,FN2B,FN2C,FN2D,FN2E
      SAVE EN1S,EN2S,FN1S,FN2S
      DATA INIT/0/

      IF(INIT.EQ.0) THEN
         INIT=1
         H=0.484375D0

         CPAI1= 3.141592653589793D0
         CPAI2= 6.283185307179586D0
         COEF0= 5.641895835477563D-1
         COEF1=-1.128379167095513D0
         COEF2=-COEF1
         COEF3=-2.256758334191025D0

         H2=H*0.5D0
         HW=H+H
         H2PAI=CPAI2/H
         HPAI=CPAI1/H
         HH=H*H

         EN11=          HH
         EN12=   4.00D0*HH
         EN13=   9.00D0*HH
         EN14=  16.00D0*HH
         EN15=  25.00D0*HH
         EN16=  36.00D0*HH
         EN17=  49.00D0*HH
         EN18=  64.00D0*HH
         EN19=  81.00D0*HH
         EN1A= 100.00D0*HH
         EN1B= 121.00D0*HH
         EN1C= 144.00D0*HH
         EN1D= 169.00D0*HH
         EN1E= 196.00D0*HH

         FN11=   0.25D0*HH
         FN12=   2.25D0*HH
         FN13=   6.25D0*HH
         FN14=  12.25D0*HH
         FN15=  20.25D0*HH
         FN16=  30.25D0*HH
         FN17=  42.25D0*HH
         FN18=  56.25D0*HH
         FN19=  72.25D0*HH
         FN1A=  90.25D0*HH
         FN1B= 110.25D0*HH
         FN1C= 132.25D0*HH
         FN1D= 156.25D0*HH
         FN1E= 182.25D0*HH

         EE1=EXP(-       HH)*HW
         EE2=EXP(- 3.0D0*HH)
         EE3=EXP(- 5.0D0*HH)
         EE4=EXP(- 7.0D0*HH)
         EE5=EXP(- 9.0D0*HH)
         EE6=EXP(-11.0D0*HH)
         EE7=EXP(-13.0D0*HH)
         EE8=EXP(-15.0D0*HH)
         EE9=EXP(-17.0D0*HH)
         EEA=EXP(-19.0D0*HH)
         EEB=EXP(-21.0D0*HH)
         EEC=EXP(-23.0D0*HH)
         EED=EXP(-25.0D0*HH)
         EEE=EXP(-27.0D0*HH)

         FE1=EXP(- .25D0*HH)*HW
         FE2=EXP(- 2.0D0*HH)
         FE3=EXP(- 4.0D0*HH)
         FE4=EXP(- 6.0D0*HH)
         FE5=EXP(- 8.0D0*HH)
         FE6=EXP(-10.0D0*HH)
         FE7=EXP(-12.0D0*HH)
         FE8=EXP(-14.0D0*HH)
         FE9=EXP(-16.0D0*HH)
         FEA=EXP(-18.0D0*HH)
         FEB=EXP(-20.0D0*HH)
         FEC=EXP(-22.0D0*HH)
         FED=EXP(-24.0D0*HH)
         FEE=EXP(-26.0D0*HH)

         EN21=EN11*EN11
         EN22=EN12*EN12
         EN23=EN13*EN13
         EN24=EN14*EN14
         EN25=EN15*EN15
         EN26=EN16*EN16
         EN27=EN17*EN17
         EN28=EN18*EN18
         EN29=EN19*EN19
         EN2A=EN1A*EN1A
         EN2B=EN1B*EN1B
         EN2C=EN1C*EN1C
         EN2D=EN1D*EN1D
         EN2E=EN1E*EN1E

         FN21=FN11*FN11
         FN22=FN12*FN12
         FN23=FN13*FN13
         FN24=FN14*FN14
         FN25=FN15*FN15
         FN26=FN16*FN16
         FN27=FN17*FN17
         FN28=FN18*FN18
         FN29=FN19*FN19
         FN2A=FN1A*FN1A
         FN2B=FN1B*FN1B
         FN2C=FN1C*FN1C
         FN2D=FN1D*FN1D
         FN2E=FN1E*FN1E

         EN1S=EN1E*EEE
         FN1S=FN1E*FEE
         EN2S=EN2E*EEE
         FN2S=FN2E*FEE
      ENDIF

      DO L=1,N
         XRE=DREAL( X(L) )
         XIM=DIMAG( X(L) )
         XXR=(XRE-XIM)*(XRE+XIM)
         XXI=XRE*XIM*2.0D0

         IF((XIM.LT.0.0D0).AND.(XXR.LT.-50.65D0)) THEN
            Z(L)   =1.D22
            DZ(L)  =1.D22
            DDZ(L) =1.D22
            DDDZ(L)=1.D22
            GOTO 1000
         ENDIF
         IF(ABS(XXI).GT.1.D-36) THEN
            SXI=XXI*XXI
         ELSE
            SXI=0.D0
         ENDIF
         IF(ABS(XIM).GT.1.D-36) THEN
            SIM=XIM*XIM*2.D0
         ELSE
            SIM=0.D0
         ENDIF
         IF(ABS(XRE).GT.1.D-36) THEN
            SRE=XRE*XRE*2.D0
         ELSE
            SRE=0.D0
         ENDIF

         ABZS2=1.0D0/(SRE+SIM)
         XREH2=ABS(XRE/H2)
         IF(XREH2.GT.1.D8) THEN
            IXRE=0
         ELSE
            IXRE=IDNINT(XRE/H2)
         ENDIF
         IF(MOD(IXRE,2).EQ.0) THEN
            D1=FN11-XXR
            D2=FN12-XXR
            D3=FN13-XXR
            D4=FN14-XXR
            D5=FN15-XXR
            D6=FN16-XXR
            D7=FN17-XXR
            D8=FN18-XXR
            D9=FN19-XXR
            DA=FN1A-XXR
            DB=FN1B-XXR
            DC=FN1C-XXR
            DD=FN1D-XXR
            DE=FN1E-XXR

            SD1=D1*D1
            SD2=D2*D2
            SD3=D3*D3
            SD4=D4*D4
            SD5=D5*D5
            SD6=D6*D6
            SD7=D7*D7
            SD8=D8*D8
            SD9=D9*D9
            SDA=DA*DA
            SDB=DB*DB
            SDC=DC*DC
            SDD=DD*DD
            SDE=DE*DE

            D11=1.0D0/(SD1+SXI)
            D12=1.0D0/(SD2+SXI)
            D13=1.0D0/(SD3+SXI)
            D14=1.0D0/(SD4+SXI)
            D15=1.0D0/(SD5+SXI)
            D16=1.0D0/(SD6+SXI)
            D17=1.0D0/(SD7+SXI)
            D18=1.0D0/(SD8+SXI)
            D19=1.0D0/(SD9+SXI)
            D1A=1.0D0/(SDA+SXI)
            D1B=1.0D0/(SDB+SXI)
            D1C=1.0D0/(SDC+SXI)
            D1D=1.0D0/(SDD+SXI)
            D1E=1.0D0/(SDE+SXI)

            D31=D11*D11*D11
            D32=D12*D12*D12
            D33=D13*D13*D13
            D34=D14*D14*D14
            D35=D15*D15*D15
            D36=D16*D16*D16
            D37=D17*D17*D17
            D38=D18*D18*D18
            D39=D19*D19*D19
            D3A=D1A*D1A*D1A
            D3B=D1B*D1B*D1B
            D3C=D1C*D1C*D1C
            D3D=D1D*D1D*D1D
            D3E=D1E*D1E*D1E

            DR1=D11*D1
            DR2=D12*D2
            DR3=D13*D3
            DR4=D14*D4
            DR5=D15*D5
            DR6=D16*D6
            DR7=D17*D7
            DR8=D18*D8
            DR9=D19*D9
            DRA=D1A*DA
            DRB=D1B*DB
            DRC=D1C*DC
            DRD=D1D*DD
            DRE=D1E*DE

            QR1=D31*(SD1-SXI*3.D0)*D1
            QR2=D32*(SD2-SXI*3.D0)*D2
            QR3=D33*(SD3-SXI*3.D0)*D3
            QR4=D34*(SD4-SXI*3.D0)*D4
            QR5=D35*(SD5-SXI*3.D0)*D5
            QR6=D36*(SD6-SXI*3.D0)*D6
            QR7=D37*(SD7-SXI*3.D0)*D7
            QR8=D38*(SD8-SXI*3.D0)*D8
            QR9=D39*(SD9-SXI*3.D0)*D9
            QRA=D3A*(SDA-SXI*3.D0)*DA
            QRB=D3B*(SDB-SXI*3.D0)*DB
            QRC=D3C*(SDC-SXI*3.D0)*DC
            QRD=D3D*(SDD-SXI*3.D0)*DD
            QRE=D3E*(SDE-SXI*3.D0)*DE

            QI1=D31*(SD1*3.D0-SXI)
            QI2=D32*(SD2*3.D0-SXI)
            QI3=D33*(SD3*3.D0-SXI)
            QI4=D34*(SD4*3.D0-SXI)
            QI5=D35*(SD5*3.D0-SXI)
            QI6=D36*(SD6*3.D0-SXI)
            QI7=D37*(SD7*3.D0-SXI)
            QI8=D38*(SD8*3.D0-SXI)
            QI9=D39*(SD9*3.D0-SXI)
            QIA=D3A*(SDA*3.D0-SXI)
            QIB=D3B*(SDB*3.D0-SXI)
            QIC=D3C*(SDC*3.D0-SXI)
            QID=D3D*(SDD*3.D0-SXI)
            QIE=D3E*(SDE*3.D0-SXI)

            SZ0R=(((((((((((((     DRE *FEE+     DRD)*FED &
              +     DRC)*FEC +     DRB)*FEB+     DRA)*FEA &
              +     DR9)*FE9 +     DR8)*FE8+     DR7)*FE7 &
              +     DR6)*FE6 +     DR5)*FE5+     DR4)*FE4 &
              +     DR3)*FE3 +     DR2)*FE2+     DR1)*FE1
            SZ0I=(((((((((((((     D1E *FEE+     D1D)*FED &
              +     D1C)*FEC +     D1B)*FEB+     D1A)*FEA &
              +     D19)*FE9 +     D18)*FE8+     D17)*FE7 &
              +     D16)*FE6 +     D15)*FE5+     D14)*FE4 &
              +     D13)*FE3 +     D12)*FE2+     D11)*FE1
            TZ1R=(((((((((((((FN1S*DRE     +FN1D*DRD)*FED &
              +FN1C*DRC)*FEC +FN1B*DRB)*FEB+FN1A*DRA)*FEA &
              +FN19*DR9)*FE9 +FN18*DR8)*FE8+FN17*DR7)*FE7 &
              +FN16*DR6)*FE6 +FN15*DR5)*FE5+FN14*DR4)*FE4 &
              +FN13*DR3)*FE3 +FN12*DR2)*FE2+FN11*DR1)*FE1
            TZ1I=(((((((((((((FN1S*D1E     +FN1D*D1D)*FED &
              +FN1C*D1C)*FEC +FN1B*D1B)*FEB+FN1A*D1A)*FEA &
              +FN19*D19)*FE9 +FN18*D18)*FE8+FN17*D17)*FE7 &
              +FN16*D16)*FE6 +FN15*D15)*FE5+FN14*D14)*FE4 &
              +FN13*D13)*FE3 +FN12*D12)*FE2+FN11*D11)*FE1*XXI
            RZ0R=(((((((((((((     QRE *FEE+     QRD)*FED &
              +     QRC)*FEC +     QRB)*FEB+     QRA)*FEA &
              +     QR9)*FE9 +     QR8)*FE8+     QR7)*FE7 &
              +     QR6)*FE6 +     QR5)*FE5+     QR4)*FE4 &
              +     QR3)*FE3 +     QR2)*FE2+     QR1)*FE1
            RZ0I=(((((((((((((     QIE *FEE+     QID)*FED &
              +     QIC)*FEC +     QIB)*FEB+     QIA)*FEA &
              +     QI9)*FE9 +     QI8)*FE8+     QI7)*FE7 &
              +     QI6)*FE6 +     QI5)*FE5+     QI4)*FE4 &
              +     QI3)*FE3 +     QI2)*FE2+     QI1)*FE1
            RZ1R=(((((((((((((FN1S*QRE     +FN1D*QRD)*FED &
              +FN1C*QRC)*FEC +FN1B*QRB)*FEB+FN1A*QRA)*FEA &
              +FN19*QR9)*FE9 +FN18*QR8)*FE8+FN17*QR7)*FE7 &
              +FN16*QR6)*FE6 +FN15*QR5)*FE5+FN14*QR4)*FE4 &
              +FN13*QR3)*FE3 +FN12*QR2)*FE2+FN11*QR1)*FE1
            RZ1I=(((((((((((((FN1S*QIE     +FN1D*QID)*FED &
              +FN1C*QIC)*FEC +FN1B*QIB)*FEB+FN1A*QIA)*FEA &
              +FN19*QI9)*FE9 +FN18*QI8)*FE8+FN17*QI7)*FE7 &
              +FN16*QI6)*FE6 +FN15*QI5)*FE5+FN14*QI4)*FE4 &
              +FN13*QI3)*FE3 +FN12*QI2)*FE2+FN11*QI1)*FE1
            RZ2R=(((((((((((((FN2S*QRE     +FN2D*QRD)*FED &
              +FN2C*QRC)*FEC +FN2B*QRB)*FEB+FN2A*QRA)*FEA &
              +FN29*QR9)*FE9 +FN28*QR8)*FE8+FN27*QR7)*FE7 &
              +FN26*QR6)*FE6 +FN25*QR5)*FE5+FN24*QR4)*FE4 &
              +FN23*QR3)*FE3 +FN22*QR2)*FE2+FN21*QR1)*FE1
            RZ2I=(((((((((((((FN2S*QIE     +FN2D*QID)*FED &
              +FN2C*QIC)*FEC +FN2B*QIB)*FEB+FN2A*QIA)*FEA &
              +FN29*QI9)*FE9 +FN28*QI8)*FE8+FN27*QI7)*FE7 &
              +FN26*QI6)*FE6 +FN25*QI5)*FE5+FN24*QI4)*FE4 &
              +FN23*QI3)*FE3 +FN22*QI2)*FE2+FN21*QI1)*FE1
            TZ0R=(SZ0R-SIM*SZ0I)*XRE
            TZ0I=(SZ0R+SRE*SZ0I)*XIM
            SZ2R=RZ1R*3.0D0+XXR*RZ0R-RZ0I*SXI
            SZ2I=RZ1I*3.0D0+    RZ0R+RZ0I*XXR
            TZ2R=(SZ2R-SIM*SZ2I)*XRE
            TZ2I=(SZ2R+SRE*SZ2I)*XIM
            TZ3R= RZ2R+(XXR*RZ1R-RZ1I*SXI)*3.0D0
            TZ3I=(RZ2I+(    RZ1R+RZ1I*XXR)*3.0D0)*XXI

            IF((XIM.GT.HPAI).OR.(XXR.GT.161.18D0)) THEN
               Z(L)   =DCMPLX(TZ0R,TZ0I)*COEF0
               DZ(L)  =DCMPLX(TZ1R,TZ1I)*COEF1
               DDZ(L) =DCMPLX(TZ2R,TZ2I)*COEF2
               DDDZ(L)=DCMPLX(TZ3R,TZ3I)*COEF3
               GOTO 1000
            ENDIF

            IF(ABS(XIM).LT.HPAI) THEN
               EXXPHR=DEXP(XIM*H2PAI)
               DC=COS(XRE*H2PAI)
               DS=SIN(XRE*H2PAI)
               CPD=CPAI2/((EXXPHR+DC*2.0D0)*EXXPHR+1.0D0)
               DP3A=CPD*EXXPHR/H
               DP3R=DP3A*DS
               DP3I=DP3A*(DC+EXXPHR)
               DI=(XIM-HPAI)*XRE*2.D0
               EXXR=EXP(-XXR)*CPD
               XN0R=(SIN(XXI)+SIN(DI)*EXXPHR)*EXXR
               XN0I=(COS(XXI)+COS(DI)*EXXPHR)*EXXR
               XN1R=XRE*XN0R-XIM*XN0I
               XN1I=XRE*XN0I+XIM*XN0R
               C1R=DP3R-XRE*2.0D0
               C1I=DP3I-XIM*2.0D0-HPAI
               HS2R=C1R*DP3R-C1I*DP3I+XXR*2.0D0-1.0D0
               HS2I=C1I*DP3R+C1R*DP3I+XXI*2.0D0
               XN2R=HS2R*XN0R-HS2I*XN0I
               XN2I=HS2R*XN0I+HS2I*XN0R
               HS3R=XRE*HS2R-XIM*HS2I+DP3R-XRE*2.0D0
               HS3I=XIM*HS2R+XRE*HS2I+DP3I-XIM*2.0D0
               XN3R=HS3R*XN0R-HS3I*XN0I
               XN3I=HS3R*XN0I+HS3I*XN0R
               Z(L)   = DCMPLX(TZ0R+XN0R,TZ0I+XN0I)*COEF0
               DZ(L)  = DCMPLX(TZ1R+XN1R,TZ1I+XN1I)*COEF1
               DDZ(L) = DCMPLX(TZ2R+XN2R,TZ2I+XN2I)*COEF2
               DDDZ(L)= DCMPLX(TZ3R+XN3R,TZ3I+XN3I)*COEF3
               GOTO 1000
            ELSE
               EXXR=EXP(-XXR)*CPAI2
               XN0R=SIN(XXI)*EXXR
               XN0I=COS(XXI)*EXXR
               XN1R=XRE*XN0R-XIM*XN0I
               XN1I=XRE*XN0I+XIM*XN0R
               X2R=XXR*2.0D0-1.0D0
               X2I=XXI*2.0D0
               XN2R=X2R*XN0R-X2I*XN0I
               XN2I=X2R*XN0I+X2I*XN0R
               X3R=(SRE-SIM*3.0D0-3.0D0)*XRE
               X3I=(SRE*3.0D0-SIM-3.0D0)*XIM
               XN3R=X3R*XN0R-X3I*XN0I
               XN3I=X3R*XN0I+X3I*XN0R
               Z(L)   = DCMPLX(TZ0R+XN0R,TZ0I+XN0I)*COEF0
               DZ(L)  = DCMPLX(TZ1R+XN1R,TZ1I+XN1I)*COEF1
               DDZ(L) = DCMPLX(TZ2R+XN2R,TZ2I+XN2I)*COEF2
               DDDZ(L)= DCMPLX(TZ3R+XN3R,TZ3I+XN3I)*COEF3
               GOTO 1000
            ENDIF
         ELSE
            D1=EN11-XXR
            D2=EN12-XXR
            D3=EN13-XXR
            D4=EN14-XXR
            D5=EN15-XXR
            D6=EN16-XXR
            D7=EN17-XXR
            D8=EN18-XXR
            D9=EN19-XXR
            DA=EN1A-XXR
            DB=EN1B-XXR
            DC=EN1C-XXR
            DD=EN1D-XXR
            DE=EN1E-XXR

            SD1=D1*D1
            SD2=D2*D2
            SD3=D3*D3
            SD4=D4*D4
            SD5=D5*D5
            SD6=D6*D6
            SD7=D7*D7
            SD8=D8*D8
            SD9=D9*D9
            SDA=DA*DA
            SDB=DB*DB
            SDC=DC*DC
            SDD=DD*DD
            SDE=DE*DE

            D11=1.0D0/(SD1+SXI)
            D12=1.0D0/(SD2+SXI)
            D13=1.0D0/(SD3+SXI)
            D14=1.0D0/(SD4+SXI)
            D15=1.0D0/(SD5+SXI)
            D16=1.0D0/(SD6+SXI)
            D17=1.0D0/(SD7+SXI)
            D18=1.0D0/(SD8+SXI)
            D19=1.0D0/(SD9+SXI)
            D1A=1.0D0/(SDA+SXI)
            D1B=1.0D0/(SDB+SXI)
            D1C=1.0D0/(SDC+SXI)
            D1D=1.0D0/(SDD+SXI)
            D1E=1.0D0/(SDE+SXI)

            D31=D11*D11*D11
            D32=D12*D12*D12
            D33=D13*D13*D13
            D34=D14*D14*D14
            D35=D15*D15*D15
            D36=D16*D16*D16
            D37=D17*D17*D17
            D38=D18*D18*D18
            D39=D19*D19*D19
            D3A=D1A*D1A*D1A
            D3B=D1B*D1B*D1B
            D3C=D1C*D1C*D1C
            D3D=D1D*D1D*D1D
            D3E=D1E*D1E*D1E

            DR1=D11*D1
            DR2=D12*D2
            DR3=D13*D3
            DR4=D14*D4
            DR5=D15*D5
            DR6=D16*D6
            DR7=D17*D7
            DR8=D18*D8
            DR9=D19*D9
            DRA=D1A*DA
            DRB=D1B*DB
            DRC=D1C*DC
            DRD=D1D*DD
            DRE=D1E*DE

            QR1=D31*(SD1-SXI*3.D0)*D1
            QR2=D32*(SD2-SXI*3.D0)*D2
            QR3=D33*(SD3-SXI*3.D0)*D3
            QR4=D34*(SD4-SXI*3.D0)*D4
            QR5=D35*(SD5-SXI*3.D0)*D5
            QR6=D36*(SD6-SXI*3.D0)*D6
            QR7=D37*(SD7-SXI*3.D0)*D7
            QR8=D38*(SD8-SXI*3.D0)*D8
            QR9=D39*(SD9-SXI*3.D0)*D9
            QRA=D3A*(SDA-SXI*3.D0)*DA
            QRB=D3B*(SDB-SXI*3.D0)*DB
            QRC=D3C*(SDC-SXI*3.D0)*DC
            QRD=D3D*(SDD-SXI*3.D0)*DD
            QRE=D3E*(SDE-SXI*3.D0)*DE

            QI1=D31*(SD1*3.D0-SXI)
            QI2=D32*(SD2*3.D0-SXI)
            QI3=D33*(SD3*3.D0-SXI)
            QI4=D34*(SD4*3.D0-SXI)
            QI5=D35*(SD5*3.D0-SXI)
            QI6=D36*(SD6*3.D0-SXI)
            QI7=D37*(SD7*3.D0-SXI)
            QI8=D38*(SD8*3.D0-SXI)
            QI9=D39*(SD9*3.D0-SXI)
            QIA=D3A*(SDA*3.D0-SXI)
            QIB=D3B*(SDB*3.D0-SXI)
            QIC=D3C*(SDC*3.D0-SXI)
            QID=D3D*(SDD*3.D0-SXI)
            QIE=D3E*(SDE*3.D0-SXI)

            SZ0R=(((((((((((((     DRE *EEE+     DRD)*EED &
              +     DRC)*EEC +     DRB)*EEB+     DRA)*EEA &
              +     DR9)*EE9 +     DR8)*EE8+     DR7)*EE7 &
              +     DR6)*EE6 +     DR5)*EE5+     DR4)*EE4 &
              +     DR3)*EE3 +     DR2)*EE2+     DR1)*EE1
            SZ0I=(((((((((((((     D1E *EEE+     D1D)*EED &
              +     D1C)*EEC +     D1B)*EEB+     D1A)*EEA &
              +     D19)*EE9 +     D18)*EE8+     D17)*EE7 &
              +     D16)*EE6 +     D15)*EE5+     D14)*EE4 &
              +     D13)*EE3 +     D12)*EE2+     D11)*EE1
            TZ1R=(((((((((((((EN1S*DRE     +EN1D*DRD)*EED &
              +EN1C*DRC)*EEC +EN1B*DRB)*EEB+EN1A*DRA)*EEA &
              +EN19*DR9)*EE9 +EN18*DR8)*EE8+EN17*DR7)*EE7 &
              +EN16*DR6)*EE6 +EN15*DR5)*EE5+EN14*DR4)*EE4 &
              +EN13*DR3)*EE3 +EN12*DR2)*EE2+EN11*DR1)*EE1
            TZ1I=(((((((((((((EN1S*D1E     +EN1D*D1D)*EED &
              +EN1C*D1C)*EEC +EN1B*D1B)*EEB+EN1A*D1A)*EEA &
              +EN19*D19)*EE9 +EN18*D18)*EE8+EN17*D17)*EE7 &
              +EN16*D16)*EE6 +EN15*D15)*EE5+EN14*D14)*EE4 &
              +EN13*D13)*EE3 +EN12*D12)*EE2+EN11*D11)*EE1*XXI
            RZ0R=(((((((((((((     QRE *EEE+     QRD)*EED &
              +     QRC)*EEC +     QRB)*EEB+     QRA)*EEA &
              +     QR9)*EE9 +     QR8)*EE8+     QR7)*EE7 &
              +     QR6)*EE6 +     QR5)*EE5+     QR4)*EE4 &
              +     QR3)*EE3 +     QR2)*EE2+     QR1)*EE1
            RZ0I=(((((((((((((     QIE *EEE+     QID)*EED &
              +     QIC)*EEC +     QIB)*EEB+     QIA)*EEA &
              +     QI9)*EE9 +     QI8)*EE8+     QI7)*EE7 &
              +     QI6)*EE6 +     QI5)*EE5+     QI4)*EE4 &
              +     QI3)*EE3 +     QI2)*EE2+     QI1)*EE1
            RZ1R=(((((((((((((EN1S*QRE     +EN1D*QRD)*EED &
              +EN1C*QRC)*EEC +EN1B*QRB)*EEB+EN1A*QRA)*EEA &
              +EN19*QR9)*EE9 +EN18*QR8)*EE8+EN17*QR7)*EE7 &
              +EN16*QR6)*EE6 +EN15*QR5)*EE5+EN14*QR4)*EE4 &
              +EN13*QR3)*EE3 +EN12*QR2)*EE2+EN11*QR1)*EE1
            RZ1I=(((((((((((((EN1S*QIE     +EN1D*QID)*EED &
              +EN1C*QIC)*EEC +EN1B*QIB)*EEB+EN1A*QIA)*EEA &
              +EN19*QI9)*EE9 +EN18*QI8)*EE8+EN17*QI7)*EE7 &
              +EN16*QI6)*EE6 +EN15*QI5)*EE5+EN14*QI4)*EE4 &
              +EN13*QI3)*EE3 +EN12*QI2)*EE2+EN11*QI1)*EE1
            RZ2R=(((((((((((((EN2S*QRE     +EN2D*QRD)*EED &
              +EN2C*QRC)*EEC +EN2B*QRB)*EEB+EN2A*QRA)*EEA &
              +EN29*QR9)*EE9 +EN28*QR8)*EE8+EN27*QR7)*EE7 &
              +EN26*QR6)*EE6 +EN25*QR5)*EE5+EN24*QR4)*EE4 &
              +EN23*QR3)*EE3 +EN22*QR2)*EE2+EN21*QR1)*EE1
            RZ2I=(((((((((((((EN2S*QIE     +EN2D*QID)*EED &
              +EN2C*QIC)*EEC +EN2B*QIB)*EEB+EN2A*QIA)*EEA &
              +EN29*QI9)*EE9 +EN28*QI8)*EE8+EN27*QI7)*EE7 &
              +EN26*QI6)*EE6 +EN25*QI5)*EE5+EN24*QI4)*EE4 &
              +EN23*QI3)*EE3 +EN22*QI2)*EE2+EN21*QI1)*EE1
            DD=HW*ABZS2
            DQ=DD*ABZS2*ABZS2*2.0D0
            TZ0R=(SZ0R-SIM*SZ0I-DD)*XRE
            TZ0I=(SZ0R+SRE*SZ0I+DD)*XIM
            SZ2R=RZ1R*3.0D0+XXR*RZ0R-RZ0I*SXI
            SZ2I=RZ1I*3.0D0+    RZ0R+RZ0I*XXR
            TZ2R=(SZ2R-SIM*SZ2I-(SRE      -SIM*3.0D0)*DQ)*XRE
            TZ2I=(SZ2R+SRE*SZ2I+(SRE*3.0D0-SIM      )*DQ)*XIM
            TZ3R= RZ2R+(XXR*RZ1R-RZ1I*SXI)*3.0D0
            TZ3I=(RZ2I+(    RZ1R+RZ1I*XXR)*3.0D0)*XXI

            IF((XIM.GT.HPAI).OR.(XXR.GT.161.18D0)) THEN
               Z(L)   =DCMPLX(TZ0R,TZ0I)*COEF0
               DZ(L)  =DCMPLX(TZ1R,TZ1I)*COEF1
               DDZ(L) =DCMPLX(TZ2R,TZ2I)*COEF2
               DDDZ(L)=DCMPLX(TZ3R,TZ3I)*COEF3
               GOTO 1000
            ENDIF

            IF(ABS(XIM).LT.HPAI) THEN
               EXXPHR=DEXP(XIM*H2PAI)
               DC=COS(XRE*H2PAI)
               DS=SIN(XRE*H2PAI)
               CPD=CPAI2/((EXXPHR-DC*2.0D0)*EXXPHR+1.0D0)
               DP3A=CPD*EXXPHR/H
               DP3R=DP3A*DS
               DP3I=DP3A*(DC-EXXPHR)
               DI=(XIM-HPAI)*XRE*2.D0
               EXXR=EXP(-XXR)*CPD
               XN0R=(SIN(XXI)-SIN(DI)*EXXPHR)*EXXR
               XN0I=(COS(XXI)-COS(DI)*EXXPHR)*EXXR
               XN1R=XRE*XN0R-XIM*XN0I
               XN1I=XRE*XN0I+XIM*XN0R
               C1R=DP3R+XRE*2.0D0
               C1I=DP3I+XIM*2.0D0+HPAI
               HS2R=C1R*DP3R-C1I*DP3I+XXR*2.0D0-1.0D0
               HS2I=C1I*DP3R+C1R*DP3I+XXI*2.0D0
               XN2R=HS2R*XN0R-HS2I*XN0I
               XN2I=HS2R*XN0I+HS2I*XN0R
               HS3R=XRE*HS2R-XIM*HS2I-DP3R-XRE*2.0D0
               HS3I=XIM*HS2R+XRE*HS2I-DP3I-XIM*2.0D0
               XN3R=HS3R*XN0R-HS3I*XN0I
               XN3I=HS3R*XN0I+HS3I*XN0R
               Z(L)   = DCMPLX(TZ0R+XN0R,TZ0I+XN0I)*COEF0
               DZ(L)  = DCMPLX(TZ1R+XN1R,TZ1I+XN1I)*COEF1
               DDZ(L) = DCMPLX(TZ2R+XN2R,TZ2I+XN2I)*COEF2
               DDDZ(L)= DCMPLX(TZ3R+XN3R,TZ3I+XN3I)*COEF3
               GOTO 1000
            ELSE
               EXXR=EXP(-XXR)*CPAI2
               XN0R=SIN(XXI)*EXXR
               XN0I=COS(XXI)*EXXR
               XN1R=XRE*XN0R-XIM*XN0I
               XN1I=XRE*XN0I+XIM*XN0R
               X2R=XXR*2.0D0-1.0D0
               X2I=XXI*2.0D0
               XN2R=X2R*XN0R-X2I*XN0I
               XN2I=X2R*XN0I+X2I*XN0R
               X3R=(SRE-SIM*3.0D0-3.0D0)*XRE
               X3I=(SRE*3.0D0-SIM-3.0D0)*XIM
               XN3R=X3R*XN0R-X3I*XN0I
               XN3I=X3R*XN0I+X3I*XN0R
               Z(L)   = DCMPLX(TZ0R+XN0R,TZ0I+XN0I)*COEF0
               DZ(L)  = DCMPLX(TZ1R+XN1R,TZ1I+XN1I)*COEF1
               DDZ(L) = DCMPLX(TZ2R+XN2R,TZ2I+XN2I)*COEF2
               DDDZ(L)= DCMPLX(TZ3R+XN3R,TZ3I+XN3I)*COEF3
               GOTO 1000
            ENDIF
         ENDIF
 1000    CONTINUE
      ENDDO
      RETURN
  END SUBROUTINE DSPFNVA

!     ***** PLASMA DISPERSION FUNCTION *****

  SUBROUTINE DSPFN(X,Z,DZ)
    IMPLICIT NONE
    COMPLEX(8),INTENT(IN):: X
    COMPLEX(8),INTENT(OUT):: Z,DZ
    COMPLEX(8):: XA(1),ZA(1),DZA(1)
    XA(1)=X
    CALL DSPFNA(1,XA,ZA,DZA)
    Z=ZA(1)
    DZ=DZA(1)
    RETURN
  END SUBROUTINE DSPFN

  SUBROUTINE DSPFNA(N,X,Z,DZ)

!      PROGRAMMED BY T. WATANABE (1991/01/09)
!      CODE IN KAKUYUUGOU-KENKYUU 

      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 X(*),Z(*),DZ(*)
      SAVE INIT,H,CPAI1,CPAI2,COEF0,COEF1
      SAVE H2,HW,H2PAI,HPAI,HH
      SAVE EN11,EN12,EN13,EN14,EN15,EN16,EN17
      SAVE EN18,EN19,EN1A,EN1B,EN1C,EN1D,EN1E
      SAVE FN11,FN12,FN13,FN14,FN15,FN16,FN17
      SAVE FN18,FN19,FN1A,FN1B,FN1C,FN1D,FN1E
      SAVE EE1,EE2,EE3,EE4,EE5,EE6,EE7
      SAVE EE8,EE9,EEA,EEB,EEC,EED,EEE
      SAVE FE1,FE2,FE3,FE4,FE5,FE6,FE7
      SAVE FE8,FE9,FEA,FEB,FEC,FED,FEE
      SAVE EN1S,FN1S
      DATA INIT/0/

      IF(INIT.EQ.0) THEN
         INIT=1
         H=0.484375D0

         CPAI1= 3.141592653589793D0
         CPAI2= 6.283185307179586D0
         COEF0= 5.641895835477563D-1
         COEF1=-1.128379167095513D0

         H2=H*0.5D0
         HW=H+H
         H2PAI=CPAI2/H
         HPAI=CPAI1/H
         HH=H*H

         EN11=          HH
         EN12=   4.00D0*HH
         EN13=   9.00D0*HH
         EN14=  16.00D0*HH
         EN15=  25.00D0*HH
         EN16=  36.00D0*HH
         EN17=  49.00D0*HH
         EN18=  64.00D0*HH
         EN19=  81.00D0*HH
         EN1A= 100.00D0*HH
         EN1B= 121.00D0*HH
         EN1C= 144.00D0*HH
         EN1D= 169.00D0*HH
         EN1E= 196.00D0*HH

         FN11=   0.25D0*HH
         FN12=   2.25D0*HH
         FN13=   6.25D0*HH
         FN14=  12.25D0*HH
         FN15=  20.25D0*HH
         FN16=  30.25D0*HH
         FN17=  42.25D0*HH
         FN18=  56.25D0*HH
         FN19=  72.25D0*HH
         FN1A=  90.25D0*HH
         FN1B= 110.25D0*HH
         FN1C= 132.25D0*HH
         FN1D= 156.25D0*HH
         FN1E= 182.25D0*HH

         EE1=EXP(-       HH)*HW
         EE2=EXP(- 3.0D0*HH)
         EE3=EXP(- 5.0D0*HH)
         EE4=EXP(- 7.0D0*HH)
         EE5=EXP(- 9.0D0*HH)
         EE6=EXP(-11.0D0*HH)
         EE7=EXP(-13.0D0*HH)
         EE8=EXP(-15.0D0*HH)
         EE9=EXP(-17.0D0*HH)
         EEA=EXP(-19.0D0*HH)
         EEB=EXP(-21.0D0*HH)
         EEC=EXP(-23.0D0*HH)
         EED=EXP(-25.0D0*HH)
         EEE=EXP(-27.0D0*HH)

         FE1=EXP(- .25D0*HH)*HW
         FE2=EXP(- 2.0D0*HH)
         FE3=EXP(- 4.0D0*HH)
         FE4=EXP(- 6.0D0*HH)
         FE5=EXP(- 8.0D0*HH)
         FE6=EXP(-10.0D0*HH)
         FE7=EXP(-12.0D0*HH)
         FE8=EXP(-14.0D0*HH)
         FE9=EXP(-16.0D0*HH)
         FEA=EXP(-18.0D0*HH)
         FEB=EXP(-20.0D0*HH)
         FEC=EXP(-22.0D0*HH)
         FED=EXP(-24.0D0*HH)
         FEE=EXP(-26.0D0*HH)

         EN1S=EN1E*EEE
         FN1S=FN1E*FEE
      ENDIF

      DO L=1,N
         XRE=DREAL( X(L) )
         XIM=DIMAG( X(L) )
         XXR=(XRE-XIM)*(XRE+XIM)
         XXI=XRE*XIM*2.0D0

         IF((XIM.LT.0.0D0).AND.(XXR.LT.-50.65D0)) THEN
            Z(L)   =1.D22
            DZ(L)  =1.D22
            GOTO 1000
         ENDIF
         IF(ABS(XXI).GT.1.D-36) THEN
            SXI=XXI*XXI
         ELSE
            SXI=0.D0
         ENDIF
         IF(ABS(XIM).GT.1.D-36) THEN
            SIM=XIM*XIM*2.D0
         ELSE
            SIM=0.D0
         ENDIF
         IF(ABS(XRE).GT.1.D-36) THEN
            SRE=XRE*XRE*2.D0
         ELSE
            SRE=0.D0
         ENDIF

         ABZS2=1.0D0/(SRE+SIM)
         XREH2=ABS(XRE/H2)
         IF(XREH2.GT.1.D8) THEN
            IXRE=0
         ELSE
            IXRE=IDNINT(XRE/H2)
         ENDIF
         IF(MOD(IXRE,2).EQ.0) THEN
            D1=FN11-XXR
            D2=FN12-XXR
            D3=FN13-XXR
            D4=FN14-XXR
            D5=FN15-XXR
            D6=FN16-XXR
            D7=FN17-XXR
            D8=FN18-XXR
            D9=FN19-XXR
            DA=FN1A-XXR
            DB=FN1B-XXR
            DC=FN1C-XXR
            DD=FN1D-XXR
            DE=FN1E-XXR

            SD1=D1*D1
            SD2=D2*D2
            SD3=D3*D3
            SD4=D4*D4
            SD5=D5*D5
            SD6=D6*D6
            SD7=D7*D7
            SD8=D8*D8
            SD9=D9*D9
            SDA=DA*DA
            SDB=DB*DB
            SDC=DC*DC
            SDD=DD*DD
            SDE=DE*DE

            D11=1.0D0/(SD1+SXI)
            D12=1.0D0/(SD2+SXI)
            D13=1.0D0/(SD3+SXI)
            D14=1.0D0/(SD4+SXI)
            D15=1.0D0/(SD5+SXI)
            D16=1.0D0/(SD6+SXI)
            D17=1.0D0/(SD7+SXI)
            D18=1.0D0/(SD8+SXI)
            D19=1.0D0/(SD9+SXI)
            D1A=1.0D0/(SDA+SXI)
            D1B=1.0D0/(SDB+SXI)
            D1C=1.0D0/(SDC+SXI)
            D1D=1.0D0/(SDD+SXI)
            D1E=1.0D0/(SDE+SXI)

            DR1=D11*D1
            DR2=D12*D2
            DR3=D13*D3
            DR4=D14*D4
            DR5=D15*D5
            DR6=D16*D6
            DR7=D17*D7
            DR8=D18*D8
            DR9=D19*D9
            DRA=D1A*DA
            DRB=D1B*DB
            DRC=D1C*DC
            DRD=D1D*DD
            DRE=D1E*DE

            SZ0R=(((((((((((((     DRE *FEE+     DRD)*FED &
              +     DRC)*FEC +     DRB)*FEB+     DRA)*FEA &
              +     DR9)*FE9 +     DR8)*FE8+     DR7)*FE7 &
              +     DR6)*FE6 +     DR5)*FE5+     DR4)*FE4 &
              +     DR3)*FE3 +     DR2)*FE2+     DR1)*FE1
            SZ0I=(((((((((((((     D1E *FEE+     D1D)*FED &
              +     D1C)*FEC +     D1B)*FEB+     D1A)*FEA &
              +     D19)*FE9 +     D18)*FE8+     D17)*FE7 &
              +     D16)*FE6 +     D15)*FE5+     D14)*FE4 &
              +     D13)*FE3 +     D12)*FE2+     D11)*FE1
            TZ1R=(((((((((((((FN1S*DRE     +FN1D*DRD)*FED &
              +FN1C*DRC)*FEC +FN1B*DRB)*FEB+FN1A*DRA)*FEA &
              +FN19*DR9)*FE9 +FN18*DR8)*FE8+FN17*DR7)*FE7 &
              +FN16*DR6)*FE6 +FN15*DR5)*FE5+FN14*DR4)*FE4 &
              +FN13*DR3)*FE3 +FN12*DR2)*FE2+FN11*DR1)*FE1
            TZ1I=(((((((((((((FN1S*D1E     +FN1D*D1D)*FED &
              +FN1C*D1C)*FEC +FN1B*D1B)*FEB+FN1A*D1A)*FEA &
              +FN19*D19)*FE9 +FN18*D18)*FE8+FN17*D17)*FE7 &
              +FN16*D16)*FE6 +FN15*D15)*FE5+FN14*D14)*FE4 &
              +FN13*D13)*FE3 +FN12*D12)*FE2+FN11*D11)*FE1*XXI
            TZ0R=(SZ0R-SIM*SZ0I)*XRE
            TZ0I=(SZ0R+SRE*SZ0I)*XIM

            IF((XIM.GT.HPAI).OR.(XXR.GT.161.18D0)) THEN
               Z(L)   =DCMPLX(TZ0R,TZ0I)*COEF0
               DZ(L)  =DCMPLX(TZ1R,TZ1I)*COEF1
               GOTO 1000
            ENDIF

            IF(ABS(XIM).LT.HPAI) THEN
               EXXPHR=DEXP(XIM*H2PAI)
               DC=COS(XRE*H2PAI)
               CPD=CPAI2/((EXXPHR+DC*2.0D0)*EXXPHR+1.0D0)
               DI=(XIM-HPAI)*XRE*2.D0
               EXXR=EXP(-XXR)*CPD
               XN0R=(SIN(XXI)+SIN(DI)*EXXPHR)*EXXR
               XN0I=(COS(XXI)+COS(DI)*EXXPHR)*EXXR
               XN1R=XRE*XN0R-XIM*XN0I
               XN1I=XRE*XN0I+XIM*XN0R
               Z(L)   = DCMPLX(TZ0R+XN0R,TZ0I+XN0I)*COEF0
               DZ(L)  = DCMPLX(TZ1R+XN1R,TZ1I+XN1I)*COEF1
               GOTO 1000
            ELSE
               EXXR=EXP(-XXR)*CPAI2
               XN0R=SIN(XXI)*EXXR
               XN0I=COS(XXI)*EXXR
               XN1R=XRE*XN0R-XIM*XN0I
               XN1I=XRE*XN0I+XIM*XN0R
               Z(L)   = DCMPLX(TZ0R+XN0R,TZ0I+XN0I)*COEF0
               DZ(L)  = DCMPLX(TZ1R+XN1R,TZ1I+XN1I)*COEF1
               GOTO 1000
            ENDIF
         ELSE
            D1=EN11-XXR
            D2=EN12-XXR
            D3=EN13-XXR
            D4=EN14-XXR
            D5=EN15-XXR
            D6=EN16-XXR
            D7=EN17-XXR
            D8=EN18-XXR
            D9=EN19-XXR
            DA=EN1A-XXR
            DB=EN1B-XXR
            DC=EN1C-XXR
            DD=EN1D-XXR
            DE=EN1E-XXR

            SD1=D1*D1
            SD2=D2*D2
            SD3=D3*D3
            SD4=D4*D4
            SD5=D5*D5
            SD6=D6*D6
            SD7=D7*D7
            SD8=D8*D8
            SD9=D9*D9
            SDA=DA*DA
            SDB=DB*DB
            SDC=DC*DC
            SDD=DD*DD
            SDE=DE*DE

            D11=1.0D0/(SD1+SXI)
            D12=1.0D0/(SD2+SXI)
            D13=1.0D0/(SD3+SXI)
            D14=1.0D0/(SD4+SXI)
            D15=1.0D0/(SD5+SXI)
            D16=1.0D0/(SD6+SXI)
            D17=1.0D0/(SD7+SXI)
            D18=1.0D0/(SD8+SXI)
            D19=1.0D0/(SD9+SXI)
            D1A=1.0D0/(SDA+SXI)
            D1B=1.0D0/(SDB+SXI)
            D1C=1.0D0/(SDC+SXI)
            D1D=1.0D0/(SDD+SXI)
            D1E=1.0D0/(SDE+SXI)

            DR1=D11*D1
            DR2=D12*D2
            DR3=D13*D3
            DR4=D14*D4
            DR5=D15*D5
            DR6=D16*D6
            DR7=D17*D7
            DR8=D18*D8
            DR9=D19*D9
            DRA=D1A*DA
            DRB=D1B*DB
            DRC=D1C*DC
            DRD=D1D*DD
            DRE=D1E*DE

            SZ0R=(((((((((((((     DRE *EEE+     DRD)*EED &
              +     DRC)*EEC +     DRB)*EEB+     DRA)*EEA &
              +     DR9)*EE9 +     DR8)*EE8+     DR7)*EE7 &
              +     DR6)*EE6 +     DR5)*EE5+     DR4)*EE4 &
              +     DR3)*EE3 +     DR2)*EE2+     DR1)*EE1
            SZ0I=(((((((((((((     D1E *EEE+     D1D)*EED &
              +     D1C)*EEC +     D1B)*EEB+     D1A)*EEA &
              +     D19)*EE9 +     D18)*EE8+     D17)*EE7 &
              +     D16)*EE6 +     D15)*EE5+     D14)*EE4 &
              +     D13)*EE3 +     D12)*EE2+     D11)*EE1
            TZ1R=(((((((((((((EN1S*DRE     +EN1D*DRD)*EED &
              +EN1C*DRC)*EEC +EN1B*DRB)*EEB+EN1A*DRA)*EEA &
              +EN19*DR9)*EE9 +EN18*DR8)*EE8+EN17*DR7)*EE7 &
              +EN16*DR6)*EE6 +EN15*DR5)*EE5+EN14*DR4)*EE4 &
              +EN13*DR3)*EE3 +EN12*DR2)*EE2+EN11*DR1)*EE1
            TZ1I=(((((((((((((EN1S*D1E     +EN1D*D1D)*EED &
              +EN1C*D1C)*EEC +EN1B*D1B)*EEB+EN1A*D1A)*EEA &
              +EN19*D19)*EE9 +EN18*D18)*EE8+EN17*D17)*EE7 &
              +EN16*D16)*EE6 +EN15*D15)*EE5+EN14*D14)*EE4 &
              +EN13*D13)*EE3 +EN12*D12)*EE2+EN11*D11)*EE1*XXI
            DD=HW*ABZS2
            TZ0R=(SZ0R-SIM*SZ0I-DD)*XRE
            TZ0I=(SZ0R+SRE*SZ0I+DD)*XIM

            IF((XIM.GT.HPAI).OR.(XXR.GT.161.18D0)) THEN
               Z(L)   =DCMPLX(TZ0R,TZ0I)*COEF0
               DZ(L)  =DCMPLX(TZ1R,TZ1I)*COEF1
               GOTO 1000
            ENDIF

            IF(ABS(XIM).LT.HPAI) THEN
               EXXPHR=DEXP(XIM*H2PAI)
               DC=COS(XRE*H2PAI)
               CPD=CPAI2/((EXXPHR-DC*2.0D0)*EXXPHR+1.0D0)
               DI=(XIM-HPAI)*XRE*2.D0
               EXXR=EXP(-XXR)*CPD
               XN0R=(SIN(XXI)-SIN(DI)*EXXPHR)*EXXR
               XN0I=(COS(XXI)-COS(DI)*EXXPHR)*EXXR
               XN1R=XRE*XN0R-XIM*XN0I
               XN1I=XRE*XN0I+XIM*XN0R
               Z(L)   = DCMPLX(TZ0R+XN0R,TZ0I+XN0I)*COEF0
               DZ(L)  = DCMPLX(TZ1R+XN1R,TZ1I+XN1I)*COEF1
               GOTO 1000
            ELSE
               EXXR=EXP(-XXR)*CPAI2
               XN0R=SIN(XXI)*EXXR
               XN0I=COS(XXI)*EXXR
               XN1R=XRE*XN0R-XIM*XN0I
               XN1I=XRE*XN0I+XIM*XN0R
               Z(L)   = DCMPLX(TZ0R+XN0R,TZ0I+XN0I)*COEF0
               DZ(L)  = DCMPLX(TZ1R+XN1R,TZ1I+XN1I)*COEF1
               GOTO 1000
            ENDIF
         ENDIF
 1000    CONTINUE
      ENDDO
      RETURN
  END SUBROUTINE DSPFNA
END MODULE libdsp
