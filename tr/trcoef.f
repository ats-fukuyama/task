C     $Id$
C     ***********************************************************
C
C           CALCULATE TRANSPORT COEFFICIENTS
C 
C     ***********************************************************
C
      SUBROUTINE TRCOEF
C
      INCLUDE 'trcomm.h'
C
      CALL TRCFDW
      CALL TRCFNC
      CALL TRCFET
      CALL TRCFAD
C
      RETURN
      END
C
C     ***********************************************************
C
      SUBROUTINE TRCFDW
C
      INCLUDE 'trcomm.h'
C
      AMD=PA(2)*AMM
      AMT=PA(3)*AMM
      AMA=PA(4)*AMM
      OMEGAD = PZ(2)*AEE*BB/AMD
C
      IF(MDLKAI.LT.10) THEN
         KGR1='/TI  vs r/'
         KGR2='/DTI  vs r/'
         KGR3='/CLN,CLT  vs r/'
         KGR4='/CLS  vs r/'
      ELSEIF(MDLKAI.LT.20) THEN
         KGR1='/LOG(DEDW),LOG(DIDW)  vs r/'
         KGR2='/ETAI,ETAC  vs r/'
         KGR3='/CLN,CLT  vs r/'
         KGR4='/CLS  vs r/'
      ELSEIF(MDLKAI.LT.30) THEN
         KGR1='/DTE,CRTCL  vs r/'
         KGR2='/ETAI,ETAC  vs r/'
         KGR3='/CLN,CLT  vs r/'
         KGR4='/CLS  vs r/'
      ELSEIF(MDLKAI.LT.40) THEN
         KGR1='/NST^2 vs r/'
         KGR2='/OmegaST vs r/'
         KGR3='@lambda vs r@'
         KGR4='@Lambda,1/(1+OmgST^2),1/(1+G*WE1^2)vs r@'
      ELSE
         KGR1='//'
         KGR2='//'
         KGR3='//'
         KGR4='//'
      ENDIF
C    
C     ***** PP   is the total pressure (nT) *****
C     ***** DPP  is the derivative of total pressure (dnT/dr) *****
C     ***** TI   is the ion temperature (Ti) *****
C     ***** DTI  is the derivative of ion temperature (dTi/dr) *****
C     ***** TE   is the electron temperature (Te) *****
C     ***** DTE  is the derivative of electron temperature (dTe/dr) *****
C     ***** ANE  is the electron density (ne) *****
C     ***** DNE  is the derivative of electron density (dne/dr) *****
C     ***** CLN  is the density scale length ne/(dne/dr) *****
C     ***** CLT  is the ion temerature scale length Ti/(dTi/dr) *****
C     ***** DQ   is the derivative of safety factor (dq/dr) *****
C     ***** CLS  is the shear length R*q**2/(r*dq/dr) *****
C
      DO 100 NR=1,NRMAX
         IF(NR.EQ.NRMAX) THEN
            ANE=PNSS(1)
            AND=PNSS(2)
            ANT=PNSS(3)
            ANA=PNSS(4)
            TE=PTS(1)
            TD=PTS(2)
            TT=PTS(3)
            TA=PTS(4)
C
            RNTP= PNSS(2)*PTS(2)
     &           +PNSS(3)*PTS(3)
     &           +PNSS(4)*PTS(4)
            RNP = PNSS(2)+PNSS(3)+PNSS(4)
            RNTM= RN(NR-1,2)*RT(NR-1,2)
     &           +RN(NR-1,3)*RT(NR-1,3)
     &           +RN(NR-1,4)*RT(NR-1,4)
     &           +RW(NR-1,1)+RW(NR-1,2)
     &           +RN(NR,  2)*RT(NR,  2)
     &           +RN(NR,  3)*RT(NR,  3)
     &           +RN(NR,  4)*RT(NR,  4)
     &           +RW(NR,  1)+RW(NR,  2)
            RNM = RN(NR-1,2)+RN(NR-1,3)+RN(NR-1,4)
     &           +RN(NR,  2)+RN(NR,  3)+RN(NR, 4)
            RPP = RNTP+PNSS(1)   *PTS(1)
            RPM = RNTM+RN(NR-1,1)*RT(NR-1,1)
     &                +RN(NR  ,1)*RT(NR  ,1)
            RPEP= PNSS(1)   *PTS(1)
            RPEM= 0.5D0*(RN(NR-1,1)*RT(NR-1,1)
     &                  +RN(NR  ,1)*RT(NR  ,1))
            RNTM= 0.5D0*RNTM
            RNM = 0.5D0*RNM
            RPM = 0.5D0*RPM
C 
            DTE = (PTS(1) -RT(NR,1))/DR
            DNE = (PNSS(1)-RN(NR,1))/DR
            ZEFFL=ZEFF(NR)
            EZOHL=EZOH(NR)
C
            DPP = (RPP-RPM)/DR
C
            TI  = RNTP/RNP
            DTI = (RNTP/RNP-RNTM/RNM)/DR
         ELSE
            ANE    = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
            AND    = 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
            ANT    = 0.5D0*(RN(NR+1,3)+RN(NR  ,3))
            ANA    = 0.5D0*(RN(NR+1,4)+RN(NR  ,4))
            TE     = 0.5D0*(RT(NR+1,1)+RT(NR  ,1))
            TD     = 0.5D0*(RT(NR+1,2)+RT(NR  ,2))
            TT     = 0.5D0*(RT(NR+1,3)+RT(NR  ,3))
            TA     = 0.5D0*(RT(NR+1,4)+RT(NR  ,4))
C
            RNTP= RN(NR+1,2)*RT(NR+1,2)
     &           +RN(NR+1,3)*RT(NR+1,3)
     &           +RN(NR+1,4)*RT(NR+1,4)
     &           +RW(NR+1,1)+RW(NR+1,2)
            RNP = RN(NR+1,2)+RN(NR+1,3)+RN(NR+1,4)
            RNTM= RN(NR,  2)*RT(NR,  2)
     &           +RN(NR,  3)*RT(NR,  3)
     &           +RN(NR,  4)*RT(NR,  4)
     &           +RW(NR,  1)+RW(NR,  2)
            RNM = RN(NR,  2)+RN(NR,  3)+RN(NR, 4)
            RPP = RNTP+RN(NR+1,1)*RT(NR+1,1)
            RPM = RNTM+RN(NR  ,1)*RT(NR  ,1)
            RPEP= RN(NR+1,1)*RT(NR+1,1)
            RPEM= RN(NR  ,1)*RT(NR  ,1)
C
            DTE = (RT(NR+1,1)-RT(NR  ,1))/DR
            DNE = (RN(NR+1,1)-RN(NR  ,1))/DR
            ZEFFL  = 0.5D0*(ZEFF(NR+1)+ZEFF(NR))
            EZOHL  = 0.5D0*(EZOH(NR+1)+EZOH(NR))
C
            DPP = (RPP-RPM)/DR
C
            TI  = 0.5D0*(RNTP/RNP+RNTM/RNM)
            DTI = (RNTP/RNP-RNTM/RNM)/DR
         ENDIF
C
         IF(NR.LE.1) THEN
            RPI4= RN(NR,  2)*RT(NR,  2)
     &           +RN(NR,  3)*RT(NR,  3)
     &           +RN(NR,  4)*RT(NR,  4)
     &           +RW(NR,  1)+RW(NR,  2)
         ELSE
            RPI4= RN(NR-1,2)*RT(NR-1,2)
     &           +RN(NR-1,3)*RT(NR-1,3)
     &           +RN(NR-1,4)*RT(NR-1,4)
     &           +RW(NR-1,1)+RW(NR-1,2)
         ENDIF
            RPI3= RN(NR,  2)*RT(NR,  2)
     &           +RN(NR,  3)*RT(NR,  3)
     &           +RN(NR,  4)*RT(NR,  4)
     &           +RW(NR,  1)+RW(NR,  2)
         IF(NR.GE.NR-1) THEN
            RPI2= RN(NR,  2)*RT(NR,  2)
     &           +RN(NR,  3)*RT(NR,  3)
     &           +RN(NR,  4)*RT(NR,  4)
     &           +RW(NR,  1)+RW(NR,  2)
         ELSE
            RPI2= RN(NR+1,2)*RT(NR+1,2)
     &           +RN(NR+1,3)*RT(NR+1,3)
     &           +RN(NR+1,4)*RT(NR+1,4)
     &           +RW(NR+1,1)+RW(NR+1,2)
         ENDIF
         IF(NR.GE.NR-2) THEN
            RPI1= RN(NR,  2)*RT(NR,  2)
     &           +RN(NR,  3)*RT(NR,  3)
     &           +RN(NR,  4)*RT(NR,  4)
     &           +RW(NR,  1)+RW(NR,  2)
         ELSE
            RPI1= RN(NR+2,2)*RT(NR+2,2)
     &           +RN(NR+2,3)*RT(NR+2,3)
     &           +RN(NR+2,4)*RT(NR+2,4)
     &           +RW(NR+2,1)+RW(NR+2,2)
         ENDIF
         RPIM=0.5D0*(RPI1+RPI2)
         RPI0=0.5D0*(RPI2+RPI3)
         RPIP=0.5D0*(RPI3+RPI4)
         DPPP=(RPIP-2*RPI0+RPIM)/(DR*DR)
C
         IF(NR.EQ.1) THEN
            DQ = (QP(NR+1)-QP(NR))/(1.5D0*DR)
         ELSEIF(NR.EQ.NRMAX) THEN
            DQ = (QP(NR)-QP(NR-1))/DR
         ELSE
            DQ = (QP(NR+1)-QP(NR-1))/(2.D0*DR)
         ENDIF
         QL = QP(NR)
         RL = RG(NR)
C
         VTE = SQRT(ABS(TE*RKEV/AME))
C
         IF(ABS(DTI).GT.1.D-32) THEN
            CLT=TI/DTI 
         ELSE
            IF(DTI.GE.0.D0) THEN
               CLT = 1.D32
            ELSE
               CLT =-1.D32
            ENDIF
         ENDIF
C
         IF(ABS(DNE).GT.1.D-32) THEN
            CLN=ANE/DNE
         ELSE
            IF(DNE.GE.0.D0) THEN
               CLN = 1.D32
            ELSE
               CLN =-1.D32
            ENDIF
         ENDIF
C
         IF(ABS(RPEP-RPEM).GT.1.D-32) THEN
            CLPE=0.5D0*DR*(RPEP+RPEM)/(RPEP-RPEM)
         ELSE
            IF(RPEP-RPEM.GE.0.D0) THEN
               CLPE = 1.D32
            ELSE
               CLPE =-1.D32
            ENDIF
         ENDIF
C
         IF(ABS(DQ).GT.1.D-32) THEN
            CLS=RR*QL*QL/(DQ*RL)
         ELSE
            IF(DQ.GE.0.D0) THEN
               CLS = 1.D32
            ELSE
               CLS =-1.D32
            ENDIF
         ENDIF
C
         COEF = 6.D0*PI*SQRT(2.D0*PI)*AEPS0**2/(1.D20*AEE**4*15.D0)
         TAUE = COEF*SQRT(AME)*(ABS(TE)*RKEV)**1.5D0/ANE
         ANYUE = 0.5D0*(1.D0+ZEFFL)/TAUE
C
         ROUS = DSQRT(ABS(TE)*RKEV/AMD)/OMEGAD
         PPK  = 0.3D0/ROUS
         EPS  = RL/RR
         OMEGAS  = PPK*TE*RKEV/(AEE*BB*ABS(CLN))
         OMEGATT = DSQRT(2.D0)*VTE/(RR*QL)
C
C        *****  0.GE.MDLKAI.LT.10 : CONSTANT COEFFICIENT MODEL *****
C        ***** 10.GE.MDLKAI.LT.20 : DRIFT WAVE (+ITG +ETG) MODEL *****
C        ***** 20.GE.MDLKAI.LT.30 : REBU-LALLA MODEL *****
C        ***** 30.GE.MDLKAI.LT.40 : CURRENT-DIFFUSIVITY DRIVEN MODEL *****
C                                                                  
         IF(MDLKAI.LT.10) THEN
C           *****  MDLKAI.EQ. 0   : CONSTANT *****
C           *****  MDLKAI.EQ. 1   : CONSTANT/(1-A*r**2) *****
C           *****  MDLKAI.EQ. 2   : CONSTANT*(dTi/dr)**B/(1-A*r**2) *****
C           *****  MDLKAI.EQ. 3   : CONSTANT*(dTi/dr)**B*Ti**C *****
C
            IF(MDLKAI.EQ.0) THEN
               AKDWL=CK0
            ELSEIF(MDLKAI.EQ.1) THEN
               AKDWL=CK0/(1.D0-CKALFA*RL**2/RA**2)
            ELSEIF(MDLKAI.EQ.2) THEN
               AKDWL=CK0/(1.D0-CKALFA*RL**2/RA**2)
     &                  *(ABS(DTI)*RA)**CKBETA
            ELSEIF(MDLKAI.EQ.3) THEN
               AKDWL=CK0*(ABS(DTI)*RA)**CKBETA*ABS(TI)**CKGUMA
            ELSE                                           
               WRITE(6,*) 'XX INVALID MDLAD : ',MDLAD
               AKDWL=0.D0
            ENDIF
            AKDW(NR,1)=AKDWL
            AKDW(NR,2)=AKDWL
            AKDW(NR,3)=AKDWL
            AKDW(NR,4)=AKDWL
C
            VGR1(NR,1)=TI
            VGR1(NR,2)=0.D0
            VGR2(NR,1)=DTI
            VGR2(NR,2)=0.D0
            VGR3(NR,1)=CLN
            VGR3(NR,2)=CLT
            VGR4(NR,1)=CLS
            VGR4(NR,2)=CLS
C
         ELSEIF(MDLKAI.LT.20) THEN
C
            ETAI=CLN/CLT
            IF(MDLKAI.EQ.10) THEN
               ETAC   = 1.D0
               RRSTAR = RR
               FDREV  = 1.D0
               HETA   = 1.D0
            ELSEIF(MDLKAI.EQ.11.OR.MDLKAI.EQ.12.OR.MDLKAI.EQ.13) THEN
               ETAC   = 1.D0
               RRSTAR = RR
               FDREV  = 1.D0
CCCCCC            ARG = 6.D0*(ETAI(NR)-ETAC(NR))*RA/CLN(NR)
               ARG = 6.D0*(ETAI-ETAC)
               IF(ARG.LE.-150.D0) THEN
                  HETA = 0.D0 
               ELSEIF(ARG.GT.150.D0) THEN
                  HETA = 1.D0
               ELSE
                  HETA = 1.D0/(1.D0+EXP(-ARG))
               ENDIF
            ELSEIF(MDLKAI.EQ.14) THEN
               IF(ABS(CLN)/RR.GE.0.2D0) THEN
                  ETAC = 1.D0+2.5D0*(ABS(CLN)/RR-0.2D0)
               ELSE
                  ETAC = 1.D0
               ENDIF
               RRSTAR = RR
               FDREV  = 1.D0
               ARG = 6.D0*(ETAI-ETAC)*RA/CLN
               IF(ARG.LE.-150.D0) THEN
                  HETA = 0.D0
               ELSEIF(ARG.GT.150.D0) THEN
                  HETA = 1.D0
               ELSE
                  HETA = 1.D0/(1.D0+EXP(-ARG))
               ENDIF
            ELSEIF(MDLKAI.EQ.15) THEN
               ETAC    = 1.D0
               RRSTAX  = ABS(RR*(1.D0-(2.D0+1.D0   /QL**2)*EPS)
     &                   /(1.D0-(2.D0+RKAP**2/QL**2)*EPS))
               RREFF = RRSTAX/(1.2D0*ABS(CLN))
               FDREV = SQRT(2.D0*PI)*RREFF**1.5D0*(RREFF-1.5D0)
     &                *EXP(-RREFF)
               RRSTAR= MIN(RRSTAX,CLS)
               FDREV = MAX(FDREV,ANYUE/(EPS*OMEGAS)) 
               ARG = 6.D0*(ETAI-ETAC)
               IF(ARG.LE.-150.D0) THEN
                  HETA = 0.D0
               ELSEIF(ARG.GT.150.D0) THEN
                  HETA = 1.D0
               ELSE
                  HETA = 1.D0/(1.D0+EXP(-ARG))
               ENDIF
            ELSE
               WRITE(6,*) 'XX INVALID MDLAD : ',MDLAD
               RRSTAR = RR 
               FDREV  = 1.D0
               HETA   = 1.D0
            ENDIF
            DEDW = CK0*2.5D0*OMEGAS/PPK**2
     &            *(SQRT(EPS)*MIN(FDREV,EPS*OMEGAS/ANYUE)
     &             +OMEGAS/OMEGATT*MAX(1.D0,ANYUE/OMEGATT))
C
            DIDW = CK0*2.5D0*OMEGAS/PPK**2*HETA
     &            *SQRT(2.D0*ABS(TI)*ABS(ETAI)*ABS(CLN)/(TE*RRSTAR))
C
            AKDW(NR,1) = CDW(1)*DEDW+CDW(2)*DIDW
            AKDW(NR,2) = CDW(3)*DEDW+CDW(4)*DIDW
            AKDW(NR,3) = CDW(5)*DEDW+CDW(6)*DIDW
            AKDW(NR,4) = CDW(7)*DEDW+CDW(8)*DIDW
            IF(MDLKAI.EQ.12) THEN
               AKDW(NR,1) = AKDW(NR,1)*QL
               AKDW(NR,2) = AKDW(NR,2)*QL
               AKDW(NR,3) = AKDW(NR,3)*QL
               AKDW(NR,4) = AKDW(NR,4)*QL
            ELSEIF(MDLKAI.EQ.13) THEN
               AKDW(NR,1) = AKDW(NR,1)*(1+QL**2)
               AKDW(NR,2) = AKDW(NR,2)*(1+QL**2)
               AKDW(NR,3) = AKDW(NR,3)*(1+QL**2)
               AKDW(NR,4) = AKDW(NR,4)*(1+QL**2)
            ENDIF
            VGR1(NR,1)=PLOG(DEDW,1.D-2,1.D2)
            VGR1(NR,2)=PLOG(DIDW,1.D-2,1.D2)
            VGR2(NR,1)=ETAI
            VGR2(NR,2)=ETAC
            VGR3(NR,1)=CLN
            VGR3(NR,2)=CLT
            VGR4(NR,1)=CLS
            VGR4(NR,2)=CLS
C
         ELSEIF(MDLKAI.LT.30) THEN
C
            IF(MDLKAI.EQ.20) THEN
               CRTCL=6.D-2*SQRT(EZOHL*BB**3/(ANE*1.D20*SQRT(TE*RKEV)))
     &              *SQRT(AEE**2/(AMYU0*SQRT(AME)))/(QL*RKEV)
               IF(ABS(DTE).LE.CRTCL.OR.
     &            DQ.LE.0.D0.OR.
     &            NR.EQ.1) THEN
                  DKAI=0.D0
               ELSE
                  DKAI=1.5D-1*(ABS(DTE)/TE
     &                 +2.D0*ABS(DNE)/ANE)
     &                 *SQRT(TE/TI)*(1.D0/EPS)
     &                 *(QL**2/(DQ*BB*SQRT(RR)))*VC**2
     &                 *SQRT(AMYU0*1.5D0*AMM)
     &                 *(1.D0-CRTCL/ABS(DTE))
               ENDIF
               AKDW(NR,1)=                  DKAI
               AKDW(NR,2)=ZEFFL*SQRT(TE/TD)*DKAI
               AKDW(NR,3)=ZEFFL*SQRT(TE/TT)*DKAI
               AKDW(NR,4)=ZEFFL*SQRT(TE/TA)*DKAI
            ELSE                                           
               WRITE(6,*) 'XX INVALID MDLAD : ',MDLAD
               AKDW(NR,1)=0.D0
               AKDW(NR,2)=0.D0
               AKDW(NR,3)=0.D0
               AKDW(NR,4)=0.D0
            ENDIF
            VGR1(NR,1)=DTE
            VGR1(NR,2)=CRTCL
            VGR2(NR,1)=ETAI
            VGR2(NR,2)=ETAC
            VGR3(NR,1)=CLN
            VGR3(NR,2)=CLT
            VGR4(NR,1)=CLS
            VGR4(NR,2)=CLS
C
         ELSEIF(MDLKAI.LT.40) THEN
C
            PNI=AND+ANT+ANA
            AMI=(AMD*AND+AMT*ANT+AMA*ANA)/PNI
C
            VA=SQRT(BB**2/(AMYU0*ANE*1.D20*AMI))
            WPE2=ANE*1.D20*AEE*AEE/(AME*AEPS0)
            S=RL*DQ/QL
            DELTA2=VC**2/WPE2
            DBDR=DPP*1.D20*RKEV*RA/(BB**2/(2*AMYU0))
            ALFA=-QL*QL*DBDR*RR/RA
            RKCV=-EPS*(1.D0-1.D0/(QL*QL))
C
            IF(MDLKAI.EQ.30) THEN
               FS=1.D0/(1.7D0+SQRT(6.D0)*S)
               AKDWEL=CK0*FS*SQRT(ABS(ALFA))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK0*FS*SQRT(ABS(ALFA))**3*DELTA2*VA/(QL*RR)
            ELSEIF(MDLKAI.EQ.31) THEN
               FS=TRCOFS(S,ALFA,RKCV)
               AKDWEL=CK0*FS*SQRT(ABS(ALFA))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK0*FS*SQRT(ABS(ALFA))**3*DELTA2*VA/(QL*RR)
            ELSEIF(MDLKAI.EQ.32) THEN
               FS=TRCOFSX(S,ALFA,RKCV,RA/RR)
               AKDWEL=CK0*FS*SQRT(ABS(ALFA))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK0*FS*SQRT(ABS(ALFA))**3*DELTA2*VA/(QL*RR)
            ELSEIF(MDLKAI.EQ.33) THEN
               FS=TRCOFS(S,ALFA,RKCV)
               DBDRR=DPPP*1.D20*RKEV*RA*RA/(BB**2/(2*AMYU0))
               DELTAE=SQRT(DELTA2)
               SL=SQRT(S**2+0.1D0**2)
               WE1=SQRT(PA(2)/PA(1))*(QL*RR*DELTAE)/(2*SL*RA*RA)*DBDRR
               RG1=10.D0
C               RG1=8.D0
               FS=FS/(1.D0+RG1*WE1*WE1)
               AKDWEL=CK0*FS*SQRT(ABS(ALFA))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK0*FS*SQRT(ABS(ALFA))**3*DELTA2*VA/(QL*RR)
C
            ELSEIF(MDLKAI.EQ.34) THEN
               AEI=(PZ(2)*AND+PZ(3)*ANT+PZ(4)*ANA)*AEE/PNI
               WCI=AEI*BB/AMI
               PTI=(TD*AND+TT*ANT+TA*ANA)/PNI
               VTI=SQRT(ABS(PTI*RKEV/AMI))
               RHOI=VTI/WCI
C
               FS=TRCOFT(S,ALFA,RKCV,RA/RR)
               SA=S-ALFA
               RNST2=0.5D0/((1.D0-2.D0*SA+3.D0*SA*SA)*FS)
               RKPP2=RNST2/(FS*ABS(ALFA)*DELTA2)
C
               SLAMDA=RKPP2*RHOI**2
               RLAMDA=RLAMBDA(SLAMDA)
               OMEGAS= SQRT(RKPP2)*TE*RKEV/(AEE*BB*ABS(CLPE))
               TAUAP=(QL*RR)/VA
               OMEGASS=(OMEGAS*TAUAP)/(RNST2*SQRT(ALFA))
               AKDWEL=CK0*FS*SQRT(ABS(ALFA))**3*DELTA2*VA/(QL*RR)
     &               /(RLAMDA*(1.D0+OMEGASS**2))
               AKDWIL=CK0*FS*SQRT(ABS(ALFA))**3*DELTA2*VA/(QL*RR)
     &               /(1.D0+OMEGASS**2)
            ELSE                                           
               WRITE(6,*) 'XX INVALID MDLAD : ',MDLAD
               FS=1.D0
               AKDWEL=CK0*FS*SQRT(ABS(ALFA))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK0*FS*SQRT(ABS(ALFA))**3*DELTA2*VA/(QL*RR)
            ENDIF
C
            AKDW(NR,1)=AKDWEL
            AKDW(NR,2)=AKDWIL
            AKDW(NR,3)=AKDWIL
            AKDW(NR,4)=AKDWIL
C
            VGR1(NR,1)=FS
            VGR1(NR,2)=S
            VGR1(NR,3)=ALFA
            VGR2(NR,1)=RNST2
            VGR2(NR,2)=OMEGASS
            VGR2(NR,3)=0.D0
            VGR3(NR,1)=SLAMDA
            VGR3(NR,2)=0.D0
            VGR3(NR,3)=0.D0
            VGR4(NR,1)=RLAMDA
            VGR4(NR,2)=1.D0/(1.D0+OMEGASS**2)
            VGR4(NR,3)=1.D0/(1.D0+RG1*WE1*WE1)
         ELSE                                           
            WRITE(6,*) 'XX INVALID MDLAD : ',MDLAD       
            AKDW(NR,1)=0.D0
            AKDW(NR,2)=0.D0
            AKDW(NR,3)=0.D0
            AKDW(NR,4)=0.D0
         ENDIF
  100 CONTINUE
      RETURN
      END
C
C     ***********************************************************
C
      SUBROUTINE TRCFNC
C
      INCLUDE 'trcomm.h'
C
C     ZEFF=1
C
      DATA RK11,RA11,RB11,RC11/1.04D0,2.01D0,1.53D0,0.89D0/
      DATA RK12,RA12,RB12,RC12/1.20D0,0.76D0,0.67D0,0.56D0/
      DATA RK22,RA22,RB22,RC22/2.55D0,0.45D0,0.43D0,0.43D0/
      DATA RK13,RA13,RB13,RC13/2.30D0,1.02D0,1.07D0,1.07D0/
      DATA RK23,RA23,RB23,RC23/4.19D0,0.57D0,0.61D0,0.61D0/
      DATA RK33,RA33,RB33,RC33/1.83D0,0.68D0,0.32D0,0.66D0/
      DATA RK2 ,RA2 ,RB2 ,RC2 /0.66D0,1.03D0,0.31D0,0.74D0/
C
      AMD=PA(2)*AMM
      AMT=PA(3)*AMM
      AMA=PA(4)*AMM
C
      DO 100 NR=1,NRMAX
         IF(NR.EQ.NRMAX) THEN
            ANE=PNSS(1)
            AND=PNSS(2)
            ANT=PNSS(3)
            ANA=PNSS(4)
            TE=PTS(1)
            TD=PTS(2)
            TT=PTS(3)
            TA=PTS(4)
C
            ZEFFL=ZEFF(NR)
C
         ELSE
            ANE    = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
            AND    = 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
            ANT    = 0.5D0*(RN(NR+1,3)+RN(NR  ,3))
            ANA    = 0.5D0*(RN(NR+1,4)+RN(NR  ,4))
            TE     = 0.5D0*(RT(NR+1,1)+RT(NR  ,1))
            TD     = 0.5D0*(RT(NR+1,2)+RT(NR  ,2))
            TT     = 0.5D0*(RT(NR+1,3)+RT(NR  ,3))
            TA     = 0.5D0*(RT(NR+1,4)+RT(NR  ,4))
C
            ZEFFL  = 0.5D0*(ZEFF(NR+1)+ZEFF(NR))
C
         ENDIF
C
         QL = QP(NR)
         RL = RG(NR)
C
         VTE = SQRT(ABS(TE*RKEV/AME))
         VTD = SQRT(ABS(TD*RKEV/AMD))
         VTT = SQRT(ABS(TT*RKEV/AMT))
         VTA = SQRT(ABS(TA*RKEV/AMA))
C
         COEF = 6.D0*PI*SQRT(2.D0*PI)*AEPS0**2/(1.D20*AEE**4*15.D0)
         TAUE = COEF*SQRT(AME)*(ABS(TE)*RKEV)**1.5D0/ANE
         TAUD = COEF*SQRT(AMD)*(ABS(TD)*RKEV)**1.5D0/ANE
         TAUT = COEF*SQRT(AMT)*(ABS(TT)*RKEV)**1.5D0/ANE
         TAUA = COEF*SQRT(AMA)*(ABS(TA)*RKEV)**1.5D0/ANE
C
         EPS  = RL/RR
C
C     ***** NEOCLASSICAL TRANSPORT (HINTON, HAZELTINE) *****
C
C
C     ***** OLD AKNC *****
C
      IF(MDLKNC.EQ.1) THEN
C
         EPS=RG(NR)/RR
         EPSS=SQRT(EPS)**3
         RHOE2=2.D0*AME*ABS(TE)*RKEV/(PZ(1)*AEE*BP(NR))**2
         RHOD2=2.D0*AMD*ABS(TD)*RKEV/(PZ(2)*AEE*BP(NR))**2
         RHOT2=2.D0*AMT*ABS(TT)*RKEV/(PZ(3)*AEE*BP(NR))**2
         RHOA2=2.D0*AMA*ABS(TA)*RKEV/(PZ(4)*AEE*BP(NR))**2
C
         COEF = 12.D0*PI*SQRT(PI)*AEPS0**2
     &         /(ANE*1.D20*ZEFFL*AEE**4*15.D0)
         TAUE = COEF*SQRT(AME)*(ABS(TE)*RKEV)**1.5D0/SQRT(2.D0)
         TAUD = COEF*SQRT(AMD)*(ABS(TD)*RKEV)**1.5D0/PZ(2)**2
         TAUT = COEF*SQRT(AMT)*(ABS(TT)*RKEV)**1.5D0/PZ(3)**2
         TAUA = COEF*SQRT(AMA)*(ABS(TA)*RKEV)**1.5D0/PZ(4)**2
C
         RNUE=ABS(QP(NR))*RR/(TAUE*VTE*EPSS)
         RNUD=ABS(QP(NR))*RR/(TAUD*VTD*EPSS)
         RNUT=ABS(QP(NR))*RR/(TAUT*VTT*EPSS)
         RNUA=ABS(QP(NR))*RR/(TAUA*VTA*EPSS)
C
C         RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)
C     &              +(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
C         RK12E=RK12*(1.D0/(1.D0+RA12*SQRT(RNUE)+RB12*RNUE)
C     &              +(EPSS*RC12)**2/RB12*RNUE/(1.D0+RC12*RNUE*EPSS))
         RK22E=RK22*(1.D0/(1.D0+RA22*SQRT(RNUE)+RB22*RNUE)
     &              +(EPSS*RC22)**2/RB22*RNUE/(1.D0+RC22*RNUE*EPSS))
C         RK13E=RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)
C     &             /(1.D0+RC13*RNUE*EPSS)
C         RK23E=RK23/(1.D0+RA23*SQRT(RNUE)+RB23*RNUE)
C     &             /(1.D0+RC23*RNUE*EPSS)
C         RK33E=RK33/(1.D0+RA33*SQRT(RNUE)+RB33*RNUE)
C     &             /(1.D0+RC33*RNUE*EPSS)
C
         RK2D =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUD)+RB2 *RNUD)
     &              +(EPSS*RC2 )**2/RB2 *RNUD/(1.D0+RC2 *RNUD*EPSS))
         RK2T =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUT)+RB2 *RNUT)
     &              +(EPSS*RC2 )**2/RB2 *RNUT/(1.D0+RC2 *RNUT*EPSS))
         RK2A =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUA)+RB2 *RNUA)
     &              +(EPSS*RC2 )**2/RB2 *RNUA/(1.D0+RC2 *RNUA*EPSS))
C         RK3D=((1.17-0.35*SQRT(RNUD))/(1.D0+0.7*SQRT(RNUD))
C     &         -2.1*(RNUD*EPSS)**2)/(1+(RNUD*EPSS)**2)
C         RK3T=((1.17-0.35*SQRT(RNUT))/(1.D0+0.7*SQRT(RNUT))
C     &         -2.1*(RNUT*EPSS)**2)/(1+(RNUT*EPSS)**2)
C         RK3A=((1.17-0.35*SQRT(RNUA))/(1.D0+0.7*SQRT(RNUA))
C     &         -2.1*(RNUA*EPSS)**2)/(1+(RNUA*EPSS)**2)
C
         AKNC(NR,1) = SQRT(EPS)*RHOE2/TAUE*RK22E
         AKNC(NR,2) = SQRT(EPS)*RHOD2/TAUD*RK2D
         AKNC(NR,3) = SQRT(EPS)*RHOT2/TAUT*RK2T
         AKNC(NR,4) = SQRT(EPS)*RHOA2/TAUA*RK2A
C
      ELSE
C
C     ***** NEW AKNC *****
C
         DELDA=0.D0
         EPS=RG(NR)/RR
         EPSS=SQRT(EPS)**3
C
         IF(NR.EQ.1) THEN
            Q0=(4.D0*QP(1)-QP(2))/3.D0
            QL=ABS(0.25D0*(3.D0*Q0+QP(NR)))
C            ZEFFL=0.5D0*(ZEFF(NR)+ZEFF(NR+1))
            ZEFFL=ZEFF(NR)
         ELSE
            QL=ABS(0.5D0*(QP(NR-1)+QP(NR)))
            ZEFFL=0.5D0*(ZEFF(NR)+ZEFF(NR-1))
         ENDIF
         RALFA=ZEFFL-1.D0
C
         rLnLamE=15.2D0-DLOG(ANE)*0.5D0+DLOG(ABS(TE))
         rLnLamD=17.3D0-DLOG(AND)*0.5D0+DLOG(ABS(TD))*1.5D0
         rLnLamT=17.3D0-DLOG(ANT)*0.5D0+DLOG(ABS(TT))*1.5D0
         rLnLamA=17.3D0-DLOG(ANA)*0.5D0+DLOG(ABS(TA))*1.5D0
C
         TAUE=6.D0*PI*SQRT(2.D0*PI)*AEPS0**2*DSQRT(AME)
     &             *(ABS(TE)*RKEV)**1.5D0
     &             /(PZ(1)**2*ANE*1.D20*AEE**4*rLnLamE)
         TAUD=12.D0*PI*SQRT(PI)*AEPS0**2*DSQRT(AMD)
     &             *(ABS(TD)*RKEV)**1.5D0
     &             /(PZ(2)**4*AND*1.D20*AEE**4*rLnLamD)
         TAUT=12.D0*PI*SQRT(PI)*AEPS0**2*DSQRT(AMT)
     &             *(ABS(TT)*RKEV)**1.5D0
     &             /(PZ(3)**4*ANT*1.D20*AEE**4*rLnLamT)
         TAUA=12.D0*PI*SQRT(PI)*AEPS0**2*DSQRT(AMA)
     &             *(ABS(TA)*RKEV)**1.5D0
     &             /(PZ(4)**4*ANA*1.D20*AEE**4*rLnLamA)
C
         RNUE=QL*RR/(TAUE*VTE*EPSS)
         RNUD=QL*RR/(TAUD*VTD*EPSS)
         RNUT=QL*RR/(TAUT*VTT*EPSS)
         RNUA=QL*RR/(TAUA*VTA*EPSS)
C
         OMEGACE=(AEE*PZ(1)*BB)/AME
         OMEGACD=(AEE*PZ(2)*BB)/AMD
         OMEGACT=(AEE*PZ(3)*BB)/AMT
         OMEGACA=(AEE*PZ(4)*BB)/AMA
C
         RHOIE=DSQRT(2.D0*ABS(TE)*RKEV/AME)/OMEGACE
         RHOID=DSQRT(2.D0*ABS(TD)*RKEV/AMD)/OMEGACD
         RHOIT=DSQRT(2.D0*ABS(TT)*RKEV/AMT)/OMEGACT
         RHOIA=DSQRT(2.D0*ABS(TA)*RKEV/AMA)/OMEGACA
C
         RMUSE=RNUE*(1.D0+1.54D0*RALFA)
         RMUSD=RNUD*(1.D0+1.54D0*RALFA)
         RMUST=RNUT*(1.D0+1.54D0*RALFA)
         RMUSA=RNUA*(1.D0+1.54D0*RALFA)
C
         F1=(1.D0+1.5D0*(EPS**2+EPS*DELDA)
     &      +3.D0/8.D0*EPS**3*DELDA)
     &      /(1.D0+0.5D0*EPS*DELDA)
         F2=DSQRT(1.D0-EPS**2)*(1.D0+0.5D0*EPS*DELDA)
     &         /(1.D0+DELDA/EPS*(DSQRT(1.D0-EPS**2)-1.D0))
C
         TERM1E=(0.66D0*(1.D0+1.54D0*RALFA)
     &       +(1.88D0*DSQRT(EPS)-1.54D0*EPS)*(1.D0+3.75D0*RALFA))*F1
     &       /(1.D0+1.03D0*DSQRT(RMUSE)+0.31D0*RMUSE)
         TERM1D=(0.66D0*(1.D0+1.54D0*RALFA)
     &       +(1.88D0*DSQRT(EPS)-1.54D0*EPS)*(1.D0+3.75D0*RALFA))*F1
     &       /(1.D0+1.03D0*DSQRT(RMUSD)+0.31D0*RMUSD)
         TERM1T=(0.66D0*(1.D0+1.54D0*RALFA)
     &       +(1.88D0*DSQRT(EPS)-1.54D0*EPS)*(1.D0+3.75D0*RALFA))*F1
     &       /(1.D0+1.03D0*DSQRT(RMUST)+0.31D0*RMUST)
         TERM1A=(0.66D0*(1.D0+1.54D0*RALFA)
     &       +(1.88D0*DSQRT(EPS)-1.54D0*EPS)*(1.D0+3.75D0*RALFA))*F1
     &       /(1.D0+1.03D0*DSQRT(RMUSA)+0.31D0*RMUSA)
C
         TERM2E=0.59D0*RMUSE*EPS/(1.D0+0.74D0*RMUSE*EPSS)
     &       *(1.D0+(1.33D0*RALFA*(1.D0+0.6D0*RALFA))
     &       /(1.D0+1.79D0*RALFA))
     &       *(F1-F2)
         TERM2D=0.59D0*RMUSD*EPS/(1.D0+0.74D0*RMUSD*EPSS)
     &       *(1.D0+(1.33D0*RALFA*(1.D0+0.6D0*RALFA))
     &       /(1.D0+1.79D0*RALFA))
     &       *(F1-F2)
         TERM2T=0.59D0*RMUST*EPS/(1.D0+0.74D0*RMUST*EPSS)
     &       *(1.D0+(1.33D0*RALFA*(1.D0+0.6D0*RALFA))
     &       /(1.D0+1.79D0*RALFA))
     &       *(F1-F2)
         TERM2A=0.59D0*RMUSA*EPS/(1.D0+0.74D0*RMUSA*EPSS)
     &       *(1.D0+(1.33D0*RALFA*(1.D0+0.6D0*RALFA))
     &       /(1.D0+1.79D0*RALFA))
     &       *(F1-F2)
C
         AKNC(NR,1)=(QL**2*RHOIE**2)/(EPSS*TAUE)*(TERM1E+TERM2E)
         AKNC(NR,2)=(QL**2*RHOID**2)/(EPSS*TAUD)*(TERM1D+TERM2D)
         AKNC(NR,3)=(QL**2*RHOIT**2)/(EPSS*TAUT)*(TERM1T+TERM2T)
         AKNC(NR,4)=(QL**2*RHOIA**2)/(EPSS*TAUA)*(TERM1A+TERM2A)
C
      ENDIF
C
         AK(NR,1) = AKDW(NR,1)+CNC*AKNC(NR,1)
         AK(NR,2) = AKDW(NR,2)+CNC*AKNC(NR,2)
         AK(NR,3) = AKDW(NR,3)+CNC*AKNC(NR,3)
         AK(NR,4) = AKDW(NR,4)+CNC*AKNC(NR,4)
C
C
  100 CONTINUE
C
      RETURN
      END
C
C     ***********************************************************
C
      SUBROUTINE TRCFET
C
      INCLUDE 'trcomm.h'
C
C     ZEFF=1
C
      DATA RK11,RA11,RB11,RC11/1.04D0,2.01D0,1.53D0,0.89D0/
      DATA RK12,RA12,RB12,RC12/1.20D0,0.76D0,0.67D0,0.56D0/
      DATA RK22,RA22,RB22,RC22/2.55D0,0.45D0,0.43D0,0.43D0/
      DATA RK13,RA13,RB13,RC13/2.30D0,1.02D0,1.07D0,1.07D0/
      DATA RK23,RA23,RB23,RC23/4.19D0,0.57D0,0.61D0,0.61D0/
      DATA RK33,RA33,RB33,RC33/1.83D0,0.68D0,0.32D0,0.66D0/
      DATA RK2 ,RA2 ,RB2 ,RC2 /0.66D0,1.03D0,0.31D0,0.74D0/
C
      DO 150 NR=1,NRMAX
C
C        ****** CLASSICAL RESISTIVITY ******
C
         ANE=RN(NR,1)
         TE =RT(NR,1)
         ZEFFL=ZEFF(NR)
C
         COEF = 12.D0*PI*SQRT(PI)*AEPS0**2
     &         /(ANE*1.D20*ZEFFL*AEE**4*15.D0)
         TAUE = COEF*SQRT(AME)*(ABS(TE)*RKEV)**1.5D0/SQRT(2.D0)
C
         ETA(NR) = AME/(ANE*1.D20*AEE*AEE*TAUE)
     &             *(0.29D0+0.46D0/(1.08D0+ZEFFL))
C
C        ****** NEOCLASSICAL RESISTIVITY ******
C
         IF(MDLETA.EQ.1) THEN
            EPS=RM(NR)/RR
            EPSS=SQRT(EPS)**3
            IF(NR.EQ.1) THEN
               QL= 0.25D0*(3.D0*Q0+QP(NR))
            ELSE
               QL= 0.5D0*(QP(NR-1)+QP(NR))
            ENDIF
            VTE=SQRT(ABS(TE)*RKEV/AME)
            RNUE=ABS(QP(NR)*RR/(TAUE*VTE*EPSS))
            RK33E=RK33/(1.D0+RA33*SQRT(RNUE)+RB33*RNUE)
     &                /(1.D0+RC33*RNUE*EPSS)
C
            FT     = 1.D0-SQRT(EPS)*RK33E
            ETA(NR)= ETA(NR)/FT
C
C        ****** NEOCLASSICAL RESISTIVITY PART II ******
C
         ELSEIF(MDLETA.EQ.2) THEN
            EPS=RM(NR)/RR
            EPSS=SQRT(EPS)**3
C
            IF(NR.EQ.1) THEN
               QL=ABS(0.25D0*(3.D0*Q0+QP(NR)))
               ZEFFL=0.5D0*(ZEFF(NR+1)+ZEFF(NR))
C               ZEFFL=ZEFF(NR)
            ELSE
               QL=ABS(0.5D0*(QP(NR-1)+QP(NR)))
               ZEFFL=0.5D0*(ZEFF(NR-1)+ZEFF(NR))
            ENDIF
C
         VTE=1.33D+7*DSQRT(ABS(TE))
         FT=1.D0-(1.D0-EPS)**2.D0
     &         /(DSQRT(1.D0-EPS**2.D0)*(1.D0+1.46D0*DSQRT(EPS)))
         rLnLam=15.2D0-DLOG(ANE)/2+DLOG(ABS(TE))
         TAUE=6.D0*PI*SQRT(2*PI)*AEPS0**2*DSQRT(AME)
     &             *(ABS(TE)*RKEV)**1.5D0/(ANE*1.D20*AEE**4*rLnLam)
         RNUSE=RR*QL/(VTE*TAUE*EPSS)
         PHI=FT/(1.D0+(0.58D0+0.20D0*ZEFFL)*RNUSE)
         ETAS=1.65D-9*rLnLam/(ABS(TE))**1.5D0
         CH=0.56D0*(3.D0-ZEFFL)/((3.D0+ZEFFL)*ZEFFL)
C
         ETA(NR)=ETAS*ZEFFL*(1.D0+0.27D0*(ZEFFL-1.D0))
     &            /(1.D0-PHI)*(1.D0-CH*PHI)*(1.D0+0.47D0*(ZEFFL-1.D0))
         ENDIF 
  150 CONTINUE
C
      RETURN
      END
C
C     ***********************************************************
C
      SUBROUTINE TRCFAD
C
      INCLUDE 'trcomm.h'
C
C     ZEFF=1
C
      DATA RK11,RA11,RB11,RC11/1.04D0,2.01D0,1.53D0,0.89D0/
      DATA RK12,RA12,RB12,RC12/1.20D0,0.76D0,0.67D0,0.56D0/
      DATA RK22,RA22,RB22,RC22/2.55D0,0.45D0,0.43D0,0.43D0/
      DATA RK13,RA13,RB13,RC13/2.30D0,1.02D0,1.07D0,1.07D0/
      DATA RK23,RA23,RB23,RC23/4.19D0,0.57D0,0.61D0,0.61D0/
      DATA RK33,RA33,RB33,RC33/1.83D0,0.68D0,0.32D0,0.66D0/
      DATA RK2 ,RA2 ,RB2 ,RC2 /0.66D0,1.03D0,0.31D0,0.74D0/
C
C        ****** AD : PARTICLE DIFFUSION ******
C        ****** AV : PARTICLE PINCH ******
C
      IF(MDLAD.EQ.0) THEN
         DO 200 NS=1,NSM
         DO 200 NR=1,NRMAX
            AD(NR,NS)=0.D0
            AV(NR,NS)=0.D0
  200    CONTINUE
      ELSEIF(MDLAD.EQ.1) THEN
         DO 210 NR=1,NRMAX
            IF(NR.EQ.NRMAX) THEN
               ANE=PNSS(1)
               AND=PNSS(2)
               ANT=PNSS(3)
               ANA=PNSS(4)
            ELSE
               ANE    = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
               AND    = 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
               ANT    = 0.5D0*(RN(NR+1,3)+RN(NR  ,3))
               ANA    = 0.5D0*(RN(NR+1,4)+RN(NR  ,4))
            ENDIF
            AD(NR,2) = PA(2)**ALP(2)*PZ(2)**ALP(3)*AD0
            AD(NR,3) = PA(3)**ALP(2)*PZ(3)**ALP(3)*AD0
            AD(NR,4) = PA(4)**ALP(2)*PZ(4)**ALP(3)*AD0
            AD(NR,1) =(PZ(2)*AND*AD(NR,2)
     &                +PZ(3)*ANT*AD(NR,3)
     &                +PZ(4)*ANA*AD(NR,4))/(AND+ANT+ANA)
C
            RX   = ALP(1)*RG(NR)/RA
            PROF0 = 1.D0-RX**PROFN1
            IF(PROF0.LE.0.D0) THEN
               PROF1=0.D0
               PROF2=0.D0
            ELSE
               PROF1=PROF0**PROFN2
               PROF2=PROFN2*PROF0**(PROFN2-1.D0)
            ENDIF
            PROF   = PROF1+PNSS(1)/(PN(1)-PNSS(1))
            DPROF  = -PROFN1*RX**(PROFN1-1.D0)*PROF2
C
            AV(NR,1) =AD(NR,1)*DPROF/PROF
            AV(NR,2) =AD(NR,2)*DPROF/PROF
            AV(NR,3) =AD(NR,3)*DPROF/PROF
            AV(NR,4) =AD(NR,4)*DPROF/PROF
  210    CONTINUE
      ELSEIF(MDLAD.EQ.2) THEN
         DO 220 NR=1,NRMAX
            IF(NR.EQ.NRMAX) THEN
               ANE=PNSS(1)
               AND=PNSS(2)
               ANT=PNSS(3)
               ANA=PNSS(4)
            ELSE
               ANE    = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
               AND    = 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
               ANT    = 0.5D0*(RN(NR+1,3)+RN(NR  ,3))
               ANA    = 0.5D0*(RN(NR+1,4)+RN(NR  ,4))
            ENDIF
            AD(NR,2) = AD0*AKDW(NR,2)
            AD(NR,3) = AD0*AKDW(NR,3)
            AD(NR,4) = AD0*AKDW(NR,4)
            AD(NR,1) =(PZ(2)*AND*AD(NR,2)
     &                +PZ(3)*ANT*AD(NR,3)
     &                +PZ(4)*ANA*AD(NR,4))/(AND+ANT+ANA)
C
            RX   = ALP(1)*RG(NR)/RA
            PROF0 = 1.D0-RX**PROFN1
            IF(PROF0.LE.0.D0) THEN
               PROF1=0.D0
               PROF2=0.D0
            ELSE
               PROF1=PROF0**PROFN2
               PROF2=PROFN2*PROF0**(PROFN2-1.D0)
            ENDIF
            PROF   = PROF1*(PN(1)-PNS(1))+PNS(1)
            DPROF  = -PROFN1*RX**(PROFN1-1.D0)*PROF2
     &                *(PN(1)-PNS(1))*ALP(1)/RA *1.5D0
C
            AV(NR,1) =AD(NR,1)*DPROF/PROF
            AV(NR,2) =AD(NR,2)*DPROF/PROF
            AV(NR,3) =AD(NR,3)*DPROF/PROF
            AV(NR,4) =AD(NR,4)*DPROF/PROF
  220    CONTINUE
      ELSEIF(MDLAD.EQ.3) THEN
         DO 230 NR=1,NRMAX
            IF(NR.EQ.NRMAX) THEN
               ANE = PNSS(1)
               TE  = PTS(1)
            ELSE
               ANE = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
               TE  = 0.5D0*(RT(NR+1,1)+RT(NR  ,1))
            ENDIF
C
            ADDW(NR,1) = AD0*AKDW(NR,1)
            ADDW(NR,2) = AD0*AKDW(NR,1)
            ADDW(NR,3) = AD0*AKDW(NR,1)
            ADDW(NR,4) = AD0*AKDW(NR,1)
C
            ZEFFL=ZEFF(NR)
            COEF = 12.D0*PI*SQRT(PI)*AEPS0**2
     &            /(ANE*1.D20*ZEFFL*AEE**4*15.D0)
            BPL   = BP(NR)
            QPL   = QP(NR)
            IF(QPL.GT.100.D0) QPL=100.D0
            EZOHL = EZOH(NR)
            EPS   = RG(NR)/RR
            EPSS  = SQRT(EPS)**3
            VTE   = SQRT(TE*RKEV/AME)
            TAUE  = COEF*SQRT(AME)*(TE*RKEV)**1.5D0/SQRT(2.D0)
            RNUE  = ABS(QPL)*RR/(TAUE*VTE*EPSS)
            RHOE2 = 2.D0*AME*TE*RKEV/(PZ(1)*AEE*BPL)**2
C
            RK13E = RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)
     &             /(1.D0+RC13*RNUE*EPSS)
C
            AV(NR,1) = -(RK13E*SQRT(EPS)*EZOHL) / BPL
            AV(NR,2) = -(RK13E*SQRT(EPS)*EZOHL) / BPL
            AV(NR,3) = -(RK13E*SQRT(EPS)*EZOHL) / BPL
            AV(NR,4) = -(RK13E*SQRT(EPS)*EZOHL) / BPL
C
            RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)
     &           +(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
C
            ADNC(NR,1) = SQRT(EPS)*RHOE2/TAUE*RK11E
            ADNC(NR,2) = SQRT(EPS)*RHOE2/TAUE*RK11E
            ADNC(NR,3) = SQRT(EPS)*RHOE2/TAUE*RK11E
            ADNC(NR,4) = SQRT(EPS)*RHOE2/TAUE*RK11E
C
            AD(NR,1) = ADDW(NR,1)+CNC*ADNC(NR,1)
            AD(NR,2) = ADDW(NR,2)+CNC*ADNC(NR,2)
            AD(NR,3) = ADDW(NR,3)+CNC*ADNC(NR,3)
            AD(NR,4) = ADDW(NR,4)+CNC*ADNC(NR,4)
C
  230    CONTINUE
      ELSEIF(MDLAD.EQ.4) THEN
         DO 240 NR=1,NRMAX
            IF(NR.EQ.NRMAX) THEN
               ANE = PNSS(1)
               TE  = PTS(1)
            ELSE
               ANE = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
               TE  = 0.5D0*(RT(NR+1,1)+RT(NR  ,1))
            ENDIF
C
            ADDW(NR,1) = AD0*AKDW(NR,1)
            ADDW(NR,2) = AD0*AKDW(NR,1)
            ADDW(NR,3) = AD0*AKDW(NR,1)
            ADDW(NR,4) = AD0*AKDW(NR,1)
C
            ZEFFL=ZEFF(NR)
            COEF = 12.D0*PI*SQRT(PI)*AEPS0**2
     &            /(ANE*1.D20*ZEFFL*AEE**4*15.D0)
            BPL   = BP(NR)
            QPL   = QP(NR)
            IF(QPL.GT.100.D0) QPL=100.D0
            EZOHL = EZOH(NR)
            EPS   = RG(NR)/RR
            EPSS  = SQRT(EPS)**3
            VTE   = SQRT(TE*RKEV/AME)
            TAUE  = COEF*SQRT(AME)*(TE*RKEV)**1.5D0/SQRT(2.D0)
            RNUE  = ABS(QPL)*RR/(TAUE*VTE*EPSS)
            RHOE2 = 2.D0*AME*TE*RKEV/(PZ(1)*AEE*BPL)**2
C
            RK13E = RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)
     &             /(1.D0+RC13*RNUE*EPSS)
C
            AVNC(NR,1) = -(RK13E*SQRT(EPS)*EZOHL) / BPL
            AVNC(NR,2) = -(RK13E*SQRT(EPS)*EZOHL) / BPL
            AVNC(NR,3) = -(RK13E*SQRT(EPS)*EZOHL) / BPL
            AVNC(NR,4) = -(RK13E*SQRT(EPS)*EZOHL) / BPL
C
            RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)
     &           +(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
C
            ADNC(NR,1) = SQRT(EPS)*RHOE2/TAUE*RK11E
            ADNC(NR,2) = SQRT(EPS)*RHOE2/TAUE*RK11E
            ADNC(NR,3) = SQRT(EPS)*RHOE2/TAUE*RK11E
            ADNC(NR,4) = SQRT(EPS)*RHOE2/TAUE*RK11E
C
            AD(NR,1) = ADDW(NR,1)+CNC*ADNC(NR,1)
            AD(NR,2) = ADDW(NR,2)+CNC*ADNC(NR,2)
            AD(NR,3) = ADDW(NR,3)+CNC*ADNC(NR,3)
            AD(NR,4) = ADDW(NR,4)+CNC*ADNC(NR,4)
C
            AVDW(NR,1) = -AV0*ADDW(NR,1)*RG(NR)/RA**2
            AVDW(NR,2) = -AV0*ADDW(NR,2)*RG(NR)/RA**2
            AVDW(NR,3) = -AV0*ADDW(NR,3)*RG(NR)/RA**2
            AVDW(NR,4) = -AV0*ADDW(NR,4)*RG(NR)/RA**2
C
            AV(NR,1) = AVDW(NR,1)+CNC*AVNC(NR,1)
            AV(NR,2) = AVDW(NR,2)+CNC*AVNC(NR,2)
            AV(NR,3) = AVDW(NR,3)+CNC*AVNC(NR,3)
            AV(NR,4) = AVDW(NR,4)+CNC*AVNC(NR,4)
  240    CONTINUE
      ELSE
         WRITE(6,*) 'XX INVALID MDLAD : ',MDLAD
         DO 250 NS=1,NSM
         DO 250 NR=1,NRMAX
            AD(NR,NS)=0.D0
            AV(NR,NS)=0.D0
  250    CONTINUE
      ENDIF
C
C        ****** AVK : HEAT PINCH ******
C
      IF(MDLAVK.EQ.0) THEN
         DO 300 NS=1,NSM
         DO 300 NR=1,NRMAX
            AVK(NR,NS)=0.D0
  300    CONTINUE
      ELSEIF(MDLAVK.EQ.1) THEN
         DO 400 NR=1,NRMAX
            AVK(NR,1)=-RG(NR)/RA*CHP
            AVK(NR,2)=-RG(NR)/RA*CHP
            AVK(NR,3)=-RG(NR)/RA*CHP
            AVK(NR,4)=-RG(NR)/RA*CHP
  400    CONTINUE
      ELSEIF(MDLAVK.EQ.2) THEN
         DO 500 NR=1,NRMAX
            AVKL=-RG(NR)/RA*(CHP*1.D6)/(ANE*1.D20*TE*RKEV)
            AVK(NR,1)=-AVKL
            AVK(NR,2)=-AVKL
            AVK(NR,3)=-AVKL
            AVK(NR,4)=-AVKL
  500    CONTINUE
      ELSE
         WRITE(6,*) 'XX INVALID MDLAVK : ',MDLAVK
         DO 600 NS=1,NSM
         DO 600 NR=1,NRMAX
            AVK(NR,NS)=0.D0
  600    CONTINUE
      ENDIF
C
      RETURN
      END
C
      FUNCTION TRCOFS(S,ALFA,RKCV)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      IF(ALFA.GE.0.D0) THEN
         SA=S-ALFA
         IF(SA.GE.0.D0) THEN
            FS1=(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0)
     &         /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA+2.0D0*SA*SA*SA))
         ELSE
            FS1=1.D0/SQRT(2.D0*(1.D0-2.D0*SA)
     &                       *(1.D0-2.D0*SA+3.D0*SA*SA))
         ENDIF
         IF(RKCV.GT.0.D0) THEN
            FS2=SQRT(RKCV)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ELSE
         SA=ALFA-S
         IF(SA.GE.0.D0) THEN
            FS1=(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0)
     &         /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA+2.0D0*SA*SA*SA))
         ELSE
            FS1=1.D0/SQRT(2.D0*(1.D0-2.D0*SA)
     &                    *(1.D0-2.D0*SA+3.D0*SA*SA))
         ENDIF
         IF(RKCV.LT.0.D0) THEN
            FS2=SQRT(-RKCV)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ENDIF
      TRCOFS=MAX(FS1,FS2)
      RETURN
      END
C
      FUNCTION TRCOFSX(S,ALFA,RKCV,EPSA)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      IF(ALFA.GE.0.D0) THEN
C        SA=S-ALFA
C        SA=S-(1.D0-(2*ALFA)/(1+3*(ALFA)**2))*ALFA
         SA=S-(1.D0-(2*ALFA)/(1+6*ALFA))*ALFA
         IF(SA.GE.0.D0) THEN
           FS1=((1.D0+RKCV)**2.5D0)*(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0)
     &         /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA*(1+RKCV)
     &                      +2.0D0*SA*SA*SA*(1.D0+RKCV)**2.5D0))
         ELSE
            FS1=((1.D0+RKCV)**2.5D0)/SQRT(2.D0*(1.D0-2.D0*SA)
     &                       *(1.D0-2.D0*SA+3.D0*SA*SA*(1+RKCV)))
         ENDIF
         IF(RKCV.GT.0.D0) THEN
            FS2=SQRT(RKCV/EPSA)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ELSE
C        SA=ALFA-S
C        SA=(1.D0-(2*ALFA)/(1+3*(ALFA)**2))*ALFA-S
         SA=(1.D0-(2*ALFA)/(1+6*ALFA))*ALFA-S
         IF(SA.GE.0.D0) THEN
            FS1=((1.D0+RKCV)**2.5D0)*(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0)
     &         /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA*(1+RKCV)
     &                      +2.0D0*SA*SA*SA*(1.D0+RKCV)**2.5D0))
         ELSE
            FS1=((1.D0+RKCV)**2.5D0)/SQRT(2.D0*(1.D0-2.D0*SA)
     &                       *(1.D0-2.D0*SA+3.D0*SA*SA*(1+RKCV)))
         ENDIF
         IF(RKCV.LT.0.D0) THEN
            FS2=SQRT(-RKCV/EPSA)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ENDIF
      TRCOFSX=MAX(FS1,FS2)
      RETURN
      END
C
      FUNCTION TRCOFT(S,ALFA,RKCV,EPSA)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      IF(ALFA.GE.0.D0) THEN
         SA=S-ALFA
         IF(SA.GE.0.D0) THEN
            FS1=(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0)
     &         /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA+2.0D0*SA*SA*SA))
         ELSE
            FS1=1.D0/SQRT(2.D0*(1.D0-2.D0*SA)
     &                       *(1.D0-2.D0*SA+3.D0*SA*SA))
         ENDIF
         IF(RKCV.GT.0.D0) THEN
            FS2=SQRT(RKCV/EPSA)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ELSE
         SA=ALFA-S
         IF(SA.GE.0.D0) THEN
            FS1=(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0)
     &         /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA+2.0D0*SA*SA*SA))
         ELSE
            FS1=1.D0/SQRT(2.D0*(1.D0-2.D0*SA)
     &                    *(1.D0-2.D0*SA+3.D0*SA*SA))
         ENDIF
         IF(RKCV.LT.0.D0) THEN
            FS2=SQRT(-RKCV/EPSA)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ENDIF
      TRCOFT=MAX(FS1,FS2)
      RETURN
      END
C
      FUNCTION RLAMBDA(X)
C
      REAL*8 RLAMBDA,X
      REAL*8 AX
      REAL*8 P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Y
      SAVE P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,
     &1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     &0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,0.2635537D-1,
     &-0.1647633D-1,0.392377D-2/
C
      AX=ABS(X)
      IF (AX.LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        RLAMBDA=(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
     &         *EXP(-AX)
      ELSE
        Y=3.75D0/AX
        RLAMBDA=(1.D0/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*
     &                           (Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
