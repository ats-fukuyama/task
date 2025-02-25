
C     $Id$
C     ***********************************************************
C
C           CALCULATE TRANSPORT COEFFICIENTS
C 
C     ***********************************************************
C
      SUBROUTINE TRCOEF
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
      INCLUDE 'trcomm.inc'
      DIMENSION S_HM(NRM)
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
         KGR1='/E$-r$=  vs r/'!'/NST$+2$= vs r/'
         KGR2='/V$-ExB$=  vs r/'!'/OmegaST vs r/'
C         KGR3='@1/(1+G*WE1$+2$=)  vs r@'!'@lambda vs r@'
         KGR3='@exp(-beta*WE1$+gamma$=)  vs r@'!'@lambda vs r@'
         KGR4='@Lambda,1/(1+OmgST$+2$=)  vs r@'
      ELSEIF(MDLKAI.LT.50) THEN
         KGR1='/NST$+2$= vs r/'
         KGR2='/OmegaST vs r/'
         KGR3='@lambda vs r@'
         KGR4='@Lambda,1/(1+OmgST$+2$=),1/(1+G*WE1$+2$=)vs r@'
      ELSEIF(MDLKAI.LT.60) THEN
         KGR1='/NST$+2$= vs r/'
         KGR2='/OmegaST vs r/'
         KGR3='@lambda vs r@'
         KGR4='@Lambda,1/(1+OmgST$+2$=),1/(1+G*WE1$+2$=)vs r@'
      ELSEIF(MDLKAI.LT.70) THEN
         KGR1='/V$-ExB$=/'
         KGR2='/E$-r$=/'
         KGR3='/GROWTH RATE/'
         KGR4='/ExB SHEARING RATE/'
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
      DO NR=1,NRMAX
         IF(SUMPBM.EQ.0.D0) THEN
            PADD(NR)=0.D0
         ELSE
            PADD(NR)=PBM(NR)*1.D-20/RKEV-RNF(NR,1)*RT(NR,2)
         ENDIF
C     Calculate ExB velocity in advance 
C                            for ExB shearing rate calculation
         VEXB(NR)= -ER(NR)/BB
      ENDDO
C
      DO NR=1,NRMAX
C     characteristic time of temporal change of transport coefficients
         TAUK(NR)=QP(NR)*RR/SQRT(RT(NR,2)*RKEV/(PA(2)*AMM))*DBLE(MDTC)
C         DRL=RJCB(NR)/DR
         DRL=1.D0/(DR*RA)
         EPS=EPSRHO(NR)
         RNTP=0.D0
         RNP =0.D0
         RNTM=0.D0
         RNM =0.D0
         IF(NR.EQ.NRMAX) THEN
            ANE =PNSS(1)
            ANDX=PNSS(2)
            ANT =PNSS(3)
            ANA =PNSS(4)
            TE=PTS(1)
            TD=PTS(2)
            TT=PTS(3)
            TA=PTS(4)
C
C     In the following, we assume that
C        1. pressures of beam and fusion at rho=1 are negligible,
C        2. calibration of Pbeam (from UFILE) is not necessary at rho=1.
C
            DO NS=2,NSM
               RNTP=RNTP+PNSS(NS)*PTS(NS)
               RNP =RNP +PNSS(NS)
               RNTM=RNTM+RN(NR-1,NS)*RT(NR-1,NS)
     &                  +RN(NR  ,NS)*RT(NR  ,NS)
               RNM =RNM +RN(NR-1,NS)+RN(NR  ,NS)
            ENDDO
            RNTM= RNTM+RW(NR-1,1)+RW(NR-1,2)
     &                +RW(NR  ,1)+RW(NR  ,2)
            RPP = RNTP+PNSS(1)   *PTS(1)
            RPM = RNTM+RN(NR-1,1)*RT(NR-1,1)
     &                +RN(NR  ,1)*RT(NR  ,1)+PADD(NR-1)+PADD(NR)
            RPEP= PNSS(1)   *PTS(1)
            RPEM= 0.5D0*(RN(NR-1,1)*RT(NR-1,1)
     &                  +RN(NR  ,1)*RT(NR  ,1))
            RNTM= 0.5D0*RNTM
            RNM = 0.5D0*RNM
            RPM = 0.5D0*RPM
C
            DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),
     &                  RHOG(NR),RHOM(NR),RHOM(NR-1))
            DNE=DERIV3P(PTS(1),RN(NR,1),RN(NR-1,1),
     &                  RHOG(NR),RHOM(NR),RHOM(NR-1))
            DPE=DERIV3P(PNSS(1)*PTS(1),
     &                  RN(NR,1)*RT(NR,1),
     &                  RN(NR-1,1)*RT(NR-1,1),
     &                  RHOG(NR),RHOM(NR),RHOM(NR-1))
            DTD=DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2),
     &                  RHOG(NR),RHOM(NR),RHOM(NR-1))
            DND=DERIV3P(PTS(2),RN(NR,2),RN(NR-1,2),
     &                  RHOG(NR),RHOM(NR),RHOM(NR-1))
            ZEFFL=ZEFF(NR)
            EZOHL=EZOH(NR)
C
            DPP = (RPP-RPM)*DRL
C
            TI  = RNTP/RNP
            DTI = (RNTP/RNP-RNTM/RNM)*DRL
         ELSE
C     density and temperature for each species on grid
            ANE    = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
            ANDX   = 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
            ANT    = 0.5D0*(RN(NR+1,3)+RN(NR  ,3))
            ANA    = 0.5D0*(RN(NR+1,4)+RN(NR  ,4))
            TE     = 0.5D0*(RT(NR+1,1)+RT(NR  ,1))
            TD     = 0.5D0*(RT(NR+1,2)+RT(NR  ,2))
            TT     = 0.5D0*(RT(NR+1,3)+RT(NR  ,3))
            TA     = 0.5D0*(RT(NR+1,4)+RT(NR  ,4))
C
C     incremental presssure and density for ion species
            DO NS=2,NSM
               RNTP=RNTP+RN(NR+1,NS)*RT(NR+1,NS)
               RNP =RNP +RN(NR+1,NS)
               RNTM=RNTM+RN(NR  ,NS)*RT(NR  ,NS)
               RNM =RNM +RN(NR  ,NS)
            ENDDO
C     incremental presssure and density for fast particles
            RNTP= RNTP+RW(NR+1,1)+RW(NR+1,2)
            RNTM= RNTM+RW(NR  ,1)+RW(NR  ,2)
C     incremental presssure and density for electron
            RPP = RNTP+RN(NR+1,1)*RT(NR+1,1)+PADD(NR+1)
            RPM = RNTM+RN(NR  ,1)*RT(NR  ,1)+PADD(NR  )
C     electron pressure
            RPEP= RN(NR+1,1)*RT(NR+1,1)
            RPEM= RN(NR  ,1)*RT(NR  ,1)
C
C     gradients of temperature, density and pressure for electron
            DTE = (RT(NR+1,1)-RT(NR  ,1))*DRL
            DNE = (RN(NR+1,1)-RN(NR  ,1))*DRL
            DPE = (RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL
C     gradients of temperature and density for ion
            DTD = (RT(NR+1,2)-RT(NR  ,2))*DRL
            DND = (RN(NR+1,2)-RN(NR  ,2))*DRL
C     effective charge number and parallel electric field on grid
            ZEFFL  = 0.5D0*(ZEFF(NR+1)+ZEFF(NR))
            EZOHL  = 0.5D0*(EZOH(NR+1)+EZOH(NR))
C
C     pressure gradient on grid
            DPP = (RPP-RPM)*DRL
C
C     effective ion temperature and its gradient on grid
            TI  = 0.5D0*(RNTP/RNP+RNTM/RNM)
            DTI = (RNTP/RNP-RNTM/RNM)*DRL
         ENDIF
C
c$$$C     second derivative of effective pressure for old version WE1
c$$$         RPI4=0.D0
c$$$         RPI3=0.D0
c$$$         RPI2=0.D0
c$$$         RPI1=0.D0
c$$$         IF(NR.LE.1) THEN
c$$$            DO NS=2,NSM
c$$$               RPI4=RPI4+RN(NR,  NS)*RT(NR,  NS)
c$$$            ENDDO
c$$$            RPI4=RPI4+RW(NR,  1)+RW(NR,  2)+PADD(NR  )
c$$$         ELSE
c$$$            DO NS=2,NSM
c$$$               RPI4=RPI4+RN(NR-1,NS)*RT(NR-1,NS)
c$$$            ENDDO
c$$$            RPI4=RPI4+RW(NR-1,1)+RW(NR-1,2)+PADD(NR-1)
c$$$         ENDIF
c$$$            DO NS=2,NSM
c$$$               RPI3=RPI3+RN(NR  ,NS)*RT(NR  ,NS)
c$$$            ENDDO
c$$$            RPI3=RPI3+RW(NR  ,1)+RW(NR  ,2)+PADD(NR  )
c$$$         IF(NR.GE.NRMAX-1) THEN
c$$$            DO NS=2,NSM
c$$$               RPI2=RPI2+RN(NR  ,NS)*RT(NR  ,NS)
c$$$            ENDDO
c$$$            RPI2=RPI2+RW(NR  ,1)+RW(NR  ,2)+PADD(NR  )
c$$$         ELSE
c$$$            DO NS=2,NSM
c$$$               RPI2=RPI2+RN(NR+1,NS)*RT(NR+1,NS)
c$$$            ENDDO
c$$$            RPI2=RPI2+RW(NR+1,1)+RW(NR+1,2)+PADD(NR+1)
c$$$         ENDIF
c$$$         IF(NR.GE.NRMAX-2) THEN
c$$$            DO NS=2,NSM
c$$$               RPI1=RPI1+RN(NR  ,NS)*RT(NR  ,NS)
c$$$            ENDDO
c$$$            RPI1=RPI1+RW(NR  ,1)+RW(NR  ,2)+PADD(NR  )
c$$$         ELSE
c$$$            DO NS=2,NSM
c$$$               RPI1=RPI1+RN(NR+1,NS)*RT(NR+1,NS)
c$$$            ENDDO
c$$$            RPI1=RPI1+RW(NR+1,1)+RW(NR+1,2)+PADD(NR+1)
c$$$         ENDIF
c$$$         RPIM=0.5D0*(RPI1+RPI2)
c$$$         RPI0=0.5D0*(RPI2+RPI3)
c$$$         RPIP=0.5D0*(RPI3+RPI4)
c$$$         DPPP=(RPIP-2*RPI0+RPIM)*DRL*DRL
C
C     safety factor and its gradient on grid
         DQ = DERIV3(NR,RHOG,QP,NRMAX,NRM,1)
         QL = QP(NR)
C
C     sound speed for electron
         VTE = SQRT(ABS(TE*RKEV/AME))
C
C     characteristic length of ion temperature
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
C     characteristic length of electron density
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
C     characteristic length of electron pressure
         IF(ABS(RPEP-RPEM).GT.1.D-32) THEN
            CLPE=0.5D0*(RPEP+RPEM)/(RPEP-RPEM)/DRL
         ELSE
            IF(RPEP-RPEM.GE.0.D0) THEN
               CLPE = 1.D32
            ELSE
               CLPE =-1.D32
            ENDIF
         ENDIF
C
C     ???
         IF(ABS(DQ).GT.1.D-32) THEN
            CLS=QL*QL/(DQ*EPS)
         ELSE
            IF(DQ.GE.0.D0) THEN
               CLS = 1.D32
            ELSE
               CLS =-1.D32
            ENDIF
         ENDIF
C
C     collision time between electrons and ions
         TAUE = FTAUE(ANE,ANDX,TE,ZEFFL)
         ANYUE = 0.5D0*(1.D0+ZEFFL)/TAUE
C
         ROUS = DSQRT(ABS(TE)*RKEV/AMD)/OMEGAD
         PPK  = 0.3D0/ROUS
         OMEGAS  = PPK*TE*RKEV/(AEE*BB*ABS(CLN))
         OMEGATT = DSQRT(2.D0)*VTE/(RR*QL)
C
C     Alfven wave velocity
         PNI=ANDX+ANT+ANA
         AMI=(AMD*ANDX+AMT*ANT+AMA*ANA)/PNI
         VA=SQRT(BB**2/(RMU0*ANE*1.D20*AMI))
C     magnetic shear
         S(NR)=RHOG(NR)/QL*DQ
C     pressure gradient for MHD instability
         ALPHA(NR)=-2.D0*RMU0*QL**2*RR/BB**2*(DPP*1.D20*RKEV)
C     magnetic curvature
         RKCV(NR)=-EPS*(1.D0-1.D0/(QL*QL))
C
C     rotational shear
C        omega(or gamma)_e=r/q d(q v_exb/r)/dr
         DVE = DERIV3(NR,RHOG,VEXB,NRMAX,NRM,1)
         WEXB(NR) = (S(NR)-1.D0)*VEXB(NR)/RHOG(NR)+DVE
C     Doppler shear
         AGMP(NR) = QP(NR)/EPS*WEXB(NR)
C
C   *************************************************************
C   ***  0.GE.MDLKAI.LT.10 : CONSTANT COEFFICIENT MODEL       ***
C   *** 10.GE.MDLKAI.LT.20 : DRIFT WAVE (+ITG +ETG) MODEL     ***
C   *** 20.GE.MDLKAI.LT.30 : REBU-LALLA MODEL                 ***
C   *** 30.GE.MDLKAI.LT.40 : CURRENT-DIFFUSIVITY DRIVEN MODEL ***
C   *** 40.GE.MDLKAI.LT.60 : DRIFT WAVE BALLOONING MODEL      ***
C   ***       MDLKAI.GE.60 : ITG(/TEM, ETG) MODEL ETC         ***
C   *************************************************************
C                                                                  
         IF(MDLKAI.LT.10) THEN
C   *********************************************************
C   ***  MDLKAI.EQ. 0   : CONSTANT                        ***
C   ***  MDLKAI.EQ. 1   : CONSTANT/(1-A*r**2)             ***
C   ***  MDLKAI.EQ. 2   : CONSTANT*(dTi/dr)**B/(1-A*r**2) ***
C   ***  MDLKAI.EQ. 3   : CONSTANT*(dTi/dr)**B*Ti**C      ***
C   ***  MDLKAI.EQ. 4   : PROP. TO CUBIC FUNC.            ***
C   *********************************************************
C
            IF(MDLKAI.EQ.0) THEN
               AKDWL=1.D0
            ELSEIF(MDLKAI.EQ.1) THEN
               AKDWL=1.D0/(1.D0-CKALFA*RHOG(NR)**2)
            ELSEIF(MDLKAI.EQ.2) THEN
               AKDWL=1.D0/(1.D0-CKALFA*RHOG(NR)**2)
     &                  *(ABS(DTI)*RA)**CKBETA
            ELSEIF(MDLKAI.EQ.3) THEN
               AKDWL=1.D0*(ABS(DTI)*RA)**CKBETA*ABS(TI)**CKGUMA
            ELSEIF(MDLKAI.EQ.4) THEN
               FSDFIX=0.1D0
               BPA=AR1RHOG(NRMAX)*RDPS/RR
               PROFDL=(PTS(1)*RKEV/(16.D0*AEE*SQRT(BB**2+BPA**2)))
     &               /FSDFIX
               AKDWL=FSDFIX*(1.D0+(PROFDL-1.D0)*(RHOG(NR)/ RA)**2)
            ELSE                                           
               WRITE(6,*) 'XX INVALID MDLKAI : ',MDLKAI
               AKDWL=0.D0
            ENDIF
            AKDW(NR,1)=CK0*AKDWL
            AKDW(NR,2)=CK1*AKDWL
            AKDW(NR,3)=CK1*AKDWL
            AKDW(NR,4)=CK1*AKDWL
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
            ELSEIF(MDLKAI.EQ.11.OR.
     &             MDLKAI.EQ.12.OR.
     &             MDLKAI.EQ.13) THEN
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
            ELSEIF(MDLKAI.EQ.15.OR.
     &             MDLKAI.EQ.16) THEN
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
               WRITE(6,*) 'XX INVALID MDLKAI : ',MDLKAI
               RRSTAR = RR 
               FDREV  = 1.D0
               HETA   = 1.D0
            ENDIF
            IF(MDLKAI.EQ.16) THEN
               DEDW = CK0*2.5D0*OMEGAS/PPK**2
     &               *(SQRT(EPS)*MIN(FDREV,EPS*OMEGAS/ANYUE)
     &                +OMEGAS/OMEGATT*MAX(1.D0,ANYUE/OMEGATT))
               RGL   = 2.5D0*OMEGAS*HETA
     &               *SQRT(2.D0*ABS(TI)*ABS(ETAI)*ABS(CLN)/(TE*RRSTAR))
               AMD=PA(2)*AMM
               TAUD = FTAUI(ANE,ANDX,TD,PZ(2),PA(2))
               RNUZ= 1.D0/(TAUD*SQRT(EPS))
               XCHI1=SQRT(RNUZ)/(SQRT(RNUZ)+SQRT(RGL))
               XXH=2.D0/(1.D0+PPK**2)
               RGLC=OMEGAS/(2*SQRT(2.D0)*XXH)
               XXA=XXH*PPK**2/4.D0
               IF(RGL.GT.RGLC) THEN
                  XCHI2=XXA*(SQRT(1.D0+(2.D0/XXA)*(RGL-RGLC)/RGL)
     &                       -1.D0)
               ELSE
                  XCHI2=0.D0
               ENDIF
               XCHI0=SQRT(XCHI1**2+XCHI2**2)
               DIDW = CK1*XCHI0*RGL/PPK**2
            ELSE
            DEDW = CK0*2.5D0*OMEGAS/PPK**2
     &            *(SQRT(EPS)*MIN(FDREV,EPS*OMEGAS/ANYUE)
     &             +OMEGAS/OMEGATT*MAX(1.D0,ANYUE/OMEGATT))
C
            DIDW = CK1*2.5D0*OMEGAS/PPK**2*HETA
     &            *SQRT(2.D0*ABS(TI)*ABS(ETAI)*ABS(CLN)/(TE*RRSTAR))
            ENDIF
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
     &              *SQRT(AEE**2/(RMU0*SQRT(AME)))/(QL*RKEV)
               IF(ABS(DTE).LE.CRTCL.OR.
     &            DQ.LE.0.D0.OR.
     &            NR.EQ.1) THEN
                  DKAI=0.D0
               ELSE
                  DKAI=1.5D-1*(ABS(DTE)/TE
     &                 +2.D0*ABS(DNE)/ANE)
     &                 *SQRT(TE/TI)*(1.D0/EPS)
     &                 *(QL**2/(DQ*BB*SQRT(RR)))*VC**2
     &                 *SQRT(RMU0*1.5D0*AMM)
     &                 *(1.D0-CRTCL/ABS(DTE))
               ENDIF
               AKDW(NR,1)=                  DKAI
               AKDW(NR,2)=ZEFFL*SQRT(TE/TD)*DKAI
               AKDW(NR,3)=ZEFFL*SQRT(TE/TT)*DKAI
               AKDW(NR,4)=ZEFFL*SQRT(TE/TA)*DKAI
            ELSE                                           
               WRITE(6,*) 'XX INVALID MDLKAI : ',MDLKAI
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
         ELSEIF(MDLKAI.LT.41) THEN
C
            WPE2=ANE*1.D20*AEE*AEE/(AME*EPS0)
            DELTA2=VC**2/WPE2
C
            RNST2=0.D0
            OMEGASS=0.D0
            SLAMDA=0.D0
            RLAMDA=0.D0
            RG1=1.D0
C
            IF(MOD(MDLKAI,2).EQ.0) THEN
               SL=SQRT(S(NR)**2+0.1D0**2)
               WE1=-QL*RR/(SL*VA)*DVE
               RG1=CWEB*FEXB(ABS(WE1),S(NR),ALPHA(NR))
C               DBDRR=DPPP*1.D20*RKEV*RA*RA/(BB**2/(2*RMU0))
C               DELTAE=SQRT(DELTA2)
C               WE1=SQRT(PA(2)/PA(1))*(QL*RR*DELTAE)/(2*SL*RA*RA)*DBDRR
            ENDIF
C
            IF(MDLKAI.EQ.30) THEN
               FS=1.D0/(1.7D0+SQRT(6.D0)*S(NR))
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            ELSEIF(MDLKAI.EQ.31) THEN
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFS(S(NR),ALPHAL,RKCV(NR))
               IF(MDCD05.NE.0) 
     &         FS=FS*(2.D0*SQRT(RKPRHO(NR))/(1.D0+RKPRHO(NR)**2))**1.5D0
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            ELSEIF(MDLKAI.EQ.32) THEN
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFS(S(NR),ALPHAL,RKCV(NR))
C               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               IF(MDCD05.NE.0) 
     &         FS=FS*(2.D0*SQRT(RKPRHO(NR))/(1.D0+RKPRHO(NR)**2))**1.5D0
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            ELSEIF(MDLKAI.EQ.33) THEN
               FS=TRCOFS(S(NR),0.D0,RKCV(NR))
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            ELSEIF(MDLKAI.EQ.34) THEN
               FS=TRCOFS(S(NR),0.D0,RKCV(NR))
C               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            ELSEIF(MDLKAI.EQ.35) THEN
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFSS(S(NR),ALPHAL)
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            ELSEIF(MDLKAI.EQ.36) THEN
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFSS(S(NR),ALPHAL)
C               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            ELSEIF(MDLKAI.EQ.37) THEN
               FS=TRCOFSS(S(NR),0.D0)
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            ELSEIF(MDLKAI.EQ.38) THEN
               FS=TRCOFSS(S(NR),0.D0)
C               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            ELSEIF(MDLKAI.EQ.39) THEN
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFSX(S(NR),ALPHAL,RKCV(NR),RA/RR)
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
C
c$$$            ELSEIF(MDLKAI.EQ.40) THEN
c$$$               AEI=(PZ(2)*ANDX+PZ(3)*ANT+PZ(4)*ANA)*AEE/PNI
c$$$               WCI=AEI*BB/AMI
c$$$               PTI=(TD*ANDX+TT*ANT+TA*ANA)/PNI
c$$$               VTI=SQRT(ABS(PTI*RKEV/AMI))
c$$$               RHOI=VTI/WCI
c$$$C
c$$$               FS=TRCOFT(S(NR),ALPHA(NR),RKCV(NR),RA/RR)
c$$$               SA=S(NR)-ALPHA(NR)
c$$$               RNST2=0.5D0/((1.D0-2.D0*SA+3.D0*SA*SA)*FS)
c$$$               RKPP2=RNST2/(FS*ABS(ALPHA(NR))*DELTA2)
c$$$C
c$$$               SLAMDA=RKPP2*RHOI**2
c$$$               RLAMDA=RLAMBDA(SLAMDA)
c$$$               OMEGAS= SQRT(RKPP2)*TE*RKEV/(AEE*BB*ABS(CLPE))
c$$$               TAUAP=(QL*RR)/VA
c$$$               OMEGASS=(OMEGAS*TAUAP)/(RNST2*SQRT(ALPHA(NR)))
c$$$C
c$$$c$$$               FS=FS/(1.D0+RG1*WE1*WE1)
c$$$C
c$$$               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
c$$$     &               /(RLAMDA*(1.D0+OMEGASS**2))
c$$$               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
c$$$     &               /(1.D0+OMEGASS**2)
            ELSE                                           
               WRITE(6,*) 'XX INVALID MDLKAI : ',MDLKAI
               FS=1.D0
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            ENDIF
            AKDW(NR,1)=AKDWEL
            AKDW(NR,2)=AKDWIL
            AKDW(NR,3)=AKDWIL
            AKDW(NR,4)=AKDWIL
C
            VGR1(NR,1)=FS
            VGR1(NR,2)=S(NR)
            VGR1(NR,3)=ALPHA(NR)
            VGR2(NR,1)=ER(NR)!RNST2
            VGR2(NR,2)=VEXB(NR)!OMEGASS
            VGR2(NR,3)=0.D0
C            VGR3(NR,1)=1.D0/(1.D0+RG1*WE1*WE1)!SLAMDA
            VGR3(NR,1)=RG1!SLAMDA
            VGR3(NR,2)=ABS(WE1)
            VGR3(NR,3)=0.D0
            VGR4(NR,1)=RLAMDA
            VGR4(NR,2)=1.D0/(1.D0+OMEGASS**2)
            VGR4(NR,3)=0.D0
C
         ELSEIF(MDLKAI.LT.51) THEN
C
            WPE2=ANE*1.D20*AEE*AEE/(AME*EPS0)
            DELTA2=VC**2/WPE2
C
            RNST2=0.D0
            OMEGASS=0.D0
            SLAMDA=0.D0
            RLAMDA=0.D0
C
            IF(MOD(MDLKAI,2).EQ.0) THEN
               RG1=CWEB*FEXB(ABS(WE1),S(NR),ALPHA(NR))
               SL=SQRT(S(NR)**2+0.1D0**2)
               WE1=-QL*RR/(SL*VA)*DVE
C               DBDRR=DPPP*1.D20*RKEV*RA*RA/(BB**2/(2*RMU0))
C               DELTAE=SQRT(DELTA2)
C               WE1=SQRT(PA(2)/PA(1))*(QL*RR*DELTAE)/(2*SL*RA*RA)*DBDRR
            ENDIF
C
            F=VTE/VA
C     
            IF(MDLKAI.EQ.40) THEN
               FS=1.D0/(1.7D0+SQRT(6.D0)*S(NR))
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            ELSEIF(MDLKAI.EQ.41) THEN
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFS(S(NR),ALPHAL,RKCV(NR))
C               IF (NR.LE.2) write(6,'(I5,4F15.10)') NR,S(NR),ALPHA(NR),RKCV(NR),FS
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            ELSEIF(MDLKAI.EQ.42) THEN
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFS(S(NR),ALPHAL,RKCV(NR))
C               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            ELSEIF(MDLKAI.EQ.43) THEN
               FS=TRCOFS(S(NR),0.D0,RKCV(NR))
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            ELSEIF(MDLKAI.EQ.44) THEN
               FS=TRCOFS(S(NR),0.D0,RKCV(NR))
C               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            ELSEIF(MDLKAI.EQ.45) THEN
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFSS(S(NR),ALPHAL)
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            ELSEIF(MDLKAI.EQ.46) THEN
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFSS(S(NR),ALPHAL)
C               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            ELSEIF(MDLKAI.EQ.47) THEN
               FS=TRCOFSS(S(NR),0.D0)
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            ELSEIF(MDLKAI.EQ.48) THEN
               FS=TRCOFSS(S(NR),0.D0)
C               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            ELSEIF(MDLKAI.EQ.49) THEN
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFSX(S(NR),ALPHAL,RKCV(NR),RA/RR)
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
C
            ELSEIF(MDLKAI.EQ.50) THEN
               AEI=(PZ(2)*ANDX+PZ(3)*ANT+PZ(4)*ANA)*AEE/PNI
               WCI=AEI*BB/AMI
               PTI=(TD*ANDX+TT*ANT+TA*ANA)/PNI
               VTI=SQRT(ABS(PTI*RKEV/AMI))
               RHOI=VTI/WCI
C
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFT(S(NR),ALPHAL,RKCV(NR),RA/RR)
               SA=S(NR)-ALPHA(NR)
               RNST2=0.5D0/((1.D0-2.D0*SA+3.D0*SA*SA)*FS)
               RKPP2=RNST2/(FS*ABS(ALPHA(NR))*DELTA2)
C
               SLAMDA=RKPP2*RHOI**2
               RLAMDA=RLAMBDA(SLAMDA)
               OMEGAS= SQRT(RKPP2)*TE*RKEV/(AEE*BB*ABS(CLPE))
               TAUAP=(QL*RR)/VA
               OMEGASS=(OMEGAS*TAUAP)/(RNST2*SQRT(ALPHA(NR)))
               AKDWEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)
     &               /(RLAMDA*(1.D0+OMEGASS**2))
               AKDWIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)
     &               /(1.D0+OMEGASS**2)
            ENDIF
C
            AKDW(NR,1)=AKDWEL
            AKDW(NR,2)=AKDWIL
            AKDW(NR,3)=AKDWIL
            AKDW(NR,4)=AKDWIL
C
            VGR1(NR,1)=FS
            VGR1(NR,2)=S(NR)
            VGR1(NR,3)=ALPHA(NR)
            VGR2(NR,1)=RNST2
            VGR2(NR,2)=OMEGASS
            VGR2(NR,3)=0.D0
            VGR3(NR,1)=SLAMDA
            VGR3(NR,2)=0.D0
            VGR3(NR,3)=0.D0
            VGR4(NR,1)=RLAMDA
            VGR4(NR,2)=1.D0/(1.D0+OMEGASS**2)
C            VGR4(NR,3)=1.D0/(1.D0+RG1*WE1*WE1)
            VGR4(NR,3)=RG1
C
         ELSEIF(MDLKAI.GE.60) THEN
C
            WPE2=ANE*1.D20*AEE*AEE/(AME*EPS0)
            DELTA2=VC**2/WPE2
C
            RNST2=0.D0
            OMEGASS=0.D0
            SLAMDA=0.D0
            RLAMDA=0.D0
C
            IF(NR.EQ.1) THEN
               DRL=RJCB(NR)/DR
               S_HM(NR) = RM(NR)*RA/(0.5D0*(QP(NR)+Q0))
     &                   *(QP(NR)-Q0)*DRL
            ELSE
               DRL=RJCB(NR)/DR
               S_HM(NR) = RM(NR)/(0.5D0*(QP(NR)+QP(NR-1)))
     &                   *(QP(NR)-QP(NR-1))/DR
            ENDIF
C
            VGR1(NR,2)=S(NR)
            VGR1(NR,3)=ALPHA(NR)
            VGR2(NR,1)=VEXB(NR)
            VGR2(NR,2)=ER(NR)
            VGR2(NR,3)=0.D0
            VGR3(NR,1)=AGMP(NR)
            VGR3(NR,2)=0.D0
            VGR3(NR,3)=0.D0
            VGR4(NR,1)=WEXB(NR)
            VGR4(NR,2)=WEXB(NR)
            VGR4(NR,3)=0.D0
            IF(MDLKAI.EQ.65) THEN
               DPERHO=DPE/RJCB(NR)
               DTERHO=DTE/RJCB(NR)
               NR08=INT(0.8D0*NRMAX)
               CHIB  = (ABS(DPERHO*1.D3)/(ANE*BB))*QL*QL
     &                *((RT(NR08,1)-RT(NRMAX,1))/RT(NRMAX,1))
               RHOS  = 1.02D-4*SQRT(PA(2)*TE*1.D3/PZ(2))/(RA*BB)
C               RHOS  = SQRT(2.D0*AMM/AEE)*SQRT(PA(2)*TE*1.D3)
C     &                /(PZ(2)*RA*BB)
               CHIGB = RHOS*ABS(DTERHO*1.D3)/BB
               CS    = SQRT(ABS(TE*RKEV/(PA(2)*AMM)))
               ALNI  = ABS(DND/ANDX)
               ALTI  = ABS(DTD/TD)
               AGITG = 0.1D0*CS/RA*SQRT(RA*ALNI+RA*ALTI)*SQRT(TD/TE)
               WEXB(NR)=0.D0
               AKDW(NR,1) = 8.D-5  *CHIB*FBHM(WEXB(NR),AGITG,S(NR))
     &                     +7.D-2  *CHIGB
               AKDW(NR,2) = 1.6D-4 *CHIB*FBHM(WEXB(NR),AGITG,S(NR))
     &                     +1.75D-2*CHIGB
               AKDW(NR,3) = AKDW(NR,2)
               AKDW(NR,4) = AKDW(NR,2)
               VGR1(NR,1) = FBHM(WEXB(NR),AGITG,S(NR))
            ENDIF
         ELSE                                           
            WRITE(6,*) 'XX INVALID MDLKAI : ',MDLKAI
            AKDW(NR,1)=0.D0
            AKDW(NR,2)=0.D0
            AKDW(NR,3)=0.D0
            AKDW(NR,4)=0.D0
         ENDIF
      ENDDO
C      STOP
C
      IF(MDLKAI.EQ.60.OR.MDLKAI.EQ.61) THEN
         CALL GLF23_DRIVER(S_HM)
      ELSEIF(MDLKAI.EQ.62) THEN
         CALL AITKEN(1.D0,RBEEDG,RM,RNF(1,1),2,NRMAX)
         RBEEDG=RBEEDG/PNSS(1)
         CALL IFSPPPL_DRIVER(NRM,NSM,NSTM,NRMAX,RN,RR,DR,RJCB,RHOG,RHOM,
     &                       QP,S,EPSRHO,RKPRHOG,RT,BB,AMM,AME,
     &                       PNSS,PTS,RNF(1,1),RBEEDG,MDLUF,NSMAX,
     &                       AR1RHOG,AR2RHOG,AKDW)
      ELSEIF(MDLKAI.EQ.63.OR.MDLKAI.EQ.64) THEN
         CALL WEILAND_DRIVER
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
      SUBROUTINE TRCFNC
C
      INCLUDE 'trcomm.inc'
C
C      DATA RK11,RA11,RB11,RC11/1.04D0,2.01D0,1.53D0,0.89D0/
C      DATA RK12,RA12,RB12,RC12/1.20D0,0.76D0,0.67D0,0.56D0/
      DATA RK22,RA22,RB22,RC22/2.55D0,0.45D0,0.43D0,0.43D0/
C      DATA RK13,RA13,RB13,RC13/2.30D0,1.02D0,1.07D0,1.07D0/
C      DATA RK23,RA23,RB23,RC23/4.19D0,0.57D0,0.61D0,0.61D0/
C      DATA RK33,RA33,RB33,RC33/1.83D0,0.68D0,0.32D0,0.66D0/
      DATA RK2 ,RA2 ,RB2 ,RC2 /0.66D0,1.03D0,0.31D0,0.74D0/
C
      SAVE CDHSV
C
      AMD=PA(2)*AMM
      AMT=PA(3)*AMM
      AMA=PA(4)*AMM
C
      DO NR=1,NRMAX
         IF(NR.EQ.NRMAX) THEN
            ANE =PNSS(1)
            ANDX=PNSS(2)
            ANT =PNSS(3)
            ANA =PNSS(4)
            TE=PTS(1)
            TD=PTS(2)
            TT=PTS(3)
            TA=PTS(4)
            ZEFFL=ZEFF(NR)
         ELSE
            ANE    = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
            ANDX   = 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
            ANT    = 0.5D0*(RN(NR+1,3)+RN(NR  ,3))
            ANA    = 0.5D0*(RN(NR+1,4)+RN(NR  ,4))
            TE     = 0.5D0*(RT(NR+1,1)+RT(NR  ,1))
            TD     = 0.5D0*(RT(NR+1,2)+RT(NR  ,2))
            TT     = 0.5D0*(RT(NR+1,3)+RT(NR  ,3))
            TA     = 0.5D0*(RT(NR+1,4)+RT(NR  ,4))
            ZEFFL  = 0.5D0*(ZEFF(NR+1)+ZEFF(NR))
         ENDIF
C
         QL = QP(NR)
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
C
         VTE = SQRT(ABS(TE*RKEV/AME))
         VTD = SQRT(ABS(TD*RKEV/AMD))
         VTT = SQRT(ABS(TT*RKEV/AMT))
         VTA = SQRT(ABS(TA*RKEV/AMA))
C
         RHOE2=2.D0*AME*ABS(TE)*RKEV/(PZ(1)*AEE*BP(NR))**2
         RHOD2=2.D0*AMD*ABS(TD)*RKEV/(PZ(2)*AEE*BP(NR))**2
         RHOT2=2.D0*AMT*ABS(TT)*RKEV/(PZ(3)*AEE*BP(NR))**2
         RHOA2=2.D0*AMA*ABS(TA)*RKEV/(PZ(4)*AEE*BP(NR))**2
C
c$$$         TAUE = FTAUE(ANE,ANDX,TE,1.D0)
c$$$         TAUD = FTAUI(ANE,ANDX,TD,1.D0,PA(2))
c$$$         TAUT = FTAUI(ANE,ANT ,TT,1.D0,PA(3))
c$$$         TAUA = FTAUI(ANE,ANA ,TA,2.D0,PA(4))
         TAUE = FTAUE(ANE,ANDX,TE,ZEFFL)
         TAUD = FTAUI(ANE,ANDX,TD,PZ(2),PA(2))
         TAUT = FTAUI(ANE,ANT ,TT,PZ(3),PA(3))
         TAUA = FTAUI(ANE,ANA ,TA,PZ(4),PA(4))
C
         RNUE=ABS(QL)*RR/(TAUE*VTE*EPSS)
         RNUD=ABS(QL)*RR/(TAUD*VTD*EPSS)
         RNUT=ABS(QL)*RR/(TAUT*VTT*EPSS)
         RNUA=ABS(QL)*RR/(TAUA*VTA*EPSS)
C
C     ***** NEOCLASSICAL TRANSPORT (HINTON, HAZELTINE) *****
C
      IF(MDLKNC.EQ.1) THEN
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
C     ***** CHANG HINTON *****
C
         DELDA=0.D0
C
         IF(MDLUF.EQ.0) THEN
            RALPHA=ZEFFL-1.D0
         ELSE
            RALPHA=PZ(3)**2*ANT/(PZ(2)**2*ANDX)
         ENDIF
C
         RMUSD=RNUD*(1.D0+1.54D0*RALPHA)
         RMUST=RNUT*(1.D0+1.54D0*RALPHA)
         RMUSA=RNUA*(1.D0+1.54D0*RALPHA)
C
         F1=(1.D0+1.5D0*(EPS**2+EPS*DELDA)
     &      +3.D0/8.D0*EPS**3*DELDA)
     &      /(1.D0+0.5D0*EPS*DELDA)
         F2=DSQRT(1.D0-EPS**2)*(1.D0+0.5D0*EPS*DELDA)
     &         /(1.D0+DELDA/EPS*(DSQRT(1.D0-EPS**2)-1.D0))
C
         TERM1D=(0.66D0*(1.D0+1.54D0*RALPHA)
     &       +(1.88D0*DSQRT(EPS)-1.54D0*EPS)*(1.D0+3.75D0*RALPHA))*F1
     &       /(1.D0+1.03D0*DSQRT(RMUSD)+0.31D0*RMUSD)
         TERM1T=(0.66D0*(1.D0+1.54D0*RALPHA)
     &       +(1.88D0*DSQRT(EPS)-1.54D0*EPS)*(1.D0+3.75D0*RALPHA))*F1
     &       /(1.D0+1.03D0*DSQRT(RMUST)+0.31D0*RMUST)
         TERM1A=(0.66D0*(1.D0+1.54D0*RALPHA)
     &       +(1.88D0*DSQRT(EPS)-1.54D0*EPS)*(1.D0+3.75D0*RALPHA))*F1
     &       /(1.D0+1.03D0*DSQRT(RMUSA)+0.31D0*RMUSA)
C
         TERM2D=0.583D0*RMUSD*EPS/(1.D0+0.74D0*RMUSD*EPSS)
     &       *(1.D0+(1.33D0*RALPHA*(1.D0+0.6D0*RALPHA))
     &       /(1.D0+1.79D0*RALPHA))
     &       *(F1-F2)
         TERM2T=0.583D0*RMUST*EPS/(1.D0+0.74D0*RMUST*EPSS)
     &       *(1.D0+(1.33D0*RALPHA*(1.D0+0.6D0*RALPHA))
     &       /(1.D0+1.79D0*RALPHA))
     &       *(F1-F2)
         TERM2A=0.583D0*RMUSA*EPS/(1.D0+0.74D0*RMUSA*EPSS)
     &       *(1.D0+(1.33D0*RALPHA*(1.D0+0.6D0*RALPHA))
     &       /(1.D0+1.79D0*RALPHA))
     &       *(F1-F2)
C
         AKNC(NR,1)=0.D0
C         AKNC(NR,2)=(QL**2*RHOD2)/(EPSS*TAUD)*(TERM1D+TERM2D)
C         AKNC(NR,3)=(QL**2*RHOT2)/(EPSS*TAUT)*(TERM1T+TERM2T)
C         AKNC(NR,4)=(QL**2*RHOA2)/(EPSS*TAUA)*(TERM1A+TERM2A)
         AKNC(NR,2)=(RHOD2*SQRT(EPS))/TAUD*(TERM1D+TERM2D)
         AKNC(NR,3)=(RHOT2*SQRT(EPS))/TAUT*(TERM1T+TERM2T)
         AKNC(NR,4)=(RHOA2*SQRT(EPS))/TAUA*(TERM1A+TERM2A)
C
      ENDIF
C
C     Limit of neoclassical diffusivity
         DO NS=1,4
            CHECK=ABS(RT(NR,NS)*RKEV/(2.D0*PZ(NS)*AEE*RR*BB))*RG(NR)*RA
            IF(AKNC(NR,NS).GT.CHECK) AKNC(NR,NS)=CHECK
         ENDDO
      ENDDO
C
      ENTRY TRCFDW_AKDW
C
      IF(MDEDGE.EQ.1) CDHSV=CDH
      DO NR=1,NRMAX
         IF(MDEDGE.EQ.1.AND.NR.GE.NREDGE) CDH=CSPRS
         DO NS=1,NSM
            AKDW(NR,NS) = CDH*AKDW(NR,NS)
            AK(NR,NS) = AKDW(NR,NS)+CNH*AKNC(NR,NS)
         ENDDO
      ENDDO
      IF(MDEDGE.EQ.1) CDH=CDHSV
C
C     ***** OFF-DIAGONAL TRANSPORT COEFFICIENTS *****
C
C     AKLP : heat flux coefficient for pressure gradient
C     AKLD : heat flux coefficient for density gradient
C
      IF(MDDIAG.EQ.1) THEN
         DO NR=1,NRMAX
            DO NS=1,NSLMAX
               DO NS1=1,NSLMAX
                  IF(NS.EQ.NS1) THEN
                     AKLP(NR,NS,NS1)= AKDW(NR,NS)
     &                               +CNH*( AKNCT(NR,NS,NS1)
     &                                     +AKNCP(NR,NS,NS1))
                     AKLD(NR,NS,NS1)=-AKDW(NR,NS)
     &                               +CNH*(-AKNCT(NR,NS,NS1))
                  ELSE
                     AKLP(NR,NS,NS1)= CNH*( AKNCT(NR,NS,NS1)
     &                                     +AKNCP(NR,NS,NS1))
                     AKLD(NR,NS,NS1)= CNH*(-AKNCT(NR,NS,NS1))
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ELSEIF(MDDIAG.EQ.2) THEN
         DO NR=1,NRMAX
            IF(MDEDGE.EQ.1.AND.NR.GE.NREDGE) CDH=CSPRS
            DO NS=1,NSLMAX
               DO NS1=1,NSLMAX
                  IF(NS.EQ.NS1) THEN
                     AKLP(NR,NS,NS1)= CDH*AKDWP(NR,NS,NS1)
     &                               +CNH*AKNC(NR,NS)
                     AKLD(NR,NS,NS1)= CDH*AKDWD(NR,NS,NS1)
                  ELSE
                     AKLP(NR,NS,NS1)= CDH*AKDWP(NR,NS,NS1)
                     AKLD(NR,NS,NS1)= CDH*AKDWD(NR,NS,NS1)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ELSEIF(MDDIAG.EQ.3) THEN
         DO NR=1,NRMAX
            IF(MDEDGE.EQ.1.AND.NR.GE.NREDGE) CDH=CSPRS
            DO NS=1,NSLMAX
               DO NS1=1,NSLMAX
                  AKLP(NR,NS,NS1)= CDH*  AKDWP(NR,NS,NS1)
     &                            +CNH*( AKNCT(NR,NS,NS1)
     &                                  +AKNCP(NR,NS,NS1))
                  AKLD(NR,NS,NS1)= CDH*  AKDWD(NR,NS,NS1)
     &                            +CNH*(-AKNCT(NR,NS,NS1))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      IF(MDEDGE.EQ.1) CDH=CDHSV
C     
      RETURN
      END
C
C     ***********************************************************
C
      SUBROUTINE TRCFET
C
      INCLUDE 'trcomm.inc'
C
      DATA RK33,RA33,RB33,RC33/1.83D0,0.68D0,0.32D0,0.66D0/
C
      DO NR=1,NRMAX
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
C
C        ****** CLASSICAL RESISTIVITY (Spitzer) from JAERI Report ******
C
         ANE=RN(NR,1)
         ANI=RN(NR,2)
         TE =RT(NR,1)
         ZEFFL=ZEFF(NR)
         TAUE = FTAUE(ANE,ANI,TE,ZEFFL)
C
         ETA(NR) = AME/(ANE*1.D20*AEE**2*TAUE)
     &        *(0.29D0+0.46D0/(1.08D0+ZEFFL))
C
C        ****** NEOCLASSICAL RESISTIVITY (Hinton, Hazeltine) ******
C
         IF(MDLETA.EQ.1) THEN
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
            H      = BB/(BB+BP(NR))
            FT     = 1.D0/H-SQRT(EPS)*RK33E
            ETA(NR)= ETA(NR)/FT
C
C        ****** NEOCLASSICAL RESISTIVITY (Hirshman, Hawryluk) ******
C
         ELSEIF(MDLETA.EQ.2) THEN
            IF(NR.EQ.1) THEN
               QL=ABS(0.25D0*(3.D0*Q0+QP(NR)))
            ELSE
               QL=ABS(0.5D0*(QP(NR-1)+QP(NR)))
            ENDIF
            ZEFFL=ZEFF(NR)
            VTE=SQRT(ABS(TE)*RKEV/AME)
            FT=FTPF(MDLTPF,EPS)
            TAUEL=6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)
     &           *(TE*RKEV)**1.5D0/(ANE*1.D20*AEE**4*COULOG(1,2,ANE,TE))
            RNUSE=RR*QL/(VTE*TAUEL*EPSS)
            PHI=FT/(1.D0+(0.58D0+0.20D0*ZEFFL)*RNUSE)
            ETAS=1.65D-9*COULOG(1,2,ANE,TE)/(ABS(TE)**1.5D0)
            CH=0.56D0*(3.D0-ZEFFL)/((3.D0+ZEFFL)*ZEFFL)
C
            ETA(NR)=ETAS*ZEFFL*(1.D0+0.27D0*(ZEFFL-1.D0))
     &           /((1.D0-PHI)*(1.D0-CH*PHI)*(1.D0+0.47D0*(ZEFFL-1.D0)))
C
C        ****** NEOCLASSICAL RESISTIVITY (Sauter)  ******
C
         ELSEIF(MDLETA.EQ.3) THEN
            IF(NR.EQ.1) THEN
               QL= 0.25D0*(3.D0*Q0+QP(NR))
            ELSE
               QL= 0.5D0*(QP(NR-1)+QP(NR))
            ENDIF
            ZEFFL=ZEFF(NR)
            rLnLame=31.3D0-LOG(SQRT(ANE*1.D20)/ABS(TE*1.D3))
            RNZ=0.58D0+0.74D0/(0.76D0+ZEFFL)
            SGMSPTZ=1.9012D4*(TE*1.D3)**1.5D0/(ZEFFL*RNZ*rLnLame)
            FT=FTPF(MDLTPF,EPS)
            RNUE=6.921D-18*ABS(QL)*RR*ANE*1.D20*ZEFFL*rLnLame
     &          /((TE*1.D3)**2*EPSS)
            F33TEFF=FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE)
     &             +0.45D0*(1.D0-FT)*RNUE/ZEFFL**1.5D0)
            ETA(NR)=1.D0/(SGMSPTZ*F33(F33TEFF,ZEFFL))
C
C        ****** NEOCLASSICAL RESISTIVITY (Hirshman, Sigmar)  ******
C
         ELSEIF(MDLETA.EQ.4) THEN
            IF(NR.EQ.1) THEN
               QL= 0.25D0*(3.D0*Q0+QP(NR))
            ELSE
               QL= 0.5D0*(QP(NR-1)+QP(NR))
            ENDIF
            ZEFFL=ZEFF(NR)
            ANE  =RN(NR,1)
            TEL  =ABS(RT(NR,1))
C
C     p1157 (7.36)
            ETAS = AME/(ANE*1.D20*AEE*AEE*TAUE)
     &           *( (1.D0+1.198D0*ZEFFL+0.222D0*ZEFFL**2)
     &             /(1.D0+2.966D0*ZEFFL+0.753D0*ZEFFL**2))
C
            FT=FTPF(MDLTPF,EPS)
            XI=0.58D0+0.2D0*ZEFFL
            CR=0.56D0/ZEFFL*(3.D0-ZEFFL)/(3.D0+ZEFFL)
C     RNUE expressions is given by the paper by Hirshman, Hawryluk.
            RNUE=SQRT(2.D0)/EPSS*RR*QL/SQRT(2.D0*TEL*RKEV/AME)/TAUE
C     p1158 (7.41)
            ETA(NR)=ETAS/(1.D0-FT/(1.D0+XI*RNUE))
     &                     /(1.D0-CR*FT/(1.D0+XI*RNUE))
         ENDIF
      ENDDO
C
C        ****** NEOCLASSICAL RESISTIVITY BY NCLASS  ******
C
      IF(NT.NE.0.AND.MDNCLS.NE.0) THEN
         DO NR=1,NRMAX
            ETA(NR)=ETANC(NR)
         ENDDO
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
      SUBROUTINE TRCFAD
C
      INCLUDE 'trcomm.inc'
      DIMENSION ACOEF(2),BCOEF(5)
C
C     ZEFF=1
C
      DATA RK11,RA11,RB11,RC11/1.04D0,2.01D0,1.53D0,0.89D0/
      DATA RK13,RA13,RB13,RC13/2.30D0,1.02D0,1.07D0,1.07D0/
      DATA RK23,RA23,RB23,RC23/4.19D0,0.57D0,0.61D0,0.61D0/
C
      SAVE CDPSV
C
C        ****** AD : PARTICLE DIFFUSION ******
C        ****** AV : PARTICLE PINCH ******
C
      IF(MDEDGE.EQ.1) CDPSV=CDP
      IF(MDNCLS.EQ.0) THEN
C     NCLASS has already calculated neoclassical particle pinch(AVNC)
C     beforehand if MDNCLS=1 so that MDLAD becomes no longer valid.
      IF(MDLAD.EQ.1) THEN
         DO NR=1,NRMAX
            IF(NR.EQ.NRMAX) THEN
               ANE =PNSS(1)
               ANDX=PNSS(2)
               ANT =PNSS(3)
               ANA =PNSS(4)
            ELSE
               ANE    = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
               ANDX   = 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
               ANT    = 0.5D0*(RN(NR+1,3)+RN(NR  ,3))
               ANA    = 0.5D0*(RN(NR+1,4)+RN(NR  ,4))
            ENDIF
            ADNC(NR,2) = PA(2)**ALP(2)*PZ(2)**ALP(3)*AD0
            ADNC(NR,3) = PA(3)**ALP(2)*PZ(3)**ALP(3)*AD0
            ADNC(NR,4) = PA(4)**ALP(2)*PZ(4)**ALP(3)*AD0
            ADNC(NR,1) =(PZ(2)*ANDX*AD(NR,2)
     &                  +PZ(3)*ANT *AD(NR,3)
     &                  +PZ(4)*ANA *AD(NR,4))/(ANDX+ANT+ANA)
C
            RX   = ALP(1)*RHOG(NR)
            PROF0 = 1.D0-RX**PROFN1
            IF(PROF0.LE.0.D0) THEN
               PROF1=0.D0
               PROF2=0.D0
            ELSE
               PROF1=PROF0**PROFN2
               PROF2=PROFN2*PROF0**(PROFN2-1.D0)
            ENDIF
            PROF   = PROF1+PNSS(1)/(PN(1)-PNSS(1))
            DPROF  =-PROFN1*RX**(PROFN1-1.D0)*PROF2
C
            DO NS=1,NSM
               AD  (NR,NS) = CNP*ADNC(NR,NS)
               AVNC(NR,NS) = AD(NR,NS)*DPROF/PROF
               AVDW(NR,NS) = 0.D0
            ENDDO
         ENDDO
      ELSEIF(MDLAD.EQ.2) THEN
         DO NR=1,NRMAX
            IF(NR.EQ.NRMAX) THEN
               ANE =PNSS(1)
               ANDX=PNSS(2)
               ANT =PNSS(3)
               ANA =PNSS(4)
            ELSE
               ANE = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
               ANDX= 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
               ANT = 0.5D0*(RN(NR+1,3)+RN(NR  ,3))
               ANA = 0.5D0*(RN(NR+1,4)+RN(NR  ,4))
            ENDIF
            ADNC(NR,2) = AD0*AKDW(NR,2)
            ADNC(NR,3) = AD0*AKDW(NR,3)
            ADNC(NR,4) = AD0*AKDW(NR,4)
            ADNC(NR,1) =(PZ(2)*ANDX*AD(NR,2)
     &                  +PZ(3)*ANT *AD(NR,3)
     &                  +PZ(4)*ANA *AD(NR,4))/(ANDX+ANT+ANA)
C
            RX   = ALP(1)*RHOG(NR)
            PROF0 = 1.D0-RX**PROFN1
            IF(PROF0.LE.0.D0) THEN
               PROF1=0.D0
               PROF2=0.D0
            ELSE
               PROF1=PROF0**PROFN2
               PROF2=PROFN2*PROF0**(PROFN2-1.D0)
            ENDIF
            PROF   = PROF1*(PN(1)-PNSS(1))+PNSS(1)
            DPROF  = -PROFN1*RX**(PROFN1-1.D0)*PROF2
     &                *(PN(1)-PNSS(1))*ALP(1)/RA *1.5D0
C
            DO NS=1,NSM
               AD  (NR,NS) = CNP*ADNC(NR,NS)
               AVNC(NR,NS) = AD(NR,NS)*DPROF/PROF
               AVDW(NR,NS) = 0.D0
            ENDDO
         ENDDO
      ELSEIF(MDLAD.EQ.3) THEN
C     *** Hinton & Hazeltine model w/o anomalous part of ***
C     *** transport effect of heat pinch ***
         DO NR=1,NRMAX
            IF(NR.EQ.NRMAX) THEN
               ANE = PNSS(1)
               ANI = PNSS(2)
               TE  = PTS(1)
            ELSE
               ANE = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
               ANI = 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
               TE  = 0.5D0*(RT(NR+1,1)+RT(NR  ,1))
            ENDIF
C
            IF(MDDW.EQ.0) THEN
               DO NS=1,NSM
                  ADDW(NR,NS) = AD0*AKDW(NR,NS)
               ENDDO
            ENDIF
C
            ZEFFL = ZEFF(NR)
            BPL   = BP(NR)
            QPL   = QP(NR)
            IF(QPL.GT.100.D0) QPL=100.D0
            EZOHL = EZOH(NR)
            EPS   = EPSRHO(NR)
            EPSS  = SQRT(EPS)**3
            VTE   = SQRT(TE*RKEV/AME)
            TAUE  = FTAUE(ANE,ANI,TE,ZEFFL)
            RNUE  = ABS(QPL)*RR/(TAUE*VTE*EPSS)
            RHOE2 = 2.D0*AME*TE*RKEV/(PZ(1)*AEE*BPL)**2
C
            RK13E = RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)
     &             /(1.D0+RC13*RNUE*EPSS)
C
            H     = BB/(BB+BPL)
C
            RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)
     &           +(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
C
            IF(MDEDGE.EQ.1.AND.NR.GE.NREDGE) CDP=CSPRS
            DO NS=1,NSM
               AVNC(NR,NS) = -(RK13E*SQRT(EPS)*EZOHL)/BPL/H
               ADNC(NR,NS) = SQRT(EPS)*RHOE2/TAUE*RK11E
               AD  (NR,NS) = CDP*ADDW(NR,NS)+CNP*ADNC(NR,NS)
               AVDW(NR,NS) = 0.D0
            ENDDO
C
         ENDDO
      ELSEIF(MDLAD.EQ.4) THEN
C     *** Hinton & Hazeltine model with anomalous part of ***
C     *** transport effect of heat pinch ***
         DO NR=1,NRMAX
            IF(NR.EQ.NRMAX) THEN
               ANE = PNSS(1)
               ANI = PNSS(2)
               TE  = PTS(1)
            ELSE
               ANE = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
               ANI = 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
               TE  = 0.5D0*(RT(NR+1,1)+RT(NR  ,1))
            ENDIF
C
            IF(MDDW.EQ.0) THEN
               DO NS=1,NSM
                  ADDW(NR,NS) = AD0*AKDW(NR,NS)
               ENDDO
            ENDIF
C
            ZEFFL = ZEFF(NR)
            BPL   = BP(NR)
            QPL   = QP(NR)
            IF(QPL.GT.100.D0) QPL=100.D0
            EZOHL = EZOH(NR)
            EPS   = EPSRHO(NR)
            EPSS  = SQRT(EPS)**3
            VTE   = SQRT(TE*RKEV/AME)
            TAUE  = FTAUE(ANE,ANI,TE,ZEFFL)
            RNUE  = ABS(QPL)*RR/(TAUE*VTE*EPSS)
            RHOE2 = 2.D0*AME*TE*RKEV/(PZ(1)*AEE*BPL)**2
C
            RK13E = RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)
     &             /(1.D0+RC13*RNUE*EPSS)
C
            H     = BB/(BB+BPL)
C
            RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)
     &           +(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
C
            IF(MDEDGE.EQ.1.AND.NR.GE.NREDGE) CDP=CSPRS
            DO NS=1,NSM
               AVNC(NR,NS) =-(RK13E*SQRT(EPS)*EZOHL)/BPL/H
               ADNC(NR,NS) = SQRT(EPS)*RHOE2/TAUE*RK11E
               AD  (NR,NS) = CDP*ADDW(NR,NS)+CNP*ADNC(NR,NS)
               AVDW(NR,NS) =-AV0*ADDW(NR,NS)*RHOG(NR)/RA
            ENDDO
         ENDDO
      ELSE
         IF(MDLAD.NE.0) WRITE(6,*) 'XX INVALID MDLAD : ',MDLAD
         DO NS=1,NSM
         DO NR=1,NRMAX
            AD  (NR,NS)=0.D0
            AV  (NR,NS)=0.D0
            ADNC(NR,NS)=0.D0
            AVNC(NR,NS)=0.D0
            AVDW(NR,NS)=0.D0
         ENDDO
         ENDDO
      ENDIF
      ENDIF
      IF(MDEDGE.EQ.1) CDP=CDPSV
C
C     ***** NET PARTICLE PINCH *****
C
      DO NS=1,NSLMAX
         DO NR=1,NRMAX
            IF(MDEDGE.EQ.1.AND.NR.GE.NREDGE) CDP=CSPRS
            AV(NR,NS)=CDP*AVDW(NR,NS)+CNP*AVNC(NR,NS)
         ENDDO
      ENDDO
      IF(MDEDGE.EQ.1) CDP=CDPSV
C
C     ***** OFF-DIAGONAL TRANSPORT COEFFICIENTS *****
C
C     ADLP : particle flux coefficient for pressure gradient
C     ADLD : particle flux coefficient for density gradient
C
      IF(MDDIAG.EQ.1) THEN
         DO NR=1,NRMAX
            IF(MDEDGE.EQ.1.AND.NR.GE.NREDGE) CDP=CSPRS
            DO NS=1,NSLMAX
               IF(MDDW.EQ.0) ADDW(NR,NS) = AD0*AKDW(NR,NS)
               DO NS1=1,NSLMAX
                  IF(NS.EQ.NS1) THEN
                     ADLD(NR,NS,NS1)= CDP*  ADDW(NR,NS)
     &                               +CNP*(-ADNCT(NR,NS,NS1))
                     ADLP(NR,NS,NS1)= CNP*( ADNCT(NR,NS,NS1)
     &                                     +ADNCP(NR,NS,NS1))
                  ELSE
                     ADLD(NR,NS,NS1)= CNP*(-ADNCT(NR,NS,NS1))
                     ADLP(NR,NS,NS1)= CNP*( ADNCT(NR,NS,NS1)
     &                                     +ADNCP(NR,NS,NS1))
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ELSEIF(MDDIAG.EQ.2) THEN
         DO NR=1,NRMAX
            IF(MDEDGE.EQ.1.AND.NR.GE.NREDGE) CDP=CSPRS
            DO NS=1,NSLMAX
               DO NS1=1,NSLMAX
                  IF(NS.EQ.NS1) THEN
                     ADLD(NR,NS,NS1)= CDP*ADDWD(NR,NS,NS1)
     &                               +CNP*ADNC(NR,NS)
                     ADLP(NR,NS,NS1)= CDP*ADDWP(NR,NS,NS1)
                  ELSE
                     ADLD(NR,NS,NS1)= CDP*ADDWD(NR,NS,NS1)
                     ADLP(NR,NS,NS1)= CDP*ADDWP(NR,NS,NS1)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ELSEIF(MDDIAG.EQ.3) THEN
         DO NR=1,NRMAX
            IF(MDEDGE.EQ.1.AND.NR.GE.NREDGE) CDP=CSPRS
            DO NS=1,NSLMAX
               DO NS1=1,NSLMAX
                  ADLD(NR,NS,NS1)= CDP*  ADDWD(NR,NS,NS1)
     &                            +CNP*(-ADNCT(NR,NS,NS1))
                  ADLP(NR,NS,NS1)= CDP*  ADDWP(NR,NS,NS1)
     &                            +CNP*( ADNCT(NR,NS,NS1)
     &                                  +ADNCP(NR,NS,NS1))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      IF(MDEDGE.EQ.1) CDP=CDPSV
C
C     /* for nuetral deuterium */
C
      ACOEF(1)=-3.231141D+1
      ACOEF(2)=-1.386002D-1
      BCOEF(1)=-3.330843D+1
      BCOEF(2)=-5.738374D-1
      BCOEF(3)=-1.028610D-1
      BCOEF(4)=-3.920980D-3
      BCOEF(5)= 5.964135D-4
      DO NR=1,NRMAX
         SUMA=0.D0
         SUMB=0.D0
         EDCM=0.5D0*(1.5D0*RT(NR,2)*1.D3)
         DO NA=1,2
            SUMA=SUMA+ACOEF(NA)*(LOG(4.D0*EDCM))**(NA-1)
         ENDDO
         DO NB=1,5
            SUMB=SUMB+BCOEF(NB)*(LOG(4.D0*EDCM))**(NB-1)
         ENDDO
         SGMNI=2.D0*EXP(SUMA)*1.D-4
         SGMNN=2.D0*EXP(SUMB)*1.D-4
C
         VNI=SQRT(RT(NR,2)*RKEV/(PA(2)*AMM))
         VNC=SQRT(0.025D-3*RKEV/(PA(7)*AMM))
         VNH=SQRT(RT(NR,2)*RKEV/(PA(8)*AMM))
         CFNCI =RN(NR,2)*1.D20*SGMNI*VNI
         CFNCNC=RN(NR,7)*1.D20*SGMNN*VNC
         CFNCNH=RN(NR,8)*1.D20*SGMNN*VNH
         CFNHI =RN(NR,2)*1.D20*SGMNI*VNI
         CFNHNC=RN(NR,7)*1.D20*SGMNN*VNH
         CFNHNH=RN(NR,8)*1.D20*SGMNN*VNH
         AD(NR,7) = CNN*VNC**2/(CFNCI+CFNCNC+CFNCNH)
         AD(NR,8) = CNN*VNH**2/(CFNHI+CFNHNH+CFNHNC)
C
         AV(NR,7) = 0.D0
         AV(NR,8) = 0.D0
      ENDDO
C
C        ****** AVK : HEAT PINCH ******
C
C     --- NEOCLASSICAL PART ---
C
      IF(MDNCLS.EQ.0) THEN
C     NCLASS has already calculated neoclassical heat pinch(AVKNC)
C     beforehand if MDNCLS=1 so that MDLAVK becomes no longer valid.
      IF(MDLAVK.EQ.1) THEN
         DO NR=1,NRMAX
            DO NS=1,NSM
               AVKNC(NR,NS) =-RHOG(NR)*CHP
               AVKDW(NR,NS) = 0.D0
            ENDDO
         ENDDO
      ELSEIF(MDLAVK.EQ.2) THEN
         DO NR=1,NRMAX
            DO NS=1,NSM
               AVKNC(NR,NS) =-RHOG(NR)*(CHP*1.D6)
     &                       /(ANE*1.D20*TE*RKEV)
               AVKDW(NR,NS) = 0.D0
            ENDDO
         ENDDO
      ELSEIF(MDLAVK.EQ.3) THEN
C     *** Hinton & Hazeltine model ***
         DO NR=1,NRMAX
            IF(NR.EQ.NRMAX) THEN
               ANE = PNSS(1)
               ANI = PNSS(2)
               TE  = PTS(1)
               TD  = PTS(2)
            ELSE
               ANE = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
               ANI = 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
               TE  = 0.5D0*(RT(NR+1,1)+RT(NR  ,1))
               TD  = 0.5D0*(RT(NR+1,2)+RT(NR  ,2))
            ENDIF
            ANED= ANE/ANI
C
            ZEFFL = ZEFF(NR)
            BPL   = BP(NR)
            QPL   = QP(NR)
            IF(QPL.GT.100.D0) QPL=100.D0
            EZOHL = EZOH(NR)
            EPS   = EPSRHO(NR)
            EPSS  = SQRT(EPS)**3
            VTE   = SQRT(TE*RKEV/AME)
            VTD   = SQRT(TD*RKEV/(PA(2)*AMM))
            TAUE  = FTAUE(ANE,ANI,TE,ZEFFL)
            TAUD  = FTAUI(ANE,ANI,TD,PZ(2),PA(2))
            RNUE  = ABS(QPL)*RR/(TAUE*VTE*EPSS)
            RNUD  = ABS(QPL)*RR/(TAUD*VTD*EPSS)
C
            RK13E = RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)
     &             /(1.D0+RC13*RNUE*EPSS)
            RK23E = RK23/(1.D0+RA23*SQRT(RNUE)+RB23*RNUE)
     &             /(1.D0+RC23*RNUE*EPSS)
            RK3D  =((1.17D0-0.35D0*SQRT(RNUD))/(1.D0+0.7D0*SQRT(RNUD))
     &             -2.1D0*(RNUD*EPSS)**2)/(1.D0+(RNUD*EPSS)**2)
C
            H     = BB/(BB+BPL)
            AVKNC(NR,1) = (-RK23E+2.5D0*RK13E)*SQRT(EPS)*EZOHL/BPL/H
            AVK(NR,1)   = CNH*AVKNC(NR,1)
            DO NS=2,NSM
               AVKNC(NR,NS) = RK3D/((1.D0+(RNUE*EPSS)**2)*PZ(2))
     &                       *RK13E*SQRT(EPS)*EZOHL/BPL/H*ANED
               AVKDW(NR,NS) = 0.D0
            ENDDO
         ENDDO
      ELSE
         IF(MDLAVK.NE.0) WRITE(6,*) 'XX INVALID MDLAVK : ',MDLAVK
         DO NS=1,NSM
         DO NR=1,NRMAX
            AVKNC(NR,NS)=0.D0
            AVKDW(NR,NS)=0.D0
            AVK  (NR,NS)=0.D0
         ENDDO
         ENDDO
      ENDIF
      ENDIF
C
C     ***** NET HEAT PINCH *****
C
      DO NR=1,NRMAX
         DO NS=1,NSLMAX
            AVKDW(NR,NS)=CDH*AVKDW(NR,NS)
            AVK(NR,NS)=CDH*AVKDW(NR,NS)+CNH*AVKNC(NR,NS)
         ENDDO
      ENDDO
C
      RETURN
      END
C
      REAL*8 FUNCTION TRCOFS(S,ALPHA,RKCV)
C
      IMPLICIT NONE
      REAL*8 S,ALPHA,RKCV
      REAL*8 SA,FS1,FS2
C
      IF(ALPHA.GE.0.D0) THEN
         SA=S-ALPHA
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
         SA=ALPHA-S
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
      REAL*8 FUNCTION TRCOFSX(S,ALPHA,RKCV,EPSA)
C
      IMPLICIT NONE
      REAL*8 S,ALPHA,RKCV,EPSA
      REAL*8 SA,FS1,FS2
C
      IF(ALPHA.GE.0.D0) THEN
C        SA=S-ALPHA
C        SA=S-(1.D0-(2.D0*ALPHA)/(1+3.D0*(ALPHA)**2))*ALPHA
         SA=S-(1.D0-(2.D0*ALPHA)/(1+6.D0*ALPHA))*ALPHA
         IF(SA.GE.0.D0) THEN
           FS1=((1.D0+RKCV)**2.5D0)*(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0)
     &         /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA*(1.D0+RKCV)
     &                      +2.0D0*SA*SA*SA*(1.D0+RKCV)**2.5D0))
         ELSE
            FS1=((1.D0+RKCV)**2.5D0)/SQRT(2.D0*(1.D0-2.D0*SA)
     &                       *(1.D0-2.D0*SA+3.D0*SA*SA*(1.D0+RKCV)))
         ENDIF
         IF(RKCV.GT.0.D0) THEN
            FS2=SQRT(RKCV/EPSA)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ELSE
C        SA=ALPHA-S
C        SA=(1.D0-(2.D0*ALPHA)/(1+3.D0*(ALPHA)**2))*ALPHA-S
         SA=(1.D0-(2.D0*ALPHA)/(1.D0+6.D0*ALPHA))*ALPHA-S
         IF(SA.GE.0.D0) THEN
            FS1=((1.D0+RKCV)**2.5D0)*(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0)
     &         /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA*(1.D0+RKCV)
     &                      +2.0D0*SA*SA*SA*(1.D0+RKCV)**2.5D0))
         ELSE
            FS1=((1.D0+RKCV)**2.5D0)/SQRT(2.D0*(1.D0-2.D0*SA)
     &                       *(1.D0-2.D0*SA+3.D0*SA*SA*(1.D0+RKCV)))
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
      REAL*8 FUNCTION TRCOFSS(S,ALPHA)
C
      IMPLICIT NONE
      REAL*8 S,ALPHA
      REAL*8 SA,FS1
C
      SA=S-ALPHA
      FS1=2.D0*SA**2/(1.D0+(2.D0/9.D0)*SQRT(ABS(SA))**5)
      TRCOFSS=FS1
      RETURN
      END
C
      REAL*8 FUNCTION TRCOFT(S,ALPHA,RKCV,EPSA)
C
      IMPLICIT NONE
      REAL*8 S,ALPHA,RKCV,EPSA
      REAL*8 SA,FS1,FS2
C
      IF(ALPHA.GE.0.D0) THEN
         SA=S-ALPHA
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
         SA=ALPHA-S
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
      REAL*8 FUNCTION RLAMBDA(X)
C
      IMPLICIT NONE
      REAL*8 X
      REAL*8 AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Y
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
C
C     *** for Mixed Bohm/gyro-Bohm model ***
C
      REAL*8 FUNCTION FBHM(WEXB,AGITG,S)
C
      IMPLICIT NONE
      REAL*8 WEXB,AGITG,S
C
      FBHM=1.D0/(1.D0+EXP(20.D0*(0.05D0+WEXB/AGITG-S)))
C
      RETURN
      END
C
C     *** ExB shearing effect for CDBM model ***
C
      REAL*8 FUNCTION FEXB(X,S,ALPHA)
C
      IMPLICIT NONE
      REAL*8 X,S,ALPHA,BETA,GAMMA,ALPHAL,A
C
      IF(ABS(ALPHA).LT.1.D-3) THEN
         ALPHAL=1.D-3
      ELSE
         ALPHAL=ABS(ALPHA)
      ENDIF
      BETA=0.5D0*ALPHAL**(-0.602D0)
     &    *(13.018D0-22.28915D0*S+17.018D0*S**2)
     &    /(1.D0-0.277584D0*S+1.42913D0*S**2)
C
      A=-10.D0/3.D0*ALPHA+16.D0/3.D0
      IF(S.LT.0.D0) THEN
         GAMMA = 1.D0/(1.1D0*SQRT(1.D0-S-2.D0*S**2-3.D0*S**3))+0.75D0
      ELSE
         GAMMA = (1.D0-0.5D0*S)/(1.1D0-2.D0*S+A*S**2+4.D0*S**3)+0.75D0
      ENDIF
      FEXB=EXP(-BETA*X**GAMMA)
C
      RETURN
      END

      
