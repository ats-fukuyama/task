!     ***********************************************************

!           CALCULATE TRANSPORT COEFFICIENTS
!
!     ***********************************************************

  SUBROUTINE TRCOEF(ierr)
    USE trcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: ns,nr
    
    DO ns=1,nsmax
       DO nr=1,nrmax
          ADTB(NR,NS)=0.D0
          AKTB(NR,NS)=0.D0
          AVPTB(NR,NS)=0.D0
          AVKTB(NR,NS)=0.D0
          ADNC(NR,NS)=0.D0
          AKNC(NR,NS)=0.D0
          AVPNC(NR,NS)=0.D0
          AVKNC(NR,NS)=0.D0
       END DO
    END DO

    IF(model_chi_nc.EQ.1.OR. &
       model_dp_nc.EQ.1.OR. &
       model_vk_nc.EQ.1.OR. &
       model_vp_nc.EQ.1.OR. &
       model_eta.EQ.1 .OR. &
       model_bs.EQ.1) THEN
       CALL TR_NCLASS(ierr)
       IF(ierr.NE.0) THEN
          WRITE(6,'(A,I6)') 'XX trcalc: tr_nclass: ierr=',ierr
          RETURN
       END IF
    END IF

    CALL TRCFDW
    CALL TRCFNC
    CALL TRCFET
    CALL TRCFAD

    RETURN
  END SUBROUTINE TRCOEF

!     ***********************************************************

  SUBROUTINE TRCFDW

    USE trcomm
    USE cdbm_mod
    USE trmodels,ONLY: mbgb_driver,mmm95_driver,mmm71_driver
    USE libitp
    USE libgrf
    USE libplog, ONLY: plog
    IMPLICIT NONE
    INTEGER:: &
         NS, NR08, NR, I
    REAL(rkind)   :: &
         AEI, AGITG, AKTBEL, AKTBIL, AKTBL, ALPHAL, ALNI, ALTI, AMA, AMD, &
         AMI, AMT, ANA, ANDX, ANE, ANT, ANYUE, ARG, CHIB, CHIGB, CLN, CLPE, &
         CLS, CLT, CRTCL, CS, DEDW, DELTA2, DIDW, DKAI, DND, DNE, &
         DPE, DPERHO, DPP, DQ, DRL, DTD, DTE, DTERHO, DTI, DVE, EPS, ETAC, &
         ETAI, EZOHL, F, FBHM, FDREV, FEXB, FS, FTAUE, FTAUI, HETA, OMEGAD, &
         OMEGAS, OMEGASS, OMEGATT, PNI, PPK, PTI, QL, RBEEDG, RG1, &
         RGL, RGLC, RHOI, RHOS, RKPP2, RLAMBDA, RLAMDA, RNM, RNP, RNST2, &
         RNTM, RNTP, RNUZ, ROUS, RPEM, RPEP, RPM, RPP, RREFF, RRSTAR, &
         RRSTAX, SA, SL, SLAMDA, TA, TAUAP, TAUD, TAUE, TD, TE, TI, TRCOFS, &
         TRCOFSS, TRCOFSX, TRCOFT, TT, VA, VTE, VTI, WCI, WE1, WPE2, XCHI0, &
         XCHI1, XCHI2, XXA, XXH, ZEFFL, FSDFIX, BPA, PROFDL
    REAL(rkind):: &
         RS,RKAPL,SHEARL,PNEL,RHONI,DPDRL,DVEXBDRL,CEXB,CKAP,chi_cdbm, &
         PAL,PZL,ADFFI,ACHIE,ACHII,ACHIEB,ACHIIB,ACHIEGB,ACHIIGB, &
         VTIL, GAMMA0, EXBfactor, SHRfactor
    REAL(rkind):: &
         CHIIW,DIFHW,CHIEW,DIMPW, &
         CHIIRB,DIFHRB,CHIERB,DIMPRB, &
         CHIIKB,DIFHKB,CHIEKB,DIMPKB, &
         ADTBHL,ADTBZL
    REAL(rkind):: &
         CHII,DIFH,CHIE,DIFZ,VIST,VISP, &
         CHIIB,DIFHB,CHIEB,CHIEG 
    INTEGER:: MODEL,ierr
    REAL(rkind),DIMENSION(NRMAX):: S_HM

      AMD=PA(2)*AMM
      AMT=PA(3)*AMM
      AMA=PA(4)*AMM
      OMEGAD = PZ(2)*AEE*BB/AMD

      select case(model_chi_tb)
      case(0:9)
         KGR1='/TI  vs r/'
         KGR2='/DTI  vs r/'
         KGR3='/CLN,CLT  vs r/'
         KGR4='/CLS  vs r/'
      case(10:19)
         KGR1='/LOG(DEDW),LOG(DIDW)  vs r/'
         KGR2='/ETAI,ETAC  vs r/'
         KGR3='/CLN,CLT  vs r/'
         KGR4='/CLS  vs r/'
      case(20:29)
         KGR1='/DTE,CRTCL  vs r/'
         KGR2='/ETAI,ETAC  vs r/'
         KGR3='/CLN,CLT  vs r/'
         KGR4='/CLS  vs r/'
      case(30:40,130:139)
         KGR1='/E$-r$=  vs r/'
         KGR2='/V$-ExB$=  vs r/'
         KGR3='@exp(-beta*WE1$+gamma$=)  vs r@'
         KGR4='@Lambda,1/(1+OmgST$+2$=)  vs r@'
      case(41:50)
         KGR1='/NST$+2$= vs r/'
         KGR2='/OmegaST vs r/'
         KGR3='@lambda vs r@'
         KGR4='@Lambda,1/(1+OmgST$+2$=),1/(1+G*WE1$+2$=)vs r@'
      case(51:59)
         KGR1='/NST$+2$= vs r/'
         KGR2='/OmegaST vs r/'
         KGR3='@lambda vs r@'
         KGR4='@Lambda,1/(1+OmgST$+2$=),1/(1+G*WE1$+2$=)vs r@'
      case(60:69)
         KGR1='/V$-ExB$=/'
         KGR2='/E$-r$=/'
         KGR3='/GROWTH RATE/'
         KGR4='/ExB SHEARING RATE/'
      case(140:149)
         KGR1='/ExB SHEARING RATE/'
         KGR2='/ExB Factor/'
         KGR3='/magnetic shear/'
         KGR4='/alpha/'
      case(150:159)
         KGR1='/chi-tb-e/'
         KGR2='/chi-tb-i/'
         KGR3='/dif-tb-h/'
         KGR4='/dif-tb-z/'
      case(160:169)
         KGR1='/chi-tb-e/'
         KGR2='/chi-tb-i/'
         KGR3='/Dif-tb-h/'
         KGR4='/Dif-tb-z/'
      case default
         KGR1='//'
         KGR2='//'
         KGR3='//'
         KGR4='//'
      end select

      DO I=1,3
         DO NR=1,NRMAX
            VGR1(NR,I)=0.D0
            VGR2(NR,I)=0.D0
            VGR3(NR,I)=0.D0
            VGR4(NR,I)=0.D0
         ENDDO
      ENDDO

!
!     ***** PP   is the total pressure (nT) *****
!     ***** DPP  is the derivative of total pressure (dnT/dr) *****
!     ***** TI   is the ion temperature (Ti) *****
!     ***** DTI  is the derivative of ion temperature (dTi/dr) *****
!     ***** TE   is the electron temperature (Te) *****
!     ***** DTE  is the derivative of electron temperature (dTe/dr) *****
!     ***** ANE  is the electron density (ne) *****
!     ***** DNE  is the derivative of electron density (dne/dr) *****
!     ***** CLN  is the density scale length ne/(dne/dr) *****
!     ***** CLT  is the ion temerature scale length Ti/(dTi/dr) *****
!     ***** DQ   is the derivative of safety factor (dq/dr) *****
!     ***** CLS  is the shear length R*q**2/(r*dq/dr) *****

      DO NR=1,NRMAX
         IF(SUMPBM.EQ.0.D0) THEN
            PADD(NR)=0.D0
         ELSE
            PADD(NR)=PBM(NR)*1.D-20/RKEV-RNF(NR,1)*RT(NR,2)
         ENDIF
!     Calculate ExB velocity in advance
!                            for ExB shearing rate calculation
         VEXB(NR)= -ER(NR)/BB
         VEXBP(NR)= -ER(NR)/(RR*BP(NR))
      ENDDO

      DO NR=1,NRMAX
!     characteristic time of temporal change of transport coefficients
         TAUK(NR)=QP(NR)*RR/SQRT(RT(NR,2)*RKEV/(PA(2)*AMM))
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

!     In the following, we assume that
!        1. pressures of beam and fusion at rho=1 are negligible,
!        2. calibration of Pbeam (from UFILE) is not necessary at rho=1.

            DO NS=2,NSM
               RNTP=RNTP+PNSS(NS)*PTS(NS)
               RNP =RNP +PNSS(NS)
               RNTM=RNTM+RN(NR-1,NS)*RT(NR-1,NS)+RN(NR  ,NS)*RT(NR  ,NS)
               RNM =RNM +RN(NR-1,NS)+RN(NR  ,NS)
            ENDDO
            RNTM= RNTM+RW(NR-1,1)+RW(NR-1,2)+RW(NR  ,1)+RW(NR  ,2)
            RPP = RNTP+PNSS(1)   *PTS(1)
            RPM = RNTM+RN(NR-1,1)*RT(NR-1,1) &
                      +RN(NR  ,1)*RT(NR  ,1)+PADD(NR-1)+PADD(NR)
            RPEP= PNSS(1)   *PTS(1)
            RPEM= 0.5D0*(RN(NR-1,1)*RT(NR-1,1)+RN(NR  ,1)*RT(NR  ,1))
            RNTM= 0.5D0*RNTM
            RNM = 0.5D0*RNM
            RPM = 0.5D0*RPM

            DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1), &
                 RHOG(NR),RHOM(NR),RHOM(NR-1))
            DNE=DERIV3P(PTS(1),RN(NR,1),RN(NR-1,1), &
                 RHOG(NR),RHOM(NR),RHOM(NR-1))
            DPE=DERIV3P(PNSS(1)*PTS(1), RN(NR,1)*RT(NR,1), &
                 RN(NR-1,1)*RT(NR-1,1), &
                 RHOG(NR),RHOM(NR),RHOM(NR-1))
            DTD=DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2), &
                 RHOG(NR),RHOM(NR),RHOM(NR-1))
            DND=DERIV3P(PTS(2),RN(NR,2),RN(NR-1,2), &
                 RHOG(NR),RHOM(NR),RHOM(NR-1))
            ZEFFL=ZEFF(NR)
            EZOHL=EZOH(NR)

            DPP = (RPP-RPM)*DRL

            TI  = RNTP/RNP
            DTI = (RNTP/RNP-RNTM/RNM)*DRL
         ELSE
!     density and temperature for each species on grid
            ANE    = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
            ANDX   = 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
            ANT    = 0.5D0*(RN(NR+1,3)+RN(NR  ,3))
            ANA    = 0.5D0*(RN(NR+1,4)+RN(NR  ,4))
            TE     = 0.5D0*(RT(NR+1,1)+RT(NR  ,1))
            TD     = 0.5D0*(RT(NR+1,2)+RT(NR  ,2))
            TT     = 0.5D0*(RT(NR+1,3)+RT(NR  ,3))
            TA     = 0.5D0*(RT(NR+1,4)+RT(NR  ,4))

!     incremental presssure and density for ion species
            DO NS=2,NSM
               RNTP=RNTP+RN(NR+1,NS)*RT(NR+1,NS)
               RNP =RNP +RN(NR+1,NS)
               RNTM=RNTM+RN(NR  ,NS)*RT(NR  ,NS)
               RNM =RNM +RN(NR  ,NS)
            ENDDO
!     incremental presssure and density for fast particles
            RNTP= RNTP+RW(NR+1,1)+RW(NR+1,2)
            RNTM= RNTM+RW(NR  ,1)+RW(NR  ,2)
!     incremental presssure and density for electron
            RPP = RNTP+RN(NR+1,1)*RT(NR+1,1)+PADD(NR+1)
            RPM = RNTM+RN(NR  ,1)*RT(NR  ,1)+PADD(NR  )
!     electron pressure
            RPEP= RN(NR+1,1)*RT(NR+1,1)
            RPEM= RN(NR  ,1)*RT(NR  ,1)

!     gradients of temperature, density and pressure for electron
            DTE = (RT(NR+1,1)-RT(NR  ,1))*DRL
            DNE = (RN(NR+1,1)-RN(NR  ,1))*DRL
            DPE = (RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL
!     gradients of temperature and density for ion
            DTD = (RT(NR+1,2)-RT(NR  ,2))*DRL
            DND = (RN(NR+1,2)-RN(NR  ,2))*DRL
!     effective charge number and parallel electric field on grid
            ZEFFL  = 0.5D0*(ZEFF(NR+1)+ZEFF(NR))
            EZOHL  = 0.5D0*(EZOH(NR+1)+EZOH(NR))

!     pressure gradient on grid
            DPP = (RPP-RPM)*DRL

!     effective ion temperature and its gradient on grid
            TI  = 0.5D0*(RNTP/RNP+RNTM/RNM)
            DTI = (RNTP/RNP-RNTM/RNM)*DRL
         ENDIF
!         WRITE(6,'(1PE12.4)') DPP

!     safety factor and its gradient on grid
!         DQ = DERIV3(NR,RHOG,QP,NRMAX,1)
         IF(NR.EQ.1) THEN
            DQ=(4.D0*QP(2)-3.D0*QP(1)-QP(3))/(2.D0*DR)
         ELSE IF(NR.EQ.NRMAX) THEN
            DQ=(3.D0*QP(NRMAX)-4.D0*QP(NRMAX-1)+QP(NRMAX-2))/(2.D0*DR)
         ELSE
            DQ=(QP(NR+1)-QP(NR-1))/(2.D0*DR)
         ENDIF
         
         QL = QP(NR)

!     sound speed for electron
         VTE = SQRT(ABS(TE*RKEV/AME))

!     characteristic length of ion temperature
         IF(ABS(DTI).GT.1.D-32) THEN
            CLT=TI/DTI
         ELSE
            IF(DTI.GE.0.D0) THEN
               CLT = 1.D32
            ELSE
               CLT =-1.D32
            ENDIF
         ENDIF

!     characteristic length of electron density
         IF(ABS(DNE).GT.1.D-32) THEN
            CLN=ANE/DNE
         ELSE
            IF(DNE.GE.0.D0) THEN
               CLN = 1.D32
            ELSE
               CLN =-1.D32
            ENDIF
         ENDIF

!     characteristic length of electron pressure
         IF(ABS(RPEP-RPEM).GT.1.D-32) THEN
            CLPE=0.5D0*(RPEP+RPEM)/(RPEP-RPEM)/DRL
         ELSE
            IF(RPEP-RPEM.GE.0.D0) THEN
               CLPE = 1.D32
            ELSE
               CLPE =-1.D32
            ENDIF
         ENDIF

!     ???
         IF(ABS(DQ).GT.1.D-32) THEN
            CLS=QL*QL/(DQ*EPS)
         ELSE
            IF(DQ.GE.0.D0) THEN
               CLS = 1.D32
            ELSE
               CLS =-1.D32
            ENDIF
         ENDIF

!     collision time between electrons and ions
         TAUE = FTAUE(ANE,ANDX,TE,ZEFFL)
         ANYUE = 0.5D0*(1.D0+ZEFFL)/TAUE

         ROUS = DSQRT(ABS(TE)*RKEV/AMD)/OMEGAD
         PPK  = 0.3D0/ROUS
         OMEGAS  = PPK*TE*RKEV/(AEE*BB*ABS(CLN))
         OMEGATT = DSQRT(2.D0)*VTE/(RR*QL)

!     Alfven wave velocity
         PNI=ANDX+ANT+ANA
         AMI=(AMD*ANDX+AMT*ANT+AMA*ANA)/PNI
         VA=SQRT(BB**2/(RMU0*ANE*1.D20*AMI))
!     magnetic shear
         S(NR)=RHOG(NR)/QL*DQ
!     pressure gradient for MHD instability
         ALPHA(NR)=-2.D0*RMU0*QL**2*RR/BB**2*(DPP*1.D20*RKEV)
!     magnetic curvature
         RKCV(NR)=-EPS*(1.D0-1.D0/(QL*QL))

!     rotational shear
!        omega(or gamma)_e=(r/q) d(q v_exb/r)/dr
         DVE = DERIV3(NR,RHOG,VEXB,NRMAX,1)
         WEXB(NR) = (S(NR)-1.D0)*VEXB(NR)/RHOG(NR)+DVE
         DVE = DERIV3(NR,RHOG,VEXBP,NRMAX,1)
         WEXBP(NR) = RR*BP(NR)*DVE/BB
!     Doppler shear
         AGMP(NR) = QP(NR)/EPS*WEXB(NR)


!   *************************************************************
!        [  0..  9] : CONSTANT COEFFICIENT MODEL
!        [ 10.. 19] : DRIFT WAVE (+ITG +ETG) MODEL
!        [ 20.. 29] : REBU-LALLA MODEL
!        [ 30.. 39] : CURRENT-DIFFUSIVITY DRIVEN MODEL
!        [ 40.. 49] : DRIFT WAVE BALLOONING MODEL
!        [ 60.. 64] : GLF23,IFS/PPPL,Weiland models
!        [130..134] : CDBM model
!        [140..143] : Mixed Boam and gyroBohm model
!        [150..151] : mmm95 model
!        [160..162] : mmm7_1 model (ETG not included)
!
         SELECT CASE(model_chi_tb)
            CASE(0:9)

!   ***  model_chi_tb.EQ. 0   : CONSTANT*(1+A*r**2)             ***
!   ***  model_chi_tb.EQ. 1   : CONSTANT/(1-A*r**2)             ***
!   ***  model_chi_tb.EQ. 2   : CONSTANT*(dTi/dr)**B/(1-A*r**2) ***
!   ***  model_chi_tb.EQ. 3   : CONSTANT*(dTi/dr)**B*Ti**C      ***
!   ***  model_chi_tb.EQ. 4   : PROP. TO CUBIC FUNC. with Bohm  ***
!   ***  model_chi_tb.EQ. 5   : PROP. TO CUBIC FUNC.            ***

            select case(model_chi_tb)
            case(0)
               AKTBL=1.D0+CKALFA*RG(NR)**2
            case(1)
               AKTBL=1.D0/(1.D0-CKALFA*RG(NR)**2)
            case(2)
               AKTBL=1.D0/(1.D0-CKALFA*RG(NR)**2)*(ABS(DTI)*RA)**CKBETA
            case(3)
               AKTBL=1.D0*(ABS(DTI)*RA)**CKBETA*ABS(TI)**CKGUMA
            case(4)
               FSDFIX=0.1D0
               BPA=AR1RHOG(NRMAX)*RDPS/RR
               PROFDL=(PTS(1)*RKEV/(16.D0*AEE*SQRT(BB**2+BPA**2)))/FSDFIX
               AKTBL=FSDFIX*(1.D0+(PROFDL-1.D0)*(RHOG(NR)/ RA)**2)
            case(5)
               FSDFIX=1.D0
               PROFDL=20.D0
               AKTBL=FSDFIX*(1.D0+(PROFDL-1.D0)*(RHOG(NR)/ RA)**2)
            case default
               WRITE(6,*) 'XX INVALID model_chi_tb : ',model_chi_tb
               AKTBL=0.D0
            end select
            AKTB(NR,1)=CK0*AKTBL
            AKTB(NR,2)=CK1*AKTBL
            AKTB(NR,3)=CK1*AKTBL
            AKTB(NR,4)=CK1*AKTBL

            VGR1(NR,1)=TI
            VGR1(NR,2)=0.D0
            VGR2(NR,1)=DTI
            VGR2(NR,2)=0.D0
            VGR3(NR,1)=CLN
            VGR3(NR,2)=CLT
            VGR4(NR,1)=CLS
            VGR4(NR,2)=CLS

          CASE(10:19)

!      ***  model_chi_tb.EQ. 10  : etac=1                                ***
!      ***  model_chi_tb.EQ. 11  : etac=1 1/(1+exp)                      ***
!      ***  model_chi_tb.EQ. 12  : etac=1 1/(1+exp) *q                   ***
!      ***  model_chi_tb.EQ. 13  : etac=1 1/(1+exp) *(1+q^2)             ***
!      ***  model_chi_tb.EQ. 14  : etac=1+2.5*(Ln/RR-0.2) 1/(1+exp)      ***
!      ***  model_chi_tb.EQ. 15  : etac=1 1/(1+exp) func(q,eps,Ln)       ***
!      ***  model_chi_tb.EQ. 16  : (model_chi_tb=15) + ZONAL FLOW        ***

             ETAI=CLN/CLT
            select case(model_chi_tb)
            case(10)
               ETAC   = 1.D0
               RRSTAR = RR
               FDREV  = 1.D0
               HETA   = 1.D0
            case(11:13)
               ETAC   = 1.D0
               RRSTAR = RR
               FDREV  = 1.D0
!!!!!!            ARG = 6.D0*(ETAI(NR)-ETAC(NR))*RA/CLN(NR)
               ARG = 6.D0*(ETAI-ETAC)
               IF(ARG.LE.-150.D0) THEN
                  HETA = 0.D0
               ELSEIF(ARG.GT.150.D0) THEN
                  HETA = 1.D0
               ELSE
                  HETA = 1.D0/(1.D0+EXP(-ARG))
               ENDIF
            case(14)
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
            case(15:16)
               ETAC    = 1.D0
               RRSTAX  = ABS(RR*(1.D0-(2.D0+1.D0   /QL**2)*EPS) &
                             /(1.D0-(2.D0+RKAP**2/QL**2)*EPS))
               RREFF = RRSTAX/(1.2D0*ABS(CLN))
               FDREV = SQRT(2.D0*PI)*RREFF**1.5D0*(RREFF-1.5D0) *EXP(-RREFF)
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
            case default
               WRITE(6,*) 'XX INVALID model_chi_tb : ',model_chi_tb
               RRSTAR = RR
               FDREV  = 1.D0
               HETA   = 1.D0
            end select
            IF(model_chi_tb.EQ.16) THEN
               DEDW = CK0*2.5D0*OMEGAS/PPK**2 &
                    *(SQRT(EPS)*MIN(FDREV,EPS*OMEGAS/ANYUE) &
                     +OMEGAS/OMEGATT*MAX(1.D0,ANYUE/OMEGATT))
               RGL  = 2.5D0*OMEGAS*HETA &
                    *SQRT(2.D0*ABS(TI)*ABS(ETAI)*ABS(CLN)/(TE*RRSTAR))
               AMD=PA(2)*AMM
               TAUD = FTAUI(ANE,ANDX,TD,PZ(2),PA(2))
               RNUZ= 1.D0/(TAUD*SQRT(EPS))
               XCHI1=SQRT(RNUZ)/(SQRT(RNUZ)+SQRT(RGL))
               XXH=2.D0/(1.D0+PPK**2)
               RGLC=OMEGAS/(2*SQRT(2.D0)*XXH)
               XXA=XXH*PPK**2/4.D0
               IF(RGL.GT.RGLC) THEN
                  XCHI2=XXA*(SQRT(1.D0+(2.D0/XXA)*(RGL-RGLC)/RGL)-1.D0)
               ELSE
                  XCHI2=0.D0
               ENDIF
               XCHI0=SQRT(XCHI1**2+XCHI2**2)
               DIDW = CK1*XCHI0*RGL/PPK**2
            ELSE
               DEDW = CK0*2.5D0*OMEGAS/PPK**2 &
                    *(SQRT(EPS)*MIN(FDREV,EPS*OMEGAS/ANYUE) &
                     +OMEGAS/OMEGATT*MAX(1.D0,ANYUE/OMEGATT))
               DIDW = CK1*2.5D0*OMEGAS/PPK**2*HETA &
                    *SQRT(2.D0*ABS(TI)*ABS(ETAI)*ABS(CLN)/(TE*RRSTAR))
            ENDIF

            AKTB(NR,1) = CDW(1)*DEDW+CDW(2)*DIDW
            AKTB(NR,2) = CDW(3)*DEDW+CDW(4)*DIDW
            AKTB(NR,3) = CDW(5)*DEDW+CDW(6)*DIDW
            AKTB(NR,4) = CDW(7)*DEDW+CDW(8)*DIDW
            IF(model_chi_tb.EQ.12) THEN
               AKTB(NR,1) = AKTB(NR,1)*QL
               AKTB(NR,2) = AKTB(NR,2)*QL
               AKTB(NR,3) = AKTB(NR,3)*QL
               AKTB(NR,4) = AKTB(NR,4)*QL
            ELSEIF(model_chi_tb.EQ.13) THEN
               AKTB(NR,1) = AKTB(NR,1)*(1+QL**2)
               AKTB(NR,2) = AKTB(NR,2)*(1+QL**2)
               AKTB(NR,3) = AKTB(NR,3)*(1+QL**2)
               AKTB(NR,4) = AKTB(NR,4)*(1+QL**2)
            ENDIF
            VGR1(NR,1)=PLOG(DEDW,1.D-2,1.D2)
            VGR1(NR,2)=PLOG(DIDW,1.D-2,1.D2)
            VGR2(NR,1)=ETAI
            VGR2(NR,2)=ETAC
            VGR3(NR,1)=CLN
            VGR3(NR,2)=CLT
            VGR4(NR,1)=CLS
            VGR4(NR,2)=CLS

         CASE(20:29)

!      ***  model_chi_tb.EQ. 20  : Rebu-Lalla model                      ***

            select case(model_chi_tb)
            case(20)
               CRTCL=6.D-2*SQRT(EZOHL*BB**3/(ANE*1.D20*SQRT(TE*RKEV))) &
                    *SQRT(AEE**2/(RMU0*SQRT(AME)))/(QL*RKEV)
               IF(ABS(DTE).LE.CRTCL .OR. DQ.LE.0.D0 .OR. NR.EQ.1) THEN
                  DKAI=0.D0
               ELSE
                  DKAI=1.5D-1*(ABS(DTE)/TE +2.D0*ABS(DNE)/ANE) &
                       *SQRT(TE/TI)*(1.D0/EPS)*(QL**2/(DQ*BB*SQRT(RR)))*VC**2 &
                       *SQRT(RMU0*1.5D0*AMM)*(1.D0-CRTCL/ABS(DTE))
               ENDIF
               AKTB(NR,1)=                  DKAI
               AKTB(NR,2)=ZEFFL*SQRT(TE/TD)*DKAI
               AKTB(NR,3)=ZEFFL*SQRT(TE/TT)*DKAI
               AKTB(NR,4)=ZEFFL*SQRT(TE/TA)*DKAI
            case default
               WRITE(6,*) 'XX INVALID model_chi_tb : ',model_chi_tb
               AKTB(NR,1)=0.D0
               AKTB(NR,2)=0.D0
               AKTB(NR,3)=0.D0
               AKTB(NR,4)=0.D0
            end select
            VGR1(NR,1)=DTE
            VGR1(NR,2)=CRTCL
            VGR2(NR,1)=ETAI
            VGR2(NR,2)=ETAC
            VGR3(NR,1)=CLN
            VGR3(NR,2)=CLT
            VGR4(NR,1)=CLS
            VGR4(NR,2)=CLS

         CASE(30:39)

!      ***  model_chi_tb.EQ. 30  : CDBM 1/(1+s)                          ***
!      ***  model_chi_tb.EQ. 31  : CDBM F(s,alpha,kappaq)                ***
!      ***  model_chi_tb.EQ. 32  : CDBM F(s,alpha,kappaq)/(1+WE1^2)      ***
!      ***  model_chi_tb.EQ. 33  : CDBM F(s,0,kappaq)                    ***
!      ***  model_chi_tb.EQ. 34  : CDBM F(s,0,kappaq)/(1+WE1^2)          ***
!      ***  model_chi_tb.EQ. 35  : CDBM (s-alpha)^2/(1+s^2.5)            ***
!      ***  model_chi_tb.EQ. 36  : CDBM (s-alpha)^2/(1+s^2.5)/(1+WE1^2)  ***
!      ***  model_chi_tb.EQ. 37  : CDBM s^2/(1+s^2.5)                    ***
!      ***  model_chi_tb.EQ. 38  : CDBM s^2/(1+s^2.5)/(1+WE1^2)          ***
!      ***  model_chi_tb.EQ. 39  : CDBM F2(s,alpha,kappaq,a/R)           ***

            WPE2=ANE*1.D20*AEE*AEE/(AME*EPS0)
            DELTA2=VC**2/WPE2

            RNST2=0.D0
            OMEGASS=0.D0
            SLAMDA=0.D0
            RLAMDA=0.D0
            RG1=1.D0
            WE1=0.D0

            IF(MOD(model_chi_tb,2).EQ.0) THEN
               SL=(S(NR)**2+0.1D0**2)
               WE1=-QL*RR/(SL*VA)*DVE
               RG1=CWEB*FEXB(ABS(WE1),S(NR),ALPHA(NR))
!               DBDRR=DPPP*1.D20*RKEV*RA*RA/(BB**2/(2*RMU0))
!               DELTAE=SQRT(DELTA2)
!               WE1=SQRT(PA(2)/PA(1))*(QL*RR*DELTAE)/(2*SL*RA*RA)*DBDRR
            ENDIF

            select case(model_chi_tb)
            case(30)
               FS=1.D0/(1.7D0+SQRT(6.D0)*S(NR))
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            case(31)
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFS(S(NR),ALPHAL,RKCV(NR))
               IF(MDCD05.NE.0) &
                    FS=FS*(2.D0*SQRT(RKPRHO(NR))/(1.D0+RKPRHO(NR)**2))**1.5D0
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            case(32)
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFS(S(NR),ALPHAL,RKCV(NR))
!               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               IF(MDCD05.NE.0) &
                    FS=FS*(2.D0*SQRT(RKPRHO(NR))/(1.D0+RKPRHO(NR)**2))**1.5D0
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            case(33)
               FS=TRCOFS(S(NR),0.D0,RKCV(NR))
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            case(34)
               FS=TRCOFS(S(NR),0.D0,RKCV(NR))
!               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            case(35)
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFSS(S(NR),ALPHAL)
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            case(36)
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFSS(S(NR),ALPHAL)
!               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            case(37)
               FS=TRCOFSS(S(NR),0.D0)
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            case(38)
               FS=TRCOFSS(S(NR),0.D0)
!               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            case(39)
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFSX(S(NR),ALPHAL,RKCV(NR),RA/RR)
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            case default
               WRITE(6,*) 'XX INVALID model_chi_tb : ',model_chi_tb
               FS=1.D0
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**3*DELTA2*VA/(QL*RR)
            end select
            AKTB(NR,1)=AKTBEL
            AKTB(NR,2)=AKTBIL
            AKTB(NR,3)=AKTBIL
            AKTB(NR,4)=AKTBIL

            VGR1(NR,1)=FS
            VGR1(NR,2)=S(NR)
            VGR1(NR,3)=ALPHA(NR)
            VGR2(NR,1)=ER(NR)!RNST2
            VGR2(NR,2)=VEXB(NR)!OMEGASS
            VGR2(NR,3)=0.D0
!            VGR3(NR,1)=1.D0/(1.D0+RG1*WE1*WE1)!SLAMDA
            VGR3(NR,1)=RG1!SLAMDA
            VGR3(NR,2)=ABS(WE1)
            VGR3(NR,3)=0.D0
            VGR4(NR,1)=RLAMDA
            VGR4(NR,2)=1.D0/(1.D0+OMEGASS**2)
            VGR4(NR,3)=0.D0

         CASE(40:50)

!        [ 40.. 50] : DRIFT WAVE BALLOONING MODEL
             
            WPE2=ANE*1.D20*AEE*AEE/(AME*EPS0)
            DELTA2=VC**2/WPE2

            RNST2=0.D0
            OMEGASS=0.D0
            SLAMDA=0.D0
            RLAMDA=0.D0

            IF(MOD(model_chi_tb,2).EQ.0) THEN
               RG1=CWEB*FEXB(ABS(WE1),S(NR),ALPHA(NR))
               SL=SQRT(S(NR)**2+0.1D0**2)
               WE1=-QL*RR/(SL*VA)*DVE
!               DBDRR=DPPP*1.D20*RKEV*RA*RA/(BB**2/(2*RMU0))
!               DELTAE=SQRT(DELTA2)
!               WE1=SQRT(PA(2)/PA(1))*(QL*RR*DELTAE)/(2*SL*RA*RA)*DBDRR
            ENDIF

            F=VTE/VA
!
            select case(model_chi_tb)
            case(40)
               FS=1.D0/(1.7D0+SQRT(6.D0)*S(NR))
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            case(41)
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFS(S(NR),ALPHAL,RKCV(NR))
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            case(42)
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFS(S(NR),ALPHAL,RKCV(NR))
!               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            case(43)
               FS=TRCOFS(S(NR),0.D0,RKCV(NR))
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            case(44)
               FS=TRCOFS(S(NR),0.D0,RKCV(NR))
!               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            case(45)
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFSS(S(NR),ALPHAL)
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            case(46)
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFSS(S(NR),ALPHAL)
!               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            case(47)
               FS=TRCOFSS(S(NR),0.D0)
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            case(48)
               FS=TRCOFSS(S(NR),0.D0)
!               FS=FS/(1.D0+RG1*WE1*WE1)
               FS=FS*RG1
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            case(49)
               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFSX(S(NR),ALPHAL,RKCV(NR),RA/RR)
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2*VA/(QL*RR)*F
            case(50)
               AEI=(PZ(2)*ANDX+PZ(3)*ANT+PZ(4)*ANA)*AEE/PNI
               WCI=AEI*BB/AMI
               PTI=(TD*ANDX+TT*ANT+TA*ANA)/PNI
               VTI=SQRT(ABS(PTI*RKEV/AMI))
               RHOI=VTI/WCI

               ALPHAL=ALPHA(NR)*CALF
               FS=TRCOFT(S(NR),ALPHAL,RKCV(NR),RA/RR)
               SA=S(NR)-ALPHA(NR)
               RNST2=0.5D0/((1.D0-2.D0*SA+3.D0*SA*SA)*FS)
               RKPP2=RNST2/(FS*ABS(ALPHA(NR))*DELTA2)

               SLAMDA=RKPP2*RHOI**2
               RLAMDA=RLAMBDA(SLAMDA)
               OMEGAS= SQRT(RKPP2)*TE*RKEV/(AEE*BB*ABS(CLPE))
               TAUAP=(QL*RR)/VA
               OMEGASS=(OMEGAS*TAUAP)/(RNST2*SQRT(ALPHA(NR)))
               AKTBEL=CK0*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2 &
                    *VA/(QL*RR)/(RLAMDA*(1.D0+OMEGASS**2))
               AKTBIL=CK1*FS*SQRT(ABS(ALPHA(NR)))**2*DELTA2 &
                    *VA/(QL*RR)/(1.D0+OMEGASS**2)
            end select

            AKTB(NR,1)=AKTBEL
            AKTB(NR,2)=AKTBIL
            AKTB(NR,3)=AKTBIL
            AKTB(NR,4)=AKTBIL

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
!            VGR4(NR,3)=1.D0/(1.D0+RG1*WE1*WE1)
            VGR4(NR,3)=RG1

         CASE(60:65)

!      ***  model_chi_tb.EQ. 60  : GLF23 model                           ***
!      ***  model_chi_tb.EQ. 61  : GLF23 (stability enhanced version)    ***
!      ***  model_chi_tb.EQ. 62  : IFS/PPPL model                        ***
!      ***  model_chi_tb.EQ. 63  : Weiland model                         ***
!      ***  model_chi_tb.EQ. 64  : Modified Weiland model                ***
!      ***  model_chi_tb.EQ. 65  : ?                                     ***

            SELECT CASE(model_chi_tb)
            CASE(60,61)
               CALL GLF23_DRIVER(S_HM)
            CASE(62)
               CALL AITKEN(1.D0,RBEEDG,RM,RNF(:,1),2,NRMAX)
               RBEEDG=RBEEDG/PNSS(1)
               CALL IFSPPPL_DRIVER(NSTM,NRMAX,RN,RR,DR,RJCB,RHOG,RHOM, &
                    QP,S,EPSRHO,RKPRHOG,RT,BB,AMM,AME,PNSS,PTS, &
                    RNF(1:NRMAX,1), RBEEDG,NSMAX,AR1RHOG,AR2RHOG,AKTB)
            CASE(63,64)
               CALL WEILAND_DRIVER
            CASE(65)
               DPERHO=DPE/RJCB(NR)
               DTERHO=DTE/RJCB(NR)
               NR08=INT(0.8D0*NRMAX)
               CHIB  =(ABS(DPERHO*1.D3)/(ANE*BB)) &
                    *QL*QL*((RT(NR08,1)-RT(NRMAX,1))/RT(NRMAX,1))
               RHOS  = 1.02D-4*SQRT(PA(2)*TE*1.D3/PZ(2))/(RA*BB)
!               RHOS  = SQRT(2.D0*AMM/AEE)*SQRT(PA(2)*TE*1.D3)/(PZ(2)*RA*BB)
               CHIGB = RHOS*ABS(DTERHO*1.D3)/BB
               CS    = SQRT(ABS(TE*RKEV/(PA(2)*AMM)))
               ALNI  = ABS(DND/ANDX)
               ALTI  = ABS(DTD/TD)
               AGITG = 0.1D0*CS/RA*SQRT(RA*ALNI+RA*ALTI)*SQRT(TD/TE)
               WEXB(NR)=0.D0
               AKTB(NR,1)=8.D-5  *CHIB*FBHM(WEXB(NR),AGITG,S(NR))+7.D-2  *CHIGB
               AKTB(NR,2)=1.6D-4 *CHIB*FBHM(WEXB(NR),AGITG,S(NR))+1.75D-2*CHIGB
               AKTB(NR,3) = AKTB(NR,2)
               AKTB(NR,4) = AKTB(NR,2)
               VGR1(NR,1) = FBHM(WEXB(NR),AGITG,S(NR))
            END SELECT

            WPE2=ANE*1.D20*AEE*AEE/(AME*EPS0)
            DELTA2=VC**2/WPE2

            RNST2=0.D0
            OMEGASS=0.D0
            SLAMDA=0.D0
            RLAMDA=0.D0

            IF(NR.EQ.1) THEN
               DRL=RJCB(NR)/DR
               S_HM(NR) = RM(NR)*RA/(0.5D0*(QP(NR)+Q0))*(QP(NR)-Q0)*DRL
            ELSE
               DRL=RJCB(NR)/DR
               S_HM(NR)=RM(NR)/(0.5D0*(QP(NR)+QP(NR-1)))*(QP(NR)-QP(NR-1))/DR
            ENDIF

            VGR1(NR,2)=S(NR)
            VGR1(NR,3)=ALPHA(NR)
            VGR2(NR,1)=VEXB(NR)
            VGR2(NR,2)=ER(NR)
            VGR2(NR,3)=0.D0
            VGR3(NR,1)=AGMP(NR)
            VGR3(NR,2)=0.D0
            VGR3(NR,3)=0.D0
            VGR4(NR,1)=WEXB(NR)
            VGR4(NR,2)=WEXBP(NR)
            VGR4(NR,3)=0.D0

         CASE(130:139)

!      ***  model_chi_tb.EQ. 130 : CDBM model                            ***
!      ***  model_chi_tb.EQ. 131 : CDBM05 model                          ***
!      ***  model_chi_tb.EQ. 132 : CDBM model with ExB shear             ***
!      ***  model_chi_tb.EQ. 134 : CDBM05 model with ExB shear           ***

            RS=RA*RG(NR)
            RKAPL=RKPRHO(NR)
            SHEARL=S(NR)
            PNEL=ANE*1.D20
            RHONI=(AMD*ANDX+AMT*ANT+AMA*ANA)*1.D20
            DPDRl=DPP*1.D20*RKEV
            DVEXBDRL=DVE/RA
            SL=(S(NR)**2+0.1D0**2)
            WE1=-QL*RR/(SL*VA)*DVE
            RG1=CWEB*FEXB(ABS(WE1),S(NR),ALPHA(NR))
            cexb=RG1
            ckap=1.d0
            MODEL=model_chi_tb-130

            CALL cdbm(BB,RR,RS,RKAPL,QL,SHEARL,PNEL,rhoni,dpdrl, &
                 &    dvexbdrl,calf,ckap,cexb,MODEL,chi_cdbm)

            AKTBEL=(CK0/12.D0)*chi_cdbm
            AKTBIL=(CK1/12.D0)*chi_cdbm

            AKTB(NR,1)=AKTBEL
            AKTB(NR,2)=AKTBIL
            AKTB(NR,3)=AKTBIL
            AKTB(NR,4)=AKTBIL

            VGR1(NR,1)=0.d0
            VGR1(NR,2)=S(NR)
            VGR1(NR,3)=ALPHA(NR)
            VGR2(NR,1)=ER(NR)!RNST2
            VGR2(NR,2)=VEXB(NR)!OMEGASS
            VGR2(NR,3)=WEXB(NR)
            VGR3(NR,1)=0.D0
            VGR3(NR,2)=0.D0
            VGR3(NR,3)=0.D0
            VGR4(NR,1)=0.D0
            VGR4(NR,2)=0.D0
            VGR4(NR,3)=0.D0

         CASE(140:149)

!      ***  model_chi_tb.EQ. 140 : mBgB (mixed Bonm and gyro-Boahm) model ***
!      ***  model_chi_tb.EQ. 141 : mBgB model suppressed by Tara          ***
!      ***  model_chi_tb.EQ. 142 : mBgB model suppressed by Pacher (EXB)  ***
!      ***  model_chi_tb.EQ. 143 : mBgB model supp. by Pacher (EXB+SHEAR) ***

            PAL=0.D0
            PZL=0.D0
            DO NS=2,NSMAX
               PAL=PAL+PA(NS)*RN(NR,NS)
               PZL=PZL+PZ(NS)*RN(NR,NS)
            ENDDO
            PAL=PAL/RN(NR,1)
            PZL=PZL/RN(NR,1)
           
            CALL mbgb_driver(RA*RG(NR),RR,ANE*1.D20,TE,TI,DTE,DNE,PAL,PZL, &
                             QP(NR),BB,WEXBP(NR),S(NR), &
                             ADFFI,ACHIE,ACHII,ACHIEB,ACHIIB,ACHIEGB,ACHIIGB, &
                             ierr)
!            WRITE(6,'(1P6E12.4)') &
!                 RA*RG(NR),RR,ANE*1.D20,TE,TI,DTE, &
!                 DNE,PAL,PZL,QP(NR),BB,WEXBP(NR), &
!                 S(NR),ADFFI, &
!                 ACHIE,ACHII,ACHIEB,ACHIIB,ACHIEGB,ACHIIGB

            AKTBEL=ACHIE
            AKTBIL=ACHII

            AKTB(NR,1)=AKTBEL
            AKTB(NR,2)=AKTBIL
            AKTB(NR,3)=AKTBIL
            AKTB(NR,4)=AKTBIL

            VGR1(NR,1)=ER(NR)
            VTIL=SQRT(2.D0*TI*RKEV/(PAL*AMM))
            GAMMA0=VTIL/(QP(NR)*RR)
            EXBfactor=1.D0/(1.D0+(WEXBP(NR)/GAMMA0)**2)
            SHRfactor=1.D0/MAX(1.D0,(S(NR)-0.5d0)**2)
            VGR2(NR,1)=EXBfactor
            VGR2(NR,2)=SHRfactor
            VGR3(NR,1)=S(NR)
            VGR3(NR,2)=0.D0
            VGR3(NR,3)=0.D0
            VGR4(NR,1)=ALPHA(NR)
            VGR4(NR,2)=0.D0
            VGR4(NR,3)=0.D0

         CASE(150:159) ! MMM95 transport model

!      ***  model_chi_tb.EQ. 150 : mmm95 (Multi-Mode tr Model) (no ExB) ***
!      ***  model_chi_tb.EQ. 151 : mmm95 (Multi-Mode tr Model) (with ExB) ***

            CALL mmm95_driver(NR,CHIIW,DIFHW,CHIEW,DIMPW, &
                                 CHIIRB,DIFHRB,CHIERB,DIMPRB, &
                                 CHIIKB,DIFHKB,CHIEKB,DIMPKB,IERR)

            AKTBEL=CHIEW+CHIERB+CHIEKB
            AKTBIL=CHIIW+CHIIRB+CHIIKB
            ADTBHL=DIFHW+DIFHRB+DIFHKB
            ADTBZL=DIMPW+DIMPRB+DIMPKB

            AKTB(NR,1)=AKTBEL         ! electron thermal diffusivity
            AKTB(NR,2:NSMAX)=AKTBIL   ! ion thermal diffusivity
            ADTB(NR,1)=ADTBHL         ! electron partilcle diffusivity (adhoc)
            DO NS=2,NSMAX
               IF(PA(NS).LE.3.D0.AND.PZ(NS).EQ.1.D0) THEN  
                  ADTB(NR,2:NSMAX)=ADTBHL   ! hydrogen ion particle diffusivity
               ELSE
                  ADTB(NR,2:NSMAX)=ADTBZL   ! impurity ion particle diffusivity
               END IF
            END DO
            
            VGR1(NR,1)=CHIEW
            VGR1(NR,2)=CHIEW+CHIERB
            VGR1(NR,3)=CHIEW+CHIERB+CHIEKB
            VGR2(NR,1)=CHIIW
            VGR2(NR,2)=CHIIW+CHIIRB
            VGR2(NR,3)=CHIIW+CHIIRB+CHIIKB
            VGR3(NR,1)=DIFHW
            VGR3(NR,2)=DIFHW+DIFHRB
            VGR3(NR,3)=DIFHW+DIFHRB+DIFHKB
            VGR4(NR,1)=DIMPW
            VGR4(NR,2)=DIMPW+DIMPRB
            VGR4(NR,3)=DIMPW+DIMPRB+DIMPKB

         CASE(160:169) ! mmm7_1 transport model

!      ***  model_chi_tb.EQ. 160 : mmm7_1 (Multi-Mode tr Model) (no ExB) ***
!      ***  model_chi_tb.EQ. 161 : mmm7_1 (Multi-Mode tr Model) (with ExB) ***

            CALL mmm71_driver(NR,CHII,DIFH,CHIE,DIFZ,VIST,VISP, &
                                 CHIIW,DIFHW,CHIEW,CHIIB,DIFHB,CHIEB,CHIEG, &
                                 IERR)

            AKTBEL=CHIE
            AKTBIL=CHII
            ADTBHL=DIFH
            ADTBZL=DIFZ

            AKTB(NR,1)=CHIE           ! electron thermal diffusivity
            AKTB(NR,2:NSMAX)=CHII     ! ion thermal diffusivity
            ADTB(NR,1)=DIFH           ! electron partilcle diffusivity (adhoc)
            DO NS=2,NSMAX
               IF(PA(NS).LE.3.D0.AND.PZ(NS).EQ.1.D0) THEN  
                  ADTB(NR,2:NSMAX)=DIFH   ! hydrogen ion particle diffusivity
               ELSE
                  ADTB(NR,2:NSMAX)=DIFZ   ! impurity ion particle diffusivity
               END IF
            END DO
            
            VGR1(NR,1)=CHIEW
            VGR1(NR,2)=CHIEW+CHIEB
            VGR1(NR,3)=CHIEW+CHIEB+CHIEG
            VGR2(NR,1)=CHIIW
            VGR2(NR,2)=CHIIW+CHIIB
            VGR2(NR,3)=CHII
            VGR3(NR,1)=DIFHW
            VGR3(NR,2)=DIFHW+DIFHB
            VGR3(NR,3)=DIFH
            VGR4(NR,1)=VIST
            VGR4(NR,2)=VISP
            VGR4(NR,3)=0.D0

         CASE DEFAULT
            WRITE(6,*) 'XX INVALID model_chi_tb : ',model_chi_tb
         END SELECT

         SELECT CASE(model_chi_tb)
         CASE(0:149) ! ADTB not defined
            SELECT CASE(model_dp_tb)
            CASE(0)
               ADTB(NR,1:4)=0.D0
            CASE(1)
               ADTB(NR,1:4)=factor_dp_tb*AKTB(NR,NS)
            END SELECT
         END SELECT

         SELECT CASE(model_vp_tb)
         CASE(0)
            AVPTB(NR,1:4)=0.D0
         CASE(1)
            AVPTB(NR,1:4)=-factor_vp_tb*RM(NR)/RA**2*ADTB(NR,1:4)
         END SELECT
         
         SELECT CASE(model_vk_tb)
         CASE(0)
            AVKTB(NR,1:4)=0.D0
         CASE(1)
            AVKTB(NR,1:4)=-factor_vk_tb*RM(NR)/RA**2*ADTB(NR,1:4)
         END SELECT
      END DO
         
!      DO NR=1,NRMAX
!         WRITE(6,'(I5,1P6E12.4)') NR,VGR1(NR,1),VGR1(NR,2),VGR1(NR,3), &
!                                     VGR2(NR,1),VGR2(NR,2),VGR2(NR,3)
!      END DO
!      DO NR=1,NRMAX
!         WRITE(6,'(I5,1P6E12.4)') NR,VGR3(NR,1),VGR3(NR,2),VGR3(NR,3), &
!                                     VGR4(NR,1),VGR4(NR,2),VGR4(NR,3)
!      END DO
         
!            CALL PAGES
!            CALL GRD1D(1,RG,VGR1,NRMAX,NRMAX,3,'CHIE vs RS')
!            CALL GRD1D(2,RG,VGR2,NRMAX,NRMAX,3,'CHII vs RS')
!            CALL GRD1D(3,RG,VGR3,NRMAX,NRMAX,3,'ADIH vs RS')
!            CALL GRD1D(4,RG,VGR4,NRMAX,NRMAX,3,'ADIZ vs RS')
!            CALL PAGEE

      RETURN
  END SUBROUTINE TRCFDW
      
!     ***********************************************************

  SUBROUTINE TRCFNC

    USE trcomm
      IMPLICIT NONE
      INTEGER:: NR, NS, NS1
      REAL(rkind)   :: AMA, AMD, AMT, ANA, ANDX, ANE, ANT, CHECK, &
           DELDA, EPS, EPSS, F1, F2, FTAUE, FTAUI, QL, RALPHA, &
           RHOA2, RHOD2, RHOE2, RHOT2, RK22E, RK2A, RK2D, RK2T, &
           RMUSA, RMUSD, RMUST, RNUA, RNUD, RNUE, RNUT, TA, TAUA, &
           TAUD, TAUE, TAUT, TD, TE, TERM1A, TERM1D, TERM1T, TERM2A, &
           TERM2D, TERM2T, TT, VTA, VTD, VTE, VTT, ZEFFL
      REAL(rkind) :: RK22=2.55D0, RA22=0.45D0, RB22=0.43D0, RC22=0.43D0
      REAL(rkind) :: RK2=0.66D0 , RA2=1.03D0 , RB2=0.31D0 , RC2=0.74D0
      REAL(rkind),SAVE :: CDHSV

!      DATA RK11,RA11,RB11,RC11/1.04D0,2.01D0,1.53D0,0.89D0/
!      DATA RK12,RA12,RB12,RC12/1.20D0,0.76D0,0.67D0,0.56D0/
!!      DATA RK22,RA22,RB22,RC22/2.55D0,0.45D0,0.43D0,0.43D0/
!      DATA RK13,RA13,RB13,RC13/2.30D0,1.02D0,1.07D0,1.07D0/
!      DATA RK23,RA23,RB23,RC23/4.19D0,0.57D0,0.61D0,0.61D0/
!      DATA RK33,RA33,RB33,RC33/1.83D0,0.68D0,0.32D0,0.66D0/
!!      DATA RK2 ,RA2 ,RB2 ,RC2 /0.66D0,1.03D0,0.31D0,0.74D0/

      IF(model_chi_nc.LE.1.AND. &
         model_dp_nc.LE.1.AND. &
         model_vk_nc.LE.1.AND. &
         model_vp_nc.LE.1) RETURN  ! all NCLASS or no NC

      AMD=PA(2)*AMM
      AMT=AMD
      AMA=AMD
      IF(NSMAX.GE.3) AMT=PA(3)*AMM
      IF(NSMAX.GE.4) AMA=PA(4)*AMM

      DO NR=1,NRMAX
         IF(NR.EQ.NRMAX) THEN
            ANE =PNSS(1)
            ANDX=PNSS(2)
            ANT=0.D0
            ANA=0.D0
            IF(NSMAX.GE.3) ANT =PNSS(3)
            IF(NSMAX.GE.4) ANA =PNSS(4)
            TE=PTS(1)
            TD=PTS(2)
            TT=TD
            TA=TD
            IF(NSMAX.GE.3) TT=PTS(3)
            IF(NSMAX.GE.4) TA=PTS(4)
            ZEFFL=ZEFF(NR)
         ELSE
            ANE    = 0.5D0*(RN(NR+1,1)+RN(NR  ,1))
            ANDX   = 0.5D0*(RN(NR+1,2)+RN(NR  ,2))
            ANT=0.D0
            ANA=0.D0
            IF(NSMAX.GE.3) ANT    = 0.5D0*(RN(NR+1,3)+RN(NR  ,3))
            IF(NSMAX.GE.4) ANA    = 0.5D0*(RN(NR+1,4)+RN(NR  ,4))
            TE     = 0.5D0*(RT(NR+1,1)+RT(NR  ,1))
            TD     = 0.5D0*(RT(NR+1,2)+RT(NR  ,2))
            TT=TD
            TA=TD
            IF(NSMAX.GE.3) TT     = 0.5D0*(RT(NR+1,3)+RT(NR  ,3))
            IF(NSMAX.GE.4) TA     = 0.5D0*(RT(NR+1,4)+RT(NR  ,4))
            ZEFFL  = 0.5D0*(ZEFF(NR+1)+ZEFF(NR))
         ENDIF

         QL = QP(NR)
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3

         VTE = SQRT(ABS(TE*RKEV/AME))
         VTD = SQRT(ABS(TD*RKEV/AMD))
         VTT = SQRT(ABS(TT*RKEV/AMT))
         VTA = SQRT(ABS(TA*RKEV/AMA))

         RHOE2=2.D0*AME*ABS(TE)*RKEV/(PZ(1)*AEE*BP(NR))**2
         RHOD2=2.D0*AMD*ABS(TD)*RKEV/(PZ(2)*AEE*BP(NR))**2
         RHOT2=RHOD2
         RHOA2=RHOD2
         IF(NSMAX.GE.3) RHOT2=2.D0*AMT*ABS(TT)*RKEV/(PZ(3)*AEE*BP(NR))**2
         IF(NSMAX.GE.4) RHOA2=2.D0*AMA*ABS(TA)*RKEV/(PZ(4)*AEE*BP(NR))**2

         TAUE = FTAUE(ANE,ANDX,TE,ZEFFL)
         TAUD = FTAUI(ANE,ANDX,TD,PZ(2),PA(2))
         TAUT=TAUD
         TAUA=TAUD
         IF(NSMAX.GE.3) TAUT = FTAUI(ANE,ANT ,TT,PZ(3),PA(3))
         IF(NSMAX.GE.4) TAUA = FTAUI(ANE,ANA ,TA,PZ(4),PA(4))

         RNUE=ABS(QL)*RR/(TAUE*VTE*EPSS)
         RNUD=ABS(QL)*RR/(TAUD*VTD*EPSS)
         RNUT=ABS(QL)*RR/(TAUT*VTT*EPSS)
         RNUA=ABS(QL)*RR/(TAUA*VTA*EPSS)

!     ***** NEOCLASSICAL TRANSPORT (HINTON, HAZELTINE) *****

         SELECT CASE(model_chi_nc)
         CASE(2)

!            RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE) &
!                       +(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
!            RK12E=RK12*(1.D0/(1.D0+RA12*SQRT(RNUE)+RB12*RNUE) &
!                       +(EPSS*RC12)**2/RB12*RNUE/(1.D0+RC12*RNUE*EPSS))
            RK22E=RK22*(1.D0/(1.D0+RA22*SQRT(RNUE)+RB22*RNUE) &
                        +(EPSS*RC22)**2/RB22*RNUE/(1.D0+RC22*RNUE*EPSS))
!         RK13E=RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)/(1.D0+RC13*RNUE*EPSS)
!         RK23E=RK23/(1.D0+RA23*SQRT(RNUE)+RB23*RNUE)/(1.D0+RC23*RNUE*EPSS)
!         RK33E=RK33/(1.D0+RA33*SQRT(RNUE)+RB33*RNUE)/(1.D0+RC33*RNUE*EPSS)

            RK2D =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUD)+RB2 *RNUD) &
                        +(EPSS*RC2 )**2/RB2 *RNUD/(1.D0+RC2 *RNUD*EPSS))
            RK2T =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUT)+RB2 *RNUT) &
                        +(EPSS*RC2 )**2/RB2 *RNUT/(1.D0+RC2 *RNUT*EPSS))
            RK2A =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUA)+RB2 *RNUA) &
                        +(EPSS*RC2 )**2/RB2 *RNUA/(1.D0+RC2 *RNUA*EPSS))
!         RK3D=((1.17-0.35*SQRT(RNUD))/(1.D0+0.7*SQRT(RNUD)) 6
!              -2.1*(RNUD*EPSS)**2)/(1+(RNUD*EPSS)**2)
!         RK3T=((1.17-0.35*SQRT(RNUT))/(1.D0+0.7*SQRT(RNUT)) &
!              -2.1*(RNUT*EPSS)**2)/(1+(RNUT*EPSS)**2)
!         RK3A=((1.17-0.35*SQRT(RNUA))/(1.D0+0.7*SQRT(RNUA)) &
!              -2.1*(RNUA*EPSS)**2)/(1+(RNUA*EPSS)**2)

         AKNC(NR,1) = SQRT(EPS)*RHOE2/TAUE*RK22E
         AKNC(NR,2) = SQRT(EPS)*RHOD2/TAUD*RK2D
         AKNC(NR,3) = SQRT(EPS)*RHOT2/TAUT*RK2T
         AKNC(NR,4) = SQRT(EPS)*RHOA2/TAUA*RK2A

      CASE(3)

!     ***** CHANG HINTON *****

         DELDA=0.D0

         RALPHA=ZEFFL-1.D0

         RMUSD=RNUD*(1.D0+1.54D0*RALPHA)
         RMUST=RNUT*(1.D0+1.54D0*RALPHA)
         RMUSA=RNUA*(1.D0+1.54D0*RALPHA)

         F1=(1.D0+1.5D0*(EPS**2+EPS*DELDA)+3.D0/8.D0*EPS**3*DELDA) &
              /(1.D0+0.5D0*EPS*DELDA)
         F2=DSQRT(1.D0-EPS**2)*(1.D0+0.5D0*EPS*DELDA) &
              /(1.D0+DELDA/EPS*(DSQRT(1.D0-EPS**2)-1.D0))

         TERM1D=(0.66D0*(1.D0+1.54D0*RALPHA) &
              +(1.88D0*DSQRT(EPS)-1.54D0*EPS)*(1.D0+3.75D0*RALPHA))*F1 &
              /(1.D0+1.03D0*DSQRT(RMUSD)+0.31D0*RMUSD)
         TERM1T=(0.66D0*(1.D0+1.54D0*RALPHA) &
              +(1.88D0*DSQRT(EPS)-1.54D0*EPS)*(1.D0+3.75D0*RALPHA))*F1 &
              /(1.D0+1.03D0*DSQRT(RMUST)+0.31D0*RMUST)
         TERM1A=(0.66D0*(1.D0+1.54D0*RALPHA) &
              +(1.88D0*DSQRT(EPS)-1.54D0*EPS)*(1.D0+3.75D0*RALPHA))*F1 &
              /(1.D0+1.03D0*DSQRT(RMUSA)+0.31D0*RMUSA)

         TERM2D=0.583D0*RMUSD*EPS/(1.D0+0.74D0*RMUSD*EPSS) &
              *(1.D0+(1.33D0*RALPHA*(1.D0+0.6D0*RALPHA)) &
              /(1.D0+1.79D0*RALPHA))*(F1-F2)
         TERM2T=0.583D0*RMUST*EPS/(1.D0+0.74D0*RMUST*EPSS) &
              *(1.D0+(1.33D0*RALPHA*(1.D0+0.6D0*RALPHA)) &
              /(1.D0+1.79D0*RALPHA))*(F1-F2)
         TERM2A=0.583D0*RMUSA*EPS/(1.D0+0.74D0*RMUSA*EPSS) &
              *(1.D0+(1.33D0*RALPHA*(1.D0+0.6D0*RALPHA)) &
              /(1.D0+1.79D0*RALPHA))*(F1-F2)

         AKNC(NR,1)=0.D0
         AKNC(NR,2)=(RHOD2*SQRT(EPS))/TAUD*(TERM1D+TERM2D)
         AKNC(NR,3)=(RHOT2*SQRT(EPS))/TAUT*(TERM1T+TERM2T)
         AKNC(NR,4)=(RHOA2*SQRT(EPS))/TAUA*(TERM1A+TERM2A)

      END SELECT
      END DO

      RETURN
    END SUBROUTINE TRCFNC

!     ***********************************************************

    SUBROUTINE TRCFET

      USE trcomm
      IMPLICIT NONE
      INTEGER:: NR
      REAL(rkind):: &
           ANE, ANI, CH, COULOG, CR, EPS, EPSS, ETAS, F33, F33TEFF, &
           FT, FTAUE, FTPF, H, PHI, QL, RK33E, RLNLAME, &
           RNUE, RNUSE, RNZ, SGMSPTZ, TAUE, TAUEL, TE, TEL, VTE, XI, ZEFFL
      REAL(rkind):: RK33=1.83D0, RA33=0.68D0, RB33=0.32D0, RC33=0.66D0

      DO NR=1,NRMAX
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3

!        ****** CLASSICAL RESISTIVITY (Spitzer) from JAERI Report ******

         ANE=RN(NR,1)
         ANI=RN(NR,2)
         TE =RT(NR,1)
         ZEFFL=ZEFF(NR)
         TAUE = FTAUE(ANE,ANI,TE,ZEFFL)

         ETA(NR) = AME/(ANE*1.D20*AEE**2*TAUE)*(0.29D0+0.46D0/(1.08D0+ZEFFL))

!        ****** NEOCLASSICAL RESISTIVITY (Hinton, Hazeltine) ******

         select case(model_eta)
         case(2)
            IF(NR.EQ.1) THEN
               QL= 0.25D0*(3.D0*Q0+QP(NR))
            ELSE
               QL= 0.5D0*(QP(NR-1)+QP(NR))
            ENDIF
            VTE=SQRT(ABS(TE)*RKEV/AME)
            RNUE=ABS(QP(NR)*RR/(TAUE*VTE*EPSS))
            RK33E=RK33/(1.D0+RA33*SQRT(RNUE)+RB33*RNUE)/(1.D0+RC33*RNUE*EPSS)

            H      = BB/SQRT(BB**2+BP(NR)**2)
            FT     = 1.D0/H-SQRT(EPS)*RK33E
            ETA(NR)= ETA(NR)/FT

!        ****** NEOCLASSICAL RESISTIVITY (Hirshman, Hawryluk) ******

         case(3)
            IF(NR.EQ.1) THEN
               QL=ABS(0.25D0*(3.D0*Q0+QP(NR)))
            ELSE
               QL=ABS(0.5D0*(QP(NR-1)+QP(NR)))
            ENDIF
            ZEFFL=ZEFF(NR)
            VTE=SQRT(ABS(TE)*RKEV/AME)
            FT=FTPF(model_tpfrac,EPS)
            TAUEL=6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)*(TE*RKEV)**1.5D0 &
                 /(ANE*1.D20*AEE**4*COULOG(1,2,ANE,TE))
            RNUSE=RR*QL/(VTE*TAUEL*EPSS)
            PHI=FT/(1.D0+(0.58D0+0.20D0*ZEFFL)*RNUSE)
            ETAS=1.65D-9*COULOG(1,2,ANE,TE)/(ABS(TE)**1.5D0)
            CH=0.56D0*(3.D0-ZEFFL)/((3.D0+ZEFFL)*ZEFFL)

            ETA(NR)=ETAS*ZEFFL*(1.D0+0.27D0*(ZEFFL-1.D0)) &
                 /((1.D0-PHI)*(1.D0-CH*PHI)*(1.D0+0.47D0*(ZEFFL-1.D0)))

!        ****** NEOCLASSICAL RESISTIVITY (Hirshman, Sigmar)  ******

         case(4)
            IF(NR.EQ.1) THEN
               QL= 0.25D0*(3.D0*Q0+QP(NR))
            ELSE
               QL= 0.5D0*(QP(NR-1)+QP(NR))
            ENDIF
            ZEFFL=ZEFF(NR)
            ANE  =RN(NR,1)
            TEL  =ABS(RT(NR,1))

!     p1157 (7.36)
            ETAS = AME/(ANE*1.D20*AEE*AEE*TAUE) &
                 *( (1.D0+1.198D0*ZEFFL+0.222D0*ZEFFL**2) &
                   /(1.D0+2.966D0*ZEFFL+0.753D0*ZEFFL**2))

            FT=FTPF(MDLTPF,EPS)
            XI=0.58D0+0.2D0*ZEFFL
            CR=0.56D0/ZEFFL*(3.D0-ZEFFL)/(3.D0+ZEFFL)
!     RNUE expressions is given by the paper by Hirshman, Hawryluk.
            RNUE=SQRT(2.D0)/EPSS*RR*QL/SQRT(2.D0*TEL*RKEV/AME)/TAUE
!     p1158 (7.41)
            ETA(NR)=ETAS/(1.D0-FT/(1.D0+XI*RNUE))/(1.D0-CR*FT/(1.D0+XI*RNUE))

!        ****** NEOCLASSICAL RESISTIVITY (Sauter)  ******

         case(5)
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
            RNUE=6.921D-18*ABS(QL)*RR*ANE*1.D20*ZEFFL*rLnLame &
                 /((TE*1.D3)**2*EPSS)
            F33TEFF=FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE) &
                            +0.45D0*(1.D0-FT)*RNUE/ZEFFL**1.5D0)
            ETA(NR)=1.D0/(SGMSPTZ*F33(F33TEFF,ZEFFL))
         END SELECT
      ENDDO

      RETURN
    END SUBROUTINE TRCFET

!     ***********************************************************

    SUBROUTINE TRCFAD

        USE trcomm
      IMPLICIT NONE
      INTEGER:: NR, NS, NS1, NA, NB
      REAL(rkind)   :: ANA, ANDX, ANE, ANED, ANI, ANT, BPL, CFNCI, CFNCNC, CFNCNH, CFNHI, CFNHNC, CFNHNH, DPROF, EDCM, EPS, EPSS,&
     &             EZOHL, FTAUE, FTAUI, H, PROF, PROF0, PROF1, PROF2, QPL, RHOE2, RK11E, RK13E, RK23E, RK3D, RNUD, RNUE, RX, &
     &             SGMNI, SGMNN, SUMA, SUMB, TAUD, TAUE, TD, TE, VNC, VNH, VNI, VTD, VTE, ZEFFL
      REAL(rkind),DIMENSION(2):: ACOEF
      REAL(rkind),DIMENSION(5):: BCOEF
      REAL(rkind) :: RK11=1.04D0, RA11=2.01D0, RB11=1.53D0, RC11=0.89D0
      REAL(rkind) :: RK13=2.30D0, RA13=1.02D0, RB13=1.07D0, RC13=1.07D0
      REAL(rkind) :: RK23=4.19D0, RA23=0.57D0, RB23=0.61D0, RC23=0.61D0
      REAL(rkind),SAVE :: CDPSV

!        ****** AD : PARTICLE DIFFUSION ******
!        ****** AVP : PARTICLE PINCH ******

      IF(MDEDGE.EQ.1) CDPSV=CDP

      select case(model_dp_tb)
      case(1)
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
            ADTB(NR,2) = factor_dp_tb*AKTB(NR,2)
            ADTB(NR,3) = factor_dp_tb*AKTB(NR,3)
            ADTB(NR,4) = factor_dp_tb*AKTB(NR,4)
            ADTB(NR,1) =(PZ(2)*ANDX*AD(NR,2) &
                        +PZ(3)*ANT *AD(NR,3) &
                        +PZ(4)*ANA *AD(NR,4))/(ANDX+ANT+ANA)
         ENDDO
      case(2)
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
            ADTB(NR,2) = PA(2)**ALP(2)*PZ(2)**ALP(3)*factor_dp_tb
            ADTB(NR,3) = PA(3)**ALP(2)*PZ(3)**ALP(3)*factor_dp_tb
            ADTB(NR,4) = PA(4)**ALP(2)*PZ(4)**ALP(3)*factor_dp_tb
            ADTB(NR,1) =(PZ(2)*ANDX*ADTB(NR,2) &
                        +PZ(3)*ANT *ADTB(NR,3) &
                        +PZ(4)*ANA *ADTB(NR,4))/(ANDX+ANT+ANA)
         ENDDO
      END select
         
      SELECT CASE(model_dp_nc)
      case(2)
!     *** Hinton & Hazeltine ***
!     *** transport effect of heat pinch ***
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

            IF(MDDW.EQ.0) ADTB(NR,1:NSM) = AD0*AKTB(NR,1:NSM)

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

            RK13E = RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)/(1.D0+RC13*RNUE*EPSS)

            H     = BB/SQRT(BB**2+BPL**2)

            RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE) &
                 +(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))

            IF(MDEDGE.EQ.1.AND.NR.GE.NREDGE) CDP=CSPRS
            ADNC(NR,1:NSM) = SQRT(EPS)*RHOE2/TAUE*RK11E
            AVPNC(NR,1:NSM) = -(RK13E*SQRT(EPS)*EZOHL)/BPL/H
         ENDDO
      END SELECT

      !     ***** OFF-DIAGONAL TRANSPORT COEFFICIENTS *****

      DO NR=1,NRMAX
         DO NS=1,NSLMAX
            DO NS1=1,NSLMAX
               IF(NS.EQ.NS1) THEN
                  ADLD(NR,NS,NS1)= ADTB(NR,NS)+ADNC(NR,NS)-ADNCT(NR,NS,NS1)
                  ADLP(NR,NS,NS1)= ADNCT(NR,NS,NS1)+ADNCP(NR,NS,NS1)
               ELSE
                  ADLD(NR,NS,NS1)=-ADNCT(NR,NS,NS1)
                  ADLP(NR,NS,NS1)= ADNCT(NR,NS,NS1)+ADNCP(NR,NS,NS1)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

!     /* for nuetral deuterium */

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

         AV(NR,7) = 0.D0
         AV(NR,8) = 0.D0
      ENDDO

!        ****** AVK : HEAT PINCH ******

!     --- NEOCLASSICAL PART ---

      IF(MDNCLS.EQ.0) THEN
!     NCLASS has already calculated neoclassical heat pinch(AVKNC)
!     beforehand if MDNCLS=1 so that MDLAVK becomes no longer valid.
      select case(MDLAVK)
      case(1)
         DO NS=1,NSM
            AVKNC(1:NRMAX,NS) =-RHOG(1:NRMAX)*CHP
            AVKDW(1:NRMAX,NS) = 0.D0
         END DO
      case(2)
         DO NS=1,NSM
            AVKNC(1:NRMAX,NS) =-RHOG(1:NRMAX)*(CHP*1.D6) /(ANE*1.D20*TE*RKEV)
            AVKDW(1:NRMAX,NS) = 0.D0
         END DO
      case(3)
!     *** Hinton & Hazeltine model ***
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

            RK13E = RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)/(1.D0+RC13*RNUE*EPSS)
            RK23E = RK23/(1.D0+RA23*SQRT(RNUE)+RB23*RNUE)/(1.D0+RC23*RNUE*EPSS)
            RK3D  =((1.17D0-0.35D0*SQRT(RNUD))/(1.D0+0.7D0*SQRT(RNUD)) &
                 & -2.1D0*(RNUD*EPSS)**2)/(1.D0+(RNUD*EPSS)**2)

            H     = BB/SQRT(BB**2+BPL**2)
            AVKNC(NR,1) = (-RK23E+2.5D0*RK13E)*SQRT(EPS)*EZOHL/BPL/H
            AVK(NR,1)   = CNH*AVKNC(NR,1)
            AVKNC(NR,2:NSM) = RK3D/((1.D0+(RNUE*EPSS)**2)*PZ(2))*RK13E &
                 & *SQRT(EPS)*EZOHL/BPL/H*ANED
            AVKDW(NR,2:NSM) = 0.D0
         ENDDO
      case default
         IF(MDLAVK.NE.0) WRITE(6,*) 'XX INVALID MDLAVK : ',MDLAVK
         AVKNC(1:NRMAX,1:NSM)=0.D0
         AVKDW(1:NRMAX,1:NSM)=0.D0
         AVK  (1:NRMAX,1:NSM)=0.D0
      end select
      ENDIF

!     ***** NET HEAT PINCH *****

      AVKDW(1:NRMAX,1:NSLMAX)=CDH*AVKDW(1:NRMAX,1:NSLMAX)
      AVK(1:NRMAX,1:NSLMAX)=CDH*AVKDW(1:NRMAX,1:NSLMAX) &
           &               +CNH*AVKNC(1:NRMAX,1:NSLMAX)

      RETURN
      END SUBROUTINE TRCFAD

      FUNCTION TRCOFS(S,ALPHA,RKCV)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind):: S, ALPHA, RKCV, TRCOFS
      REAL(rkind):: FS1,FS2, SA

      IF(ALPHA.GE.0.D0) THEN
         SA=S-ALPHA
         IF(SA.GE.0.D0) THEN
            FS1=(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0) &
                 & /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA+2.0D0*SA*SA*SA))
         ELSE
            FS1=1.D0/SQRT(2.D0*(1.D0-2.D0*SA)*(1.D0-2.D0*SA+3.D0*SA*SA))
         ENDIF
         IF(RKCV.GT.0.D0) THEN
            FS2=SQRT(RKCV)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ELSE
         SA=ALPHA-S
         IF(SA.GE.0.D0) THEN
            FS1=(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0) &
                 & /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA+2.0D0*SA*SA*SA))
         ELSE
            FS1=1.D0/SQRT(2.D0*(1.D0-2.D0*SA)*(1.D0-2.D0*SA+3.D0*SA*SA))
         ENDIF
         IF(RKCV.LT.0.D0) THEN
            FS2=SQRT(-RKCV)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ENDIF
      TRCOFS=MAX(FS1,FS2)
      RETURN
      END FUNCTION TRCOFS

      FUNCTION TRCOFSX(S,ALPHA,RKCV,EPSA)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind):: S, ALPHA, RKCV, EPSA, TRCOFSX
      REAL(rkind):: FS1, FS2, SA

      IF(ALPHA.GE.0.D0) THEN
!        SA=S-ALPHA
!        SA=S-(1.D0-(2.D0*ALPHA)/(1+3.D0*(ALPHA)**2))*ALPHA
         SA=S-(1.D0-(2.D0*ALPHA)/(1+6.D0*ALPHA))*ALPHA
         IF(SA.GE.0.D0) THEN
           FS1=((1.D0+RKCV)**2.5D0)*(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0) &
     &         /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA*(1.D0+RKCV) &
     &                       +2.0D0*SA*SA*SA*(1.D0+RKCV)**2.5D0))
         ELSE
            FS1=((1.D0+RKCV)**2.5D0)/SQRT(2.D0*(1.D0-2.D0*SA) &
     &          *(1.D0-2.D0*SA+3.D0*SA*SA*(1.D0+RKCV)))
         ENDIF
         IF(RKCV.GT.0.D0) THEN
            FS2=SQRT(RKCV/EPSA)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ELSE
!        SA=ALPHA-S
!        SA=(1.D0-(2.D0*ALPHA)/(1+3.D0*(ALPHA)**2))*ALPHA-S
         SA=(1.D0-(2.D0*ALPHA)/(1.D0+6.D0*ALPHA))*ALPHA-S
         IF(SA.GE.0.D0) THEN
            FS1=((1.D0+RKCV)**2.5D0)*(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0) &
     &         /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA*(1.D0+RKCV) &
     &                       +2.0D0*SA*SA*SA*(1.D0+RKCV)**2.5D0))
         ELSE
            FS1=((1.D0+RKCV)**2.5D0) &
            & /SQRT(2.D0*(1.D0-2.D0*SA)*(1.D0-2.D0*SA+3.D0*SA*SA*(1.D0+RKCV)))
         ENDIF
         IF(RKCV.LT.0.D0) THEN
            FS2=SQRT(-RKCV/EPSA)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ENDIF
      TRCOFSX=MAX(FS1,FS2)
      RETURN
      END FUNCTION TRCOFSX

      FUNCTION TRCOFSS(S,ALPHA)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind):: S,ALPHA,TRCOFSS
      REAL(rkind):: SA,FS1

      SA=S-ALPHA
      FS1=2.D0*SA**2/(1.D0+(2.D0/9.D0)*SQRT(ABS(SA))**5)
      TRCOFSS=FS1
      RETURN
      END FUNCTION TRCOFSS

      FUNCTION TRCOFT(S,ALPHA,RKCV,EPSA)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind):: S,ALPHA,RKCV,EPSA,TRCOFT
      REAL(rkind):: FS1,FS2, SA

      IF(ALPHA.GE.0.D0) THEN
         SA=S-ALPHA
         IF(SA.GE.0.D0) THEN
            FS1=(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0) &
                 & /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA+2.0D0*SA*SA*SA))
         ELSE
            FS1=1.D0/SQRT(2.D0*(1.D0-2.D0*SA)*(1.D0-2.D0*SA+3.D0*SA*SA))
         ENDIF
         IF(RKCV.GT.0.D0) THEN
            FS2=SQRT(RKCV/EPSA)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ELSE
         SA=ALPHA-S
         IF(SA.GE.0.D0) THEN
            FS1=(1.D0+9.0D0*SQRT(2.D0)*SA**2.5D0) &
                 & /(SQRT(2.D0)*(1.D0-2.D0*SA+3.D0*SA*SA+2.0D0*SA*SA*SA))
         ELSE
            FS1=1.D0/SQRT(2.D0*(1.D0-2.D0*SA)*(1.D0-2.D0*SA+3.D0*SA*SA))
         ENDIF
         IF(RKCV.LT.0.D0) THEN
            FS2=SQRT(-RKCV/EPSA)**3/(S*S)
         ELSE
            FS2=0.D0
         ENDIF
      ENDIF
      TRCOFT=MAX(FS1,FS2)
      RETURN
      END FUNCTION TRCOFT

      FUNCTION RLAMBDA(X)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind):: X,RLAMBDA
      REAL(rkind):: AX, Y
      REAL(rkind),SAVE:: P1=1.0D0, P2=3.5156229D0, P3=3.0899424D0, &
           &         P4=1.2067492D0, P5=0.2659732D0, P6=0.360768D-1, &
           &         P7=0.45813D-2
      REAL(rkind),SAVE:: Q1=0.39894228D0,  Q2=0.1328592D-1, Q3=0.225319D-2, &
           &         Q4=-0.157565D-2, Q5=0.916281D-2, Q6=-0.2057706D-1, &
           &         Q7=0.2635537D-1, Q8=-0.1647633D-1,Q9=0.392377D-2


      AX=ABS(X)
      IF (AX.LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        RLAMBDA=(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))) *EXP(-AX)
      ELSE
        Y=3.75D0/AX
        RLAMBDA=(1.D0/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5 &
        &                          +Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END FUNCTION RLAMBDA

!     *** for Mixed Bohm/gyro-Bohm model ***

      FUNCTION FBHM(WEXB,AGITG,S)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind):: AGITG,S, WEXB,FBHM

      FBHM=1.D0/(1.D0+EXP(20.D0*(0.05D0+WEXB/AGITG-S)))

      RETURN
      END FUNCTION FBHM

!     *** ExB shearing effect for CDBM model ***

      FUNCTION FEXB(X,S,ALPHA)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind):: A, ALPHA, ALPHAL, BETA, GAMMA, S, X,FEXB

      IF(ABS(ALPHA).LT.1.D-3) THEN
         ALPHAL=1.D-3
      ELSE
         ALPHAL=ABS(ALPHA)
      ENDIF
      BETA=0.5D0*ALPHAL**(-0.602D0)*(13.018D0-22.28915D0*S+17.018D0*S**2) &
           & /(1.D0-0.277584D0*S+1.42913D0*S**2)

      A=-10.D0/3.D0*ALPHA+16.D0/3.D0
      IF(S.LT.0.D0) THEN
         GAMMA = 1.D0/(1.1D0*SQRT(1.D0-S-2.D0*S**2-3.D0*S**3))+0.75D0
      ELSE
         GAMMA = (1.D0-0.5D0*S)/(1.1D0-2.D0*S+A*S**2+4.D0*S**3)+0.75D0
      ENDIF
      FEXB=EXP(-BETA*X**GAMMA)

      RETURN
      END FUNCTION FEXB
