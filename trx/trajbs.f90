!     ***********************************************

!         BOOTSTRAP CURRENT (NCLASS)

!     ***********************************************

      SUBROUTINE TRAJBS_NCLASS

      USE TRCOMM, ONLY : AJBS, AJBSNC, BB, CJBSP, CJBST, DR, NRMAX, NSMAX, PBSCD, PNSS, PTS, RG, RM, RN, RT, rkind
      USE libitp  
      IMPLICIT NONE
      INTEGER:: NR, NS, NSW
      REAL(rkind)   :: DRPNW, DRTNW, RPNW, RTNW, SUML
      REAL(rkind),DIMENSION(NRMAX)::  AJBSL


      IF(PBSCD.LE.0.D0) RETURN

      NSW=1
      IF(NSW.EQ.0) THEN
         AJBSL(1:NRMAX)=PBSCD*AJBSNC(1:NRMAX)
      ELSE
         DO NR=1,NRMAX-1
            SUML=0.D0
            DO NS=1,NSMAX
               RTNW=0.5D0*(RT(NR+1,NS)+RT(NR,NS))
               RPNW=0.5D0*(RN(NR+1,NS)*RT(NR+1,NS)+RN(NR  ,NS)*RT(NR  ,NS))
               DRTNW=(RT(NR+1,NS)-RT(NR,NS))/DR
               DRPNW=(RN(NR+1,NS)*RT(NR+1,NS)-RN(NR,NS)*RT(NR,NS))/DR
               SUML=SUML+CJBST(NR,NS)*DRTNW/RTNW +CJBSP(NR,NS)*DRPNW/RPNW
            ENDDO
            AJBSL(NR)=-PBSCD*SUML/BB
         ENDDO

         NR=NRMAX
         SUML=0.D0
         DO NS=1,NSMAX
            RTNW=PTS(NS)
            RPNW=PNSS(NS)*PTS(NS)
            DRTNW=DERIV3P(PTS(NS),RT(NR,NS),RT(NR-1,NS),RG(NR),RM(NR),RM(NR-1))
            DRPNW=DERIV3P(PNSS(NS)*PTS(NS),RN(NR  ,NS)*RT(NR  ,NS),RN(NR-1,NS)*RT(NR-1,NS),RG(NR),RM(NR),RM(NR-1))
            SUML=SUML+CJBST(NR,NS)*DRTNW/RTNW +CJBSP(NR,NS)*DRPNW/RPNW
!CC            if(ns.eq.3) write(6,*) NR,PNSS(NS),PTS(NS)
         ENDDO
         AJBSL(NR)=-PBSCD*SUML/BB
      ENDIF
!
      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO

      RETURN
      END SUBROUTINE TRAJBS_NCLASS

!     ***********************************************

!         BOOTSTRAP CURRENT (O. SAUTER)

!     ***********************************************

      SUBROUTINE TRAJBSSAUTER

      USE TRCOMM
      USE libitp
      IMPLICIT NONE
      INTEGER:: NR, NS, NF
      REAL(rkind)   :: ANE, DPE, DPI, DRL, DTE, DTI, EPS, EPSS,&
           & F31TEFF, F32EETEFF, F32EITEFF, F34TEFF, FT, FTPF, PE,&
           & PPI, QL, RL31, RL32, RL34, RLNLAME, RLNLAMII, RNM, RNP,&
           & RNTM, RNTP, RNUE, RNUI, RPIM, RPIP, SALFA, SALFA0, TE,&
           & TI, ZEFFL
      REAL(rkind),DIMENSION(NRMAX):: AJBSL, ANI


      IF(PBSCD.LE.0.D0) RETURN

      DO NR=1,NRMAX-1

         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         QL=ABS(QP(NR))
         ZEFFL=0.5D0*(ZEFF(NR)+ZEFF(NR+1))
         DRL=1.D0/DR

         RNTP=0.D0
         RNP =0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+RN(NR+1,NS)*RT(NR+1,NS)
            RNP =RNP +RN(NR+1,NS)
            RNTM=RNTM+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR  ,NS)
         ENDDO
         DO NF=1,NFMAX
            RNTP=RNTP+RW(NR+1,NF)
            RNTM=RNTM+RW(NR  ,NF)
         END DO
         RPIP=RNTP+PADD(NR+1)
         RPIM=RNTM+PADD(NR  )

!     ****** ION PARAMETER ******

!     *** ANI  is the the ion density (ni) ***
!     *** TI   is the ion temperature (Ti) ***
!     *** DTI  is the derivative of ion temperature (dTi/dr) ***
!     *** PPI  is the ion pressure (Pi) ***
!     *** DPI  is the derivative of ion pressure (dPi/dr) ***

         ANI(NR)=0.D0
         DO NS=2,NSMAX
            ANI(NR)=ANI(NR)+0.5D0*PZ(NS)*(RN(NR+1,NS)+RN(NR,NS))
         ENDDO
         PPI=0.5D0*(RPIP+RPIM)
         DPI=(RPIP-RPIM)*DRL
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         DTI=(RNTP/RNP-RNTM/RNM)*DRL

         rLnLamii=30.D0-LOG(PZ(2)**3*SQRT(ANI(NR)*1.D20)/(ABS(TI*1.D3)**1.5D0))
         RNUI=4.90D-18*QL*RR*ANI(NR)*1.D20*PZ(2)**4*rLnLamii /(ABS(TI*1.D3)**2*EPSS)
!
!     ****** ELECTORON PARAMETER ******

!     *** ANE  is the the electron density (ne) ***
!     *** TE   is the electron temperature (Te) ***
!     *** PE   is the electron pressure (Pe) ***
!     *** DTE  is the derivative of electron temperature (dTe/dr) ***
!     *** DPE  is the derivative of electron pressure (dPe/dr) ***

         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
         TE= 0.5D0*(RT(NR+1,1)+RT(NR,1))
         PE= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
         DTE=(RT(NR+1,1)-RT(NR,1))*DRL
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL

         rLnLame=31.3D0-LOG(SQRT(ANE*1.D20)/ABS(TE*1.D3))
         RNUE=6.921D-18*QL*RR*ANE*1.D20*ZEFFL*rLnLame /(ABS(TE*1.D3)**2*EPSS)

         RPE=PE/(PE+PPI)
         FT=FTPF(MDLTPF,EPS)

!         F33TEFF=FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE)+0.45D0*(1.D0-FT)*RNUE/ZEFFL**1.5)
         F31TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)   +0.5D0*(1.D0-FT)*RNUE/ZEFFL)
         F32EETEFF=FT/(1.D0+0.26D0*(1.D0-FT)*SQRT(RNUE)+0.18D0*(1.D0-0.37D0*FT)*RNUE/SQRT(ZEFFL))
         F32EITEFF=FT/(1.D0+(1.D0+0.6D0*FT)*SQRT(RNUE) +0.85D0*(1.D0-0.37D0*FT)*RNUE*(1.D0+ZEFFL))
         F34TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)   +0.5D0*(1.D0-0.5D0*FT)*RNUE/ZEFFL)

         SALFA0=-1.17D0*(1.D0-FT)/(1.D0-0.22D0*FT-0.19D0*FT**2)
         SALFA=((SALFA0+0.25D0*(1.D0-FT**2)*SQRT(RNUI))/(1.D0+0.5D0*SQRT(RNUI))+0.315D0*RNUI**2*FT**6) &
     &        /(1.D0+0.15D0*RNUI**2*FT**6)

!         RNZ=0.58D0+0.74D0/(0.76D0+ZEFFL)
!         SGMSPTZ=1.9012D4*(TE*1.D3)**1.5/(ZEFFL*RNZ*rLnLame)
!         SGMNEO=SGMSPTZ*F33(F33TEFF,ZEFFL)

         RL31=F31(F31TEFF,ZEFFL)
         RL32=F32EE(F32EETEFF,ZEFFL)+F32EI(F32EITEFF,ZEFFL)
         RL34=F31(F34TEFF,ZEFFL)

         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV *( RL31*(DPE/PE+DPI/PE)+RL32*DTE/TE &
     &             +RL34*SALFA*(1.D0-RPE)/RPE*DTI/TI)/RDP(NR)/BB
      ENDDO

      NR=NRMAX
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         QL=ABS(QP(NR))
         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
         DRL=1.D0/DR

!     In the following, we assume that
!        1. pressures of beam and fusion at rho=1 are negligible,
!        2. calibration of Pbeam (from UFILE) is not necessary at rho=1.

         RNTP=0.D0
         RNP =0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+PNSS(NS)*PTS(NS)
            RNP =RNP +PNSS(NS)
            RNTM=RNTM+RN(NR-1,NS)*RT(NR-1,NS)+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR-1,NS)+RN(NR  ,NS)
         ENDDO
         DO NF=1,NFMAX
            RNTM=RNTM+RW(NR-1,NF)+RW(NR  ,NF)
         END DO
         RPIP=RNTP
         RPIM=RNTM+PADD(NR-1)+PADD(NR  )
         RNTM=0.5D0*RNTM
         RNM =0.5D0*RNM
         RPIM=0.5D0*RPIM

!     ****** ION PARAMETER ******

!     *** ANI  is the the ion density (ni) ***
!     *** TI   is the ion temperature (Ti) ***
!     *** DTI  is the derivative of ion temperature (dTi/dr) ***
!     *** PPI  is the ion pressure (Pi) ***
!     *** DPI  is the derivative of ion pressure (dPi/dr) ***

         ANI(NR)=0.D0
         DO NS=2,NSMAX
            ANI(NR)=ANI(NR)+PZ(NS)*PNSS(NS)
         ENDDO
         PPI=RPIP
         DPI=(RPIP-RPIM)*DRL
         TI =RNTP/RNP
         DTI=(RNTP/RNP-RNTM/RNM)*DRL

!         WRITE(6,'(A,I6,4ES12.4)') 'rLnLamii:',NR,PZ(2),ANI(NR),TI,QL
!         WRITE(6,'(A,I6,4ES12.4)') 'RNUI:    ',NR,RR,ANI(NR),rLnLamii,EPSS
         rLnLamii=30.D0-LOG(PZ(2)**3*SQRT(ANI(NR)*1.D20)/(ABS(TI*1.D3)**1.5D0))
         RNUI=4.90D-18*QL*RR*ANI(NR)*1.D20*PZ(2)**4*rLnLamii /(ABS(TI*1.D3)**2*EPSS)
!
!     ****** ELECTORON PARAMETER ******

!     *** ANE  is the the electron density (ne) ***
!     *** TE   is the electron temperature (Te) ***
!     *** PE   is the electron pressure (Pe) ***
!     *** DTE  is the derivative of electron temperature (dTe/dr) ***
!     *** DPE  is the derivative of electron pressure (dPe/dr) ***

         ANE=PNSS(1)
         TE =PTS(1)
         PE =PNSS(1)*PTS(1)
         DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPE=DERIV3P(PNSS(1)*PTS(1),RN(NR,1)*RT(NR,1),RN(NR-1,1)*RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
!
         rLnLame=31.3D0-LOG(SQRT(ANE*1.D20)/ABS(TE*1.D3))
         RNUE=6.921D-18*QL*RR*ANE*1.D20*ZEFFL*rLnLame /(ABS(TE*1.D3)**2*EPSS)
!
         RPE=PE/(PE+PPI)
         FT=FTPF(MDLTPF,EPS)

!     F33TEFF=FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE)+0.45D0*(1.D0-FT)*RNUE/ZEFFL**1.5)
         F31TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)   +0.5D0*(1.D0-FT)*RNUE/ZEFFL)
         F32EETEFF=FT/(1.D0+0.26D0*(1.D0-FT)*SQRT(RNUE)+0.18D0*(1.D0-0.37D0*FT)*RNUE/SQRT(ZEFFL))
         F32EITEFF=FT/(1.D0+(1.D0+0.6D0*FT)*SQRT(RNUE) +0.85D0*(1.D0-0.37D0*FT)*RNUE*(1.D0+ZEFFL))
         F34TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)   +0.5D0*(1.D0-0.5D0*FT)*RNUE/ZEFFL)

         SALFA0=-1.17D0*(1.D0-FT)/(1.D0-0.22D0*FT-0.19D0*FT**2)

!         WRITE(6,'(A,I6,3ES12.4)') 'SALFA:',NR,SALFA0,FT,RNUI
         SALFA=((SALFA0+0.25D0*(1.D0-FT**2)*SQRT(RNUI)) /(1.D0+0.5D0*SQRT(RNUI))+0.315D0*RNUI**2*FT**6) &
     &        /(1.D0+0.15D0*RNUI**2*FT**6)

!     RNZ=0.58D0+0.74D0/(0.76D0+ZEFFL)
!     SGMSPTZ=1.9012D4*(TE*1.D3)**1.5/(ZEFFL*RNZ*rLnLame)
!     SGMNEO=SGMSPTZ*F33(F33TEFF,ZEFFL)
!
         RL31=F31(F31TEFF,ZEFFL)
         RL32=F32EE(F32EETEFF,ZEFFL)+F32EI(F32EITEFF,ZEFFL)
         RL34=F31(F34TEFF,ZEFFL)

         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV*( RL31*(DPE/PE+DPI/PE)+RL32*DTE/TE &
     &             +RL34*SALFA*(1.D0-RPE)/RPE*DTI/TI)/RDP(NR)/BB

      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO

      RETURN
      END SUBROUTINE TRAJBSSAUTER

!     *********************
!     *  Fitting Function *
!     *********************

      FUNCTION F33(X,Z)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind):: X,Z,F33

      F33 = 1.D0-(1.D0+0.36D0/Z)*X+0.59D0/Z*X**2-0.23D0/Z*X**3

      RETURN
      END FUNCTION F33

      FUNCTION F31(X,Z)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind) X,Z,F31

      F31 = (1.D0+1.4D0/(Z+1.D0))*X &
           -1.9D0/(Z+1.D0)*X**2 &
           +0.3D0/(Z+1.D0)*X**3 &
           +0.2D0/(Z+1.D0)*X**4

      RETURN
      END FUNCTION F31

      FUNCTION F32EE(X,Z)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind) X,Z,F32EE

      F32EE = (0.05D0+0.62D0*Z)/(Z*(1.D0+0.44D0*Z))*(X-X**4) &
           +1.D0/(1.D0+0.22D0*Z)*(X**2-X**4-1.2D0*(X**3-X**4)) &
           +1.2D0/(1.D0+0.5D0*Z)*X**4

      RETURN
      END FUNCTION F32EE

      FUNCTION F32EI(X,Z)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind) X,Z,F32EI

      F32EI =-(0.56D0+1.93D0*Z)/(Z*(1.D0+0.44D0*Z))*(X-X**4) &
           +4.95D0/(1.D0+2.48D0*Z)*(X**2-X**4-0.55D0*(X**3-X**4)) &
           -1.2D0/(1.D0+0.5D0*Z)*X**4

      RETURN
      END FUNCTION F32EI

!     ************************************************

!         BOOTSTRAP CURRENT (Hirshman (Wilson))

!     ************************************************

      SUBROUTINE TRAJBSNEW

      USE TRCOMM
      USE libitp
      IMPLICIT NONE
      INTEGER:: NR, NS, NF
      REAL(rkind) :: DDD, DDX, DPE, DPI, DRL, DTE, DTI, EPS, FT, PE, PPI, RL31, RL32, RNM, RNP, RNTM, RNTP, RPIM, RPIP, TE, TI
      REAL(rkind),DIMENSION(NRMAX)::  AJBSL, ANI


      IF(PBSCD.LE.0.D0) RETURN

      DO NR=1,NRMAX-1

         EPS=EPSRHO(NR)
!         EPSS=SQRT(EPS)**3
!         QL=ABS(QP(NR))
!         ZEFFL=0.5D0*(ZEFF(NR)+ZEFF(NR+1))
         DRL=1.D0/DR

         RNTP=0.D0
         RNP= 0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+RN(NR+1,NS)*RT(NR+1,NS)
            RNP =RNP +RN(NR+1,NS)
            RNTM=RNTM+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR  ,NS)
         ENDDO
         DO NF=1,NFMAX
            RNTP=RNTP+RW(NR+1,NF)
            RNTM=RNTM+RW(NR  ,NF)
         END DO
         RPIP=RNTP+PADD(NR+1)
         RPIM=RNTM+PADD(NR  )

!     ****** ION PARAMETER ******

!     ***** ANI  is the the ion density (ni) *****
!     ***** TI   is the ion temperature (Ti) *****
!     ***** DTI  is the derivative of ion temperature (dTi/dr) *****
!     ***** PPI  is the ion pressure (Pi) *****
!     ***** DPI  is the derivative of ion pressure (dPi/dr) *****
!     ***** VTI  is the ion velocity (VTi) *****

         ANI(NR)=0.D0
         DO NS=2,NSMAX
            ANI(NR)=ANI(NR)+0.5D0*PZ(NS)*(RN(NR+1,NS)+RN(NR,NS))
         ENDDO
         PPI=0.5D0*(RPIP+RPIM)
         DPI=(RPIP-RPIM)*DRL
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         DTI=(RNTP/RNP-RNTM/RNM)*DRL
!         VTI=SQRT(ABS(TI)*RKEV/AMP)

!         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
!         rLnLam=17.3D0-DLOG(ANE)*0.5D0+DLOG(ABS(TI))*1.5D0
!         TAUI=12.D0*PI*SQRT(PI)*EPS0**2*SQRT(AMP)
!     &             *(ABS(TI)*RKEV)**1.5D0/(ANI(NR)*1.D20
!     &             *ZEFFL**4*AEE**4*rLnLam)

!         RNUI=QL*RR/(TAUI*VTI*EPSS)

!     ****** ELECTORON PARAMETER ******

!     ***** ANE  is the the electron density (ne) *****
!     ***** TE   is the electron temperature (Te) *****
!     ***** DTE  is the derivative of electron temperature (dTe/dr) ****
!     ***** PE   is the electron pressure (Pe) *****
!     ***** DPE  is the derivative of electron pressure (dPe/dr) *****
!     ***** VTE  is the electron velocity (VTe) *****

         TE= 0.5D0*(RT(NR+1,1)+RT(NR,1))
         PE= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
         DTE=(RT(NR+1,1)-RT(NR,1))*DRL
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL
!         VTE=SQRT(ABS(TE)*RKEV/AME)

!         rLnLam=15.2D0-DLOG(ANE)*0.5D0+DLOG(ABS(TE))
!         TAUE=6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)*(ABS(TE)*RKEV)**1.5D0/(ANI*1.D20*ZEFFL**2*AEE**4*rLnLam)

!         RNUE=QL*RR/(TAUE*VTE*EPSS)

!         FT=FTPF(MDLTPF,EPS)
         FT=(1.46D0*SQRT(EPS)+2.4D0*EPS)/(1.D0-EPS)**1.5D0

!         DDX=2.4D0+5.4D0*FT+2.6D0*FT**2
!         DDD=-1.17D0/(1.D0+0.46D0*FT)
!         C1=(4.D0+2.6D0*FT)/((1.D0+1.02D0*DSQRT(RNUE)+1.07D0*RNUE)*(1.D0+1.07D0*EPSS*RNUE))
!         C2=C1*TI/TE
!         C3=(7.D0+6.5D0*FT)/((1.D0+0.57D0*DSQRT(RNUE)+0.61D0*RNUE)*(1.D0+0.61D0*EPSS*RNUE))-2.5D0*C1
!         C4=((DDD+0.35D0*DSQRT(RNUI))/(1.D0+0.7D0*DSQRT(RNUI))
!     &      +2.1D0*EPS**3*RNUI**2)*C2/((1.D0-EPS**3*RNUI**2)*(1.D0+EPS**3*RNUE**2))

!         AJBSL(NR)=-PBSCD*FT*PE*1.D20*RKEV/(DDX*BP(NR))*(C1*(DPE/PE)
!     &         +C2*(DPI/PPI)+C3*(DTE/TE)+C4*(DTI/TI))

!     *** S. P. Hirshman, Phys Fluids 31, 1988 3150 ***
!     *** cited (H.R. WILSON, Nucl.Fusion 32,no.2,1992 259-263) ***

         DDX=1.414D0*PZ(2)+PZ(2)**2+FT*(0.754D0+2.657D0*PZ(2) &
     &        +2.D0*PZ(2)**2)+FT**2*(0.348D0+1.243D0*PZ(2)+PZ(2)**2)
         RL31= FT*(0.754D0+2.210D0*PZ(2)+PZ(2)**2 +FT*(0.348D0+1.243D0*PZ(2)+PZ(2)**2))/DDX
         RL32=-FT*(0.884D0+2.074D0*PZ(2))/DDX
         DDD=-1.172D0/(1.D0+0.462D0*FT)

         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV*(RL31*((DPE/PE)+(TI/(PZ(2)*TE)) &
     &        *((DPI/PPI)+DDD*(DTI/TI)))+RL32*(DTE/TE))/RDP(NR)/BB
      ENDDO

      NR=NRMAX
         EPS=EPSRHO(NR)
!         EPSS=SQRT(EPS)**3
!         QL=ABS(QP(NR))
!         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
         DRL=1.D0/DR

!     In the following, we assume that
!        1. pressures of beam and fusion at rho=1 are negligible,
!        2. calibration of Pbeam (from UFILE) is not necessary at rho=1.

         RNTP=0.D0
         RNP= 0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+PNSS(NS)*PTS(NS)
            RNP =RNP +PNSS(NS)
            RNTM=RNTM+RN(NR-1,NS)*RT(NR-1,NS)+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR-1,NS)+RN(NR  ,NS)
         ENDDO
         RNTM=RNTP+RW(NR-1,1)+RW(NR-1,2)+RW(NR  ,1)+RW(NR  ,2)
         RPIP=RNTP
         RPIM=RNTM+PADD(NR-1)+PADD(NR  )
         RNTM=0.5D0*RNTM
         RNM =0.5D0*RNM
         RPIM=0.5D0*RPIM

!     ****** ION PARAMETER ******

!     ***** ANI  is the the ion density (ni) *****
!     ***** TI   is the ion temperature (Ti) *****
!     ***** DTI  is the derivative of ion temperature (dTi/dr) *****
!     ***** PPI  is the ion pressure (Pi) *****
!     ***** DPI  is the derivative of ion pressure (dPi/dr) *****
!     ***** VTI  is the ion velocity (VTi) *****

         ANI(NR)=0.D0
         DO NS=2,NSMAX
            ANI(NR)=ANI(NR)+PZ(NS)*PNSS(NS)
         ENDDO
         PPI=0.5D0*(RPIP+RPIM)
         DPI=(RPIP-RPIM)*DRL
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         DTI=(RNTP/RNP-RNTM/RNM)*DRL
!         VTI=SQRT(ABS(TI)*RKEV/AMP)

!         ANE=PNSS(1)
!         rLnLam=17.3D0-DLOG(ANE)*0.5D0+DLOG(ABS(TI))*1.5D0
!         TAUI=12.D0*PI*SQRT(PI)*EPS0**2*SQRT(AMP)
!     &             *(ABS(TI)*RKEV)**1.5D0/(ANI(NR)*1.D20
!     &             *ZEFFL**4*AEE**4*rLnLam)

!         RNUI=QL*RR/(TAUI*VTI*EPSS)

!     ****** ELECTORON PARAMETER ******

!     ***** ANE  is the the electron density (ne) *****
!     ***** TE   is the electron temperature (Te) *****
!     ***** DTE  is the derivative of electron temperature (dTe/dr) ****
!     ***** PE   is the electron pressure (Pe) *****
!     ***** DPE  is the derivative of electron pressure (dPe/dr) *****
!     ***** VTE  is the electron velocity (VTe) *****

         TE =PTS(1)
         PE =PNSS(1)*PTS(1)
         DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPE=DERIV3P(PNSS(1)*PTS(1),RN(NR,1)*RT(NR,1),RN(NR-1,1)*RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
!         VTE=SQRT(ABS(TE)*RKEV/AME)

!         rLnLam=15.2D0-DLOG(ANE)*0.5D0+DLOG(ABS(TE))
!         TAUE=6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)
!     &             *(ABS(TE)*RKEV)**1.5D0/(ANI*1.D20*ZEFFL**2*AEE**4*rLnLam)

!         RNUE=QL*RR/(TAUE*VTE*EPSS)

!         FT=FTPF(MDLTPF,EPS)
         FT=(1.46D0*SQRT(EPS)+2.4D0*EPS)/(1.D0-EPS)**1.5D0

!         DDX=2.4D0+5.4D0*FT+2.6D0*FT**2
!         DDD=-1.17D0/(1.D0+0.46D0*FT)
!         C1=(4.D0+2.6D0*FT)/((1.D0+1.02D0*DSQRT(RNUE)+1.07D0*RNUE)*(1.D0+1.07D0*EPSS*RNUE))
!         C2=C1*TI/TE
!         C3=(7.D0+6.5D0*FT)/((1.D0+0.57D0*DSQRT(RNUE)+0.61D0*RNUE)*(1.D0+0.61D0*EPSS*RNUE))
!     &      -2.5D0*C1
!         C4=((DDD+0.35D0*DSQRT(RNUI))/(1.D0+0.7D0*DSQRT(RNUI))+2.1D0*EPS**3*RNUI**2)*C2
!     &      /((1.D0-EPS**3*RNUI**2)*(1.D0+EPS**3*RNUE**2))

!         AJBSL(NR)=-PBSCD*FT*PE*1.D20*RKEV/(DDX*BP(NR))*(C1*(DPE/PE)
!     &         +C2*(DPI/PPI)+C3*(DTE/TE)+C4*(DTI/TI))

!     *** S. P. Hirshman, Phys Fluids 31, 1988 3150 ***
!     *** cited (H.R. WILSON, Nucl.Fusion 32,no.2,1992 259-263) ***

         DDX=1.414D0*PZ(2)+PZ(2)**2+FT*(0.754D0+2.657D0*PZ(2) &
     &        +2.D0*PZ(2)**2)+FT**2*(0.348D0+1.243D0*PZ(2)+PZ(2)**2)
         RL31=FT*( 0.754D0+2.21D0*PZ(2)+PZ(2)**2+FT*(0.348D0+1.243D0 &
     &            *PZ(2)+PZ(2)**2))/DDX
         RL32=-FT*(0.884D0+2.074D0*PZ(2))/DDX
         DDD=-1.172D0/(1.D0+0.462D0*FT)

         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV &
     &        *(RL31*((DPE/PE)+(TI/(PZ(2)*TE))*((DPI/PPI)+DDD*(DTI/TI)))+RL32*(DTE/TE))/RDP(NR)/BB

      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO

      RETURN
      END SUBROUTINE TRAJBSNEW

!     ***********************************************************

!           BOOTSTRAP CURRENT (Hinton & Hazeltine)

!     ***********************************************************

      SUBROUTINE TRAJBS

      USE TRCOMM
      USE libitp
      IMPLICIT NONE
      INTEGER:: NR
      REAL(rkind)   :: A, ANDX, ANE, ANT, ANA, BPL, DPA, DPD, DPE, DPT, DRL, DTA, DTD, DTE, DTT, EPS, EPSS, FACT, &
     &             FTAUE, FTAUI, H, PAL, PDL, PEL, PTL, RK13E, RK23E, RK3A, RK3D, RK3T, RNUA, RNUD, RNUE, RNUT, TAL, TAUA,   &
     &             TAUD, TAUE, TAUT, TDL, TEL, TTL, VTA, VTD, VTE, VTT, ZEFFL
      REAL(rkind),DIMENSION(NRMAX):: AJBSL

      REAL(rkind):: RK13=2.30D0, RA13=1.02D0, RB13=1.07D0, RC13=1.07D0, RK23=4.19D0, RA23=0.57D0, RB23=0.61D0, RC23=0.61D0

!     ZEFF=1

!      DATA RK11,RA11,RB11,RC11/1.04D0,2.01D0,1.53D0,0.89D0/
!      DATA RK12,RA12,RB12,RC12/1.20D0,0.76D0,0.67D0,0.56D0/
!      DATA RK22,RA22,RB22,RC22/2.55D0,0.45D0,0.43D0,0.43D0/
!!      DATA RK13,RA13,RB13,RC13/2.30D0,1.02D0,1.07D0,1.07D0/
!!      DATA RK23,RA23,RB23,RC23/4.19D0,0.57D0,0.61D0,0.61D0/
!      DATA RK33,RA33,RB33,RC33/1.83D0,0.68D0,0.32D0,0.66D0/
!      DATA RK2 ,RA2 ,RB2 ,RC2 /0.66D0,1.03D0,0.31D0,0.74D0/

      IF(PBSCD.LE.0.D0) RETURN

      DO NR=1,NRMAX-1

         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
         ANDX=0.5D0*(RN(NR+1,2)+RN(NR,2))
         ANT =0.5D0*(RN(NR+1,3)+RN(NR,3))
         ANA =0.5D0*(RN(NR+1,4)+RN(NR,4))
         TEL=ABS(0.5D0*(RT(NR+1,1)+RT(NR,1)))
         TDL=ABS(0.5D0*(RT(NR+1,2)+RT(NR,2)))
         TTL=ABS(0.5D0*(RT(NR+1,3)+RT(NR,3)))
         TAL=ABS(0.5D0*(RT(NR+1,4)+RT(NR,4)))
         PEL=0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
         PDL=0.5D0*(RN(NR+1,2)*RT(NR+1,2)+RN(NR,2)*RT(NR,2))
         PTL=0.5D0*(RN(NR+1,3)*RT(NR+1,3)+RN(NR,3)*RT(NR,3))
         PAL=0.5D0*(RN(NR+1,4)*RT(NR+1,4)+RN(NR,4)*RT(NR,4))
         ZEFFL=0.5D0*(ZEFF(NR+1)+ZEFF(NR))

         TAUE = FTAUE(ANE,ANDX,TEL,ZEFFL)
         TAUD = FTAUI(ANE,ANDX,TDL,PZ(2),PM(2))
         TAUT = FTAUI(ANE,ANT ,TTL,PZ(3),PM(3))
         TAUA = FTAUI(ANE,ANA ,TAL,PZ(4),PM(4))

         VTE=SQRT(TEL*RKEV/AME)
         VTD=SQRT(TDL*RKEV/AMD)
         VTT=SQRT(TTL*RKEV/AMT)
         VTA=SQRT(TAL*RKEV/AMA)

         RNUE=ABS(QP(NR))*RR/(TAUE*VTE*EPSS)
         RNUD=ABS(QP(NR))*RR/(TAUD*VTD*EPSS)
         RNUT=ABS(QP(NR))*RR/(TAUT*VTT*EPSS)
         RNUA=ABS(QP(NR))*RR/(TAUA*VTA*EPSS)

!         RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)+(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
!         RK12E=RK12*(1.D0/(1.D0+RA12*SQRT(RNUE)+RB12*RNUE)+(EPSS*RC12)**2/RB12*RNUE/(1.D0+RC12*RNUE*EPSS))
!         RK22E=RK22*(1.D0/(1.D0+RA22*SQRT(RNUE)+RB22*RNUE)+(EPSS*RC22)**2/RB22*RNUE/(1.D0+RC22*RNUE*EPSS))
         RK13E=RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)/(1.D0+RC13*RNUE*EPSS)
         RK23E=RK23/(1.D0+RA23*SQRT(RNUE)+RB23*RNUE)/(1.D0+RC23*RNUE*EPSS)
!         RK33E=RK33/(1.D0+RA33*SQRT(RNUE)+RB33*RNUE)/(1.D0+RC33*RNUE*EPSS)

!         RK2D =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUD)+RB2 *RNUD)+(EPSS*RC2 )**2/RB2 *RNUD/(1.D0+RC2 *RNUD*EPSS))
!         RK2T =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUT)+RB2 *RNUT)+(EPSS*RC2 )**2/RB2 *RNUT/(1.D0+RC2 *RNUT*EPSS))
!         RK2A =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUA)+RB2 *RNUA)+(EPSS*RC2 )**2/RB2 *RNUA/(1.D0+RC2 *RNUA*EPSS))
         RK3D=((1.17D0-0.35D0*SQRT(RNUD))/(1.D0+0.7D0*SQRT(RNUD))-2.1D0*(RNUD*EPSS)**2)/(1.D0+(RNUD*EPSS)**2)
         RK3T=((1.17D0-0.35D0*SQRT(RNUT))/(1.D0+0.7D0*SQRT(RNUT))-2.1D0*(RNUT*EPSS)**2)/(1.D0+(RNUT*EPSS)**2)
         RK3A=((1.17D0-0.35D0*SQRT(RNUA))/(1.D0+0.7D0*SQRT(RNUA))-2.1D0*(RNUA*EPSS)**2)/(1.D0+(RNUA*EPSS)**2)

         DRL=RJCB(NR)/DR
         DTE=(RT(NR+1,1)-RT(NR,1))*DRL
         DTD=(RT(NR+1,2)-RT(NR,2))*DRL
         DTT=(RT(NR+1,3)-RT(NR,3))*DRL
         DTA=(RT(NR+1,4)-RT(NR,4))*DRL
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL
         DPD=(RN(NR+1,2)*RT(NR+1,2)-RN(NR,2)*RT(NR,2))*DRL
         DPT=(RN(NR+1,3)*RT(NR+1,3)-RN(NR,3)*RT(NR,3))*DRL
         DPA=(RN(NR+1,4)*RT(NR+1,4)-RN(NR,4)*RT(NR,4))*DRL
         BPL=BP(NR)

         FACT=1.D0/(1.D0+(RNUE*EPS)**2)
         IF(NSMAX.EQ.2) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
         ELSEIF(NSMAX.EQ.3) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL &
     &       +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)+(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)
         ELSEIF(NSMAX.EQ.4) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL) &
     &       +(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)+(PAL/PEL)*(DPA/PAL-RK3A*FACT*DTA/TAL)
         ENDIF
         H=BB/(BB+BP(NR))
         AJBSL(NR)=-PBSCD*H*SQRT(EPS)*ANE*1.D20*TEL*RKEV/BPL *(RK13E*A+RK23E*DTE/TEL)
!
!        WRITE(6,'(I3,1P6E12.4)') NR,RNUE,RK13E,RK23E,A,BPL,AJBSL(NR)
!         IF(NR.EQ.NRMAX) THEN
!            WRITE(6,'(1P6E12.4)') DN,DTE,ANE,PNSS(1),TE,PTS(1)
!            WRITE(6,'(1P6E12.4)') RK13E*(RN(NR,1)*RT(NR,1)
!     &                   +RN(NR,2)*RT(NR,2)
!     &                   +RN(NR,3)*RT(NR,3)
!     &                   +RN(NR,4)*RT(NR,4)
!     &                   +RW(NR,1)
!     &                   +RW(NR,2)         )*DN/ANE,
!     &            +(RK23E-1.5D0*RK13E)*RN(NR,1)*DTE,
!     &            +RK13D*(RK23D-1.5D0)*RN(NR,2)*DTD,
!     &            BPL,AJBS(NR-1),AJBS(NR)
!         ENDIF
      ENDDO

      NR=NRMAX
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         ANE=PNSS(1)
         ANDX=PNSS(2)
!         ANT =PNSS(3)
!         ANA =PNSS(4)
         TEL=ABS(PTS(1))
         TDL=ABS(PTS(2))
         TTL=ABS(PTS(3))
         TAL=ABS(PTS(4))
         PEL=PNSS(1)*PTS(1)
         PDL=PNSS(2)*PTS(2)
         PTL=PNSS(3)*PTS(3)
         PAL=PNSS(4)*PTS(4)
         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)

         TAUE = FTAUE(ANE,ANDX,TEL,ZEFFL)
         TAUD = FTAUI(ANE,ANDX,TDL,PZ(2),PM(2))
         TAUT = FTAUI(ANE,ANT ,TTL,PZ(3),PM(3))
         TAUA = FTAUI(ANE,ANA ,TAL,PZ(4),PM(4))

         VTE=SQRT(TEL*RKEV/AME)
         VTD=SQRT(TDL*RKEV/AMD)
         VTT=SQRT(TTL*RKEV/AMT)
         VTA=SQRT(TAL*RKEV/AMA)

         RNUE=ABS(QP(NR))*RR/(TAUE*VTE*EPSS)
         RNUD=ABS(QP(NR))*RR/(TAUD*VTD*EPSS)
         RNUT=ABS(QP(NR))*RR/(TAUT*VTT*EPSS)
         RNUA=ABS(QP(NR))*RR/(TAUA*VTA*EPSS)

!         RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)+(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
!         RK12E=RK12*(1.D0/(1.D0+RA12*SQRT(RNUE)+RB12*RNUE)+(EPSS*RC12)**2/RB12*RNUE/(1.D0+RC12*RNUE*EPSS))
!         RK22E=RK22*(1.D0/(1.D0+RA22*SQRT(RNUE)+RB22*RNUE)+(EPSS*RC22)**2/RB22*RNUE/(1.D0+RC22*RNUE*EPSS))
         RK13E=RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE) /(1.D0+RC13*RNUE*EPSS)
         RK23E=RK23/(1.D0+RA23*SQRT(RNUE)+RB23*RNUE) /(1.D0+RC23*RNUE*EPSS)
!         RK33E=RK33/(1.D0+RA33*SQRT(RNUE)+RB33*RNUE)/(1.D0+RC33*RNUE*EPSS)

!         RK2D =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUD)+RB2 *RNUD)+(EPSS*RC2 )**2/RB2 *RNUD/(1.D0+RC2 *RNUD*EPSS))
!         RK2T =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUT)+RB2 *RNUT)+(EPSS*RC2 )**2/RB2 *RNUT/(1.D0+RC2 *RNUT*EPSS))
!         RK2A =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUA)+RB2 *RNUA)+(EPSS*RC2 )**2/RB2 *RNUA/(1.D0+RC2 *RNUA*EPSS))
         RK3D=((1.17D0-0.35D0*SQRT(RNUD))/(1.D0+0.7D0*SQRT(RNUD))-2.1D0*(RNUD*EPSS)**2)/(1.D0+(RNUD*EPSS)**2)
         RK3T=((1.17D0-0.35D0*SQRT(RNUT))/(1.D0+0.7D0*SQRT(RNUT))-2.1D0*(RNUT*EPSS)**2)/(1.D0+(RNUT*EPSS)**2)
         RK3A=((1.17D0-0.35D0*SQRT(RNUA))/(1.D0+0.7D0*SQRT(RNUA))-2.1D0*(RNUA*EPSS)**2)/(1.D0+(RNUA*EPSS)**2)

         DRL=RJCB(NR)/DR
         DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DTD=DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DTT=DERIV3P(PTS(3),RT(NR,3),RT(NR-1,3),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DTA=DERIV3P(PTS(4),RT(NR,4),RT(NR-1,4),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPE=DERIV3P(PNSS(1)*PTS(1),RN(NR,1)*RT(NR,1),RN(NR-1,1)*RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPD=DERIV3P(PNSS(2)*PTS(2),RN(NR,2)*RT(NR,2),RN(NR-1,2)*RT(NR-1,2),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPT=DERIV3P(PNSS(3)*PTS(3),RN(NR,3)*RT(NR,3),RN(NR-1,3)*RT(NR-1,3),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPA=DERIV3P(PNSS(4)*PTS(4),RN(NR,4)*RT(NR,4),RN(NR-1,4)*RT(NR-1,4),RHOG(NR),RHOM(NR),RHOM(NR-1))
         BPL=BP(NR)

         FACT=1.D0/(1.D0+(RNUE*EPS)**2)
         IF(NSMAX.EQ.2) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
         ELSEIF(NSMAX.EQ.3) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL &
     &       +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)+(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)
         ELSEIF(NSMAX.EQ.4) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL) &
     &       +(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)+(PAL/PEL)*(DPA/PAL-RK3A*FACT*DTA/TAL)
         ENDIF
         H=BB/(BB+BP(NR))
         AJBSL(NR)=-PBSCD*H*SQRT(EPS)*ANE*1.D20*TEL*RKEV/BPL*(RK13E*A+RK23E*DTE/TEL)
!
!        WRITE(6,'(I3,1P6E12.4)') NR,RNUE,RK13E,RK23E,A,BPL,AJBSL(NR)
!         IF(NR.EQ.NRMAX) THEN
!            WRITE(6,'(1P6E12.4)') DN,DTE,ANE,PNSS(1),TE,PTS(1)
!            WRITE(6,'(1P6E12.4)') RK13E*(RN(NR,1)*RT(NR,1)
!     &                   +RN(NR,2)*RT(NR,2)
!     &                   +RN(NR,3)*RT(NR,3)
!     &                   +RN(NR,4)*RT(NR,4)
!     &                   +RW(NR,1)
!     &                   +RW(NR,2)         )*DN/ANE,
!     &            +(RK23E-1.5D0*RK13E)*RN(NR,1)*DTE,
!     &            +RK13D*(RK23D-1.5D0)*RN(NR,2)*DTD,
!     &            BPL,AJBS(NR-1),AJBS(NR)
!         ENDIF

      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO

      RETURN
      END SUBROUTINE TRAJBS
