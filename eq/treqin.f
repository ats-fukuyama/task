C
C   *** The variables(data) necessary for this calculation ***
C
C     All of these variables below must be given by TASK/TR.
C
C     NTRMAX        : Maximum array number
C     PSILL(NTRMAX) : Poloidal flux                  (m^2kgs^-2A^-2)
C     HJ(NTRMAX)    : Plasma current density                    (MA)
C     HJTOT(NTRMAX) : Integrated plasma current density
C     OMG(NTRMAX)   : Toroidal angular speed                  (s^-1)
C     HN(NTRMAX)    : Number density                          (m^-3)
C     T(NTRMAX)     : Temperature                              (keV)
C
C   *** Complement ***
C
C     P(PSIN)    : Plasma pressure (P=BLTZ*HN*T; T(K))          (Pa)
C     OMGT(PSIN) : OMG^2/T (OMGT=OMG**2/T; T(K))          (K^-1s^-2)
C
C   ***************************************************************
C
C   ************************************************  
C   **    Caluculate initial PSI from TASK/TR     **
C   ************************************************
C
      SUBROUTINE TRPSIN(IERR)
C              
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      DIMENSION DERIV(NPSM)
C
      NMAX = 32
C
C   *** from PSILL(NTR) to PSI(NTG,NSG) ***
C
      PSIMIN=0.D0
      DO NTR=1,NTRMAX
         IF (PSILL(NTR).LE.PSIMIN) PSIMIN=PSILL(NTR)
      ENDDO
      PSI0 = PSIMIN
C
      CALL SPL1D(PSILL,PSIT,DERIV,UPSIT,NMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRPSIN: SPL1D ERROR : IERR=',IERR
C
      DO NC=1,NMAX
      CALL SPL1DF(NC,PSIT,PSILL,UPSIT,NMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRPSIN: SPL1DF ERROR : IERR=',IERR
      ENDDO
C
      RAXIS=RR
      ZAXIS=0.0D0
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         PSI(NTG,NSG) = PSIT(NTG)
         DELPSI(NTG,NSG) = 0.D0
      ENDDO
      ENDDO
C
C   *** Interpolation ***    
C
C     *** Calculate Pressure and OMEGA**2/TEMP ***
C
C     A unit of T is still [keV], so P[keV/m^3] and OMGT[s^-2keV^-1].
C
      DO NTR=1,NTRMAX
         P(NTR) = HN(NTR)*T(NTR)
         OMGT(NTR) = OMG(NTR)**2/T(NTR)
      ENDDO
C
C     *** Calculate SPLINE coefficients ***
C
      CALL SPL1D(PSILL,HJ,DERIV,UHJP,NTRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRRHSV: SPL1D ERROR1 : IERR=',IERR
      CALL SPL1D(PSILL,P,DERIV,UPP,NTRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRRHSV: SPL1D ERROR2 : IERR=',IERR
      CALL SPL1D(PSILL,OMGT,DERIV,UOMGTP,NTRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRRHSV: SPL1D ERROR3 : IERR=',IERR
      CALL SPL1D(PSILL,T,DERIV,UTP,NTRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRRHSV: SPL1D ERROR4 : IERR=',IERR
      CALL SPL1D(PSILL,HJTOT,DERIV,UHJTOT,NTRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRRHSV: SPL1D ERROR5 : IERR=',IERR
C
      RETURN
      END
C
C   ********************************************************************
C
C   ************************************************  
C   **             Instead of EQRHSV              **
C   ************************************************
C
      SUBROUTINE TRRHSV(IERR)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      DIMENSION PSIRG(NRGM,NZGM),PSIZG(NRGM,NZGM),PSIRZG(NRGM,NZGM)
      EXTERNAL EQPSID
C
      IERR=0
      CALL EQTORZ
      CALL SPL2D(RG,ZG,PSIRZ,PSIRG,PSIZG,PSIRZG,URZ,
     &           NRGM,NRGMAX,NZGMAX,0,0,IERR)
C
      DELT=1.D-8
      EPS=1.D-8
      ILMAX=40
      LIST=0
      RINIT=RAXIS
      ZINIT=ZAXIS
      CALL NEWTN(EQPSID,RINIT,ZINIT,RAXIS,ZAXIS,
     &     DELT,EPS,ILMAX,LIST,IER)
      IF(IER.NE.0) THEN
         WRITE(6,'(A,I5)') 'XX EQRHSV: NEWTN ERROR: IER=',IER
         IERR=101
         RETURN
      ENDIF
      IF(RAXIS.LE.RR+RB) THEN
         PSI0=PSIG(RAXIS,ZAXIS)
      ELSE
         WRITE(6,'(A)') 'XX EQRHSV: AXIS OUT OF PLASMA:'
         IERR=102
         RETURN
      ENDIF
C      PSIMIN=PSI(1,1)
C      PSIMAX=PSI(1,1)
C      DO NSG=1,NSGMAX
C      DO NTG=1,NTGMAX
C         IF(PSI(NTG,NSG).LT.PSIMIN) PSIMIN=PSI(NTG,NSG)
C         IF(PSI(NTG,NSG).GT.PSIMAX) PSIMAX=PSI(NTG,NSG)
C      ENDDO
C      ENDDO
C      PSIMIN=PSIMIN/PSI0
C      PSIMAX=PSIMAX/PSI0
C      IF(MAX(ABS(PSIMIN),ABS(PSIMAX)).GT.3.D0*ABS(PSI0)) THEN
C         WRITE(6,'(A)') 'XX EQRHSV: PSI OUT OF RANGE:'
C         WRITE(6,'(A,1P3E12.4)') 
C     &        '  PSIMIN,PSIMAX,PSI0=',PSIMIN,PSIMAX,PSI0
C         IERR=103
C         RETURN
C      ENDIF
C
C     *** Calculate HJT  ***
C
C     ----- positive current density, jp.gt.0-----
      RRC=RR-RA
C     ----- quasi-symmetric current density, jp:anti-symmetric -----
C      RRC=RR
C
      FJP=0.D0
      FJT=0.D0
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         PSIN=PSI(NTG,NSG)/PSI0
         PP(NTG,NSG)=PPSI(PSIN)
         RMM(NTG,NSG)=RR+SIGM(NSG)*RHOM(NTG)*COS(THGM(NTG))
         ZMM(NTG,NSG)=SIGM(NSG)*RHOM(NTG)*SIN(THGM(NTG))
         HJP1(NTG,NSG)=RMM(NTG,NSG)*DP(PSIN)
         HJT1(NTG,NSG)=HJ(PSIN)
         HJP2A=EXP(RMM(NTG,NSG)**2*OMGT(PSIN)/(2.D0*RGAS))
         HJP2B=EXP(RRC**2*OMGT(PSIN)/(2.D0*RGAS))
         HJP2C=HJP2A-(RRC**2/RMM(NTG,NSG)**2)*HJP2B
         HJP2D=HJP2A-(RRC**4/RMM(NTG,NSG)**4)*HJP2B
         HJP2E=0.5D0*P(PSIN)*RMM(NTG,NSG)**3
         HJP2F=(1/RGAS)*DOMGT(PSIN)
         HJP2(NTG,NSG)=HJP2C*HJP1(NTG,NSG)+HJP2D*HJP2E*HJP2F
         HJT2(NTG,NSG)=(RRC/RMM(NTG,NSG))*HJT1(NTG,NSG)
         DVOL=SIGM(NSG)*RHOM(NTG)*RHOM(NTG)*DSG*DTG
         FJP=FJP+HJP2(NTG,NSG)*DVOL
         FJT=FJT+HJT2(NTG,NSG)*DVOL
      ENDDO
      ENDDO
C     
      TJ=(-RIP*1.D6-FJP)/FJT
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         HJT(NTG,NSG)=HJP2(NTG,NSG)+TJ*HJT2(NTG,NSG)
         PSIN=PSI(NTG,NSG)/PSI0
         TT(NTG,NSG)=SQRT(BB**2*RR**2
     &        +2.D0*RMU0*RRC
     &        *(TJ*HJTOT(PSIN)-RRC*P(PSIN)
     &        *EXP(RRC**2*OMGT(PSIN)/(2.D0*RGAS))))
         RHO(NTG,NSG)=(P(PSIN)/(RGAS*(T(PSIN))))
     &        *EXP(RMM(NTG,NSG)**2*OMGT(PSIN)
     &        /(2.D0*RGAS))
      ENDDO
      ENDDO
C     
      RETURN
      END
C
C   ********************************************************************
C
C   ************************************************
C   **      CALCULATE pp,tt,temp,omega,rho        **
C   ************************************************
C
      SUBROUTINE TREQCP
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      DPS=PSI0/(NPSMAX-1)
      DO NPS=1,NPSMAX
         PSIPS(NPS)=DPS*(NPS-1)
         PSIN=PSIPS(NPS)/PSI0
         PPPS(NPS)=P(PSIN)
         TTPS(NPS)=SQRT(BB**2*RR**2
     &        +2.D0*RMU0*RRC
     &        *(TJ*HJTOT(PSIN)-RRC*P(PSIN)
     &        *EXP(RRC**2*OMGT(PSIN)/(2.D0*RGAS))))
         TEPS(NPS)=T(PSIN)*BLTZ/(AEE*1.D3)
         OMPS(NPS)=SQRT(OMGT(PSIN)*T(PSIN))
      ENDDO
      RETURN
      END
C
C   ***************************************************************
C
C   ***********************************************
C   *                  FUNCTIONS                  *
C   ***********************************************
C
C   *** HJ & HJTOT ***
C
      REAL*8 FUNCTION HJ(PSIN1)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      CALL SPL1DD(PSIN,HJ,DHJ,PSILL,UHJP,NTRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRRHSV: SPL1DD HJ : IERR=',IERR
      HJ=HJ*1.D6
C
      RETURN
      END
C
      REAL*8 FUNCTION HJTOT(PSIN1)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      CALL SPL1DF(PSIN,HJTOT,PSILL,UHJTOT,NTRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQIN: SPL1DF HJTOT : IERR=',IERR
      HJTOT=HJTOT*1.D6
C
      RETURN
      END
C
C   *** P & DP ***
C
      REAL*8 FUNCTION P(PSIN1)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      CALL SPL1DD(PSIN,P,DP,PSILL,UPP,NTRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRRHSV: SPL1DD P : IERR=',IERR
C
C   Now, a unit of P changes [keV/m^3] into [Pa]=[kg m^-1s^-2].
C
      P=P*1.D3*AEE
C
      RETURN
      END
C
C
      REAL*8 FUNCTION DP(PSIN1)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      CALL SPL1DD(PSIN,P,DP,PSILL,UPP,NTRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRRHSV: SPL1DD DP : IERR=',IERR
      DP=DP*1.D3*AEE/PSI0
C
      RETURN
      END
C
C   *** OMGT & DOMGT ***
C
      REAL*8 FUNCTION OMGT(PSIN1)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      CALL SPL1DD(PSIN,OMGT,DOMGT,PSILL,UOMGTP,NTRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRRHSV: SPL1DD OMGT : IERR=',IERR
C
C   Now, a unit of OMGT changes [s^-2keV^-1] into [K^-1s^-2].
C
      OMGT=OMGT*(BLTZ/(1.D3*AEE))
C
      RETURN
      END
C
      REAL*8 FUNCTION DOMGT(PSIN1)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      CALL SPL1DD(PSIN,OMGT,DOMGT,PSILL,UOMGTP,NTRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRRHSV: SPL1DD DOMGT : IERR=',IERR
      DOMGT=DOMGT*(BLTZ/(1.D3*AEE))
C
      RETURN
      END
C
C   *** T ***
C
      REAL*8 FUNCTION T(PSIN1)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      CALL SPL1DD(PSIN,T,DT,PSILL,UTP,NTRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRRHSV: SPL1DD T : IERR=',IERR
C
C   Now, a unit of T changes [keV] into [K].
C
      T=T*1.D3*AEE/BLTZ
C
      RETURN
      END
