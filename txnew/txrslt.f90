!     $Id: txrslt.f90,v 1.44 2010/04/06 05:54:57 honda Exp $
!***********************************************************
!
!           CALCULATE GLOBAL QUANTITIES
!
!***********************************************************

SUBROUTINE TXGLOB

  use tx_commons
  use tx_interface, only : dfdx
  use tx_core_module, only : intg_vol, intg_area, intg_vol_p, sub_intg_vol
  implicit none

  INTEGER(4) :: I, NR!, NS, NF
  real(8), parameter :: toMega = 1.d-6
  REAL(8) :: RNINT, RPINT, RPEINT, RPIINT, ANFINT, RWINT, &
       &     PAI, BPave, denomINT, volume!, FKAP = 1.d0
  REAL(8) :: SUMdenom, SUMPNiV!, SUMM, SUMP
  REAL(8), DIMENSION(1:NRMAX) :: BETA, BETAP, BETAL, BETAPL, BETAQ
  real(8), dimension(0:NRMAX) :: BthV2, PNdiff
  real(8), dimension(:), allocatable :: denom, dPNV, RSmb
  real(8) :: suml, sum1, sum2
  real(8) :: DERIV4, FCTR

  !     Volume-Averaged Density and Temperature
  !     Core Density and Temperature
  !     Electoron and ion Work Quantities

  !  Plasma volume
  volume = vlt(nrmax)

  PNdiff(:) =  achg(2) * Var(:,2)%n + achgb * PNbV(:) &
       &           + achgb * rip_rat(:) * PNbrpV(:) - Var(:,1)%n
  VOLAVN =  intg_vol(PNdiff) / volume

  ! Electron
  RNINT  = intg_vol(Var(:,1)%n)
  RPEINT = intg_vol(Var(:,1)%p)
  ANSAV(1) = RNINT / volume
  ANS0(1)  = Var(0,1)%n
  IF(RNINT > 0.D0) THEN
     TSAV(1) = RPEINT/RNINT ! [keV]
  ELSE
     TSAV(1)=0.D0
  END IF
  TS0(1) = Var(0,1)%T
  WST(1) = 1.5D0*RPEINT*rKeV*1.D20*toMega ! [MJ]

  ! Ion
  RNINT  = intg_vol(Var(:,2)%n)
  RPIINT = intg_vol(Var(:,2)%p)
  ANSAV(2) = RNINT / volume
  ANS0(2)  = Var(0,2)%n
  IF(RNINT > 0.D0) THEN
     TSAV(2) = RPIINT/RNINT ! [keV]
  ELSE
     TSAV(2)=0.D0
  END IF
  TS0(2) = Var(0,2)%T
  WST(2) = 1.5D0*RPIINT*rKeV*1.D20*toMega ! [MJ]

  ! Beam ion
  ANFINT = intg_vol(PNbV) + intg_vol(PNbrpV)
  IF(ANFINT > 0.D0) THEN
     RWINT  = intg_vol(SNB/rNuB)
!     WFT(1) = 0.5D0*AMb*amp*RWINT**2*1.5D0*2.D0*PI*RR*2.D0*PI*elip(NRA)*1.D-6!rKeV*1.D14
     WFT(1) = 0.5d0*RWINT*Eb*rKeV*1.d20*toMega ! [MJ]
     ANFAV(1) = ANFINT / volume 
     ANF0(1)  = PNbV(0) + PNbrpV(0)
     TFAV(1)  = 0.5d0*RWINT*Eb/ANFINT ! [keV]
     if(ANF0(1) > 0.D0) TF0(1) = SNB(0)/ANF0(1)
  ELSE
     WFT(1)   = 0.d0
     ANFAV(1) = 0.d0
     ANF0(1)  = 0.d0
     TFAV(1)  = 0.d0
     TF0(1)   = 0.D0
  END IF

  ! Num. of neutrals

  PnumN0 = intg_vol(PN01V+PN02V+PN03V)*1.d20

  !     Input powers

  POHT = intg_vol(POH)*toMega ! [MW]
  PNBT = intg_vol(PNB)*toMega ! [MW]
  PNFT =(intg_vol(PALFe)+intg_vol(PALFi))*toMega ! [MW]

  PRFTe   = intg_vol(PRFe)*toMega
  PRFTi   = intg_vol(PRFi)*toMega
  PRFT    = PRFTe+PRFTi

  !      PFINT = intg_vol(PFIN)*toMega
  !      DO NS=1,NSM
  !        PFCLT(NS) = intg_vol(PFCL(:,NS))*toMega
  !      END DO

  !     Output powers

  SIE(:) = Var(:,1)%n*rNuION(:)*1.D20
  PIE(:) = SIE(:)*EION*AEE
  SCX(:) = FSCX * SiVcxA(:) * Var(:,2)%n * PN01V(:) * 1.D40
  PCX(:) = 1.5D0*Var(:,2)%n*rNuiCXT(:)*1.D20*Var(:,2)%T*rKeV

  !    PRLT = intg_vol(PRL)*toMega
  PCXT = intg_vol(PCX)*toMega
  PIET = intg_vol(PIE)*toMega

  !     Currents

  AJT   = intg_area(AJ  )*toMega
  AJOHT = intg_area(AJOH)*toMega
  AJNBT = intg_area(AJNB)*toMega
  AJBST = intg_area(AJBS)*toMega

  !      DRH=0.5D0*DR
  !      DO NS=1,NSM
  !         VNP=AV(NRMAX,NS)
  !         DNP=AD(NRMAX,NS)
  !
  !         VTP=(AVK(NRMAX,NS)/1.5D0+AV(NRMAX,NS))
  !         DTP= AK(NRMAX,NS)/(1.5D0)
  !
  !         VXP=0.D0
  !         DXP=(1.5D0*AD(NRMAX,NS)-AK(NRMAX,NS))*PTS(NS)
  !
  !         SLT(NS) =((     DNP/DRH)*RN(NRMAX,NS)
  !     &            +( VNP-DNP/DRH)*PNSS(NS))
  !     &            *2.D0*PI*RR*2.D0*PI*RA*FKAP
  !
  !         PLT(NS) =((     DXP/DRH)*RN(NRMAX,NS)
  !     &            +(     DTP/DRH)*RN(NRMAX,NS)*RT(NRMAX,NS)*1.5D0
  !     &            +( VXP-DXP/DRH)*PNSS(NS)
  !     &            +( VTP-DTP/DRH)*PNSS(NS)*PTS(NS)*1.5D0)
  !     &            *2.D0*PI*RR*2.D0*PI*RA*FKAP*rKeV*1.D14
  !      END DO
  !

  SIET = intg_vol(SIE)
  !      SNFT = intg_vol(SNF)
  SNBT = intg_vol(SNB)

  !      DO NS=1,NSM
  !         SPET(NS) = intg_vol(SPE)
  !      END DO

  !     Torque deposition

  allocate(RSmb(0:NRMAX))
  RSmb(:) = fipol(:) / bbt(:) * BSmb(:) * 1.d20
  TNBcol = intg_vol(RSmb)
  deallocate(RSmb)
  TTqt = intg_vol(Tqt)

  WBULKT = SUM(WST(1:2))
  WTAILT = SUM(WFT(1:1))

  WPT  = WBULKT+WTAILT
  PINT = PNBT+PRFT+POHT
  POUT = PCXT+PIET
  SINT = SIET+SNBT
  !    SOUT = SLST

  IF(ABS(T_TX - TPRE) <= epsilon(1.d0)) THEN
     WPDOT = 0.D0
  ELSE
     WPDOT = (WPT - WPPRE)/(T_TX - TPRE)
  END IF
  WPPRE = WPT
  TPRE  = T_TX

  IF(PINT <= 0.D0) THEN
     TAUE1 = 0.D0
     TAUE2 = 0.D0
  ELSE
     TAUE1 = WPT/PINT
     TAUE2 = WPT/(PINT - WPDOT)
  END IF

!     *** Local beta ***
!        BETA  : volume-averaged toroidal beta
!        BETAL : toroidal beta
!        BETAP : volume-averaged poloidal beta
!        BETAPL: poloidal beta
!        BETAQ : toroidal beta for reaction rate
!               (ref. TOKAMAKS 3rd, p115)

  ! Fast ion components are not included currently.

  sum1 = 0.d0 ; sum2 = 0.d0
  DO NR = 1, NRMAX
     suml = sum(Var(NR,1:NSM)%p)*rKeV*1.D20
     do i = 1, NSM
        RPINT = intg_vol_p(Var(:,i)%p,NR)*rKeV*1.D20
        sum1  = sum1 + RPINT
        sum2  = sum2 + RPINT*RPINT
     end do
     BETA  (NR) = 2.D0*rMU0*sum1      /(     vlt(nr) *bb*bb)
     BETAL (NR) = 2.D0*rMU0*suml      /(              bb*bb)
     BETAP (NR) = 2.D0*rMU0*sum1      /(     vlt(nr) *BthV(NRMAX)*BthV(NRMAX))
     BETAPL(NR) = 2.D0*rMU0*suml      /(              BthV(NRMAX)*BthV(NRMAX))
     BETAQ (NR) = 2.D0*rMU0*sqrt(sum2)/(sqrt(vlt(nr))*BthV(NR)*BthV(NR))
  END DO

  BETA0  = FCTR(R(1),R(2),BETA  (1),BETA  (2))
  BETAP0 = FCTR(R(1),R(2),BETAPL(1),BETAPL(2))
  BETAQ0 = FCTR(R(1),R(2),BETAQ (1),BETAQ(2))

!     *** Global beta ***
!        BETAPA: poloidal beta at separatrix
!        BETAA : toroidal beta at separatrix
!        BETAN : normalized toroidal beta (Troyon beta)

  BETAPA = BETAP(NRMAX)
  BETAA  = BETA (NRMAX)
  BETAN  = BETAA / (rIP / (RA * BB)) * 100.D0

  ! Internal inductance li(3)
  ! (In our case, li(1)=li(2)=li(3) due to the circular cross-section.)
  BthV2(:) = BthV(:)*BthV(:)
  BPave = intg_vol(BthV2) / volume
  ALI   = BPave / ((rMU0 * rIp * 1.D6)**2 * RR / (2.D0 * volume))
  VLOOP = 2.d0*Pi*PsidotV(NRMAX)

  !  PAI=(PA(2)*PN(2)+PA(3)*PN(3)+PA(4)*PN(4))/(PN(2)+PN(3)+PN(4))
  PAI=amas(2)

  IF(PINT <= 0.D0) THEN
     TAUEP=0.D0
     TAUEH=0.D0
     QF=0.D0
  ELSE
     TAUEP=4.8D-2*(rIP**0.85D0) &
          &      *(RR**1.2D0) &
          &      *(RA**0.3D0) &
          &      *(elip(NRA)**0.5D0) &
          &      *(ANSAV(1)**0.1D0) &
          &      *(BB**0.2D0) &
          &      *(PAI**0.5D0) &
          &      *(PINT**(-0.5D0))
     TAUEH=0.0562D0*(rIP**0.93D0) &
          &        *(BB**0.15D0) &
          &        *(PINT**(-0.69D0)) &
          &        *(ANSAV(1)**0.41D0) &
          &        *(PAI**0.19D0) &
          &        *(RR**1.97D0) &
          &        *((RA/RR)**0.58D0) &
          &        *(elip(NRA)**0.78D0)
     QF=5.D0*PNFT/PINT
  END IF

  IF(Q(0) >= 1.D0) THEN
     RQ1=0.D0
  ELSE
     OUTER:DO 
        IF(Q(1) > 1.D0) THEN
           RQ1=SQRT( (1.D0-Q(0) )/(Q(1)-Q(0)) )*(R(2)-R(1))
           EXIT OUTER
        END IF
        DO I=2,NRMAX
           IF(Q(I) > 1.D0) THEN
              RQ1=(R(I)-R(I-1))*(1.D0-Q(I-1))/(Q(I)-Q(I-1))+R(I-1)
              EXIT OUTER
           END IF
        END DO
        RQ1=RA
        EXIT OUTER
     END DO OUTER
  END IF

  !  ZEFF0=FCTR(R(1),R(2),ZEFF(1),ZEFF(2))
  ZEFF0=Zeff

!!$  !  Evaluate effective particle diffusivity from particle flux minus Ware pinch
!!$  !  Gamma = n v - Deff dn/dr  --> Deff = (n v - ft n Eph / Bth) / (- dn/dr)
!!$
!!$  Deff(0) = 0.d0
!!$  do nr = 1, nrmax
!!$     ! Including Ware pinch
!!$!     Deff(nr) = (Var(NR,2)%n*Var(NR,2)%Ur-ft(NR)*Var(NR,2)%n*Etor(NR)/BthV(NR)) &
!!$!          & /(-DERIV4(NR,R,Var(:,2)%n,NRMAX,0))
!!$     ! Excluding Ware pinch
!!$     Deff(nr) = -Var(NR,2)%n*Var(NR,2)%Ur/DERIV4(NR,R,Var(:,2)%n,NRMAX,0)
!!$     if(Deff(nr) < 0.d0) Deff(nr) = 0.d0
!!$!     write(6,*) 'nr=',nr,'Deff=',deff(nr)
!!$  end do

  ! Effective ion particle diffusivity
  ! n_s <u_s.grad V> = - Deff <|grad.V|^2> dn_s/dV

  allocate(dPNV(0:NRMAX))
  Deff(0) = 0.d0
  dPNV(:) = dfdx(vv,Var(:,1)%n,NRMAX,0)
  do NR = 1, NRMAX
     Deff(NR) = - Var(NR,1)%n * Var(NR,1)%UrV / (sst(NR) * dPNV(NR))
  end do

  ! *** Particle confinement time; TAUP ***

  allocate(denom(0:NRMAX))
  denom(:) = Var(:,2)%n * rNuION(:)
  denomINT = intg_vol(denom)
  if(denomINT < epsilon(1.d0)) then
     TAUP = 0.d0
  else
     TAUP  = intg_vol(Var(:,2)%n) / denomINT
  end if
  deallocate(dPNV,denom)

  allocate(denom(0:NRA))
  denom(0:NRA) = Var(0:NRA,2)%n * rNuION(0:NRA)
  call sub_intg_vol(Var(:,2)%n,NRA,SUMPNiV)
  call sub_intg_vol(denom,NRA,SUMdenom)
  if(SUMdenom < epsilon(1.d0)) then
     TAUPA = 0.d0
  else
     TAUPA = SUMPNiV / SUMdenom
  end if
  deallocate(denom)

  ! *** Ion outflux through the separatrix ***

  call cal_flux

  RETURN
END SUBROUTINE TXGLOB

!************************************************************
!
!  Outflux of ions through the separatrix in case of no NBI
!
!    Output : Gamma_a [m^-2s^-1]
!
!************************************************************

subroutine cal_flux

  use tx_commons, only : NRA, PI, rNuION, achg, DT, Gamma_a, surflcfs, Var
  use tx_core_module, only : sub_intg_vol
  implicit none

  real(8) :: total_ion, SizINT, SumPNiV
  real(8), save :: total_ion_save = 0.d0
  real(8), dimension(:), allocatable :: Siz

  !  Number of ions inside the separatrix [particles]
  call sub_intg_vol(Var(:,2)%n,NRA,SumPNiV)
  total_ion = SumPNiV*1.D20

  !  Rate of generation of ions inside the separatix
  allocate(Siz(0:NRA))
  Siz(0:NRA) = rNuION(0:NRA)*Var(0:NRA,1)%n/achg(2)*1.D20 ! [m^-3s^-1]
  call sub_intg_vol(Siz,NRA,SizINT)
  SizINT = SizINT ! [s^-1]
  deallocate(Siz)

  !  Outflux of ions through the separatix [m^-2s^-1]
  Gamma_a = (SizINT - (total_ion - total_ion_save) / DT) / surflcfs

  total_ion_save = total_ion

end subroutine cal_flux
