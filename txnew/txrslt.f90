module tx_glob
  implicit none
  public

contains
!***********************************************************
!
!           CALCULATE GLOBAL QUANTITIES
!
!***********************************************************

  subroutine TXGLOB

    use tx_commons, only : Pi, rKeV, AEE, rMU0, EPS0, amp, amb, T_TX, TPRE, NRMAX, NRA, NSM &
         & , Eb, EION, FSCX, rIP, RA, RR, bb &
         & , array_init_NR, rpt, vv, Var, vlt, achg, achgb, amas, rip_rat &
         & , PNbV, PNbrpV, RUbphV, SNB, rNuB, BSmb &
         & , PN01V, PN02V, PN03V, PN0zV, PnumN0, PnumN0z & 
         & , POH, PNB, PALFe, PALFi, PRFe, PRFi, rNuION, SiVcxA, rNuiCXT &
         & , PsidotV, PsitdotV, BthV, Q, AJ, AJOH, AJNB, AJBS, Tqt, Deff, qneut &
         & , fipol, bbt, sst, sdt, elip, Zeff, bthco, gtti, aat, ckt, hdt &
         & , irestart, iPoyntpol, iPoynttor, VLOOPpol, PoyntS, PoyntR, PoyntI, Vpoynt &
         & , CPsi, CPsi_old, CEjima &
         & , VOLAVN, ANSAV, ANS0, TSAV, TS0, WST, WFT, ANFAV, ANF0, TFAV, TF0 &
         & , POHT, PNBT, PNFT, PRFTe, PRFTi, PRFT &
         & , SIE, PIE, SCX, PCX, PCXT, PIET, AJT, AJOHT, AJNBT, AJBST &
         & , SIET, SNBT, TNBcol, TTqt, WBULKT, WTAILT, WPT &
         & , PINT, POUT, SINT, BETA0, BETAP0, BETAQ0, BETAPA, BETAA, BETAN &
         & , TAUP, TAUPA, totmnRV, WPDOT, WPPRE &
         & , TAUE1, TAUE2, TAUEP, TAUEH, RQ1, ZEFF0, ALI, totjxB, QF, VLOOP
    use tx_interface, only : dfdx
    use tx_core_module, only : intg_vol, intg_area, intg_vol_p, sub_intg_vol
    use libitp, only : FCTR
    implicit none

    integer(4) :: I, NR!, NS, NF
    real(8), parameter :: toMega = 1.d-6, fourPisq = 4.d0 * Pi * Pi
    real(8) :: RNINT, RPINT, ANFINT, RWINT, &
         &     PAI, BPave, denomINT, volume!, FKAP = 1.d0
    real(8) :: SUMdenom, SUMPNiV, suml, sum1, sum2
    real(8), dimension(1:NRMAX) :: BETA, BETAP, BETAL, BETAPL, BETAQ
    real(8), dimension(0:NRMAX) :: BthV2
    real(8), dimension(:), allocatable :: denom, dPNV, RSmb
    ! ---------------------------------

    !     Volume-Averaged Density and Temperature
    !     Core Density and Temperature
    !     Electron and ion Work Quantities

    !  Plasma volume
    volume = vlt(NRMAX)

    do NR = 0, NRMAX
       qneut(NR) = sum(achg(2:NSM) * Var(NR,2:NSM)%n) - Var(NR,1)%n &
            &     + achgb * PNbV(NR) + achgb * rip_rat(NR) * PNbrpV(NR)
    end do
    VOLAVN =  intg_vol(qneut) / volume

    !  Volume-averaged density & temperature, stored energy for thermal species
    do i = 1, NSM
       RNINT = intg_vol(Var(:,i)%n)
       RPINT = intg_vol(Var(:,i)%p)
       ANSAV(i) = RNINT / volume
       ANS0(i)  = Var(0,i)%n
       if(RNINT > 0.d0) then
          TSAV(i) = RPINT/RNINT ! [keV]
       else
          TSAV(i) = 0.d0
       end if
       TS0(i) = Var(0,i)%T
       WST(i) = 1.5d0*RPINT*rKeV*1.d20*toMega ! [MJ]
    end do

    ! Beam ion
    ANFINT = intg_vol(PNbV) + intg_vol(PNbrpV)
    if(ANFINT > 0.d0) then
       RWINT  = intg_vol(SNB/rNuB)
       !     WFT(1) = 0.5d0*AMb*amp*RWINT**2*1.5d0*2.d0*PI*RR*2.d0*PI*elip(NRA)*1.D-6!rKeV*1.D14
       WFT(1) = 0.5d0*RWINT*Eb*rKeV*1.d20*toMega ! [MJ]
       ANFAV(1) = ANFINT / volume 
       ANF0(1)  = PNbV(0) + PNbrpV(0)
       TFAV(1)  = 0.5d0*RWINT*Eb/ANFINT ! [keV]
       if(ANF0(1) > 0.d0) TF0(1) = SNB(0)/ANF0(1)
    else
       WFT(1)   = 0.d0
       ANFAV(1) = 0.d0
       ANF0(1)  = 0.d0
       TFAV(1)  = 0.d0
       TF0(1)   = 0.d0
    end if

    ! Num. of neutrals

    PnumN0(1) = intg_vol(PN01V)*1.d20
    PnumN0(2) = intg_vol(PN02V)*1.d20
    PnumN0(3) = intg_vol(PN03V)*1.d20
    PnumN0(0) = sum(PnumN0(1:3)) 

    PnumN0z   = intg_vol(PN0zV)*1.d20

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
    PCX(:) = 1.5d0*Var(:,2)%n*rNuiCXT(:)*1.D20*Var(:,2)%T*rKeV

    !    PRLT = intg_vol(PRL)*toMega
    PCXT = intg_vol(PCX)*toMega
    PIET = intg_vol(PIE)*toMega

    !     Currents

    AJT   = intg_area(AJ  )*toMega
    AJOHT = intg_area(AJOH)*toMega
    AJNBT = intg_area(AJNB)*toMega
    AJBST = intg_area(AJBS)*toMega

    !      DRH=0.5d0*DR
    !      DO NS=1,NSM
    !         VNP=AV(NRMAX,NS)
    !         DNP=AD(NRMAX,NS)
    !
    !         VTP=(AVK(NRMAX,NS)/1.5d0+AV(NRMAX,NS))
    !         DTP= AK(NRMAX,NS)/(1.5d0)
    !
    !         VXP=0.d0
    !         DXP=(1.5d0*AD(NRMAX,NS)-AK(NRMAX,NS))*PTS(NS)
    !
    !         SLT(NS) =((     DNP/DRH)*RN(NRMAX,NS)
    !     &            +( VNP-DNP/DRH)*PNSS(NS))
    !     &            *2.d0*PI*RR*2.d0*PI*RA*FKAP
    !
    !         PLT(NS) =((     DXP/DRH)*RN(NRMAX,NS)
    !     &            +(     DTP/DRH)*RN(NRMAX,NS)*RT(NRMAX,NS)*1.5d0
    !     &            +( VXP-DXP/DRH)*PNSS(NS)
    !     &            +( VTP-DTP/DRH)*PNSS(NS)*PTS(NS)*1.5d0)
    !     &            *2.d0*PI*RR*2.d0*PI*RA*FKAP*rKeV*1.D14
    !      END DO
    !

    SIET = intg_vol(SIE)
    !      SNFT = intg_vol(SNF)
    SNBT = intg_vol(SNB)

    !      DO NS=1,NSM
    !         SPET(NS) = intg_vol(SPE)
    !      END DO

    !     Torque deposition

    allocate(RSmb, mold=fipol)
    RSmb(:) = fipol(:) / bbt(:) * BSmb(:) * 1.d20
    TNBcol = intg_vol(RSmb)
    deallocate(RSmb)
    TTqt = intg_vol(Tqt)

    WBULKT = sum(WST(:))
    WTAILT = sum(WFT(:))

    WPT  = WBULKT+WTAILT
    PINT = PNBT+PRFT+POHT
    POUT = PCXT+PIET
    SINT = SIET+SNBT
    !    SOUT = SLST

    if(abs(T_TX - TPRE) <= epsilon(1.d0)) then
       WPDOT = 0.d0
    else
       WPDOT = (WPT - WPPRE)/(T_TX - TPRE)
    end if
    WPPRE = WPT
    TPRE  = T_TX

    if(PINT <= 0.d0) then
       TAUE1 = 0.d0
       TAUE2 = 0.d0
    else
       TAUE1 = WPT/PINT
       TAUE2 = WPT/(PINT - WPDOT)
    end if

    !     *** Local beta ***
    !        BETA  : volume-averaged toroidal beta
    !        BETAL : toroidal beta
    !        BETAP : volume-averaged poloidal beta
    !        BETAPL: poloidal beta
    !        BETAQ : toroidal beta for reaction rate
    !               (ref. TOKAMAKS 3rd, p115)

    ! Fast ion components are not included currently.

    sum1 = 0.d0 ; sum2 = 0.d0
    do NR = 1, NRMAX
       suml = sum(Var(NR,1:NSM)%p)*rKeV*1.D20
       do i = 1, NSM
          RPINT = intg_vol_p(Var(:,i)%p,NR)*rKeV*1.D20
          sum1  = sum1 + RPINT
          sum2  = sum2 + RPINT*RPINT
       end do
       BETA  (NR) = 2.d0*rMU0*sum1      /(     vlt(nr) *bb*bb)
       BETAL (NR) = 2.d0*rMU0*suml      /(              bb*bb)
       BETAP (NR) = 2.d0*rMU0*sum1      /(     vlt(nr) *BthV(NRMAX)*BthV(NRMAX))
       BETAPL(NR) = 2.d0*rMU0*suml      /(              BthV(NRMAX)*BthV(NRMAX))
       BETAQ (NR) = 2.d0*rMU0*sqrt(sum2)/(sqrt(vlt(nr))*BthV(NR)*BthV(NR))
    end do

    BETA0  = FCTR(rpt(1),rpt(2),BETA  (1),BETA  (2))
    BETAP0 = FCTR(rpt(1),rpt(2),BETAPL(1),BETAPL(2))
    BETAQ0 = FCTR(rpt(1),rpt(2),BETAQ (1),BETAQ(2))

    !     *** Global beta ***
    !        BETAPA: poloidal beta at separatrix
    !        BETAA : toroidal beta at separatrix
    !        BETAN : normalized toroidal beta (Troyon beta)

    BETAPA = BETAP(NRMAX)
    BETAA  = BETA (NRMAX)
    BETAN  = BETAA / (rIP / (RA * BB)) * 100.d0

    ! Internal inductance li(3)
    ! (li(1)=li(2)=li(3) in the case of the circular cross-section.)
    BthV2(:) = BthV(:)*BthV(:)
    BPave = intg_vol(BthV2) / volume
    ALI   = BPave / ((rMU0 * rIp * 1.d6)**2 * RR / (2.d0 * volume))
    VLOOP = 2.d0*Pi*PsidotV(NRMAX)

    !  PAI=(PA(2)*PN(2)+PA(3)*PN(3)+PA(4)*PN(4))/(PN(2)+PN(3)+PN(4))
    PAI=amas(2)

    if(PINT <= 0.d0) then
       TAUEP=0.d0
       TAUEH=0.d0
       QF=0.d0
    else
       TAUEP=4.8D-2*(rIP**0.85d0) &
            &      *(RR**1.2d0) &
            &      *(RA**0.3d0) &
            &      *(elip(NRA)**0.5d0) &
            &      *(ANSAV(1)**0.1d0) &
            &      *(BB**0.2d0) &
            &      *(PAI**0.5d0) &
            &      *(PINT**(-0.5d0))
       TAUEH=0.0562d0*(rIP**0.93d0) &
            &        *(BB**0.15d0) &
            &        *(PINT**(-0.69d0)) &
            &        *(ANSAV(1)**0.41d0) &
            &        *(PAI**0.19d0) &
            &        *(RR**1.97d0) &
            &        *((RA/RR)**0.58d0) &
            &        *(elip(NRA)**0.78d0)
       QF=5.d0*PNFT/PINT
    end if

    if(Q(0) >= 1.d0) then
       RQ1=0.d0
    else
       OUTER:do 
          if(Q(1) > 1.d0) then
             RQ1=sqrt( (1.d0-Q(0) )/(Q(1)-Q(0)) )*(rpt(2)-rpt(1))
             exit OUTER
          end if
          do I=2,NRMAX
             if(Q(I) > 1.d0) then
                RQ1=(rpt(I)-rpt(I-1))*(1.d0-Q(I-1))/(Q(I)-Q(I-1))+rpt(I-1)
                exit OUTER
             end if
          end do
          RQ1=RA
          exit OUTER
       end do OUTER
    end if

    ZEFF0=Zeff(0)

    ! Effective ion particle diffusivity
    ! n_s <u_s.grad V> = - Deff <|grad.V|^2> dn_s/dV

    allocate(dPNV, mold=vv)
    Deff(0,:) = 0.d0
    do i = 1, NSM
       dPNV(:) = dfdx(vv,Var(:,i)%n,NRMAX,0)
       do NR = 1, NRMAX
          Deff(NR,i) = - Var(NR,i)%n * Var(NR,i)%UrV / (sst(NR) * dPNV(NR))
       end do
    end do

    ! *** Particle confinement time; TAUP ***

    allocate(denom, mold=rNuION)
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

    ! *** Conservation of toroidal angular momentum ***

    totmnRV = 0.d0
    do i = 1, NSM
       ! thermal species
       totmnRV = totmnRV + amas(i) * amp * intg_vol(Var(:,i)%n * Var(:,i)%RUph) * 1.d20
    end do
    ! fast particles
    totmnRV = totmnRV + amb * amp * intg_vol(PNbV(:) * RUbphV(:)) * 1.d20

    ! *** Volume integrated jxB torque ***

    totjxB = 0.d0
    do i = 1, NSM
       totjxB = totjxB + aee * 1.d20 &
            & * intg_vol(sdt(:) * achg(i) * Var(:,i)%n * Var(:,i)%UrV) ! <j.grad psi>
    end do

    ! ----------

    block
      ! --- Poynting flux calculation ---
      !  << Details are shown in [M. Honda, Nucl. Fusion 58 (2018) 026006] >>
      !    iPoyntpol and iPoynttor are flags to choose Poynting flux equations:
      !    iPoyntpol = 1 and iPoynttor = 0  ==>  eq. (35) in the reference
      !    iPoyntpol = 0 and iPoynttor = 1  ==>  eq. (36) in the reference
      !    iPoyntpol = 1 and iPoynttor = 1  ==>  eq. (25) in the reference
      real(8), dimension(:), allocatable :: EngyEp, EngyEt, EngyMp, EngyMt, EngyJP, EngyJT
      real(8) :: WJengy, WFEpengy, WFEtengy, WFEengy, WFMpengy, WFMtengy, WFMengy, &
           &     jEpol, jEtor, zIp, zIpol, zvar, zvar_old
      real(8), save :: T_TX_old, zIp_old, WJengy_old, WFEengy_old, WFMengy_old, &
           &           WFMpengy_old, WFMtengy_old, VLOOP_old, PS2_old, zIpol_old, VLOOPpol_old

      ! *** Poynting flux and flux consumption ***

      zIp = rIp * 1.d6 ! [MA] -> [A]

      !   --- Poloidal current around the plasma surface ---

      zIpol = 2.d0 * Pi * fipol(NRMAX) / rMU0

      !   --- Poloidal (toroidal flux) contribution to loop voltage ---

      VLOOPpol =-2.d0*Pi*PsitdotV(NRMAX) * iPoynttor

      !   --- Poynting flux (LHS) ---

      PoyntS(1) = - VLOOP    * zIp

      PoyntS(2) =   VLOOPpol * zIpol

      !   --- Energy of joule heating : WJengy ---

      allocate(EngyJP, EngyJT, EngyEp, EngyEt, EngyMp, EngyMt, mold=array_init_NR)
      do nr = 0, nrmax
         jEpol = 0.d0
         jEtor = 0.d0
         do i = 1, NSM
            jEpol = jEpol + achg(i) * Var(nr,i)%n * 1.d20 * Var(nr,i)%Uthhat
            jEtor = jEtor + achg(i) * Var(nr,i)%n * 1.d20 * Var(nr,i)%UphR
         end do
         EngyJP(nr) = - PsitdotV(nr) * bthco(nr) * aee * jEpol * iPoyntpol
         EngyJT(nr) =   PsidotV (nr)             * aee * jEtor * iPoynttor
      end do
      ! Poynting flux from toroidal component of particle acceleration (joule dissipation)
      PoyntR(1) = - intg_vol(EngyJT)
      ! Poynting flux from poloidal component of particle acceleration (joule dissipation)
      PoyntR(2) = - intg_vol(EngyJP)
      WJengy = - sum(PoyntR(1:2))

      !   --- Energy of electromagnetic fields : WFengy ---

      do nr = 0, nrmax
         ! EngyE : induced electric energy density
         EngyEp(nr) = 0.5d0 * EPS0 * gtti(nr) * PsitdotV(nr)**2 ! poloidal
         EngyEt(nr) = 0.5d0 * EPS0 * aat (nr) * PsidotV (nr)**2 ! toroidal
         ! EngyM : magnetic energy density
         EngyMp(nr) = 0.5d0 / rMU0 * ckt(NR) * sdt(NR)*sdt(NR)                     ! poloidal
         EngyMt(nr) = 0.5d0 / rMU0 * fourPisq*fourPisq * hdt(NR)*hdt(NR) / aat(NR) ! toroidal
      end do
      WFEpengy   = intg_vol(EngyEp)
      WFEtengy   = intg_vol(EngyEt)
      WFEengy    =(WFEpengy + WFEtengy) * iPoynttor ! Electric field energy
      WFMpengy   = intg_vol(EngyMp)     * iPoynttor
      WFMtengy   = intg_vol(EngyMt)     * iPoyntpol
      WFMengy    = WFMpengy + WFMtengy              ! Magnetic field energy
      
      !   --- The inductive and resistive parts of the magnetic flux change at the plasma surface ---

      if( T_TX /= 0.d0 .and. irestart == 0 ) then
         ! Poynting flux from electric energy
         PoyntI(1) = - (WFEengy - WFEengy_old) / (T_TX - T_TX_old)
         ! Poynting flux from magnetic energy : PoyntI(2) = PoyntI(3) + PoyntI(4)
         PoyntI(2) = - (WFMengy - WFMengy_old) / (T_TX - T_TX_old)
         ! Poynting flux from poloidal magnetic energy
         PoyntI(3) = - (WFMpengy - WFMpengy_old) / (T_TX - T_TX_old)
         ! Poynting flux from toroidal magnetic energy
         PoyntI(4) = - (WFMtengy - WFMtengy_old) / (T_TX - T_TX_old)

         ! Loop voltage
         VPoynt(0) = VLOOP
         !  ** Integrate "Vloop" w.r.t. time (CPsi : Capital Psi) ; eq. (28) on LHS
         if( iPoyntpol == 1 .and. iPoynttor == 0 ) then
            zvar     = VLOOPpol
            zvar_old = VLOOPpol_old
         else
            zvar     = VLOOP
            zvar_old = VLOOP_old
         end if
         CPsi(0) = CPsi_old(0)  &
              &   + 0.5d0 * (zvar + zvar_old) * (T_TX - T_TX_old)

         ! Inductive flux
         if( iPoyntpol == 1 .and. iPoynttor == 0 ) then
            zvar     = zIpol
            zvar_old = zIpol_old
         else
            zvar     = zIp
            zvar_old = zIp_old
         end if
         VPoynt(1) = ((WFEengy + WFMengy) - (WFEengy_old + WFMengy_old)) / (T_TX - T_TX_old) &
              &  / ( 0.5d0 * (zvar + zvar_old) )
         !          &  / zvar
         !  ** Integrate "(dWF/dt)/Ip" w.r.t. time ; 1st term on RHS of eq. (28) or of eq. (40)
         CPsi(1) = CPsi_old(1) &
              &   + 0.5d0 * (1.d0 / zvar + 1.d0 / zvar_old) &
              &   * ((WFEengy + WFMengy) - (WFEengy_old + WFMengy_old))
         ! Resistive flux
         VPoynt(2) = WJengy / zvar
         !  ** Integrate "WJ/Ip" w.r.t. time ; 2nd term on RHS of eq. (28) or of eq. (40)
         if( iPoyntpol == 1 .and. iPoynttor == 0 ) then
            zvar     = WJengy     / zIpol
            zvar_old = WJengy_old / zIpol_old
         else
            zvar     = WJengy     / zIp
            zvar_old = WJengy_old / zIp_old
         end if
         CPsi(2) = CPsi_old(2) &
              &   + 0.5d0 * (zvar + zvar_old) * (T_TX - T_TX_old)

         ! Poloidal loop voltage
         VPoynt(3) = VLOOPpol * (zIpol / zIp)
         !  ** Integrate "VLOOPpol*Ipol/Ip" w.r.t. time ; 3rd term on RHS of eq. (28)
         CPsi(3) = CPsi_old(3) &
              &   + 0.5d0 * (PoyntS(2) / zIp + PS2_old / zIp_old) &
              &   * (T_TX - T_TX_old)

         CPsi_old(:) = CPsi(:)
      else if( irestart == 1 ) then
         ! Right after restart
         !     CPsi_old(:) = CPsi(:)
         CPsi(:)     = 0.d0
         CPsi_old(:) = 0.d0
      else
         ! At initial
         PoyntI(:)   = 0.d0
         CPsi(:)     = 0.d0
         CPsi_old(:) = 0.d0
      end if
      T_TX_old     = T_TX
      zIp_old      = zIp
      VLOOP_old    = VLOOP
      zIpol_old    = zIpol
      VLOOPpol_old = VLOOPpol
      PS2_old      = PoyntS(2)
      WFEengy_old  = WFEengy
      WFMengy_old  = WFMengy
      WJengy_old   = WJengy
      WFMpengy_old = WFMpengy
      WFMtengy_old = WFMtengy

      !   --- Ejima coefficient ---

      CEjima = CPsi(2) / (rMU0 * RR * zIp)

      !  write(6,*) intg_vol(dfdx(vv,fipol,nrmax,0)*psitdotv(nrmax))/(fipol(nrmax)*psitdotv(nrmax))
      !  write(6,'(6ES15.7)') PoyntS(2),PoyntI(4),PoyntS(1)-PoyntI(1)-PoyntI(3)-sum(PoyntR(1:2)),&
      !       &    sum(PoyntS(1:2))-sum(PoyntI(1:2))-sum(PoyntR(1:2)),PoyntS(2),PoyntI(4)+PoyntR(2)

      deallocate(EngyJP,EngyJT,EngyEp,EngyEt,EngyMp,EngyMt)
      !  write(6,'(6ES15.7)') t_tx,CPsi(1),CPsi(2),WFengy,WJengy,CEjima
      !  write(6,'(6ES15.7)') t_tx,VPoynt(1)+VPoynt(2),Vloop
      !  write(6,'(6ES15.7)') t_tx,sum(PoyntS),sum(PoyntI)+sum(PoyntR),PoyntI(1),PoyntI(2),sum(PoyntR)
    end block

!!$!  Beam-driven particle flux estimated by the model [cf. Honda NF 2012]
!!$  do nr = 0, nrmax
!!$     write(6,'(F8.5,3ES15.7)') rho(nr) &
!!$          & , xmu(NR,1,1,1)*(PNbV(NR)*BUbparV(NR)/(Var(NR,1)%n*bbrt(NR)))/aee*rr/sdt(NR)*1.d-20 &
!!$          & , fipol(NR)/(aee*bbt(NR))*(achgb**2*PNbV(NR)*BUbparV(NR)/(Zeff(NR)*Var(NR,1)%n)) &
!!$          &  *xmu(NR,1,1,1)*(1.d0-linliu_L31(ft(NR),Zeff(NR)))/sdt(NR)*1.d-20 &
!!$          & , fipol(NR)/(aee*bbt(NR))*(achgb**2*PNbV(NR)*BUbparV(NR)/(Zeff(NR)*Var(NR,1)%n)) &
!!$          &  *xmu(NR,1,1,1)*(1.d0-cL31(NR,1))/sdt(NR)*1.d-20
!!$  end do

  end subroutine TXGLOB

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

end module tx_glob

!!$real(8) function linliu_L31(xgt,zeff)
!!$  implicit none
!!$  real(8) :: xgt, zeff
!!$  
!!$  linliu_L31 = xgt * ( ( 0.754d0 + 2.21d0  * zeff + zeff**2 ) &
!!$       &     + xgt *   ( 0.348d0 + 1.243d0 * zeff + zeff**2 )) &
!!$       &  / (                       1.414d0 * zeff +        zeff**2 &
!!$       &     + xgt    * ( 0.754d0 + 2.657d0 * zeff + 2.d0 * zeff**2) &
!!$       &     + xgt**2 * ( 0.348d0 + 1.243d0 * zeff +        zeff**2 ))
!!$
!!$end function linliu_L31
