!***************************************************************
!
!   Write Data
!
!***************************************************************

subroutine TXWDAT
  use tx_commons, only : PI, RB, PN01V, NRMAX, WPT, Var
  use tx_interface, only : rLINEAVE
  use tx_core_module, only : intg_area

  implicit none
  real(8) :: rNbar

  !     ***** Volume-averaged density *****

  rNbar = intg_area(Var(:,1)%n) * 1.D20
  rNbar = rNbar / (PI * RB**2)

  write(6,'((1X,A," =",ES9.2,3(2X,A,"=",ES9.2)))') &
       &     'Ne(0)',    Var(0,1)%n,  &
       &     'UePhi(0)', Var(0,1)%Uph / 1.D3,  &
       &     'UiPhi(0)', Var(0,2)%Uph / 1.D3,  &
       &     'N0(RB)',   PN01V(NRMAX) * 1.D20, &
       &     'NB(0)',    rLINEAVE(0.D0)   / 1.D20,  &
       &     'NB(0.24)', rLINEAVE(0.24D0) / 1.D20,  &
       &     'NB(0.60)', rLINEAVE(0.6D0)  / 1.D20,  &
       &     'PF    ',   Var(0,1)%n * 1.D20 / rNbar
  write(6,'(1X,A," =",ES9.2,2(2X,A,"=",ES9.2))') &
       &     'Te(0)',    Var(0,1)%T,  &
       &     'Ti(0)   ', Var(0,2)%T,  &
       &     'Wst     ', WPT

end subroutine TXWDAT

!***************************************************************
!
!   Write Data 2
!
!***************************************************************

subroutine TXWDAT2

  use tx_graphic, only : GTY, NGT
  implicit none
  real(8) :: PNeMIN, PNeMAX, NB0MIN, NB0MAX, UiphMIN, UiphMAX

  PNeMIN  = minval(real(GTY(0:NGT+1, 1),8))
  PNeMAX  = maxval(real(GTY(0:NGT+1, 1),8))
  NB0MIN  = minval(real(GTY(0:NGT+1, 2),8))
  NB0MAX  = maxval(real(GTY(0:NGT+1, 2),8))
  UiphMIN = minval(real(GTY(0:NGT+1,11),8))
  UiphMAX = maxval(real(GTY(0:NGT+1,11),8))

  write(6,'((1X,A,"=",ES9.2,2(2X,A,"=",ES9.2)))') &
       &     'MAX(Ne(0))',     PNeMAX / 1.E20,  &
       &     'MAX(NB(0))',     NB0MAX / 1.E20,  &
       &     'MAX(UiPhi(0))', UiphMAX / 1.E3,   &
       &     'MIN(Ne(0))',     PNeMIN / 1.E20,  &
       &     'MIN(NB(0))',     NB0MIN / 1.E20,  &
       &     'MIN(UiPhi(0))', UiphMIN / 1.E3

end subroutine TXWDAT2

!***************************************************************
!
!   Write Statistic Data
!
!***************************************************************

subroutine TXSTAT
  use tx_commons, only : VOLAVN, ALI, VLOOP, TAUE1, TAUE2, TAUEP, TAUEH, BETAA, &
       &                 BETAPA, BETAN, Q, ANSAV, rIp, PI, RA, NRA, NRMAX, &
       &                 rMus, Chis, TSAV, Gamma_a, TAUPA, achg, &
       &                 rKeV, amas, amp, RR, rNuei, rho, Var, vlt
  use tx_glob, only : cal_flux
  use libitp, only: aitken2p
  implicit none
  integer(4) :: NR, NRL1, NRL2
  real(8) :: rhol1, rhol2, rmuil, chiil, uiphl, PTeVL, WDe, rNueiL
  character(len=21) :: cfmt = '(1X,2(A27,ES10.3,3X))'

  rhol1 = 0.3d0 ; rhol2 = 0.5d0
  do NR = 0, NRMAX-1
     if(rho(NR) <= rhol1 .and. rho(NR+1) >= rhol1) then
        NRL1 = NR
     else if(rho(NR) <= rhol2 .and. rho(NR+1) >= rhol2) then
        NRL2 = NR
        exit
     end if
  end do

  call cal_flux

  rmuil = aitken2p(rhol1,rmus(NRL1,2),rmus(NRL1+1,2),rmus(NRL1+2,2),rho(NRL1),rho(NRL1+1),rho(NRL1+2))
  chiil = aitken2p(rhol1,chis(NRL1,2),chis(NRL1+1,2),chis(NRL1+2,2),rho(NRL1),rho(NRL1+1),rho(NRL1+2))
  uiphl = aitken2p(rhol1,Var(NRL1,2)%Uph,Var(NRL1+1,2)%Uph,Var(NRL1+2,2)%Uph,rho(NRL1),rho(NRL1+1),rho(NRL1+2))

  ! For effective collision frequency at rho = 0.5

  PTeVL = aitken2p(rhol2,Var(NRL2,1)%T,Var(NRL2+1,1)%T,Var(NRL2+2,1)%T,rho(NRL2),rho(NRL2+1),rho(NRL2+2))
  ! WDe : curvature drift frequency
  WDe   = 2.d0 * sqrt(0.1d0) * sqrt(PTeVL * rKeV / (amas(2)*amp)) / RR
  rNueiL = aitken2p(rhol2,rNuei(NRL2),rNuei(NRL2+1),rNuei(NRL2+2),rho(NRL2),rho(NRL2+1),rho(NRL2+2))

  write(6,cfmt) "Vol. ave. of neutrality  = ", VOLAVN
  write(6,cfmt) "Inductance               = ", ALI, &
       &        "Loop voltage             = ", VLOOP
  write(6,cfmt) "Confinement time 1       = ", TAUE1, &
       &        "Confinement time 2       = ", TAUE2
  write(6,cfmt) "L-mode scaling time      = ", TAUEP, &
       &        "IPB98(y,2) scaling time  = ", TAUEH
  write(6,cfmt) "Beta                     = ", BETAA, &
       &        "Poloidal beta            = ", BETAPA
  write(6,cfmt) "Normalized beta          = ", BETAN, &
       &        "tau_p inside sep.        = ", TAUPA
  write(6,cfmt) "Ion flux thru sep.       = ", achg(2)*Var(NRA,2)%n*Var(NRA,2)%Ur*1.D20, &
       &        "Ion flux thru sep. (est) = ", Gamma_a
  write(6,cfmt) "Vol. averaged e density  = ", ANSAV(1), &
       &        "Greenwald density        = ", rIp / (PI * RA**2)
  write(6,cfmt) "Vol. averaged e temp.    = ", TSAV(1), &
       &        "Vol. averaged i temp.    = ", TSAV(2)
  write(6,cfmt) "Safety factor on axis    = ", Q(0), &
       &        "Safety factor at sep.    = ", Q(NRA)
  write(6,cfmt) "Ion Prandtl num. at 0.3  = ", rmuil/chiil, &
       &        "Ion tor. velocity at 0.3 = ", uiphl
  write(6,cfmt) "Eff .col. freq. at 0.5   = ", rNueiL/WDe, &
       &        "Plasma volume inside sep.= ", vlt(NRA)

end subroutine TXSTAT

!***************************************************************
!
!   Steady state check
!
!***************************************************************

subroutine steady_check
  use tx_commons, only : NRMAX, Var, T_TX, rho
  use libitp, only: aitken2p

  implicit none
  integer(4) :: nr
  integer(4), save :: nrl = 0
  real(8) :: rhol, uiphl, pnevl, pevl
  real(8), save :: uiphl_old = 0.d0, pnevl_old = 0.d0, pevl_old = 0.d0

  ! Seel a grid number "nrl" nearest rho=0.3
  if(nrl == 0) then
     rhol = 0.3d0
     do NR = 0, NRMAX-1
        if(rho(NR) <= rhol .and. rho(NR+1) >= rhol) then
           NRL = NR
           exit
        end if
     end do
  end if
  
  ! Ion toroidal velocity
  uiphl = aitken2p(rhol,Var(nrl,2)%Uph,Var(nrl+1,2)%Uph,Var(nrl+2,2)%Uph,rho(nrl),rho(nrl+1),rho(nrl+2))

  ! Electron density
  PNeVl = aitken2p(rhol,Var(nrl,1)%n,  Var(nrl+1,1)%n,  Var(nrl+2,1)%n,  rho(nrl),rho(nrl+1),rho(nrl+2))

  ! Electron pressure
  PeVl  = aitken2p(rhol,Var(nrl,1)%p,  Var(nrl+1,1)%p,  Var(nrl+2,1)%p,  rho(nrl),rho(nrl+1),rho(nrl+2))

  ! Dispaly
  if(uiphl /= 0.d0) then
     write(6,*) real(t_tx),real(abs(uiphl - uiphl_old)/uiphl), &
          &                real(abs(pnevl - pnevl_old)/pnevl), &
          &                real(abs(pevl  - pevl_old )/pevl ), real(uiphl)
  end if
  uiphl_old = uiphl
  pnevl_old = pnevl
  pevl_old  = pevl

end subroutine steady_check

!***********************************************************
!
!  LINE AVERAGE OF rN
!
!***********************************************************

REAL(8) FUNCTION rLINEAVE(Rho)

  use tx_commons, only : RA, NRMAX, Var
  implicit none
  real(8), intent(IN) :: Rho
  integer(4) :: I, IR, NY = 100
  real(8) :: D, DY, Y, RL, SUML

  SUML = 0.D0
  D = Rho * RA
  DY = sqrt(RA*RA - D*D) / NY
  do I = 0, NY
     Y = DY * I
     RL = sqrt(Y*Y + D*D)
     IR = nint(RL * NRMAX / RA)
     SUML = SUML + Var(IR,1)%n * 1.D20 * DY
  end do
  rLINEAVE = SUML / sqrt(RA**2 - D**2)

end FUNCTION rLINEAVE

!***************************************************************
!
!   Save transport data
!
!***************************************************************

subroutine TXSAVE
  use tx_commons, only : &
       & SLID,RA,rhob,rhoaccum,RR,BB,rbvt,ravl,rbvl,amas,achg,amb,achgb,rIPs,rIPe, &
       & PN0,PNa,PTe0,PTea,PTi0,PTia,PTz0,PTza, &
       & PROFJ,PROFN1,PROFN2,PROFT1,PROFT2,Uiph0,PROFD,PROFD1,PROFD2,PROFDB,PROFM,PROFM1,PROFMB,PROFC,PROFC1,PROFCB, &
       & Dfs0,VWpch0,rMus0,WPM0,Chis0,ChiNC,FSDFIX,FSANOM,FSCBKP,FSCBSH,rG1, &
       & FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV,FSCX,FSLC,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSD01,FSD02,FSD03, &
       & FSRP,FSNF,FSADV,FSADVB,FSUG,FSHL,MDLC,rLn,rLT,Ebmax,esps,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,&
       & PNBHT1,PNBHT2,rNRFe,RRFew,RRFe0,PRFHe,rNRFi,RRFiw,RRFi0,PRFHi,PNBCD, &
       & PN0s,V0,rGamm0,rGASPF,PNsDIV,PTsDIV, &
       & DT,EPS,ADV,tiny_cap,CMESH0,WMESH0,CMESH,WMESH, &
       & ICMAX,NRMAX,NTMAX,NTSTEP,T_TX,TMAX,NT,NTCUM,NQMAX,IERR,X, &
       & NLCMAX,NCM,NTCOIL,DltRPn,m_pol,n_tor, &
       & MODEAV,IDIAG,MDLPCK,MODECV,oldmix,iSUPG3,iSUPG6,iSUPG8,SUPGstab, &
       & IGBDF,MDFIXT,MDBEAM,MDOSQZ,MDOSQZN,MDLETA,MDLNEO,MDANOM, &
       & MDLNBD,PNBMPD,PNBPTC,thrp,kappa,CPsi,VPoynt,PoyntI
  use tx_graphic, only : NGYTM,NGYVM,MODEG,MODEGL,NGT,NGVV,NGRSTP,NGTSTP,NGVSTP,GTY,GVY,GQY,GTX,GVX
  use tx_interface, only : TOUPPER

  implicit none
  integer(4) :: IST, NQ, NR, NC, I, IGYT, IGYV
  character(len=100) :: TXFNAM, RCSId
  character(len=1) :: STR
  logical :: LEX

  RCSId = ' '

  do 
     write(6,*) '# INPUT : SAVE FILE NAME'
     flush(6)
     read(*,'(A100)',iostat=IST) TXFNAM
     if (IST > 0) then
        cycle
     else if (IST < 0) then
        return
     end if
     inquire(file=TXFNAM,exist=LEX)
     if (LEX) then
        write(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',  &
             &              'ARE YOU SURE {Y/N} ?'
        flush(6)
        read(*,'(A1)') STR
        call TOUPPER(STR)
        if (STR == 'Y') then
           open(21,file=TXFNAM,iostat=IST,status='old',form='unformatted')
           if (IST == 0) then
              write(6,*) '# OLD FILE ( ', TXFNAM(1:len_trim(TXFNAM)), &
                   &     ' ) IS ASSIGNED FOR OUTPUT.'
              exit
           else if (IST > 0) then
              write(6,*) 'XX  OLD FILE open ERROR !, iostat = ', IST
           end if
        end if
     else
        open(21,file=TXFNAM,iostat=IST,status='new',form='unformatted')
        if (IST == 0) then
           write(6,*) '# NEW FILE ( ', TXFNAM(1:len_trim(TXFNAM)), &
                &     ' ) IS CREATED FOR OUTPUT.'
           exit
        else if (IST > 0) then
           write(6,*) 'XX  NEW FILE open ERROR !, iostat = ', IST
        end if
     end if
  end do

  ! *** Variables defined in tx_commons but not included in the following ***
  !
  !   VWpch0, Tqt0, Tqp0, NEMAX, NRA, NRC, DelRho, DelN,
  !   EpsH, Q0, QA, NCph, NCth, DMAG0, RMAGMN, RMAGMX,
  !   MDITSN, MDITST, MDINTN, MDINTT, MDINTC
  !
  ! *************************************************************************

  write(21) SLID
  write(21) RCSId

  write(21) RA,rhob,rhoaccum,RR,BB,rbvt,ravl,rbvl
  write(21) amas,achg,amb,achgb,rIPs,rIPe
  write(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PTz0,PTza,PROFJ,PROFN1,PROFN2,PROFT1,PROFT2,Uiph0
  write(21) PROFD,PROFD1,PROFD2,PROFDB,PROFM,PROFM1,PROFMB,PROFC,PROFC1,PROFCB
  write(21) Dfs0,VWpch0,rMus0,WPM0,Chis0,ChiNC
  write(21) FSDFIX,FSANOM,FSCBKP,FSCBSH,rG1,FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV
  write(21) FSCX,FSLC,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSD01,FSD02,FSD03,FSRP,FSNF,FSADV,FSADVB,FSUG,FSHL
  write(21) rLn,rLT
  write(21) Ebmax,esps,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,PNBHT1,PNBHT2
  write(21) rNRFe,RRFew,RRFe0,PRFHe,rNRFi,RRFiw,RRFi0,PRFHi,PNBCD,PNBMPD,PNBPTC
  write(21) PN0s,V0,rGamm0,rGASPF,PNsDIV,PTsDIV
  write(21) DltRPn,kappa
  write(21) DT,EPS,ADV,tiny_cap,CMESH0,WMESH0,CMESH,WMESH
  write(21) ICMAX,NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
  write(21) MODEG,MODEAV,MODEGL,IDIAG,MDLPCK,MODECV,oldmix,iSUPG3,iSUPG6,iSUPG8,SUPGstab,IGBDF
  write(21) MDFIXT,MDBEAM,MDOSQZ,MDOSQZN,MDLETA,MDLNEO,MDANOM
  write(21) MDLNBD,MDLC,NTCOIL,m_pol,n_tor

  write(21) T_TX,TMAX,NT,NTCUM,NQMAX,IERR
  write(21) ((X(NR,NQ), NR=0, NRMAX), NQ=1, NQMAX)

  write(21) NGT
  write(21) (GTX(I), I=0, NGT)
  write(21) (GVX(I), I=0, NGVV)
  write(21) ((GTY(I,IGYT), I=0, NGT),  IGYT =1, NGYTM)
  write(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
  write(21) (NLCMAX(NQ), NQ=1,NQMAX)
  write(21) (((GQY(NR,NC,NQ), NR=0, NRMAX), NC=1, NCM), NQ=1, NQMAX)
  write(21) (PoyntI(I), I=1, 2)
  write(21) (CPsi(I), I=0, 3)
  write(21) (VPoynt(I), I=0, 3)
  write(21) (thrp(I), I=1, 2*NRMAX)
  close(21)
  write(6,*) '# DATA WAS SUCCESSFULLY SAVED IN THE FILE.'

end SUBROUTINE TXSAVE

!***************************************************************
!
!   Load transport data
!
!***************************************************************

subroutine TXLOAD(IST)
  use tx_commons, only : &
       & allocate_txcomm, deallocate_txcomm, &
       & RA,RB,rhob,rhoaccum,RR,BB,rbvt,ravl,rbvl,amas,achg,amb,achgb,rIPs,rIPe, &
       & PN0,PNa,PTe0,PTea,PTi0,PTia,PTz0,PTza, &
       & PROFJ,PROFN1,PROFN2,PROFT1,PROFT2,Uiph0,PROFD,PROFD1,PROFD2,PROFDB,PROFM,PROFM1,PROFMB,PROFC,PROFC1,PROFCB, &
       & Dfs0,VWpch0,rMus0,WPM0,Chis0,ChiNC,FSDFIX,FSANOM,FSCBKP,FSCBSH,rG1, &
       & FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV,FSCX,FSLC,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSD01,FSD02,FSD03, &
       & FSRP,FSNF,FSADV,FSADVB,FSUG,FSHL,MDLC,rLn,rLT,Ebmax,esps,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP, &
       & PNBHT1,PNBHT2,rNRFe,RRFew,RRFe0,PRFHe,rNRFi,RRFiw,RRFi0,PRFHi,PNBCD, &
       & PN0s,V0,rGamm0,rGASPF,PNsDIV,PTsDIV, &
       & DT,EPS,ADV,tiny_cap,CMESH0,WMESH0,CMESH,WMESH, &
       & ICMAX,NRMAX,NTMAX,NTSTEP,T_TX,TMAX,NT,NTCUM,NQMAX,IERR,X, &
       & NLCMAX,NCM,NTCOIL,DltRPn,m_pol,n_tor, &
       & MODEAV,IDIAG,MDLPCK,MODECV,oldmix,iSUPG3,iSUPG6,iSUPG8,SUPGstab, &
       & IGBDF,MDFIXT,MDBEAM,MDOSQZ,MDOSQZN,MDLETA,MDLNEO,MDANOM, &
       & MDLNBD,PNBMPD,PNBPTC,rIP,thrp,kappa,CPsi,VPoynt,PoyntI, &
       & ErV,PTsV_FIX,PNsV_FIX,ErV_FIX, &
       & rMU0,rMUb1,rMUb2,NEMAX,ICONT,TAUE2,LQb1,Var,irestart, NSM
  use tx_graphic, only : allocate_txgraf,deallocate_txgraf, &
       &                 MODEG,MODEGL,NGYTM,NGYVM,NGR,NGT,NGVV,NGRSTP,NGTSTP,NGVSTP, &
       &                 GTY,GVY,GQY,GTX,GVX
  use tx_variables
  use tx_coefficients, only : TXCALA
  use tx_parameter_control, only : TXPARM_CHECK
  use tx_glob, only : TXGLOB
  use mod_eqneo, only : wrap_eqneo
  use mod_cross_section, only : spline_table_carbon_rate_coef_adas, spline_table_beam_rate_coef
  implicit none
  integer(4), intent(out) :: IST
  integer(4) :: NQ, NR, NC, I, IGYT, IGYV
  character(len=100) ::  TXFNAM, RCSId
  character(len=8) :: LOADSLID
  logical :: LEX

  do 
     write(6,*) '# INPUT : LOAD FILE NAME'
     flush(6)
     read(*,'(A100)',iostat=IST) TXFNAM
     if (IST > 0) then
        cycle
     else if (IST < 0) then
        return
     end if
     inquire(file=TXFNAM,exist=LEX)
     if (LEX) then
        open(21,file=TXFNAM,iostat=IST,status='old',form='unformatted')
        if (IST == 0) then
           write(6,*) '# OLD FILE ( ', TXFNAM(1:len_trim(TXFNAM)),  &
                &     ' ) IS ASSIGNED FOR INPUT.'
           exit
        else if (IST > 0) then
           write(6,*) 'XX  OLD FILE open ERROR !, iostat = ', IST
        end if
     else
        write(6,*) 'XX  FILE ( ', TXFNAM(1:len_trim(TXFNAM)), ' ) DOES NOT exist !'
     end if
  end do

  read(21,iostat=IST) LOADSLID
  if (IST > 0) then
     write(6,*) 'XX read ERROR in TXLOAD !'
     close(21)
     return
  end if
  read(21) RCSId

  read(21) RA,rhob,rhoaccum,RR,BB,rbvt,ravl,rbvl
  read(21) amas,achg,amb,achgb,rIPs,rIPe
  read(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PTz0,PTza,PROFJ,PROFN1,PROFN2,PROFT1,PROFT2,Uiph0
  read(21) PROFD,PROFD1,PROFD2,PROFDB,PROFM,PROFM1,PROFMB,PROFC,PROFC1,PROFCB
  read(21) Dfs0,VWpch0,rMus0,WPM0,Chis0,ChiNC
  read(21) FSDFIX,FSANOM,FSCBKP,FSCBSH,rG1,FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV
  read(21) FSCX,FSLC,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSD01,FSD02,FSD03,FSRP,FSNF,FSADV,FSADVB,FSUG,FSHL
  read(21) rLn,rLT
  read(21) Ebmax,esps,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,PNBHT1,PNBHT2
  read(21) rNRFe,RRFew,RRFe0,PRFHe,rNRFi,RRFiw,RRFi0,PRFHi,PNBCD,PNBMPD,PNBPTC
  read(21) PN0s,V0,rGamm0,rGASPF,PNsDIV,PTsDIV
  read(21) DltRPn,kappa
  read(21) DT,EPS,ADV,tiny_cap,CMESH0,WMESH0,CMESH,WMESH
  read(21) ICMAX,NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
  read(21) MODEG,MODEAV,MODEGL,IDIAG,MDLPCK,MODECV,oldmix,iSUPG3,iSUPG6,iSUPG8,SUPGstab,IGBDF
  read(21) MDFIXT,MDBEAM,MDOSQZ,MDOSQZN,MDLETA,MDLNEO,MDANOM
  read(21) MDLNBD,MDLC,NTCOIL,m_pol,n_tor

  read(21) T_TX,TMAX,NT,NTCUM,NQMAX,IERR

  call allocate_txcomm(ierr)
  call allocate_txgraf(ierr)
  if(ierr /= 0) then
     call deallocate_txcomm
     call deallocate_txgraf
     write(6,*) "XX Allocation error : TXLOAD"
  end if

  read(21) ((X(NR,NQ), NR=0, NRMAX), NQ=1, NQMAX)

  read(21) NGT
  read(21) (GTX(I), I=0, NGT)
  read(21) (GVX(I), I=0, NGVV)
  read(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYTM)
  read(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
  read(21) (NLCMAX(NQ), NQ=1,NQMAX)
  read(21) (((GQY(NR,NC,NQ), NR=0, NRMAX), NC=1, NCM), NQ=1, NQMAX)
  read(21) (PoyntI(I), I=1, 2)
  read(21) (CPsi(I), I=0, 3)
  read(21) (VPoynt(I), I=0, 3)
  read(21) (thrp(I), I=1, 2*NRMAX)
  close(21)
  write(6,'(2A)') '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE. ID = ',LOADSLID

  rb = rhob * ra

  rIP=rIPs
  ICONT = 1
  irestart = 1

  call TXPARM_CHECK

  !  Read spline table for rate coefficients
  call spline_table_carbon_rate_coef_adas
  call spline_table_beam_rate_coef

  NEMAX = NRMAX
  ! Create the mesh and read an equilibrium
  call TXCALM

!!  if(rMUb1 == rMU0 .and. (PNBHT1 /= 0.D0 .OR. PNBHT2 /= 0.D0 .OR. PNBHP /= 0.D0)) then
  if(rMUb1 == rMU0 .and. (maxval(X(:,LQb1)) > epsilon(1.d0))) then
     rMUb1 = 1.D0
     rMUb2 = rMU0
  end if

  call TXCALV(X,0)

  do i = 1, NSM
     PNsV_FIX(:,i) = Var(:,i)%n
     PTsV_FIX(:,i) = Var(:,i)%T
  end do
  ErV_FIX (:) = ErV (:)

  call wrap_eqneo

  call TXCALC(0)
  call TXCALA

  ! TAUE2 uses data one step before the data was stored.
  ! Then TAUE2 is reconstituted by using the graphic data of TAUE2.
  TAUE2 = real(GTY(NGT,34),8)

  ! Reset start point of graphics
  NGT=-1
  NGR=-1
  NGVV=-1

  call TXGLOB
  call TXWDAT
  call TXWDAT2

  irestart = 0

end subroutine TXLOAD

!***************************************************************
!
!   Save graphic data
!
!***************************************************************

subroutine TXGSAV

  use tx_commons, only : &
       & SLID,RA,rhob,RR,BB,rbvt,ravl,rbvl,amas,achg,amb,achgb,PTe0,PTea,PTi0,PTia,PTz0,PTza, &
       & Dfs0,rMus0,WPM0,Chis0,FSDFIX,FSANOM,FSCBKP,FSCBSH, &
       & FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV,PROFD,PROFC, &
       & FSCX,FSLC,FSRP,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION, &
       & FSD01,FSD02,FSD03,Ebmax,esps,PNBH,PNBHP,PNBHT1,PNBHT2,PRFHe,PRFHi,PNBCD,PNBMPD,PNBPTC, &
       & V0,rGamm0,rGASPF,PNsDIV,PTsDIV, &
       & DT,NRMAX,NTMAX,NTSTEP,rG1,T_TX,TMAX,NT,NTCUM,NQMAX,IERR,X, &
       & NLCMAX,NCM,DltRPn,thrp,kappa,CPsi,VPoynt,PoyntI
  use tx_graphic, only : MODEG,MODEGL,NGYRM,NGYTM,NGYVM,NGR,NGT,NGVV,NGRSTP,NGTSTP,NGVSTP, &
       &                 GTY,GVY,GQY,GY,GYT,GTX,GVX,GT
  implicit none
  integer(4) :: IST, NQ, NR, NC, IGR, I, IGYR, IGYT, IGYV
  character(len=100) :: TXFNAM, RCSId
  character(len=1) :: STR
  logical :: LEX

  RCSId = ' '

  do
     write(6,*) '# INPUT : SAVE FILE NAME'
     flush(6)
     read(*,'(A100)',iostat=IST) TXFNAM
     if (IST > 0) then
        cycle
     else if(IST < 0) then
        return
     end if
     inquire(file=TXFNAM,exist=LEX)
     if (LEX) then
        write(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',  &
             &              'ARE YOU SURE {Y/N} ?'
        flush(6)
        read(*,'(A1)') STR
        if (STR /= 'Y' .and. STR /= 'y') then
        else
           open(21,file=TXFNAM,iostat=IST,status='old',form='unformatted')
           if (IST == 0) then
              write(6,*) '# OLD FILE ( ', TXFNAM(1:len_trim(TXFNAM)), &
                   &     ' ) IS ASSIGNED FOR OUTPUT.'
              exit
           else if (IST > 0) then
              write(6,*) 'XX  OLD FILE open ERROR !, iostat = ', IST
           end if
        end if
     else
        open(21,file=TXFNAM,iostat=IST,status='new',form='unformatted')
        if (IST == 0) then
           write(6,*) '# NEW FILE ( ', TXFNAM(1:len_trim(TXFNAM)), &
                &     ' ) IS CREATED FOR OUTPUT.'
           exit
        else if (IST > 0) then
           write(6,*) 'XX  NEW FILE open ERROR !, iostat = ', IST
        end if
     end if
  end do

    write(21) SLID
!!$    write(21) RCSId

  write(21) RA,rhob,RR,BB,rbvt,ravl,rbvl
  write(21) amas,achg,amb,achgb
  write(21) PTe0,PTea,PTi0,PTia,PTz0,PTza
  write(21) Dfs0,rMus0,WPM0,Chis0
  write(21) FSDFIX,FSANOM,FSCBKP,FSCBSH,rG1,FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV
  write(21) PROFD,PROFC
  write(21) FSCX,FSLC,FSRP,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSD01,FSD02,FSD03
  write(21) Ebmax,esps,PNBH,PNBHP,PNBHT1,PNBHT2
  write(21) PRFHe,PRFHi,PNBCD,PNBMPD,PNBPTC
  write(21) V0,rGamm0,rGASPF,PNsDIV,PTsDIV
  write(21) DltRPn,kappa
  write(21) DT
  write(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
  write(21) MODEG,MODEGL

  write(21) T_TX,TMAX,NT,NTCUM,NQMAX,IERR
  write(21) ((X(NR,NQ), NR=0, NRMAX), NQ=1, NQMAX)

  write(21) NGR,NGT,NGVV
  write(21) (GT(IGR), IGR=0, NGR)
  write(21)(((GY%v(I,IGR,IGYR), I=0,NRMAX), IGR=0,NGR), IGYR=1,NGYRM)
  write(21) (GTX(I), I=0, NGT)
  write(21) (GVX(I), I=0, NGVV)
  write(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYTM)
  write(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
  write(21) (NLCMAX(NQ), NQ=1,NQMAX)
  write(21) (((GQY(NR,NC,NQ), NR=0, NRMAX), NC=1, NCM), NQ=1, NQMAX)
  write(21) (((GYT%v(NR,I,IGYR), NR=0,NRMAX), I=0,NGT), IGYR=1,NGYRM)
  write(21) (PoyntI(I), I=1, 2)
  write(21) (CPsi(I), I=0, 3)
  write(21) (VPoynt(I), I=0, 3)
  write(21) (thrp(I), I=1, 2*NRMAX)
  close(21)
  write(6,*) '# DATA WAS SUCCESSFULLY SAVED IN THE FILE.'

end subroutine TXGSAV

!***************************************************************
!
!   Load graphic data
!
!***************************************************************

subroutine TXGLOD(IST)

  use tx_commons, only : &
       & allocate_txcomm, deallocate_txcomm, &
       & RA,RB,rhob,RR,BB,rbvt,ravl,rbvl,amas,achg,amb,achgb,PTe0,PTea,PTi0,PTia,PTz0,PTza, &
       & Dfs0,rMus0,WPM0,Chis0,FSDFIX,FSANOM,FSCBKP,FSCBSH, &
       & FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV,PROFD,PROFC, &
       & FSCX,FSLC,FSRP,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION, &
       & FSD01,FSD02,FSD03,Ebmax,esps,PNBH,PNBHP,PNBHT1,PNBHT2,PRFHe,PRFHi,PNBCD,PNBMPD,PNBPTC, &
       & V0,rGamm0,rGASPF,PNsDIV,PTsDIV, &
       & DT,NRMAX,NTMAX,NTSTEP,rG1,T_TX,TMAX,NT,NTCUM,NQMAX,IERR,X, &
       & NLCMAX,NCM,DltRPn,thrp,kappa,CPsi,VPoynt,PoyntI!,rho
  use tx_graphic, only : allocate_txgraf,deallocate_txgraf, &
       &                 MODEG,MODEGL,NGYRM,NGYTM,NGYVM,NGR,NGT,NGVV,NGRSTP,NGTSTP,NGVSTP, &
       &                 GTY,GVY,GQY,GY,GYT,GTX,GVX,GT
  use tx_variables

  implicit none
  integer(4), intent(out) :: IST
  integer(4) :: NQ, NR, NC, IGR, I, IGYR, IGYT, IGYV
  character(len=100) :: TXFNAM
!!  character(len=100) :: RCSId
  character(len=8) :: LOADSLID
  logical :: LEX

  do 
     write(6,*) '# INPUT : LOAD FILE NAME'
     flush(6)
     read(*,'(A100)',iostat=IST) TXFNAM
     if (IST > 0) then
        cycle
     else if(IST < 0) then
        return
     end if
     inquire(file=TXFNAM,exist=LEX)
     if (LEX) then
        open(21,file=TXFNAM,iostat=IST,status='old',form='unformatted')
        if (IST == 0) then
           write(6,*) '# OLD FILE ( ', TXFNAM(1:len_trim(TXFNAM)),&
                &     ' ) IS ASSIGNED FOR INPUT.'
           exit
        else if (IST > 0) then
           write(6,*) 'XX  OLD FILE open ERROR !, iostat = ', IST
        end if
     else
        write(6,*) 'XX  FILE ( ', TXFNAM(1:len_trim(TXFNAM)), ' ) DOES NOT exist !'
     end if
  end do

  read(21,iostat=IST) LOADSLID
  if (IST > 0) then
     write(6,*) 'XX read ERROR in TXGLOD !'
     close(21)
     return
  end if
!!$    !  if(LOADSLID(1:5) == 'tx459') then
!!$    read(21) RCSId

  read(21) RA,rhob,RR,BB,rbvt,ravl,rbvl
  read(21) amas,achg,amb,achgb
  read(21) PTe0,PTea,PTi0,PTia,PTz0,PTza
  read(21) Dfs0,rMus0,WPM0,Chis0
  read(21) FSDFIX,FSANOM,FSCBKP,FSCBSH,rG1,FSBOHM,FSPCL,FSVAHL,FSMPCH,FSPARV
  read(21) PROFD,PROFC
  read(21) FSCX,FSLC,FSRP,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,FSD01,FSD02,FSD03
  read(21) Ebmax,esps,PNBH,PNBHP,PNBHT1,PNBHT2
  read(21) PRFHe,PRFHi,PNBCD,PNBMPD,PNBPTC
  read(21) V0,rGamm0,rGASPF,PNsDIV,PTsDIV
  read(21) DltRPn,kappa
  read(21) DT
  read(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
  read(21) MODEG,MODEGL

  read(21) T_TX,TMAX,NT,NTCUM,NQMAX,IERR

  call allocate_txcomm(ierr)
  call allocate_txgraf(ierr)
  if(ierr /= 0) then
     call deallocate_txcomm
     call deallocate_txgraf
     write(6,*) "XX Allocation error : TXGLOD"
  end if

  read(21) ((X(NR,NQ), NR=0, NRMAX), NQ=1, NQMAX)

  read(21) NGR,NGT,NGVV
  read(21) (GT(IGR), IGR=0, NGR)
  read(21) (((GY%v(I,IGR,IGYR), I=0, NRMAX), IGR=0, NGR), IGYR=1, NGYRM)
  read(21) (GTX(I), I=0, NGT)
  read(21) (GVX(I), I=0, NGVV)
  read(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYTM)
  read(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
  read(21) (NLCMAX(NQ), NQ=1,NQMAX)
  read(21) (((GQY(NR,NC,NQ), NR=0, NRMAX), NC=1, NCM), NQ=1, NQMAX)
  read(21) (((GYT%v(NR,I,IGYR), NR=0,NRMAX), I=0,NGT), IGYR=1,NGYRM)
  read(21) (PoyntI(I), I=1, 2)
  read(21) (CPsi(I), I=0, 3)
  read(21) (VPoynt(I), I=0, 3)
  read(21) (thrp(I), I=1, 2*NRMAX)
  !  end if
  close(21)
  write(6,'(2A)') '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE. ID = ',LOADSLID

  rb = rhob * ra

  call TXCALM
!!$    call TXCALV(X)
!!$    call TXCALC(0)
!!$    call TXGLOB
!!$    call TXWDAT
!!$    call TXWDAT2

!!$  do nr = 0, nrmax
!!$     write(6,*) real(RHO(NR)), GY%v(NR,NGR,1), GY%v(NR,NGR,37), GY%v(NR,NGR,14)*1.e3, GY%v(NR,NGR,15)*1.e3, GY%v(NR,NGR,35), GY%v(NR,NGR,36)
!!$  end do

  RETURN
end SUBROUTINE TXGLOD

!***************************************************************
!
!   Read ASCII data file to substitute it into array
!
!***************************************************************

subroutine ascii_input

  use tx_commons, only : infiles, nmax_file, n_infiles, iflag_file, datatype
  use tx_interface, only : KSPLIT_TX
  implicit none
  integer(4) :: IST, i, j, k, nrho, i_start, jshift
  character(len=100) :: TXFNAM
  character(len=180) :: kline, kline1, kline2
  character(len=20)  :: kmesh, kdata
  integer(4) :: nol, itype
  LOGICAL :: LEX

  ! Check a mode of input file
  do 
     write(6,*) '# NBI input from OFMC (1)(2) or arbitrary input (3) ?'
     flush(6)
     read(*,'(I1)',iostat=ist) iflag_file
     if(ist > 0) then
        cycle
     else if(ist < 0) then
        return
     end if
     if(iflag_file >= 1 .and. iflag_file <= 3) then
        exit
     else
        cycle
     end if
  end do

  ! *** Pre-defined input (OFMC: OrbitEffectDist.dat) ***************
  if(iflag_file == 1) then
     TXFNAM = 'OrbitEffectDist.dat'
     inquire(file=TXFNAM,exist=LEX)
     if (LEX) then
     else
        write(6,*) 'XX  FILE ( ', TXFNAM(1:len_trim(TXFNAM)), ' ) DOES NOT exist !'
        return
     end if

     n_infiles = 8

  ! *** Pre-defined input (OFMC: Torque.txt) ************************
  else if(iflag_file == 2) then
     TXFNAM = 'Torque.txt'
     inquire(file=TXFNAM,exist=LEX)
     if (LEX) then
     else
        write(6,*) 'XX  FILE ( ', TXFNAM(1:len_trim(TXFNAM)), ' ) DOES NOT exist !'
        return
     end if

     n_infiles = 1 ! only the TOTAL torque

  ! *** Arbitrary input *********************************************
  else if(iflag_file == 3) then
     ! Check ASCII file
     do
        write(6,*) '# INPUT : LOAD ASCII FILE NAME'
        flush(6)
        read(*,'(A100)',iostat=IST) TXFNAM
        if (IST > 0) then
           cycle
        else if (IST < 0) then
           return
        end if
        inquire(file=TXFNAM,exist=LEX)
        if (LEX) then
           exit
        else
           write(6,*) 'XX  FILE ( ', TXFNAM(1:len_trim(TXFNAM)), ' ) DOES NOT exist !'
        end if
     end do

     ! Number of data
     do
        write(6,*) '# Number of data which you would like to use as inputs ?'
        flush(6)
        read(*,'(I2)',iostat=ist) n_infiles
        if(ist > 0) then
           cycle
        else if(ist < 0) then
           return
        else
           exit   
        end if
     end do
  end if

  ! Deallocate if already allocated
  if(allocated(infiles)) deallocate(infiles)

  ! Allocate derived type
  allocate(infiles(1:n_infiles))

  ! *** Pre-defined input (OFMC: OrbitEffectDist.dat) ***************
  if(iflag_file == 1) then
     ! Read data
     open(21,file=TXFNAM,iostat=IST,status='old',form='formatted')
     if (IST == 0) then
        write(6,*) '# ASCII FILE ( ', TXFNAM(1:len_trim(TXFNAM)),  &
             &     ' ) IS ASSIGNED FOR INPUT.'
     else if (IST > 0) then
        write(6,*) 'XX  ASCII FILE open ERROR !, iostat = ', IST
     else if (IST < 0) then
        deallocate(infiles)
        return
     end if

     ! Name
     !   In case of "NBI input from OFMC (1)", sequence of data is already
     !   defined as follows:
     !     (1) S_birth_ele,  (2) S_birth_tot, (3) S_birth_trap, (4) S_birth_pass,
     !     (5) S_birth_loss, (6) S_orbit_tot, (7) S_orbit_trap, (8) S_orbit_pass
     !   Therefore name of the data is not necessary; hence null is substituted.
     !   Note: The following relations are satisfied:
     !     (1) = (3) + (4) + (5), (2) = (3) + (4), (6) = (7) + (8)
     !     (3) /= (7), (4) /= (8), (2) = (6)
     infiles(1:n_infiles)%name = ' '

     ! Read data from the file
     i    = 0  ! outermost loop count
     nol  = 0  ! number of lines
     nrho = 0  ! number of times that "RHO" appears in kline, i.e. separator
     k    = 0  ! data is taken during k/=0
     do 
        read(21,'(A)',iostat=IST) kline
        if(ist > 0) then
           stop 'Read file error'
        else if(ist < 0) then
           exit ! detect the end of the file
        end if
        i = i + 1
        if(k == 0) then
           if(index(kline,"RHO") /= 0) then ! detect the start position of the data chunk
              nrho = nrho + 1
              k = 1
              i_start = i + 1

              jshift = 2 ! denotes the index "I" and the radial coordinates "RHO"
              if(index(kline,"PSI") /= 0) then ! detect whether the file is new or old
                 !== jshift ================================================================
                 !  Old type "OrbitEffctDist.dat" includes only the "RHO" coordinate.
                 !  New type "OrbitEffctDist.dat" includes the "PSI" and "RHO" coordinates.
                 !==========================================================================
                 jshift = 3 ! denotes the radial coordinates "PSI", which will be discarded
                            !   as well as "I" and "RHO"
              end if
           end if
           cycle
        end if

        if(nrho == 1) then ! === Get "TOTAL" and "TOTAL POWER" data ===
           if(index(kline,"sec") /= 0) then
              ! Total number of ions per second
              j = 0
              do
                 kline = trim(adjustl(kline))
                 call ksplit_tx(kline,' ',kline1,kline2)
                 if(index(kline1,"E+") == 0) then
                    kline = kline2
                    cycle
                 else
                    j = j + 1
                    kdata = trim(kline1)
                    read(kdata,'(E15.7)') infiles(j)%totS
                    kline = kline2
                 end if
                 if(len_trim(kline2) == 0) exit ! all the data has been already taken.
              end do

           else if(index(kline,"W") /= 0) then
              ! Total power of ions
              j = 0
              do
                 kline = trim(adjustl(kline))
                 call ksplit_tx(kline,' ',kline1,kline2)
                 if(index(kline1,"E+") == 0) then
                    kline = kline2
                    cycle
                 else
                    j = j + 1
                    kdata = trim(kline1)
                    read(kdata,'(E15.7)') infiles(j)%totP
                    kline = kline2
                 end if
                 if(len_trim(kline2) == 0) then
                    exit ! all the data has been already taken.
                 end if
              end do
              
              k = 0 ! flag for the next data chunk
           end if

       else if(nrho == 2) then ! === Get "SMOOTHING" data ===
           if(len_trim(kline) == 0) then ! detect the end of the data chunk
              k = 0
           end if
           nol = i - i_start + 1

           do j = 1, n_infiles + jshift
              ! jshift denotes the index number and the mesh data columns
              kline = trim(adjustl(kline))
              call ksplit_tx(kline,' ',kline1,kline2)
              if(j == 1) then
                 ! First line (index I) is discarded.
                 kline = kline2
                 cycle
              else if(j == jshift) then
                 ! Second (jshift == 0) or third (jshift == 1) line is defined as a mesh data, rho.
                 kmesh = trim(kline1)
                 kline = kline2
              else if(j > jshift) then
                 ! Third to eighth lines are defined as data, in case of jshift == 0
                 ! Fourth to ninth lines are defined as data, in case of jshift == 1
                 ! Please see "Name" shown above
                 kdata = trim(kline1)
                 read(kmesh,'(E15.7)') infiles(j-jshift)%r(nol)
                 read(kdata,'(E15.7)') infiles(j-jshift)%data(nol)
                 kline = kline2
              else
                 kline = kline2
              end if
           end do

        else if(nrho == 4) then ! === Get "V//" data ===
           if(len_trim(kline) == 0) then ! detect the end of the data chunk
              k = 0
           else
              nol = i - i_start + 1
              j = 5 ! Vb data exist for (6)(7)(8) only
              infiles(1:j)%vb(nol) = 0.d0
              do
                 kline = trim(adjustl(kline))
                 call ksplit_tx(kline,' ',kline1,kline2)
                 if(index(kline1,"E+") == 0) then
                    kline = kline2
                    cycle
                 else
                    j = j + 1
                    ! Parallel velocity
                    kdata = trim(kline1)
                    read(kdata,'(E15.7)') infiles(j)%vb(nol)
                    kline = kline2
                 end if
                 if(len_trim(kline2) == 0) exit ! all the data has been already taken.
              end do
           end if

        else
           ! detect the end of the data chunk
           if(len_trim(kline) == 0) k = 0

        end if

     end do

     if(nol > nmax_file) then
        write(6,*) 'XX The number of lines in the file exceeds ',nmax_file
        write(6,*) '   Program is terminated.'
        stop
     else
        ! All the data should have same mesh data.
        infiles(:)%nol = nol
     end if

     close(21)

!!$     write(6,*) "===== data ====="
!!$     do i = 1, nol
!!$        write(6,'(I3,F6.3,8(ES12.4))') i,infiles(1)%r(i),(infiles(j)%data(i),j=1,8)
!!$     end do
!!$     write(6,*) "===== V// ====="
!!$     do i = 1, nol
!!$        write(6,'(I3,F6.3,8(ES12.4))') i,infiles(1)%r(i),(infiles(j)%vb(i),j=1,8)
!!$     end do
!!$     write(6,*) "===== total S ====="
!!$     write(6,'(9X,8(ES12.4))') (infiles(i)%totS, i = 1, 8)
!!$     write(6,*) "===== total W ====="
!!$     write(6,'(9X,8(ES12.4))') (infiles(i)%totP, i = 1, 8)
!!$     stop

  ! *** Pre-defined input (OFMC: Torque.txt) ************************
  else if(iflag_file == 2) then
     ! Read data
     open(21,file=TXFNAM,iostat=IST,status='old',form='formatted',position='rewind')
     if (IST == 0) then
        write(6,*) '# ASCII FILE ( ', TXFNAM(1:len_trim(TXFNAM)),  &
             &     ' ) IS ASSIGNED FOR INPUT.'
     else if (IST > 0) then
        write(6,*) 'XX  ASCII FILE open ERROR !, iostat = ', IST
     else if (IST < 0) then
        deallocate(infiles)
        return
     end if

     infiles(1:n_infiles)%name = 'TOTAL'

     k  = 0  ! data is taken during k/=0
     do 
        read(21,'(A)',iostat=IST) kline
        if(ist > 0) then
           stop 'Read file error'
        else if(ist < 0) then
           exit ! detect the end of the file
        end if
        
        if(k == 0) then
           if(index(kline,"I,") /= 0) then ! detect the start position of the data chunk
              k = 1
              nol = 0
           end if
           cycle
        end if

        j = 0
        do
           kline = trim(adjustl(kline))
           if(kline == ' ') then
              k = 0
              exit
           end if
           call ksplit_tx(kline,',',kline1,kline2)
           if(index(kline1,"E") == 0) then
              kdata = trim(kline1)
              if(kdata == ' ' .or. kdata == "TOTAL") then
                 nol = 0
              else
                 nol = nol + 1 ! Not use "I" index
!                 read(kdata,'(I3)') nol ! Use "I" index
              end if
              kline = kline2
              cycle
           else
              j = j + 1
              kdata = trim(kline1)
!              write(6,*) "nol",nol,kdata,"@",kline2,"@"
              if(j == 1 .and. nol /= 0) then
                 infiles(1)%nol = nol
!                 write(6,*) nol,kdata
                 if(nol /= 0) read(kdata,'(E15.7)') infiles(1)%r(nol)
!                 write(6,*) "rho",nol,infiles(1)%r(nol)
                 kline = kline2
              else
                 if(index(kline2,"E") /= 0) then
                    kdata = trim(kline2)
                    if(nol /= 0) then
                       read(kdata,'(E15.7)') infiles(1)%data(nol)
!                       write(6,*) "data",nol,infiles(1)%data(nol)
                    else
                       read(kdata,'(E15.7)') infiles(1)%totS
!                       write(6,*) "total",nol,infiles(1)%totS
                    end if
                    kline = kline2
                 else
                    exit
                 end if
              end if
!              pause
              cycle
           end if
        end do
     end do

!!$     write(6,*) infiles(1)%nol
!!$     do nol = 1, infiles(1)%nol
!!$        write(6,*) nol,infiles(1)%r(nol),infiles(1)%data(nol)
!!$     end do
!!$     write(6,*) infiles(1)%totS

     close(21)

  ! *** Arbitrary input *********************************************
  else if(iflag_file == 3) then

     ! Read data
     do i = 1, n_infiles
        open(21,file=TXFNAM,iostat=IST,status='old',form='formatted')
        if (IST == 0) then
           write(6,*) '# ASCII FILE ( ', TXFNAM(1:len_trim(TXFNAM)),  &
                &     ' ) IS ASSIGNED FOR INPUT.'
        else if (IST > 0) then
           write(6,*) 'XX  ASCII FILE open ERROR !, iostat = ', IST
        else if (IST < 0) then
           deallocate(infiles)
           return
        end if

        ! Select a type of input data from ascii files
        do 
           write(6,'(A,I2,A)') '# Select a type of input data from ascii file for ',i,' data.'
           write(6,*) '# 1: PNBP, 2: PNBT1, 3: PNBT2, 4: PRF, 5: LQe4'
           flush(6)
           read(*,'(I2)',iostat=ist) itype
           if(ist > 0) then
              cycle
           else if(ist < 0) then
              deallocate(infiles)
              return
           else
              if(itype <= 0 .or. itype > 5) cycle
              exit
           end if
        end do
        infiles(i)%name = datatype(itype)

        ! Which column indicates the mesh data
        do
           write(6,*) '# Which column indicates the mesh data in ', &
                &     TXFNAM(1:len_trim(TXFNAM)), ' ?'
           read(*,'(I2)', iostat=ist) infiles(i)%ncol_mesh
           if(infiles(i)%ncol_mesh == 0) cycle
           if (ist > 0) then
              cycle
           else if(ist < 0) then
              deallocate(infiles)
              return
           else
              exit
           end if
        end do

        ! Which column indicates the data which you would like to use
        do
           write(6,*) '# Which column indicates the data which you would like to use in ', &
                &     TXFNAM(1:len_trim(TXFNAM)), ' ?'
           read(*,'(I2)', iostat=ist) infiles(i)%ncol_data
           if(infiles(i)%ncol_data == 0) cycle
           if (ist > 0) then
              cycle
           else if(ist < 0) then
              deallocate(infiles)
              return
           else
              exit
           end if
        end do

        ! Read data from the file
        nol = 0 ! number of lines
        do 
           read(21,'(A)',iostat=IST) kline
           if(ist < 0) exit ! detect the end of the file
           nol = nol + 1

           j = 0
           do
              j = j + 1
              kline = trim(adjustl(kline))
              call ksplit_tx(kline,' ',kline1,kline2)
              if(j == infiles(i)%ncol_mesh) then
                 kmesh = trim(kline1)
              else if(j == infiles(i)%ncol_data) then
                 kdata = trim(kline1)
              end if
              kline = kline2
              if(j >= infiles(i)%ncol_mesh .and. j >= infiles(i)%ncol_data) exit
           end do

           read(kmesh,'(E15.7)') infiles(i)%r(nol)
           read(kdata,'(E15.7)') infiles(i)%data(nol)

!           write(6,*) infiles(i)%r(nol),infiles(i)%data(nol)
        end do
        if(nol > nmax_file) then
           write(6,*) 'XX The number of lines in the file exceeds ',nmax_file
           write(6,*) '   Program is terminated.'
           stop
        else
           infiles(i)%nol = nol
        end if

        close(21)

     end do

  end if

end subroutine ascii_input

!***************************************************************
!
!   Return the number of the datatype in "infiles"
!
!***************************************************************

integer(4) function detect_datatype(kchar)

  use tx_commons, only : infiles, n_infiles
  use libchar, only : kmatch
  implicit none
  character(len=*), intent(in) :: kchar
  integer(4) :: i

  do i = 1, n_infiles
     if(kmatch(infiles(i)%name,kchar)) then
        detect_datatype = i
        return
     end if
  end do

  detect_datatype = 0

end function detect_datatype

!***************************************************************
!
!   Control routine for writing ASCII data file
!
!      Just choose "for_ntmain" or "for_ofmc" routine
!
!***************************************************************

subroutine outfile
  implicit none
  integer(4) :: n, ist

  do
     write(6,'(1X,A,1X,I1)') '## outfile: input ?'
     write(6,'(1X,A)') '##          1: profiles, 2: 0d quantities, 3: 0d variables, ' &
          & // '4: for OFMC, 5: for TOPICS/NTMAIN, 9: exit'
     read(*,'(I1)',iostat=ist) n
     if (ist > 0) then
        cycle
     else if (ist < 0) then
        return
     end if

     select case(n)
     case(1)
        call write_profs
     case(2)
        call write_quantities
     case(3)
        call write_vars
     case(4)
        call for_ofmc
     case(5)
        call for_ntmain
     case default
        exit
     end select
  end do

end subroutine outfile

!***************************************************************
!
!   Write ASCII data file in terms of profiles
!
!***************************************************************

subroutine write_profs
  use tx_commons, only : NRMAX
  use tx_graphic, only : GT, GY, NGR, NGYRM
  use libfio, only : FWopen
  implicit none
  integer(4) :: NG, NR, NGYR, IERR

  call FWopen(21,'tx_profs.dat',1,1,'TX',IERR)
  if(IERR /= 0) then
     write(6,*) 'XX write_profs: FWopen: IERR=', IERR
     return
  end if

  do NG = 0, NGR
     do NR = 1, NRMAX
        write(21,'(256ES15.7)')  GT(NG), (GY%v(NR,NG,NGYR), NGYR=1,NGYRM)
     end do
     write(21,*) 
  end do

  close(21)

end subroutine write_profs

!***************************************************************
!
!   Write ASCII data file in terms of 0d quantities
!
!***************************************************************

subroutine write_quantities
  use tx_graphic, only : GTX, GTY, NGT, NGYTM
  use libfio, only : FWopen
  implicit none
  integer(4) :: NG, NGYT, IERR

  call FWopen(21,'tx_quantities.dat',1,1,'TX',IERR)
  if(IERR /= 0) then
     write(6,*) 'XX write_quantities: FWopen: IERR=', IERR
     return
  end if

  do NG = 0, NGT
     write(21,'(256ES15.7)')  GTX(NG), (GTY(NG,NGYT), NGYT=1,NGYTM)
  end do

  close(21)

end subroutine write_quantities

!***************************************************************
!
!   Write ASCII data file in terms of 0d variables
!
!***************************************************************

subroutine write_vars
  use tx_graphic, only : GVX, GVY, NGVV, NGYVM
  use libfio, only : FWopen
  implicit none
  integer(4) :: NG, NGYV, IERR

  call FWopen(21,'tx_vars.dat',1,1,'TX',IERR)
  if(IERR /= 0) then
     write(6,*) 'XX write_vars: FWopen: IERR=', IERR
     return
  end if

  do NG = 0, NGVV
     write(21,'(256ES15.7)')  GVX(NG), (GVY(NG,NGYV), NGYV=1,NGYVM)
  end do

  close(21)

end subroutine write_vars

!***************************************************************
!
!   Write ASCII data file for TOPICS/NTMAIN
!
!***************************************************************

subroutine for_ntmain
  use tx_commons, only : NRMAX, Rho, PN01V, PN02V, Var
  use libfio, only : FWopen
  implicit none
  integer(4) :: IERR, NR

  call FWopen(21,'tx_ntmain.dat',1,1,'TX',IERR)
  if(IERR /= 0) then
     write(6,*) 'XX for_ntmain: FWopen: IERR=', IERR
     return
  end if

  write(21,'(4X,I3)') NRMAX+1
  do nr = 0, nrmax
     write(21,'(7ES17.9)') Rho(NR), Var(NR,1)%n*1.D20, Var(NR,2)%n*1.D20, &
          & Var(NR,1)%T*1.D3, Var(NR,2)%T*1.D3, PN01V(NR)*1.D20, PN02V(NR)*1.D20
  end do

  close(21)

end subroutine for_ntmain

!***************************************************************
!
!   Write ASCII data file for OFMC
!
!      Prepare three profiles of electron density and
!         electron and ion temperatures from TASK/TX
!
!      For OFMC code, the ele. density should be in m^{-3} and
!                     the temperatures in eV.
!
!      This routine was originally built as task/data_spline.
!
!***************************************************************

subroutine for_ofmc
  use tx_commons, only : NRMAX, NRA, RR, RA, BB, rIp, psiV, Var, array_init_NR
  use libfio, only : FWopen
  use libspl1d, only : spl1d, spl1df
  implicit none

  integer(4) :: nmax, ist, i, j, NR, ierr
  real(8) :: dpsi, psirho_a
  real(8), dimension(:), allocatable :: psirho, psi_out, data, data_out, deriv
  real(8), dimension(:,:), allocatable :: u

  do
     write(6,'(1X,A)') '## for_ofmc: number of grids ? [Default=30]'
     read(*,'(I3)',iostat=ist) nmax
     if (ist > 0) then
        cycle
     else if (ist < 0) then
        return
     end if
     if(nmax == 0) then
        nmax = 30
     else if(nmax < 0) then
        cycle
     end if
     exit
  end do

  ! Output file open
  call FWopen(21,'tx_ofmc.dat',1,1,'TX',IERR)
  if(IERR /= 0) then
     write(6,*) 'XX for_ofmc: FWopen: IERR=', IERR
     return
  end if

  write(21,'(A)') 'PARM'
  write(21,'(1X,A5,1X,I5,2(1X,A5))') '*****',-1,('*****', j = 1, 2)
  write(21,'(1X,I5,11(1X,A5))') nmax,('*****', j = 1, 11)
  write(21,'(4(ES12.5))') RR,RA,BB,rIp
  write(21,'(3(1X,A11),2(ES12.5))') ('************', j = 1, 3), 1.d0, 0.d0
  write(21,'(5(1X,A11))') ('************', j = 1, 5)

  allocate(psirho, data, deriv, source=array_init_NR)
  allocate(u(1:4,0:NRMAX))
  allocate(psi_out(1:nmax),data_out(1:nmax+1))

  ! normalized psi with respect to rho
  do NR = 1, NRMAX
     psirho(NR) = psiV(NR) - psiV(0) !- RR * (AphV(NR) - AphV(0))
  end do
  psirho_a = psirho(NRA)
  psirho(:) = psirho(:) / psirho_a

  ! Equally-spaced psi on a half mesh
  dpsi = 1.d0 / real(nmax,8)
  do i = 1, nmax
     psi_out(i) = i * dpsi - 0.5d0 * dpsi
  end do

  do i = 1, 3
     if     (i == 1) then ! Ne
        data(:) = Var(:,1)%n * 1.D20 ! in m^{-3}
     else if(i == 2) then ! Te
        data(:) = Var(:,1)%T * 1.D3  ! in eV
     else if(i == 3) then ! Ti
        data(:) = Var(:,2)%T * 1.D3  ! in eV
     end if

     ! Output data on the magnetic axis
     write(21,'(2(ES12.5))') data(0), -1.d0

     ! Spline data for equally-spaced psi
     call spl1d(psirho,data,deriv,u,NRMAX+1,0,ierr)
     if(ierr /= 0) stop 'Error at spl1d'
     do j = 1, nmax
        call spl1df(psi_out(j),data_out(j),psirho,u,NRMAX+1,ierr)
        if(ierr /= 0) stop 'Error at spl1df'
     end do

     ! Add data on the separatrix to output data array
     data_out(nmax+1) = data(NRA)

     ! Output data
     write(21,'(4(ES18.10))') (data_out(j), j = 1, nmax+1)
  end do

  deallocate(psirho,data,deriv,u,psi_out,data_out)

  close(21)

end subroutine for_ofmc

!***************************************************************
!
!   Read ASCII data file for initial profiles
!
!      Prepare profiles of 
!         PNeV          : electron density [/m^3], 
!         PTeV, PTiV    : electron and ion temperatures [eV], and
!         AEE*PNeV*UephV: toroidal current density [A/m^2]
!         for TASK/TX
!
!      List structure has been used.
!
!***************************************************************

subroutine initprof_input(nr, idx, out)

  use tx_commons, only : Rho, AEE
  use libitp, only : aitken2p, fctr4pt
  use libfio, only : fropen
  use libspl1d, only : spl1d, spl1df

  integer(4), optional, intent(in) :: nr, idx
  real(8), optional, intent(out) :: out
  integer(4) :: nintin, ier, ist, k, iflag, j
  integer(4), save :: nrinmax
  real(8), dimension(:), allocatable, save :: rho_in, deriv
  real(8), dimension(:,:), allocatable, save :: prof_in, u1, u2, u3, u4

  type unit
     integer(4)          :: l_p
     real(8)             :: rho_p, prof1_p, prof2_p, prof3_p, prof4_p
     type(unit), pointer :: next
  end type unit
  type(unit), pointer :: ent, new, p
  integer(4) :: l
  real(8)    :: rho_r, prof1_r, prof2_r, prof3_r, prof4_r

  if((present(nr) .eqv. .false.) .and. (present(idx) .eqv. .false.)) then
     allocate(ent); nullify(ent%next)

     nintin = 24
     call FRopen(nintin,'initprof.dat',1,0,'INIT PROF',ier)
     if(ier /= 0) return

     do
!!!        read(nintin,'(1X,I3,5(E9.3))',iostat=ist) l, rho_r, prof1_r, prof2_r, prof3_r, prof4_r
        read(nintin,'(1X,I3,5(E10.3))',iostat=ist) l, rho_r, prof1_r, prof2_r, prof3_r, prof4_r
        if(ist > 0) then
           write(6,*) 'XX file format is invalid. Aborting.'
           close(nintin)
           stop
        else if (ist < 0) then
           exit
        end if
        call rearrange
     end do

     iflag = 0
     p => ent%next
     k = 0
     do
        if(associated(p)) then 
!           write(6,*) p%l_p, p%rho_p, p%prof1_p, p%prof2_p, p%prof3_p, p%prof4_p
           if(k == 0 .and. p%rho_p /= 1.d0) iflag = iflag + 1
           if(p%l_p == 1 .and. p%rho_p /= 0.d0) iflag = iflag + 2
           k = k + 1
           p => p%next
        else
           exit
        end if
     end do
     if(iflag == 0) then
        nrinmax = k
     else if(iflag == 1 .or. iflag == 2) then
        nrinmax = k + 1
     else
        nrinmax = k + 2
     end if

     allocate(rho_in(1:nrinmax),prof_in(1:nrinmax,1:4))
     p => ent%next
     do
        if(associated(p)) then
           rho_in (p%l_p)   = p%rho_p
           prof_in(p%l_p,1) = p%prof1_p
           prof_in(p%l_p,2) = p%prof2_p
           prof_in(p%l_p,3) = p%prof3_p
           prof_in(p%l_p,4) = p%prof4_p
           p => p%next
        else
           exit
        end if
     end do

     if(iflag == 1) then ! Extrapolate values at r/a = 1.
        rho_in(k+1) = 1.d0
        do j = 1, 4
           prof_in(nrinmax,j) = AITKEN2P(rho_in(nrinmax),prof_in(nrinmax-1,j), &
                & prof_in(nrinmax-2,j),prof_in(nrinmax-3,j),rho_in(nrinmax-1), &
                & rho_in(nrinmax-2),rho_in(nrinmax-3))
        end do
     else if(iflag == 2) then ! Extrapolate values at r/a = 0.
        rho_in(2:k+1) = rho_in(1:k)
        rho_in(1) = 0.d0
        prof_in(2:k+1,1:4) = prof_in(1:k,1:4)
        do j = 1, 4
           prof_in(1,j) = FCTR4pt(rho_in(2),rho_in(3),rho_in(4), &
                &                 prof_in(2,j),prof_in(3,j),prof_in(4,j))
        end do
     else if(iflag == 3) then ! Extrapolate values at r/a = 0 and 1.
        rho_in(2:k+1) = rho_in(1:k)
        rho_in(1) = 0.d0
        rho_in(k+2) = 1.d0
        prof_in(2:k+1,1:4) = prof_in(1:k,1:4)
        do j = 1, 4
           prof_in(1,j) = FCTR4pt(rho_in(2),rho_in(3),rho_in(4), &
                &                 prof_in(2,j),prof_in(3,j),prof_in(4,j))
           prof_in(nrinmax,j) = AITKEN2P(rho_in(nrinmax),prof_in(nrinmax-1,j), &
                & prof_in(nrinmax-2,j),prof_in(nrinmax-3,j),rho_in(nrinmax-1), &
                & rho_in(nrinmax-2),rho_in(nrinmax-3))
        end do
     end if

!!$     do k = 1, nrinmax
!!$        write(6,'(1X,I3,5ES10.3)') k, rho_in(k),(prof_in(k,j),j=1,4)
!!$     end do
!!$
     close(nintin)
     if(associated(ent)) deallocate(ent)
     if(associated(new)) deallocate(new)

     allocate(deriv(1:nrinmax),u1(1:4,1:nrinmax),u2(1:4,1:nrinmax),u3(1:4,1:nrinmax),u4(1:4,1:nrinmax))
     call spl1d(rho_in,prof_in(1:nrinmax,1),deriv,u1,nrinmax,0,ier)
     if(ier /= 0) stop 'Error at spl1d in initprof_input: prof_in(1)'
     call spl1d(rho_in,prof_in(1:nrinmax,2),deriv,u2,nrinmax,0,ier)
     if(ier /= 0) stop 'Error at spl1d in initprof_input: prof_in(2)'
     call spl1d(rho_in,prof_in(1:nrinmax,3),deriv,u3,nrinmax,0,ier)
     if(ier /= 0) stop 'Error at spl1d in initprof_input: prof_in(3)'
     call spl1d(rho_in,prof_in(1:nrinmax,4),deriv,u4,nrinmax,0,ier)
     if(ier /= 0) stop 'Error at spl1d in initprof_input: prof_in(4)'

  else if((present(nr) .eqv. .true.) .and. (present(idx) .eqv. .true.)) then
     if(idx == 1) then ! Electron density
        call spl1df(Rho(nr),out,rho_in,u1,nrinmax,ier)
        if(ier /= 0) stop 'Error at spl1df in initprof_input: idx = 1'
        out = out * 1.d-20
     else if(idx == 2) then ! Electron temperature
        call spl1df(Rho(nr),out,rho_in,u2,nrinmax,ier)
        if(ier /= 0) stop 'Error at spl1df in initprof_input: idx = 2'
        out = out * 1.d-3
     else if(idx == 3) then ! Ion temperature
        call spl1df(Rho(nr),out,rho_in,u3,nrinmax,ier)
        if(ier /= 0) stop 'Error at spl1df in initprof_input: idx = 3'
        out = out * 1.d-3
     else if(idx == 4) then ! Toroidal current
        call spl1df(Rho(nr),out,rho_in,u4,nrinmax,ier)
        if(ier /= 0) stop 'Error at spl1df in initprof_input: idx = 4'
        out = out
     else 
        stop 'initprof_input: wrong input!'
     end if

  else if((present(nr) .eqv. .false.) .and. (present(idx) .eqv. .true.)) then
     if(idx == 0) then
        deallocate(rho_in,prof_in,deriv,u1,u2,u3,u4)
     else
        stop 'initprof_input: wrong input!'
     end if
  else
     stop 'initprof_input: wrong input!'
  end if

contains
  subroutine rearrange
    allocate(new)
    new = unit(l, rho_r, prof1_r, prof2_r, prof3_r, prof4_r, ent%next)
    ent%next => new
  end subroutine rearrange

end subroutine initprof_input
