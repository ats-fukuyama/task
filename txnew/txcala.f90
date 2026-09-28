module tx_coefficients
  use tx_commons
  use tx_core_module
  implicit none
  private
  real(8), parameter :: fourPisq = 4.d0 * Pi * Pi
  real(8), dimension(:,:,:,:), allocatable :: ELM
  real(8), dimension(:), allocatable ::  &
       & DTf
  real(8), dimension(:), allocatable :: &
       & rrtinv, aatinv, bthcoinv, ParatoZeta, ZetatoPara, &
       & aatq, PsitdotVBth, fipolinv, fipolsdt, dZetatoPara, ribi, ribsdt, &
       & RatCXSNBe, RatCXSNBi, Ratdenhyd
!!rp_conv       &, rNubLL
  real(8), dimension(:,:), allocatable :: UsrgV, fipolsdtPNsV, BUsparVbbt, &
       & dlnPNsV, BpBVsdiag, tormflux, dtormflux, ribsdtNs!, RatPNbinv,
  real(8), dimension(:,:,:), allocatable :: Ratden, zChi
  real(8) :: DTt, invDT, BeamSW, RpplSW, ThntSW, &
       &     fact = 1.d0 ! <= SOL loss accelerator
  integer(4), save :: ICALA = 0!, ICALA2 = 0
  integer(4) :: L3, L6, L8
  ! Density equation of each thermal species, indexed the same way as Var(:,i).
  integer(4), parameter :: LQsn(1:NSM) = [ LQe1, LQi1, LQz1 ]
  ! NOTE:
  !   These workspaces are module-wide shared arrays.
  !   Current TXLOOP usage is serial-safe, but OpenMP parallel calls to TXCALA
  !   would require thread-private workspaces.
  
  public :: TXCALA

contains

!***************************************************************
!
!   Calculate ALC, BLC, CLC, PLC
!
!***************************************************************

  subroutine TXCALA

    integer(4) :: NC, NQ, iHvsLC, NSMP

    !*** Nodal Equation *****************************************!
    !   ALC  : Coefficient matrix on NR+1                        !
    !   BLC  : Coefficient matrix on NR                          !
    !   CLC  : Coefficient matrix on NR-1                        !
    !   NLCR : Matrix of identifying variables with NR index     !
    !   NLC  : Matrix of identifying variables                   !
    !   PLC  : Coefficient matrix for non-variable terms         !
    !   NLCMAX : Number of right-hand-side terms                 !
    !************************************************************!

    !*** Elemental Equation *************************************!
    !   ELM  : Elemental matrix of one term, equation, element   !
    !************************************************************!

    NSMP = NSM + 1
    call ensure_txcala_workspace(NSMP)

    ! Replace per-iteration allocate(...,source=0) with explicit clear.
    ELM = 0.d0
    rrtinv = 0.d0
    aatinv = 0.d0
    bthcoinv = 0.d0
    ParatoZeta = 0.d0
    ZetatoPara = 0.d0
    aatq = 0.d0
    PsitdotVBth = 0.d0
    fipolinv = 0.d0
    fipolsdt = 0.d0
    dZetatoPara = 0.d0
    ribi = 0.d0
    ribsdt = 0.d0
    RatCXSNBe = 0.d0
    RatCXSNBi = 0.d0
    Ratdenhyd = 0.d0
    UsrgV = 0.d0
    fipolsdtPNsV = 0.d0
    BUsparVbbt = 0.d0
    dlnPNsV = 0.d0
    BpBVsdiag = 0.d0
    tormflux = 0.d0
    dtormflux = 0.d0
    ribsdtNs = 0.d0
    Ratden = 0.d0
    zChi = 0.d0

    !     Preconditioning

    invDT = 1.d0 / DT
!!$    if(DT <= 2.D-5) then
!!$       DTt          = 1.5d1
!!$       DTf(1)       = 1.d0
!!$       DTf(2:NQMAX) = 1.5d1
!!$    else
    DTt          = 1.d0
    DTf(1:NQMAX) = 1.d0
!!$    end if

    ! In case of ohmic heating (i.e. no NBI), largeness of the term related to beam
    ! components comparable to that of other terms related to finite variables
    ! sometimes induces spurious values in beam components even if no NBI is activated.
    ! To avoid this, setting these terms, especially related to electromagnetic potentials,
    ! to zero is desirable during no NBI.

    if(PNBH == 0.d0 .and. PNbV(0) == 0.d0) then
       BeamSW = 0.d0
    else
       BeamSW = 1.d0
    end if
    if(FSRP == 1.d0) then
       RpplSW = BeamSW
    else
       RpplSW = 0.d0
    end if

    ! If ThntSW = 0, there is no particle source from PN02V.
    ThntSW = 1.d0

    ! SUPG
    L3 = abs(iSUPG3) * 100
    L6 =     iSUPG6  * 100
    L8 =     iSUPG8  * 100

    ! Initialize arrays
    ALC  = 0.d0
    BLC  = 0.d0
    CLC  = 0.d0
    NLCR = 0
    NLC  = 0
    PLC  = 0.d0
    NLCMAX = 0

    call LQcoef

    !     Maxwell

    call LQm1CC
    call LQm2CC
    call LQm3CC
    call LQm4CC
    call LQm5CC

    !     Electron

    call LQe1CC
    call LQe2CC
    call LQe3CC
    call LQe4CC
    call LQe5CC
    call LQe6CC ! Heat flux
    call LQe7CC ! Rotation frequency
    call LQe8CC ! Diamagnetic particle flow

    !     Ion

    call LQi1CC
    call LQi2CC
    call LQi3CC
    call LQi4CC
    call LQi5CC
    call LQi6CC ! Heat flux
    call LQi7CC ! Rotation frequency
    call LQi8CC ! Diamagnetic particle flow

    !     Impurity

    call LQz1CC
    call LQz2CC
    call LQz3CC
    call LQz4CC
    call LQz5CC
    call LQz6CC ! Heat flux
    call LQz7CC ! Rotation frequency
    call LQz8CC ! Diamagnetic particle flow

    !     Beam

    call LQb1CC
    call LQb2CC
    call LQb3CC
    call LQb4CC
    call LQb7CC
    call LQb8CC

    !     Neutral

    call LQn1CC
    call LQn2CC
    call LQn3CC
    call LQnzCC

    !     Ripple trapped beam

    call LQr1CC

    !     Elemental equations -> Nodal equations

    do NQ = 1, NQMAX
       NC = 0
          BLC(0:NRMAX-1,NC,NQ) = BLC(0:NRMAX-1,NC,NQ) + ELM(1:NEMAX,1,NC,NQ) * DTt
          ALC(0:NRMAX-1,NC,NQ) = ALC(0:NRMAX-1,NC,NQ) + ELM(1:NEMAX,2,NC,NQ) * DTt
          CLC(1:NRMAX  ,NC,NQ) = CLC(1:NRMAX  ,NC,NQ) + ELM(1:NEMAX,3,NC,NQ) * DTt
          BLC(1:NRMAX  ,NC,NQ) = BLC(1:NRMAX  ,NC,NQ) + ELM(1:NEMAX,4,NC,NQ) * DTt

       do NC = 1, NLCMAX(NQ)
          ! iHvsLC : Heaviside function
          !          iHvsLC = 1 when NLC = 0 ; otherwise iHvsLC = 0.
          iHvsLC = real(NQMAX-NLC(NC,NQ))/NQMAX
          BLC(0:NRMAX-1,NC,NQ) = BLC(0:NRMAX-1,NC,NQ) + ELM(1:NEMAX,1,NC,NQ)  * DTf(NQ)
          ALC(0:NRMAX-1,NC,NQ) = ALC(0:NRMAX-1,NC,NQ) + ELM(1:NEMAX,2,NC,NQ)  * DTf(NQ)
          PLC(0:NRMAX-1,NC,NQ) = PLC(0:NRMAX-1,NC,NQ) +(ELM(1:NEMAX,1,NC,NQ) &
               &                                      + ELM(1:NEMAX,2,NC,NQ)) * DTf(NQ) * iHvsLC
          CLC(1:NRMAX  ,NC,NQ) = CLC(1:NRMAX  ,NC,NQ) + ELM(1:NEMAX,3,NC,NQ)  * DTf(NQ)
          BLC(1:NRMAX  ,NC,NQ) = BLC(1:NRMAX  ,NC,NQ) + ELM(1:NEMAX,4,NC,NQ)  * DTf(NQ)
          PLC(1:NRMAX  ,NC,NQ) = PLC(1:NRMAX  ,NC,NQ) +(ELM(1:NEMAX,3,NC,NQ) &
               &                                      + ELM(1:NEMAX,4,NC,NQ)) * DTf(NQ) * iHvsLC
       end do
    end do

    NLCR(0:NCM,1:NQMAX,0) = NLC(0:NCM,1:NQMAX)
    NLCR(0:NCM,1:NQMAX,1) = NLC(0:NCM,1:NQMAX)

    !     Dirichlet condition

!    if(NTCUM >= 3) then
!       if(NTCUM == 3 .and. ICALA2 == 0) then
!          ICALA = 0 ; ICALA2 = 1
!       end if
    call BOUNDARY(0    ,LQe2,0)
    call BOUNDARY(0    ,LQi2,0)
    call BOUNDARY(0    ,LQz2,0)
    ! --- Excess B.C. for numerical stability ---
    !   NOTE: This condition originally should be satisfied by the radial force balance
    !         equation. However, due to numerical accuracy and complex nonlinearlity of
    !         the system, they are not always satisfied. Suppressing the flapping of the
    !         ion parallel flow (LQi3) at the magnetic axis, imposing the following condition
    !         is favorable. It shouldn't be set if the perpendicular viscosity is included
    !         in the parallel flow equations (LQe3, LQi3).
!!unstable    call BOUNDARY(0    ,LQe3,0,Var(0,1)%RUph*bbt(0)/fipol(0))
    if(FSPARV(2) == 0.d0) call BOUNDARY(0    ,LQi3,0,Var(0,2)%RUph*bbt(0)/fipol(0))
    if(FSPARV(3) == 0.d0) call BOUNDARY(0    ,LQz3,0,Var(0,3)%RUph*bbt(0)/fipol(0))
    ! -----------------------------------------------------------
!    end if
    call BOUNDARY(NRMAX,LQm1,0)
    call BOUNDARY(0    ,LQm2,0)
!!!    call BOUNDARY(NRMAX,LQe2,0) ! Violates quasi-neutrality and breaks down calculation
!    call BOUNDARY(NRMAX,LQe3,0)   ! OK
!!!    call BOUNDARY(NRMAX,LQe4,0) ! Not impose this B.C. which violates quasi-neutrality and 
                                   ! break down calculation, because this eq. yields particle flux.
!!!    call BOUNDARY(NRMAX,LQi2,0) ! Violates quasi-neutrality and breaks down calculation
!    call BOUNDARY(NRMAX,LQi3,0)   ! OK
!!!    call BOUNDARY(NRMAX,LQi4,0) ! Not impose this B.C. which violates quasi-neutrality and 
                                   ! break down calculation, because this eq. yields particle flux.
!!!    call BOUNDARY(NRMAX,LQz2,0) ! Violates quasi-neutrality and breaks down calculation
!    call BOUNDARY(NRMAX,LQz3,0)   ! OK
!!!    call BOUNDARY(NRMAX,LQz4,0) ! Not impose this B.C. which violates quasi-neutrality and 
                                   ! break down calculation, because this eq. yields particle flux.

    call BOUNDARY(NRMAX,LQn2,0)
    call BOUNDARY(NRMAX,LQn3,0)

    call BOUNDARY(0    ,LQb2,0)
    ! No radial beam flux at rho = rhob.  With MDBEAM = 1, LQb2 is the Lagrange
    ! multiplier of the beam radial force balance and nothing ties its wall
    ! value to n_b; the fast-ion orbit width exceeds the SOL gradient length,
    ! so that value is no outflow law.  The beam has no surface loss channel
    ! in this model (outside rho = 1 it is lost by the volume loss nu_bL).
    ! To be reconsidered with a finite-orbit-width or orbit-loss model.
    if(MDBEAM == 1) call BOUNDARY(NRMAX,LQb2,0)

    ! Neumann condition of the continuity equation at the boundary is naturally
    ! imposed through the diffusion term in the pressure equation. However, this
    ! condition is very weak because it comes from nondiagonal term. Then by setting
    ! the density constant in the last element, we can force the density gradient
    ! to be nought at the boundary explicitly.

!    call BOUNDARY(NRMAX,LQe1,0,Var(NRMAX-1,1)%n)
!    call BOUNDARY(NRMAX,LQi1,0,Var(NRMAX-1,2)%n)
!    call BOUNDARY(NRMAX,LQz1,0,Var(NRMAX-1,3)%n)

    ! MDLWBC = 1: the zero density gradient of the table in tx560.tex, imposed
    ! as a constraint n(NRMAX) = n(NRMAX-1) on each charged density in
    ! constrained-DOF form (see BOUNDARY_ZEROGRAD).  With MDLWBC = 0 (V5.52)
    ! the density has no condition at the wall, and the radial particle flux
    ! through rho = rhob then depends on the last element.
    if( MDLWBC == 1 ) then
       call BOUNDARY_ZEROGRAD(LQe1)
       call BOUNDARY_ZEROGRAD(LQi1)
       call BOUNDARY_ZEROGRAD(LQz1)
    end if
!!    call BOUNDARY(0,LQe3,0,Var(1,1)%BUpar)
!!    call BOUNDARY(0,LQi3,0,Var(1,2)%BUpar)
!!    call BOUNDARY(0,LQz3,0,Var(1,3)%BUpar)

    !     Integral term stemming from integration by parts in the diffusion term

    !  NOTE: ckt(NRMAX)*fipol(NRMAX)/(fourPisq*rMU0) is equivalent to
    !        fipol(NRMAX)*rIp*1.D6/(2.d0*Pi*sdt(NRMAX)).
    !        See "sdtvac" in TXCALV.
    call BOUNDARY(NRMAX,LQm2,1, ckt(NRMAX)*fipol(NRMAX)/(fourPisq*rMU0))
    call BOUNDARY(NRMAX,LQm3,1, 2.d0*Pi*rMUb1*rIp*1.D6)
    call BOUNDARY(NRMAX,LQn1,1, suft(NRMAX)*rGASPF)
    call BOUNDARY(NRMAX,LQnz,1, suft(NRMAX)*rGASPFz)

    if(ICALA == 0) ICALA = 1

    call dump_axis_LQs8
    call dump_axis_row(LQz1,'LQz1    ')
    call dump_axis_row(LQi1,'LQi1    ')
    call dump_axis_coefs

  end subroutine TXCALA

  subroutine ensure_txcala_workspace(NSMP)

    integer(4), intent(in) :: NSMP
    logical :: need_elm, need_dtf, need_nr, need_nrns, need_ratden, need_zchi

    need_elm = .not. allocated(ELM)
    if(.not. need_elm) then
       need_elm = (lbound(ELM,1) /= 1) .or. (ubound(ELM,1) /= NEMAX) .or. &
            &     (lbound(ELM,2) /= 1) .or. (ubound(ELM,2) /= 4)     .or. &
            &     (lbound(ELM,3) /= 0) .or. (ubound(ELM,3) /= NCM)   .or. &
            &     (lbound(ELM,4) /= 1) .or. (ubound(ELM,4) /= NQMAX)
    end if
    if(need_elm) then
      if(allocated(ELM)) deallocate(ELM)
      allocate(ELM(1:NEMAX,1:4,0:NCM,1:NQMAX))
    end if

    need_dtf = .not. allocated(DTf)
    if(.not. need_dtf) then
       need_dtf = (lbound(DTf,1) /= 1) .or. (ubound(DTf,1) /= NQMAX)
    end if
    if(need_dtf) then
      if(allocated(DTf)) deallocate(DTf)
      allocate(DTf(1:NQMAX))
    end if

    need_nr = .not. allocated(rrtinv)
    if(.not. need_nr) then
       need_nr = (lbound(rrtinv,1) /= 0) .or. (ubound(rrtinv,1) /= NRMAX)
    end if
    if(need_nr) then
      if(allocated(rrtinv)) then
        deallocate(rrtinv, aatinv, bthcoinv, ParatoZeta, ZetatoPara, &
             &     aatq, PsitdotVBth, fipolinv, fipolsdt, dZetatoPara, ribi, ribsdt, &
             &     RatCXSNBe, RatCXSNBi, Ratdenhyd)
      end if
      allocate(rrtinv(0:NRMAX), aatinv(0:NRMAX), bthcoinv(0:NRMAX), ParatoZeta(0:NRMAX), ZetatoPara(0:NRMAX), &
           &   aatq(0:NRMAX), PsitdotVBth(0:NRMAX), fipolinv(0:NRMAX), fipolsdt(0:NRMAX), dZetatoPara(0:NRMAX), &
           &   ribi(0:NRMAX), ribsdt(0:NRMAX), RatCXSNBe(0:NRMAX), RatCXSNBi(0:NRMAX), Ratdenhyd(0:NRMAX))
    end if

    need_nrns = .not. allocated(UsrgV)
    if(.not. need_nrns) then
       need_nrns = (lbound(UsrgV,1) /= 0) .or. (ubound(UsrgV,1) /= NRMAX) .or. &
            &      (lbound(UsrgV,2) /= 1) .or. (ubound(UsrgV,2) /= NSM)
    end if
    if(need_nrns) then
      if(allocated(UsrgV)) then
        deallocate(UsrgV, fipolsdtPNsV, BUsparVbbt, dlnPNsV, BpBVsdiag, tormflux, dtormflux, ribsdtNs)
      end if
      allocate(UsrgV(0:NRMAX,1:NSM), fipolsdtPNsV(0:NRMAX,1:NSM), BUsparVbbt(0:NRMAX,1:NSM), &
           &   dlnPNsV(0:NRMAX,1:NSM), BpBVsdiag(0:NRMAX,1:NSM), tormflux(0:NRMAX,1:NSM), &
           &   dtormflux(0:NRMAX,1:NSM), ribsdtNs(0:NRMAX,1:NSM))
    end if

    need_ratden = .not. allocated(Ratden)
    if(.not. need_ratden) then
       need_ratden = (lbound(Ratden,1) /= 0) .or. (ubound(Ratden,1) /= NRMAX) .or. &
            &        (lbound(Ratden,2) /= 1) .or. (ubound(Ratden,2) /= NSMP) .or. &
            &        (lbound(Ratden,3) /= 1) .or. (ubound(Ratden,3) /= NSMP)
    end if
    if(need_ratden) then
      if(allocated(Ratden)) deallocate(Ratden)
      allocate(Ratden(0:NRMAX,1:NSMP,1:NSMP))
    end if

    need_zchi = .not. allocated(zChi)
    if(.not. need_zchi) then
       need_zchi = (lbound(zChi,1) /= 0) .or. (ubound(zChi,1) /= NRMAX) .or. &
            &      (lbound(zChi,2) /= 1) .or. (ubound(zChi,2) /= NSM)   .or. &
            &      (lbound(zChi,3) /= 1) .or. (ubound(zChi,3) /= NSM)
    end if
    if(need_zchi) then
      if(allocated(zChi)) deallocate(zChi)
      allocate(zChi(0:NRMAX,1:NSM,1:NSM))
    end if

  end subroutine ensure_txcala_workspace

!***************************************************************
!
!   Coefficients for Equations
!
!**************************************************************

  subroutine LQcoef

    use tx_interface, only : dfdx
    integer(4) :: i


    do i = 1, NSM
       zChi(:,i,1) = Chis(:,i) + ChiNCt(:,i) + ChiNCp(:,i)
       zChi(:,i,2) = Chis(:,i) + ChiNCt(:,i)
    end do

    Ratden(:,1,2) = Var(:,1)%n / Var(:,2)%n ! Ne/Ni
    Ratden(:,1,3) = Var(:,1)%n / Var(:,3)%n ! Ne/Nz
    Ratden(:,2,1) = Var(:,2)%n / Var(:,1)%n ! Ni/Ne
    Ratden(:,2,3) = Var(:,2)%n / Var(:,3)%n ! Ni/Nz
    Ratden(:,3,1) = Var(:,3)%n / Var(:,1)%n ! Nz/Ne
    Ratden(:,3,2) = Var(:,3)%n / Var(:,2)%n ! Nz/Ni
    do i = 1, NSM
       Ratden(:,i,4) = Var(:,i)%n * PNbVinv(:) ! Ns/Nb
    end do
    Ratdenhyd(:)   = (Var(:,2)%n + PNbV(:)) / Var(:,1)%n ! (Ni+Nb)/Ne

    fipolinv(:)    = 1.d0 / fipol(:)
    fipolsdt(:)    = fipol(:) / sdt(:)
    bthcoinv(:)    = 1.d0 / bthco(:)
    PsitdotVBth(:) = fourPisq * sdt(:) * PsitdotV(:)
    aatq(:)        = aat(:) / Q(:)

    rrtinv(:)      = 1.d0 / rrt(:) ! 1/<R^2>
    aatinv(:)      = 1.d0 / aat(:) ! 1/<R^-2>
    ParatoZeta(:)  = fipol(:) / bbt(:) ! I/<B^2> : parallel to toroidal
    ZetatoPara(:)  = bbt(:) / fipol(:) ! <B^2>/I : toroidal to parallel
    dZetatoPara(:) = dfdx(vv,ZetatoPara,NRMAX,0) ! d(<B^2>/I)/dV
    RatCXSNBe(:)   = (1.d0 - RatCX(:)) * SNBe(:)
    RatCXSNBi(:)   = RatCX(:) * SNBi(:)
    ribi(:)        = bri(:) / (fipol(:) * bbt(:)) ! (<R^2>-I^2/<B^2>)/I
    ribsdt(:)      = ribi(:) * fipolsdt(:) ! (<B^2><R^2>-I^2)/<B^2>/(dpsi/dV)

    do i = 1, NSM
       UsrgV(:,i)        = FSADV * Var(:,i)%UrV - UgV(:)
       fipolsdtPNsV(:,i) = fipolsdt(:) / Var(:,i)%n
       BUsparVbbt(:,i)   = Var(:,i)%BUpar / bbt(:)
       BpBVsdiag(:,i)    = ribi(:) * BVsdiag(:,i,1)
       tormflux(:,i)     = suft(:) * Vmps(:,i) * Var(:,i)%n * BpBVsdiag(:,i) &
            &            - sst(:)  * rMus(:,i) * Var(:,i)%n * dfdx(vv,BpBVsdiag(:,i),NRMAX,0)
       dtormflux(:,i)    = dfdx(vv,tormflux(:,i),NRMAX,0)
       ribsdtNs(:,i)     = ribsdt(:) * Var(:,i)%n
!       dlnPNsV(:,i)      = dfdx(vv,Var(:,i)%n,NRMAX,0) / Var(:,i)%n
!       RatPNbinv(:,i)    = PNbV(:) / Var(:,i)%n ! Nb/Ns
    end do

!    rGASPFA(0:NRMAX-1) = 0.d0
!    rGASPFA(NRMAX) = rGASPF

!!rp_conv    rNubLL(:) = rNubL(:) * rip_rat(:)

  end subroutine LQcoef

!***************************************************************
!
!   Poisson Equation: phi
!
!**************************************************************

  subroutine LQm1CC

    integer(4) :: NEQ = LQm1
    real(8) :: sqeps0inv = 1.d0 / sqeps0

    ! phi'(0) : 0

    ELM(:,:,1,NEQ) =   sqeps0 / AEE * 1.D-20 * fem_int(12,sst)
    NLC(1,NEQ) = LQm1

    ELM(:,:,2,NEQ) = - achg(1) * fem_int(1) * sqeps0inv
    NLC(2,NEQ) = LQe1

    ELM(:,:,3,NEQ) = - achg(2) * fem_int(1) * sqeps0inv
    NLC(3,NEQ) = LQi1

    ELM(:,:,4,NEQ) = - achg(3) * fem_int(1) * sqeps0inv
    NLC(4,NEQ) = LQz1

    ELM(:,:,5,NEQ) = - achgb   * fem_int(1) * sqeps0inv * BeamSW
    NLC(5,NEQ) = LQb1

    ELM(:,:,6,NEQ) = - achgb   * fem_int(2,rip_rat) * sqeps0inv * RpplSW
    NLC(6,NEQ) = LQr1

    ! phi(b) : 0

    NLCMAX(NEQ) = 6

  end subroutine LQm1CC

!***************************************************************
!
!   Ampere's Law: psit'
!
!***************************************************************

  subroutine LQm2CC

    integer(4) :: NEQ = LQm2

    ! psit'(0) : 0

    ELM(:,:,0,NEQ) =   fem_int(1) * EPS0 * invDT
    NLC(0,NEQ) = LQm2

    ! rot Bphi

    ELM(:,:,1,NEQ) = - fem_int(27,aatinv,ckt) - fem_int(25,aatinv,ckt)
    NLC(1,NEQ) = LQm5

    ! Electron current

    ELM(:,:,2,NEQ) =   achg(1) * AEE * 1.D20 * fem_int(20,bthcoinv,Var(:,1)%n)
    NLC(2,NEQ) = LQe3

    ELM(:,:,3,NEQ) = - achg(1) * AEE * 1.D20 * fem_int(20,bthcoinv,fipol)
    NLC(3,NEQ) = LQe7

    ! Ion current

    ELM(:,:,4,NEQ) =   achg(2) * AEE * 1.D20 * fem_int(20,bthcoinv,Var(:,2)%n)
    NLC(4,NEQ) = LQi3
 
    ELM(:,:,5,NEQ) = - achg(2) * AEE * 1.D20 * fem_int(20,bthcoinv,fipol)
    NLC(5,NEQ) = LQi7

    ! Impurity current

    ELM(:,:,6,NEQ) =   achg(3) * AEE * 1.D20 * fem_int(20,bthcoinv,Var(:,3)%n)
    NLC(6,NEQ) = LQz3
 
    ELM(:,:,7,NEQ) = - achg(3) * AEE * 1.D20 * fem_int(20,bthcoinv,fipol)
    NLC(7,NEQ) = LQz7

    if( MDBEAM /= 2 ) then

       ! Beam ion current

       ELM(:,:,8,NEQ) =   achgb   * AEE * 1.D20 * fem_int(2,bthcoinv) * BeamSW
       NLC(8,NEQ) = LQb3

       ELM(:,:,9,NEQ) = - achgb   * AEE * 1.D20 * fem_int(20,bthcoinv,fipol) * BeamSW
       NLC(9,NEQ) = LQb7

    end if

    ! (NRMAX) : ckt(NRMAX)*Ivac/(4 pi^2)

    NLCMAX(NEQ) = 9

  end subroutine LQm2CC

!***************************************************************
!
!   Ampere's Law: psi'
!
!***************************************************************

  subroutine LQm3CC

    integer(4) :: NEQ = LQm3

    ! psi'(0) : 0

    ELM(:,:,0,NEQ) =   fem_int(2,aat) * EPS0 * rMUb1 * invDT
    NLC(0,NEQ) = LQm3

    ! rot Btheta

    ELM(:,:,1,NEQ) = - fem_int(12,ckt)
    NLC(1,NEQ) = LQm4

    ! Electron current

    ELM(:,:,2,NEQ) = - rMUb1 * achg(1) * AEE * 1.D20 * fem_int(1)
    NLC(2,NEQ) = LQe7

    ! Ion current

    ELM(:,:,3,NEQ) = - rMUb1 * achg(2) * AEE * 1.D20 * fem_int(1)
    NLC(3,NEQ) = LQi7

    ! Impurity current

    ELM(:,:,4,NEQ) = - rMUb1 * achg(3) * AEE * 1.D20 * fem_int(1)
    NLC(4,NEQ) = LQz7

    if( MDBEAM /= 2 ) then

       ! Beam ion current

       ELM(:,:,5,NEQ) = - rMUb1 * achgb   * AEE * 1.D20 * fem_int(1) * BeamSW
       NLC(5,NEQ) = LQb7

    end if

!    ! Virtual current for helical system
!
!    ELM(:,:,6,NEQ) = - rMUb1 * fem_int(-1,AJV) / rr
!    NLC(6,NEQ) = 0

    ! psi'(NRMAX) : Ip

    NLCMAX(NEQ) = 5

  end subroutine LQm3CC

!**************************************************************
!
!   Faraday's Law : psi
!
!***************************************************************

  subroutine LQm4CC

    integer(4) :: NEQ = LQm4

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQm4

    ! psi'

    ELM(:,:,1,NEQ) = fem_int(1) / rMUb2
    NLC(1,NEQ) = LQm3

    NLCMAX(NEQ) = 1

  end subroutine LQm4CC

!***************************************************************
!
!   Faraday's Law : psit
!
!***************************************************************

  subroutine LQm5CC

    integer(4) :: NEQ = LQm5

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQm5

    ! psit'

    ELM(:,:,1,NEQ) = fem_int(1) / rMU0
    NLC(1,NEQ) = LQm2

    NLCMAX(NEQ) = 1

  end subroutine LQm5CC

!***************************************************************
!
!   Electron Density Equation
!
!***************************************************************

  subroutine LQe1CC

    integer(4) :: NEQ = LQe1, i = 1
    real(8), dimension(1:NEMAX,1:4) :: int_tmp

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQe1

    ! Divergence

    ELM(:,:,1,NEQ) = - fem_int(4)
    NLC(1,NEQ) = LQe2

    ! Divergence + Grid velocity

    ELM(:,:,2,NEQ) =   fem_int(5,UgV)
    NLC(2,NEQ) = LQe1

    ! Ionization of n01, n02 and n03 by electron impact

    int_tmp = fem_int(20,SiVizA,Var(:,i)%n)
    ELM(:,:,3,NEQ) =   FSION * 1.D20 * int_tmp
    NLC(3,NEQ) = LQn1

    ELM(:,:,4,NEQ) =   FSION * 1.D20 * int_tmp * ThntSW
    NLC(4,NEQ) = LQn2

    ELM(:,:,5,NEQ) =   FSION * 1.D20 * int_tmp * BeamSW
    NLC(5,NEQ) = LQn3

    ! Increment of electrons due to effective ionization of impurities

    ELM(:,:,6,NEQ) =   FSION * 1.D20 * fem_int(20,SiVsefA,Var(:,1)%n) * achg(3) 
    NLC(6,NEQ) = LQnz

    ! Loss of electrons due to effective recombination of impurities

    ELM(:,:,7,NEQ) = - FSCX  * 1.D20 * fem_int(20,SiVa6A, Var(:,1)%n) * achg(3)
    NLC(7,NEQ) = LQz1

    ! Loss to divertor plasma region

    ELM(:,:,8,NEQ) = -             fem_int( 2,rNuL)
    NLC(8,NEQ) = LQe1

    ELM(:,:,9,NEQ) =   PNsSUP(1) * fem_int(-1,rNuL)
    NLC(9,NEQ) = 0

    ! Generated by NBI (Ionization)

!    ELM(:,:,10,NEQ) =   (1.d0 - RatCX) * fem_int(-1,SNBe)
    ELM(:,:,10,NEQ) =   fem_int(-1,RatCXSNBe)
    NLC(10,NEQ) = 0

!    !  Diffusion of electrons (***AF 2008-06-08)
!
!    ELM(:,:,11,NEQ) = - fem_int(27,sst,DMAGe)
!    NLC(11,NEQ) = LQe1

    NLCMAX(NEQ) = 10

  end subroutine LQe1CC

!***************************************************************
!
!   Electron Radial Flow
!
!***************************************************************

  subroutine LQe2CC

    integer(4) :: NEQ = LQe2, i = 1
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(i) * amp)

    ! Diamagnetic force

    ELM(:,:,1,NEQ) =   achg(i) * aee * amasinv * fem_int(30,bri,Var(:,i)%n,fipolinv)
    NLC(1,NEQ) = LQe8

    ! v x B force

    ELM(:,:,2,NEQ) =   achg(i) * aee * amasinv * fem_int(20,fipol,Var(:,i)%n)
    NLC(2,NEQ) = LQe3

    ELM(:,:,3,NEQ) = - achg(i) * aee * amasinv * fem_int( 2,bbt)
    NLC(3,NEQ) = LQe4

    NLCMAX(NEQ) = 3

  end subroutine LQe2CC

!***************************************************************
!
!   Electron Parallel Flow
!
!***************************************************************

  subroutine LQe3CC

    integer(4) :: NEQ = LQe3, i = 1
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(i) * amp)

    ELM(:,:, 0,NEQ) =   fem_int(2+L3,Var(:,i)%n) * invDT
    NLC( 0,NEQ) = LQe3

    ! Neoclassical viscosity force

    ELM(:,:, 1,NEQ) = - amasinv * 1.D-20 * fem_int(2+L3,xmu(:,i,1,1))
    NLC( 1,NEQ) = LQe3

    ELM(:,:, 2,NEQ) = - amasinv * 1.D-20 * fem_int(2+L3,xmu(:,i,1,2))
    NLC( 2,NEQ) = LQe6

    ! Diamagnetic forces (particle)

    ELM(:,:, 3,NEQ) =   amasinv * 1.D-20 * fem_int(2+L3,xmu(:,i,1,1))
    NLC( 3,NEQ) = LQe8

    ! Diamagnetic forces (heat)

    ELM(:,:, 4,NEQ) = - rKilo / achg(i) * amasinv * 1.D-20 * fem_int(22+L3,fipolsdtPNsV(:,i),xmu(:,i,1,2))
    NLC( 4,NEQ) = LQe5

    ELM(:,:, 5,NEQ) =   rKilo / achg(i) * amasinv * 1.D-20 * fem_int(32+L3,fipolsdtPNsV(:,i),Var(:,i)%T,xmu(:,i,1,2))
    NLC( 5,NEQ) = LQe1
!    ELM(:,:, 5,NEQ) =   rKilo / achg(i) * amasinv * 1.D-20 * fem_int(30+L3,fipolsdtPNsV(:,i),dlnPNsV(:,i),xmu(:,i,1,2))
!    NLC( 5,NEQ) = LQe5

    ! Collisional friction force with electrons (self collision)

    ELM(:,:, 6,NEQ) =   fem_int(20+L3,lab(:,i,1,1,1),Var(:,i)%n)
    NLC( 6,NEQ) = LQe3

    ! Collisional friction force with ions

    ELM(:,:, 7,NEQ) =   fem_int(20+L3,lab(:,i,2,1,1),Var(:,i)%n)
    NLC( 7,NEQ) = LQi3

    ! Collisional friction force with impurities

    ELM(:,:, 8,NEQ) =   fem_int(20+L3,lab(:,i,3,1,1),Var(:,i)%n)
    NLC( 8,NEQ) = LQz3

    ! Collisional friction force with beam ions

    ELM(:,:, 9,NEQ) =   fem_int(20+L3,laf(:,i,1,1),  Ratden(:,i,4))
    NLC( 9,NEQ) = LQb3

    ! Collisional heat friction force with electrons (self collision)

    ELM(:,:,10,NEQ) = - fem_int(20+L3,lab(:,i,1,1,2),Var(:,i)%n)
    NLC(10,NEQ) = LQe6

    ! Collisional heat friction force with ions

    ELM(:,:,11,NEQ) = - fem_int(20+L3,lab(:,i,2,1,2),Var(:,i)%n)
    NLC(11,NEQ) = LQi6

    ! Collisional heat friction force with impurities

    ELM(:,:,12,NEQ) = - fem_int(20+L3,lab(:,i,3,1,2),Var(:,i)%n)
    NLC(12,NEQ) = LQz6

    ! Electric field forces

    ELM(:,:,13,NEQ) = - achg(i) * aee * amasinv * fem_int(30+L3,fipol,Var(:,i)%n,aatq)
    NLC(13,NEQ) = LQm2

    ELM(:,:,14,NEQ) =   achg(i) * aee * amasinv * fem_int(30+L3,fipol,Var(:,i)%n,aat)
    NLC(14,NEQ) = LQm3

    ! Loss to divertor plasma region

    ELM(:,:,15,NEQ) = - fem_int(20+L3,rNuL,Var(:,i)%n) * fact
    NLC(15,NEQ) = LQe3

    ! Collisional friction force with neutrals

    ELM(:,:,16,NEQ) = - fem_int(20+L3,rNu0s(:,i),Var(:,i)%n)
    NLC(16,NEQ) = LQe3

    ! ---

    ! Perpendicular viscosity force: diffusion

    ELM(:,:,17,NEQ) = - FSPARV(i) * fem_int(38,sst*rMus(:,i)*Var(:,i)%n, ZetatoPara,ParatoZeta)
    NLC(17,NEQ) = LQe3

    ELM(:,:,18,NEQ) = - FSPARV(i) * fem_int(34,sst*rMus(:,i)*Var(:,i)%n,dZetatoPara,ParatoZeta)
    NLC(18,NEQ) = LQe3

    ! Perpendicular viscosity force: pinch

    ELM(:,:,19,NEQ) = - FSPARV(i) * fem_int(24,suft*Vmps(:,i),Var(:,i)%n)
    NLC(19,NEQ) = LQe3

    ELM(:,:,20,NEQ) =   FSPARV(i) * fem_int(30,suft*Vmps(:,i)*Var(:,i)%n,dZetatoPara,ParatoZeta)
    NLC(20,NEQ) = LQe3

    ! Perpendicular viscosity force: residual stress

    ELM(:,:,21,NEQ) = - amasinv * fem_int(-5,ZetatoPara,PiRess(:,i))
    NLC(21,NEQ) = 0

    ! Perpendicular viscosity force: terms associated with diamagnetic flow

    ELM(:,:,22,NEQ) = - FSPARV(i) * fem_int(-2,ZetatoPara,dtormflux(:,i))
    NLC(22,NEQ) = 0

!!$    ! Alternative implementation of perpendicular viscosity force: terms associated with diamagnetic flow
!!$
!!$    ELM(:,:,22,NEQ) = - FSPARV(i) * fem_int(38,sst *rMus(:,i)*Var(:,i)%n, ZetatoPara,ribi)
!!$    NLC(22,NEQ) = LQe8
!!$
!!$    ELM(:,:,23,NEQ) = - FSPARV(i) * fem_int(34,sst *rMus(:,i)*Var(:,i)%n,dZetatoPara,ribi)
!!$    NLC(23,NEQ) = LQe8
!!$
!!$    ELM(:,:,24,NEQ) = - FSPARV(i) * fem_int(24,suft*Vmps(:,i)*Var(:,i)%n, ZetatoPara*ribi)
!!$    NLC(24,NEQ) = LQe8
!!$
!!$    ELM(:,:,25,NEQ) =   FSPARV(i) * fem_int(30,suft*Vmps(:,i)*Var(:,i)%n,dZetatoPara,ribi)
!!$    NLC(25,NEQ) = LQe8

    ! ---

!!$    ! Perpendicular viscosity (off-diagonal)
!!$
!!$    ELM(:,:,17,NEQ) = - FSPARV(i) * fem_int(37,ZetatoPara,sst,rMus(:,i))
!!$    NLC(17,NEQ) = LQe4
!!$
!!$    ELM(:,:,18,NEQ) =   FSPARV(i) * fem_int(37,ZetatoPara,sst*rMus(:,i),Var(:,i)%RUph)
!!$    NLC(18,NEQ) = LQe1
!!$  
!!$    ELM(:,:,19,NEQ) = - FSPARV(i) * fem_int(32,dZetatoPara,sst,rMus(:,i))
!!$    NLC(19,NEQ) = LQe4
!!$
!!$    ELM(:,:,20,NEQ) =   FSPARV(i) * fem_int(32,dZetatoPara,sst*rMus(:,i),Var(:,i)%RUph)
!!$    NLC(20,NEQ) = LQe1
!!$
!!$    ! Momentum pinch (off-diagonal)
!!$
!!$    ELM(:,:,21,NEQ) = - FSPARV(i) * fem_int(24,suft*Vmps(:,i),ZetatoPara)
!!$    NLC(21,NEQ) = LQe4
!!$
!!$    ELM(:,:,22,NEQ) =   FSPARV(i) * fem_int(30,suft,Vmps(:,i),dZetatoPara)
!!$    NLC(22,NEQ) = LQe4

!    ! Diffusion of electrons (***AF 2008-06-08) (37+L3 is not possible because it is a diffusion term)
!
!    ELM(:,:,23,NEQ) = - fem_int(37,sst,DMAGe,Var(:,i)%n)
!    NLC(23,NEQ) = LQe3

!    ! Additional torque input
!
!    ELM(:,:,24,NEQ) = - amasinv * 1.D-20 * fem_int(-(2+L3),Tqt,ZetatoPara)
!    NLC(24,NEQ) = 0

    NLCMAX(NEQ) = 22

    if( MDFIXT /= 0 ) then ! Overwrite Pressure gradient force
       ELM(:,:, 4,NEQ) = - rKilo / achg(i) * amasinv * 1.D-20 * fem_int(-(22+L3),fipolsdt,xmu(:,i,1,2),Var(:,i)%T)
       NLC( 4,NEQ) = 0

       ELM(:,:, 5,NEQ) =   0.d0
       NLC( 5,NEQ) = 0
    end if

  end subroutine LQe3CC

!***************************************************************
!
!   Electron Toroidal Flow
!
!***************************************************************

  subroutine LQe4CC

    integer(4) :: NEQ = LQe4, i = 1, j
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(i) * amp)

    ! Uephi(0)' : 0

    ELM(:,:, 0,NEQ) = fem_int(1) * invDT
    NLC( 0,NEQ) = LQe4

    ! Advection + Grid velocity

    ELM(:,:, 1,NEQ) = - fem_int(5,UsrgV(:,i))
    NLC( 1,NEQ) = LQe4

    ELM(:,:, 2,NEQ) = - FSADV * fem_int(6,Var(:,i)%UrV)
    NLC( 2,NEQ) = LQe4

    ! Momentum pinch

    ELM(:,:, 3,NEQ) = - fem_int(24,suft,Vmps(:,i))
    NLC( 3,NEQ) = LQe4

    ! Viscosity force

    ELM(:,:, 4,NEQ) = - fem_int(27,sst,rMus(:,i))
    NLC( 4,NEQ) = LQe4

    ELM(:,:, 5,NEQ) =   fem_int(37,sst,rMus(:,i),Var(:,i)%RUph)
    NLC( 5,NEQ) = LQe1
!!$    ELM(:,:, 5,NEQ) = - fem_int(24,sst*rMus(:,i),dlnPNsV(:,i))
!!$    NLC( 5,NEQ) = LQe4

    ! Residual stress

    ELM(:,:, 6,NEQ) = - amasinv * fem_int(-4,PiRess(:,i))
    NLC( 6,NEQ) = 0

    ! Collisional friction force with electrons (self collision)

    ELM(:,:, 7,NEQ) =   fem_int( 2,lab(:,i,1,1,1))
    NLC( 7,NEQ) = LQe4

    ! Collisional friction force with bulk ions

    ELM(:,:, 8,NEQ) =   fem_int(20,lab(:,i,2,1,1),Ratden(:,i,2))
    NLC( 8,NEQ) = LQi4

    ! Collisional friction force with impurities

    ELM(:,:, 9,NEQ) =   fem_int(20,lab(:,i,3,1,1),Ratden(:,i,3))
    NLC( 9,NEQ) = LQz4

    ! Collisional friction force with beam ions

!    ELM(:,:,10,NEQ) =   fem_int(30,ParatoZeta,laf(:,i,1,1),Ratden(:,i,4))
!    NLC(10,NEQ) = LQb3
    ELM(:,:,10,NEQ) =   fem_int(20,laf(:,i,1,1),Ratden(:,i,4))
    NLC(10,NEQ) = LQb4

    ! Collisional heat friction force with electrons (self collision)

    j = 1
    ELM(:,:,11,NEQ) =   rKilo / achg(j) * fem_int(-22,lab(:,i,j,1,2),ribsdtNs(:,i),Var(:,j)%T)
    NLC(11,NEQ) = 0
!!$    ELM(:,:,11,NEQ) =   rKilo / achg(j) * fem_int(31,lab(:,i,j,1,2),ribsdt,Var(:,j)%T)
!!$    NLC(11,NEQ) = LQe1

    ELM(:,:,12,NEQ) = -                   fem_int(30,lab(:,i,j,1,2),ParatoZeta,Var(:,i)%n)
    NLC(12,NEQ) = LQe6

    ! Collisional heat friction force with ions

    j = 2
    ELM(:,:,13,NEQ) =   rKilo / achg(j) * fem_int(-22,lab(:,i,j,1,2),ribsdtNs(:,i),Var(:,j)%T)
    NLC(13,NEQ) = 0
!!$    ELM(:,:,13,NEQ) =   rKilo / achg(j) * fem_int(31,lab(:,i,j,1,2),ribsdt,Var(:,j)%T)
!!$    NLC(13,NEQ) = LQe1

    ELM(:,:,14,NEQ) = -                   fem_int(30,lab(:,i,j,1,2),ParatoZeta,Var(:,i)%n)
    NLC(14,NEQ) = LQi6

    ! Collisional heat friction force with impurities

    j = 3
    ELM(:,:,15,NEQ) =   rKilo / achg(j) * fem_int(-22,lab(:,i,j,1,2),ribsdtNs(:,i),Var(:,j)%T)
    NLC(15,NEQ) = 0
!!$    ELM(:,:,15,NEQ) =   rKilo / achg(j) * fem_int(31,lab(:,i,j,1,2),ribsdt,Var(:,j)%T)
!!$    NLC(15,NEQ) = LQe1

    ELM(:,:,16,NEQ) = -                   fem_int(30,lab(:,i,j,1,2),ParatoZeta,Var(:,i)%n)
    NLC(16,NEQ) = LQz6

    ! Toroidal E force

    ELM(:,:,17,NEQ) =   achg(i) * aee * amasinv * fem_int(2,Var(:,i)%n)
    NLC(17,NEQ) = LQm3

    ! v x B force

    ELM(:,:,18,NEQ) =   achg(i) * aee * amasinv * fem_int(6,PsiV)
    NLC(18,NEQ) = LQe2

    ! Turbulent particle transport driver (electron driven)

    if( FSTPTM(i) == 0 ) then
       ELM(:,:,19,NEQ) =   achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Var(:,1)%T) &
            &                           * (1.d0 - FSVAHL)
       NLC(19,NEQ) = LQe1

       ELM(:,:,20,NEQ) =   achg(i) * rKeV * amasinv * fem_int( 5,FQLcoef) &
            &                           * FSVAHL
       NLC(20,NEQ) = LQe5

    else if( FSTPTM(i) == 1 ) then
       ELM(:,:,19,NEQ) =   achg(i) * aee * amasinv * fem_int(20,FQLcoef2,ZetatoPara)
       NLC(19,NEQ) = LQe4

       ELM(:,:,20,NEQ) = - achg(i) * aee * amasinv * fem_int(20,FQLcoef2,Var(:,i)%n)
       NLC(20,NEQ) = LQe3

       ELM(:,:,21,NEQ) =   achg(i) * rKeV * amasinv * fem_int( 5,FQLcoef) &
            &                           * (FSVAHL - 1.d0)
       NLC(21,NEQ) = LQe5

       ELM(:,:,22,NEQ) = - achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Var(:,i)%T) &
            &                           * (FSVAHL - 1.d0)
       NLC(22,NEQ) = LQe1

       ELM(:,:,23,NEQ) =   achg(i) * aee * amasinv * fem_int(22,FQLcoef,Var(:,i)%n)
       NLC(23,NEQ) = LQm1

    else
       ELM(:,:,19,NEQ) =   achg(i) * aee * amasinv * fem_int(20,FQLcoef1,Var(:,i)%n)
       NLC(19,NEQ) = LQe8

       ELM(:,:,20,NEQ) =   achg(i) * aee * amasinv * fem_int(22,FQLcoef,Var(:,i)%n)
       NLC(20,NEQ) = LQm1

       ELM(:,:,21,NEQ) =   achg(i) * rKeV * amasinv * fem_int( 5,FQLcoef) &
            &                           * (FSVAHL - 1.d0)
       NLC(21,NEQ) = LQe5

       ELM(:,:,22,NEQ) = - achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Var(:,i)%T) &
            &                           * (FSVAHL - 1.d0)
       NLC(22,NEQ) = LQe1

    end if

    !  -- Turbulent pinch term

    ELM(:,:,24,NEQ) =   achg(i) * aee * amasinv * fem_int(20,FQLcoef,FVpch)
    NLC(24,NEQ) = LQe1

    ! Loss to divertor plasma region

    ELM(:,:,25,NEQ) = - fem_int(2,rNuL) * fact
    NLC(25,NEQ) = LQe4

    ! Collisional friction force with neutrals

    ELM(:,:,26,NEQ) = - fem_int(2,rNu0s(:,i))
    NLC(26,NEQ) = LQe4

!    ! Diffusion of electrons (***AF 2008-06-08)
!
!    ELM(:,:,27,NEQ) = - fem_int(27,sst,DMAGe)
!    NLC(27,NEQ) = LQe4

!    ! Additional torque input
!
!    ELM(:,:,28,NEQ) = - amasinv * 1.D-20 * fem_int(-1,Tqt) 
!    NLC(28,NEQ) = 0

    NLCMAX(NEQ) = 26

    if( MDFIXT /= 0 ) then ! Overwrite Pressure gradient force
       if( FSTPTM(i) == 0 ) then
          ELM(:,:,19,NEQ) =   achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Var(:,1)%T)
          NLC(19,NEQ) = LQe1

          ELM(:,:,20,NEQ) =   achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Var(:,i)%n) &
               &                           * FSVAHL
          NLC(20,NEQ) = LQe5

       else
          ELM(:,:,21,NEQ) =   achg(i) * rKeV * amasinv * fem_int(21,FQLcoef,Var(:,i)%T) &
               &                           * (FSVAHL - 1.d0)
          NLC(21,NEQ) = LQe1

          ELM(:,:,22,NEQ) = 0.d0
          NLC(22,NEQ) = 0

       end if
    end if

  end subroutine LQe4CC

!***************************************************************
!
!   Electron Energy Transport: Te
!
!***************************************************************

  subroutine LQe5CC

    integer(4) :: NEQ = LQe5, i = 1
    real(8), dimension(1:NEMAX,1:4) :: int_tmp

    ! Temperature evolution
    
    if( MDFIXT == 0 ) then
       ELM(:,:, 0,NEQ) =   1.5d0 * fem_int(1) * invDT
       NLC( 0,NEQ) = LQe5

       ! Advection + Grid velocity

       ELM(:,:, 1,NEQ) = - 2.5d0 * fem_int(5,UsrgV(:,i))
       NLC( 1,NEQ) = LQe5

       ELM(:,:, 2,NEQ) = - 2.5d0 * FSADV * fem_int(6,Var(:,i)%UrV)
       NLC( 2,NEQ) = LQe5

       ! Heat pinch

       ELM(:,:, 3,NEQ) = - fem_int(24,Vhps(:,i),suft)
       NLC( 3,NEQ) = LQe5

       ! Conduction transport

       ELM(:,:, 4,NEQ) = - fem_int(27,sst,zChi(:,i,1))
       NLC( 4,NEQ) = LQe5

       ELM(:,:, 5,NEQ) =   fem_int(37,sst,zChi(:,i,2),Var(:,i)%T)
       NLC( 5,NEQ) = LQe1

       ! Redundant heat advection

       ELM(:,:, 6,NEQ) = - fem_int(5,UgV)
       NLC( 6,NEQ) = LQe5

       ELM(:,:, 7,NEQ) = - fem_int(5,Var(:,2)%UrV)
       NLC( 7,NEQ) = LQi5

       ELM(:,:, 8,NEQ) = - fem_int(5,Var(:,3)%UrV)
       NLC( 8,NEQ) = LQz5

       ! Viscous heating

       ELM(:,:, 9,NEQ) = - 1.d-20 / rKeV * fem_int(-2,Var(:,2)%Uthhat,BnablaPi(:,2))
       NLC( 9,NEQ) = 0

       ELM(:,:,10,NEQ) = - 1.d-20 / rKeV * fem_int(-2,Var(:,3)%Uthhat,BnablaPi(:,3))
       NLC(10,NEQ) = 0

       ! Collisional transfer with ions (Energy equilibration)

       ELM(:,:,11,NEQ) = - 1.5d0 * fem_int( 2,rNuTei)
       NLC(11,NEQ) = LQe5

       ELM(:,:,12,NEQ) =   1.5d0 * fem_int(20,rNuTei,Ratden(:,i,2))
       NLC(12,NEQ) = LQi5

       ! Collisional transfer with impurities (Energy equilibration)

       ELM(:,:,13,NEQ) = - 1.5d0 * fem_int( 2,rNuTez)
       NLC(13,NEQ) = LQe5

       ELM(:,:,14,NEQ) =   1.5d0 * fem_int(20,rNuTez,Ratden(:,i,3))
       NLC(14,NEQ) = LQz5

       ! Joule heating

       ELM(:,:,15,NEQ) = - achg(i) * 1.d-3 * fem_int(-20,PsitdotVBth,Var(:,i)%n,Var(:,i)%Uthhat)
       NLC(15,NEQ) = 0

       ELM(:,:,16,NEQ) = - achg(2) * 1.d-3 * fem_int(-20,PsitdotVBth,Var(:,2)%n,Var(:,2)%Uthhat)
       NLC(16,NEQ) = 0

       ELM(:,:,17,NEQ) = - achg(3) * 1.d-3 * fem_int(-20,PsitdotVBth,Var(:,3)%n,Var(:,3)%Uthhat)
       NLC(17,NEQ) = 0

       int_tmp = fem_int(2,PsidotV)
       ELM(:,:,18,NEQ) =   achg(i) * 1.d-3 * int_tmp
       NLC(18,NEQ) = LQe7

       ELM(:,:,19,NEQ) =   achg(2) * 1.d-3 * int_tmp
       NLC(19,NEQ) = LQi7

       ELM(:,:,20,NEQ) =   achg(3) * 1.d-3 * int_tmp
       NLC(20,NEQ) = LQz7
!       ELM(:,:,18,NEQ) =   achg(i) * 1.d-3 * fem_int(-20,PsidotV,Var(:,i)%n,Var(:,i)%UphR)
!       NLC(18,NEQ) = 0
!
!       ELM(:,:,19,NEQ) =   achg(2) * 1.d-3 * fem_int(-20,PsidotV,Var(:,2)%n,Var(:,2)%UphR)
!       NLC(19,NEQ) = 0
!
!       ELM(:,:,20,NEQ) =   achg(3) * 1.d-3 * fem_int(-20,PsidotV,Var(:,3)%n,Var(:,3)%UphR)
!       NLC(20,NEQ) = 0

       ! Loss to divertor plasma region

       ELM(:,:,21,NEQ) = -                     fem_int( 2,rNuL)
       NLC(21,NEQ) = LQe5

       ELM(:,:,22,NEQ) =           PNsSUP(i) * fem_int(-2,rNuL,Var(:,i)%T)
       NLC(22,NEQ) = 0

       ELM(:,:,23,NEQ) = - 1.5d0             * fem_int( 2,rNuLTs(:,i))
       NLC(23,NEQ) = LQe5

       ELM(:,:,24,NEQ) =   1.5d0 * PTsSUP(i) * fem_int(-2,rNuLTs(:,i),Var(:,i)%n)
       NLC(24,NEQ) = 0
!!$       ELM(:,:,24,NEQ) =   1.5d0 * PTsSUP(i) * fem_int( 2,rNuLTs(:,i))
!!$       NLC(24,NEQ) = LQe1

       ! Ionization loss of n01, n02 and n03 by electron impact

       int_tmp = fem_int(20,SiVizA,Var(:,i)%n)
       ELM(:,:,25,NEQ) = - (EION * 1.D-3) * FSION * 1.D20 * int_tmp
       NLC(25,NEQ) = LQn1

       ELM(:,:,26,NEQ) = - (EION * 1.D-3) * FSION * 1.D20 * int_tmp
       NLC(26,NEQ) = LQn2

       ELM(:,:,27,NEQ) = - (EION * 1.D-3) * FSION * 1.D20 * int_tmp * BeamSW
       NLC(27,NEQ) = LQn3

       ! Collisional NBI heating (Perp + Tan)

       ELM(:,:,28,NEQ) = Eb * fem_int(-2,SNB,PNBcol_e)
       NLC(28,NEQ) = 0

       ! Simpified Alpha heating

       ELM(:,:,29,NEQ) = 1.D-20 / rKeV * fem_int(-1,PALFe)
       NLC(29,NEQ) = 0
       
       ! Direct heating (RF)

       ELM(:,:,30,NEQ) =   1.D-20 / rKeV * fem_int(-1,PRFe)
       NLC(30,NEQ) = 0

       ! Radiation loss (Bremsstrahlung)

       ELM(:,:,31,NEQ) = - 1.D-20 / rKeV * fem_int(-1,PBr)
       NLC(31,NEQ) = 0

!       ! Diffusion of electrons (***AF 2008-06-08)
!
!       ELM(:,:,32,NEQ) = - fem_int(27,sst,DMAGe)
!       NLC(32,NEQ) = LQe5

!       ! Collisional heating with beam
!
!       ELM(:,:,33,NEQ) = achgb * amqp * 1.d-3 * fem_int(-20,MNB,BUsparVbbt(:,i),PNbVinv)
!       NLC(33,NEQ) = LQb3

       NLCMAX(NEQ) = 31
    else

       ! Fixed temperature profile

       ELM(:,:,0,NEQ) = fem_int(1) * invDT
       NLC(0,NEQ) = LQe5

       NLCMAX(NEQ) = 0
    end if

  end subroutine LQe5CC

!***************************************************************
!
!   Electron Heat Parallel Flow
!
!***************************************************************

  subroutine LQe6CC

    integer(4) :: NEQ = LQe6, i = 1
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(i) * amp)

    ELM(:,:, 0,NEQ) =   2.5d0 * fem_int(2+L6,Var(:,i)%n) * invDT
    NLC( 0,NEQ) = LQe6

    ! Neoclassical viscosity force

    ELM(:,:, 1,NEQ) = - amasinv * 1.D-20 * fem_int(2+L6,xmu(:,i,2,1))
    NLC( 1,NEQ) = LQe3

    ELM(:,:, 2,NEQ) = - amasinv * 1.D-20 * fem_int(2+L6,xmu(:,i,2,2))
    NLC( 2,NEQ) = LQe6

    ! Diamagnetic forces

    ELM(:,:, 3,NEQ) =   amasinv * 1.D-20 * fem_int(2+L6,xmu(:,i,2,1))
    NLC( 3,NEQ) = LQe8

    ELM(:,:, 4,NEQ) = - rKilo / achg(i) * amasinv * 1.D-20 * fem_int(22+L6,fipolsdtPNsV(:,i),xmu(:,i,2,2))
    NLC( 4,NEQ) = LQe5

    ELM(:,:, 5,NEQ) =   rKilo / achg(i) * amasinv * 1.D-20 * fem_int(32+L6,fipolsdtPNsV(:,i),Var(:,i)%T,xmu(:,i,2,2))
    NLC( 5,NEQ) = LQe1

    ! Collisional friction force with electrons

    ELM(:,:, 6,NEQ) = - fem_int(20+L6,lab(:,i,1,2,1),Var(:,i)%n)
    NLC( 6,NEQ) = LQe3

    ! Collisional friction force with ions

    ELM(:,:, 7,NEQ) = - fem_int(20+L6,lab(:,i,2,2,1),Var(:,i)%n)
    NLC( 7,NEQ) = LQi3

    ! Collisional friction force with impurities

    ELM(:,:, 8,NEQ) = - fem_int(20+L6,lab(:,i,3,2,1),Var(:,i)%n)
    NLC( 8,NEQ) = LQz3

    ! Collisional friction force with beam ions (No counterpart; LQb6)
    ! It may be retained because it is required to yield an adequate shielding factor.

    ELM(:,:, 9,NEQ) = - fem_int(20+L6,laf(:,i,2,1),  Ratden(:,i,4))
    NLC( 9,NEQ) = LQb3

    ! Collisional heat friction force with electrons (self collision)

    ELM(:,:,10,NEQ) =   fem_int(20+L6,lab(:,i,1,2,2),Var(:,i)%n)
    NLC(10,NEQ) = LQe6

    ! Collisional heat friction force with ions

    ELM(:,:,11,NEQ) =   fem_int(20+L6,lab(:,i,2,2,2),Var(:,i)%n)
    NLC(11,NEQ) = LQi6

    ! Collisional heat friction force with impurities

    ELM(:,:,12,NEQ) =   fem_int(20+L6,lab(:,i,3,2,2),Var(:,i)%n)
    NLC(12,NEQ) = LQz6

    ! Loss to divertor plasma region

    ELM(:,:,13,NEQ) = - fem_int(20+L6,rNuL,Var(:,i)%n) * fact
    NLC(13,NEQ) = LQe6

    ! Collisional friction force with neutrals

    ELM(:,:,14,NEQ) = - fem_int(20+L6,rNu0s(:,i),Var(:,i)%n)
    NLC(14,NEQ) = LQe6

    NLCMAX(NEQ) = 14

    if( MDFIXT /= 0 ) then ! Overwrite Pressure gradient force
       ELM(:,:, 4,NEQ) = - rKilo / achg(i) * amasinv * 1.D-20 * fem_int(-(22+L6),fipolsdt,xmu(:,i,2,2),Var(:,i)%T)
       NLC( 4,NEQ) = 0

       ELM(:,:, 5,NEQ) =   0.d0
       NLC( 5,NEQ) = 0
    end if

  end subroutine LQe6CC

!***************************************************************
!
!   Electron rotation frequency: Ne <UePhi/R>
!
!***************************************************************

  subroutine LQe7CC

    integer(4) :: NEQ = LQe7, i = 1

    ! Electron rotation frequency

    ELM(:,:,1,NEQ) = - fem_int(1)
    NLC(1,NEQ) = LQe7

    ! Electron rotation

    ELM(:,:,2,NEQ) =   fem_int(2,rrtinv)
    NLC(2,NEQ) = LQe4

    ELM(:,:,3,NEQ) = - fem_int(20,Fqhatsq,rrtinv)
    NLC(3,NEQ) = LQe4

    ! Electron parallel velocity

    ELM(:,:,4,NEQ) =   fem_int(30,Fqhatsq,fipolinv,Var(:,i)%n)
    NLC(4,NEQ) = LQe3

    NLCMAX(NEQ) = 4

  end subroutine LQe7CC

!***************************************************************
!
!   Electron diamagnetic flow: BV1
!
!***************************************************************

  subroutine LQe8CC

    integer(4) :: NEQ = LQe8, i = 1

    ! Diamagnetic flow

    ELM(:,:,1,NEQ) = - fem_int(2+L8,Var(:,i)%n)
    NLC(1,NEQ) = LQe8

    ! Pressure gradient

    ELM(:,:,2,NEQ) = - rKilo / achg(i) * fem_int(5+L8,fipolsdt)
    NLC(2,NEQ) = LQe5

    ! Electrostatic potential gradient

    ELM(:,:,3,NEQ) = - fem_int(22+L8,fipolsdt,Var(:,i)%n)
    NLC(3,NEQ) = LQm1

    NLCMAX(NEQ) = 3

    if( MDFIXT /= 0 ) then ! Overwrite Pressure gradient force
       ELM(:,:,2,NEQ) = - rKilo / achg(i) &
            &                  * (  fem_int(21+L8,fipolsdt,Var(:,i)%T) &
            &                     + fem_int(22+L8,fipolsdt,Var(:,i)%T))
       NLC(2,NEQ) = LQe1
    end if

  end subroutine LQe8CC

!***************************************************************
!
!   Ion Density Equation
!
!***************************************************************

  subroutine LQi1CC

    integer(4) :: NEQ = LQi1, i = 2
    real(8), dimension(1:NEMAX,1:4) :: int_tmp

    ELM(:,:, 0,NEQ) = fem_int(1) * invDT
    NLC( 0,NEQ) = LQi1

    ! Divergence

    ELM(:,:, 1,NEQ) = - fem_int(4)
    NLC( 1,NEQ) = LQi2

    ! Divergence + Grid velocity

    ELM(:,:, 2,NEQ) =   fem_int(5,UgV)
    NLC( 2,NEQ) = LQi1

    ! Ionization of n01, n02 and n03 by electron impact

    int_tmp = fem_int(20,SiVizA,Var(:,1)%n)
    ELM(:,:, 3,NEQ) =     FSION * 1.D20 * int_tmp
    NLC( 3,NEQ) = LQn1

    ELM(:,:, 4,NEQ) =     FSION * 1.D20 * int_tmp * ThntSW
    NLC( 4,NEQ) = LQn2

    ELM(:,:, 5,NEQ) =     FSION * 1.D20 * int_tmp * BeamSW
    NLC( 5,NEQ) = LQn3

    ! Loss to divertor plasma region

    ELM(:,:, 6,NEQ) = -             fem_int( 2,rNuL)
    NLC( 6,NEQ) = LQi1

    ELM(:,:, 7,NEQ) =   PNsSUP(i) * fem_int(-1,rNuL)
    NLC( 7,NEQ) = 0

    ! Particle source from beam ion

    ELM(:,:, 8,NEQ) =   fem_int(2,rNuB)
    NLC( 8,NEQ) = LQb1

    ! Particle source from ripple trapped beam ions

    ELM(:,:, 9,NEQ) =   fem_int(20,rNuB,rip_rat) * RpplSW
    NLC( 9,NEQ) = LQr1

    ! NBI kick up ions (Charge exchange)

!    ELM(:,:,10,NEQ) = - RatCX * fem_int(-1,SNBi)
    ELM(:,:,10,NEQ) = - fem_int(-1,RatCXSNBi)
    NLC(10,NEQ) = 0

    ! Loss cone loss

    ELM(:,:,11,NEQ) =   fem_int(-1,SiLC)
    NLC(11,NEQ) = 0
 
    ! Ion orbit loss

    ELM(:,:,12,NEQ) = - fem_int(2,rNuOL)
    NLC(12,NEQ) = LQi1

!    ! Diffusion of ions (***AF 2008-06-08)
!
!    ELM(:,:,13,NEQ) = - fem_int(27,sst,DMAGi)
!    NLC(13,NEQ) = LQi1

!    ! Parallel loss reduction due to the potential
!    ! induced by the parallel loss of the beam ions
!
!    ELM(:,:,14,NEQ) =   fem_int(2,rNuLB) * BeamSW
!    NLC(14,NEQ) = LQb1

    NLCMAX(NEQ) = 12

  end subroutine LQi1CC

!***************************************************************
!
!   Ion Radial Flow
!
!***************************************************************
  
  subroutine LQi2CC

    integer(4) :: NEQ = LQi2, i = 2
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(i) * amp)

    ! Diamagnetic force

    ELM(:,:,1,NEQ) =   achg(i) * aee * amasinv * fem_int(30,bri,Var(:,i)%n,fipolinv)
    NLC(1,NEQ) = LQi8

    ! v x B force

    ELM(:,:,2,NEQ) =   achg(i) * aee * amasinv * fem_int(20,fipol,Var(:,i)%n)
    NLC(2,NEQ) = LQi3

    ELM(:,:,3,NEQ) = - achg(i) * aee * amasinv * fem_int( 2,bbt)
    NLC(3,NEQ) = LQi4

    NLCMAX(NEQ) = 3

  end subroutine LQi2CC

!***************************************************************
!
!   Ion Parallel Flow
!
!***************************************************************

  subroutine LQi3CC

    integer(4) :: NEQ = LQi3, i = 2
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(i) * amp)

    ELM(:,:, 0,NEQ) =   fem_int(2+L3,Var(:,i)%n) * invDT
    NLC( 0,NEQ) = LQi3

    ! Neoclassical viscosity force

    ELM(:,:, 1,NEQ) = - amasinv * 1.D-20 * fem_int(2+L3,xmu(:,i,1,1))
    NLC( 1,NEQ) = LQi3

    ELM(:,:, 2,NEQ) = - amasinv * 1.D-20 * fem_int(2+L3,xmu(:,i,1,2))
    NLC( 2,NEQ) = LQi6

    ! Diamagnetic forces (particle)

    ELM(:,:, 3,NEQ) =   amasinv * 1.D-20 * fem_int(2+L3,xmu(:,i,1,1))
    NLC( 3,NEQ) = LQi8

    ! Diamagnetic forces (heat)

    ELM(:,:, 4,NEQ) = - rKilo / achg(i) * amasinv * 1.D-20 * fem_int(22+L3,fipolsdtPNsV(:,i),xmu(:,i,1,2))
    NLC( 4,NEQ) = LQi5

    ELM(:,:, 5,NEQ) =   rKilo / achg(i) * amasinv * 1.D-20 * fem_int(32+L3,fipolsdtPNsV(:,i),Var(:,i)%T,xmu(:,i,1,2))
    NLC( 5,NEQ) = LQi1
!    ELM(:,:, 5,NEQ) =   rKilo / achg(i) * amasinv * 1.D-20 * fem_int(30+L3,fipolsdtPNsV(:,i),dlnPNsV(:,i),xmu(:,i,1,2))
!    NLC( 5,NEQ) = LQi5

    ! Collisional friction force with electrons

    ELM(:,:, 6,NEQ) =   fem_int(20+L3,lab(:,i,1,1,1),Var(:,i)%n)
    NLC( 6,NEQ) = LQe3

    ! Collisional friction force with ions (self collision)

    ELM(:,:, 7,NEQ) =   fem_int(20+L3,lab(:,i,2,1,1),Var(:,i)%n)
    NLC( 7,NEQ) = LQi3

    ! Collisional friction force with impurities

    ELM(:,:, 8,NEQ) =   fem_int(20+L3,lab(:,i,3,1,1),Var(:,i)%n)
    NLC( 8,NEQ) = LQz3

    ! Collisional friction force with beam ions

    ELM(:,:, 9,NEQ) =   fem_int(20+L3,laf(:,i,1,1),  Ratden(:,i,4))
    NLC( 9,NEQ) = LQb3

    ! Collisional heat friction force with electrons

    ELM(:,:,10,NEQ) = - fem_int(20+L3,lab(:,i,1,1,2),Var(:,i)%n)
    NLC(10,NEQ) = LQe6

    ! Collisional heat friction force with ions (self collision)

    ELM(:,:,11,NEQ) = - fem_int(20+L3,lab(:,i,2,1,2),Var(:,i)%n)
    NLC(11,NEQ) = LQi6

    ! Collisional heat friction force with impurities

    ELM(:,:,12,NEQ) = - fem_int(20+L3,lab(:,i,3,1,2),Var(:,i)%n)
    NLC(12,NEQ) = LQz6

    ! Electric field forces

    ELM(:,:,13,NEQ) = - achg(i) * aee * amasinv * fem_int(30+L3,fipol,Var(:,i)%n,aatq)
    NLC(13,NEQ) = LQm2

    ELM(:,:,14,NEQ) =   achg(i) * aee * amasinv * fem_int(30+L3,fipol,Var(:,i)%n,aat)
    NLC(14,NEQ) = LQm3

    ! Loss to divertor plasma region

    ELM(:,:,15,NEQ) = - fem_int(20+L3,rNuL,Var(:,i)%n) * fact
    NLC(15,NEQ) = LQi3

    ! Increase in momentum due to thermalization of fast ions

    ELM(:,:,16,NEQ) =   amb / amas(i) * fem_int(2+L3,rNuB) * BeamSW
    NLC(16,NEQ) = LQb3

    ! Collisional friction force with neutrals

    ELM(:,:,17,NEQ) = - fem_int(20+L3,rNu0s(:,i),Var(:,i)%n)
    NLC(17,NEQ) = LQi3

    ! Charge exchange force

    ELM(:,:,18,NEQ) = - fem_int(20+L3,rNuiCX,Var(:,i)%n)
    NLC(18,NEQ) = LQi3

    ! Loss cone loss

    ELM(:,:,19,NEQ) =   fem_int(-(2+L3),SiLCB,Var(:,i)%n)
    NLC(19,NEQ) = 0

    ! Ion orbit loss

    ELM(:,:,20,NEQ) = - fem_int(20+L3,rNuOL,Var(:,i)%n)
    NLC(20,NEQ) = LQi3

    ! Additional torque input (sheared flow)

    ELM(:,:,21,NEQ) =   amasinv * 1.D-20 * fem_int(-(1+L3),Tqp)
    NLC(21,NEQ) = 0

    ! Additional torque input

    ELM(:,:,22,NEQ) =   amasinv * 1.D-20 * fem_int(-(2+L3),Tqt,ZetatoPara)
    NLC(22,NEQ) = 0

    ! ---

    ! Perpendicular viscosity force: diffusion

    ELM(:,:,23,NEQ) = - FSPARV(i) * fem_int(38,sst*rMus(:,i)*Var(:,i)%n,ZetatoPara,ParatoZeta)
    NLC(23,NEQ) = LQi3

    ELM(:,:,24,NEQ) = - FSPARV(i) * fem_int(34,sst*rMus(:,i)*Var(:,i)%n,dZetatoPara,ParatoZeta)
    NLC(24,NEQ) = LQi3

    ! Perpendicular viscosity force: pinch

    ELM(:,:,25,NEQ) = - FSPARV(i) * fem_int(24,suft*Vmps(:,i),Var(:,i)%n)
    NLC(25,NEQ) = LQi3

    ELM(:,:,26,NEQ) =   FSPARV(i) * fem_int(30,suft*Vmps(:,i)*Var(:,i)%n,dZetatoPara,ParatoZeta)
    NLC(26,NEQ) = LQi3

    ! Perpendicular viscosity force: residual stress

    ELM(:,:,27,NEQ) = - amasinv * fem_int(-5,ZetatoPara,PiRess(:,i))
    NLC(27,NEQ) = 0

    ! Perpendicular viscosity force: terms associated with diamagnetic flow

    ELM(:,:,28,NEQ) = - FSPARV(i) * fem_int(-2,ZetatoPara,dtormflux(:,i))
    NLC(28,NEQ) = 0

!!$    ! Alternative implementation of perpendicular viscosity force: terms associated with diamagnetic flow
!!$
!!$    ELM(:,:,28,NEQ) = - FSPARV(i) * fem_int(38,sst *rMus(:,i)*Var(:,i)%n, ZetatoPara,ribi)
!!$    NLC(28,NEQ) = LQi8
!!$
!!$    ELM(:,:,29,NEQ) = - FSPARV(i) * fem_int(34,sst *rMus(:,i)*Var(:,i)%n,dZetatoPara,ribi)
!!$    NLC(29,NEQ) = LQi8
!!$
!!$    ELM(:,:,30,NEQ) = - FSPARV(i) * fem_int(24,suft*Vmps(:,i)*Var(:,i)%n, ZetatoPara*ribi)
!!$    NLC(30,NEQ) = LQi8
!!$
!!$    ELM(:,:,31,NEQ) =   FSPARV(i) * fem_int(30,suft*Vmps(:,i)*Var(:,i)%n,dZetatoPara,ribi)
!!$    NLC(31,NEQ) = LQi8

    ! ---

!!$    ! Perpendicular viscosity (off-diagonal)
!!$
!!$    ELM(:,:,23,NEQ) = - FSPARV(i) * fem_int(37,ZetatoPara,sst,rMus(:,i))
!!$    NLC(23,NEQ) = LQi4
!!$
!!$    ELM(:,:,24,NEQ) =   FSPARV(i) * fem_int(37,ZetatoPara,sst*rMus(:,i),Var(:,i)%RUph)
!!$    NLC(24,NEQ) = LQi1
!!$  
!!$    ELM(:,:,25,NEQ) = - FSPARV(i) * fem_int(32,dZetatoPara,sst,rMus(:,i))
!!$    NLC(25,NEQ) = LQi4
!!$
!!$    ELM(:,:,26,NEQ) =   FSPARV(i) * fem_int(32,dZetatoPara,sst*rMus(:,i),Var(:,i)%RUph)
!!$    NLC(26,NEQ) = LQi1
!!$
!!$    ! Momentum pinch (off-diagonal)
!!$
!!$    ELM(:,:,27,NEQ) = - FSPARV(i) * fem_int(24,suft*Vmps(:,i),ZetatoPara)
!!$    NLC(27,NEQ) = LQi4
!!$
!!$    ELM(:,:,28,NEQ) =   FSPARV(i) * fem_int(30,suft,Vmps(:,i),dZetatoPara)
!!$    NLC(28,NEQ) = LQi4

!    ! Diffusion of ions (***AF 2008-06-08) (37+L3 is not possible because it is a diffusion term)
!
!    ELM(:,:,29,NEQ) = - fem_int(37,sst,DMAGi,Var(:,i)%n)
!    NLC(29,NEQ) = LQi3

    NLCMAX(NEQ) = 28

    if( MDFIXT /= 0 ) then ! Overwrite Pressure gradient force
       ELM(:,:, 4,NEQ) = - rKilo / achg(i) * amasinv * 1.D-20 * fem_int(-(22+L3),fipolsdt,xmu(:,i,1,2),Var(:,i)%T)
       NLC( 4,NEQ) = 0

       ELM(:,:, 5,NEQ) =   0.d0
       NLC( 5,NEQ) = 0
    end if

  end subroutine LQi3CC

!***************************************************************
!
!   Ion Toroidal Flow
!
!***************************************************************

  subroutine LQi4CC

    integer(4) :: NEQ = LQi4, i = 2, j, NCcomp
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(i) * amp)

    ! Uiphi'(0) : 0

    ELM(:,:, 0,NEQ) = fem_int(1) * invDT
    NLC( 0,NEQ) = LQi4

    ! Advection + Grid velocity

    ELM(:,:, 1,NEQ) = - fem_int(5,UsrgV(:,i))
    NLC( 1,NEQ) = LQi4

    ELM(:,:, 2,NEQ) = - FSADV * fem_int(6,Var(:,i)%UrV)
    NLC( 2,NEQ) = LQi4

    ! Momentum pinch

    ELM(:,:, 3,NEQ) = - fem_int(24,suft,Vmps(:,i))
    NLC( 3,NEQ) = LQi4

    ! Viscosity force

    ELM(:,:, 4,NEQ) = - fem_int(27,sst,rMus(:,i))
    NLC( 4,NEQ) = LQi4

    ELM(:,:, 5,NEQ) =   fem_int(37,sst,rMus(:,i),Var(:,i)%RUph)
    NLC( 5,NEQ) = LQi1
!!$    ELM(:,:, 5,NEQ) = - fem_int(24,sst*rMus(:,i),dlnPNsV(:,i))
!!$    NLC( 5,NEQ) = LQi4

    ! Residual stress

    ELM(:,:, 6,NEQ) = - amasinv * fem_int(-4,PiRess(:,i))
    NLC( 6,NEQ) = 0

    ! Collisional friction force with electrons

    ELM(:,:, 7,NEQ) =   fem_int(20,lab(:,i,1,1,1),Ratden(:,i,1))
    NLC( 7,NEQ) = LQe4

    ! Collisional friction force with bulk ions (self collision)

    ELM(:,:, 8,NEQ) =   fem_int( 2,lab(:,i,2,1,1))
    NLC( 8,NEQ) = LQi4

    ! Collisional friction force with impurities

    ELM(:,:, 9,NEQ) =   fem_int(20,lab(:,i,3,1,1),Ratden(:,i,3))
    NLC( 9,NEQ) = LQz4

    ! Collisional friction force with beam ions

!    ELM(:,:,10,NEQ) =   fem_int(30,ParatoZeta,laf(:,i,1,1),Ratden(:,i,4))
!    NLC(10,NEQ) = LQb3
    ELM(:,:,10,NEQ) =   fem_int(20,laf(:,i,1,1),Ratden(:,i,4))
    NLC(10,NEQ) = LQb4

    ! Collisional heat friction force with electrons

    j = 1
    ELM(:,:,11,NEQ) =   rKilo / achg(j) * fem_int(-22,lab(:,i,j,1,2),ribsdtNs(:,i),Var(:,j)%T)
    NLC(11,NEQ) = 0
!!$    ELM(:,:,11,NEQ) =   rKilo / achg(j) * fem_int(31,lab(:,i,j,1,2),ribsdt,Var(:,j)%T)
!!$    NLC(11,NEQ) = LQi1

    ELM(:,:,12,NEQ) = -                   fem_int(30,lab(:,i,j,1,2),ParatoZeta,Var(:,i)%n)
    NLC(12,NEQ) = LQe6

    ! Collisional heat friction force with ions (self collision)

    j = 2
    ELM(:,:,13,NEQ) =   rKilo / achg(j) * fem_int(-22,lab(:,i,j,1,2),ribsdtNs(:,i),Var(:,j)%T)
    NLC(13,NEQ) = 0
!!$    ELM(:,:,13,NEQ) =   rKilo / achg(j) * fem_int(31,lab(:,i,j,1,2),ribsdt,Var(:,j)%T)
!!$    NLC(13,NEQ) = LQi1

    ELM(:,:,14,NEQ) = -                   fem_int(30,lab(:,i,j,1,2),ParatoZeta,Var(:,i)%n)
    NLC(14,NEQ) = LQi6

    ! Collisional heat friction force with impurities

    j = 3
    ELM(:,:,15,NEQ) =   rKilo / achg(j) * fem_int(-22,lab(:,i,j,1,2),ribsdtNs(:,i),Var(:,j)%T)
    NLC(15,NEQ) = 0
!!$    ELM(:,:,15,NEQ) =   rKilo / achg(j) * fem_int(31,lab(:,i,j,1,2),ribsdt,Var(:,j)%T)
!!$    NLC(15,NEQ) = LQi1

    ELM(:,:,16,NEQ) = -                   fem_int(30,lab(:,i,j,1,2),ParatoZeta,Var(:,i)%n)
    NLC(16,NEQ) = LQz6

    ! Toroidal E force

    ELM(:,:,17,NEQ) =   achg(i) * aee * amasinv * fem_int(2,Var(:,i)%n)
    NLC(17,NEQ) = LQm3

    ! v x B force

    ELM(:,:,18,NEQ) =   achg(i) * aee * amasinv * fem_int(6,PsiV)
    NLC(18,NEQ) = LQi2

    ! Turbulent particle transport driver (electron driven)

    if( FSTPTM(i) == 0 ) then
       ELM(:,:,19,NEQ) =   achg(i) * rKeV * amasinv * fem_int(32,FQLcoef,Var(:,1)%T,Ratdenhyd) &
            &                           * (1.d0 - FSVAHL)
       NLC(19,NEQ) = LQe1

       ELM(:,:,20,NEQ) =   achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Ratdenhyd) &
            &                           * FSVAHL
       NLC(20,NEQ) = LQe5

    else if( FSTPTM(i) == 1 ) then
       ELM(:,:,19,NEQ) =   achg(i) * aee * amasinv * fem_int(30,FQLcoef2,ZetatoPara,Ratdenhyd)
       NLC(19,NEQ) = LQe4

       ELM(:,:,20,NEQ) = - achg(i) * aee * amasinv * fem_int(20,FQLcoef2,Var(:,i)%n+PNbV(:))
       NLC(20,NEQ) = LQe3

       ELM(:,:,21,NEQ) =   achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Ratdenhyd) &
            &                           * (FSVAHL - 1.d0)
       NLC(21,NEQ) = LQe5

       ELM(:,:,22,NEQ) = - achg(i) * rKeV * amasinv * fem_int(32,FQLcoef,Var(:,1)%T,Ratdenhyd) &
            &                           * (FSVAHL - 1.d0)
       NLC(22,NEQ) = LQe1

       ELM(:,:,23,NEQ) =   achg(i) * aee * amasinv * fem_int(22,FQLcoef,Var(:,i)%n+PNbV(:))
       NLC(23,NEQ) = LQm1

    else
       ELM(:,:,19,NEQ) =   achg(i) * aee * amasinv * fem_int(20,FQLcoef1,Var(:,i)%n+PNbV(:))
       NLC(19,NEQ) = LQe8

       ELM(:,:,20,NEQ) =   achg(i) * aee * amasinv * fem_int(22,FQLcoef,Var(:,i)%n+PNbV(:))
       NLC(20,NEQ) = LQm1

       ELM(:,:,21,NEQ) =   achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Ratdenhyd) &
            &                           * (FSVAHL - 1.d0)
       NLC(21,NEQ) = LQe5

       ELM(:,:,22,NEQ) = - achg(i) * rKeV * amasinv * fem_int(32,FQLcoef,Var(:,1)%T,Ratdenhyd) &
            &                           * (FSVAHL - 1.d0)
       NLC(22,NEQ) = LQe1

    end if

    !  -- Turbulent pinch term

    ELM(:,:,24,NEQ) =   achg(i) * aee * amasinv * fem_int(20,FQLcoef,FVpch)
    NLC(24,NEQ) = LQi1

    ELM(:,:,25,NEQ) =   achg(i) * aee * amasinv * fem_int(20,FQLcoef,FVpch)
    NLC(25,NEQ) = LQb1

    ! Loss to divertor plasma region

    ELM(:,:,26,NEQ) = - fem_int(2,rNuL) * fact
    NLC(26,NEQ) = LQi4

    ! Increase in momentum due to thermalization of fast ions

!    ELM(:,:,27,NEQ) =   fem_int(20,ParatoZeta,rNuB) * BeamSW
!    NLC(27,NEQ) = LQb3

    ELM(:,:,27,NEQ) =   amb / amas(i) * fem_int(2,rNuB) * BeamSW
    NLC(27,NEQ) = LQb4

    ! Collisional friction force with neutrals

    ELM(:,:,28,NEQ) = - fem_int(2,rNu0s(:,i))
    NLC(28,NEQ) = LQi4

    ! Charge exchange force

    ELM(:,:,29,NEQ) = - fem_int(2,rNuiCX)
    NLC(29,NEQ) = LQi4

    ! Loss cone loss

    ELM(:,:,30,NEQ) =   fem_int(-1,SiLCph)
    NLC(30,NEQ) = 0

    ! Ion orbit loss

    ELM(:,:,31,NEQ) = - fem_int(2,rNuOL)
    NLC(31,NEQ) = LQi4

    !  Additional torque input

    ELM(:,:,32,NEQ) =   amasinv * 1.D-20 * fem_int(-1,Tqt) 
    NLC(32,NEQ) = 0

!    !  Diffusion of ions (***AF 2008-06-08)
!
!    ELM(:,:,33,NEQ) = - fem_int(27,sst,DMAGi)
!    NLC(33,NEQ) = LQi4

    NLCMAX(NEQ) = 32

    ! Multi-ion composition channel (FSTPTM = 3).
    !
    !   Delta K_a = sum_{b /= a, thermal ions} L_ab ( X_a - X_b ),
    !   X_a = (1/n_a) dn_a/dV,
    !   L_ab = L_ba = e (dpsi/dV) <|grad V|^2> D_e^QL (Z_a n_a)(Z_b n_b)/n_e.
    !
    ! Each pair contributes equal and opposite forces, so sum_a Delta K_a = 0
    ! holds element matrix by element matrix over whatever set of thermal ions
    ! is present: no net toroidal momentum is added, and since
    ! Gamma_a = -<R K_a1zeta>/(Z_a e dpsi/dV) that is the same statement as
    ! sum_a Z_a Delta Gamma_a = 0.  Beam ions are left out - they are too fast
    ! to be mixed by the background fluctuations - and their absence from the
    ! sum does not spoil the cancellation.
    !
    ! The 1/n_a of X_a cancels against the Z_a n_a of L_ab, so nothing here
    ! divides by a trace density.  With a single thermal ion the sum is empty
    ! and FSTPTM = 3 reduces exactly to FSTPTM = 2.
    if( FSTPTM(i) == 3.d0 ) then
       NCcomp = NLCMAX(NEQ)
       do j = 2, NSM
          if( j == i ) cycle
          NCcomp = NCcomp + 1
          ELM(:,:,NCcomp,NEQ) =   achg(i) * achg(j) * aee * amasinv &
               &                * fem_int(22,FQCcoef,Ratden(:,j,1))
          NLC(NCcomp,NEQ) = LQsn(i)
          NCcomp = NCcomp + 1
          ELM(:,:,NCcomp,NEQ) = - achg(i) * achg(j) * aee * amasinv &
               &                * fem_int(22,FQCcoef,Ratden(:,i,1))
          NLC(NCcomp,NEQ) = LQsn(j)
       end do
       NLCMAX(NEQ) = NCcomp
    end if


    if( MDFIXT /= 0 ) then ! Overwrite Pressure gradient force
       if( FSTPTM(i) == 0 ) then
          ELM(:,:,19,NEQ) =   achg(i) * rKeV * amasinv * fem_int(32,FQLcoef,Var(:,1)%T,Ratdenhyd)
          NLC(19,NEQ) = LQe1

          ELM(:,:,20,NEQ) =   achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Var(:,i)%n+PNbV(:)) &
               &                           * FSVAHL
          NLC(20,NEQ) = LQe5

       else
          ELM(:,:,21,NEQ) =   achg(i) * rKeV * amasinv * fem_int(21,FQLcoef,Var(:,1)%T) &
               &                           * (FSVAHL - 1.d0)
          NLC(21,NEQ) = LQi1

          ELM(:,:,22,NEQ) =   achg(i) * rKeV * amasinv * fem_int(21,FQLcoef,Var(:,1)%T) &
               &                           * (FSVAHL - 1.d0) * BeamSW
          NLC(22,NEQ) = LQb1

       end if
    end if

  end subroutine LQi4CC

!***************************************************************
!
!  Ion Energy Transport: Ti
!
!***************************************************************

  subroutine LQi5CC

    integer(4) :: NEQ = LQi5, i = 2

    ! Temperature evolution

    if( MDFIXT == 0 ) then
       ELM(:,:, 0,NEQ) =   1.5d0 * fem_int(1) * invDT
       NLC( 0,NEQ) = LQi5

       ! Advection + Grid velocity

       ELM(:,:, 1,NEQ) = - 2.5d0 * fem_int(5,UsrgV(:,i))
       NLC( 1,NEQ) = LQi5

       ELM(:,:, 2,NEQ) = - 2.5d0 * FSADV * fem_int(6,Var(:,i)%UrV)
       NLC( 2,NEQ) = LQi5

       ! Heat pinch

       ELM(:,:, 3,NEQ) = - fem_int(24,Vhps(:,i),suft)
       NLC( 3,NEQ) = LQi5

       ! Conduction transport

       ELM(:,:, 4,NEQ) = - fem_int(27,sst,zChi(:,i,1))
       NLC( 4,NEQ) = LQi5

       ELM(:,:, 5,NEQ) =   fem_int(37,sst,zChi(:,i,2),Var(:,i)%T)
       NLC( 5,NEQ) = LQi1

       ! Redundant heat convection term

       ELM(:,:, 6,NEQ) =   fem_int(5,UsrgV(:,i))
       NLC( 6,NEQ) = LQi5

       ! Viscous heating

       ELM(:,:, 7,NEQ) =   1.D-20 / rKeV * fem_int(-2,Var(:,i)%Uthhat,BnablaPi(:,i))
       NLC( 7,NEQ) = 0

       ! Collisional transfer with electrons (Energy equilibration)

       ELM(:,:, 8,NEQ) = - 1.5d0 * fem_int(20,rNuTei,Ratden(:,1,i))
       NLC( 8,NEQ) = LQi5

       ELM(:,:, 9,NEQ) =   1.5d0 * fem_int( 2,rNuTei)
       NLC( 9,NEQ) = LQe5

       ! Collisional transfer with impurities (Energy equilibration)

       ELM(:,:,10,NEQ) = - 1.5d0 * fem_int( 2,rNuTiz)
       NLC(10,NEQ) = LQi5

       ELM(:,:,11,NEQ) =   1.5d0 * fem_int(20,rNuTiz,Ratden(:,i,3))
       NLC(11,NEQ) = LQz5

       ! Loss to divertor plasma region

       ELM(:,:,12,NEQ) = -             fem_int( 2,rNuL)
       NLC(12,NEQ) = LQi5

       ELM(:,:,13,NEQ) =   PNsSUP(i) * fem_int(-2,rNuL,Var(:,i)%T)
       NLC(13,NEQ) = 0

       ELM(:,:,14,NEQ) = - 1.5d0            * fem_int( 2,rNuLTs(:,i))
       NLC(14,NEQ) = LQi5

       ELM(:,:,15,NEQ) =   1.5d0 * PTsSUP(i) * fem_int(-2,rNuLTs(:,i),Var(:,i)%n)
       NLC(15,NEQ) = 0
!       ELM(:,:,15,NEQ) =   1.5d0 * PTsSUP(i)    * fem_int( 2,rNuLTs(:,i))
!       NLC(15,NEQ) = LQi1

       ! Ionization heating of n01, n02 and n03 by electron impact

       ELM(:,:,16,NEQ) = 1.5d0 * FSION * 1.D20 * fem_int(30,SiVizA,Var(:,1)%n,PT01V)
       NLC(16,NEQ) = LQn1

       ELM(:,:,17,NEQ) = 1.5d0 * FSION * 1.D20 * fem_int(30,SiVizA,Var(:,1)%n,PT02V)
       NLC(17,NEQ) = LQn2

       ELM(:,:,18,NEQ) = 1.5d0 * FSION * 1.D20 * fem_int(30,SiVizA,Var(:,1)%n,PT03V) * BeamSW
       NLC(18,NEQ) = LQn3

       ! Charge exchange loss due to slow neutrals
       !   (Thermal neutrals are assumed to have same temperature with ions.)

       ELM(:,:,19,NEQ) = - 1.5d0 * fem_int( 2,rNuiCXT)
       NLC(19,NEQ) = LQi5

       ELM(:,:,20,NEQ) =   1.5d0 * fem_int(-20,rNuiCXT(:),Var(:,i)%n,PT01V(:))
       NLC(20,NEQ) = 0

       ! Collisional NBI heating (Perp + Tan)

       ELM(:,:,21,NEQ) = Eb * fem_int(-2,SNB,PNBcol_i)
       NLC(21,NEQ) = 0

       ! Heating due to beam momentum deposition

       ELM(:,:,22,NEQ) = achgb * amqp * 1.d-3 * fem_int(-20,MNB,BUsparVbbt(:,i),PNbVinv)
       NLC(22,NEQ) = LQb3

       ! Simplified Alpha heating

       ELM(:,:,23,NEQ) = 1.D-20 / rKeV * fem_int(-1,PALFi)
       NLC(23,NEQ) = 0

       ! Direct heating (RF)

       ELM(:,:,24,NEQ) = 1.D-20 / rKeV * fem_int(-1,PRFi)
       NLC(24,NEQ) = 0

!       !  Diffusion of ions (***AF 2008-06-08)
!
!       ELM(:,:,25,NEQ) = - fem_int(27,sst,DMAGi)
!       NLC(25,NEQ) = LQi5

       NLCMAX(NEQ) = 24
    else

       ! Fixed temperature profile

       ELM(:,:,0,NEQ) = fem_int(1) * invDT
       NLC(0,NEQ) = LQi5

       NLCMAX(NEQ) = 0
    end if

  end subroutine LQi5CC

!***************************************************************
!
!   Ion Heat Parallel Flow
!
!***************************************************************

  subroutine LQi6CC

    integer(4) :: NEQ = LQi6, i = 2
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(i) * amp)

    ELM(:,:, 0,NEQ) =   2.5d0 * fem_int(2+L6,Var(:,i)%n) * invDT
    NLC( 0,NEQ) = LQi6

    ! Neoclassical viscosity force

    ELM(:,:, 1,NEQ) = - amasinv * 1.D-20 * fem_int(2+L6,xmu(:,i,2,1))
    NLC( 1,NEQ) = LQi3

    ELM(:,:, 2,NEQ) = - amasinv * 1.D-20 * fem_int(2+L6,xmu(:,i,2,2))
    NLC( 2,NEQ) = LQi6

    ! Diamagnetic forces

    ELM(:,:, 3,NEQ) =   amasinv * 1.D-20 * fem_int(2+L6,xmu(:,i,2,1))
    NLC( 3,NEQ) = LQi8

    ELM(:,:, 4,NEQ) = - rKilo / achg(i) * amasinv * 1.D-20 * fem_int(22+L6,fipolsdtPNsV(:,i),xmu(:,i,2,2))
    NLC( 4,NEQ) = LQi5

    ELM(:,:, 5,NEQ) =   rKilo / achg(i) * amasinv * 1.D-20 * fem_int(32+L6,fipolsdtPNsV(:,i),Var(:,i)%T,xmu(:,i,2,2))
    NLC( 5,NEQ) = LQi1

    ! Collisional friction force with electrons

    ELM(:,:, 6,NEQ) = - fem_int(20+L6,lab(:,i,1,2,1),Var(:,i)%n)
    NLC( 6,NEQ) = LQe3

    ! Collisional friction force with ions (self collision)

    ELM(:,:, 7,NEQ) = - fem_int(20+L6,lab(:,i,2,2,1),Var(:,i)%n)
    NLC( 7,NEQ) = LQi3

    ! Collisional friction force with impurities

    ELM(:,:, 8,NEQ) = - fem_int(20+L6,lab(:,i,3,2,1),Var(:,i)%n)
    NLC( 8,NEQ) = LQz3

    ! Collisional friction force with beam ions (No counterpart; LQb6)
    ! It may be retained because it is required to yield an adequate shielding factor.

    ELM(:,:, 9,NEQ) = - fem_int(20+L6,laf(:,i,2,1),  Ratden(:,i,4))
    NLC( 9,NEQ) = LQb3

    ! Collisional heat friction force with electrons

    ELM(:,:,10,NEQ) =   fem_int(20+L6,lab(:,i,1,2,2),Var(:,i)%n)
    NLC(10,NEQ) = LQe6

    ! Collisional heat friction force with ions (self collision)

    ELM(:,:,11,NEQ) =   fem_int(20+L6,lab(:,i,2,2,2),Var(:,i)%n)
    NLC(11,NEQ) = LQi6

    ! Collisional heat friction force with impurities

    ELM(:,:,12,NEQ) =   fem_int(20+L6,lab(:,i,3,2,2),Var(:,i)%n)
    NLC(12,NEQ) = LQz6

    ! Loss to divertor plasma region

    ELM(:,:,13,NEQ) = - fem_int(20+L6,rNuL,Var(:,i)%n) * fact
    NLC(13,NEQ) = LQi6

    ! Collisional friction force with neutrals

    ELM(:,:,14,NEQ) = - fem_int(20+L6,rNu0s(:,i),Var(:,i)%n)
    NLC(14,NEQ) = LQi6

!    ! Charge exchange force
!
!    ELM(:,:,15,NEQ) = - fem_int(20+L6,rNuiCX,Var(:,i)%n)
!    NLC(15,NEQ) = LQi6

!    ! Loss cone loss
!
!    ELM(:,:,16,NEQ) =   fem_int(-(2+L6),SiLCB,Var(:,i)%n)
!    NLC(16,NEQ) = 0

!    ! Ion orbit loss
!
!    ELM(:,:,17,NEQ) = - fem_int(20+L6,rNuOL,Var(:,i)%n)
!    NLC(17,NEQ) = LQi6

    NLCMAX(NEQ) = 14

    if( MDFIXT /= 0 ) then ! Overwrite Pressure gradient force
       ELM(:,:, 4,NEQ) = - rKilo / achg(i) * amasinv * 1.D-20 * fem_int(-(22+L6),fipolsdt,xmu(:,i,2,2),Var(:,i)%T)
       NLC( 4,NEQ) = 0

       ELM(:,:, 5,NEQ) =   0.d0
       NLC( 5,NEQ) = 0
    end if

  end subroutine LQi6CC

!***************************************************************
!
!   Ion rotation frequency: Ni <UiPhi/R>
!
!***************************************************************

  subroutine LQi7CC

    integer(4) :: NEQ = LQi7, i = 2

    ! Ion rotation frequency

    ELM(:,:,1,NEQ) = - fem_int(1)
    NLC(1,NEQ) = LQi7

    ! Ion rotation

    ELM(:,:,2,NEQ) =   fem_int(2,rrtinv)
    NLC(2,NEQ) = LQi4

    ELM(:,:,3,NEQ) = - fem_int(20,Fqhatsq,rrtinv)
    NLC(3,NEQ) = LQi4

    ! Ion parallel velocity

    ELM(:,:,4,NEQ) =   fem_int(30,Fqhatsq,fipolinv,Var(:,i)%n)
    NLC(4,NEQ) = LQi3

    NLCMAX(NEQ) = 4

  end subroutine LQi7CC

!***************************************************************
!
!   Ion diamagnetic flow: BV1
!
!***************************************************************

  subroutine LQi8CC

    integer(4) :: NEQ = LQi8, i = 2

    ! Diamagnetic flow

    ELM(:,:,1,NEQ) = - fem_int(2+L8,Var(:,i)%n)
    NLC(1,NEQ) = LQi8

    ! Pressure gradient

    ELM(:,:,2,NEQ) = - rKilo / achg(i) * fem_int(5+L8,fipolsdt)
    NLC(2,NEQ) = LQi5

    ! Electrostatic potential gradient

    ELM(:,:,3,NEQ) = - fem_int(22+L8,fipolsdt,Var(:,i)%n)
    NLC(3,NEQ) = LQm1

    NLCMAX(NEQ) = 3

    if( MDFIXT /= 0 ) then ! Overwrite Pressure gradient force
       ELM(:,:,2,NEQ) = - rKilo / achg(i) &
            &                  * (  fem_int(21+L8,fipolsdt,Var(:,i)%T) &
            &                     + fem_int(22+L8,fipolsdt,Var(:,i)%T))
       NLC(2,NEQ) = LQi1
    end if

  end subroutine LQi8CC

!***************************************************************
!
!   Impurity Density Equation
!
!***************************************************************

  subroutine LQz1CC

    integer(4) :: NEQ = LQz1, i = 3

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQz1

    ! Divergence

    ELM(:,:,1,NEQ) = - fem_int(4)
    NLC(1,NEQ) = LQz2

    ! Divergence + Grid velocity

    ELM(:,:,2,NEQ) =   fem_int(5,UgV)
    NLC(2,NEQ) = LQz1

    ! Ionization by effective ionization

    ELM(:,:,3,NEQ) =   FSION * 1.D20 * fem_int(20,SiVsefA,Var(:,1)%n)
    NLC(3,NEQ) = LQnz

    ! Generation of thermal neutrals by effective recombination

    ELM(:,:,4,NEQ) = - FSCX  * 1.D20 * fem_int(20,SiVa6A, Var(:,1)%n)
    NLC(4,NEQ) = LQz1

    ! Loss to divertor plasma region

    ELM(:,:,5,NEQ) = -             fem_int( 2,rNuL)
    NLC(5,NEQ) = LQz1

    ELM(:,:,6,NEQ) =   PNsSUP(i) * fem_int(-1,rNuL)
    NLC(6,NEQ) = 0

    ! Loss cone loss

    ELM(:,:,7,NEQ) =   fem_int(-1,SiLC)
    NLC(7,NEQ) = 0
 
    ! Ion orbit loss

    ELM(:,:,8,NEQ) = - fem_int(2,rNuOL)
    NLC(8,NEQ) = LQz1

    ! Diagnostic only (Dnfac /= 0): conservative particle diffusion with the
    ! same operator and geometry the momentum equation uses.  This asks which
    ! branch the zero-diffusion equations select as Dnfac -> 0+; it is not a
    ! model term and is off by default.

    if(Dnfac /= 0.d0) then
       ELM(:,:,9,NEQ) = - Dnfac * fem_int(27,sst,rMus(:,i))
       NLC(9,NEQ) = LQz1
       NLCMAX(NEQ) = 9
    else
       NLCMAX(NEQ) = 8
    end if

  end subroutine LQz1CC

!***************************************************************
!
!   Impurity Radial Flow
!
!***************************************************************
  
  subroutine LQz2CC

    integer(4) :: NEQ = LQz2, i = 3
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(i) * amp)

    ! Diamagnetic force

    ELM(:,:,1,NEQ) =   achg(i) * aee * amasinv * fem_int(30,bri,Var(:,i)%n,fipolinv)
    NLC(1,NEQ) = LQz8

    ! v x B force

    ELM(:,:,2,NEQ) =   achg(i) * aee * amasinv * fem_int(20,fipol,Var(:,i)%n)
    NLC(2,NEQ) = LQz3

    ELM(:,:,3,NEQ) = - achg(i) * aee * amasinv * fem_int( 2,bbt)
    NLC(3,NEQ) = LQz4

    NLCMAX(NEQ) = 3

  end subroutine LQz2CC

!***************************************************************
!
!   Impurity Parallel Flow
!
!***************************************************************

  subroutine LQz3CC

    integer(4) :: NEQ = LQz3, i = 3
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(i) * amp)

    ELM(:,:, 0,NEQ) =   fem_int(2+L3,Var(:,i)%n) * invDT
    NLC( 0,NEQ) = LQz3

    ! Neoclassical viscosity force

    ELM(:,:, 1,NEQ) = - amasinv * 1.D-20 * fem_int(2+L3,xmu(:,i,1,1))
    NLC( 1,NEQ) = LQz3

    ELM(:,:, 2,NEQ) = - amasinv * 1.D-20 * fem_int(2+L3,xmu(:,i,1,2))
    NLC( 2,NEQ) = LQz6

    ! Diamagnetic forces (particle)

    ELM(:,:, 3,NEQ) =   amasinv * 1.D-20 * fem_int(2+L3,xmu(:,i,1,1))
    NLC( 3,NEQ) = LQz8

    ! Diamagnetic forces (heat)

    ELM(:,:, 4,NEQ) = - rKilo / achg(i) * amasinv * 1.D-20 * fem_int(22+L3,fipolsdtPNsV(:,i),xmu(:,i,1,2))
    NLC( 4,NEQ) = LQz5

    ELM(:,:, 5,NEQ) =   rKilo / achg(i) * amasinv * 1.D-20 * fem_int(32+L3,fipolsdtPNsV(:,i),Var(:,i)%T,xmu(:,i,1,2))
    NLC( 5,NEQ) = LQz1
!    ELM(:,:, 5,NEQ) =   rKilo / achg(i) * amasinv * 1.D-20 * fem_int(30+L3,fipolsdtPNsV(:,i),dlnPNsV(:,i),xmu(:,i,1,2))
!    NLC( 5,NEQ) = LQz5

    ! Collisional friction force with electrons

    ELM(:,:, 6,NEQ) =   fem_int(20+L3,lab(:,i,1,1,1),Var(:,i)%n)
    NLC( 6,NEQ) = LQe3

    ! Collisional friction force with ions

    ELM(:,:, 7,NEQ) =   fem_int(20+L3,lab(:,i,2,1,1),Var(:,i)%n)
    NLC( 7,NEQ) = LQi3

    ! Collisional friction force with impurities (self collision)

    ELM(:,:, 8,NEQ) =   fem_int(20+L3,lab(:,i,3,1,1),Var(:,i)%n)
    NLC( 8,NEQ) = LQz3

    ! Collisional friction force with beam ions

    ELM(:,:, 9,NEQ) =   fem_int(20+L3,laf(:,i,1,1),  Ratden(:,i,4))
    NLC( 9,NEQ) = LQb3

    ! Collisional heat friction force with electrons

    ELM(:,:,10,NEQ) = - fem_int(20+L3,lab(:,i,1,1,2),Var(:,i)%n)
    NLC(10,NEQ) = LQe6

    ! Collisional heat friction force with ions

    ELM(:,:,11,NEQ) = - fem_int(20+L3,lab(:,i,2,1,2),Var(:,i)%n)
    NLC(11,NEQ) = LQi6

    ! Collisional heat friction force with impurities (self collision)

    ELM(:,:,12,NEQ) = - fem_int(20+L3,lab(:,i,3,1,2),Var(:,i)%n)
    NLC(12,NEQ) = LQz6

    ! Electric field forces

    ELM(:,:,13,NEQ) = - achg(i) * aee * amasinv * fem_int(30+L3,fipol,Var(:,i)%n,aatq)
    NLC(13,NEQ) = LQm2

    ELM(:,:,14,NEQ) =   achg(i) * aee * amasinv * fem_int(30+L3,fipol,Var(:,i)%n,aat)
    NLC(14,NEQ) = LQm3

    ! Loss to divertor plasma region

    ELM(:,:,15,NEQ) = - fem_int(20+L3,rNuL,Var(:,i)%n) * fact
    NLC(15,NEQ) = LQz3

    ! Collisional friction force with neutrals

    ELM(:,:,16,NEQ) = - fem_int(20+L3,rNu0s(:,i),Var(:,i)%n)
    NLC(16,NEQ) = LQz3

    ! Loss cone loss

    ELM(:,:,17,NEQ) =   fem_int(-(2+L3),SiLCB,Var(:,i)%n)
    NLC(17,NEQ) = 0

    ! Ion orbit loss

    ELM(:,:,18,NEQ) = - fem_int(20+L3,rNuOL,Var(:,i)%n)
    NLC(18,NEQ) = LQz3

    !---

    ! Perpendicular viscosity force: diffusion

    ELM(:,:,19,NEQ) = - FSPARV(i) * fem_int(38,sst*rMus(:,i)*Var(:,i)%n, ZetatoPara,ParatoZeta)
    NLC(19,NEQ) = LQz3

    ELM(:,:,20,NEQ) = - FSPARV(i) * fem_int(34,sst*rMus(:,i)*Var(:,i)%n,dZetatoPara,ParatoZeta)
    NLC(20,NEQ) = LQz3

    ! Perpendicular viscosity force: pinch

    ELM(:,:,21,NEQ) = - FSPARV(i) * fem_int(24,suft*Vmps(:,i),Var(:,i)%n)
    NLC(21,NEQ) = LQz3

    ELM(:,:,22,NEQ) =   FSPARV(i) * fem_int(30,suft*Vmps(:,i)*Var(:,i)%n,dZetatoPara,ParatoZeta)
    NLC(22,NEQ) = LQz3

    ! Perpendicular viscosity force: residual stress

    ELM(:,:,23,NEQ) = - amasinv * fem_int(-5,ZetatoPara,PiRess(:,i))
    NLC(23,NEQ) = 0

    ! Perpendicular viscosity force: terms associated with diamagnetic flow

    ELM(:,:,24,NEQ) = - FSPARV(i) * fem_int(-2,ZetatoPara,dtormflux(:,i))
    NLC(24,NEQ) = 0

!!$    ! Alternative implementation of perpendicular viscosity force: terms associated with diamagnetic flow
!!$
!!$    ELM(:,:,24,NEQ) = - FSPARV(i) * fem_int(38,sst *rMus(:,i)*Var(:,i)%n, ZetatoPara,ribi)
!!$    NLC(24,NEQ) = LQz8
!!$
!!$    ELM(:,:,25,NEQ) = - FSPARV(i) * fem_int(34,sst *rMus(:,i)*Var(:,i)%n,dZetatoPara,ribi)
!!$    NLC(25,NEQ) = LQz8
!!$
!!$    ELM(:,:,26,NEQ) = - FSPARV(i) * fem_int(24,suft*Vmps(:,i)*Var(:,i)%n, ZetatoPara*ribi)
!!$    NLC(26,NEQ) = LQz8
!!$
!!$    ELM(:,:,27,NEQ) =   FSPARV(i) * fem_int(30,suft*Vmps(:,i)*Var(:,i)%n,dZetatoPara,ribi)
!!$    NLC(27,NEQ) = LQz8

    !---

!!$    ! Perpendicular viscosity (off-diagonal)
!!$
!!$    ELM(:,:,19,NEQ) = - FSPARV(i) * fem_int(37,ZetatoPara,sst,rMus(:,i))
!!$    NLC(19,NEQ) = LQz4
!!$
!!$    ELM(:,:,20,NEQ) =   FSPARV(i) * fem_int(37,ZetatoPara,sst*rMus(:,i),Var(:,i)%RUph)
!!$    NLC(20,NEQ) = LQz1
!!$  
!!$    ELM(:,:,21,NEQ) = - FSPARV(i) * fem_int(32,dZetatoPara,sst,rMus(:,i))
!!$    NLC(21,NEQ) = LQz4
!!$
!!$    ELM(:,:,22,NEQ) =   FSPARV(i) * fem_int(32,dZetatoPara,sst*rMus(:,i),Var(:,i)%RUph)
!!$    NLC(22,NEQ) = LQz1
!!$
!!$    ! Momentum pinch (off-diagonal)
!!$
!!$    ELM(:,:,23,NEQ) = - FSPARV(i) * fem_int(24,suft*Vmps(:,i),ZetatoPara)
!!$    NLC(23,NEQ) = LQz4
!!$
!!$    ELM(:,:,24,NEQ) =   FSPARV(i) * fem_int(30,suft,Vmps(:,i),dZetatoPara)
!!$    NLC(24,NEQ) = LQz4

    NLCMAX(NEQ) = 24

    if( MDFIXT /= 0 ) then ! Overwrite Pressure gradient force
       ELM(:,:, 4,NEQ) = - rKilo / achg(i) * amasinv * 1.D-20 * fem_int(-(22+L3),fipolsdt,xmu(:,i,1,2),Var(:,i)%T)
       NLC( 4,NEQ) = 0

       ELM(:,:, 5,NEQ) =   0.d0
       NLC( 5,NEQ) = 0
    end if

  end subroutine LQz3CC

!***************************************************************
!
!   Impurity Toroidal Flow
!
!***************************************************************

  subroutine LQz4CC

    integer(4) :: NEQ = LQz4, i = 3, j, NCcomp
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(i) * amp)

    ! Uiphi'(0) : 0

    ELM(:,:, 0,NEQ) = fem_int(1) * invDT
    NLC( 0,NEQ) = LQz4

    ! Advection + Grid velocity

    ELM(:,:, 1,NEQ) = - fem_int(5,UsrgV(:,i))
    NLC( 1,NEQ) = LQz4

    ELM(:,:, 2,NEQ) = - FSADV * fem_int(6,Var(:,i)%UrV)
    NLC( 2,NEQ) = LQz4

    ! Momentum pinch

    ELM(:,:, 3,NEQ) = - fem_int(24,suft,Vmps(:,i))
    NLC( 3,NEQ) = LQz4

    ! Viscosity force

    ELM(:,:, 4,NEQ) = - fem_int(27,sst,rMus(:,i))
    NLC( 4,NEQ) = LQz4

    ELM(:,:, 5,NEQ) =   fem_int(37,sst,rMus(:,i),Var(:,i)%RUph)
    NLC( 5,NEQ) = LQz1
!!$    ELM(:,:, 5,NEQ) = - fem_int(24,sst*rMus(:,i),dlnPNsV(:,i))
!!$    NLC( 5,NEQ) = LQz4

    ! Residual stress

    ELM(:,:, 6,NEQ) = - amasinv * fem_int(-4,PiRess(:,i))
    NLC( 6,NEQ) = 0

    ! Collisional friction force with electrons

    ELM(:,:, 7,NEQ) =   fem_int(20,lab(:,i,1,1,1),Ratden(:,i,1))
    NLC( 7,NEQ) = LQe4

    ! Collisional friction force with bulk ions

    ELM(:,:, 8,NEQ) =   fem_int(20,lab(:,i,2,1,1),Ratden(:,i,2))
    NLC( 8,NEQ) = LQi4

    ! Collisional friction force with impurities (self collision)

    ELM(:,:, 9,NEQ) =   fem_int( 2,lab(:,i,3,1,1))
    NLC( 9,NEQ) = LQz4

    ! Collisional friction force with beam ions

!    ELM(:,:,10,NEQ) =   fem_int(30,ParatoZeta,laf(:,i,1,1),Ratden(:,i,4))
!    NLC(10,NEQ) = LQb3
    ELM(:,:,10,NEQ) =   fem_int(20,laf(:,i,1,1),Ratden(:,i,4))
    NLC(10,NEQ) = LQb4

    ! Collisional heat friction force with electrons

    j = 1
    ELM(:,:,11,NEQ) =   rKilo / achg(j) * fem_int(-22,lab(:,i,j,1,2),ribsdtNs(:,i),Var(:,j)%T)
    NLC(11,NEQ) = 0
!!$    ELM(:,:,11,NEQ) =   rKilo / achg(j) * fem_int(31,lab(:,i,j,1,2),ribsdt,Var(:,j)%T)
!!$    NLC(11,NEQ) = LQz1

    ELM(:,:,12,NEQ) = -                   fem_int(30,lab(:,i,j,1,2),ParatoZeta,Var(:,i)%n)
    NLC(12,NEQ) = LQe6

    ! Collisional heat friction force with ions

    j = 2
    ELM(:,:,13,NEQ) =   rKilo / achg(j) * fem_int(-22,lab(:,i,j,1,2),ribsdtNs(:,i),Var(:,j)%T)
    NLC(13,NEQ) = 0
!!$    ELM(:,:,13,NEQ) =   rKilo / achg(j) * fem_int(31,lab(:,i,j,1,2),ribsdt,Var(:,j)%T)
!!$    NLC(13,NEQ) = LQz1

    ELM(:,:,14,NEQ) = -                   fem_int(30,lab(:,i,j,1,2),ParatoZeta,Var(:,i)%n)
    NLC(14,NEQ) = LQi6

    ! Collisional heat friction force with impurities (self collision)

    j = 3
    ELM(:,:,15,NEQ) =   rKilo / achg(j) * fem_int(-22,lab(:,i,j,1,2),ribsdtNs(:,i),Var(:,j)%T)
    NLC(15,NEQ) = 0
!!$    ELM(:,:,15,NEQ) =   rKilo / achg(j) * fem_int(31,lab(:,i,j,1,2),ribsdt,Var(:,j)%T)
!!$    NLC(15,NEQ) = LQz1

    ELM(:,:,16,NEQ) = -                   fem_int(30,lab(:,i,j,1,2),ParatoZeta,Var(:,i)%n)
    NLC(16,NEQ) = LQz6

    ! Toroidal E force

    ELM(:,:,17,NEQ) =   achg(i) * aee * amasinv * fem_int(2,Var(:,i)%n)
    NLC(17,NEQ) = LQm3

    ! v x B force

    ELM(:,:,18,NEQ) =   achg(i) * aee * amasinv * fem_int(6,PsiV)
    NLC(18,NEQ) = LQz2

    ! Turbulent particle transport driver (electron driven)

    if( FSTPTM(i) == 0 ) then
       ELM(:,:,19,NEQ) =   achg(i) * rKeV * amasinv * fem_int(32,FQLcoef,Var(:,1)%T,Ratden(:,i,1)) &
            &                           * (1.d0 - FSVAHL)
       NLC(19,NEQ) = LQe1

       ELM(:,:,20,NEQ) =   achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Ratden(:,i,1)) &
            &                           * FSVAHL
       NLC(20,NEQ) = LQe5

    else if( FSTPTM(i) == 1 ) then
       ELM(:,:,19,NEQ) =   achg(i) * aee * amasinv * fem_int(30,FQLcoef2,ZetatoPara,Ratden(:,i,1))
       NLC(19,NEQ) = LQe4

       ELM(:,:,20,NEQ) = - achg(i) * aee * amasinv * fem_int(20,FQLcoef2,Var(:,i)%n)
       NLC(20,NEQ) = LQe3

       ELM(:,:,21,NEQ) =   achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Ratden(:,i,1)) &
            &                           * (FSVAHL - 1.d0)
       NLC(21,NEQ) = LQe5

       ELM(:,:,22,NEQ) = - achg(i) * rKeV * amasinv * fem_int(32,FQLcoef,Var(:,1)%T,Ratden(:,i,1)) &
            &                           * (FSVAHL - 1.d0)
       NLC(22,NEQ) = LQe1

       ELM(:,:,23,NEQ) =   achg(i) * aee * amasinv * fem_int(22,FQLcoef,Var(:,i)%n)
       NLC(23,NEQ) = LQm1

    else
       ELM(:,:,19,NEQ) =   achg(i) * aee * amasinv * fem_int(20,FQLcoef1,Var(:,i)%n)
       NLC(19,NEQ) = LQe8

       ELM(:,:,20,NEQ) =   achg(i) * aee * amasinv * fem_int(22,FQLcoef,Var(:,i)%n)
       NLC(20,NEQ) = LQm1

       ELM(:,:,21,NEQ) =   achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Ratden(:,i,1)) &
            &                           * (FSVAHL - 1.d0)
       NLC(21,NEQ) = LQe5

       ELM(:,:,22,NEQ) = - achg(i) * rKeV * amasinv * fem_int(32,FQLcoef,Var(:,1)%T,Ratden(:,i,1)) &
            &                           * (FSVAHL - 1.d0)
       NLC(22,NEQ) = LQe1

    end if

    !  -- Turbulent pinch term

    ELM(:,:,23,NEQ) =   achg(i) * aee * amasinv * fem_int(20,FQLcoef,FVpch)
    NLC(23,NEQ) = LQz1

    ! Loss to divertor plasma region

    ELM(:,:,24,NEQ) = - fem_int(2,rNuL) * fact
    NLC(24,NEQ) = LQz4

    ! Collisional friction force with neutrals

    ELM(:,:,25,NEQ) = - fem_int(2,rNu0s(:,i))
    NLC(25,NEQ) = LQz4

    ! Loss cone loss

    ELM(:,:,26,NEQ) =   fem_int(-1,SiLCph)
    NLC(26,NEQ) = 0

    ! Ion orbit loss

    ELM(:,:,27,NEQ) = - fem_int(2,rNuOL)
    NLC(27,NEQ) = LQz4

    NLCMAX(NEQ) = 27

    ! Multi-ion composition channel (FSTPTM = 3).
    !
    !   Delta K_a = sum_{b /= a, thermal ions} L_ab ( X_a - X_b ),
    !   X_a = (1/n_a) dn_a/dV,
    !   L_ab = L_ba = e (dpsi/dV) <|grad V|^2> D_e^QL (Z_a n_a)(Z_b n_b)/n_e.
    !
    ! Each pair contributes equal and opposite forces, so sum_a Delta K_a = 0
    ! holds element matrix by element matrix over whatever set of thermal ions
    ! is present: no net toroidal momentum is added, and since
    ! Gamma_a = -<R K_a1zeta>/(Z_a e dpsi/dV) that is the same statement as
    ! sum_a Z_a Delta Gamma_a = 0.  Beam ions are left out - they are too fast
    ! to be mixed by the background fluctuations - and their absence from the
    ! sum does not spoil the cancellation.
    !
    ! The 1/n_a of X_a cancels against the Z_a n_a of L_ab, so nothing here
    ! divides by a trace density.  With a single thermal ion the sum is empty
    ! and FSTPTM = 3 reduces exactly to FSTPTM = 2.
    if( FSTPTM(i) == 3.d0 ) then
       NCcomp = NLCMAX(NEQ)
       do j = 2, NSM
          if( j == i ) cycle
          NCcomp = NCcomp + 1
          ELM(:,:,NCcomp,NEQ) =   achg(i) * achg(j) * aee * amasinv &
               &                * fem_int(22,FQCcoef,Ratden(:,j,1))
          NLC(NCcomp,NEQ) = LQsn(i)
          NCcomp = NCcomp + 1
          ELM(:,:,NCcomp,NEQ) = - achg(i) * achg(j) * aee * amasinv &
               &                * fem_int(22,FQCcoef,Ratden(:,i,1))
          NLC(NCcomp,NEQ) = LQsn(j)
       end do
       NLCMAX(NEQ) = NCcomp
    end if


    if( MDFIXT /= 0 ) then ! Overwrite Pressure gradient force
       if( FSTPTM(i) == 0 ) then
          ELM(:,:,19,NEQ) =   achg(i) * rKeV * amasinv * fem_int(32,FQLcoef,Var(:,1)%T,Ratden(:,i,1))
          NLC(19,NEQ) = LQe1

          ELM(:,:,20,NEQ) =   achg(i) * rKeV * amasinv * fem_int(22,FQLcoef,Var(:,i)%n) &
               &                           * FSVAHL
          NLC(20,NEQ) = LQe5

       else
          ELM(:,:,21,NEQ) =   achg(i) * rKeV * amasinv * fem_int(21,FQLcoef,Var(:,1)%T) &
               &                           * (FSVAHL - 1.d0)
          NLC(21,NEQ) = LQz1

          ELM(:,:,22,NEQ) = 0.d0
          NLC(22,NEQ) = 0

       end if
    end if

  end subroutine LQz4CC

!***************************************************************
!
!  Impurity Energy Transport: Tz
!
!***************************************************************

  subroutine LQz5CC

    integer(4) :: NEQ = LQz5, i = 3

    ! Temperature evolution

    if( MDFIXT == 0 ) then
       ELM(:,:, 0,NEQ) =   1.5d0 * fem_int(1) * invDT
       NLC( 0,NEQ) = LQz5

       ! Advection + Grid velocity

       ELM(:,:, 1,NEQ) = - 2.5d0 * fem_int(5,UsrgV(:,i))
       NLC( 1,NEQ) = LQz5

       ELM(:,:, 2,NEQ) = - 2.5d0 * FSADV * fem_int(6,Var(:,i)%UrV)
       NLC( 2,NEQ) = LQz5

       ! Heat pinch

       ELM(:,:, 3,NEQ) = - fem_int(24,Vhps(:,i),suft)
       NLC( 3,NEQ) = LQz5

       ! Conduction transport

       ELM(:,:, 4,NEQ) = - fem_int(27,sst,zChi(:,i,1))
       NLC( 4,NEQ) = LQz5

       ELM(:,:, 5,NEQ) =   fem_int(37,sst,zChi(:,i,2),Var(:,i)%T)
       NLC( 5,NEQ) = LQz1

       ! Redundant heat convection term

       ELM(:,:, 6,NEQ) =   fem_int(5,UsrgV(:,i))
       NLC( 6,NEQ) = LQz5

       ! Viscous heating

       ELM(:,:, 7,NEQ) =   1.D-20 / rKeV * fem_int(-2,Var(:,i)%Uthhat,BnablaPi(:,i))
       NLC( 7,NEQ) = 0

       ! Collisional transfer with electrons (Energy equilibration)

       ELM(:,:, 8,NEQ) = - 1.5d0 * fem_int(20,rNuTez,Ratden(:,1,i))
       NLC( 8,NEQ) = LQz5

       ELM(:,:, 9,NEQ) =   1.5d0 * fem_int( 2,rNuTez)
       NLC( 9,NEQ) = LQe5

       ! Collisional transfer with ions (Energy equilibration)

       ELM(:,:,10,NEQ) = - 1.5d0 * fem_int(20,rNuTiz,Ratden(:,2,i))
       NLC(10,NEQ) = LQz5

       ELM(:,:,11,NEQ) =   1.5d0 * fem_int( 2,rNuTiz)
       NLC(11,NEQ) = LQi5

       ! Loss to divertor plasma region

       ELM(:,:,12,NEQ) = -             fem_int( 2,rNuL)
       NLC(12,NEQ) = LQz5

       ELM(:,:,13,NEQ) =   PNsSUP(i) * fem_int(-2,rNuL,Var(:,i)%T)
       NLC(13,NEQ) = 0

       ELM(:,:,14,NEQ) = - 1.5d0            * fem_int( 2,rNuLTs(:,i))
       NLC(14,NEQ) = LQz5

       ELM(:,:,15,NEQ) =   1.5d0 * PTsSUP(i) * fem_int(-2,rNuLTs(:,i),Var(:,i)%n)
       NLC(15,NEQ) = 0
!       ELM(:,:,15,NEQ) =   1.5d0 * PTsSUP(i)    * fem_int( 2,rNuLTs(:,i))
!       NLC(15,NEQ) = LQz1

       ! Simplified Alpha heating

       ELM(:,:,16,NEQ) = 1.D-20 / rKeV * fem_int(-1,PALFz)
       NLC(16,NEQ) = 0

       ! Direct heating (RF)

       ELM(:,:,17,NEQ) = 1.D-20 / rKeV * fem_int(-1,PRFz)
       NLC(17,NEQ) = 0

       NLCMAX(NEQ) = 17
    else

       ! Fixed temperature profile

       ELM(:,:,0,NEQ) = fem_int(1) * invDT
       NLC(0,NEQ) = LQz5

       NLCMAX(NEQ) = 0
    end if

  end subroutine LQz5CC

!***************************************************************
!
!   Impurity Heat Parallel Flow
!
!***************************************************************

  subroutine LQz6CC

    integer(4) :: NEQ = LQz6, i = 3
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(i) * amp)

    ELM(:,:, 0,NEQ) =   2.5d0 * fem_int(2+L6,Var(:,i)%n) * invDT
    NLC( 0,NEQ) = LQz6

    ! Neoclassical viscosity force

    ELM(:,:, 1,NEQ) = - amasinv * 1.D-20 * fem_int(2+L6,xmu(:,i,2,1))
    NLC( 1,NEQ) = LQz3

    ELM(:,:, 2,NEQ) = - amasinv * 1.D-20 * fem_int(2+L6,xmu(:,i,2,2))
    NLC( 2,NEQ) = LQz6

    ! Diamagnetic forces

    ELM(:,:, 3,NEQ) =   amasinv * 1.D-20 * fem_int(2+L6,xmu(:,i,2,1))
    NLC( 3,NEQ) = LQz8

    ELM(:,:, 4,NEQ) = - rKilo / achg(i) * amasinv * 1.D-20 * fem_int(22+L6,fipolsdtPNsV(:,i),xmu(:,i,2,2))
    NLC( 4,NEQ) = LQz5

    ELM(:,:, 5,NEQ) =   rKilo / achg(i) * amasinv * 1.D-20 * fem_int(32+L6,fipolsdtPNsV(:,i),Var(:,i)%T,xmu(:,i,2,2))
    NLC( 5,NEQ) = LQz1

    ! Collisional friction force with electrons

    ELM(:,:, 6,NEQ) = - fem_int(20+L6,lab(:,i,1,2,1),Var(:,i)%n)
    NLC( 6,NEQ) = LQe3

    ! Collisional friction force with ions

    ELM(:,:, 7,NEQ) = - fem_int(20+L6,lab(:,i,2,2,1),Var(:,i)%n)
    NLC( 7,NEQ) = LQi3

    ! Collisional friction force with impurities (self collision)

    ELM(:,:, 8,NEQ) = - fem_int(20+L6,lab(:,i,3,2,1),Var(:,i)%n)
    NLC( 8,NEQ) = LQz3

    ! Collisional friction force with beam ions (No counterpart; LQb6)
    ! It may be retained because it is required to yield an adequate shielding factor.

    ELM(:,:, 9,NEQ) = - fem_int(20+L6,laf(:,i,2,1),  Ratden(:,i,4))
    NLC( 9,NEQ) = LQb3

    ! Collisional heat friction force with electrons

    ELM(:,:,10,NEQ) =   fem_int(20+L6,lab(:,i,1,2,2),Var(:,i)%n)
    NLC(10,NEQ) = LQe6

    ! Collisional heat friction force with ions

    ELM(:,:,11,NEQ) =   fem_int(20+L6,lab(:,i,2,2,2),Var(:,i)%n)
    NLC(11,NEQ) = LQi6

    ! Collisional heat friction force with impurities (self collision)

    ELM(:,:,12,NEQ) =   fem_int(20+L6,lab(:,i,3,2,2),Var(:,i)%n)
    NLC(12,NEQ) = LQz6

    ! Loss to divertor plasma region

    ELM(:,:,13,NEQ) = - fem_int(20+L6,rNuL,Var(:,i)%n) * fact
    NLC(13,NEQ) = LQz6

    ! Collisional friction force with neutrals

    ELM(:,:,14,NEQ) = - fem_int(20+L6,rNu0s(:,i),Var(:,i)%n)
    NLC(14,NEQ) = LQz6

!    ! Loss cone loss
!
!    ELM(:,:,15,NEQ) =   fem_int(-(2+L6),SiLCB,Var(:,i)%n)
!    NLC(15,NEQ) = 0

!    ! Ion orbit loss
!
!    ELM(:,:,16,NEQ) = - fem_int(20+L6,rNuOL,Var(:,i)%n)
!    NLC(16,NEQ) = LQz6

    NLCMAX(NEQ) = 14

    if( MDFIXT /= 0 ) then ! Overwrite Pressure gradient force
       ELM(:,:, 4,NEQ) = - rKilo / achg(i) * amasinv * 1.D-20 * fem_int(-(22+L6),fipolsdt,xmu(:,i,2,2),Var(:,i)%T)
       NLC( 4,NEQ) = 0

       ELM(:,:, 5,NEQ) =   0.d0
       NLC( 5,NEQ) = 0
    end if

  end subroutine LQz6CC

!***************************************************************
!
!   Impurity rotation frequency: Nz <UzPhi/R>
!
!***************************************************************

  subroutine LQz7CC

    integer(4) :: NEQ = LQz7, i = 3

    ! Ion rotation frequency

    ELM(:,:,1,NEQ) = - fem_int(1)
    NLC(1,NEQ) = LQz7

    ! Ion rotation

    ELM(:,:,2,NEQ) =   fem_int(2,rrtinv)
    NLC(2,NEQ) = LQz4

    ELM(:,:,3,NEQ) = - fem_int(20,Fqhatsq,rrtinv)
    NLC(3,NEQ) = LQz4

    ! Ion parallel velocity

    ELM(:,:,4,NEQ) =   fem_int(30,Fqhatsq,fipolinv,Var(:,i)%n)
    NLC(4,NEQ) = LQz3

    NLCMAX(NEQ) = 4

  end subroutine LQz7CC

!***************************************************************
!
!   Impurity diamagnetic flow: BV1
!
!***************************************************************

  subroutine LQz8CC

    integer(4) :: NEQ = LQz8, i = 3

    ! Diamagnetic flow

    ELM(:,:,1,NEQ) = - fem_int(2+L8,Var(:,i)%n)
    NLC(1,NEQ) = LQz8

    ! Pressure gradient

    ELM(:,:,2,NEQ) = - rKilo / achg(i) * fem_int(5+L8,fipolsdt)
    NLC(2,NEQ) = LQz5

    ! Electrostatic potential gradient

    ELM(:,:,3,NEQ) = - fem_int(22+L8,fipolsdt,Var(:,i)%n)
    NLC(3,NEQ) = LQm1

    NLCMAX(NEQ) = 3

    if( MDFIXT /= 0 ) then ! Overwrite Pressure gradient force
       ELM(:,:,2,NEQ) = - rKilo / achg(i) &
            &                  * (  fem_int(21+L8,fipolsdt,Var(:,i)%T) &
            &                     + fem_int(22+L8,fipolsdt,Var(:,i)%T))
       NLC(2,NEQ) = LQz1
    end if

  end subroutine LQz8CC

!***************************************************************
!
!   Beam Ion Density
!
!***************************************************************

  subroutine LQb1CC

    integer(4) :: NEQ = LQb1

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQb1

    if( MDBEAM /= 1 ) then

       ! Grid velocity

       ELM(:,:,1,NEQ) = - fem_int(6,UgV)
       NLC(1,NEQ) = LQb1

    else

       ! Divergence

       ELM(:,:,1,NEQ) = - fem_int(4)
       NLC(1,NEQ) = LQb2

       ! Grid velocity

       ELM(:,:,2,NEQ) =   fem_int(5,UgV)
       NLC(2,NEQ) = LQb1

    end if

    ! NBI particle source (Both charge exchange and ionization)

    ELM(:,:,3,NEQ) =   fem_int(-1,SNBb)
    NLC(3,NEQ) = 0

    ! Extracted NBI perpendicular component fallen into the ripple well region
    ! (only when the ripple trapped ions are solved, FSRP = 1; LQr1 receives it)

    ELM(:,:,4,NEQ) = - fem_int(-2,SNBPDi,rip_rat) * RpplSW
    NLC(4,NEQ) = 0

    ! Relaxation to thermal ions

    ELM(:,:,5,NEQ) = - fem_int(2,rNuB)
    NLC(5,NEQ) = LQb1

    ! Loss to divertor plasma region

    ELM(:,:,6,NEQ) = - fem_int(2,rNuLB) * BeamSW
    NLC(6,NEQ) = LQb1

    ! Ripple trapped beam ions collision with otherwise beam ions

    ELM(:,:,7,NEQ) = - fem_int(20,rNubrp2,rip_rat) * RpplSW
    NLC(7,NEQ) = LQb1

    ELM(:,:,8,NEQ) =   fem_int(20,rNubrp1,rip_rat) * RpplSW
    NLC(8,NEQ) = LQr1

!!rp_conv    ELM(:,:,8,NEQ) =  fem_int(-2,PNbrpLV,rNubLL)
!!rp_conv    NLC(8,NEQ) = 0

    ! Ripple diffusion

    ELM(:,:,9,NEQ) = - fem_int(37,sst,Dbrp,ft)
    NLC(9,NEQ) = LQb1

    NLCMAX(NEQ) = 9

  end subroutine LQb1CC

!***************************************************************
!
!   Beam Ion Radial Flow
!
!***************************************************************
  
  subroutine LQb2CC

    integer(4) :: NEQ = LQb2
    real(8) :: amasinv

    amasinv = 1.d0 / (amb * amp)

    if( MDBEAM /= 1 ) then

       ! Dummy term to avoid null when not solving LQb2

       ELM(:,:,0,NEQ) = fem_int(1) * invDT
       NLC(0,NEQ) = LQb2

    else

       ! Diamagnetic force

       ELM(:,:,1,NEQ) =   achgb * aee * amasinv * fem_int(20,bri,fipolinv)
       NLC(1,NEQ) = LQb8

       ! v x B force

       ELM(:,:,2,NEQ) =   achgb * aee * amasinv * fem_int( 2,fipol)
       NLC(2,NEQ) = LQb3

       ELM(:,:,3,NEQ) = - achgb * aee * amasinv * fem_int( 2,bbt)
       NLC(3,NEQ) = LQb4

    end if

    NLCMAX(NEQ) = 3

  end subroutine LQb2CC

!***************************************************************
!
!   Beam Ion Parallel flow
!
!***************************************************************

  subroutine LQb3CC

    integer(4) :: NEQ = LQb3
    real(8) :: amasinv

    amasinv = 1.d0 / (amb * amp)

!    if( iSUPG3 > 0 ) L3 = 0

    ELM(:,:, 0,NEQ) = fem_int(1+L3) * invDT
    NLC( 0,NEQ) = LQb3

    if( MDBEAM /= 2 ) then

       ! Neoclassical viscosity

       ELM(:,:, 1,NEQ) = - amasinv * 1.d-20 * fem_int(20+L3,PNbVinv,xmuf(:,1))
       NLC( 1,NEQ) = LQb3

       ! Diamagnetic forces

       ELM(:,:, 2,NEQ) =   amasinv * 1.d-20 * fem_int(20+L3,PNbVinv,xmuf(:,1))
       NLC( 2,NEQ) = LQb8

       ! Collisional friction force with beam ions (self collision)

       ELM(:,:, 3,NEQ) =   (amas(1) / amb) * fem_int(20+L3,lff(:,1,1),Ratden(:,1,4))
       NLC( 3,NEQ) = LQb3

       ! Collisional friction force with electrons

       ELM(:,:, 4,NEQ) =   (amas(1) / amb) * fem_int(20+L3,lfb(:,1,1,1),Var(:,1)%n)
       NLC( 4,NEQ) = LQe3

       ! Collisional friction force with ions

       ELM(:,:, 5,NEQ) =   (amas(1) / amb) * fem_int(20+L3,lfb(:,2,1,1),Var(:,1)%n)
       NLC( 5,NEQ) = LQi3

       ! Collisional friction force with impurities

       ELM(:,:, 6,NEQ) =   (amas(1) / amb) * fem_int(20+L3,lfb(:,3,1,1),Var(:,1)%n)
       NLC( 6,NEQ) = LQz3

       ! Electric field forces

       ELM(:,:, 7,NEQ) = - achgb * aee * amasinv * fem_int(30+L3,PNbV,fipol,aatq)
       NLC( 7,NEQ) = LQm2

       ELM(:,:, 8,NEQ) =   achgb * aee * amasinv * fem_int(30+L3,PNbV,fipol,aat)
       NLC( 8,NEQ) = LQm3

       ! Loss to divertor plasma regino

       ELM(:,:, 9,NEQ) = - fem_int(2+L3,rNuLB) * fact! * BeamSW
       NLC( 9,NEQ) = LQb3

       ! Momentum loss due to thermalization

       ELM(:,:,10,NEQ) = - fem_int(2+L3,rNuB)! * BeamSW
       NLC(10,NEQ) = LQb3

       ! Collisional friction force with neutrals

       ELM(:,:,11,NEQ) = - fem_int(2+L3,rNu0b)! * BeamSW
       NLC(11,NEQ) = LQb3

       ! Charge exchange force

       ELM(:,:,12,NEQ) = - fem_int(2+L3,rNubCX)! * BeamSW
       NLC(12,NEQ) = LQb3

       ! NBI momentum source

       ELM(:,:,13,NEQ) =  amasinv * fem_int(-(1+L3),BSmb)
       NLC(13,NEQ) = 0

       ! Ripple trapping and ripple diffusion carry no beam momentum at the
       ! order kept here: the ions scattered into the ripple wells have
       ! |v_par|/v < sqrt(delta), and the banana ions that diffuse have zero
       ! bounce-averaged parallel and toroidal momentum up to finite orbit
       ! width, which this beam fluid (no beam viscosity) does not keep.
       ! Neither LQb3 nor LQb4 has a ripple term (see doc/tx561.tex).

    end if

    NLCMAX(NEQ) = 13

  end subroutine LQb3CC

!***************************************************************
!
!   Beam Ion Toroidal Flow
!
!***************************************************************

  subroutine LQb4CC

    integer(4) :: NEQ = LQb4
    real(8) :: amasinv

    amasinv = 1.d0 / (amb * amp)

    ! Ubphi'(0) : 0

    ELM(:,:, 0,NEQ) =   fem_int(1) * invDT
    NLC( 0,NEQ) = LQb4

    if( MDBEAM /= 2 ) then

       ! Advection + Grid velocity

       ELM(:,:, 1,NEQ) =   fem_int(5,UgV)
       NLC( 1,NEQ) = LQb4

       if( MDBEAM == 1 ) then

          ELM(:,:, 2,NEQ) = - FSADVB * fem_int(3,RUbphV)
          NLC( 2,NEQ) = LQb2

       end if

       ! Collisional friction force with electrons

       ELM(:,:, 3,NEQ) =   (amas(1) / amb) * fem_int(2,lfb(:,1,1,1))
       NLC( 3,NEQ) = LQe4

       ! Collisional friction force with bulk ions

       ELM(:,:, 4,NEQ) =   (amas(1) / amb) * fem_int(20,lfb(:,2,1,1),Ratden(:,1,2))
       NLC( 4,NEQ) = LQi4

       ! Collisional friction force with impurities

       ELM(:,:, 5,NEQ) =   (amas(1) / amb) * fem_int(20,lfb(:,3,1,1),Ratden(:,1,3))
       NLC( 5,NEQ) = LQz4

       ! Collisional friction force with beam ions (self collision)

       ELM(:,:, 6,NEQ) =   (amas(1) / amb) * fem_int(20,lff(:,1,1),Ratden(:,1,4))
       NLC( 6,NEQ) = LQb4

       ! Toroidal E force

       ELM(:,:, 7,NEQ) =   achgb * aee * amasinv * fem_int(2,PNbV)
       NLC( 7,NEQ) = LQm3

       if( MDBEAM == 1 ) then

          ! v x B force

          ELM(:,:, 8,NEQ) =   achgb * aee * amasinv * fem_int(6,PsiV)
          NLC( 8,NEQ) = LQb2

       end if

       ! Loss to divertor plasma region

       ELM(:,:, 9,NEQ) = - fem_int(2,rNuLB) * fact
       NLC( 9,NEQ) = LQb4

       ! Momentum loss due to thermalization

       ELM(:,:,10,NEQ) = - fem_int(2,rNuB)
       NLC(10,NEQ) = LQb4

       ! Collisional friction force with neutrals

       ELM(:,:,11,NEQ) = - fem_int(2,rNu0b)
       NLC(11,NEQ) = LQb4

       ! Charge exchange force

       ELM(:,:,12,NEQ) = - fem_int(2,rNubCX)
       NLC(12,NEQ) = LQb4

       ! NBI momentum source : (I / <B^2>) m_b dot{n}_b <B v_//0>

       ELM(:,:,13,NEQ) =   amasinv * fem_int(-2,ParatoZeta,BSmb)
       NLC(13,NEQ) = 0

!       ! Turbulent particle transport driver (electron driven)
!
!       ELM(:,:,14,NEQ) =   achgb * aee * amasinv * fem_int(30,FQLcoef2,ZetatoPara,RatPNbinv(:,1))
!       NLC(14,NEQ) = LQe4
!
!       ELM(:,:,15,NEQ) = - achgb * aee * amasinv * fem_int(20,FQLcoef2,PNbV)
!       NLC(15,NEQ) = LQe3
!
!       ELM(:,:,16,NEQ) =   achgb * rKeV * amasinv * fem_int(22,FQLcoef,RatPNbinv(:,1)) &
!            &                           * (FSVAHL - 1.d0)
!       NLC(16,NEQ) = LQe5
!
!       ELM(:,:,17,NEQ) = - achgb * rKeV * amasinv * fem_int(32,FQLcoef,Var(:,1)%T,RatPNbinv(:,1)) &
!            &                           * (FSVAHL - 1.d0)
!       NLC(17,NEQ) = LQe1
!
!       ELM(:,:,18,NEQ) =   achgb * aee * amasinv * fem_int(22,FQLcoef,PNbV)
!       NLC(18,NEQ) = LQm1
!
!       !  -- Turbulent pinch term
!
!       ELM(:,:,19,NEQ) =   achgb * aee * amasinv * fem_int(20,FQLcoef,FVpch)
!       NLC(19,NEQ) = LQb1

    end if

    NLCMAX(NEQ) = 13

  end subroutine LQb4CC

!***************************************************************
!
!  Beam Ion rotation frequency: Nb <UbPhi/R>
!
!***************************************************************

  subroutine LQb7CC

    integer(4) :: NEQ = LQb7

    ! Beam Ion rotation frequency

    ELM(:,:,1,NEQ) = - fem_int(1)
    NLC(1,NEQ) = LQb7

    ! Beam Ion rotation

    ELM(:,:,2,NEQ) =   fem_int(2,rrtinv)
    NLC(2,NEQ) = LQb4

    ELM(:,:,3,NEQ) = - fem_int(20,Fqhatsq,rrtinv)
    NLC(3,NEQ) = LQb4

    ! Beam Ion parallel velocity

    ELM(:,:,4,NEQ) =   fem_int(20,Fqhatsq,fipolinv)
    NLC(4,NEQ) = LQb3

    NLCMAX(NEQ) = 4

  end subroutine LQb7CC

!***************************************************************
!
!   Beam Ion diamagnetic flow: Nb BV1
!
!***************************************************************

  subroutine LQb8CC

    integer(4) :: NEQ = LQb8

    if( MDBEAM /= 1 ) then

       ! Dummy term to avoid null when not solving LQb2

       ELM(:,:,0,NEQ) = fem_int(1) * invDT
       NLC(0,NEQ) = LQb8

    else

       ! Diamagnetic flow

       ELM(:,:,1,NEQ) = - fem_int(1+L8)
       NLC(1,NEQ) = LQb8

       ! Pressure gradient

       ELM(:,:,2,NEQ) = - rKilo / achgb * fem_int(21+L8,fipolsdt,PTbV)
       NLC(2,NEQ) = LQb1

       ELM(:,:,3,NEQ) = - rKilo / achgb * fem_int(22+L8,fipolsdt,PTbV)
       NLC(3,NEQ) = LQb1

       ! Electrostatic potential gradient

       ELM(:,:,4,NEQ) = - fem_int(22+L8,fipolsdt,PNbV)
       NLC(4,NEQ) = LQm1

    end if

    NLCMAX(NEQ) = 4

  end subroutine LQb8CC

!***************************************************************
!
!   Slow Neutral Transport: n01
!
!***************************************************************

  subroutine LQn1CC

    integer(4) :: NEQ = LQn1

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQn1

    ! Grid velocity

    ELM(:,:,1,NEQ) = - fem_int(6,UgV)
    NLC(1,NEQ) = LQn1

    ! Diffusion of neutrals

    if( MDLN0FL == 2 ) then
       ! Larsen limiter evaluated element by element, as for the C0 part of the
       ! carbon bundle (LQnzCC): from the element mean nbar of n01 and the
       ! element difference dn of the current iterate,
       !   D_e = Dbar / (1 + Dbar |dn| sqrt(<|grad V|^2>) / (h FLlim vbar nbar)).
       ! |dn|/nbar <= 2 on an element, so the factor is bounded by
       ! 1 + 2 lambda/(3 FLlim h) whatever the local profile.  The limited
       ! ionization front is exponential with the decay length FLlim lambda
       ! (5 mm at rho ~ 0.87 with FLlim = 1/4), shorter than an element at
       ! NRMAX = 60; the nodal form of MDLN0FL = 1 evaluates |grad n|/n on that
       ! unresolved front and in the tiny-n01 tail inside it, where it is large
       ! and erratic, and the Picard iteration does not settle.  This form
       ! needs no relaxation.  D01 itself is left unlimited.
       block
         integer(4) :: NE
         real(8) :: V0aveL, Del, sstl, dn, an
         V0aveL = sqrt(4.d0 * V0*V0 / Pi)
         do NE = 1, NEMAX
            Del  = 0.5d0 * (D01(NE-1) + D01(NE))
            sstl = 0.5d0 * (sst(NE-1) + sst(NE))
            dn   = PN01V(NE) - PN01V(NE-1)
            an   = 0.5d0 * (PN01V(NE-1) + PN01V(NE))
            if( an > 0.d0 ) Del = Del / (1.d0 + Del * abs(dn) / hv(NE) * sqrt(sstl) &
                 &                                    / (FLlim * V0aveL * an))
            ELM(NE,1,2,NEQ) = - sstl * Del / hv(NE)
            ELM(NE,2,2,NEQ) =   sstl * Del / hv(NE)
            ELM(NE,3,2,NEQ) =   sstl * Del / hv(NE)
            ELM(NE,4,2,NEQ) = - sstl * Del / hv(NE)
         end do
       end block
    else
       ELM(:,:,2,NEQ) = - fem_int(27,sst,D01)
    end if
    NLC(2,NEQ) = LQn1

    ! Ionization by electron impact
    !   cf. [C.E. Singer et al., Comput Phys. Commun. 49 (1988) 275, (2.9.1a)]

    ELM(:,:,3,NEQ) = - FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n)
    NLC(3,NEQ) = LQn1

    ! Generation of thermal neutrals by charge exchange

    ELM(:,:,4,NEQ) = - FSCX * 1.D20 * fem_int(20,SiVcxA,Var(:,2)%n)
    NLC(4,NEQ) = LQn1

    ! Recycling in divertor plasma region

    ELM(:,:,5,NEQ) =   rGamm0             * fem_int( 2,rNuL)
    NLC(5,NEQ) = LQi1

    ELM(:,:,6,NEQ) = - rGamm0 * PNsSUP(2) * fem_int(-1,rNuL)
    NLC(6,NEQ) = 0

!    ! Gas puff (alternative way to impose B.C.;
!                Deactivate call boundary(lqn1) if this term is activated)
!
!    ELM(:,:,7,NEQ) = fem_int(-3,suft,rGASPFA)
!    NLC(7,NEQ) = 0

    NLCMAX(NEQ) = 6

  end subroutine LQn1CC

!***************************************************************
!
!   Thermal Neutral Transport: n02
!
!***************************************************************

  subroutine LQn2CC

    integer(4) :: NEQ = LQn2

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQn2

    ! Grid velocity

    ELM(:,:,1,NEQ) = - fem_int(6,UgV)
    NLC(1,NEQ) = LQn2

    ! Diffusion of neutrals

    ELM(:,:,2,NEQ) = - fem_int(27,sst,D02)
    NLC(2,NEQ) = LQn2

    ! Ionization by electron impact

    ELM(:,:,3,NEQ) = - FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n)
    NLC(3,NEQ) = LQn2

    ! Generation of thermal neutrals by charge exchange

    ELM(:,:,4,NEQ) = FSCX * 1.D20 * fem_int(20,SiVcxA,Var(:,2)%n)
    NLC(4,NEQ) = LQn1

!    ! NBI particle source (Charge exchange), replaced by LQn3CC
!
!!    ELM(:,:,5,NEQ) = RatCX * fem_int(-1,SNBi)
!    ELM(:,:,5,NEQ) = fem_int(-1,RatCXSNBi)
!    NLC(5,NEQ) = 0

    NLCMAX(NEQ) = 4

  end subroutine LQn2CC

!***************************************************************
!
!   Halo Neutral Transport: n03
!
!***************************************************************

  subroutine LQn3CC

    integer(4) :: NEQ = LQn3

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQn3

    ! Grid velocity

    ELM(:,:,1,NEQ) = - fem_int(6,UgV)
    NLC(1,NEQ) = LQn3

    ! Diffusion of neutrals

    ELM(:,:,2,NEQ) = - fem_int(27,sst,D03)
    NLC(2,NEQ) = LQn3

    ! Ionization by electron impact

    ELM(:,:,3,NEQ) = - FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n)
    NLC(3,NEQ) = LQn3

    ! NBI particle source (Charge exchange)

!    ELM(:,:,4,NEQ) = RatCX * fem_int(-1,SNBi)
    ELM(:,:,4,NEQ) = fem_int(-1,RatCXSNBi)
    NLC(4,NEQ) = 0

    NLCMAX(NEQ) = 4

  end subroutine LQn3CC

!***************************************************************
!
!   Superstaged impurity Neutral Transport: n0z
!
!***************************************************************

  subroutine LQnzCC

    integer(4) :: NEQ = LQnz, i = 3, j = 2
    integer(4) :: NE
    real(8) :: V0zaveL, Del, sstl, dnC0, anC0, Fel

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQnz

    ! Grid velocity

    ELM(:,:,1,NEQ) = - fem_int(6,UgV)
    NLC(1,NEQ) = LQnz

    ! Diffusion: Dnzb = D01 (MDLNZD = 0) or Dfs (MDLNZD = 1), set in TXCALC.

    ELM(:,:,2,NEQ) = - fem_int(27,sst,Dnzb)
    NLC(2,NEQ) = LQnz

    ! Effective ionization generating C^{6+}

    ELM(:,:,3,NEQ) = - FSION * 1.D20 * fem_int(20,SiVsefA,Var(:,1)%n)
    NLC(3,NEQ) = LQnz

    ! Effective recombination generating C^{0-5+}

    ELM(:,:,4,NEQ) =   FSCX  * 1.D20 * fem_int(20,SiVa6A, Var(:,1)%n)
    NLC(4,NEQ) = LQz1

    ! Recycling due to self sputtering in divertor plasma region

    ELM(:,:,5,NEQ) =   rGamm0zz             * fem_int( 2,rNuL)
    NLC(5,NEQ) = LQz1

    ELM(:,:,6,NEQ) = - rGamm0zz * PNsSUP(i) * fem_int(-1,rNuL)
    NLC(6,NEQ) = 0

    ! Recycling due to sputtering by bulk ions in divertor plasma region

    ELM(:,:,7,NEQ) =   rGamm0iz             * fem_int( 2,rNuL)
    NLC(7,NEQ) = LQi1

    ELM(:,:,8,NEQ) = - rGamm0iz * PNsSUP(j) * fem_int(-1,rNuL)
    NLC(8,NEQ) = 0

    NLCMAX(NEQ) = 8

    ! MDLNZD = 1: the C0 part of the bundle, n_C0, is a known field slaved
    ! to the sputtering source (see TXCALC).  It adds
    !   - the C0 flux in excess of the D_c term, -(D_C0 - D_c) dn_C0/dV, and
    !   - the parallel loss of the charged part n_b - n_C0 in the SOL, of which
    !     the fraction rGamm0zz comes back as sputtered C0 (as for C6+ above).
    !
    ! The C0 flux is assembled element by element, with the Larsen limiter
    ! evaluated on each element from the element mean of n_C0 and the element
    ! gradient: at a node, |grad n_C0|/n_C0 is 0/0 at the separatrix, where the
    ! source and hence n_C0 vanish, and whether rounding left rNuL(NRA) at 0 or
    ! at 1e-25 decided whether that node got the unlimited D_C0.  The limiter
    ! acts on the whole D_C0, since it bounds the whole C0 flux by free
    ! streaming; D_c, whose part of the C0 flux is already in the D_c term,
    ! is subtracted afterwards, with its sign: the bundle flux is then exactly
    ! -D_c grad(n_b - n_C0) - D_C0 grad n_C0.  Delta D < 0, which occurs in
    ! the first element outside the separatrix on fine meshes, where the
    ! limited D_C0 falls to about alpha vbar_0z h / 2, only removes the D_c
    ! part of the C0 flux; n_C0 is a known field, so the matrix is unchanged.  The element flux F_e enters the two rows of
    ! the element as +F_e and -F_e, so the term is exactly conservative.
    ! fem_int cannot be used here: its coefficients are nodal, while the
    ! limited diffusivity is defined per element.

    if( MDLNZD == 1 ) then
       V0zaveL = sqrt(4.d0 * V0z*V0z / Pi)
       ELM(:,:,9,NEQ) = 0.d0
       do NE = 1, NEMAX
          Del  = 0.5d0 * (DC0z(NE-1) + DC0z(NE))
          sstl = 0.5d0 * (sst(NE-1)  + sst(NE))
          dnC0 = nC0z(NE) - nC0z(NE-1)
          anC0 = 0.5d0 * (nC0z(NE-1) + nC0z(NE))
          if( anC0 > 0.d0 ) Del = Del / (1.d0 + Del * abs(dnC0) / hv(NE) * sqrt(sstl) &
               &                                    / (FLlim * V0zaveL * anC0))
          Del  = Del - 0.5d0 * (Dnzb(NE-1) + Dnzb(NE))
          Fel = sstl * Del * dnC0 / hv(NE)
          ELM(NE,1,9,NEQ) =   Fel
          ELM(NE,3,9,NEQ) = - Fel
       end do
       NLC(9,NEQ) = 0

       ELM(:,:,10,NEQ) = - (1.d0 - rGamm0zz) * fem_int( 2,rNuL)
       NLC(10,NEQ) = LQnz

       ELM(:,:,11,NEQ) =   (1.d0 - rGamm0zz) * fem_int(-2,rNuL,nC0z)
       NLC(11,NEQ) = 0
       NLCMAX(NEQ) = 11
    end if

  end subroutine LQnzCC

!***************************************************************
!
!   Ripple Trapped Beam Ion Density
!
!     X(:,LQr1) = n_r, the density of the ripple-trapped beam ions inside
!     the ripple wells, which occupy the fraction g = rip_rat of a flux
!     surface.  The inventory per unit volume is N_r = g n_r, and the
!     equation is the conservation law of N_r:
!
!       d(g n_r)/dt + d(g U n_r)/dV - d/dV(g K dn_r/dV)
!         = g S_perp + g nu_2 n_b - g nu_1 n_r - g nu_B n_r,
!
!     U = <|grad V|> U_brp (vertical grad-B drift, V space), K = <|grad V|^2> D_brp.
!     g stays inside every divergence: the flux vanishes where the wells end
!     (g = 0), and no 1/g appears.  The exchange terms use the element
!     matrices of LQb1, so n_b + g n_r is conserved exactly.  The advection is
!     in strong form, so the ions leave at the wall with the flux g U n_r
!     (the ripple loss); the diffusion has the natural condition.  Nodes with
!     no well on either side carry the identity dn_r/dt = 0.
!
!***************************************************************

  subroutine LQr1CC

    integer(4) :: NEQ = LQr1, nr
    real(8), dimension(0:NRMAX) :: gU

    if(RpplSW == 0.d0) then
       ELM(:,:,0,NEQ) = fem_int(1) * invDT
       NLC(0,NEQ) = LQr1
       NLCMAX(NEQ) = 0
       return
    end if

    gU(0) = 0.d0
    gU(1:NRMAX) = rip_rat(1:NRMAX) * suft(1:NRMAX) * Ubrp(1:NRMAX)

    ELM(:,:,0,NEQ) = fem_int(2,rip_rat) * invDT
    do nr = 0, NRMAX
       if(rip_rat(nr) /= 0.d0) cycle
       if(nr > 0) then
          if(rip_rat(nr-1) /= 0.d0) cycle
       end if
       if(nr < NRMAX) then
          if(rip_rat(nr+1) /= 0.d0) cycle
          ELM(nr+1,1,0,NEQ) = ELM(nr+1,1,0,NEQ) + hv(nr+1) * invDT
       else
          ELM(nr,4,0,NEQ) = ELM(nr,4,0,NEQ) + hv(nr) * invDT
       end if
    end do
    NLC(0,NEQ) = LQr1

    ! NBI perpendicular particle source falling into the ripple wells
    ! (M.H. Redi, et al., NF 35 (1995) 1191, p.1201)

    ELM(:,:,1,NEQ) =   fem_int(-2,SNBPDi,rip_rat)
    NLC(1,NEQ) = 0

    ! Trapping of and detrapping to the ripple-untrapped beam ions

    ELM(:,:,2,NEQ) =   fem_int(20,rNubrp2,rip_rat)
    NLC(2,NEQ) = LQb1

    ELM(:,:,3,NEQ) = - fem_int(20,rNubrp1,rip_rat)
    NLC(3,NEQ) = LQr1

    ! Relaxation to thermal ions

    ELM(:,:,4,NEQ) = - fem_int(20,rNuB,rip_rat)
    NLC(4,NEQ) = LQr1

    ! Ripple loss transport (convective)

    ELM(:,:,5,NEQ) = - fem_int(3,gU)
    NLC(5,NEQ) = LQr1

    ! Ripple loss transport (diffusive)

    ELM(:,:,6,NEQ) = - fem_int(37,sst,rip_rat,Dbrp)
    NLC(6,NEQ) = LQr1

    NLCMAX(NEQ) = 6

  end subroutine LQr1CC

!***************************************************************
!
!   Boundary condition
!
!      NLCR(NC,NQ,0) : inner boundary (magnetic axis)
!      NLCR(NC,NQ,1) : outer boundary (virtual wall)
!
!***************************************************************

  subroutine dump_axis_LQs8

!   Everything that enters the LQ*8 rows at the first few nodes, so that the
!   near-axis behaviour of BV_1s can be attributed.  Written once per time step
!   under IDIAG ones-digit >= 7; nothing here runs otherwise.
!
!   The axis row takes contributions from element 1 only.  With f = fipolsdt,
!   n = n_s and h = hv(1), and the stencils of fem_int ids 2, 5 and 22, it reads
!
!     -h/12 [ (3n0+n1) BV0 + (n0+n1) BV1 ]
!       - (rKilo/Z)(2f0+f1)/6 (p1-p0)
!       - (3f0n0+f0n1+f1n0+f1n1)/12 (Phi1-Phi0)  =  0
!
!   which expands to the continuum relation with an O(h) truncation, so the row
!   itself converges.  What is dumped here is each of those three pieces as the
!   code actually forms them, plus the coefficients they are built from.

    use tx_commons, only : IDIAG, NT, NTMAX, NRMAX, NEMAX, hv, vv, sdt, fipol, X, Var, &
         &                 LQm1, LQe1, LQi1, LQz1, LQe5, LQi5, LQz5, LQe8, LQi8, LQz8, &
         &                 achg, rKilo

    integer(4), parameter :: NRD = 3
    integer(4) :: nr, i, NEQ, LQs5, LQs8
    integer(4), save :: iua, nt_done = -1
    logical, save :: openeda = .false.
    real(8) :: rowmass, rowpres, rowphi, rowsum

    if(mod(IDIAG,10) < 7) return
    if(NT == nt_done) return
    nt_done = NT

    if(.not. openeda) then
       openeda = .true.
       open(newunit=iua, file='tx_axis.dat', status='replace', action='write')
       write(iua,'(A)') '# LQ*8 near-axis breakdown.  NRMAX, NT, species, NR, then:'
       write(iua,'(A)') '#   vv  hv  sdt  fipol  fipolsdt  n_s  p_s  Phi'
       write(iua,'(A)') '#   BLC(NR,1) ALC(NR,1) CLC(NR,1)   [mass, on LQ*8]'
       write(iua,'(A)') '#   BLC(NR,2) ALC(NR,2) CLC(NR,2)   [pressure, on LQ*5]'
       write(iua,'(A)') '#   BLC(NR,3) ALC(NR,3) CLC(NR,3)   [potential, on LQm1]'
       write(iua,'(A)') '#   mass.BV  pres.p  phi.Phi  sum   BV_1s'
    end if

    do i = 1, 3
       select case(i)
       case(1) ; NEQ = LQe8 ; LQs5 = LQe5
       case(2) ; NEQ = LQi8 ; LQs5 = LQi5
       case(3) ; NEQ = LQz8 ; LQs5 = LQz5
       end select
       do nr = 0, NRD
          ! the row as assembled, applied to the current X
          rowmass = BLC(nr,1,NEQ)*X(nr,NEQ)
          rowpres = BLC(nr,2,NEQ)*X(nr,LQs5)
          rowphi  = BLC(nr,3,NEQ)*X(nr,LQm1)
          if(nr < NRMAX) then
             rowmass = rowmass + ALC(nr,1,NEQ)*X(nr+1,NEQ)
             rowpres = rowpres + ALC(nr,2,NEQ)*X(nr+1,LQs5)
             rowphi  = rowphi  + ALC(nr,3,NEQ)*X(nr+1,LQm1)
          end if
          if(nr > 0) then
             rowmass = rowmass + CLC(nr,1,NEQ)*X(nr-1,NEQ)
             rowpres = rowpres + CLC(nr,2,NEQ)*X(nr-1,LQs5)
             rowphi  = rowphi  + CLC(nr,3,NEQ)*X(nr-1,LQm1)
          end if
          rowsum = rowmass + rowpres + rowphi
          write(iua,'(4I6,22ES16.8)') NRMAX, NT, i, nr, &
               & vv(nr), hv(max(nr,1)), sdt(nr), fipol(nr), fipol(nr)/sdt(nr), &
               & Var(nr,i)%n, X(nr,LQs5), X(nr,LQm1), &
               & BLC(nr,1,NEQ), ALC(nr,1,NEQ), CLC(nr,1,NEQ), &
               & BLC(nr,2,NEQ), ALC(nr,2,NEQ), CLC(nr,2,NEQ), &
               & BLC(nr,3,NEQ), ALC(nr,3,NEQ), CLC(nr,3,NEQ), &
               & rowmass, rowpres, rowphi, rowsum, X(nr,NEQ)
       end do
    end do
    flush(iua)

  end subroutine dump_axis_LQs8

  subroutine dump_axis_row(NEQ,tagc)

!   Every term of one equation's row at the first few nodes, as assembled and
!   applied to the current X, so that a non-smooth forcing can be attributed.
!   Written under IDIAG ones-digit >= 7 alongside dump_axis_LQs8.
!
!   For LQz1 this matters because that row has no second-order term at all --
!   the impurity continuity is first-order advection and reaction -- so n_z is
!   not smoothed by its own equation and inherits the regularity of whatever
!   forces it, term by term.

    use tx_commons, only : IDIAG, NT, NRMAX, X, NQMAX

    integer(4), intent(in) :: NEQ
    character(len=*), intent(in) :: tagc
    integer(4), parameter :: NRD = 5
    integer(4) :: nr, nc, nv
    integer(4), save :: iur, nt_done = -1
    logical, save :: openedr = .false.
    real(8) :: contrib(0:NCM), tot

    if(mod(IDIAG,10) < 7) return

    if(.not. openedr) then
       openedr = .true.
       open(newunit=iur, file='tx_axis_row.dat', status='replace', action='write')
       write(iur,'(A)') '# per-term contribution to one row.  eq NRMAX NT NR nterm  then' // &
            & ' (NLC, contribution) per term, then the total'
    end if

    do nr = 0, min(NRD,NRMAX)
       do nc = 0, NLCMAX(NEQ)
          nv = NLC(nc,NEQ)
          if(nv == 0) then
             contrib(nc) = PLC(nr,nc,NEQ)
          else
             contrib(nc) = BLC(nr,nc,NEQ) * X(nr,nv)
             if(nr < NRMAX) contrib(nc) = contrib(nc) + ALC(nr,nc,NEQ) * X(nr+1,nv)
             if(nr > 0)     contrib(nc) = contrib(nc) + CLC(nr,nc,NEQ) * X(nr-1,nv)
          end if
       end do
       tot = sum(contrib(0:NLCMAX(NEQ)))
       write(iur,'(A8,4I6,50(I4,ES15.7))') tagc, NRMAX, NT, nr, NLCMAX(NEQ), &
            & (NLC(nc,NEQ), contrib(nc), nc = 0, NLCMAX(NEQ))
       write(iur,'(A8,4I6,A,ES15.7)') tagc, NRMAX, NT, nr, -1, '   total', tot
    end do

  end subroutine dump_axis_row

  subroutine dump_axis_coefs

!   The coefficient arrays that force LQz1 and LQi1, at the first few nodes.

    use tx_commons, only : IDIAG, NT, NRMAX, vv, UgV, SiVsefA, SiVa6A, SiVizA, &
         &                 rNuL, rNuOL, SiLC, Var, PN01V, PN02V, PN0zV

    integer(4), parameter :: NRD = 5
    integer(4) :: nr
    integer(4), save :: iuc2
    logical, save :: openedc2 = .false.

    if(mod(IDIAG,10) < 7) return
    if(.not. openedc2) then
       openedc2 = .true.
       open(newunit=iuc2, file='tx_axis_coef.dat', status='replace', action='write')
       write(iuc2,'(A)') '# NRMAX NT NR  vv  UgV  SiVsefA  SiVa6A  SiVizA  rNuL  rNuOL' // &
            & '  SiLC  Te  ne  nz  n0z'
    end if
    do nr = 0, min(NRD,NRMAX)
       write(iuc2,'(3I6,12ES15.7)') NRMAX, NT, nr, vv(nr), UgV(nr), SiVsefA(nr), &
            & SiVa6A(nr), SiVizA(nr), rNuL(nr), rNuOL(nr), SiLC(nr), &
            & Var(nr,1)%T, Var(nr,1)%n, Var(nr,3)%n, PN0zV(nr)
    end do

  end subroutine dump_axis_coefs

  subroutine BOUNDARY(NR,LQ,ID,VAL)

    integer(4), intent(in) :: NR, LQ, ID
    real(8), intent(in), optional :: VAL
    integer(4), parameter :: NQMAXL = 50
    integer(4) :: NQ, NC, I, NB
    integer(4), save :: IMAX(1:NQM)
    type list
       integer(4) :: IDXNC
       integer(4) :: IDXNQ
    end type list
    type(list), save :: IDX(1:50,1:NQM,0:1) ! IDX(num. of B.C., equation, axis or bndry)

    ! NB = 0 at axis, = 1 at the boundary
    NB = NR / NRMAX

    ! === Dirichelt boundary condition ===
    if(ID == 0) then
       ! Initialize ALC, BLC and CLC at NR
       ALC(NR,:,LQ) = 0.d0
       BLC(NR,:,LQ) = 0.d0
       CLC(NR,:,LQ) = 0.d0
       PLC(NR,:,LQ) = 0.d0

       ! The boundary condition matrices IDX and IMAX are created if ICALA=0
       !   to store the information of the Dirichlet boundary conditions.
       ! ** IDX stores both the equation and term numbers, NQ and NC, where the variable LQ appears.
       ! ** IMAX stores the maximum number of IDX for each LQ.
       if(ICALA == 0) then
          I = 0
          do NQ = 1, NQMAX
             do NC = 1, NLCMAX(NQ)
                if(NLCR(NC,NQ,NB) == LQ) THEN
                   I = I + 1
                   IDX(I,LQ,NB)%IDXNC = NC
                   IDX(I,LQ,NB)%IDXNQ = NQ
                end if
             end do
          end do
          IMAX(LQ) = I
       else
       ! Read indices from the boundary condition matrices created at initial (ICALA=1)
          if(present(VAL)) then
             if(NR == 0) then ! axis
                do I = 1, IMAX(LQ)
                   NC = IDX(I,LQ,NB)%IDXNC
                   NQ = IDX(I,LQ,NB)%IDXNQ
                   PLC(NR  ,NC,NQ) = PLC(NR,NC,NQ)   + BLC(NR,NC,NQ)   * VAL
                   PLC(NR+1,NC,NQ) = PLC(NR+1,NC,NQ) + CLC(NR+1,NC,NQ) * VAL
                end do
             else ! wall
                do I = 1, IMAX(LQ)
                   NC = IDX(I,LQ,NB)%IDXNC
                   NQ = IDX(I,LQ,NB)%IDXNQ
                   PLC(NR  ,NC,NQ) = PLC(NR,NC,NQ)   + BLC(NR,NC,NQ)   * VAL
                   PLC(NR-1,NC,NQ) = PLC(NR-1,NC,NQ) + ALC(NR-1,NC,NQ) * VAL
                end do
             end if
          end if

          if(NR == 0) then ! axis
             do I = 1, IMAX(LQ)
                NC = IDX(I,LQ,NB)%IDXNC
                NQ = IDX(I,LQ,NB)%IDXNQ
                BLC(NR  ,NC,NQ) = 0.d0
                CLC(NR+1,NC,NQ) = 0.d0
             end do
          else ! wall
             do I = 1, IMAX(LQ)
                NC = IDX(I,LQ,NB)%IDXNC
                NQ = IDX(I,LQ,NB)%IDXNQ
                BLC(NR  ,NC,NQ) = 0.d0
                ALC(NR-1,NC,NQ) = 0.d0
             end do
          end if
       end if

       ! Diagonal term on handled variable
       BLC (NR,1,LQ) = 1.d0
       NLCR(1,LQ,NB) = LQ

       if(present(VAL)) PLC(NR,1,LQ) = PLC(NR,1,LQ) - VAL
    else
    ! === Neumann boundary condition ===
       NLCMAX(LQ) = NLCMAX(LQ) + 1
       PLC(NR,NLCMAX(LQ),LQ) = VAL * DTf(LQ)
       NLCR(NLCMAX(LQ),LQ,NB) = 0
    end if

  end subroutine BOUNDARY

!***************************************************************
!
!   Zero density gradient at the wall in constrained-DOF form (MDLWBC = 1)
!
!   For a linear system A x = b and the constraint x_N = x_(N-1), write
!   x = T y with T copying y_(N-1) into x_N; the constrained Galerkin system
!   is  T^T A T y = T^T b.  T^T adds row N to row N-1, A T adds column N to
!   column N-1.  In P1 finite elements this is the basis reduction
!   phi_(N-1) -> phi_(N-1) + phi_N, i.e. functions constant on the last
!   element, which is the P1 form of dn/drho = 0.
!
!   The system is kept square: row N-1 <- row N-1 + row N, and row N is
!   replaced by the constraint x_N - x_(N-1) = 0 (zero right-hand side).
!   x_N stays an unknown and the other equations keep their coupling to it,
!   so the column sum is carried by the constraint.  The time term and its
!   right-hand side are built from the same ALC/BLC/CLC arrays, and PLC is
!   folded with them, so the wall half-cell balance is not dropped: the sum
!   of the rows of the continuity equation, i.e. the global particle
!   balance, is unchanged.
!
!   Band storage: in row NR, ALC, BLC and CLC multiply X(NR+1,v), X(NR,v)
!   and X(NR-1,v), v = NLC(NC,LQ) (NLCR(NC,LQ,1) in row NRMAX).
!
!***************************************************************

  subroutine BOUNDARY_ZEROGRAD(LQ)

    integer(4), intent(in) :: LQ
    integer(4) :: NC

    do NC = 0, NLCMAX(LQ)
       if( NLCR(NC,LQ,1) /= NLC(NC,LQ) ) stop 'BOUNDARY_ZEROGRAD: NLCR /= NLC'
       BLC(NRMAX-1,NC,LQ) = BLC(NRMAX-1,NC,LQ) + CLC(NRMAX,NC,LQ)
       ALC(NRMAX-1,NC,LQ) = ALC(NRMAX-1,NC,LQ) + BLC(NRMAX,NC,LQ)
       if( NC >= 1 ) PLC(NRMAX-1,NC,LQ) = PLC(NRMAX-1,NC,LQ) + PLC(NRMAX,NC,LQ)
    end do

    ALC(NRMAX,:,LQ) = 0.d0
    BLC(NRMAX,:,LQ) = 0.d0
    CLC(NRMAX,:,LQ) = 0.d0
    PLC(NRMAX,:,LQ) = 0.d0

    NLCR(1,LQ,1) = LQ
    BLC(NRMAX,1,LQ) =   1.d0
    CLC(NRMAX,1,LQ) = - 1.d0

  end subroutine BOUNDARY_ZEROGRAD

end module tx_coefficients
