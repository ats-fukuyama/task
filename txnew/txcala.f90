module tx_coefficients
  use tx_commons
  use tx_core_module
  implicit none
  private
  real(8), parameter :: fourPisq = 4.d0 * Pi * Pi
  real(8), dimension(:,:,:,:), allocatable :: ELM
  real(8), dimension(:), allocatable ::  &
       & DTf, Dbrpft
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

    allocate(ELM(1:NEMAX,1:4,0:NCM,1:NQMAX), source=0.d0)
    allocate(DTf(1:NQMAX))

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

    !     Coefficients

    NSMP = NSM + 1

    allocate(Dbrpft, source=array_init_NR)
    allocate(rrtinv, aatinv, bthcoinv, ParatoZeta, ZetatoPara, &
       &     aatq, PsitdotVBth, fipolinv, fipolsdt, dZetatoPara, ribi, ribsdt, &
       &     RatCXSNBe, RatCXSNBi, Ratdenhyd, source=array_init_NR)
    allocate(UsrgV, fipolsdtPNsV, BUsparVbbt, dlnPNsV, BpBVsdiag, tormflux, dtormflux, &
       &     ribsdtNs, source=array_init_NRNS)
    !allocate(RatPNbinv, source=array_init_NRNS)
    ! not working with ifort15, but working with ifort18
    !    allocate(Ratden(0:NRMAX,NSMP,NSMP), zChi(0:NRMAX,NSM,NSM),source=0.d0)
    allocate(Ratden(0:NRMAX,NSMP,NSMP), zChi(0:NRMAX,NSM,NSM))
    Ratden = 0.d0 ; zChi = 0.d0

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
    ! When ripple effect is on (FSRP /= 0), ripple diffusion term will be activated.
    ! Then we must impose two boundary conditions at each equation.
!    call BOUNDARY(0    ,LQb3,0)
!    call BOUNDARY(NRMAX,LQb3,0)
!    IF(FSRP /= 0.d0) call BOUNDARY(0,LQr1,0)

    ! Neumann condition of the continuity equation at the boundary is naturally
    ! imposed through the diffusion term in the pressure equation. However, this
    ! condition is very weak because it comes from nondiagonal term. Then by setting
    ! the density constant in the last element, we can force the density gradient
    ! to be nought at the boundary explicitly.

!    call BOUNDARY(NRMAX,LQe1,0,Var(NRMAX-1,1)%n)
!    call BOUNDARY(NRMAX,LQi1,0,Var(NRMAX-1,2)%n)
!    call BOUNDARY(NRMAX,LQz1,0,Var(NRMAX-1,3)%n)
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

    deallocate(ELM, DTf)
    deallocate(Dbrpft)
    deallocate(rrtinv, aatinv, bthcoinv, ParatoZeta, ZetatoPara, &
         &     aatq, PsitdotVBth, fipolinv, fipolsdt, fipolsdtPNsV, dZetatoPara, ribi, ribsdt, &
         &     RatCXSNBe, RatCXSNBi, Ratdenhyd)
    deallocate(UsrgV, BUsparVbbt, dlnPNsV, BpBVsdiag, tormflux, dtormflux, ribsdtNs)
    deallocate(Ratden, zChi)
    !deallocate(RatPNbinv)

    if(ICALA == 0) ICALA = 1

  end subroutine TXCALA

!***************************************************************
!
!   Coefficients for Equations
!
!**************************************************************

  subroutine LQcoef

    use tx_interface, only : dfdx
    integer(4) :: i

    Dbrpft(:) = Dbrp(:) * ft(:)

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

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQe1

    ! Divergence

    ELM(:,:,1,NEQ) = - fem_int(4)
    NLC(1,NEQ) = LQe2

    ! Divergence + Grid velocity

    ELM(:,:,2,NEQ) =   fem_int(5,UgV)
    NLC(2,NEQ) = LQe1

    ! Ionization of n01, n02 and n03 by electron impact

    ELM(:,:,3,NEQ) =   FSION * 1.D20 * fem_int(20,SiVizA,Var(:,i)%n)
    NLC(3,NEQ) = LQn1

    ELM(:,:,4,NEQ) =   FSION * 1.D20 * fem_int(20,SiVizA,Var(:,i)%n) * ThntSW
    NLC(4,NEQ) = LQn2

    ELM(:,:,5,NEQ) =   FSION * 1.D20 * fem_int(20,SiVizA,Var(:,i)%n) * BeamSW
    NLC(5,NEQ) = LQn3

    ! Increment of electrons due to effective ionization of impurities

    ELM(:,:,6,NEQ) =   FSION * 1.D20 * fem_int(20,SiVsefA,Var(:,1)%n) * achg(3) 
    NLC(6,NEQ) = LQnz

    ! Loss of electrons due to effective recombination of impurities

    ELM(:,:,7,NEQ) = - FSCX  * 1.D20 * fem_int(20,SiVa6A, Var(:,1)%n) * achg(3)
    NLC(7,NEQ) = LQz1

    ! Loss to divertor

    ELM(:,:,8,NEQ) = -             fem_int( 2,rNuL)
    NLC(8,NEQ) = LQe1

    ELM(:,:,9,NEQ) =   PNsDIV(1) * fem_int(-1,rNuL)
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

    ! Loss to divertor

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

    ! Loss to divertor

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

       ELM(:,:,18,NEQ) =   achg(i) * 1.d-3 * fem_int(  2,PsidotV)
       NLC(18,NEQ) = LQe7

       ELM(:,:,19,NEQ) =   achg(2) * 1.d-3 * fem_int(  2,PsidotV)
       NLC(19,NEQ) = LQi7

       ELM(:,:,20,NEQ) =   achg(3) * 1.d-3 * fem_int(  2,PsidotV)
       NLC(20,NEQ) = LQz7
!       ELM(:,:,18,NEQ) =   achg(i) * 1.d-3 * fem_int(-20,PsidotV,Var(:,i)%n,Var(:,i)%UphR)
!       NLC(18,NEQ) = 0
!
!       ELM(:,:,19,NEQ) =   achg(2) * 1.d-3 * fem_int(-20,PsidotV,Var(:,2)%n,Var(:,2)%UphR)
!       NLC(19,NEQ) = 0
!
!       ELM(:,:,20,NEQ) =   achg(3) * 1.d-3 * fem_int(-20,PsidotV,Var(:,3)%n,Var(:,3)%UphR)
!       NLC(20,NEQ) = 0

       ! Loss to divertor

       ELM(:,:,21,NEQ) = -                     fem_int( 2,rNuL)
       NLC(21,NEQ) = LQe5

       ELM(:,:,22,NEQ) =           PNsDIV(i) * fem_int(-2,rNuL,Var(:,i)%T)
       NLC(22,NEQ) = 0

       ELM(:,:,23,NEQ) = - 1.5d0             * fem_int( 2,rNuLTs(:,i))
       NLC(23,NEQ) = LQe5

       ELM(:,:,24,NEQ) =   1.5d0 * PTsDIV(i) * fem_int(-2,rNuLTs(:,i),Var(:,i)%n)
       NLC(24,NEQ) = 0
!!$       ELM(:,:,24,NEQ) =   1.5d0 * PTsDIV(i) * fem_int( 2,rNuLTs(:,i))
!!$       NLC(24,NEQ) = LQe1

       ! Ionization loss of n01, n02 and n03 by electron impact

       ELM(:,:,25,NEQ) = - (EION * 1.D-3) * FSION * 1.D20 * fem_int(20,SiVizA,Var(:,i)%n)
       NLC(25,NEQ) = LQn1

       ELM(:,:,26,NEQ) = - (EION * 1.D-3) * FSION * 1.D20 * fem_int(20,SiVizA,Var(:,i)%n)
       NLC(26,NEQ) = LQn2

       ELM(:,:,27,NEQ) = - (EION * 1.D-3) * FSION * 1.D20 * fem_int(20,SiVizA,Var(:,i)%n) * BeamSW
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

    ! Loss to divertor

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

    ELM(:,:, 0,NEQ) = fem_int(1) * invDT
    NLC( 0,NEQ) = LQi1

    ! Divergence

    ELM(:,:, 1,NEQ) = - fem_int(4)
    NLC( 1,NEQ) = LQi2

    ! Divergence + Grid velocity

    ELM(:,:, 2,NEQ) =   fem_int(5,UgV)
    NLC( 2,NEQ) = LQi1

    ! Ionization of n01, n02 and n03 by electron impact

    ELM(:,:, 3,NEQ) =     FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n)
    NLC( 3,NEQ) = LQn1

    ELM(:,:, 4,NEQ) =     FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n) * ThntSW
    NLC( 4,NEQ) = LQn2

    ELM(:,:, 5,NEQ) =     FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n) * BeamSW
    NLC( 5,NEQ) = LQn3

    ! Loss to divertor

    ELM(:,:, 6,NEQ) = -             fem_int( 2,rNuL)
    NLC( 6,NEQ) = LQi1

    ELM(:,:, 7,NEQ) =   PNsDIV(i) * fem_int(-1,rNuL)
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

    ! Loss to divertor

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

    integer(4) :: NEQ = LQi4, i = 2, j
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

    ! Loss to divertor

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

       ! Loss to divertor

       ELM(:,:,12,NEQ) = -             fem_int( 2,rNuL)
       NLC(12,NEQ) = LQi5

       ELM(:,:,13,NEQ) =   PNsDIV(i) * fem_int(-2,rNuL,Var(:,i)%T)
       NLC(13,NEQ) = 0

       ELM(:,:,14,NEQ) = - 1.5d0            * fem_int( 2,rNuLTs(:,i))
       NLC(14,NEQ) = LQi5

       ELM(:,:,15,NEQ) =   1.5d0 * PTsDIV(i) * fem_int(-2,rNuLTs(:,i),Var(:,i)%n)
       NLC(15,NEQ) = 0
!       ELM(:,:,15,NEQ) =   1.5d0 * PTsDIV(i)    * fem_int( 2,rNuLTs(:,i))
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

    ! Loss to divertor

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

    ! Loss to divertor

    ELM(:,:,5,NEQ) = -             fem_int( 2,rNuL)
    NLC(5,NEQ) = LQz1

    ELM(:,:,6,NEQ) =   PNsDIV(i) * fem_int(-1,rNuL)
    NLC(6,NEQ) = 0

    ! Loss cone loss

    ELM(:,:,7,NEQ) =   fem_int(-1,SiLC)
    NLC(7,NEQ) = 0
 
    ! Ion orbit loss

    ELM(:,:,8,NEQ) = - fem_int(2,rNuOL)
    NLC(8,NEQ) = LQz1

    NLCMAX(NEQ) = 8

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

    ! Loss to divertor

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

    integer(4) :: NEQ = LQz4, i = 3, j
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

    ! Loss to divertor

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

       ! Loss to divertor

       ELM(:,:,12,NEQ) = -             fem_int( 2,rNuL)
       NLC(12,NEQ) = LQz5

       ELM(:,:,13,NEQ) =   PNsDIV(i) * fem_int(-2,rNuL,Var(:,i)%T)
       NLC(13,NEQ) = 0

       ELM(:,:,14,NEQ) = - 1.5d0            * fem_int( 2,rNuLTs(:,i))
       NLC(14,NEQ) = LQz5

       ELM(:,:,15,NEQ) =   1.5d0 * PTsDIV(i) * fem_int(-2,rNuLTs(:,i),Var(:,i)%n)
       NLC(15,NEQ) = 0
!       ELM(:,:,15,NEQ) =   1.5d0 * PTsDIV(i)    * fem_int( 2,rNuLTs(:,i))
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

    ! Loss to divertor

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

    ELM(:,:,4,NEQ) = - fem_int(-2,SNBPDi,rip_rat)
    NLC(4,NEQ) = 0

    ! Relaxation to thermal ions

    ELM(:,:,5,NEQ) = - fem_int(2,rNuB)
    NLC(5,NEQ) = LQb1

    ! Loss to divertor

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

       ! Loss to divertor

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

       ! Momentum loss due to collisional ripple trapping

       ELM(:,:,14,NEQ) = - fem_int(20+L3,rNubrp2,rip_rat) * RpplSW
       NLC(14,NEQ) = LQb3

       ! Momentum diffusion arising from beam ion convective flux due to ripple

       ELM(:,:,15,NEQ) = - fem_int(37,sst,BUbparV,Dbrpft)
       NLC(15,NEQ) = LQb1

    end if

    NLCMAX(NEQ) = 15

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

       ! Loss to divertor

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

    ELM(:,:,2,NEQ) = - fem_int(27,sst,D01)
    NLC(2,NEQ) = LQn1

    ! Ionization by electron impact
    !   cf. [C.E. Singer et al., Comput Phys. Commun. 49 (1988) 275, (2.9.1a)]

    ELM(:,:,3,NEQ) = - FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n)
    NLC(3,NEQ) = LQn1

    ! Generation of thermal neutrals by charge exchange

    ELM(:,:,4,NEQ) = - FSCX * 1.D20 * fem_int(20,SiVcxA,Var(:,2)%n)
    NLC(4,NEQ) = LQn1

    ! Recycling from divertor

    ELM(:,:,5,NEQ) =   rGamm0             * fem_int( 2,rNuL)
    NLC(5,NEQ) = LQi1

    ELM(:,:,6,NEQ) = - rGamm0 * PNsDIV(2) * fem_int(-1,rNuL)
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

    integer(4) :: NEQ = LQnz, i = 3

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQnz

    ! Grid velocity

    ELM(:,:,1,NEQ) = - fem_int(6,UgV)
    NLC(1,NEQ) = LQnz

    ! Diffusion of neutrals

    ELM(:,:,2,NEQ) = - fem_int(27,sst,D01)
    NLC(2,NEQ) = LQnz

    ! Effective ionization generating C^{6+}

    ELM(:,:,3,NEQ) = - FSION * 1.D20 * fem_int(20,SiVsefA,Var(:,1)%n)
    NLC(3,NEQ) = LQnz

    ! Effective recombination generating C^{0-5+}

    ELM(:,:,4,NEQ) =   FSCX  * 1.D20 * fem_int(20,SiVa6A, Var(:,1)%n)
    NLC(4,NEQ) = LQz1

    ! Recycling from divertor

    ELM(:,:,5,NEQ) =   rGamm0z             * fem_int( 2,rNuL)
    NLC(5,NEQ) = LQz1

    ELM(:,:,6,NEQ) = - rGamm0z * PNsDIV(i) * fem_int(-1,rNuL)
    NLC(6,NEQ) = 0

    NLCMAX(NEQ) = 6

  end subroutine LQnzCC

!***************************************************************
!
!   Ripple Trapped Beam Ion Density (SUPG)
!
!***************************************************************

  subroutine LQr1CC

    integer(4) :: ne, NEQ = LQr1
    real(8) :: RUbrpl, SDbrpl, peclet, coef

    do ne = 1, nemax
       RUbrpl = RUbrp(ne-1) + RUbrp(ne)
       SDbrpl = Dbrp(ne-1)*vv(ne-1) + Dbrp(ne)*vv(ne)
       if (RUbrpl == 0.d0) then
          coef = 0.d0
       elseif (SDbrpl == 0.d0) then
          coef = 1.d0 / sqrt(15.d0) * hv(ne)
       else
          peclet = 0.5d0 * RUbrpl * hv(ne) / SDbrpl
          coef = 0.5d0 * langevin(peclet) * hv(ne)
       end if
       coef=0.d0 !!! no SUPG

       ELM(NE,1:4,0,NEQ) =   fem_int_point(1,NE) * invDT &
            &              + fem_int_point(8,NE) * coef * invDT
       NLC(0,NEQ) = LQr1

       ! NBI perpendicular particle source (Both charge exchange and ionization)
       ! (Beam ions by perpendicular NBI have few parallel velocity, hence
       !  they can be easily trapped by ripple wells if they go into the ripple well region.)
       !  (M.H. Redi, et al., NF 35 (1995) 1191, p.1201 sixth line from the bottom
       !   at right-hand-side column)
       
       ELM(NE,1:4,1,NEQ) =(  fem_int_point(-1,NE,SNBPDi) &
            &              + fem_int_point(-8,NE,SNBPDi) * coef) * RpplSW
       NLC(1,NEQ) = 0

       ! Ripple trapped beam ions collision with otherwise beam ions

       ELM(NE,1:4,2,NEQ) =(  fem_int_point(2,NE,rNubrp2) &
            &              + fem_int_point(9,NE,rNubrp2) * coef) * RpplSW
       NLC(2,NEQ) = LQb1

       ELM(NE,1:4,3,NEQ) =(- fem_int_point(2,NE,rNubrp1) &
            &              - fem_int_point(9,NE,rNubrp1) * coef) * RpplSW
       NLC(3,NEQ) = LQr1

       ! Relaxation to thermal ions

       ELM(NE,1:4,4,NEQ) =(- fem_int_point(2,NE,rNuB) &
            &              - fem_int_point(9,NE,rNuB) * coef) * RpplSW
       NLC(4,NEQ) = LQr1

       ! Ripple loss transport (convective)

       ELM(NE,1:4,5,NEQ) =(- fem_int_point( 3,NE,RUbrp) &
            &              - fem_int_point(10,NE,RUbrp) * coef) * RpplSW
       NLC(5,NEQ) = LQr1

!!rp_conv       ELM(NE,1:4,5,NEQ) = - fem_int_point(2,NE,rNubL)
!!rp_conv       NLC(5,NEQ) = LQr1
!!rp_conv       ELM(NE,1:4,5,NEQ) = - fem_int_point(-2,NE,rNubL,PNbrpV)
!!rp_conv       NLC(5,NEQ) = 0

       ! Ripple loss transport (diffusive)

       ELM(NE,1:4,6,NEQ) = - fem_int_point(27,NE,sst,Dbrp) * RpplSW
       NLC(6,NEQ) = LQr1
    end do

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
!   Approximate Langevin function
!
!***************************************************************

  real(8) function langevin(x) result(y)

    real(8), intent(in) :: x

    if (x < -3.d0) then
       y = - 1.d0 - 1.d0 / x
    elseif (x > 3.d0) then
       y =   1.d0 - 1.d0 / x
    else
       y = x / 3.d0 * (1.d0 - abs(x) / 9.d0)
    end if

  end function langevin

end module tx_coefficients
