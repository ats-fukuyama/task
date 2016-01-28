!     $ID: txcala.f90,v 1.88 2011/05/02 06:44:16 fukuyama Exp $
module tx_coefficients
  use tx_commons
  use tx_core_module
  implicit none
  private
  real(8), dimension(:,:,:,:), allocatable :: ELM
  real(8), dimension(:), allocatable ::  &
       & DTf, Dbrpft, &
       & Chie1, Chie2, Chii1, Chii2
!       & rGASPFA, rNueHLphthVR, rNuiHLphthVR   ! for UHphSwitch == 1 miki_m 2010/8/26
  real(8), dimension(:), allocatable :: &
       & rrtinv, aatinv, bthcoinv, FqhatsqI, ParatoZeta, ZetatoPara, &
       & aatq, Bpsqbbt, PsitdotVBth, &
       & fipolsdt, fipolsdtPNbV, dlogZetatoPara
!!rp_conv       &, rNubLL
  real(8), dimension(:,:), allocatable :: UsrgV, RatPN, RatPNb, &!RatPNbinv, &
       & fipolsdtPNsV, BUsparVbbt, dlnPNsV
  real(8) :: DTt, invDT, BeamSW, RpplSW, ThntSW, &
       &     fact = 1.d0 ! <= SOL loss accelerator
  integer(4), save :: ICALA = 0!, ICALA2 = 0
  integer(4) :: L2, L3, L6
  
  public :: TXCALA

contains

!***************************************************************
!
!   Calculate ALC, BLC, CLC, PLC
!
!***************************************************************

  SUBROUTINE TXCALA

    INTEGER(4) :: NE, NR, NC, NQ, N, iHvsLC

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

    allocate(ELM(1:NEMAX,1:4,0:NCM,1:NQMAX),DTf(1:NQMAX))

    !     Preconditioning

    invDT = 1.d0 / DT
!!$    IF(DT <= 2.D-5) THEN
!!$       DTt          = 1.5d1
!!$       DTf(1)       = 1.d0
!!$       DTf(2:NQMAX) = 1.5d1
!!$    ELSE
       DTt          = 1.d0
       DTf(1:NQMAX) = 1.d0
!!$    END IF

    ! In case of ohmic heating (i.e. no NBI), largeness of the term related to beam
    ! components comparable to that of other terms related to finite variables
    ! sometimes induces spurious values in beam components even if no NBI is activated.
    ! To avoid this, setting these terms, especially related to electromagnetic potentials,
    ! to zero is desirable during no NBI.

    IF(PNBH == 0.D0 .and. PNbV(0) == 0.D0) THEN
       BeamSW = 0.D0
    ELSE
       BeamSW = 1.D0
    END IF
    IF(FSRP == 1.D0) THEN
       RpplSW = BeamSW
    ELSE
       RpplSW = 0.D0
    END IF

    ! If ThntSW = 0, there is no particle source from PN02V.
    ThntSW = 1.D0

    ! SUPG
    L2 =     iSUPG2  * 100
    L3 = abs(iSUPG3) * 100
    L6 =     iSUPG6  * 100

    !     Coefficients

    N = NRMAX

    allocate(Dbrpft(0:N), Chie1(0:N), Chie2(0:N), Chii1(0:N), Chii2(0:N))
!       &     rGASPFA(0:N), rNueHLphthVR(0:N), rNuiHLphthVR(0:N))
    allocate(UsrgV(0:N,NSM), rrtinv(0:N), aatinv(0:N), bthcoinv(0:N), &
       &     FqhatsqI(0:N), ParatoZeta(0:N), ZetatoPara(0:N), &
       &     aatq(0:N), Bpsqbbt(0:N), PsitdotVBth(0:N), &
       &     fipolsdt(0:N), fipolsdtPNbV(0:N), dlogZetatoPara(0:N))
    allocate(RatPN(0:N,NSM), RatPNb(0:N,NSM), &!RatPNbinv(0:N,NSM), &
         &   fipolsdtPNsV(0:N,NSM), BUsparVbbt(0:N,NSM), dlnPNsV(0:N,NSM))

!$omp parallel
!$omp workshare
    ! Initialize arrays
    ELM = 0.d0
    ALC = 0.d0
    BLC = 0.d0
    CLC = 0.d0
    NLCR = 0
    NLC = 0
    PLC = 0.d0
    NLCMAX = 0
!$omp end workshare

    CALL LQCOEF

!$omp sections
!$omp section
    !     Maxwell

    CALL LQm1CC
    CALL LQm2CC
    CALL LQm3CC
    CALL LQm4CC
    CALL LQm5CC

    !     Electron

    CALL LQe1CC
    CALL LQe2CC
!$omp section
    CALL LQe3CC
    CALL LQe4CC
    CALL LQe5CC
    CALL LQe6CC ! Heat flux
    CALL LQe7CC ! Rotation frequency

    !     Ion

    CALL LQi1CC
    CALL LQi2CC
!$omp section
    CALL LQi3CC
    CALL LQi4CC
    CALL LQi5CC
    CALL LQi6CC ! Heat flux
    CALL LQi7CC ! Rotation frequency

    !     Beam

    CALL LQb1CC
    CALL LQb2CC
!$omp section
    CALL LQb3CC
    CALL LQb4CC
    CALL LQb7CC

    !     Neutral

    CALL LQn1CC
    CALL LQn2CC
    CALL LQn3CC

    !     Ripple trapped beam

    CALL LQr1CC
!$omp end sections

    !     Elemental equations -> Nodal equations

!$omp do private(NC,iHvsLC)
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
!$omp end do
!$omp end parallel

    NLCR(0:NCM,1:NQMAX,0) = NLC(0:NCM,1:NQMAX)
    NLCR(0:NCM,1:NQMAX,1) = NLC(0:NCM,1:NQMAX)

    !     Dirichlet condition

!    if(NTCUM >= 3) then
!       if(NTCUM == 3 .and. ICALA2 == 0) then
!          ICALA = 0 ; ICALA2 = 1
!       end if
    CALL BOUNDARY(0    ,LQe2,0)
    CALL BOUNDARY(0    ,LQi2,0)
    ! --- Excess B.C. for numerical stability ---
    !   NOTE: This condition originally should be satisfied by the radial force balance
    !         equation. However, due to numerical accuracy and complex nonlinearlity of
    !         the system, they are not always satisfied. Suppressing the flapping of the
    !         ion parallel flow (LQi3) at the magnetic axis, imposing the following condition
    !         is favorable. It shouldn't be set if the perpendicular viscosity is included
    !         in the parallel flow equations (LQe3, LQi3).
!!unstable    CALL BOUNDARY(0    ,LQe3,0,Var(0,1)%RUph*bbt(0)/fipol(0))
    if(FSPARV(2) == 0.d0) CALL BOUNDARY(0    ,LQi3,0,Var(0,2)%RUph*bbt(0)/fipol(0))
    ! -----------------------------------------------------------
!    end if
    CALL BOUNDARY(NRMAX,LQm1,0)
    CALL BOUNDARY(0    ,LQm2,0)
!!!    CALL BOUNDARY(NRMAX,LQe2,0) ! Violates quasi-neutrality and breaks down calculation
!    CALL BOUNDARY(NRMAX,LQe3,0)   ! OK
!!!    CALL BOUNDARY(NRMAX,LQe4,0) ! Not impose this B.C. which violates quasi-neutrality and 
                                   ! break down calculation, because this eq. yields particle flux.
!!!    CALL BOUNDARY(NRMAX,LQi2,0) ! Violates quasi-neutrality and breaks down calculation
!    CALL BOUNDARY(NRMAX,LQi3,0)   ! OK
!!!    CALL BOUNDARY(NRMAX,LQi4,0) ! Not impose this B.C. which violates quasi-neutrality and 
                                   ! break down calculation, because this eq. yields particle flux.

    CALL BOUNDARY(NRMAX,LQn2,0)
    CALL BOUNDARY(NRMAX,LQn3,0)

    CALL BOUNDARY(0    ,LQb2,0)
    ! When ripple effect is on (FSRP /= 0), ripple diffusion term will be activated.
    ! Then we must impose two boundary conditions at each equation.
!    CALL BOUNDARY(0    ,LQb3,0)
!    CALL BOUNDARY(NRMAX,LQb3,0)
!    IF(FSRP /= 0.D0) CALL BOUNDARY(0,LQr1,0)

    ! Neumann condition of the continuity equation at the boundary is naturally
    ! imposed through the diffusion term in the pressure equation. However, this
    ! condition is very weak because it comes from nondiagonal term. Then by setting
    ! the density constant in the last element, we can force the density gradient
    ! to be nought at the boundary explicitly.

!    CALL BOUNDARY(NRMAX,LQe1,0,Var(NRMAX-1,1)%n)
!    CALL BOUNDARY(NRMAX,LQi1,0,Var(NRMAX-1,2)%n)
!!    CALL BOUNDARY(0,LQe3,0,Var(1,1)%BUpar)
!!    CALL BOUNDARY(0,LQi3,0,Var(1,2)%BUpar)

    !     Integral term stemming from integration by parts in the diffusion term

    CALL BOUNDARY(NRMAX,LQm2,1, ckt(NRMAX)*fipol(NRMAX)/(4.d0*Pisq*rMU0))
    CALL BOUNDARY(NRMAX,LQm3,1, 2.D0*Pi*rMUb1*rIp*1.D6)
    CALL BOUNDARY(NRMAX,LQn1,1, suft(NRMAX)*rGASPF)

    deallocate(ELM, DTf)
    deallocate(Dbrpft, Chie1, Chie2, Chii1, Chii2)
!       &       rGASPFA, rNueHLphthVR, rNuiHLphthVR)
    deallocate(UsrgV, rrtinv, aatinv, bthcoinv, FqhatsqI, &
       &       ParatoZeta, ZetatoPara, RatPN, RatPNb, &!RatPNbinv, &
       &       aatq, Bpsqbbt, PsitdotVBth, &
       &       fipolsdt, fipolsdtPNsV, dlnPNsV, fipolsdtPNbV, &
       &       BUsparVbbt, dlogZetatoPara)

    IF(ICALA == 0) ICALA = 1

    RETURN
  END SUBROUTINE TXCALA

!***************************************************************
!
!   Coefficients for Equations
!
!**************************************************************

  SUBROUTINE LQCOEF

    use tx_interface, only : dfdx
    INTEGER(4) :: NR, i

    fipolsdt(:)     = fipol(:) / sdt(:)
!$omp workshare
    fipolsdtPNbV(:) = fipolsdt(:) * PNbVinv(:) * 1.d-20

    Dbrpft(:) = Dbrp(:) * ft(:)

    Chie1(:) = Chie(:) + ChiNCTe(:) + ChiNCpe(:)
    Chie2(:) = Chie(:) + ChiNCTe(:)
    Chii1(:) = Chii(:) + ChiNCTi(:) + ChiNCpi(:)
    Chii2(:) = Chii(:) + ChiNCTi(:)

    RatPN(:,1) = Var(:,1)%n / Var(:,2)%n ! Ne/Ni
    RatPN(:,2) = Var(:,2)%n / Var(:,1)%n ! Ni/Ne

    bthcoinv(:)    = 1.d0 / bthco(:)
    Bpsqbbt(:)     = Bpsq(:) / bbt(:)
    PsitdotVBth(:) = 4.d0 * Pisq * sdt(:) * PsitdotV(:)
    aatq(:)        = aat(:) / Q(:)

    UsrgV(:,1)    = FSADV * Var(:,1)%UrV - UgV(:)
    UsrgV(:,2)    = FSADV * Var(:,2)%UrV - UgV(:)
    rrtinv(:)     = 1.d0 / rrt(:) ! 1/<R^2>
    aatinv(:)     = 1.d0 / aat(:) ! 1/<R^-2>
    FqhatsqI(:)   = Fqhatsq(:) / fipol(:)
    ParatoZeta(:) = fipol(:) / bbt(:) ! I/<B^2> : parallel to toroidal
    ZetatoPara(:) = bbt(:) / fipol(:) ! <B^2>/I : toroidal to parallel
!$omp end workshare
    dlogZetatoPara(:) = dfdx(vv,log(ZetatoPara),NRMAX,0)   ! dln(<B^2>/I)/dV

    do i = 1, NSM
       fipolsdtPNsV(:,i) = fipolsdt(:) / (Var(:,i)%n * 1.d20)
       BUsparVbbt(:,i)   = Var(:,i)%BUpar / bbt(:)
       dlnPNsV(:,i)      = dfdx(vv,Var(:,i)%n,NRMAX,0) / Var(:,i)%n
       RatPNb(:,i)       = Var(:,i)%n * PNbVinv(:)
!       RatPNbinv(:,i)    = PNbV(:) / Var(:,i)%n
    end do

!    rGASPFA(0:NRMAX-1) = 0.D0
!    rGASPFA(NRMAX) = rGASPF

!!rp_conv    rNubLL(:) = rNubL(:) * rip_rat(:)

!    rNueHLphthVR(1:NRMAX)   = rNueHLphth(1:NRMAX) / R(1:NRMAX)
!    rNueHLphthVR(0)         = 0.D0 ! Any value is OK. (Never affect the result.)
!    rNuiHLphthVR(1:NRMAX)   = rNuiHLphth(1:NRMAX) / R(1:NRMAX)
!    rNuiHLphthVR(0)         = 0.D0 ! Any value is OK. (Never affect the result.)

  END SUBROUTINE LQCOEF

!***************************************************************
!
!   Poisson Equation: phi
!
!**************************************************************

  SUBROUTINE LQm1CC

    integer(4) :: NEQ = LQm1
    real(8) :: sqeps0inv = 1.d0 / sqeps0

    ! phi'(0) : 0

    ELM(:,:,1,NEQ) =   sqeps0 / AEE * 1.D-20 * fem_int(12,sst)
    NLC(1,NEQ) = LQm1

    ELM(:,:,2,NEQ) = - achg(1) * fem_int(1) * sqeps0inv
    NLC(2,NEQ) = LQe1

    ELM(:,:,3,NEQ) = - achg(2) * fem_int(1) * sqeps0inv
    NLC(3,NEQ) = LQi1

    ELM(:,:,4,NEQ) = - achgb   * fem_int(1) * sqeps0inv * BeamSW
    NLC(4,NEQ) = LQb1

    ELM(:,:,5,NEQ) = - achgb   * fem_int(2,rip_rat) * sqeps0inv * RpplSW
    NLC(5,NEQ) = LQr1

    ! phi(b) : 0

    NLCMAX(NEQ) = 5
    RETURN
  END SUBROUTINE LQm1CC

!***************************************************************
!
!   Ampere's Law: psit'
!
!***************************************************************

  SUBROUTINE LQm2CC

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

   ! Beam ion current

    ELM(:,:,6,NEQ) =   achgb   * AEE * 1.D20 * fem_int(2,bthcoinv) * BeamSW
    NLC(6,NEQ) = LQb3

    ELM(:,:,7,NEQ) = - achgb   * AEE * 1.D20 * fem_int(20,bthcoinv,fipol) * BeamSW
    NLC(7,NEQ) = LQb7

    ! (NRMAX) : ckt(NRMAX)*Ivac/(4 pi^2)

    NLCMAX(NEQ) = 7
    RETURN
  END SUBROUTINE LQm2CC

!***************************************************************
!
!   Ampere's Law: psi'
!
!***************************************************************

  SUBROUTINE LQm3CC

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

    ! Beam ion current

    ELM(:,:,4,NEQ) = - rMUb1 * achgb   * AEE * 1.D20 * fem_int(1) * BeamSW
    NLC(4,NEQ) = LQb7

!    ! Virtual current for helical system
!
!    ELM(:,:,5,NEQ) = - rMUb1 * fem_int(-1,AJV) / rr
!    NLC(5,NEQ) = 0

    ! psi'(NRMAX) : Ip

    NLCMAX(NEQ) = 4
    RETURN
  END SUBROUTINE LQm3CC

!**************************************************************
!
!   Faraday's Law : psi
!
!***************************************************************

  SUBROUTINE LQm4CC

    integer(4) :: NEQ = LQm4

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQm4

    ! psi'

    ELM(:,:,1,NEQ) = fem_int(1) / rMUb2
    NLC(1,NEQ) = LQm3

    NLCMAX(NEQ) = 1
    RETURN
  END SUBROUTINE LQm4CC

!***************************************************************
!
!   Faraday's Law : psit
!
!***************************************************************

  SUBROUTINE LQm5CC

    integer(4) :: NEQ = LQm5

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQm5

    ! psit'

    ELM(:,:,1,NEQ) = fem_int(1) / rMU0
    NLC(1,NEQ) = LQm2

    NLCMAX(NEQ) = 1
    RETURN
  END SUBROUTINE LQm5CC

!***************************************************************
!
!   Electron Density Equation
!
!***************************************************************

  SUBROUTINE LQe1CC

    integer(4) :: NEQ = LQe1

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQe1

    ! Divergence

    ELM(:,:,1,NEQ) = - fem_int(4)
    NLC(1,NEQ) = LQe2

    ! Divergence + Grid velocity

    ELM(:,:,2,NEQ) =   fem_int(5,UgV)
    NLC(2,NEQ) = LQe1

    ! Ionization of n01, n02 and n03

    ELM(:,:,3,NEQ) =   FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n)
    NLC(3,NEQ) = LQn1

    ELM(:,:,4,NEQ) =   FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n) * ThntSW
    NLC(4,NEQ) = LQn2

    ELM(:,:,5,NEQ) =   FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n) * BeamSW
    NLC(5,NEQ) = LQn3

    ! Loss to divertor

    ELM(:,:,6,NEQ) = -          fem_int( 2,rNuL)
    NLC(6,NEQ) = LQe1

    ELM(:,:,7,NEQ) =   PNeDIV * fem_int(-1,rNuL)
    NLC(7,NEQ) = 0

    ! Generated by NBI (Ionization)

    ELM(:,:,8,NEQ) =   (1.D0 - RatCX) * fem_int(-1,SNBe)
    NLC(8,NEQ) = 0

    !  Diffusion of electrons (***AF 2008-06-08)

    ELM(:,:,9,NEQ) = - fem_int(27,sst,DMAGe)
    NLC(9,NEQ) = LQe1

    NLCMAX(NEQ) = 8
    RETURN
  END SUBROUTINE LQe1CC

!***************************************************************
!
!   Electron Radial Flow
!
!***************************************************************

  SUBROUTINE LQe2CC

    integer(4) :: NEQ = LQe2
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(1) * amp)

    ! Pressure gradient force

    ELM(:,:,1,NEQ) = - rKilo   * aee * amasinv * fem_int(5+L2,bri)
    NLC(1,NEQ) = LQe5

    ! Radial E force

    ELM(:,:,2,NEQ) = - achg(1) * aee * amasinv * fem_int(22+L2,bri,Var(:,1)%n)
    NLC(2,NEQ) = LQm1

    ! v x B force

    ELM(:,:,3,NEQ) =   achg(1) * aee * amasinv * fem_int(30+L2,sdt,fipol,Var(:,1)%n)
    NLC(3,NEQ) = LQe3

    ELM(:,:,4,NEQ) = - achg(1) * aee * amasinv * fem_int(20+L2,sdt,bbt)
    NLC(4,NEQ) = LQe4

    NLCMAX(NEQ) = 4

    IF(MDFIXT /= 0) THEN ! Overwrite Pressure gradient force
       ELM(:,:,1,NEQ) = - rKilo * aee * amasinv &
            &                 *(  fem_int(21+L2,bri,Var(:,1)%T) &
            &                   + fem_int(22+L2,bri,Var(:,1)%T))
       NLC(1,NEQ) = LQe1
    END IF

    RETURN
  END SUBROUTINE LQe2CC

!***************************************************************
!
!   Electron Parallel Flow
!
!***************************************************************

  SUBROUTINE LQe3CC

    integer(4) :: NEQ = LQe3
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(1) * amp)

    ELM(:,:, 0,NEQ) =   fem_int(2+L3,Var(:,1)%n) * invDT
    NLC( 0,NEQ) = LQe3

    ! Neoclassical viscosity force

    ELM(:,:, 1,NEQ) = - amasinv * 1.D-20 * fem_int(2+L3,xmu(:,1,1,1))
    NLC( 1,NEQ) = LQe3

    ELM(:,:, 2,NEQ) = - amasinv * 1.D-20 * fem_int(2+L3,xmu(:,1,1,2))
    NLC( 2,NEQ) = LQe6

    ! Diamagnetic forces (particle)

    ELM(:,:, 3,NEQ) = - rKilo / achg(1) * amasinv * fem_int(22+L3,fipolsdtPNsV(:,1),xmu(:,1,1,1))
    NLC( 3,NEQ) = LQe5

    ELM(:,:, 4,NEQ) = - amasinv * 1.D-20 * fem_int(22+L3,fipolsdt,xmu(:,1,1,1))
    NLC( 4,NEQ) = LQm1

    ! Diamagnetic forces (heat)

    ELM(:,:, 5,NEQ) = - rKilo / achg(1) * amasinv * fem_int(22+L3,fipolsdtPNsV(:,1),xmu(:,1,1,2))
    NLC( 5,NEQ) = LQe5

    ELM(:,:, 6,NEQ) =   rKilo / achg(1) * amasinv * fem_int(32+L3,fipolsdtPNsV(:,1),Var(:,1)%T,xmu(:,1,1,2))
    NLC( 6,NEQ) = LQe1
!    ELM(:,:, 6,NEQ) =   rKilo / achg(1) * amasinv * fem_int(30+L3,fipolsdtPNsV(:,1),dlnPNsV(:,1),xmu(:,1,1,2))
!    NLC( 6,NEQ) = LQe5

!!$    ELM(:,:, 5,NEQ) = - rKilo / achg(1) * amasinv * 1.D-20 * fem_int(-(22+L3),fipolsdt,xmu(:,1,1,2),Var(:,1)%T)
!!$!!    ELM(:,:, 5,NEQ) =   amasinv * 1.D-20 * fem_int(-(2+L3),xmu(:,1,1,2),BVsdiag(:,1,2))
!!$    NLC( 5,NEQ) = 0

    ! Collisional friction force with electrons (self collision)

    ELM(:,:, 7,NEQ) =   fem_int(20+L3,lab(:,1,1,1,1),Var(:,1)%n)
    NLC( 7,NEQ) = LQe3

    ! Collisional friction force with ions

    ELM(:,:, 8,NEQ) =   fem_int(20+L3,lab(:,1,2,1,1),Var(:,1)%n)
    NLC( 8,NEQ) = LQi3

    ! Collisional friction force with beam ions

    ELM(:,:, 9,NEQ) =   fem_int(20+L3,laf(:,1,1,1),  RatPNb(:,1))
    NLC( 9,NEQ) = LQb3

    ! Collisional heat friction force with electrons

    ELM(:,:,10,NEQ) = - fem_int(20+L3,lab(:,1,1,1,2),Var(:,1)%n)
    NLC(10,NEQ) = LQe6

    ! Collisional heat friction force with ions

    ELM(:,:,11,NEQ) = - fem_int(20+L3,lab(:,1,2,1,2),Var(:,1)%n)
    NLC(11,NEQ) = LQi6

    ! Electric field forces

    ELM(:,:,12,NEQ) = - achg(1) * aee * amasinv * fem_int(30+L3,fipol,Var(:,1)%n,aatq)
    NLC(12,NEQ) = LQm2

    ELM(:,:,13,NEQ) =   achg(1) * aee * amasinv * fem_int(30+L3,fipol,Var(:,1)%n,aat)
    NLC(13,NEQ) = LQm3

    ! Loss to divertor

    ELM(:,:,14,NEQ) = - fem_int(20+L3,rNuL,Var(:,1)%n) * fact
    NLC(14,NEQ) = LQe3

    ! Collisional friction force with neutrals

    ELM(:,:,15,NEQ) = - fem_int(20+L3,rNu0e,Var(:,1)%n)
    NLC(15,NEQ) = LQe3

!!$    ! Perpendicular viscosity (off-diagonal)
!!$
!!$    ELM(:,:,16,NEQ) = - FSPARV(1) * fem_int(37,ParatoZeta,sst,rMue)
!!$    NLC(16,NEQ) = LQe4
!!$
!!$    ELM(:,:,17,NEQ) =   FSPARV(1) * fem_int(37,ParatoZeta,sst*rMue,Var(:,1)%RUph)
!!$    NLC(17,NEQ) = LQe1
!!$  
!!$    ELM(:,:,18,NEQ) = - FSPARV(1) * fem_int(25,sst*rMue,ParatoZeta)
!!$    NLC(18,NEQ) = LQe4
!!$
!!$    ELM(:,:,19,NEQ) =   FSPARV(1) * fem_int(25,sst*rMue*Var(:,1)%RUph,ParatoZeta)
!!$    NLC(19,NEQ) = LQe1
!!$
!!$    ! Momentum pinch (off-diagonal)
!!$
!!$    ELM(:,:,20,NEQ) = - FSPARV(1) * fem_int(24,suft*Vmps(:,1),ZetatoPara)
!!$    NLC(20,NEQ) = LQe4
!!$
!!$    ELM(:,:,21,NEQ) =   FSPARV(1) * fem_int(21,suft*Vmps(:,1),ZetatoPara)
!!$    NLC(21,NEQ) = LQe4

    ! Perpendicular viscosity (diagonal)

    ELM(:,:,16,NEQ) = - FSPARV(1) * fem_int(37,sst,rMue,Var(:,1)%n)
    NLC(16,NEQ) = LQe3

    ELM(:,:,17,NEQ) =   FSPARV(1) * fem_int(24,sst*rMue*Var(:,1)%n,dlogZetatoPara)
    NLC(17,NEQ) = LQe3
  
    ELM(:,:,18,NEQ) = - FSPARV(1) * fem_int(32,sst*rMue,Var(:,1)%n,dlogZetatoPara)
    NLC(18,NEQ) = LQe3

    ELM(:,:,19,NEQ) =   FSPARV(1) * fem_int(30,sst*rMue,Var(:,1)%n,dlogZetatoPara**2)
    NLC(19,NEQ) = LQe3

    ! Momentum pinch (diagonal)

    ELM(:,:,20,NEQ) = - FSPARV(1) * fem_int(24,suft*Vmps(:,1),Var(:,1)%n)
    NLC(20,NEQ) = LQe3

    ELM(:,:,21,NEQ) =   FSPARV(1) * fem_int(30,suft*Vmps(:,1),Var(:,1)%n,dlogZetatoPara)
    NLC(21,NEQ) = LQe3

    !  Diffusion of electrons (***AF 2008-06-08) (37+L3 is not possible because it is a diffusion term)

    ELM(:,:,22,NEQ) = - fem_int(37,sst,DMAGe,Var(:,1)%n)
    NLC(22,NEQ) = LQe3

    !  Additional torque input

!    ELM(:,:,17,NEQ) = - amasinv * 1.D-20 * fem_int(-(2+L3),Tqt,ZetatoPara)
!    NLC(17,NEQ) = 0

!    ! Helical neoclassical viscosity force (***AF 2008-06-08)
!
!    ELM(:,:,17,NEQ) = - fem_int(20+L3,rNueHLthth,Var(:,1)%n)
!    NLC(17,NEQ) = LQe3
!
!!---- 09/11/26 AF
!    ! To deal with the existence of (m=0, n>0) component   miki_m  2010/8/26
!!    ELM(:,:,18,NEQ) =   fem_int(15,rNueHLthph)
!    IF(UHphSwitch == 0) THEN
!       ELM(:,:,18,NEQ) =   fem_int(15,rNueHLthph)
!    else
!       ELM(:,:,18,NEQ) =   fem_int(22,rNueHLthph)
!    endif
!    NLC(18,NEQ) = LQe4

    NLCMAX(NEQ) = 21

    IF(MDFIXT /= 0) THEN ! Overwrite Pressure gradient force
       ELM(:,:, 3,NEQ) = - rKilo / achg(1) * amasinv &
            &                  *(  fem_int(31+L3,fipolsdtPNsV(:,1),xmu(:,1,1,1),Var(:,1)%T) &
            &                    + fem_int(32+L3,fipolsdtPNsV(:,1),xmu(:,1,1,1),Var(:,1)%T))
       NLC( 3,NEQ) = LQe1

       ELM(:,:, 5,NEQ) =   amasinv * 1.D-20 * fem_int(-(2+L3),xmu(:,1,1,2),BVsdiag(:,1,2))
       NLC( 5,NEQ) = 0

       ELM(:,:, 6,NEQ) =   0.d0
       NLC( 6,NEQ) = 0
    END IF

    RETURN
  END SUBROUTINE LQe3CC

!***************************************************************
!
!   Electron Toroidal Flow
!
!***************************************************************

  SUBROUTINE LQe4CC

    integer(4) :: NEQ = LQe4
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(1) * amp)

    ! Uephi(0)' : 0

    ELM(:,:, 0,NEQ) = fem_int(1) * invDT
    NLC( 0,NEQ) = LQe4

    ! Advection + Grid velocity

    ELM(:,:, 1,NEQ) = - fem_int(5,UsrgV(:,1))
    NLC( 1,NEQ) = LQe4

    ELM(:,:, 2,NEQ) = - FSADV * fem_int(6,Var(:,1)%UrV)
    NLC( 2,NEQ) = LQe4

    ! Momentum pinch

    ELM(:,:, 3,NEQ) = - fem_int(24,suft,Vmps(:,1))
    NLC( 3,NEQ) = LQe4

    ! Viscosity force

    ELM(:,:, 4,NEQ) = - fem_int(27,sst,rMue)
    NLC( 4,NEQ) = LQe4

    ELM(:,:, 5,NEQ) =   fem_int(37,sst,rMue,Var(:,1)%RUph)
    NLC( 5,NEQ) = LQe1

    ! Residual stress

    ELM(:,:, 6,NEQ) = - amasinv * fem_int(-4,PiRess(:,1))
    NLC( 6,NEQ) = 0

    ! Collisional friction force with electrons (self collision)

    ELM(:,:, 7,NEQ) =   fem_int(2,lab(:,1,1,1,1))
    NLC( 7,NEQ) = LQe4

    ! Collisional friction force with bulk ions

    ELM(:,:, 8,NEQ) =   fem_int(20,lab(:,1,2,1,1),RatPN(:,1))
    NLC( 8,NEQ) = LQi4

    ! Collisional friction force with beam ions

!    ELM(:,:, 9,NEQ) =   fem_int(30,ParatoZeta,laf(:,1,1,1),RatPNb(:,1))
!    NLC( 9,NEQ) = LQb3

    ELM(:,:, 9,NEQ) =   fem_int(20,laf(:,1,1,1),RatPNb(:,1))
    NLC( 9,NEQ) = LQb4

    ! Collisional heat friction force with electrons

    ELM(:,:,10,NEQ) = - fem_int(30,ParatoZeta,lab(:,1,1,1,2),Var(:,1)%n)
    NLC(10,NEQ) = LQe6

    ! Collisional heat friction force with ions

    ELM(:,:,11,NEQ) = - fem_int(30,ParatoZeta,lab(:,1,2,1,2),Var(:,1)%n)
    NLC(11,NEQ) = LQi6

    ! Toroidal E force

    ELM(:,:,12,NEQ) =   achg(1) * aee * amasinv * fem_int(2,Var(:,1)%n)
    NLC(12,NEQ) = LQm3

    ! v x B force

    ELM(:,:,13,NEQ) =   achg(1) * aee * amasinv * fem_int(6,PsiV)
    NLC(13,NEQ) = LQe2

    ! Turbulent particle transport driver

    ELM(:,:,14,NEQ) =   achg(1) * aee * amasinv * fem_int( 2,FQLcoef1)
    NLC(14,NEQ) = LQe4

    ELM(:,:,15,NEQ) = - achg(1) * aee * amasinv * fem_int(20,FQLcoef2,Var(:,1)%n)
    NLC(15,NEQ) = LQe3

    ELM(:,:,16,NEQ) =   achg(1) * rKeV * amasinv * fem_int( 5,FQLcoef) &
         &                           * (FSVAHL - 1.d0)
    NLC(16,NEQ) = LQe5

    ELM(:,:,17,NEQ) = - achg(1) * rKeV * amasinv * fem_int(22,FQLcoef,Var(:,1)%T) &
         &                           * (FSVAHL - 1.d0)
    NLC(17,NEQ) = LQe1

    ELM(:,:,18,NEQ) =   achg(1) * aee * amasinv * fem_int(22,FQLcoef,Var(:,1)%n)
    NLC(18,NEQ) = LQm1

    !  -- Turbulent pinch term

    ELM(:,:,19,NEQ) =   achg(1) * aee * amasinv * fem_int(20,FQLcoef,FVpch)
    NLC(19,NEQ) = LQe1

    ! Loss to divertor

    ELM(:,:,20,NEQ) = - fem_int(2,rNuL) * fact
    NLC(20,NEQ) = LQe4

    ! Collisional friction force with neutrals

    ELM(:,:,21,NEQ) = - fem_int(2,rNu0e)
    NLC(21,NEQ) = LQe4

    ! Diffusion of electrons (***AF 2008-06-08)

    ELM(:,:,22,NEQ) = - fem_int(27,sst,DMAGe)
    NLC(22,NEQ) = LQe4

    ! Additional torque input

!    ELM(:,:,23,NEQ) = - amasinv * 1.D-20 * fem_int(-1,Tqt) 
!    NLC(23,NEQ) = 0

!    ! Helical neoclassical viscosity force (***AF 2008-06-08)
!
!    ! To deal with the existence of (m=0, n>0) component   miki_m  2010/8/26
!!    ELM(:,:,24,NEQ) =   fem_int(2,rNueHLphth)
!    IF(UHphSwitch == 0) THEN
!       ELM(:,:,24,NEQ) =   fem_int(2,rNueHLphth) * rr
!    ELSE
!       ELM(:,:,24,NEQ) =   fem_int(2,rNueHLphthVR) * rr
!       ELM(1      ,1:4,24,NEQ) =   fem_int_point(2,0,rNueHLphthVR) * rr
!    ENDIF
!    NLC(24,NEQ) = LQe3
!
!!    ELM(:,:,25,NEQ) = - fem_int(15,rNueHLphph)
!    IF(UHphSwitch == 0) THEN
!       ELM(:,:,25,NEQ) = - fem_int(15,rNueHLphph) 
!    else
!       ELM(:,:,25,NEQ) = - fem_int(2,rNueHLphph)
!    endif
!    NLC(25,NEQ) = LQe4
 
    NLCMAX(NEQ) = 21

    IF(MDFIXT /= 0) THEN ! Overwrite Pressure gradient force
       ELM(:,:,16,NEQ) =   achg(1) * rKeV * amasinv * fem_int(21,FQLcoef,Var(:,1)%T) &
         &                           * (FSVAHL - 1.d0)
       NLC(16,NEQ) = LQe1

       ELM(:,:,17,NEQ) = 0.d0
       NLC(17,NEQ) = 0
    END IF

    RETURN
  END SUBROUTINE LQe4CC

!***************************************************************
!
!   Electron Energy Transport: Te
!
!***************************************************************

  SUBROUTINE LQe5CC

    integer(4) :: NEQ = LQe5

    ! Temperature evolution
    
    IF(MDFIXT == 0) THEN
       ELM(:,:, 0,NEQ) =   1.5D0 * fem_int(1) * invDT
       NLC( 0,NEQ) = LQe5

       ! Advection + Grid velocity

       ELM(:,:, 1,NEQ) = - 2.5d0 * fem_int(5,UsrgV(:,1))
       NLC( 1,NEQ) = LQe5

       ELM(:,:, 2,NEQ) = - 2.5d0 * FSADV * fem_int(6,Var(:,1)%UrV)
       NLC( 2,NEQ) = LQe5

       ! Heat pinch

       ELM(:,:, 3,NEQ) = - fem_int(24,Vhps(:,1),suft)
       NLC( 3,NEQ) = LQe5

       ! Conduction transport

       ELM(:,:, 4,NEQ) = - fem_int(27,sst,Chie1)
       NLC( 4,NEQ) = LQe5

       ELM(:,:, 5,NEQ) =   fem_int(37,sst,Chie2,Var(:,1)%T)
       NLC( 5,NEQ) = LQe1

       ! Redundant heat advection

       ELM(:,:, 6,NEQ) = - fem_int(5,UgV)
       NLC( 6,NEQ) = LQe5

       ELM(:,:, 7,NEQ) = - fem_int(5,Var(:,2)%UrV)
       NLC( 7,NEQ) = LQi5

       ! Viscous heating

       ELM(:,:, 8,NEQ) = - 1.D-20 / rKeV * fem_int(-2,Var(:,2)%Uthhat,BnablaPi(:,2))
       NLC( 8,NEQ) = 0

       ! Collisional transfer with ions (Energy equilibration)

       ELM(:,:, 9,NEQ) = - 1.5d0 * fem_int( 2,rNuTei)
       NLC( 9,NEQ) = LQe5

       ELM(:,:,10,NEQ) =   1.5d0 * fem_int(20,rNuTei,RatPN(:,1))
       NLC(10,NEQ) = LQi5

       ! Joule heating

       ELM(:,:,11,NEQ) = - achg(1) * 1.d-3 * fem_int(-20,PsitdotVBth,Var(:,1)%n,Var(:,1)%Uthhat)
       NLC(11,NEQ) = 0

       ELM(:,:,12,NEQ) = - achg(2) * 1.d-3 * fem_int(-20,PsitdotVBth,Var(:,2)%n,Var(:,2)%Uthhat)
       NLC(12,NEQ) = 0

       ELM(:,:,13,NEQ) =   achg(1) * 1.d-3 * fem_int(  2,PsidotV)
       NLC(13,NEQ) = LQe7

       ELM(:,:,14,NEQ) =   achg(2) * 1.d-3 * fem_int(  2,PsidotV)
       NLC(14,NEQ) = LQi7
!       ELM(:,:,13,NEQ) =   achg(1) * 1.d-3 * fem_int(-20,PsidotV,Var(:,1)%n,Var(:,1)%UphR)
!       NLC(13,NEQ) = 0
!
!       ELM(:,:,14,NEQ) =   achg(2) * 1.d-3 * fem_int(-20,PsidotV,Var(:,2)%n,Var(:,2)%UphR)
!       NLC(14,NEQ) = 0

       ! Loss to diverter

       ELM(:,:,15,NEQ) = -                  fem_int( 2,rNuL)
       NLC(15,NEQ) = LQe5

       ELM(:,:,16,NEQ) =           PNeDIV * fem_int(-2,rNuL,Var(:,1)%T)
       NLC(16,NEQ) = 0

       ELM(:,:,17,NEQ) = - 1.5D0          * fem_int( 2,rNuLTe)
       NLC(17,NEQ) = LQe5

       ELM(:,:,18,NEQ) =   1.5D0 * PTeDIV * fem_int(-2,rNuLTe,Var(:,1)%n)
       NLC(18,NEQ) = 0
!!$       ELM(:,:,18,NEQ) =   1.5D0 * PTeDIV * fem_int( 2,rNuLTe)
!!$       NLC(18,NEQ) = LQe1

       ! Ionization loss of n01, n02 and n03

       ELM(:,:,19,NEQ) = - (EION * 1.D-3) * FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n)
       NLC(19,NEQ) = LQn1

       ELM(:,:,20,NEQ) = - (EION * 1.D-3) * FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n)
       NLC(20,NEQ) = LQn2

       ELM(:,:,21,NEQ) = - (EION * 1.D-3) * FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n) * BeamSW
       NLC(21,NEQ) = LQn3

       ! Collisional NBI heating (Perp + Tan)

       ELM(:,:,22,NEQ) = Eb * fem_int(-2,SNB,PNBcol_e)
       NLC(22,NEQ) = 0

       ! Simpified Alpha heating

       ELM(:,:,23,NEQ) = 1.D-20 / rKeV * fem_int(-1,PALFe)
       NLC(23,NEQ) = 0
       
       ! Direct heating (RF)

       ELM(:,:,24,NEQ) =   1.D-20 / rKeV * fem_int(-1,PRFe)
       NLC(24,NEQ) = 0

       ! Radiation loss (Bremsstrahlung)

       ELM(:,:,25,NEQ) = - 1.D-20 / rKeV * fem_int(-1,PBr)
       NLC(25,NEQ) = 0

       !  Diffusion of electrons (***AF 2008-06-08)

       ELM(:,:,26,NEQ) = - fem_int(27,sst,DMAGe)
       NLC(26,NEQ) = LQe5

!       ! Collisional heating with beam
!
!       ELM(:,:,27,NEQ) = achgb * amqp * 1.d-3 * fem_int(-20,MNB,BUsparVbbt(:,1),PNbVinv)
!       NLC(27,NEQ) = LQb3

       NLCMAX(NEQ) = 25
    ELSE

       !  Fixed temperature profile

       ELM(:,:,0,NEQ) = fem_int(1) * invDT
       NLC(0,NEQ) = LQe5

       NLCMAX(NEQ) = 0
    END IF

    RETURN
  END SUBROUTINE LQe5CC

!***************************************************************
!
!   Electron Heat Parallel Flow
!
!***************************************************************

  SUBROUTINE LQe6CC

    integer(4) :: NEQ = LQe6
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(1) * amp)

    ELM(:,:, 0,NEQ) =   2.5d0 * fem_int(2+L6,Var(:,1)%n) * invDT
    NLC( 0,NEQ) = LQe6

    ! Neoclassical viscosity force

    ELM(:,:, 1,NEQ) = - amasinv * 1.D-20 * fem_int(2+L6,xmu(:,1,2,1))
    NLC( 1,NEQ) = LQe3

    ELM(:,:, 2,NEQ) = - amasinv * 1.D-20 * fem_int(2+L6,xmu(:,1,2,2))
    NLC( 2,NEQ) = LQe6

    ! Diamagnetic forces

    ELM(:,:, 3,NEQ) = - rKilo / achg(1) * amasinv * fem_int(22+L6,fipolsdtPNsV(:,1),xmu(:,1,2,1))
    NLC( 3,NEQ) = LQe5

    ELM(:,:, 4,NEQ) = - amasinv * 1.D-20 * fem_int(22+L6,fipolsdt,xmu(:,1,2,1))
    NLC( 4,NEQ) = LQm1

    ELM(:,:, 5,NEQ) = - rKilo / achg(1) * amasinv * fem_int(22+L6,fipolsdtPNsV(:,1),xmu(:,1,2,2))
    NLC( 5,NEQ) = LQe5

    ELM(:,:, 6,NEQ) =   rKilo / achg(1) * amasinv * fem_int(32+L6,fipolsdtPNsV(:,1),Var(:,1)%T,xmu(:,1,2,2))
    NLC( 6,NEQ) = LQe1

    ! Collisional friction force with electrons

    ELM(:,:, 7,NEQ) = - fem_int(20+L6,lab(:,1,1,2,1),Var(:,1)%n)
    NLC( 7,NEQ) = LQe3

    ! Collisional friction force with ions

    ELM(:,:, 8,NEQ) = - fem_int(20+L6,lab(:,1,2,2,1),Var(:,1)%n)
    NLC( 8,NEQ) = LQi3

    ! Collisional friction force with beam ions

    ELM(:,:, 9,NEQ) = - fem_int(20+L6,laf(:,1,2,1),  RatPNb(:,1))
    NLC( 9,NEQ) = LQb3

    ! Collisional heat friction force with electrons (self collision)

    ELM(:,:,10,NEQ) =   fem_int(20+L6,lab(:,1,1,2,2),Var(:,1)%n)
    NLC(10,NEQ) = LQe6

    ! Collisional heat friction force with ions

    ELM(:,:,11,NEQ) =   fem_int(20+L6,lab(:,1,2,2,2),Var(:,1)%n)
    NLC(11,NEQ) = LQi6

    ! Loss to divertor

    ELM(:,:,12,NEQ) = - fem_int(20+L6,rNuL,Var(:,1)%n) * fact
    NLC(12,NEQ) = LQe6

    ! Collisional friction force with neutrals

    ELM(:,:,13,NEQ) = - fem_int(20+L6,rNu0e,Var(:,1)%n)
    NLC(13,NEQ) = LQe6

    NLCMAX(NEQ) = 13

    IF(MDFIXT /= 0) THEN ! Overwrite Pressure gradient force
       ELM(:,:, 3,NEQ) = - rKilo / achg(1) * amasinv &
            &                  *(  fem_int(31+L6,fipolsdtPNsV(:,1),xmu(:,1,2,1),Var(:,1)%T) &
            &                    + fem_int(22+L6,fipolsdtPNsV(:,1),xmu(:,1,2,1),Var(:,1)%T))
       NLC( 3,NEQ) = LQe1

       ELM(:,:, 5,NEQ) =   amasinv * 1.D-20 * fem_int(-(2+L6),xmu(:,1,2,2),BVsdiag(:,1,2))
       NLC( 5,NEQ) = 0

       ELM(:,:, 6,NEQ) =   0.d0
       NLC( 6,NEQ) = 0
    END IF

    RETURN
  END SUBROUTINE LQe6CC

!***************************************************************
!
!   Electron rotation frequency: Ne <UePhi/R>
!
!***************************************************************

  SUBROUTINE LQe7CC

    integer(4) :: NEQ = LQe7

    !  Electron rotation frequency

    ELM(:,:,1,NEQ) = - fem_int(1)
    NLC(1,NEQ) = LQe7

    !  Electron rotation

    ELM(:,:,2,NEQ) =   fem_int(2,rrtinv)
    NLC(2,NEQ) = LQe4

    ELM(:,:,3,NEQ) = - fem_int(20,Fqhatsq,rrtinv)
    NLC(3,NEQ) = LQe4

    !  Electron parallel velocity

    ELM(:,:,4,NEQ) =   fem_int(20,FqhatsqI,Var(:,1)%n)
    NLC(4,NEQ) = LQe3

    NLCMAX(NEQ) = 4
    RETURN
  END SUBROUTINE LQe7CC

!***************************************************************
!
!   Ion Density Equation
!
!***************************************************************

  SUBROUTINE LQi1CC

    integer(4) :: NEQ = LQi1

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQi1

    ! Divergence

    ELM(:,:,1,NEQ) = - fem_int(4)
    NLC(1,NEQ) = LQi2

    ! Divergence + Grid velocity

    ELM(:,:,2,NEQ) =   fem_int(5,UgV)
    NLC(2,NEQ) = LQi1

    ! Ionization of n01, n02 and n03

    ELM(:,:,3,NEQ) =     FSION * 1.D20 * fem_int(20,SiVizA,Var(:,2)%n)
    NLC(3,NEQ) = LQn1

    ELM(:,:,4,NEQ) =     FSION * 1.D20 * fem_int(20,SiVizA,Var(:,2)%n) * ThntSW
    NLC(4,NEQ) = LQn2

    ELM(:,:,5,NEQ) =     FSION * 1.D20 * fem_int(20,SiVizA,Var(:,2)%n) * BeamSW
    NLC(5,NEQ) = LQn3

!!$    ELM(:,:,3,NEQ) =     1.D0 / achg(2)  * FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n)
!!$    NLC(3,NEQ) = LQn1
!!$
!!$    ELM(:,:,4,NEQ) =     1.D0 / achg(2)  * FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n) * ThntSW
!!$    NLC(4,NEQ) = LQn2
!!$
!!$    ELM(:,:,5,NEQ) =     1.D0 / achg(2)  * FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n) * BeamSW
!!$    NLC(5,NEQ) = LQn3

    ! Loss to divertor

    ELM(:,:,6,NEQ) = -          fem_int( 2,rNuL)
    NLC(6,NEQ) = LQi1

    ELM(:,:,7,NEQ) =   PNiDIV * fem_int(-1,rNuL)
    NLC(7,NEQ) = 0

!!$    ELM(:,:,6,NEQ) = -   1.D0 / achg(2) * fem_int( 2,rNuL)
!!$    NLC(6,NEQ) = LQe1
!!$
!!$    ELM(:,:,7,NEQ) =   PNeDIV / achg(2) * fem_int(-1,rNuL)
!!$    NLC(7,NEQ) = 0

    ! Particle source from beam ion

    ELM(:,:,8,NEQ) =   fem_int(2,rNuB)
    NLC(8,NEQ) = LQb1

    ! Particle source from ripple trapped beam ions

    ELM(:,:,9,NEQ) =   fem_int(20,rNuB,rip_rat) * RpplSW
    NLC(9,NEQ) = LQr1

    ! NBI kick up ions (Charge exchange)

    ELM(:,:,10,NEQ) = - RatCX * fem_int(-1,SNBi)
    NLC(10,NEQ) = 0

    ! Loss cone loss

    ELM(:,:,11,NEQ) =   fem_int(-1,SiLC)
    NLC(11,NEQ) = 0
 
    ! Ion orbit loss

    ELM(:,:,12,NEQ) = - fem_int(2,rNuOL)
    NLC(12,NEQ) = LQi1

    !  Diffusion of ions (***AF 2008-06-08)

    ELM(:,:,13,NEQ) = - fem_int(27,sst,DMAGi)
    NLC(13,NEQ) = LQi1

!    ! Parallel loss reduction due to the potential
!    ! induced by the parallel loss of the beam ions
!
!    ELM(:,:,14,NEQ) =   fem_int(2,rNuLB) * BeamSW
!    NLC(14,NEQ) = LQb1

    NLCMAX(NEQ) = 12
    RETURN
  END SUBROUTINE LQi1CC

!***************************************************************
!
!   Ion Radial Flow
!
!     Note: Ion pressure gradient term may destabilize the neutrality at the edge at initial.
!
!***************************************************************
  
  SUBROUTINE LQi2CC

    integer(4) :: NEQ = LQi2
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(2) * amp)

    ! Pressure gradient force

    ELM(:,:,1,NEQ) = - rKilo   * aee * amasinv * fem_int(5+L2,bri)
    NLC(1,NEQ) = LQi5

    ! Radial E force

    ELM(:,:,2,NEQ) = - achg(2) * aee * amasinv * fem_int(22+L2,bri,Var(:,2)%n)
    NLC(2,NEQ) = LQm1

    ! v x B force

    ELM(:,:,3,NEQ) =   achg(2) * aee * amasinv * fem_int(30+L2,sdt,fipol,Var(:,2)%n)
    NLC(3,NEQ) = LQi3

    ELM(:,:,4,NEQ) = - achg(2) * aee * amasinv * fem_int(20+L2,sdt,bbt)
    NLC(4,NEQ) = LQi4

    NLCMAX(NEQ) = 4

    IF(MDFIXT /= 0) THEN ! Overwrite Pressure gradient force
       ELM(:,:,1,NEQ) = - rKilo * aee * amasinv &
            &                 *(  fem_int(21+L2,bri,Var(:,2)%T) &
            &                   + fem_int(22+L2,bri,Var(:,2)%T))
       NLC(1,NEQ) = LQi1
    END IF

    RETURN
  END SUBROUTINE LQi2CC

!***************************************************************
!
!   Ion Parallel Flow
!
!***************************************************************

  SUBROUTINE LQi3CC

    integer(4) :: NEQ = LQi3
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(2) * amp)

    ELM(:,:, 0,NEQ) =   fem_int(2+L3,Var(:,2)%n) * invDT
    NLC( 0,NEQ) = LQi3

    ! Neoclassical viscosity force

    ELM(:,:, 1,NEQ) = - amasinv * 1.D-20 * fem_int(2+L3,xmu(:,2,1,1))
    NLC( 1,NEQ) = LQi3

    ELM(:,:, 2,NEQ) = - amasinv * 1.D-20 * fem_int(2+L3,xmu(:,2,1,2))
    NLC( 2,NEQ) = LQi6

    ! Diamagnetic forces (particle)

    ELM(:,:, 3,NEQ) = - rKilo / achg(2) * amasinv * fem_int(22+L3,fipolsdtPNsV(:,2),xmu(:,2,1,1))
    NLC( 3,NEQ) = LQi5

    ELM(:,:, 4,NEQ) = - amasinv * 1.D-20 * fem_int(22+L3,fipolsdt,xmu(:,2,1,1))
    NLC( 4,NEQ) = LQm1

    ! Diamagnetic forces (heat)

    ELM(:,:, 5,NEQ) = - rKilo / achg(2) * amasinv * fem_int(22+L3,fipolsdtPNsV(:,2),xmu(:,2,1,2))
    NLC( 5,NEQ) = LQi5

    ELM(:,:, 6,NEQ) =   rKilo / achg(2) * amasinv * fem_int(32+L3,fipolsdtPNsV(:,2),Var(:,2)%T,xmu(:,2,1,2))
    NLC( 6,NEQ) = LQi1
!    ELM(:,:, 6,NEQ) =   rKilo / achg(2) * amasinv * fem_int(30+L3,fipolsdtPNsV(:,2),dlnPNsV(:,2),xmu(:,2,1,2))
!    NLC( 6,NEQ) = LQi5

!!$    ELM(:,:, 5,NEQ) = - rKilo / achg(2) * amasinv * 1.D-20 * fem_int(-(22+L3),fipolsdt,xmu(:,2,1,2),Var(:,2)%T)
!!$!!    ELM(:,:, 5,NEQ) =   amasinv * 1.D-20 * fem_int(-(2+L3),xmu(:,2,1,2),BVsdiag(:,2,2))
!!$    NLC( 5,NEQ) = 0

    ! Collisional friction force with electrons

    ELM(:,:, 7,NEQ) =   fem_int(20+L3,lab(:,2,1,1,1),Var(:,2)%n)
    NLC( 7,NEQ) = LQe3

    ! Collisional friction force with ions (self collision)

    ELM(:,:, 8,NEQ) =   fem_int(20+L3,lab(:,2,2,1,1),Var(:,2)%n)
    NLC( 8,NEQ) = LQi3

    ! Collisional friction force with beam ions

    ELM(:,:, 9,NEQ) =   fem_int(20+L3,laf(:,2,1,1),  RatPNb(:,2))
    NLC( 9,NEQ) = LQb3

    ! Collisional heat friction force with electrons

    ELM(:,:,10,NEQ) = - fem_int(20+L3,lab(:,2,1,1,2),Var(:,2)%n)
    NLC(10,NEQ) = LQe6

    ! Collisional heat friction force with ions

    ELM(:,:,11,NEQ) = - fem_int(20+L3,lab(:,2,2,1,2),Var(:,2)%n)
    NLC(11,NEQ) = LQi6

    ! Electric field forces

    ELM(:,:,12,NEQ) = - achg(2) * aee * amasinv * fem_int(30+L3,fipol,Var(:,2)%n,aatq)
    NLC(12,NEQ) = LQm2

    ELM(:,:,13,NEQ) =   achg(2) * aee * amasinv * fem_int(30+L3,fipol,Var(:,2)%n,aat)
    NLC(13,NEQ) = LQm3

    ! Loss to divertor

    ELM(:,:,14,NEQ) = - fem_int(20+L3,rNuL,Var(:,2)%n) * fact
    NLC(14,NEQ) = LQi3

    ! Increase in momentum due to thermalization of fast ions

    ELM(:,:,15,NEQ) =   amb / amas(2) * fem_int(2+L3,rNuB) * BeamSW
    NLC(15,NEQ) = LQb3

    ! Collisional friction force with neutrals

    ELM(:,:,16,NEQ) = - fem_int(20+L3,rNu0i,Var(:,2)%n)
    NLC(16,NEQ) = LQi3

    ! Charge exchange force

    ELM(:,:,17,NEQ) = - fem_int(20+L3,rNuiCX,Var(:,2)%n)
    NLC(17,NEQ) = LQi3

    ! Loss cone loss

    ELM(:,:,18,NEQ) =   fem_int(-(2+L3),SiLCth,Var(:,2)%n)
    NLC(18,NEQ) = 0

    ! Ion orbit loss

    ELM(:,:,19,NEQ) = - fem_int(20+L3,rNuOL,Var(:,2)%n)
    NLC(19,NEQ) = LQi3

    !  Additional torque input (sheared flow)

    ELM(:,:,20,NEQ) =   amasinv * 1.D-20 * fem_int(-(1+L3),Tqp)
    NLC(20,NEQ) = 0

    !  Additional torque input

    ELM(:,:,21,NEQ) =   amasinv * 1.D-20 * fem_int(-(2+L3),Tqt,ZetatoPara)
    NLC(21,NEQ) = 0

!!$    ! Perpendicular viscosity (off-diagonal)
!!$
!!$    ELM(:,:,22,NEQ) = - FSPARV(2) * fem_int(37,ParatoZeta,sst,rMui)
!!$    NLC(22,NEQ) = LQi4
!!$
!!$    ELM(:,:,23,NEQ) =   FSPARV(2) * fem_int(37,ParatoZeta,sst*rMui,Var(:,2)%RUph)
!!$    NLC(23,NEQ) = LQi1
!!$  
!!$    ELM(:,:,24,NEQ) = - FSPARV(2) * fem_int(25,sst*rMui,ParatoZeta)
!!$    NLC(24,NEQ) = LQi4
!!$
!!$    ELM(:,:,25,NEQ) =   FSPARV(2) * fem_int(25,sst*rMui*Var(:,2)%RUph,ParatoZeta)
!!$    NLC(25,NEQ) = LQi1
!!$
!!$    ! Momentum pinch (off-diagonal)
!!$
!!$    ELM(:,:,26,NEQ) = - FSPARV(2) * fem_int(24,suft*Vmps(:,2),ZetatoPara)
!!$    NLC(26,NEQ) = LQi4
!!$
!!$    ELM(:,:,27,NEQ) =   FSPARV(2) * fem_int(21,suft*Vmps(:,2),ZetatoPara)
!!$    NLC(27,NEQ) = LQi4

    ! Perpendicular viscosity (diagonal)

    ELM(:,:,22,NEQ) = - FSPARV(2) * fem_int(37,sst,rMui,Var(:,2)%n)
    NLC(22,NEQ) = LQi3

    ELM(:,:,23,NEQ) =   FSPARV(2) * fem_int(24,sst*rMui*Var(:,2)%n,dlogZetatoPara)
    NLC(23,NEQ) = LQi3
  
    ELM(:,:,24,NEQ) = - FSPARV(2) * fem_int(32,sst*rMui,Var(:,2)%n,dlogZetatoPara)
    NLC(24,NEQ) = LQi3

    ELM(:,:,25,NEQ) =   FSPARV(2) * fem_int(30,sst*rMui,Var(:,2)%n,dlogZetatoPara**2)
    NLC(25,NEQ) = LQi3

    ! Momentum pinch (diagonal)

    ELM(:,:,26,NEQ) = - FSPARV(2) * fem_int(24,suft*Vmps(:,2),Var(:,2)%n)
    NLC(26,NEQ) = LQi3

    ELM(:,:,27,NEQ) =   FSPARV(2) * fem_int(30,suft*Vmps(:,2),Var(:,2)%n,dlogZetatoPara)
    NLC(27,NEQ) = LQi3

    !  Diffusion of ions (***AF 2008-06-08) (37+L3 is not possible because it is a diffusion term)

    ELM(:,:,28,NEQ) = - fem_int(37,sst,DMAGi,Var(:,2)%n)
    NLC(28,NEQ) = LQi3

!    ! Helical Neoclassical viscosity force
!
!    ELM(:,:,23,NEQ) = - fem_int(20+L3,rNuiHLthth,Var(:,2)%n)
!    NLC(23,NEQ) = LQi3
!
!    ! To deal with the existence of (m=0, n>0) component   miki_m  2010/8/26
!!    ELM(:,:,24,NEQ) =   fem_int(15,rNuiHLthph)
!    IF(UHphSwitch == 0) THEN
!       ELM(:,:,24,NEQ) =   fem_int(15,rNuiHLthph)
!    else
!       ELM(:,:,24,NEQ) =   fem_int(22,rNuiHLthph)
!    endif
!    NLC(24,NEQ) = LQi4

    NLCMAX(NEQ) = 27

    IF(MDFIXT /= 0) THEN ! Overwrite Pressure gradient force
       ELM(:,:, 3,NEQ) = - rKilo / achg(2) * amasinv &
            &                  *(  fem_int(31+L3,fipolsdtPNsV(:,2),xmu(:,2,1,1),Var(:,2)%T) &
            &                    + fem_int(32+L3,fipolsdtPNsV(:,2),xmu(:,2,1,1),Var(:,2)%T))
       NLC( 3,NEQ) = LQi1

       ELM(:,:, 5,NEQ) =   amasinv * 1.D-20 * fem_int(-(2+L3),xmu(:,2,1,2),BVsdiag(:,2,2))
       NLC( 5,NEQ) = 0

       ELM(:,:, 6,NEQ) =   0.d0
       NLC( 6,NEQ) = 0
    END IF

    RETURN
  END SUBROUTINE LQi3CC

!***************************************************************
!
!   Ion Toroidal Flow
!
!***************************************************************

  SUBROUTINE LQi4CC

    integer(4) :: NEQ = LQi4
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(2) * amp)

    ! Uiphi'(0) : 0

    ELM(:,:, 0,NEQ) = fem_int(1) * invDT
    NLC( 0,NEQ) = LQi4

    ! Advection + Grid velocity

    ELM(:,:, 1,NEQ) = - fem_int(5,UsrgV(:,2))
    NLC( 1,NEQ) = LQi4

    ELM(:,:, 2,NEQ) = - FSADV * fem_int(6,Var(:,2)%UrV)
    NLC( 2,NEQ) = LQi4

    ! Momentum pinch

    ELM(:,:, 3,NEQ) = - fem_int(24,suft,Vmps(:,2))
    NLC( 3,NEQ) = LQi4

    ! Viscosity force

    ELM(:,:, 4,NEQ) = - fem_int(27,sst,rMui)
    NLC( 4,NEQ) = LQi4

    ELM(:,:, 5,NEQ) =   fem_int(37,sst,rMui,Var(:,2)%RUph)
    NLC( 5,NEQ) = LQi1

    ! Residual stress

    ELM(:,:, 6,NEQ) =   amasinv * fem_int(-4,PiRess(:,2))
    NLC( 6,NEQ) = 0

    ! Collisional friction force with electrons

    ELM(:,:, 7,NEQ) =   fem_int(20,lab(:,2,1,1,1),RatPN(:,2))
    NLC( 7,NEQ) = LQe4

    ! Collisional friction force with bulk ions (self collision)

    ELM(:,:, 8,NEQ) =   fem_int(2,lab(:,2,2,1,1))
    NLC( 8,NEQ) = LQi4

    ! Collisional friction force with beam ions

!    ELM(:,:, 9,NEQ) =   fem_int(30,ParatoZeta,laf(:,2,1,1),RatPNb(:,2))
!    NLC( 9,NEQ) = LQb3

    ELM(:,:, 9,NEQ) =   fem_int(20,laf(:,2,1,1),RatPNb(:,2))
    NLC( 9,NEQ) = LQb4

    ! Collisional heat friction force with electrons

    ELM(:,:,10,NEQ) = - fem_int(30,ParatoZeta,lab(:,2,1,1,2),Var(:,2)%n)
    NLC(10,NEQ) = LQe6

    ! Collisional heat friction force with ions

    ELM(:,:,11,NEQ) = - fem_int(30,ParatoZeta,lab(:,2,2,1,2),Var(:,2)%n)
    NLC(11,NEQ) = LQi6

    ! Toroidal E force

    ELM(:,:,12,NEQ) =   achg(2) * aee * amasinv * fem_int(2,Var(:,2)%n)
    NLC(12,NEQ) = LQm3

    ! v x B force

    ELM(:,:,13,NEQ) =   achg(2) * aee * amasinv * fem_int(6,PsiV)
    NLC(13,NEQ) = LQi2

    ! Turbulent particle transport driver (electron driven)

    ELM(:,:,14,NEQ) =   achg(2) * aee * amasinv * fem_int(20,FQLcoef1,RatPN(:,2))
    NLC(14,NEQ) = LQe4

    ELM(:,:,15,NEQ) = - achg(2) * aee * amasinv * fem_int(20,FQLcoef2,Var(:,2)%n)
    NLC(15,NEQ) = LQe3

    ELM(:,:,16,NEQ) =   achg(2) * rKeV * amasinv * fem_int(22,FQLcoef,RatPN(:,2)) &
         &                           * (FSVAHL - 1.d0)
    NLC(16,NEQ) = LQe5

    ELM(:,:,17,NEQ) = - achg(2) * rKeV * amasinv * fem_int(32,FQLcoef,Var(:,1)%T,RatPN(:,2)) &
         &                           * (FSVAHL - 1.d0)
    NLC(17,NEQ) = LQe1

    ELM(:,:,18,NEQ) =   achg(2) * aee * amasinv * fem_int(22,FQLcoef,Var(:,2)%n)
    NLC(18,NEQ) = LQm1

    !  -- Turbulent pinch term

    ELM(:,:,19,NEQ) =   achg(2) * aee * amasinv * fem_int(20,FQLcoef,FVpch)
    NLC(19,NEQ) = LQi1

    ! Loss to divertor

    ELM(:,:,20,NEQ) = - fem_int(2,rNuL) * fact
    NLC(20,NEQ) = LQi4

    ! Increase in momentum due to thermalization of fast ions

!    ELM(:,:,21,NEQ) =   fem_int(20,ParatoZeta,rNuB) * BeamSW
!    NLC(21,NEQ) = LQb3

    ELM(:,:,21,NEQ) =   amb / amas(2) * fem_int(2,rNuB) * BeamSW
    NLC(21,NEQ) = LQb3

    ! Collisional friction force with neutrals

    ELM(:,:,22,NEQ) = - fem_int(2,rNu0i)
    NLC(22,NEQ) = LQi4

    ! Charge exchange force

    ELM(:,:,23,NEQ) = - fem_int(2,rNuiCX)
    NLC(23,NEQ) = LQi4

    ! Loss cone loss

    ELM(:,:,24,NEQ) =   fem_int(-1,SiLCph) * rr
    NLC(24,NEQ) = 0

    ! Ion orbit loss

    ELM(:,:,25,NEQ) = - fem_int(2,rNuOL)
    NLC(25,NEQ) = LQi4

    !  Additional torque input

    ELM(:,:,26,NEQ) =   amasinv * 1.D-20 * fem_int(-1,Tqt) 
    NLC(26,NEQ) = 0

    !  Diffusion of ions (***AF 2008-06-08)

    ELM(:,:,27,NEQ) = - fem_int(27,sst,DMAGi)
    NLC(27,NEQ) = LQi4

!   ! Helical Neoclassical viscosity force
!   ! To deal with the existence of (m=0, n>0) component   miki_m  2010/8/26

!!    ELM(:,:,27,NEQ) =   fem_int(2,rNuiHLphth) * rr
!    IF(UHphSwitch == 0) THEN
!       ELM(:,:,27,NEQ) =   fem_int(2,rNuiHLphth) * rr
!    ELSE
!       ELM(:,:,27,NEQ) =   fem_int(2,rNuiHLphthVR) * rr
!       ELM(1      ,1:4,27,NEQ) =   fem_int_point(2,0,rNuiHLphthVR) * rr
!    ENDIF
!    NLC(27,NEQ) = LQi3

!!    ELM(:,:,28,NEQ) = - fem_int(15,rNuiHLphph)
!    IF(UHphSwitch == 0) THEN
!       ELM(:,:,28,NEQ) = - fem_int(15,rNuiHLphph)
!    else
!       ELM(:,:,28,NEQ) = - fem_int(2,rNuiHLphph)
!    endif
!    NLC(28,NEQ) = LQi4

    NLCMAX(NEQ) = 26

    IF(MDFIXT /= 0) THEN ! Overwrite Pressure gradient force
       ELM(:,:,16,NEQ) = - achg(1) * rKeV * amasinv * fem_int(21,FQLcoef,Var(:,2)%T) &
         &                           * (FSVAHL - 1.d0)
       NLC(16,NEQ) = LQe1

       ELM(:,:,17,NEQ) = 0.D0
       NLC(17,NEQ) = 0
    END IF

    RETURN
  END SUBROUTINE LQi4CC

!***************************************************************
!
!  Ion Energy Transport: Ti
!
!***************************************************************

  SUBROUTINE LQi5CC

    integer(4) :: NEQ = LQi5

    ! Temperature evolution

    IF(MDFIXT == 0) THEN
       ELM(:,:, 0,NEQ) =   1.5D0 * fem_int(1) * invDT
       NLC( 0,NEQ) = LQi5

       ! Advection + Grid velocity

       ELM(:,:, 1,NEQ) = - 2.5d0 * fem_int(5,UsrgV(:,2))
       NLC( 1,NEQ) = LQi5

       ELM(:,:, 2,NEQ) = - 2.5d0 * FSADV * fem_int(6,Var(:,1)%UrV)
       NLC( 2,NEQ) = LQi5

       ! Heat pinch

       ELM(:,:, 3,NEQ) = - fem_int(24,Vhps(:,2),suft)
       NLC( 3,NEQ) = LQi5

       ! Conduction transport

       ELM(:,:, 4,NEQ) = - fem_int(27,sst,Chii1)
       NLC( 4,NEQ) = LQi5

       ELM(:,:, 5,NEQ) =   fem_int(37,sst,Chii2,Var(:,2)%T)
       NLC( 5,NEQ) = LQi1

       ! Redundant heat convection term

       ELM(:,:, 6,NEQ) =   fem_int(5,UsrgV(:,2))
       NLC( 6,NEQ) = LQi5

       ! Viscous heating

       ELM(:,:, 7,NEQ) =   1.D-20 / rKeV * fem_int(-2,Var(:,2)%Uthhat,BnablaPi(:,2))
       NLC( 7,NEQ) = 0

       ! Collisional transfer with electrons (Energy equilibration)

       ELM(:,:, 8,NEQ) = - 1.5d0 * fem_int(20,rNuTei,RatPN(:,1))
       NLC( 8,NEQ) = LQi5

       ELM(:,:, 9,NEQ) =   1.5d0 * fem_int( 2,rNuTei)
       NLC( 9,NEQ) = LQe5

       ! Loss to diverter

       ELM(:,:,10,NEQ) = -          fem_int( 2,rNuL)
       NLC(10,NEQ) = LQi5

       ELM(:,:,11,NEQ) =   PNiDIV * fem_int(-2,rNuL,Var(:,2)%T)
       NLC(11,NEQ) = 0

!!$       ELM(:,:,10,NEQ) = - 1.D0   / achg(2) * fem_int(20,rNuL,RatPN(:,2))
!!$       NLC(10,NEQ) = LQi5
!!$
!!$       ELM(:,:,11,NEQ) =   PNeDIV / achg(2) * fem_int(-2,rNuL,Var(:,2)%T)
!!$       NLC(11,NEQ) = 0

       ELM(:,:,12,NEQ) = - 1.5D0            * fem_int( 2,rNuLTi)
       NLC(12,NEQ) = LQi5

       ELM(:,:,13,NEQ) =   1.5D0 * PTiDIV   * fem_int(-2,rNuLTi,Var(:,2)%n)
       NLC(13,NEQ) = 0
!!$       ELM(:,:,13,NEQ) =   1.5D0 * PTiDIV      * fem_int( 2,rNuLTi)
!!$       NLC(13,NEQ) = LQi1

       ! Ionization heating of n01, n02 and n03

       ELM(:,:,14,NEQ) = 1.5D0 * FSION * 1.D20 * fem_int(30,SiVizA,Var(:,2)%n,PT01V)
       NLC(14,NEQ) = LQn1

       ELM(:,:,15,NEQ) = 1.5D0 * FSION * 1.D20 * fem_int(30,SiVizA,Var(:,2)%n,PT02V)
       NLC(15,NEQ) = LQn2

       ELM(:,:,16,NEQ) = 1.5D0 * FSION * 1.D20 * fem_int(30,SiVizA,Var(:,2)%n,PT03V) * BeamSW
       NLC(16,NEQ) = LQn3

!!$       ELM(:,:,14,NEQ) = 1.5D0 / achg(2) * FSION * 1.D20 * fem_int(30,SiVizA,Var(:,1)%n,PT01V)
!!$       NLC(14,NEQ) = LQn1
!!$
!!$       ELM(:,:,15,NEQ) = 1.5D0 / achg(2) * FSION * 1.D20 * fem_int(30,SiVizA,Var(:,1)%n,PT02V)
!!$       NLC(15,NEQ) = LQn2
!!$
!!$       ELM(:,:,16,NEQ) = 1.5D0 / achg(2) * FSION * 1.D20 * fem_int(30,SiVizA,Var(:,1)%n,PT03V) * BeamSW
!!$       NLC(16,NEQ) = LQn3

       ! Charge exchange loss due to slow neutrals
       !   (Thermal neutrals are assumed to have same temperature with ions.)

       ELM(:,:,17,NEQ) = - 1.5D0 * fem_int( 2,rNuiCXT)
       NLC(17,NEQ) = LQi5

       ELM(:,:,18,NEQ) =   1.5D0 * fem_int(-20,rNuiCXT(:),Var(:,2)%n,PT01V(:))
       NLC(18,NEQ) = 0

       ! Collisional NBI heating (Perp + Tan)

       ELM(:,:,19,NEQ) = Eb * fem_int(-2,SNB,PNBcol_i)
       NLC(19,NEQ) = 0

       ! Heating due to beam momentum deposition

       ELM(:,:,20,NEQ) = achgb * amqp * 1.d-3 * fem_int(-20,MNB,BUsparVbbt(:,2),PNbVinv)
       NLC(20,NEQ) = LQb3

       ! Simplified Alpha heating

       ELM(:,:,21,NEQ) = 1.D-20 / rKeV * fem_int(-1,PALFi)
       NLC(21,NEQ) = 0

       ! Direct heating (RF)

       ELM(:,:,22,NEQ) = 1.D-20 / rKeV * fem_int(-1,PRFi)
       NLC(22,NEQ) = 0

       !  Diffusion of ions (***AF 2008-06-08)

       ELM(:,:,23,NEQ) = - fem_int(27,sst,DMAGi)
       NLC(23,NEQ) = LQi5

       NLCMAX(NEQ) = 22
    ELSE

       !  Fixed temperature profile

       ELM(:,:,0,NEQ) = fem_int(1) * invDT
       NLC(0,NEQ) = LQi5

       NLCMAX(NEQ) = 0
    END IF

    RETURN
  END SUBROUTINE LQi5CC

!***************************************************************
!
!   Ion Heat Parallel Flow
!
!***************************************************************

  SUBROUTINE LQi6CC

    integer(4) :: NEQ = LQi6
    real(8) :: amasinv

    amasinv = 1.d0 / (amas(2) * amp)

    ELM(:,:, 0,NEQ) =   2.5d0 * fem_int(2+L6,Var(:,2)%n) * invDT
    NLC( 0,NEQ) = LQi6

    ! Neoclassical viscosity force

    ELM(:,:, 1,NEQ) = - amasinv * 1.D-20 * fem_int(2+L6,xmu(:,2,2,1))
    NLC( 1,NEQ) = LQi3

    ELM(:,:, 2,NEQ) = - amasinv * 1.D-20 * fem_int(2+L6,xmu(:,2,2,2))
    NLC( 2,NEQ) = LQi6

    ! Diamagnetic forces

    ELM(:,:, 3,NEQ) = - rKilo / achg(2) * amasinv * fem_int(22+L6,fipolsdtPNsV(:,2),xmu(:,2,2,1))
    NLC( 3,NEQ) = LQi5

    ELM(:,:, 4,NEQ) = - amasinv * 1.D-20 * fem_int(22+L6,fipolsdt,xmu(:,2,2,1))
    NLC( 4,NEQ) = LQm1

    ELM(:,:, 5,NEQ) = - rKilo / achg(2) * amasinv * fem_int(22+L6,fipolsdtPNsV(:,2),xmu(:,2,2,2))
    NLC( 5,NEQ) = LQi5

    ELM(:,:, 6,NEQ) =   rKilo / achg(2) * amasinv * fem_int(32+L6,fipolsdtPNsV(:,2),Var(:,2)%T,xmu(:,2,2,2))
    NLC( 6,NEQ) = LQi1

    ! Collisional friction force with electrons

    ELM(:,:, 7,NEQ) = - fem_int(20+L6,lab(:,2,1,2,1),Var(:,2)%n)
    NLC( 7,NEQ) = LQe3

    ! Collisional friction force with ions

    ELM(:,:, 8,NEQ) = - fem_int(20+L6,lab(:,2,2,2,1),Var(:,2)%n)
    NLC( 8,NEQ) = LQi3

    ! Collisional friction force with beam ions

    ELM(:,:, 9,NEQ) = - fem_int(20+L6,laf(:,2,2,1),  RatPNb(:,2))
    NLC( 9,NEQ) = LQb3

    ! Collisional heat friction force with electrons

    ELM(:,:,10,NEQ) =   fem_int(20+L6,lab(:,2,1,2,2),Var(:,2)%n)
    NLC(10,NEQ) = LQe6

    ! Collisional heat friction force with ions (self collision)

    ELM(:,:,11,NEQ) =   fem_int(20+L6,lab(:,2,2,2,2),Var(:,2)%n)
    NLC(11,NEQ) = LQi6

    ! Loss to divertor

    ELM(:,:,12,NEQ) = - fem_int(20+L6,rNuL,Var(:,2)%n) * fact
    NLC(12,NEQ) = LQi6

    ! Collisional friction force with neutrals

    ELM(:,:,13,NEQ) = - fem_int(20+L6,rNu0i,Var(:,2)%n)
    NLC(13,NEQ) = LQi6

    ! Charge exchange force

    ELM(:,:,14,NEQ) = - fem_int(20+L6,rNuiCX,Var(:,2)%n)
    NLC(14,NEQ) = LQi6

    ! Loss cone loss

    ELM(:,:,15,NEQ) =   fem_int(-(2+L6),SiLCth,Var(:,2)%n)
    NLC(15,NEQ) = 0

    ! Ion orbit loss

    ELM(:,:,16,NEQ) = - fem_int(20+L6,rNuOL,Var(:,2)%n)
    NLC(16,NEQ) = LQi6

    NLCMAX(NEQ) = 13

    IF(MDFIXT /= 0) THEN ! Overwrite Pressure gradient force
       ELM(:,:, 3,NEQ) = - rKilo / achg(2) * amasinv &
            &                  *(  fem_int(31+L6,fipolsdtPNsV(:,2),xmu(:,2,2,1),Var(:,2)%T) &
            &                    + fem_int(32+L6,fipolsdtPNsV(:,2),xmu(:,2,2,1),Var(:,2)%T))
       NLC( 3,NEQ) = LQi1

       ELM(:,:, 5,NEQ) =   amasinv * 1.D-20 * fem_int(-(2+L6),xmu(:,2,2,2),BVsdiag(:,2,2))
       NLC( 5,NEQ) = 0

       ELM(:,:, 6,NEQ) =   0.d0
       NLC( 6,NEQ) = 0
    END IF

    RETURN
  END SUBROUTINE LQi6CC

!***************************************************************
!
!   Ion rotation frequency: Ni <UiPhi/R>
!
!***************************************************************

  SUBROUTINE LQi7CC

    integer(4) :: NEQ = LQi7

    !  Ion rotation frequency

    ELM(:,:,1,NEQ) = - fem_int(1)
    NLC(1,NEQ) = LQi7

    !  Ion rotation

    ELM(:,:,2,NEQ) =   fem_int(2,rrtinv)
    NLC(2,NEQ) = LQi4

    ELM(:,:,3,NEQ) = - fem_int(20,Fqhatsq,rrtinv)
    NLC(3,NEQ) = LQi4

    !  Ion parallel velocity

    ELM(:,:,4,NEQ) =   fem_int(20,FqhatsqI,Var(:,2)%n)
    NLC(4,NEQ) = LQi3

    NLCMAX(NEQ) = 4
    RETURN
  END SUBROUTINE LQi7CC

!***************************************************************
!
!   Beam Ion Density
!
!***************************************************************

  SUBROUTINE LQb1CC

    integer(4) :: NEQ = LQb1

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQb1

    if(MDBEAM == 0) then

    !  Grid velocity

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
    RETURN
  END SUBROUTINE LQb1CC

!***************************************************************
!
!   Beam Ion Radial Flow
!
!***************************************************************
  
  SUBROUTINE LQb2CC

    integer(4) :: NEQ = LQb2
    real(8) :: amasinv

    amasinv = 1.d0 / (amb * amp)

    if( MDBEAM == 0 ) then

    ! Dummy term to avoid null when not solving LQb2

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQb2

    else

    ! Pressure gradient force

    ELM(:,:,1,NEQ) = - rKilo * aee * amasinv &
         &                 *(  fem_int(21+L2,bri,PTbV) &
         &                   + fem_int(22+L2,bri,PTbV))
    NLC(1,NEQ) = LQb1

    ! Radial E force

    ELM(:,:,2,NEQ) = - achgb * aee * amasinv * fem_int(22+L2,bri,PNbV)
    NLC(2,NEQ) = LQm1

    ! v x B force

    ELM(:,:,3,NEQ) =   achgb * aee * amasinv * fem_int(20+L2,sdt,fipol)
    NLC(3,NEQ) = LQb3

    ELM(:,:,4,NEQ) = - achgb * aee * amasinv * fem_int(20+L2,sdt,bbt)
    NLC(4,NEQ) = LQb4

    end if

    NLCMAX(NEQ) = 4
    RETURN
  END SUBROUTINE LQb2CC

!***************************************************************
!
!   Beam Ion Parallel flow
!
!***************************************************************

  SUBROUTINE LQb3CC

    integer(4) :: NEQ = LQb3
    real(8) :: amasinv

    amasinv = 1.d0 / (amb * amp)

!!$    if( iSUPG3 > 0 ) L3 = 0

    ELM(:,:,0,NEQ) = fem_int(1+L3) * invDT
    NLC(0,NEQ) = LQb3

    ! Neoclassical viscosity

    ELM(:,:,1,NEQ) = - amasinv * 1.d-20 * fem_int(20+L3,PNbVinv,xmuf(:,1))
    NLC(1,NEQ) = LQb3

    ! Diamagnetic forces

!    ELM(:,:,2,NEQ) = - rKilo / achgb * amasinv &
!            &                  *(  fem_int(31+L3,fipolsdtPNbV,xmuf(:,1),PTbV) &
!            &                    + fem_int(32+L3,fipolsdtPNbV,xmuf(:,1),PTbV))
    ELM(:,:,2,NEQ) = - rKilo / achgb * amasinv &
            &                  * fem_int(32+L3,fipolsdtPNbV,xmuf(:,1),PTbV)
    NLC(2,NEQ) = LQb1

    ELM(:,:,3,NEQ) = - rKilo / achgb * amasinv * 1.D-20 &
            &                  * fem_int(-(22+L3),fipolsdt,xmuf(:,1),PTbV)
    NLC(3,NEQ) = 0

    ELM(:,:,4,NEQ) = - amasinv * 1.D-20 * fem_int(22+L3,fipolsdt,xmuf(:,1))
    NLC(4,NEQ) = LQm1

    ! Collisional friction force with beam ions (self collision)

    ELM(:,:,5,NEQ) =   (amas(1) / amb) * fem_int(20+L3,lff(:,1,1),RatPNb(:,1))
    NLC(5,NEQ) = LQb3

    ! Collisional friction force with electrons

    ELM(:,:,6,NEQ) =   (amas(1) / amb) * fem_int(20+L3,lfb(:,1,1,1),Var(:,1)%n)
    NLC(6,NEQ) = LQe3

    ! Collisional friction force with ions

    ELM(:,:,7,NEQ) =   (amas(1) / amb) * fem_int(20+L3,lfb(:,2,1,1),Var(:,1)%n)
    NLC(7,NEQ) = LQi3

    ! Electric field forces

    ELM(:,:,8,NEQ) = - achgb * aee * amasinv * fem_int(30+L3,PNbV,fipol,aatq)
    NLC(8,NEQ) = LQm2

    ELM(:,:,9,NEQ) =   achgb * aee * amasinv * fem_int(30+L3,PNbV,fipol,aat)
    NLC(9,NEQ) = LQm3

    ! Loss to divertor

    ELM(:,:,10,NEQ) = - fem_int(2+L3,rNuLB) * fact! * BeamSW
    NLC(10,NEQ) = LQb3

    ! Momentum loss due to thermalization

    ELM(:,:,11,NEQ) = - fem_int(2+L3,rNuB)! * BeamSW
    NLC(11,NEQ) = LQb3

    ! Collisional friction force with neutrals

    ELM(:,:,12,NEQ) = - fem_int(2+L3,rNu0b)! * BeamSW
    NLC(12,NEQ) = LQb3

    ! Charge exchange force

    ELM(:,:,13,NEQ) = - fem_int(2+L3,rNubCX)! * BeamSW
    NLC(13,NEQ) = LQb3

    ! NBI momentum source

    ELM(:,:,14,NEQ) =  amasinv * fem_int(-(1+L3),BSmb)
    NLC(14,NEQ) = 0

    ! Momentum loss due to collisional ripple trapping

    ELM(:,:,15,NEQ) = - fem_int(20+L3,rNubrp2,rip_rat) * RpplSW
    NLC(15,NEQ) = LQb3

    ! Momentum diffusion arising from beam ion convective flux due to ripple

    ELM(:,:,16,NEQ) = - fem_int(37,sst,BUbparV,Dbrpft)
    NLC(16,NEQ) = LQb1

    NLCMAX(NEQ) = 16
    RETURN
  END SUBROUTINE LQb3CC

!***************************************************************
!
!   Beam Ion Toroidal Flow
!
!***************************************************************

  SUBROUTINE LQb4CC

    integer(4) :: NEQ = LQb4
    real(8) :: amasinv

    amasinv = 1.d0 / (amb * amp)

    ! Ubphi'(0) : 0

    ELM(:,:, 0,NEQ) =   fem_int(1) * invDT
    NLC( 0,NEQ) = LQb4

    ! Advection + Grid velocity

    ELM(:,:, 1,NEQ) =   fem_int(5,UgV)
    NLC( 1,NEQ) = LQb4

    if( MDBEAM /= 0 ) then

    ELM(:,:, 2,NEQ) = - FSADVB * fem_int(3,RUbphV)
    NLC( 2,NEQ) = LQb2

    end if

    ! Collisional friction force with electrons

    ELM(:,:, 3,NEQ) =   (amas(1) / amb) * fem_int(2,lfb(:,1,1,1))
    NLC( 3,NEQ) = LQe4

    ! Collisional friction force with bulk ions

    ELM(:,:, 4,NEQ) =   (amas(1) / amb) * fem_int(20,lfb(:,2,1,1),RatPN(:,1))
    NLC( 4,NEQ) = LQi4

    ! Collisional friction force with beam ions (self collision)

    ELM(:,:, 5,NEQ) =   (amas(1) / amb) * fem_int(20,lff(:,1,1),RatPNb(:,1))
    NLC( 5,NEQ) = LQb4

    ! Toroidal E force

    ELM(:,:, 6,NEQ) =   achgb * aee * amasinv * fem_int(2,PNbV)
    NLC( 6,NEQ) = LQm3

    if( MDBEAM /= 0 ) then

    ! v x B force

    ELM(:,:, 7,NEQ) =   achgb * aee * amasinv * fem_int(6,PsiV)
    NLC( 7,NEQ) = LQb2

    end if

    ! Loss to divertor

    ELM(:,:, 8,NEQ) = - fem_int(2,rNuLB) * fact
    NLC( 8,NEQ) = LQb4

    ! Momentum loss due to thermalization

    ELM(:,:, 9,NEQ) = - fem_int(2,rNuB)
    NLC( 9,NEQ) = LQb4

    ! Collisional friction force with neutrals

    ELM(:,:,10,NEQ) = - fem_int(2,rNu0b)
    NLC(10,NEQ) = LQb4

    ! Charge exchange force

    ELM(:,:,11,NEQ) = - fem_int(2,rNubCX)
    NLC(11,NEQ) = LQb4

    ! NBI momentum source : (I / <B^2>) m_b dot{n}_b <B v_//0>

    ELM(:,:,12,NEQ) =   amasinv * fem_int(-2,ParatoZeta,BSmb)
    NLC(12,NEQ) = 0

!    ! Turbulent particle transport driver (electron driven)
!
!    ELM(:,:,13,NEQ) =   achgb * aee * amasinv * fem_int(20,FQLcoef1,RatPNbinv(:,1))
!    NLC(13,NEQ) = LQe4
!
!    ELM(:,:,14,NEQ) = - achgb * aee * amasinv * fem_int(20,FQLcoef2,PNbV)
!    NLC(14,NEQ) = LQe3
!
!    ELM(:,:,15,NEQ) =   achgb * rKeV * amasinv * fem_int(22,FQLcoef,RatPNbinv(:,1)) &
!         &                           * (FSVAHL - 1.d0)
!    NLC(15,NEQ) = LQe5
!
!    ELM(:,:,16,NEQ) = - achgb * rKeV * amasinv * fem_int(32,FQLcoef,Var(:,1)%T,RatPNbinv(:,1)) &
!         &                           * (FSVAHL - 1.d0)
!    NLC(16,NEQ) = LQe1
!
!    ELM(:,:,17,NEQ) =   achgb * aee * amasinv * fem_int(22,FQLcoef,PNbV)
!    NLC(17,NEQ) = LQm1
!
!    !  -- Turbulent pinch term
!
!    ELM(:,:,18,NEQ) =   achgb * aee * amasinv * fem_int(20,FQLcoef,FVpch)
!    NLC(18,NEQ) = LQb1

!    ! Loss cone loss
!
!    ELM(:,:,13,NEQ) =   fem_int(-(1),SiLCph) * rr
!    NLC(13,NEQ) = 0
!
!    ! Ion orbit loss
!
!    ELM(:,:,14,NEQ) = - fem_int(2,rNuOL)
!    NLC(14,NEQ) = LQi4

    NLCMAX(NEQ) = 12
    RETURN
  END SUBROUTINE LQb4CC

!***************************************************************
!
!  Beam Ion rotation frequency: Nb <UbPhi/R>
!
!***************************************************************

  SUBROUTINE LQb7CC

    integer(4) :: NEQ = LQb7

    !  Beam Ion rotation frequency

    ELM(:,:,1,NEQ) = - fem_int(1)
    NLC(1,NEQ) = LQb7

    !  Beam Ion rotation

    ELM(:,:,2,NEQ) =   fem_int(2,rrtinv)
    NLC(2,NEQ) = LQb4

    ELM(:,:,3,NEQ) = - fem_int(20,Fqhatsq,rrtinv)
    NLC(3,NEQ) = LQb4

    !  Beam Ion parallel velocity

    ELM(:,:,4,NEQ) =   fem_int(2,FqhatsqI)
    NLC(4,NEQ) = LQb3

    NLCMAX(NEQ) = 4
    RETURN
  END SUBROUTINE LQb7CC

!***************************************************************
!
!   Slow Neutral Transport: n01
!
!***************************************************************

  SUBROUTINE LQn1CC

    integer(4) :: NEQ = LQn1

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQn1

    !  Grid velocity

    ELM(:,:,1,NEQ) = - fem_int(6,UgV)
    NLC(1,NEQ) = LQn1

    !  Diffusion of neutrals

    ELM(:,:,2,NEQ) = - fem_int(27,sst,D01)
    NLC(2,NEQ) = LQn1

    ! Ionization

    ELM(:,:,3,NEQ) = - FSION * 1.D20 * fem_int(20,SiVizA,Var(:,2)%n)
    NLC(3,NEQ) = LQn1

!!$    ELM(:,:,3,NEQ) = - 1.D0 / achg(2) * FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n)
!!$    NLC(3,NEQ) = LQn1

    ! Generation of thermal neutrals by charge exchange

    ELM(:,:,4,NEQ) = - FSCX * 1.D20 * fem_int(20,SiVcxA,Var(:,2)%n)
    NLC(4,NEQ) = LQn1

    ! Recycling from divertor

    ELM(:,:,5,NEQ) =   rGamm0          * fem_int( 2,rNuL)
    NLC(5,NEQ) = LQi1

    ELM(:,:,6,NEQ) = - rGamm0 * PNiDIV * fem_int(-1,rNuL)
    NLC(6,NEQ) = 0

!!$    ELM(:,:,5,NEQ) =   rGamm0          / achg(2) * fem_int( 2,rNuL)
!!$    NLC(5,NEQ) = LQe1
!!$
!!$    ELM(:,:,6,NEQ) = - rGamm0 * PNeDIV / achg(2) * fem_int(-1,rNuL)
!!$    NLC(6,NEQ) = 0

!    ! Gas puff (alternative way to impose B.C.;
!                Deactivate call boundary(lqn1) if this term is activated)
!
!    ELM(:,:,7,NEQ) = fem_int(-3,suft,rGASPFA)
!    NLC(7,NEQ) = 0

    NLCMAX(NEQ) = 6
    RETURN
  END SUBROUTINE LQn1CC

!***************************************************************
!
!   Thermal Neutral Transport: n02
!
!***************************************************************

  SUBROUTINE LQn2CC

    integer(4) :: NEQ = LQn2

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQn2

    !  Grid velocity

    ELM(:,:,1,NEQ) = - fem_int(6,UgV)
    NLC(1,NEQ) = LQn2

    !  Diffusion of neutrals

    ELM(:,:,2,NEQ) = - fem_int(27,sst,D02)
    NLC(2,NEQ) = LQn2

    ! Ionization

    ELM(:,:,3,NEQ) = - FSION * 1.D20 * fem_int(20,SiVizA,Var(:,2)%n)
    NLC(3,NEQ) = LQn2

!!$    ELM(:,:,3,NEQ) = - 1.D0 / achg(2) * FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n)
!!$    NLC(3,NEQ) = LQn2

    ! Generation of thermal neutrals by charge exchange

    ELM(:,:,4,NEQ) = FSCX * 1.D20 * fem_int(20,SiVcxA,Var(:,2)%n)
    NLC(4,NEQ) = LQn1

!    ! NBI particle source (Charge exchange), replaced by LQn3CC
!
!    ELM(:,:,5,NEQ) = RatCX * fem_int(-1,SNBi)
!    NLC(5,NEQ) = 0

    NLCMAX(NEQ) = 4
    RETURN
  END SUBROUTINE LQn2CC

!***************************************************************
!
!   Halo Neutral Transport: n03
!
!***************************************************************

  SUBROUTINE LQn3CC

    integer(4) :: NEQ = LQn3

    ELM(:,:,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQn3

    !  Grid velocity

    ELM(:,:,1,NEQ) = - fem_int(6,UgV)
    NLC(1,NEQ) = LQn3

    !  Diffusion of neutrals

    ELM(:,:,2,NEQ) = - fem_int(27,sst,D03)
    NLC(2,NEQ) = LQn3

    ! Ionization

    ELM(:,:,3,NEQ) = - FSION * 1.D20 * fem_int(20,SiVizA,Var(:,2)%n)
    NLC(3,NEQ) = LQn3

!!$    ELM(:,:,3,NEQ) = - 1.D0 / achg(2) * FSION * 1.D20 * fem_int(20,SiVizA,Var(:,1)%n)
!!$    NLC(3,NEQ) = LQn3

    ! NBI particle source (Charge exchange)

    ELM(:,:,4,NEQ) = RatCX * fem_int(-1,SNBi)
    NLC(4,NEQ) = 0

    NLCMAX(NEQ) = 4
    RETURN
  END SUBROUTINE LQn3CC

!***************************************************************
!
!   Ripple Trapped Beam Ion Density (SUPG)
!
!***************************************************************

  SUBROUTINE LQr1CC

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
    RETURN
  END SUBROUTINE LQr1CC

!***************************************************************
!
!   Boundary condition
!
!      NLCR(NC,NQ,0) : inner boundary (magnetic axis)
!      NLCR(NC,NQ,1) : outer boundary (virtual wall)
!
!***************************************************************

  SUBROUTINE BOUNDARY(NR,LQ,ID,VAL)

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
    IF(ID == 0) THEN
       ! Initialize ALC, BLC and CLC at NR
       ALC(NR,:,LQ) = 0.D0
       BLC(NR,:,LQ) = 0.D0
       CLC(NR,:,LQ) = 0.D0
       PLC(NR,:,LQ) = 0.D0

       ! The boundary condition matrices IDX and IMAX are created if ICALA=0
       !   to store the information of the Dirichlet boundary conditions.
       ! ** IDX stores both the equation and term numbers, NQ and NC, where the variable LQ appears.
       ! ** IMAX stores the maximum number of IDX for each LQ.
       IF(ICALA == 0) THEN
          I = 0
          DO NQ = 1, NQMAX
             DO NC = 1, NLCMAX(NQ)
                IF(NLCR(NC,NQ,NB) == LQ) THEN
                   I = I + 1
                   IDX(I,LQ,NB)%IDXNC = NC
                   IDX(I,LQ,NB)%IDXNQ = NQ
                END IF
             END DO
          END DO
          IMAX(LQ) = I
       ELSE
       ! Read indices from the boundary condition matrices created at initial (ICALA=1)
          IF(PRESENT(VAL)) THEN
             IF(NR == 0) THEN ! axis
                DO I = 1, IMAX(LQ)
                   NC = IDX(I,LQ,NB)%IDXNC
                   NQ = IDX(I,LQ,NB)%IDXNQ
                   PLC(NR  ,NC,NQ) = PLC(NR,NC,NQ)   + BLC(NR,NC,NQ)   * VAL
                   PLC(NR+1,NC,NQ) = PLC(NR+1,NC,NQ) + CLC(NR+1,NC,NQ) * VAL
                END DO
             ELSE ! wall
                DO I = 1, IMAX(LQ)
                   NC = IDX(I,LQ,NB)%IDXNC
                   NQ = IDX(I,LQ,NB)%IDXNQ
                   PLC(NR  ,NC,NQ) = PLC(NR,NC,NQ)   + BLC(NR,NC,NQ)   * VAL
                   PLC(NR-1,NC,NQ) = PLC(NR-1,NC,NQ) + ALC(NR-1,NC,NQ) * VAL
                END DO
             END IF
          END IF

          IF(NR == 0) THEN ! axis
             DO I = 1, IMAX(LQ)
                NC = IDX(I,LQ,NB)%IDXNC
                NQ = IDX(I,LQ,NB)%IDXNQ
                BLC(NR  ,NC,NQ) = 0.D0
                CLC(NR+1,NC,NQ) = 0.D0
             END DO
          ELSE ! wall
             DO I = 1, IMAX(LQ)
                NC = IDX(I,LQ,NB)%IDXNC
                NQ = IDX(I,LQ,NB)%IDXNQ
                BLC(NR  ,NC,NQ) = 0.D0
                ALC(NR-1,NC,NQ) = 0.D0
             END DO
          END IF
       END IF

       ! Diagonal term on handled variable
       BLC (NR,1,LQ) = 1.D0
       NLCR(1,LQ,NB) = LQ

       IF(PRESENT(VAL)) PLC(NR,1,LQ) = PLC(NR,1,LQ) - VAL
    ELSE
    ! === Neumann boundary condition ===
       NLCMAX(LQ) = NLCMAX(LQ) + 1
       PLC(NR,NLCMAX(LQ),LQ) = VAL * DTf(LQ)
       NLCR(NLCMAX(LQ),LQ,NB) = 0
    END IF

  END SUBROUTINE BOUNDARY

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
