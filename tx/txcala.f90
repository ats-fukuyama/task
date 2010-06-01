!     $Id$
module tx_coefficients
  use tx_commons
  use tx_core_module
  implicit none
  private
  real(8), dimension(:,:,:,:), allocatable :: ELM
  real(8), dimension(:), allocatable ::  &
       & rNuIN0, rNuCXN0, rNubeBE, rNubiBI, rNuTeiEI,&
       & rNuCXN1, rMueNe, rMuiNi, & !dPNeV, dPNiV, &
       & RUbthV, UethVR, UithVR, EthVR, RUerV, RUirV, UerVR, UirVR, &
       & FWpheBB, FWpheBB2, FWthpheB, FWpheB, FWphe, dAphV, &
!!ion       & FWphiBB, FWphiBB2, FWthphiB, FWphiB, &
       & BphBNi, BthBNi, Dbrpft, &
       & rNuei1EI, rNuei2BthEI, rNuei3EI, &
       & rNube1BE, rNube2BthBE, rNube3BE, &
       & Vbparaph, RVbparath, &
       & Chie1, Chie2, Chii1, Chii2, &
       & FVpchph, rGASPFA, RNeUer, RNiUir, &
       & rNue2NCN, rNui2NCN, rat_ei
!!rp_conv       &, rNubLL
  real(8), dimension(:), allocatable :: UNITY
  real(8) :: DTt, DTf(1:NQM), invDT, BeamSW, RpplSW, ThntSW, FSVAHLT, FSVAHLE, &
       &     fact = 1.d0 ! <= SOL loss accelerator
  integer(4), save :: ICALA = 0
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

    allocate(ELM(1:NEMAX,1:4,0:NCM,1:NQMAX))

    ELM(1:NEMAX,1:4,0:NCM,1:NQMAX) = 0.D0

    ALC(0:NCM,1:NQMAX,0:NRMAX) = 0.D0
    BLC(0:NCM,1:NQMAX,0:NRMAX) = 0.D0
    CLC(0:NCM,1:NQMAX,0:NRMAX) = 0.D0
    NLCR(0:NCM,1:NQMAX,0:NRMAX) = 0
    NLC(0:NCM,1:NQMAX) = 0
    PLC(1:NCM,1:NQMAX,0:NRMAX) = 0.D0
    NLCMAX(1:NQMAX) = 0

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
       RpplSW = 1.D0
    ELSE
       RpplSW = 0.D0
    END IF

    ! If ThntSW = 0, there is no particle source from PN02V.
    ThntSW = 1.D0

    ! Annihilator of inherent particle pinch induced by FWth terms
    IF(MDVAHL == 1) THEN
       FSVAHLT = FSVAHL
       FSVAHLE = 0.D0
    ELSE IF(MDVAHL >= 2) THEN
       IF(MDVAHL == 2) THEN
          FSVAHLT = FSVAHL
       ELSE
          FSVAHLT = 0.D0
       END IF
       FSVAHLE = 1.D0
    ELSE ! MDVAHL == 0
       FSVAHLT = 1.D0
       FSVAHLE = 0.D0
    END IF

    !     Coefficients

    N = NRMAX

    allocate(rNuIN0(0:N), rNuCXN0(0:N), rNubeBE(0:N), rNubiBI(0:N), rNuTeiEI(0:N), &
       &     rNuCXN1(0:N), rMueNe(0:N), rMuiNi(0:N), &!dPNeV(0:N), dPNiV(0:N), &
       &     RUbthV(0:N), UethVR(0:N), UithVR(0:N), &
       &     EthVR(0:N), RUerV(0:N), RUirV(0:N), UerVR(0:N), UirVR(0:N), &
       &     FWpheBB(0:N), FWpheBB2(0:N), FWthpheB(0:N), FWpheB(0:N), FWphe(0:N), &
       &     dAphV(0:N), BphBNi(0:N), BthBNi(0:N), Dbrpft(0:N), &
       &     rNuei1EI(0:N), rNuei2BthEI(0:N), rNuei3EI(0:N), &
       &     rNube1BE(0:N), rNube2BthBE(0:N), rNube3BE(0:N), &
       &     Vbparaph(0:N), RVbparath(0:N), &
       &     Chie1(0:N), Chie2(0:N), Chii1(0:N), Chii2(0:N), &
       &     FVpchph(0:N), rGASPFA(0:N), RNeUer(0:N), RNiUir(0:N))!, &
!       &     rNue2NCN(0:N), rNui2NCN(0:N), rat_ei(0:N))
!!ion    allocate(FWphiBB(0:N), FWphiBB2(0:N), FWthphiB(0:N), FWphiB(0:N))
    allocate(UNITY(0:N))
    UNITY(0:N) = 1.D0

    CALL LQCOEF

    !     Maxwell

    CALL LQm1CC
    CALL LQm2CC
    CALL LQm3CC
    CALL LQm4CC
    CALL LQm5CC

    !     Electron

    CALL LQe1CC
    CALL LQe2CC
    CALL LQe3CC
    CALL LQe4CC
    CALL LQe5CC

    !     Ion

    CALL LQi1CC
    CALL LQi2CC
    CALL LQi3CC
    CALL LQi4CC
    CALL LQi5CC

    !     Beam

    CALL LQb1CC
    CALL LQb3CC
    CALL LQb4CC

    !     Neutral

    CALL LQn1CC
    CALL LQn2CC
    CALL LQn3CC

    !     Ripple trapped beam

    CALL LQr1CC

    !     Elemental equations -> Nodal equations

    DO NE = 1, NEMAX
       NR = NE - 1
       DO NQ = 1, NQMAX
          NC = 0
             BLC(NC,NQ,NR)   = BLC(NC,NQ,NR)   + ELM(NE,1,NC,NQ) * DTt
             ALC(NC,NQ,NR)   = ALC(NC,NQ,NR)   + ELM(NE,2,NC,NQ) * DTt
             CLC(NC,NQ,NR+1) = CLC(NC,NQ,NR+1) + ELM(NE,3,NC,NQ) * DTt
             BLC(NC,NQ,NR+1) = BLC(NC,NQ,NR+1) + ELM(NE,4,NC,NQ) * DTt

          DO NC = 1, NLCMAX(NQ)
             ! iHvsLC : Heaviside function
             !          iHvsLC = 1 when NLC = 0 ; otherwise iHvsLC = 0.
             iHvsLC = 1 - ceiling(real(NLC(NC,NQ))/(NQMAX+1))
             BLC(NC,NQ,NR)   = BLC(NC,NQ,NR)   + ELM(NE,1,NC,NQ) * DTf(NQ)
             ALC(NC,NQ,NR)   = ALC(NC,NQ,NR)   + ELM(NE,2,NC,NQ) * DTf(NQ)
             PLC(NC,NQ,NR)   = PLC(NC,NQ,NR)   + sum(ELM(NE,1:2,NC,NQ)) * DTf(NQ) * iHvsLC
             CLC(NC,NQ,NR+1) = CLC(NC,NQ,NR+1) + ELM(NE,3,NC,NQ) * DTf(NQ)
             BLC(NC,NQ,NR+1) = BLC(NC,NQ,NR+1) + ELM(NE,4,NC,NQ) * DTf(NQ)
             PLC(NC,NQ,NR+1) = PLC(NC,NQ,NR+1) + sum(ELM(NE,3:4,NC,NQ)) * DTf(NQ) * iHvsLC
          END DO
       END DO
    END DO

    DO NR = 0, NRMAX
       NLCR(0:NCM,1:NQMAX,NR) = NLC(0:NCM,1:NQMAX)
    END DO

    !     Dirichlet condition

    CALL BOUNDARY(NRMAX,LQm1,0)
    CALL BOUNDARY(0    ,LQm2,0)
    CALL BOUNDARY(0    ,LQe2,0)
!    CALL BOUNDARY(NRMAX,LQe2,0)
    CALL BOUNDARY(0    ,LQe3,0)
    CALL BOUNDARY(NRMAX,LQe3,0)
    CALL BOUNDARY(NRMAX,LQe4,0)
    CALL BOUNDARY(0    ,LQi2,0)
!    CALL BOUNDARY(NRMAX,LQi2,0)
    CALL BOUNDARY(0    ,LQi3,0)
    CALL BOUNDARY(NRMAX,LQi3,0)
    CALL BOUNDARY(NRMAX,LQi4,0)
    CALL BOUNDARY(NRMAX,LQn2,0)
    CALL BOUNDARY(NRMAX,LQn3,0)
    ! When ripple effect is on (FSRP /= 0), ripple diffusion term will be activated.
    ! Then we must impose two boundary conditions at each equation.
    CALL BOUNDARY(0    ,LQb3,0)
    CALL BOUNDARY(NRMAX,LQb3,0)
    CALL BOUNDARY(NRMAX,LQb4,0)
!    IF(FSRP /= 0.D0) CALL BOUNDARY(0,LQr1,0)

    ! Neumann condition of the continuity equation at the boundary is naturally
    ! imposed through the diffusion term in the pressure equation. However, this
    ! condition is very weak because it comes from nondiagonal term. Then by setting
    ! the density constant in the last element, we can force the density gradient
    ! to be nought at the boundary explicitly.

    CALL BOUNDARY(NRMAX,LQe1,0,PNeV(NRMAX-1))
    CALL BOUNDARY(NRMAX,LQi1,0,PNiV(NRMAX-1))

    !     Integral term stemming from integration by parts in the diffusion term

    CALL BOUNDARY(NRMAX,LQm2,1, 2.D0*PSI(NRMAX)*BB/rMU0)
    CALL BOUNDARY(NRMAX,LQm3,1,-2.D0*R(NRMAX)*Bthb/rMUb2)
!    CALL BOUNDARY(NRMAX,LQn1,1, 2.D0*R(NRMAX)*rGASPF)

    deallocate(ELM)
    deallocate(rNuIN0, rNuCXN0, rNubeBE, rNubiBI, rNuTeiEI,&
       &       rNuCXN1, rMueNe, rMuiNi, & !dPNeV, dPNiV, &
       &       RUbthV, UethVR, UithVR, EthVR, RUerV, RUirV, UerVR, UirVR, &
       &       FWpheBB, FWpheBB2, FWthpheB, FWpheB, FWphe, dAphV, &
       &       BphBNi, BthBNi, Dbrpft, &
       &       rNuei1EI, rNuei2BthEI, rNuei3EI, &
       &       rNube1BE, rNube2BthBE, rNube3BE, &
       &       Vbparaph, RVbparath, &
       &       Chie1, Chie2, Chii1, Chii2, &
       &       FVpchph, rGASPFA, RNeUer, RNiUir)!, &
!       &       rNue2NCN, rNui2NCN, rat_ei)
!!ion    deallocate(FWphiBB, FWphiBB2, FWthphiB, FWphiB)
    deallocate(UNITY)

    IF(ICALA ==0) ICALA = 1

    RETURN
  END SUBROUTINE TXCALA

!***************************************************************
!
!   Coefficients for Equations
!
!**************************************************************

  SUBROUTINE LQCOEF

    use tx_interface, only : dfdx!, derivs
    INTEGER(4) :: NR
    REAL(8) :: BBL

!!$    rNuIN0(0:NRMAX)   = rNuION(0:NRMAX) * PNeV(0:NRMAX) &
!!$         &            / (PN01V(0:NRMAX) + PN02V(0:NRMAX) + PN03V(0:NRMAX))
!!$    rNuCXN0(0:NRMAX)  = rNuiCX(0:NRMAX) * PNiV(0:NRMAX) &
!!$         &            / (PN01V(0:NRMAX) + PN02V(0:NRMAX) + PN03V(0:NRMAX))
    rNuIN0(0:NRMAX)   = FSION * SiVizA(0:NRMAX) * PNeV(0:NRMAX) * 1.D20
    rNuCXN0(0:NRMAX)  = FSCX  * SiVcxA(0:NRMAX) * PNiV(0:NRMAX) * 1.D20
    rNuCXN1(0:NRMAX)  = rNuiCXT(0:NRMAX) * PNiV(0:NRMAX) * PT01V(0:NRMAX)
    rNuTeiEI(0:NRMAX) = rNuTei(0:NRMAX) * PNeV(0:NRMAX) / PNiV(0:NRMAX)
!!oldNBI    rNubeBE(0:NRMAX)  = rNube(0:NRMAX)  * PNbV(0:NRMAX) / PNeV(0:NRMAX)
    rNubiBI(0:NRMAX)  = rNubi(0:NRMAX)  * PNbV(0:NRMAX) / PNiV(0:NRMAX)
!    rMueNe(0:NRMAX)   = rMue(0:NRMAX)   / PNeV(0:NRMAX)
!    rMuiNi(0:NRMAX)   = rMui(0:NRMAX)   / PNiV(0:NRMAX)

    RUerV(0:NRMAX)    = R(0:NRMAX)      * UerV(0:NRMAX)
    RUirV(0:NRMAX)    = R(0:NRMAX)      * UirV(0:NRMAX)
    UerVR(1:NRMAX)    = UerV(1:NRMAX) / R(1:NRMAX)
    UerVR(0)          = 0.D0 ! Any value is OK. (Never affect the result.)
    UirVR(1:NRMAX)    = UirV(1:NRMAX) / R(1:NRMAX)
    UirVR(0)          = 0.D0 ! Any value is OK. (Never affect the result.)
    UethVR(1:NRMAX)   = UethV(1:NRMAX) / R(1:NRMAX)
    UethVR(0)         = 0.D0 ! Any value is OK. (Never affect the result.)
    UithVR(1:NRMAX)   = UithV(1:NRMAX) / R(1:NRMAX)
    UithVR(0)         = 0.D0 ! Any value is OK. (Never affect the result.)
    EthVR(1:NRMAX)    = EthV(1:NRMAX) / R(1:NRMAX)
    EthVR(0)          = 0.D0 ! Any value is OK. (Never affect the result.)
!    CALL DERIVS(PSI,X,LQe1,NQMAX,NRMAX,dPNeV)
!    CALL DERIVS(PSI,X,LQi1,NQMAX,NRMAX,dPNiV)
    dAphV(0:NRMAX)    = dfdx(PSI,AphV,NRMAX,0)
    FWphe(0:NRMAX)    =- 2.D0 * dAphV(0:NRMAX) * FWe(0:NRMAX)
    FWpheBB(0:NRMAX)  =- 2.D0 * dAphV(0:NRMAX) * FWthphe(0:NRMAX)
    FWpheBB2(0:NRMAX) = FWe(0:NRMAX) * BthV(0:NRMAX)**2
!!ion    FWphi(0:NRMAX)    =- 2.D0 * dAphV(0:NRMAX) * FWi(0:NRMAX)
!!ion    FWphiBB(0:NRMAX)  =- 2.D0 * dAphV(0:NRMAX) * FWthphi(0:NRMAX)
!!ion    FWphiBB2(0:NRMAX) = FWthi(0:NRMAX) * BthV(0:NRMAX)**2

    DO NR = 0, NRMAX
       BBL = SQRT(BphV(NR)**2 + BthV(NR)**2)
       BphBNi(NR) = BphV(NR) / (BBL * PNiV(NR))
       IF(NR /= 0) BthBNi(NR) = BthV(NR) / (BBL * R(NR) * PNiV(NR))
       if(MDLMOM == 0) then
          Vbparaph(NR)  = Vbpara(NR)
          RVbparath(NR) = 0.d0
       else
          Vbparaph(NR)  = Vbpara(NR) * (BphV(NR) / BBL)
          RVbparath(NR) = Vbpara(NR) * (BthV(NR) / BBL) * R(NR)
       end if

       FWthpheB(NR) = FWthphe(NR) * BBL
       FWpheB  (NR) = FWphe  (NR) * BBL
!!ion       FWthphiB(NR) = FWthphi(NR) * BBL
!!ion       FWphiB  (NR) = FWphi  (NR) * BBL
    END DO
    BthBNi(0) = 0.D0 ! Any value is OK. (Never affect the result.)

!!rp_conv    rNubLL(0:NRMAX) = rNubL(0:NRMAX) * rip_rat(0:NRMAX)

    RUbthV(0:NRMAX) = R(0:NRMAX) * UbthV(0:NRMAX)
    Dbrpft(0:NRMAX) = Dbrp(0:NRMAX) * ft(0:NRMAX)

    rNuei1EI   (0:NRMAX)  = rNuei1   (0:NRMAX)  * PNeV(0:NRMAX) / PNiV(0:NRMAX)
    rNuei2BthEI(0:NRMAX)  = rNuei2Bth(0:NRMAX)  * PNeV(0:NRMAX) / PNiV(0:NRMAX)
    rNuei3EI   (0:NRMAX)  = rNuei3   (0:NRMAX)  * PNeV(0:NRMAX) / PNiV(0:NRMAX)
    rNube1BE   (0:NRMAX)  = rNube1   (0:NRMAX)  * PNbV(0:NRMAX) / PNeV(0:NRMAX)
    rNube2BthBE(0:NRMAX)  = rNube2Bth(0:NRMAX)  * PNbV(0:NRMAX) / PNeV(0:NRMAX)
    rNube3BE   (0:NRMAX)  = rNube3   (0:NRMAX)  * PNbV(0:NRMAX) / PNeV(0:NRMAX)

    Chie1(0:NRMAX) = Chie(0:NRMAX) + ChiNCTe(0:NRMAX) + ChiNCpe(0:NRMAX)
    Chie2(0:NRMAX) = Chie(0:NRMAX) + ChiNCTe(0:NRMAX)
    Chii1(0:NRMAX) = Chii(0:NRMAX) + ChiNCTi(0:NRMAX) + ChiNCpi(0:NRMAX)
    Chii2(0:NRMAX) = Chii(0:NRMAX) + ChiNCTi(0:NRMAX)

    FVpchph (0:NRMAX) = FVpch(0:NRMAX)  / BphV(0:NRMAX)

    rGASPFA(0:NRMAX-1) = 0.D0
    rGASPFA(NRMAX) = rGASPF
    RNeUer(0:NRMAX) = R(0:NRMAX) * PNeV(0:NRMAX) * UerV(0:NRMAX)
    RNiUir(0:NRMAX) = R(0:NRMAX) * PNiV(0:NRMAX) * UirV(0:NRMAX)

!    rNue2NCN(0:NRMAX) = rNue2NC(0:NRMAX) / PNeV(0:NRMAX)
!    rNui2NCN(0:NRMAX) = rNui2NC(0:NRMAX) / PNiV(0:NRMAX)
!    rat_ei(0:NRMAX) = PNeV(0:NRMAX) / PNiV(0:NRMAX)

  END SUBROUTINE LQCOEF

!***************************************************************
!
!   Poisson Equation: phi
!
!**************************************************************

  SUBROUTINE LQm1CC

    integer(4) :: NEQ = LQm1

    ! phi'(0) : 0

    ELM(1:NEMAX,1:4,1,NEQ) =   4.D0 * sqeps0 / (AEE * 1.D20) * fem_int(18,UNITY)
    NLC(1,NEQ) = LQm1

    ELM(1:NEMAX,1:4,2,NEQ) =        fem_int(1) / sqeps0
    NLC(2,NEQ) = LQe1

    ELM(1:NEMAX,1:4,3,NEQ) = - PZ * fem_int(1) / sqeps0
    NLC(3,NEQ) = LQi1

    ELM(1:NEMAX,1:4,4,NEQ) = - PZ * fem_int(1) / sqeps0 * BeamSW
    NLC(4,NEQ) = LQb1

    ELM(1:NEMAX,1:4,5,NEQ) = - PZ * fem_int(2,rip_rat) / sqeps0 * BeamSW * RpplSW
    NLC(5,NEQ) = LQr1

    ! phi(b) : 0

    NLCMAX(NEQ) = 5
    RETURN
  END SUBROUTINE LQm1CC

!***************************************************************
!
!   Ampere's Law: Atheta'
!
!***************************************************************

  SUBROUTINE LQm2CC

    integer(4) :: NEQ = LQm2

    ! (r*Atheta)'(0) : 0

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * EPS0 * invDT
    NLC(0,NEQ) = LQm2

    ! rot Bphi

    ELM(1:NEMAX,1:4,1,NEQ) = - 4.D0 * fem_int(18,UNITY) - 4.D0 * fem_int(4)
    NLC(1,NEQ) = LQm5

    ! Electron current

    ELM(1:NEMAX,1:4,2,NEQ) = -      AEE * 1.D20 * fem_int(1)
    NLC(2,NEQ) = LQe3

    ! Ion current

    ELM(1:NEMAX,1:4,3,NEQ) =   PZ * AEE * 1.D20 * fem_int(1)
    NLC(3,NEQ) = LQi3

    ! Beam ion current

    ELM(1:NEMAX,1:4,4,NEQ) =   PZ * AEE * 1.D20 * fem_int(1) * BeamSW
    NLC(4,NEQ) = LQb3

    ! (1/r)*(r*Atheta)'(NRMAX) : BB

    NLCMAX(NEQ) = 4
    RETURN
  END SUBROUTINE LQm2CC

!***************************************************************
!
!   Ampere's Law: Aphi'
!
!***************************************************************

  SUBROUTINE LQm3CC

    integer(4) :: NEQ = LQm3

    ! Aphi'(0) : 0

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * EPS0 * rMUb1 * invDT
    NLC(0,NEQ) = LQm3

    ! rot Btheta

    ELM(1:NEMAX,1:4,1,NEQ) = - 4.D0 * fem_int(18,UNITY)
    NLC(1,NEQ) = LQm4

    ! Electron current

    ELM(1:NEMAX,1:4,2,NEQ) = - rMUb1      * AEE * 1.D20 * fem_int(1)
    NLC(2,NEQ) = LQe4

    ! Ion current

    ELM(1:NEMAX,1:4,3,NEQ) =   rMUb1 * PZ * AEE * 1.D20 * fem_int(1)
    NLC(3,NEQ) = LQi4

    ! Beam ion current

    ELM(1:NEMAX,1:4,4,NEQ) =   rMUb1 * PZ * AEE * 1.D20 * fem_int(1) * BeamSW
    NLC(4,NEQ) = LQb4

!   ! Virtual current for helical system

    ELM(1:NEMAX,1:4,5,NEQ) =   rMUb1 * fem_int(-1,AJV)
    NLC(5,NEQ) = 0

    ! Aphi'(NRMAX) : -Bthb

    NLCMAX(NEQ) = 5
    RETURN
  END SUBROUTINE LQm3CC

!**************************************************************
!
!   Faraday's Law : Aphi
!
!***************************************************************

  SUBROUTINE LQm4CC

    integer(4) :: NEQ = LQm4

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQm4

    ! Aphi'

    ELM(1:NEMAX,1:4,1,NEQ) = fem_int(1) / rMUb2
    NLC(1,NEQ) = LQm3

    NLCMAX(NEQ) = 1
    RETURN
  END SUBROUTINE LQm4CC

!***************************************************************
!
!   Faraday's Law : Atheta
!
!***************************************************************

  SUBROUTINE LQm5CC

    integer(4) :: NEQ = LQm5

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQm5

    ! r * Atheta'

    ELM(1:NEMAX,1:4,1,NEQ) = fem_int(1) / rMU0
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

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQe1

    ! Divergence

    ELM(1:NEMAX,1:4,1,NEQ) = - 2.D0 * fem_int(4)
    NLC(1,NEQ) = LQe2

    ! Ionization of n01, n02 and n03

    ELM(1:NEMAX,1:4,2,NEQ) =   fem_int(2,rNuIN0)
    NLC(2,NEQ) = LQn1

    ELM(1:NEMAX,1:4,3,NEQ) =   fem_int(2,rNuIN0) * ThntSW
    NLC(3,NEQ) = LQn2

    ELM(1:NEMAX,1:4,4,NEQ) =   fem_int(2,rNuIN0) * BeamSW
    NLC(4,NEQ) = LQn3

    ! Loss to divertor

    ELM(1:NEMAX,1:4,5,NEQ) = - fem_int(2,rNuL)
    NLC(5,NEQ) = LQe1

    ELM(1:NEMAX,1:4,6,NEQ) =  PNeDIV * fem_int(-1,rNuL)
    NLC(6,NEQ) = 0

    ! Generated by NBI (Ionization)

    ELM(1:NEMAX,1:4,7,NEQ) =   (1.D0 - RatCX) * fem_int(-1,SNBe)
    NLC(7,NEQ) = 0

    !  Diffusion of electrons (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,8,NEQ) = - 4.D0 * fem_int(18,DMAGe)
    NLC(8,NEQ) = LQe1

    NLCMAX(NEQ) = 8
    RETURN
  END SUBROUTINE LQe1CC

!***************************************************************
!
!   Electron Radial Flow (SUPG)
!
!***************************************************************

  SUBROUTINE LQe2CC

    integer(4) :: NEQ = LQe2

    ! Ns*Usr(0) : fixed

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT &
         &                 + fem_int(8) * fem_int(0) * invDT
    NLC(0,NEQ) = LQe2

    ! Advection
    
    ELM(1:NEMAX,1:4,1,NEQ) = - 2.D0 * fem_int( 3,RUerV) + fem_int(2,UerVR) &
         &                 +(- 2.D0 * fem_int(10,RUerV) + fem_int(9,UerVR)) * fem_int(0)
    ELM(1      ,1:4,1,NEQ) = - 2.D0 * fem_int_point( 3,0,RUerV) + fem_int_point(2,0,UerVR) &
         &                 +(- 2.D0 * fem_int_point(10,0,RUerV) + fem_int_point(9,0,UerVR)) &
         &                 * fem_int_point(0,1)
    NLC(1,NEQ) = LQe2

    ! Centrifugal force

    ELM(1:NEMAX,1:4,2,NEQ) = fem_int(2,UethVR) &
         &                 + fem_int(9,UethVR) * fem_int(0)
    ELM(1      ,1:4,2,NEQ) = fem_int_point(2,0,UethVR) &
         &                 + fem_int_point(9,0,UethVR) * fem_int_point(0,1)
    NLC(2,NEQ) = LQe3

    ! Pressure gradient force

    ELM(1:NEMAX,1:4,3,NEQ) = - 2.D0 * rKeV / AME * fem_int(17,UNITY) &
         &                   - 2.D0 * rKeV / AME * fem_int(18,UNITY) * fem_int(0)
    NLC(3,NEQ) = LQe5

    ! Radial E force

    ELM(1:NEMAX,1:4,4,NEQ) =   2.D0 * (AEE / AME) * fem_int(17,PNeV) &
         &                 +   2.D0 * (AEE / AME) * fem_int(18,PNeV) * fem_int(0)
    NLC(4,NEQ) = LQm1

    ! v x B force

    ELM(1:NEMAX,1:4,5,NEQ) = -        (AEE / AME) * fem_int( 2,BphV) &
         &                   -        (AEE / AME) * fem_int( 9,BphV) * fem_int(0)
    NLC(5,NEQ) = LQe3

    ELM(1:NEMAX,1:4,6,NEQ) = - 2.D0 * (AEE / AME) * fem_int(16,AphV) &
         &                   - 2.D0 * (AEE / AME) * fem_int(20,UNITY,AphV) * fem_int(0)
    NLC(6,NEQ) = LQe4

    NLCMAX(NEQ) = 6

    IF(MDFIXT /= 0) THEN
       ELM(1:NEMAX,1:4,3,NEQ) = - 2.D0 * rKeV / AME &
            &                 *(  fem_int(17,PTeV) + fem_int(18,PTeV) * fem_int(0) &
            &                   + fem_int(16,PTeV) + fem_int(34,PSI,PTeV) * fem_int(0))
       NLC(3,NEQ) = LQe1
    END IF

    RETURN
  END SUBROUTINE LQe2CC

!***************************************************************
!
!   Electron Poloidal Flow
!
!***************************************************************

  SUBROUTINE LQe3CC

    integer(4) :: NEQ = LQe3

    ! Ns*UsTheta(0) : 0

    ELM(1:NEMAX,1:4, 0,NEQ) = fem_int(1) * invDT
    NLC( 0,NEQ) = LQe3

    ! Advection

    ELM(1:NEMAX,1:4, 1,NEQ) = - 2.D0 * fem_int(3,RUerV)
    NLC( 1,NEQ) = LQe3

    ! Viscosity force

!    ELM(1:NEMAX,1:4, 2,NEQ) = - 4.D0 * fem_int(18,rMue) &
!         &                    - 4.D0 * fem_int(39,rMueNe,dPNeV) &
!         &                    - 4.D0 * fem_int( 3,rMue)
    ELM(1:NEMAX,1:4, 2,NEQ) = - 4.D0 * fem_int(18,rMue) &
         &                    - 4.D0 * fem_int( 3,rMue)
    NLC( 2,NEQ) = LQe3

    ELM(1:NEMAX,1:4, 3,NEQ) =   4.D0 * fem_int(41,rMue,RUethV)
    NLC( 3,NEQ) = LQe1

    ! Poloidal E force

    ELM(1:NEMAX,1:4, 4,NEQ) =   (AEE / AME) * fem_int(2,PNeV)
    NLC( 4,NEQ) = LQm2

    ! v x B force

    ELM(1:NEMAX,1:4, 5,NEQ) =   (AEE / AME) * fem_int(2,BphV)
    NLC( 5,NEQ) = LQe2

    ! Neoclassical viscosity force

    ELM(1:NEMAX,1:4, 6,NEQ) = - fem_int(2,rNueNC)
    NLC( 6,NEQ) = LQe3

    ! Collisional friction force with ions

    ELM(1:NEMAX,1:4, 7,NEQ) = - fem_int(2,rNuei1)
    NLC( 7,NEQ) = LQe3

    ELM(1:NEMAX,1:4, 8,NEQ) =   fem_int(2,rNuei1EI)
    NLC( 8,NEQ) = LQi3

    ELM(1:NEMAX,1:4, 9,NEQ) =   2.D0 * fem_int(37,rNuei2Bth,AphV)
    NLC( 9,NEQ) = LQe4

    ELM(1:NEMAX,1:4,10,NEQ) = - 2.D0 * fem_int(37,rNuei2BthEI,AphV)
    NLC(10,NEQ) = LQi4

    ! Collisional friction with beam ions

    ELM(1:NEMAX,1:4,11,NEQ) = - (AMB / AME) * fem_int(2,rNube1BE)
    NLC(11,NEQ) = LQe3

    ELM(1:NEMAX,1:4,12,NEQ) =   (AMB / AME) * fem_int(2,rNube1) * BeamSW
    NLC(12,NEQ) = LQb3

    ELM(1:NEMAX,1:4,13,NEQ) =   2.D0 * (AMB / AME) * fem_int(37,rNube2BthBE,AphV)
    NLC(13,NEQ) = LQe4

    ELM(1:NEMAX,1:4,14,NEQ) = - 2.D0 * (AMB / AME) * fem_int(37,rNube2Bth,AphV) * BeamSW
    NLC(14,NEQ) = LQb4

    ! Turbulent particle transport driver

!!$    ELM(1:NEMAX,1:4,15,NEQ) =   2.D0 * AEE / AME * fem_int(38,De,BphV)
!!$    NLC(15,NEQ) = LQe1
!!$
!!$    ELM(1:NEMAX,1:4,16,NEQ) =   2.D0 * AEE / AME * fem_int(38,De,BphV/PTeV) * FSVAHLT
!!$    NLC(16,NEQ) = LQe5
!!$
!!$    ELM(1:NEMAX,1:4,17,NEQ) = - 2.D0 * AEE / AME * fem_int(38,De,BphV) * FSVAHLT
!!$    NLC(17,NEQ) = LQe1
!!$
    ELM(1:NEMAX,1:4,15,NEQ) = - 1.D0 / AME * fem_int(2,FWthe)
    NLC(15,NEQ) = LQe3

!    ELM(1:NEMAX,1:4,16,NEQ) =   1.D0 / AME * fem_int(28,FWthe,rat_ei)
!    NLC(16,NEQ) = LQi3

    ELM(1:NEMAX,1:4,17,NEQ) =   1.D0 / AME * fem_int(15,FWpheBB)
    NLC(17,NEQ) = LQe4

!    ELM(1:NEMAX,1:4,18,NEQ) = - 1.D0 / AME * fem_int(44,FWpheBB,rat_ei)
!    NLC(18,NEQ) = LQi4

    ! Wave interaction force (electron driven)

    ELM(1:NEMAX,1:4,19,NEQ) =   1.D0 / AME * fem_int(44,FWthpheB,WPM)
    NLC(19,NEQ) = LQe1

    ! Contribution of off-diagonal term due to heat flux

    ELM(1:NEMAX,1:4,20,NEQ) = - 2.D0 * (rKeV / (AME * AEE)) * fem_int(17,FWthphe) &
         &                           * (1.D0 - FSVAHLT)
    NLC(20,NEQ) = LQe5

    ELM(1:NEMAX,1:4,21,NEQ) = + 2.D0 * (rKeV / (AME * AEE)) * fem_int(38,FWthphe,PTeV) &
         &                           * (1.D0 - FSVAHLT)
    NLC(21,NEQ) = LQe1

    ELM(1:NEMAX,1:4,22,NEQ) = + 2.D0 * (1.D0 /  AME       ) * fem_int(37,FWthphe,PhiV) * FSVAHLE
    NLC(22,NEQ) = LQe1

    ! Ad hoc turbulent pinch term

    ELM(1:NEMAX,1:4,23,NEQ) = - 1.D0 / AME * fem_int(15,FVpch)
    NLC(23,NEQ) = LQe1

!!ion       ! Turbulent particle transport driver (ion driven)
!!ion
!!ion       ELM(1:NEMAX,1:4,15,NEQ) =   1.D0 / AME * fem_int(2,FWthi)
!!ion       NLC(15,NEQ) = LQi3
!!ion
!!ion       ELM(1:NEMAX,1:4,16,NEQ) = - 1.D0 / AME * fem_int(2,FWthi)
!!ion       NLC(16,NEQ) = LQe3
!!ion
!!ion       ELM(1:NEMAX,1:4,17,NEQ) = - 1.D0 / AME * fem_int(15,FWphiBB)
!!ion       NLC(17,NEQ) = LQi4
!!ion
!!ion       ELM(1:NEMAX,1:4,18,NEQ) =   1.D0 / AME * fem_int(15,FWphiBB)
!!ion       NLC(18,NEQ) = LQe4
!!ion
!!ion       ELM(1:NEMAX,1:4,19,NEQ) =   2.D0 * (rKeV / (AME * AEE)) * fem_int(37,FWthphi,PTiV) &
!!ion                                        * (1.D0 - FSVAHLT) &
!!ion            &                    - 2.D0 * (1.D0 /  AME       ) * fem_int(37,FWthphi,PhiV) * FSVAHLE
!!ion       NLC(19,NEQ) = LQi1
!!ion
!!ion       ELM(1:NEMAX,1:4,20,NEQ) =   1.D0 / AME * fem_int(44,FWthphiB,WPM)
!!ion       NLC(20,NEQ) = LQi1

    ! Loss to divertor

    ELM(1:NEMAX,1:4,24,NEQ) = - 2.D0 * fem_int(2,rNuL) * fact
    NLC(24,NEQ) = LQe3

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,25,NEQ) = - fem_int(2,rNu0e)
    NLC(25,NEQ) = LQe3

    ! Helical neoclassical viscosity force (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,26,NEQ) = - fem_int(2,rNueHLthth)
    NLC(26,NEQ) = LQe3

!---- 09/11/26 AF
    ELM(1:NEMAX,1:4,27,NEQ) =   fem_int(15,rNueHLthph)
    NLC(27,NEQ) = LQe4

    !  Diffusion of electrons (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,28,NEQ) = - 4.D0 * fem_int(18,DMAGe)
    NLC(28,NEQ) = LQe3

    !  Poloidal torque due to neoclassical heat flux

    ELM(1:NEMAX,1:4,29,NEQ) =   FSNCPL / 1.D20 * fem_int(-2,R,FQeth)
    NLC(29,NEQ) = 0

!    ELM(1:NEMAX,1:4,29,NEQ) =   FSNCPL * 2.D0 * rKeV / AEE * fem_int(17,rNue2NC)
!    NLC(29,NEQ) = LQe5
!
!    ELM(1:NEMAX,1:4,30,NEQ) = - FSNCPL * 2.D0 * rKeV / AEE * fem_int(37,rNue2NCN,PNeV)
!    NLC(30,NEQ) = LQe5

    ! Ns*UsTheta(NRMAX) : 0

    NLCMAX(NEQ) = 29

    IF(MDFIXT /= 0) THEN
       ELM(1:NEMAX,1:4,20,NEQ) = - 2.D0 * (rKeV / (AME * AEE)) * fem_int(37,FWthphe,PTeV) &
            &                           * (1.D0 - FSVAHLT)
       NLC(20,NEQ) = LQe1

       ELM(1:NEMAX,1:4,21,NEQ) = 0.D0
       NLC(21,NEQ) = 0

!       ELM(1:NEMAX,1:4,29,NEQ) =   FSNCPL * 2.D0 * rKeV / AEE * fem_int(37,rNue2NC,PTeV)
!       NLC(29,NEQ) = LQe1
!
!       ELM(1:NEMAX,1:4,30,NEQ) = 0.D0
!       NLC(30,NEQ) = 0
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

    ! Uephi(0)' : 0

    ELM(1:NEMAX,1:4, 0,NEQ) = fem_int(1) * invDT
    NLC( 0,NEQ) = LQe4

    ! Advection

    ELM(1:NEMAX,1:4, 1,NEQ) = - 2.D0 * fem_int(3,RUerV)
    NLC( 1,NEQ) = LQe4
    
    ! Viscosity force
    
!    ELM(1:NEMAX,1:4, 2,NEQ) = - 4.D0 *(fem_int(18,rMue) + fem_int(39,rMueNe,dPNeV))
    ELM(1:NEMAX,1:4, 2,NEQ) = - 4.D0 * fem_int(18,rMue)
    NLC( 2,NEQ) = LQe4

    ELM(1:NEMAX,1:4, 3,NEQ) =   4.D0 * fem_int(41,rMue,UephV)
    NLC( 3,NEQ) = LQe1

    ! Toroidal E force

    ELM(1:NEMAX,1:4, 4,NEQ) =   (AEE / AME) * fem_int(2,PNeV)
    NLC( 4,NEQ) = LQm3

    ! v x B force

    ELM(1:NEMAX,1:4, 5,NEQ) =   2.D0 * (AEE / AME) * fem_int(6,AphV)
    NLC( 5,NEQ) = LQe2

    ! Collisional friction with bulk ions

    ELM(1:NEMAX,1:4, 6,NEQ) = - fem_int(2,rNuei3)
    NLC( 6,NEQ) = LQe4

    ELM(1:NEMAX,1:4, 7,NEQ) =   fem_int(2,rNuei3EI)
    NLC( 7,NEQ) = LQi4

    ELM(1:NEMAX,1:4, 8,NEQ) =   2.D0 * fem_int(29,rNuei2Bth,AphV)
    NLC( 8,NEQ) = LQe3

    ELM(1:NEMAX,1:4, 9,NEQ) = - 2.D0 * fem_int(29,rNuei2BthEI,AphV)
    NLC( 9,NEQ) = LQi3

    ! Collisional friction with beam ions

    ELM(1:NEMAX,1:4,10,NEQ) = - (AMB / AME) * fem_int(2,rNube3BE)
    NLC(10,NEQ) = LQe4

    ELM(1:NEMAX,1:4,11,NEQ) =   (AMB / AME) * fem_int(2,rNube3) * BeamSW
    NLC(11,NEQ) = LQb4

    ELM(1:NEMAX,1:4,12,NEQ) =   2.D0 * (AMB / AME) * fem_int(29,rNube2BthBE,AphV)
    NLC(12,NEQ) = LQe3

    ELM(1:NEMAX,1:4,13,NEQ) = - 2.D0 * (AMB / AME) * fem_int(29,rNube2Bth,AphV) * BeamSW
    NLC(13,NEQ) = LQb3

    ! Turbulent particle transport driver

!!$    ELM(1:NEMAX,1:4,14,NEQ) = - 2.D0 * AEE / AME * fem_int(38,De,BthV)
!!$    NLC(14,NEQ) = LQe1
!!$
!!$    ELM(1:NEMAX,1:4,15,NEQ) = - 2.D0 * AEE / AME * fem_int(38,De,BthV/PTeV) * FSVAHLT
!!$    NLC(15,NEQ) = LQe5
!!$
!!$    ELM(1:NEMAX,1:4,16,NEQ) =   2.D0 * AEE / AME * fem_int(38,De,BthV) * FSVAHLT
!!$    NLC(16,NEQ) = LQe1
!!$
    ELM(1:NEMAX,1:4,14,NEQ) =   1.D0 / AME * fem_int(2,FWpheBB)
    NLC(14,NEQ) = LQe3

!    ELM(1:NEMAX,1:4,15,NEQ) = - 1.D0 / AME * fem_int(28,FWpheBB,rat_ei)
!    NLC(15,NEQ) = LQi3

    ELM(1:NEMAX,1:4,16,NEQ) = - 1.D0 / AME * fem_int(2,FWpheBB2)
    NLC(16,NEQ) = LQe4

!    ELM(1:NEMAX,1:4,17,NEQ) =   1.D0 / AME * fem_int(28,FWpheBB2,rat_ei)
!    NLC(17,NEQ) = LQi4

    ! Wave interaction force (electron driven)

    ELM(1:NEMAX,1:4,18,NEQ) = - 1.D0 / AME * fem_int(44,FWpheB,WPM)
    NLC(18,NEQ) = LQe1

    ! Contribution of off-diagonal term due to heat flux

    ELM(1:NEMAX,1:4,19,NEQ) =   2.D0 * (rKeV / (AME * AEE)) * fem_int(17,FWphe) &
         &                           * (1.D0 - FSVAHLT)
    NLC(19,NEQ) = LQe5

    ELM(1:NEMAX,1:4,20,NEQ) = - 2.D0 * (rKeV / (AME * AEE)) * fem_int(38,FWphe,PTeV) &
         &                           * (1.D0 - FSVAHLT)
    NLC(20,NEQ) = LQe1

    ELM(1:NEMAX,1:4,21,NEQ) = - 2.D0 * (1.D0 /  AME       ) * fem_int(37,FWphe,PhiV) * FSVAHLE
    NLC(21,NEQ) = LQe1

    ! Ad hoc turbulent pinch term

    ELM(1:NEMAX,1:4,22,NEQ) = - 2.D0 / AME * fem_int(36,AphV,FVpchph)
    NLC(22,NEQ) = LQe1

!!ion       ! Turbulent particle transport driver (ion driven)
!!ion
!!ion       ELM(1:NEMAX,1:4,14,NEQ) = - 1.D0 / AME * fem_int(2,FWphiBB)
!!ion       NLC(14,NEQ) = LQi3
!!ion
!!ion       ELM(1:NEMAX,1:4,15,NEQ) =   1.D0 / AME * fem_int(2,FWphiBB)
!!ion       NLC(15,NEQ) = LQe3
!!ion
!!ion       ELM(1:NEMAX,1:4,16,NEQ) =   1.D0 / AME * fem_int(2,FWphiBB2)
!!ion       NLC(16,NEQ) = LQi4
!!ion
!!ion       ELM(1:NEMAX,1:4,17,NEQ) = - 1.D0 / AME * fem_int(2,FWphiBB2)
!!ion       NLC(17,NEQ) = LQe4
!!ion
!!ion       ELM(1:NEMAX,1:4,18,NEQ) = - 2.D0 * (rKeV / (AME * AEE)) * fem_int(37,FWphiBB,PTiV) &
!!ion            &                           * (1.D0 - FSVAHLT) &
!!ion            &                    + 2.D0 * (1.D0 /  AME       ) * fem_int(37,FWphiBB,PhiV) * FSVAHLE
!!ion       NLC(18,NEQ) = LQi1
!!ion
!!ion       ELM(1:NEMAX,1:4,19,NEQ) = - 1.D0 / AME * fem_int(44,FWphiB,WPM)
!!ion       NLC(19,NEQ) = LQi1

    ! Loss to divertor

    ELM(1:NEMAX,1:4,23,NEQ) = - 2.D0 * fem_int(2,rNuL) * fact
    NLC(23,NEQ) = LQe4

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,24,NEQ) = - fem_int(2,rNu0e)
    NLC(24,NEQ) = LQe4

    ! Helical neoclassical viscosity force (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,25,NEQ) =   fem_int(2,rNueHLphth)
    NLC(25,NEQ) = LQe3

    ELM(1:NEMAX,1:4,26,NEQ) = - fem_int(15,rNueHLphph)
    NLC(26,NEQ) = LQe4

    !  Diffusion of electrons (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,27,NEQ) = - 4.D0 * fem_int(18,DMAGe)
    NLC(27,NEQ) = LQe4

    NLCMAX(NEQ) = 27

    IF(MDFIXT /= 0) THEN
       ELM(1:NEMAX,1:4,19,NEQ) =   2.D0 * (rKeV / (AME * AEE)) * fem_int(37,FWphe,PTeV) &
            &                           * (1.D0 - FSVAHLT)
       NLC(19,NEQ) = LQe1

       ELM(1:NEMAX,1:4,20,NEQ) = 0.D0
       NLC(20,NEQ) = 0
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
       ELM(1:NEMAX,1:4, 0,NEQ) =   1.5D0 * fem_int(1) * invDT
       NLC( 0,NEQ) = LQe5

       ! Advection

       ELM(1:NEMAX,1:4, 1,NEQ) = - 5.D0 * fem_int(3,RUerV)
       NLC( 1,NEQ) = LQe5

       ! Conduction transport

       ELM(1:NEMAX,1:4, 2,NEQ) = - 4.D0 * fem_int(18,Chie1)
       NLC( 2,NEQ) = LQe5

       ELM(1:NEMAX,1:4, 3,NEQ) =   4.D0 * fem_int(41,Chie2,PTeV)
       NLC( 3,NEQ) = LQe1

       ! Redundant heat advection

       ELM(1:NEMAX,1:4, 4,NEQ) = 2.D0 * fem_int(5,RUerV)
       NLC( 4,NEQ) = LQe5

       ! Joule heating

       ELM(1:NEMAX,1:4, 5,NEQ) = - AEE / rKeV * fem_int(2,EthVR)
       ELM(1      ,1:4, 5,NEQ) = - AEE / rKeV * fem_int_point(2,0,EthVR)
       NLC( 5,NEQ) = LQe3

       ELM(1:NEMAX,1:4, 6,NEQ) = - AEE / rKeV * fem_int(2,EphV)
       NLC( 6,NEQ) = LQe4

       ! Collisional transfer with ions (Energy equilibration)

       ELM(1:NEMAX,1:4, 7,NEQ) = - 1.5d0 * fem_int(2,rNuTei)
       NLC( 7,NEQ) = LQe5

       ELM(1:NEMAX,1:4, 8,NEQ) =   1.5d0 * fem_int(2,rNuTeiEI)
       NLC( 8,NEQ) = LQi5

       ! Loss to diverter

       ELM(1:NEMAX,1:4, 9,NEQ) = -                  fem_int( 2,rNuL)
       NLC( 9,NEQ) = LQe5

       ELM(1:NEMAX,1:4,10,NEQ) =           PNeDIV * fem_int(-2,rNuL,PTeV)
       NLC(10,NEQ) = 0

       ELM(1:NEMAX,1:4,11,NEQ) = - 1.5D0          * fem_int( 2,rNuLTe)
       NLC(11,NEQ) = LQe5

!!unstable       ELM(1:NEMAX,1:4,12,NEQ) =  1.5D0 * PTeDIV * fem_int(-2,rNuLTe,PNeV)
!!unstable       NLC(12,NEQ) = 0
       ELM(1:NEMAX,1:4,12,NEQ) =   1.5D0 * PTeDIV * fem_int( 2,rNuLTe)
       NLC(12,NEQ) = LQe1

       ! Ionization loss of n01, n02 and n03

       ELM(1:NEMAX,1:4,13,NEQ) = - (EION * 1.D-3) * fem_int(2,rNuIN0)
       NLC(13,NEQ) = LQn1

       ELM(1:NEMAX,1:4,14,NEQ) = - (EION * 1.D-3) * fem_int(2,rNuIN0)
       NLC(14,NEQ) = LQn2

       ELM(1:NEMAX,1:4,15,NEQ) = - (EION * 1.D-3) * fem_int(2,rNuIN0) * BeamSW
       NLC(15,NEQ) = LQn3

       ! Collisional NBI heating (Perp + Tan)

       ELM(1:NEMAX,1:4,16,NEQ) = Eb * fem_int(-2,SNB,PNBcol_e)
       NLC(16,NEQ) = 0

       ! Collisional heating with beam

!!oldNBI      ELM(1:NEMAX,1:4,21,NEQ) = - 0.5D0 * AMB * (PNBCD * Vb) / rKeV &
!!oldNBI           &                * fem_int(2,rNubeBE)
!!oldNBI      NLC(21,NEQ) = LQe4
!!oldNBI
!!oldNBI      ELM(1:NEMAX,1:4,22,NEQ) =   0.5D0 * AMB * (PNBCD * Vb) / rKeV &
!!oldNBI           &                * fem_int(2,rNube) * BeamSW
!!oldNBI      NLC(22,NEQ) = LQb4

       ! Simpified Alpha heating

       ELM(1:NEMAX,1:4,17,NEQ) = 1.D0 / (1.D20 * rKeV) * fem_int(-1,PALFe)
       NLC(17,NEQ) = 0
       
       ! Direct heating (RF)

       ELM(1:NEMAX,1:4,18,NEQ) =   1.D0 / (1.D20 * rKeV) * fem_int(-1,PRFe)
       NLC(18,NEQ) = 0

       ! Radiation loss (Bremsstrahlung)

       ELM(1:NEMAX,1:4,19,NEQ) = - 1.D0 / (1.D20 * rKeV) * fem_int(-1,PBr)
       NLC(19,NEQ) = 0

       !  Diffusion of electrons (***AF 2008-06-08)

       ELM(1:NEMAX,1:4,20,NEQ) = - 4.D0 * fem_int(18,DMAGe)
       NLC(20,NEQ) = LQe5

       NLCMAX(NEQ) = 20
    ELSE

       !  Fixed temperature profile

       ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT
       NLC(0,NEQ) = LQe5

       NLCMAX(NEQ) = 0
    END IF

    RETURN
  END SUBROUTINE LQe5CC

!***************************************************************
!
!   Ion Density Equation
!
!***************************************************************

  SUBROUTINE LQi1CC

    integer(4) :: NEQ = LQi1

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQi1

    ! Divergence

    ELM(1:NEMAX,1:4,1,NEQ) = - 2.D0 * fem_int(4)
    NLC(1,NEQ) = LQi2

    ! Ionization of n01, n02 and n03

    ELM(1:NEMAX,1:4,2,NEQ) =     1.D0 / PZ * fem_int(2,rNuIN0)
    NLC(2,NEQ) = LQn1

    ELM(1:NEMAX,1:4,3,NEQ) =     1.D0 / PZ * fem_int(2,rNuIN0) * ThntSW
    NLC(3,NEQ) = LQn2

    ELM(1:NEMAX,1:4,4,NEQ) =     1.D0 / PZ * fem_int(2,rNuIN0) * BeamSW
    NLC(4,NEQ) = LQn3

    ! Loss to divertor

    ELM(1:NEMAX,1:4,5,NEQ) = -   1.D0 / PZ * fem_int(2,rNuL)
    NLC(5,NEQ) = LQe1

    ELM(1:NEMAX,1:4,6,NEQ) =  PNeDIV / PZ * fem_int(-1,rNuL)
    NLC(6,NEQ) = 0

    ! Particle source from beam ion

    ELM(1:NEMAX,1:4,7,NEQ) =   fem_int(2,rNuB)
    NLC(7,NEQ) = LQb1

    ! Particle source from ripple trapped beam ions

    ELM(1:NEMAX,1:4,8,NEQ) =   fem_int(28,rNuB,rip_rat) * RpplSW
    NLC(8,NEQ) = LQr1

    ! NBI kick up ions (Charge exchange)

    ELM(1:NEMAX,1:4,9,NEQ) = - RatCX * fem_int(-1,SNBi)
    NLC(9,NEQ) = 0

    ! Parallel loss reduction due to the potential
    ! induced by the parallel loss of the beam ions

    ELM(1:NEMAX,1:4,10,NEQ) =   fem_int(2,rNuLB) * BeamSW
    NLC(10,NEQ) = LQb1

    ! Loss cone loss

    ELM(1:NEMAX,1:4,11,NEQ) =   fem_int(-1,SiLC)
    NLC(11,NEQ) = 0
 
    ! Ion orbit loss

    ELM(1:NEMAX,1:4,12,NEQ) = - fem_int(2,rNuOL)
    NLC(12,NEQ) = LQi1

    !  Diffusion of ions (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,13,NEQ) = - 4.D0 * fem_int(18,DMAGi)
    NLC(13,NEQ) = LQi1

    NLCMAX(NEQ) = 13
    RETURN
  END SUBROUTINE LQi1CC

!***************************************************************
!
!   Ion Radial Flow (SUPG)
!
!***************************************************************
  
  SUBROUTINE LQi2CC

    integer(4) :: NEQ = LQi2

    ! Ns*Usr(0) : fixed

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT &
         &                 + fem_int(8) * fem_int(0) * invDT
    NLC(0,NEQ) = LQi2

    ! Advection

    ELM(1:NEMAX,1:4,1,NEQ) = - 2.D0 * fem_int( 3,RUirV) + fem_int(2,UirVR) &
         &                 +(- 2.D0 * fem_int(10,RUirV) + fem_int(9,UirVR)) * fem_int(0)
    ELM(1      ,1:4,1,NEQ) = - 2.D0 * fem_int_point( 3,0,RUirV) + fem_int_point(2,0,UirVR) &
         &                 +(- 2.D0 * fem_int_point(10,0,RUirV) + fem_int_point(9,0,UirVR)) &
         &                 * fem_int_point(0,1)
    NLC(1,NEQ) = LQi2

    ! Centrifugal force

    ELM(1:NEMAX,1:4,2,NEQ) = fem_int(2,UithVR) &
         &                 + fem_int(9,UithVR) * fem_int(0)
    ELM(1      ,1:4,2,NEQ) = fem_int_point(2,0,UithVR) &
         &                 + fem_int_point(9,0,UithVR) * fem_int_point(0,1)
    NLC(2,NEQ) = LQi3

    ! Pressure gradient force

    ELM(1:NEMAX,1:4,3,NEQ) = - 2.D0 * rKeV / AMI * fem_int(17,UNITY) &
         &                   - 2.D0 * rKeV / AMI * fem_int(18,UNITY) * fem_int(0)
    NLC(3,NEQ) = LQi5

    ! Radial E force

    ELM(1:NEMAX,1:4,4,NEQ) = - 2.D0 * (PZ * AEE / AMI) * fem_int(17,PNiV) &
         &                   - 2.D0 * (PZ * AEE / AMI) * fem_int(18,PNiV) * fem_int(0)
    NLC(4,NEQ) = LQm1

    ! v x B force

    ELM(1:NEMAX,1:4,5,NEQ) =          (PZ * AEE / AMI) * fem_int( 2,BphV) &
         &                 +          (PZ * AEE / AMI) * fem_int( 9,BphV) * fem_int(0)
    NLC(5,NEQ) = LQi3

    ELM(1:NEMAX,1:4,6,NEQ) =   2.D0 * (PZ * AEE / AMI) * fem_int(16,AphV) &
         &                 +   2.D0 * (PZ * AEE / AMI) * fem_int(20,UNITY,AphV) * fem_int(0)
    NLC(6,NEQ) = LQi4

    NLCMAX(NEQ) = 6

    IF(MDFIXT /= 0) THEN
       ELM(1:NEMAX,1:4,3,NEQ) = - 2.D0 * rKeV / AMI &
            &                 *(  fem_int(17,PTiV) + fem_int(18,PTiV) * fem_int(0) &
            &                   + fem_int(16,PTiV) + fem_int(34,PSI,PTiV) * fem_int(0))
       NLC(3,NEQ) = LQi1
    END IF

    RETURN
  END SUBROUTINE LQi2CC

!***************************************************************
!
!   Ion Poloidal Flow
!
!***************************************************************

  SUBROUTINE LQi3CC

    integer(4) :: NEQ = LQi3

    ! Ni*UiTheta(0) : 0

    ELM(1:NEMAX,1:4, 0,NEQ) = fem_int(1) * invDT
    NLC( 0,NEQ) = LQi3

    ! Advection

    ELM(1:NEMAX,1:4, 1,NEQ) = - 2.D0 * fem_int(3,RUirV)
    NLC( 1,NEQ) = LQi3

    ! Viscosity force

!    ELM(1:NEMAX,1:4, 2,NEQ) = - 4.D0 * fem_int(18,rMui) &
!         &                    - 4.D0 * fem_int(39,rMuiNi,dPNiV) &
!         &                    - 4.D0 * fem_int( 3,rMui)
    ELM(1:NEMAX,1:4, 2,NEQ) = - 4.D0 * fem_int(18,rMui) &
         &                    - 4.D0 * fem_int( 3,rMui)
    NLC( 2,NEQ) = LQi3

    ELM(1:NEMAX,1:4, 3,NEQ) =   4.D0 * fem_int(41,rMui,RUithV)
    NLC( 3,NEQ) = LQi1

    ! Poloidal E force

    ELM(1:NEMAX,1:4, 4,NEQ) = - (PZ * AEE / AMI) * fem_int(2,PNiV)
    NLC( 4,NEQ) = LQm2

    ! v x B force

    ELM(1:NEMAX,1:4, 5,NEQ) = - (PZ * AEE / AMI) * fem_int(2,BphV)
    NLC( 5,NEQ) = LQi2

    ! Neoclassical viscosity force

    ELM(1:NEMAX,1:4, 6,NEQ) = - fem_int(2,rNuiNC)
    NLC( 6,NEQ) = LQi3

    ! Collisional friction force

    ELM(1:NEMAX,1:4, 7,NEQ) = - (AME / AMI) * fem_int(2,rNuei1EI)
    NLC( 7,NEQ) = LQi3

    ELM(1:NEMAX,1:4, 8,NEQ) =   (AME / AMI) * fem_int(2,rNuei1)
    NLC( 8,NEQ) = LQe3

    ELM(1:NEMAX,1:4, 9,NEQ) =   2.D0 * (AME / AMI) * fem_int(37,rNuei2BthEI,AphV)
    NLC( 9,NEQ) = LQi4

    ELM(1:NEMAX,1:4,10,NEQ) = - 2.D0 * (AME / AMI) * fem_int(37,rNuei2Bth,AphV)
    NLC(10,NEQ) = LQe4

    ! Collisional friction with beam ions

    ELM(1:NEMAX,1:4,11,NEQ) = - (AMB / AMI) * fem_int(2,rNubiBI)
    NLC(11,NEQ) = LQi3

    ELM(1:NEMAX,1:4,12,NEQ) =   (AMB / AMI) * fem_int(2,rNubi) * BeamSW
    NLC(12,NEQ) = LQb3

    ! Turbulent particle transport driver

!!$    ELM(1:NEMAX,1:4,13,NEQ) = - 2.D0 * AEE / AMI * fem_int(38,De,BphV)
!!$    NLC(13,NEQ) = LQe1
!!$
!!$    ELM(1:NEMAX,1:4,14,NEQ) = - 2.D0 * AEE / AMI * fem_int(38,De,BphV/PTeV) * FSVAHLT
!!$    NLC(14,NEQ) = LQe5
!!$
!!$    ELM(1:NEMAX,1:4,15,NEQ) =   2.D0 * AEE / AMI * fem_int(38,De,BphV) * FSVAHLT
!!$    NLC(15,NEQ) = LQe1
!!$
    ELM(1:NEMAX,1:4,13,NEQ) =   1.D0 / AMI * fem_int(2,FWthe)
    NLC(13,NEQ) = LQe3

!    ELM(1:NEMAX,1:4,14,NEQ) = - 1.D0 / AMI * fem_int(28,FWthe,rat_ei)
!    NLC(14,NEQ) = LQi3

    ELM(1:NEMAX,1:4,15,NEQ) = - 1.D0 / AMI * fem_int(15,FWpheBB)
    NLC(15,NEQ) = LQe4

!    ELM(1:NEMAX,1:4,16,NEQ) =   1.D0 / AMI * fem_int(44,FWpheBB,rat_ei)
!    NLC(16,NEQ) = LQi4

    ! Wave interaction force (electron driven)

    ELM(1:NEMAX,1:4,17,NEQ) = - 1.D0 / AMI * fem_int(44,FWthpheB,WPM)
    NLC(17,NEQ) = LQe1

    ! Contribution of off-diagonal term due to heat flux

    ELM(1:NEMAX,1:4,18,NEQ) =   2.D0 * (rKeV / (AMI * AEE)) * fem_int(17,FWthphe) &
         &                           * (1.D0 - FSVAHLT)
    NLC(18,NEQ) = LQe5

    ELM(1:NEMAX,1:4,19,NEQ) = - 2.D0 * (rKeV / (AMI * AEE)) * fem_int(38,FWthphe,PTeV) &
         &                           * (1.D0 - FSVAHLT)
    NLC(19,NEQ) = LQe1

    ELM(1:NEMAX,1:4,20,NEQ) = - 2.D0 * (1.D0 /  AMI       ) * fem_int(37,FWthphe,PhiV) * FSVAHLE
    NLC(20,NEQ) = LQe1

    ! Ad hoc turbulent pinch term

    ELM(1:NEMAX,1:4,21,NEQ) =   1.D0 / AMI * fem_int(15,FVpch)
    NLC(21,NEQ) = LQe1

!!ion       ! Turbulent particle transport driver (ion driven)
!!ion
!!ion       ELM(1:NEMAX,1:4,13,NEQ) = - 1.D0 / AMI * fem_int(2,FWthi)
!!ion       NLC(13,NEQ) = LQi3
!!ion
!!ion       ELM(1:NEMAX,1:4,14,NEQ) =   1.D0 / AMI * fem_int(2,FWthi)
!!ion       NLC(14,NEQ) = LQe3
!!ion
!!ion       ELM(1:NEMAX,1:4,15,NEQ) =   1.D0 / AMI * fem_int(15,FWphiBB)
!!ion       NLC(15,NEQ) = LQi4
!!ion
!!ion       ELM(1:NEMAX,1:4,16,NEQ) = - 1.D0 / AMI * fem_int(15,FWphiBB)
!!ion       NLC(16,NEQ) = LQe4
!!ion
!!ion       ELM(1:NEMAX,1:4,17,NEQ) = - 2.D0 * (rKeV / (AMI * AEE)) * fem_int(37,FWthphi,PTiV) &
!!ion            &                           * (1.D0 - FSVAHLT) &
!!ion            &                    + 2.D0 * (1.D0 /  AMI       ) * fem_int(37,FWthphi,PhiV) * FSVAHLE
!!ion       NLC(17,NEQ) = LQi1
!!ion
!!ion       ELM(1:NEMAX,1:4,18,NEQ) = - 1.D0 / AMI * fem_int(44,FWthphiB,WPM)
!!ion       NLC(18,NEQ) = LQi1

    ! Loss to divertor

    ELM(1:NEMAX,1:4,22,NEQ) = - 2.D0 * fem_int(2,rNuL) * fact
    NLC(22,NEQ) = LQi3

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,23,NEQ) = - fem_int(2,rNu0i)
    NLC(23,NEQ) = LQi3

    ! Charge exchange force

    ELM(1:NEMAX,1:4,24,NEQ) = - fem_int(2,rNuiCX)
    NLC(24,NEQ) = LQi3

    ! Loss cone loss

    ELM(1:NEMAX,1:4,25,NEQ) = fem_int(-1,SiLCth)
    NLC(25,NEQ) = 0

    ! Ion orbit loss

    ELM(1:NEMAX,1:4,26,NEQ) = - fem_int(2,rNuOL)
    NLC(26,NEQ) = LQi3

    ! Helical Neoclassical viscosity force

    ELM(1:NEMAX,1:4,27,NEQ) = - fem_int(2,rNuiHLthth)
    NLC(27,NEQ) = LQi3

    ELM(1:NEMAX,1:4,28,NEQ) =   fem_int(15,rNuiHLthph)
    NLC(28,NEQ) = LQi4

    !  Diffusion of ions (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,29,NEQ) = - 4.D0 * fem_int(18,DMAGi)
    NLC(29,NEQ) = LQi3

    !  Poloidal torque due to neoclassical heat flux

    ELM(1:NEMAX,1:4,30,NEQ) =   FSNCPL / 1.D20 * fem_int(-2,R,FQith)
    NLC(30,NEQ) = 0

!    ELM(1:NEMAX,1:4,30,NEQ) = - FSNCPL * 2.D0 * rKeV / (PZ * AEE) * fem_int(17,rNui2NC)
!    NLC(30,NEQ) = LQi5
!
!    ELM(1:NEMAX,1:4,31,NEQ) =   FSNCPL * 2.D0 * rKeV / (PZ * AEE) * fem_int(37,rNui2NCN,PNiV)
!    NLC(31,NEQ) = LQi5
    
    !  Virtual torque input

    ELM(1:NEMAX,1:4,32,NEQ) =   1.D0 / (AMI * 1.D20) * fem_int(-1,Tqp)
    NLC(32,NEQ) = 0

    ! Ns*UsTheta(NRMAX) : 0

    NLCMAX(NEQ) = 32

    IF(MDFIXT /= 0) THEN
       ELM(1:NEMAX,1:4,18,NEQ) =   2.D0 * (rKeV / (AMI * AEE)) * fem_int(37,FWthphe,PTeV) &
            &                           * (1.D0 - FSVAHLT)
       NLC(18,NEQ) = LQe1

       ELM(1:NEMAX,1:4,19,NEQ) = 0.D0
       NLC(19,NEQ) = 0

!       ELM(1:NEMAX,1:4,30,NEQ) = - FSNCPL * 2.D0 * rKeV / (PZ * AEE) * fem_int(37,rNui2NC,PTiV)
!       NLC(30,NEQ) = LQi1
!
!       ELM(1:NEMAX,1:4,31,NEQ) = 0.D0
!       NLC(31,NEQ) = 0
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

    ! Uiphi'(0) : 0

    ELM(1:NEMAX,1:4, 0,NEQ) = fem_int(1) * invDT
    NLC( 0,NEQ) = LQi4

    ! Advection

    ELM(1:NEMAX,1:4, 1,NEQ) = - 2.D0 * fem_int(3,RUirV)
    NLC( 1,NEQ) = LQi4

    ! Viscosity force

!    ELM(1:NEMAX,1:4, 2,NEQ) = - 4.D0 *(fem_int(18,rMui) + fem_int(39,rMuiNi,dPNiV))
    ELM(1:NEMAX,1:4, 2,NEQ) = - 4.D0 * fem_int(18,rMui)
    NLC( 2,NEQ) = LQi4

    ELM(1:NEMAX,1:4, 3,NEQ) =   4.D0 * fem_int(41,rMui,UiphV)
    NLC( 3,NEQ) = LQi1

    ! Toroidal E force

    ELM(1:NEMAX,1:4, 4,NEQ) = - (PZ * AEE / AMI) * fem_int(2,PNiV)
    NLC( 4,NEQ) = LQm3

    ! v x B force

    ELM(1:NEMAX,1:4, 5,NEQ) = - 2.D0 * (PZ * AEE / AMI) * fem_int(6,AphV)
    NLC( 5,NEQ) = LQi2

    ! Collisional friction with bulk ions

    ELM(1:NEMAX,1:4, 6,NEQ) = - (AME / AMI) * fem_int(2,rNuei3EI)
    NLC( 6,NEQ) = LQi4

    ELM(1:NEMAX,1:4, 7,NEQ) =   (AME / AMI) * fem_int(2,rNuei3)
    NLC( 7,NEQ) = LQe4

    ELM(1:NEMAX,1:4, 8,NEQ) =   2.D0 * (AME / AMI) * fem_int(29,rNuei2BthEI,AphV)
    NLC( 8,NEQ) = LQi3

    ELM(1:NEMAX,1:4, 9,NEQ) = - 2.D0 * (AME / AMI) * fem_int(29,rNuei2Bth,AphV)
    NLC( 9,NEQ) = LQe3

    ! Collisional friction with beam ions

    ELM(1:NEMAX,1:4,10,NEQ) = - (AMB / AMI) * fem_int(2,rNubiBI)
    NLC(10,NEQ) = LQi4

    ELM(1:NEMAX,1:4,11,NEQ) =   (AMB / AMI) * fem_int(2,rNubi) * BeamSW
    NLC(11,NEQ) = LQb4

    ! Turbulent particle transport driver (electron driven)

!!$    ELM(1:NEMAX,1:4,12,NEQ) =   2.D0 * AEE / AMI * fem_int(38,De,BthV)
!!$    NLC(12,NEQ) = LQe1
!!$
!!$    ELM(1:NEMAX,1:4,13,NEQ) =   2.D0 * AEE / AMI * fem_int(38,De,BthV/PTeV) * FSVAHLT
!!$    NLC(13,NEQ) = LQe5
!!$
!!$    ELM(1:NEMAX,1:4,14,NEQ) = - 2.D0 * AEE / AMI * fem_int(38,De,BthV) * FSVAHLT
!!$    NLC(14,NEQ) = LQe1
!!$
    ELM(1:NEMAX,1:4,12,NEQ) = - 1.D0 / AMI * fem_int(2,FWpheBB)
    NLC(12,NEQ) = LQe3

!    ELM(1:NEMAX,1:4,13,NEQ) =   1.D0 / AMI * fem_int(28,FWpheBB,rat_ei)
!    NLC(13,NEQ) = LQi3

    ELM(1:NEMAX,1:4,14,NEQ) =   1.D0 / AMI * fem_int(2,FWpheBB2)
    NLC(14,NEQ) = LQe4

!    ELM(1:NEMAX,1:4,15,NEQ) = - 1.D0 / AMI * fem_int(28,FWpheBB2,rat_ei) 
!    NLC(15,NEQ) = LQi4

    ! Wave interaction force (electron driven)

    ELM(1:NEMAX,1:4,16,NEQ) =   1.D0 / AMI * fem_int(44,FWpheB,WPM)
    NLC(16,NEQ) = LQe1

    ! Contribution of off-diagonal term due to heat flux

    ELM(1:NEMAX,1:4,17,NEQ) = - 2.D0 * (rKeV / (AMI * AEE)) * fem_int(17,FWphe) &
         &                           * (1.D0 - FSVAHLT)
    NLC(17,NEQ) = LQe5

    ELM(1:NEMAX,1:4,18,NEQ) = + 2.D0 * (rKeV / (AMI * AEE)) * fem_int(38,FWphe,PTeV) &
         &                           * (1.D0 - FSVAHLT)
    NLC(18,NEQ) = LQe1

    ELM(1:NEMAX,1:4,19,NEQ) = + 2.D0 * (1.D0 /  AMI       ) * fem_int(37,FWphe,PhiV) * FSVAHLE
    NLC(19,NEQ) = LQe1

    ! Ad hoc turbulent pinch term

    ELM(1:NEMAX,1:4,20,NEQ) =   2.D0 / AMI * fem_int(36,AphV,FVpchph)
    NLC(20,NEQ) = LQe1

!!ion       ! Turbulent particle transport driver (ion driven)
!!ion
!!ion       ELM(1:NEMAX,1:4,12,NEQ) =   1.D0 / AMI * fem_int(2,FWphiBB)
!!ion       NLC(12,NEQ) = LQi3
!!ion
!!ion       ELM(1:NEMAX,1:4,13,NEQ) = - 1.D0 / AMI * fem_int(2,FWphiBB)
!!ion       NLC(13,NEQ) = LQe3
!!ion
!!ion       ELM(1:NEMAX,1:4,14,NEQ) = - 1.D0 / AMI * fem_int(2,FWphiBB2)
!!ion       NLC(14,NEQ) = LQi4
!!ion
!!ion       ELM(1:NEMAX,1:4,15,NEQ) =   1.D0 / AMI * fem_int(2,FWphiBB2)
!!ion       NLC(15,NEQ) = LQe4
!!ion
!!ion       ELM(1:NEMAX,1:4,16,NEQ) =   2.D0 * (rKeV / (AMI * AEE)) * fem_int(37,FWphiBB,PTiV) &
!!ion            &                           * (1.D0 - FSVAHLT) &
!!ion            &                    - 2.D0 * (1.D0 /  AMI       ) * fem_int(37,FWphiBB,PhiV) * FSVAHLE
!!ion       NLC(16,NEQ) = LQi1
!!ion
!!ion       ELM(1:NEMAX,1:4,17,NEQ) =   1.D0 / AMI * fem_int(44,FWphiB,WPM)
!!ion       NLC(17,NEQ) = LQi1

    ! Loss to divertor

    ELM(1:NEMAX,1:4,21,NEQ) = - 2.D0 * fem_int(2,rNuL) * fact
    NLC(21,NEQ) = LQi4

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,22,NEQ) = - fem_int(2,rNu0i)
    NLC(22,NEQ) = LQi4

    ! Charge exchange force

    ELM(1:NEMAX,1:4,23,NEQ) = - fem_int(2,rNuiCX)
    NLC(23,NEQ) = LQi4

    ! Loss cone loss

    ELM(1:NEMAX,1:4,24,NEQ) =   fem_int(-1,SiLCph)
    NLC(24,NEQ) = 0

    ! Ion orbit loss

    ELM(1:NEMAX,1:4,25,NEQ) = - fem_int(2,rNuOL)
    NLC(25,NEQ) = LQi4

   ! Helical Neoclassical viscosity force

    ELM(1:NEMAX,1:4,26,NEQ) =   fem_int(2,rNuiHLphth)
    NLC(26,NEQ) = LQi3

    ELM(1:NEMAX,1:4,27,NEQ) = - fem_int(15,rNuiHLphph)
    NLC(27,NEQ) = LQi4

    !  Diffusion of ions (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,28,NEQ) = - 4.D0 * fem_int(18,DMAGi)
    NLC(28,NEQ) = LQi4

    !  Virtual torque input

    ELM(1:NEMAX,1:4,29,NEQ) =   1.D0 / (AMI * RR * 1.D20) * fem_int(-1,Tqt)
    NLC(29,NEQ) = 0

    NLCMAX(NEQ) = 29

    IF(MDFIXT /= 0) THEN
       ELM(1:NEMAX,1:4,17,NEQ) = - 2.D0 * (rKeV / (AMI * AEE)) * fem_int(37,FWphe,PTeV) &
            &                           * (1.D0 - FSVAHLT)
       NLC(17,NEQ) = LQe1

       ELM(1:NEMAX,1:4,18,NEQ) = 0.D0
       NLC(18,NEQ) = 0
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
       ELM(1:NEMAX,1:4, 0,NEQ) =   1.5D0 * fem_int(1) * invDT
       NLC( 0,NEQ) = LQi5

       ! Advection

       ELM(1:NEMAX,1:4, 1,NEQ) = - 5.D0 * fem_int(3,RUirV)
       NLC( 1,NEQ) = LQi5

       ! Conduction transport

       ELM(1:NEMAX,1:4, 2,NEQ) = - 4.D0 * fem_int(18,Chii1)
       NLC( 2,NEQ) = LQi5

       ELM(1:NEMAX,1:4, 3,NEQ) =   4.D0 * fem_int(41,Chii2,PTiV)
       NLC( 3,NEQ) = LQi1

       ! Redundant heat convection term

       ELM(1:NEMAX,1:4,4,NEQ) = 2.D0 * fem_int(5,RUirV)
       NLC(4,NEQ) = LQi5

       ! Joule heating

       ELM(1:NEMAX,1:4, 5,NEQ) =   PZ * AEE / rKeV * fem_int(2,EthVR)
       ELM(1      ,1:4, 5,NEQ) =   PZ * AEE / rKeV * fem_int_point(2,0,EthVR)
       NLC( 5,NEQ) = LQi3

       ELM(1:NEMAX,1:4, 6,NEQ) =   PZ * AEE / rKeV * fem_int(2,EphV)
       NLC( 6,NEQ) = LQi4

       ! Collisional transfer with electrons (Energy equilibration)

       ELM(1:NEMAX,1:4, 7,NEQ) = - 1.5d0 * fem_int(2,rNuTeiEI)
       NLC( 7,NEQ) = LQi5

       ELM(1:NEMAX,1:4, 8,NEQ) =   1.5d0 * fem_int(2,rNuTei)
       NLC( 8,NEQ) = LQe5

       ! Heating due to beam momentum deposition

!!oldNBI      ELM(1:NEMAX,1:4, 9,NEQ) = - 0.5D0 * AMB * (PNBCD * Vb) / rKeV &
!!oldNBI           &                * fem_int(2,rNubiBI)
!!oldNBI      NLC( 9,NEQ) = LQi4
!!oldNBI
!!oldNBI      ELM(1:NEMAX,1:4,10,NEQ) =   0.5D0 * AMB * (PNBCD * Vb) / rKeV &
!!oldNBI           &                * fem_int(2,rNubi) * BeamSW
!!oldNBI      NLC(10,NEQ) = LQb4

       ELM(1:NEMAX,1:4, 9,NEQ) = AMB * Vb / rKeV * fem_int(28,BthBNi,MNB)
       ELM(1      ,1:4, 9,NEQ) = AMB * Vb / rKeV * fem_int_point(28,0,BthBNi,MNB)
       NLC( 9,NEQ) = LQi3

       ELM(1:NEMAX,1:4,10,NEQ) = AMB * Vb / rKeV * fem_int(28,BphBNi,MNB)
       NLC(10,NEQ) = LQi4

       ! Loss to diverter

       ELM(1:NEMAX,1:4,11,NEQ) = -         1.D0   / PZ * fem_int( 2,rNuL)
       NLC(11,NEQ) = LQi5

       ELM(1:NEMAX,1:4,12,NEQ) =           PNeDIV / PZ * fem_int(-2,rNuL,PTiV)
       NLC(12,NEQ) = 0

       ELM(1:NEMAX,1:4,13,NEQ) = - 1.5D0               * fem_int( 2,rNuLTi)
       NLC(13,NEQ) = LQi5

       ELM(1:NEMAX,1:4,14,NEQ) =   1.5D0 * PTiDIV      * fem_int(-2,rNuLTi,PNiV)
       NLC(14,NEQ) = 0

       ! Ionization heating of n01, n02 and n03

       ELM(1:NEMAX,1:4,15,NEQ) = 1.5D0 / PZ * fem_int(28,rNuIN0,PT01V)
       NLC(15,NEQ) = LQn1

       ELM(1:NEMAX,1:4,16,NEQ) = 1.5D0 / PZ * fem_int(28,rNuIN0,PT02V)
       NLC(16,NEQ) = LQn2

       ELM(1:NEMAX,1:4,17,NEQ) = 1.5D0 / PZ * fem_int(28,rNuIN0,PT03V) * BeamSW
       NLC(17,NEQ) = LQn3

       ! Charge exchange loss due to slow neutrals
       !   (Thermal neutrals are assumed to have same temperature with ions.)

       ELM(1:NEMAX,1:4,18,NEQ) = - 1.5D0 * fem_int( 2,rNuiCXT)
       NLC(18,NEQ) = LQi5

       ELM(1:NEMAX,1:4,19,NEQ) =   1.5D0 * fem_int(-1,rNuCXN1)
       NLC(19,NEQ) = 0

       ! Collisional NBI heating (Perp + Tan)

       ELM(1:NEMAX,1:4,20,NEQ) = Eb * fem_int(-2,SNB,PNBcol_i)
       NLC(20,NEQ) = 0

       ! Simplified Alpha heating

       ELM(1:NEMAX,1:4,21,NEQ) = 1.D0 / (1.D20 * rKeV) * fem_int(-1,PALFi)
       NLC(21,NEQ) = 0

       ! Direct heating (RF)

       ELM(1:NEMAX,1:4,22,NEQ) = 1.D0 / (1.D20 * rKeV) * fem_int(-1,PRFi)
       NLC(22,NEQ) = 0

       !  Diffusion of ions (***AF 2008-06-08)

       ELM(1:NEMAX,1:4,23,NEQ) = - 4.D0 * fem_int(18,DMAGi)
       NLC(23,NEQ) = LQi5

       NLCMAX(NEQ) = 23
    ELSE

       !  Fixed temperature profile

       ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT
       NLC(0,NEQ) = LQi5

       NLCMAX(NEQ) = 0
    END IF

    RETURN
  END SUBROUTINE LQi5CC

!***************************************************************
!
!   Beam Ion Density
!
!***************************************************************

  SUBROUTINE LQb1CC

    integer(4) :: NEQ = LQb1

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQb1

    ! NBI particle source (Both charge exchange and ionization)

    ELM(1:NEMAX,1:4,1,NEQ) =   fem_int(-1,SNBb)
    NLC(1,NEQ) = 0

    ! Extracted NBI perpendicular component fallen into the ripple well region

    ELM(1:NEMAX,1:4,2,NEQ) = - fem_int(-2,SNBPDi,rip_rat)
    NLC(2,NEQ) = 0

    ! Relaxation to thermal ions

    ELM(1:NEMAX,1:4,3,NEQ) = - fem_int(2,rNuB)
    NLC(3,NEQ) = LQb1

    ! Loss to divertor

    ELM(1:NEMAX,1:4,4,NEQ) = - fem_int(2,rNuLB) * BeamSW
    NLC(4,NEQ) = LQb1

    ! Ripple trapped beam ions collision with otherwise beam ions

    ELM(1:NEMAX,1:4,5,NEQ) = - fem_int(28,rNubrp2,rip_rat) * RpplSW
    NLC(5,NEQ) = LQb1

    ELM(1:NEMAX,1:4,6,NEQ) =   fem_int(28,rNubrp1,rip_rat) * RpplSW
    NLC(6,NEQ) = LQr1

!!rp_conv    ELM(1:NEMAX,1:4,6,NEQ) =  fem_int(-2,PNbrpLV,rNubLL)
!!rp_conv    NLC(6,NEQ) = 0

    ! Ripple diffusion

    ELM(1:NEMAX,1:4,7,NEQ) = - 4.d0 * fem_int(18,Dbrpft)
    NLC(7,NEQ) = LQb1

    NLCMAX(NEQ) = 7
    RETURN
  END SUBROUTINE LQb1CC

!***************************************************************
!
!   Beam Ion Poloidal Flow
!
!***************************************************************

  SUBROUTINE LQb3CC

    integer(4) :: NEQ = LQb3

    ! Ubth(0) : 0

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQb3

    ! Poloidal E force

    ELM(1:NEMAX,1:4,1,NEQ) = - (PZ * AEE / AMB) * fem_int(2,PNbV)
    NLC(1,NEQ) = LQm2

    ! Collisional friction force with electrons

    ELM(1:NEMAX,1:4,2,NEQ) = - fem_int(2,rNube1)
    NLC(2,NEQ) = LQb3

    ELM(1:NEMAX,1:4,3,NEQ) =   fem_int(2,rNube1BE)
    NLC(3,NEQ) = LQe3

    ELM(1:NEMAX,1:4,4,NEQ) =   2.D0 * fem_int(37,rNube2Bth,AphV) * BeamSW
    NLC(4,NEQ) = LQb4

    ELM(1:NEMAX,1:4,5,NEQ) = - 2.D0 * fem_int(37,rNube2BthBE,AphV)
    NLC(5,NEQ) = LQe4

    ! Collisional friction force with ions

    ELM(1:NEMAX,1:4,6,NEQ) = - fem_int(2,rNubi)
    NLC(6,NEQ) = LQb3

    ELM(1:NEMAX,1:4,7,NEQ) =   fem_int(2,rNubiBI)
    NLC(7,NEQ) = LQi3

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,8,NEQ) = - fem_int(2,rNu0b) * BeamSW
    NLC(8,NEQ) = LQb3

    ! Charge exchange force

    ELM(1:NEMAX,1:4,9,NEQ) = - fem_int(2,rNubCX) * BeamSW
    NLC(9,NEQ) = LQb3

    ! NBI momentum source

    ELM(1:NEMAX,1:4,10,NEQ) =  fem_int(-2,RVbparath,MNB)
    NLC(10,NEQ) = 0

    ! Loss to divertor

    ELM(1:NEMAX,1:4,11,NEQ) = - fem_int(2,rNuLB) * BeamSW * fact
    NLC(11,NEQ) = LQb3

    ! Momentum loss due to collisional ripple trapping

    ELM(1:NEMAX,1:4,12,NEQ) = - fem_int(28,rNubrp2,rip_rat) * RpplSW
    NLC(12,NEQ) = LQb3

    ! Momentum diffusion arising from beam ion convective flux due to ripple

    ELM(1:NEMAX,1:4,13,NEQ) = - 4.D0 * fem_int(18,Dbrpft)
    NLC(13,NEQ) = LQb3

    ELM(1:NEMAX,1:4,14,NEQ) = - 4.D0 * fem_int(45,Dbrpft,RUbthV)
    NLC(14,NEQ) = LQb1

!!neo    ! Neoclassical viscosity force
!!neo
!!neo    ELM(1:NEMAX,1:4,15,NEQ) = - fem_int(2,rNuiNC)
!!neo    NLC(15,NEQ) = LQb3

    ! Ubth(NRMAX) : 0

    NLCMAX(NEQ) = 14
    RETURN
  END SUBROUTINE LQb3CC

!***************************************************************
!
!   Beam Ion Toroidal Flow
!
!***************************************************************
 
 SUBROUTINE LQb4CC

    integer(4) :: NEQ = LQb4

    ! - UbPhi(0)' : 0

    ELM(1:NEMAX,1:4,0,NEQ) =   fem_int(1) * invDT
    NLC(0,NEQ) = LQb4

    ! Toroidal E force

    ELM(1:NEMAX,1:4,1,NEQ) = - (PZ * AEE / AMB) * fem_int(2,PNbV)
    NLC(1,NEQ) = LQm3

    ! Collisional friction with electrons

    ELM(1:NEMAX,1:4,2,NEQ) = - fem_int(2,rNube3)
    NLC(2,NEQ) = LQb4

    ELM(1:NEMAX,1:4,3,NEQ) =   fem_int(2,rNube3BE)
    NLC(3,NEQ) = LQe4

    ELM(1:NEMAX,1:4,4,NEQ) =   2.D0 * fem_int(29,rNube2Bth,AphV)
    NLC(4,NEQ) = LQb3

    ELM(1:NEMAX,1:4,5,NEQ) = - 2.D0 * fem_int(29,rNube2BthBE,AphV)
    NLC(5,NEQ) = LQe3

    ! Collisional friction with ions

    ELM(1:NEMAX,1:4,6,NEQ) = - fem_int(2,rNubi)
    NLC(6,NEQ) = LQb4

    ELM(1:NEMAX,1:4,7,NEQ) =   fem_int(2,rNubiBI)
    NLC(7,NEQ) = LQi4

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,8,NEQ) = - fem_int(2,rNu0b) * BeamSW
    NLC(8,NEQ) = LQb4

    ! Charge exchange force

    ELM(1:NEMAX,1:4,9,NEQ) = - fem_int(2,rNubCX) * BeamSW
    NLC(9,NEQ) = LQb4

    ! NBI momentum source

    ELM(1:NEMAX,1:4,10,NEQ) =  fem_int(-2,Vbparaph,MNB)
    NLC(10,NEQ) = 0

    ! Loss to divertor

    ELM(1:NEMAX,1:4,11,NEQ) = - fem_int(2,rNuLB) * BeamSW * fact
    NLC(11,NEQ) = LQb4

    ! Momentum loss due to collisional ripple trapping

    ELM(1:NEMAX,1:4,12,NEQ) = - fem_int(28,rNubrp2,rip_rat) * RpplSW
    NLC(12,NEQ) = LQb4

    ! Momentum diffusion arising from beam ion convective flux due to ripple

    ELM(1:NEMAX,1:4,13,NEQ) = - 4.D0 * fem_int(18,Dbrpft)
    NLC(13,NEQ) = LQb4

    ELM(1:NEMAX,1:4,14,NEQ) = - 4.D0 * fem_int(45,Dbrpft,UbphV)
    NLC(14,NEQ) = LQb1

    ! Ubphi(NRMAX) : 0

    NLCMAX(NEQ) = 14
    RETURN
  END SUBROUTINE LQb4CC

!***************************************************************
!
!   Slow Neutral Transport: n01
!
!***************************************************************

  SUBROUTINE LQn1CC

    integer(4) :: NEQ = LQn1

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQn1

    !  Diffusion of neutrals

    ELM(1:NEMAX,1:4,1,NEQ) = - 4.D0 * fem_int(18,D01)
    NLC(1,NEQ) = LQn1

    ! Ionization

    ELM(1:NEMAX,1:4,2,NEQ) = - 1.D0 / PZ * fem_int(2,rNuIN0)
    NLC(2,NEQ) = LQn1

    ! Generation of thermal neutrals by charge exchange

    ELM(1:NEMAX,1:4,3,NEQ) = - fem_int(2,rNuCXN0)
    NLC(3,NEQ) = LQn1

    ! Recycling from divertor

    ELM(1:NEMAX,1:4,4,NEQ) =   rGamm0 / PZ * fem_int(2,rNuL)
    NLC(4,NEQ) = LQe1

    ELM(1:NEMAX,1:4,5,NEQ) = - rGamm0 * PNeDIV / PZ * fem_int(-1,rNuL)
    NLC(5,NEQ) = 0

    ! Gas puff

    ELM(1:NEMAX,1:4,6,NEQ) = 2.D0 * fem_int(-3,R,rGASPFA)
    NLC(6,NEQ) = 0

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

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQn2

    !  Diffusion of neutrals

    ELM(1:NEMAX,1:4,1,NEQ) = - 4.D0 * fem_int(18,D02)
    NLC(1,NEQ) = LQn2

    ! Ionization

    ELM(1:NEMAX,1:4,2,NEQ) = - 1.D0 / PZ * fem_int(2,rNuIN0)
    NLC(2,NEQ) = LQn2

    ! Generation of thermal neutrals by charge exchange

    ELM(1:NEMAX,1:4,3,NEQ) = fem_int(2,rNuCXN0)
    NLC(3,NEQ) = LQn1

!    ! NBI particle source (Charge exchange)
!
!    ELM(1:NEMAX,1:4,4,NEQ) = RatCX * fem_int(-1,SNBi)
!    NLC(4,NEQ) = 0

    NLCMAX(NEQ) = 3
    RETURN
  END SUBROUTINE LQn2CC

!***************************************************************
!
!   Halo Neutral Transport: n03
!
!***************************************************************

  SUBROUTINE LQn3CC

    integer(4) :: NEQ = LQn3

    ELM(1:NEMAX,1:4,0,NEQ) = fem_int(1) * invDT
    NLC(0,NEQ) = LQn3

    !  Diffusion of neutrals

    ELM(1:NEMAX,1:4,1,NEQ) = - 4.D0 * fem_int(18,D03)
    NLC(1,NEQ) = LQn3

    ! Ionization

    ELM(1:NEMAX,1:4,2,NEQ) = - 1.D0 / PZ * fem_int(2,rNuIN0)
    NLC(2,NEQ) = LQn3

    ! NBI particle source (Charge exchange)

    ELM(1:NEMAX,1:4,3,NEQ) = RatCX * fem_int(-1,SNBi)
    NLC(3,NEQ) = 0

    NLCMAX(NEQ) = 3
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
       SDbrpl = Dbrp(ne-1)*PSI(ne-1) + Dbrp(ne)*PSI(ne)
       if (RUbrpl == 0.d0) then
          coef = 0.d0
       elseif (SDbrpl == 0.d0) then
          coef = 1.d0 / sqrt(15.d0) * hpsi(ne)
       else
          peclet = 0.5d0 * RUbrpl * hpsi(ne) / SDbrpl
          coef = 0.5d0 * langevin(peclet) * hpsi(ne)
       end if
       coef=0.d0 !!! no SUPG

       ELM(NE,1:4,0,NEQ) =   fem_int_point(1,NE) * invDT &
            &               + fem_int_point(8,NE) * coef * invDT
       NLC(0,NEQ) = LQr1

       ! NBI perpendicular particle source (Both charge exchange and ionization)
       ! (Beam ions by perpendicular NBI have few parallel velocity, hence
       !  they can be easily trapped by ripple wells if they go into the ripple well region.)
       !  (M.H. Redi, et al., NF 35 (1995) 1191, p.1201 sixth line from the bottom
       !   at right-hand-side column)
       
       ELM(NE,1:4,1,NEQ) =(  fem_int_point(-1,NE,SNBPDi) &
            &                + fem_int_point(-8,NE,SNBPDi) * coef) * RpplSW
       NLC(1,NEQ) = 0

       ! Ripple trapped beam ions collision with otherwise beam ions

       ELM(NE,1:4,2,NEQ) =(  fem_int_point(2,NE,rNubrp2) &
            &               + fem_int_point(9,NE,rNubrp2) * coef) * RpplSW
       NLC(2,NEQ) = LQb1

       ELM(NE,1:4,3,NEQ) =(- fem_int_point(2,NE,rNubrp1) &
            &               - fem_int_point(9,NE,rNubrp1) * coef) * RpplSW
       NLC(3,NEQ) = LQr1

       ! Relaxation to thermal ions

       ELM(NE,1:4,4,NEQ) =(- fem_int_point(2,NE,rNuB) &
            &               - fem_int_point(9,NE,rNuB) * coef) * RpplSW
       NLC(4,NEQ) = LQr1

       ! Ripple loss transport (convective)

       ELM(NE,1:4,5,NEQ) =(- 2.d0 * fem_int_point( 3,NE,RUbrp) &
            &               - 2.d0 * fem_int_point(10,NE,RUbrp) * coef) * RpplSW
       NLC(5,NEQ) = LQr1

!!rp_conv       ELM(NE,1:4,5,NEQ) = - fem_int_point(2,NE,rNubL)
!!rp_conv       NLC(5,NEQ) = LQr1
!!rp_conv       ELM(NE,1:4,5,NEQ) = - fem_int_point(-2,NE,rNubL,PNbrpV)
!!rp_conv       NLC(5,NEQ) = 0

       ! Ripple loss transport (diffusive)

       ELM(NE,1:4,6,NEQ) = - 4.d0 * fem_int_point(18,NE,Dbrp) * RpplSW
       NLC(6,NEQ) = LQr1
    end do

    NLCMAX(NEQ) = 6
    RETURN
  END SUBROUTINE LQr1CC

!***************************************************************
!
!   Dirichlet condition
!
!***************************************************************

  SUBROUTINE BOUNDARY(NR,LQ,ID,VAL)

    integer(4), intent(in) :: NR, LQ, ID
    real(8), intent(in), optional :: VAL
    integer(4) :: NQ, NC, I
    integer(4), save :: IMAX(1:NQM)
    type list
       integer(4) :: IDXNC
       integer(4) :: IDXNQ
    end type list
    type(list), save :: IDX(1:NQM,0:1,1:50)

    IF(ID == 0) THEN
       ! Initialize ALC, BLC and CLC at NR
       ALC(0:NCM,LQ,NR) = 0.D0
       BLC(0:NCM,LQ,NR) = 0.D0
       CLC(0:NCM,LQ,NR) = 0.D0
       PLC(1:NCM,LQ,NR) = 0.D0

       IF(ICALA == 0) THEN
          I = 0
          DO NQ = 1, NQMAX
             DO NC = 1, NLCMAX(NQ)
                IF(NLCR(NC,NQ,NR) == LQ) THEN
                   I = I + 1
                   IDX(LQ,MOD(NR,NRMAX),I)%IDXNC = NC
                   IDX(LQ,MOD(NR,NRMAX),I)%IDXNQ = NQ
                END IF
             END DO
          END DO
          IMAX(LQ) = I
       ELSE
          IF(PRESENT(VAL)) THEN
             IF(NR == 0) THEN
                DO I = 1, IMAX(LQ)
                   NC = IDX(LQ,MOD(NR,NRMAX),I)%IDXNC
                   NQ = IDX(LQ,MOD(NR,NRMAX),I)%IDXNQ
                   PLC(NC,NQ,NR)   = PLC(NC,NQ,NR)   + BLC(NC,NQ,NR)   * VAL
                   PLC(NC,NQ,NR+1) = PLC(NC,NQ,NR+1) + CLC(NC,NQ,NR+1) * VAL
                END DO
             ELSE
                DO I = 1, IMAX(LQ)
                   NC = IDX(LQ,MOD(NR,NRMAX),I)%IDXNC
                   NQ = IDX(LQ,MOD(NR,NRMAX),I)%IDXNQ
                   PLC(NC,NQ,NR)   = PLC(NC,NQ,NR)   + BLC(NC,NQ,NR)   * VAL
                   PLC(NC,NQ,NR-1) = PLC(NC,NQ,NR-1) + ALC(NC,NQ,NR-1) * VAL
                END DO
             END IF
          END IF

          IF(NR == 0) THEN
             DO I = 1, IMAX(LQ)
                NC = IDX(LQ,MOD(NR,NRMAX),I)%IDXNC
                NQ = IDX(LQ,MOD(NR,NRMAX),I)%IDXNQ
                BLC(NC,NQ,NR)   = 0.D0
                CLC(NC,NQ,NR+1) = 0.D0
             END DO
          ELSE
             DO I = 1, IMAX(LQ)
                NC = IDX(LQ,MOD(NR,NRMAX),I)%IDXNC
                NQ = IDX(LQ,MOD(NR,NRMAX),I)%IDXNQ
                BLC(NC,NQ,NR)   = 0.D0
                ALC(NC,NQ,NR-1) = 0.D0
             END DO
          END IF
       END IF

       ! Diagonal term on handled variable
       BLC(1,LQ,NR) = 1.D0
       NLCR(1,LQ,NR) = LQ

       IF(PRESENT(VAL)) PLC(1,LQ,NR) = PLC(1,LQ,NR) - VAL
    ELSE
       NLCMAX(LQ) = NLCMAX(LQ) + 1
       PLC(NLCMAX(LQ),LQ,NR) = VAL * DTf(LQ)
       NLCR(NLCMAX(LQ),LQ,NR) = 0
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
