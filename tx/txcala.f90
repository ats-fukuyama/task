!     $Id$
module tx_coefficients
  use tx_commons
  use tx_core_module
  implicit none
  private
  real(8), dimension(:,:,:,:), allocatable :: ELM, PELM
  real(8), dimension(:), allocatable ::  &
       & rNuIN0, rNuCXN0, rNubeBE, rNubiBI, rNuTeiEI,&
       & rNuCXN1, rMueNe, rMuiNi, dPNeV, dPNiV, &
       & RUbthV, UethVR, UithVR, EthVR, RUerV, RUirV, UerVR, UirVR, &
       & FWpheBB, FWphiBB, dAphV, FWpheBB2, FWphiBB2, &
       & BphBNi, BthBNi, Dbrpft, &
       & rNuei1EI, rNuei2BthEI, rNuei3EI, &
       & rNube1BE, rNube2BthBE, rNube3BE, &
       & Vbparaph, RVbparath, &
       & Chie1, Chie2, Chii1, Chii2
!!sqeps       &, sqeps_perp, sqeps_perp_inv
!!rp_conv       &, rNubLL
  real(8), dimension(:), allocatable :: UNITY
  real(8) :: DTt, DTf(1:NQM), invDT, BeamSW, RpplSW, &
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

    INTEGER(4) :: NE, NR, NC, NQ, N

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
    !   PELM : Elemental matrix of one term, equation, element   !
    !          for non-variable terms                            !
    !************************************************************!

    allocate(ELM(1:NEMAX,1:4,0:NCM,1:NQMAX))
    allocate(PELM(1:NEMAX,1:4,1:NCM,1:NQMAX))

    ELM(1:NEMAX,1:4,0:NCM,1:NQMAX) = 0.D0
    PELM(1:NEMAX,1:4,1:NCM,1:NQMAX) = 0.D0

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

    !     Coefficients

    N = NRMAX

    allocate(rNuIN0(0:N), rNuCXN0(0:N), rNubeBE(0:N), rNubiBI(0:N), rNuTeiEI(0:N), &
       &     rNuCXN1(0:N), rMueNe(0:N), rMuiNi(0:N), dPNeV(0:N), dPNiV(0:N), &
       &     RUbthV(0:N), UethVR(0:N), UithVR(0:N), &
       &     EthVR(0:N), RUerV(0:N), RUirV(0:N), UerVR(0:N), UirVR(0:N), &
       &     FWpheBB(0:N), FWphiBB(0:N), dAphV(0:N), FWpheBB2(0:N), FWphiBB2(0:N), &
       &     BphBNi(0:N), BthBNi(0:N), Dbrpft(0:N), &
       &     rNuei1EI(0:N), rNuei2BthEI(0:N), rNuei3EI(0:N), &
       &     rNube1BE(0:N), rNube2BthBE(0:N), rNube3BE(0:N), &
       &     Vbparaph(0:N), RVbparath(0:N), &
       &     Chie1(0:N), Chie2(0:N), Chii1(0:N), Chii2(0:N))
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

    !     Ripple trapped beam

    CALL LQr1CC

    !     Elemental equations -> Nodal equations

    DO NE = 1, NEMAX
       NR = NE - 1
       DO NQ = 1, NQMAX
          NC = 0
             BLC(NC,NQ,NR) = BLC(NC,NQ,NR) + ELM(NE,1,NC,NQ) * DTt
             ALC(NC,NQ,NR) = ALC(NC,NQ,NR) + ELM(NE,2,NC,NQ) * DTt
          DO NC = 1, NLCMAX(NQ)
             BLC(NC,NQ,NR) = BLC(NC,NQ,NR) + ELM(NE,1,NC,NQ) * DTf(NQ)
             ALC(NC,NQ,NR) = ALC(NC,NQ,NR) + ELM(NE,2,NC,NQ) * DTf(NQ)
             PLC(NC,NQ,NR) = PLC(NC,NQ,NR) +(PELM(NE,1,NC,NQ) + PELM(NE,2,NC,NQ)) * DTf(NQ)
          END DO
       END DO
       NR = NE
       DO NQ = 1, NQMAX
          NC = 0
             CLC(NC,NQ,NR) = CLC(NC,NQ,NR) + ELM(NE,3,NC,NQ) * DTt
             BLC(NC,NQ,NR) = BLC(NC,NQ,NR) + ELM(NE,4,NC,NQ) * DTt
          DO NC = 1, NLCMAX(NQ)
             CLC(NC,NQ,NR) = CLC(NC,NQ,NR) + ELM(NE,3,NC,NQ) * DTf(NQ)
             BLC(NC,NQ,NR) = BLC(NC,NQ,NR) + ELM(NE,4,NC,NQ) * DTf(NQ)
             PLC(NC,NQ,NR) = PLC(NC,NQ,NR) +(PELM(NE,3,NC,NQ) + PELM(NE,4,NC,NQ)) * DTf(NQ)
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
    ! When ripple effect is on (FSRP /= 0), ripple diffusion term will be activated.
    ! Then we must impose two boundary conditions at each equation.
    IF(FSRP /= 0.D0) THEN
       CALL BOUNDARY(0    ,LQb3,0)
       CALL BOUNDARY(NRMAX,LQb3,0)
       CALL BOUNDARY(NRMAX,LQb4,0)
!       CALL BOUNDARY(0    ,LQr1,0)
    END IF

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
    CALL BOUNDARY(NRMAX,LQn1,1, 2.D0*R(NRMAX)*rGASPF)

    deallocate(ELM,PELM)
    deallocate(rNuIN0, rNuCXN0, rNubeBE, rNubiBI, rNuTeiEI,&
       &       rNuCXN1, rMueNe, rMuiNi, dPNeV, dPNiV, &
       &       RUbthV, UethVR, UithVR, EthVR, RUerV, RUirV, UerVR, UirVR, &
       &       FWpheBB, FWphiBB, dAphV, FWpheBB2, FWphiBB2, &
       &       BphBNi, BthBNi, Dbrpft, &
       &       rNuei1EI, rNuei2BthEI, rNuei3EI, &
       &       rNube1BE, rNube2BthBE, rNube3BE, &
       &       Vbparaph, RVbparath, &
       &       Chie1, Chie2, Chii1, Chii2)
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

    use tx_interface, only : DERIVS
    INTEGER(4) :: NR
    REAL(8) :: BBL

    rNuIN0(0:NRMAX)   = rNuION(0:NRMAX) * PNeV(0:NRMAX) &
         &            / (PN01V(0:NRMAX) + PN02V(0:NRMAX))
    rNuCXN0(0:NRMAX)  = rNuiCX(0:NRMAX) * PNiV(0:NRMAX) &
         &            / (PN01V(0:NRMAX) + PN02V(0:NRMAX))
    rNuCXN1(0:NRMAX)  = rNuiCX(0:NRMAX) * PNiV(0:NRMAX) * PT01V(0:NRMAX)
    rNuTeiEI(0:NRMAX) = rNuTei(0:NRMAX) * PNeV(0:NRMAX) / PNiV(0:NRMAX)
    rNubeBE(0:NRMAX)  = rNube(0:NRMAX)  * PNbV(0:NRMAX) / PNeV(0:NRMAX)
    rNubiBI(0:NRMAX)  = rNubi(0:NRMAX)  * PNbV(0:NRMAX) / PNiV(0:NRMAX)
    rMueNe(0:NRMAX)   = rMue(0:NRMAX)   / PNeV(0:NRMAX)
    rMuiNi(0:NRMAX)   = rMui(0:NRMAX)   / PNiV(0:NRMAX)
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
    CALL DERIVS(PSI,X,LQe1,NQMAX,NRMAX,dPNeV)
    CALL DERIVS(PSI,X,LQi1,NQMAX,NRMAX,dPNiV)
    CALL DERIVS(PSI,X,LQm4,NQMAX,NRMAX,dAphV)
    FWpheBB(0:NRMAX)  =- 2.D0 * (dAphV(0:NRMAX) * rMUb2) * FWthphe(0:NRMAX)
    FWphiBB(0:NRMAX)  =- 2.D0 * (dAphV(0:NRMAX) * rMUb2) * FWthphi(0:NRMAX)
    FWpheBB2(0:NRMAX) =(BthV(0:NRMAX) / BphV(0:NRMAX))**2 * FWthe(0:NRMAX)
    FWphiBB2(0:NRMAX) =(BthV(0:NRMAX) / BphV(0:NRMAX))**2 * FWthi(0:NRMAX)

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

!!sqeps    sqeps_perp(0:NRMAX) = SQRT(PNiV(0:NRMAX)*1.D20*AMI/(BphV(0:NRMAX)**2+BthV(0:NRMAX)**2))
!!sqeps    sqeps_perp_inv(0:NRMAX) = 1.D0 / sqeps_perp(0:NRMAX)

  END SUBROUTINE LQCOEF

!***************************************************************
!
!   Poisson Equation: phi
!
!**************************************************************

  SUBROUTINE LQm1CC

    ! phi'(0) : 0

    ELM(1:NEMAX,1:4,1,LQm1) =   4.D0 * sqeps0 / (AEE * 1.D20) * fem_int(18,UNITY)
    NLC(1,LQm1) = LQm1

    ELM(1:NEMAX,1:4,2,LQm1) =        fem_int(1) / sqeps0
    NLC(2,LQm1) = LQe1

    ELM(1:NEMAX,1:4,3,LQm1) = - PZ * fem_int(1) / sqeps0
    NLC(3,LQm1) = LQi1

    ELM(1:NEMAX,1:4,4,LQm1) = - PZ * fem_int(1) / sqeps0 * BeamSW
    NLC(4,LQm1) = LQb1

    ELM(1:NEMAX,1:4,5,LQm1) = - PZ * fem_int(2,rip_rat) / sqeps0 * BeamSW * RpplSW
    NLC(5,LQm1) = LQr1
!!sqeps    ELM(1:NEMAX,1:4,1,LQm1) =   4.D0 / (AEE * 1.D20) * fem_int(18,sqeps_perp)
!!sqeps    NLC(1,LQm1) = LQm1
!!sqeps
!!sqeps    ELM(1:NEMAX,1:4,2,LQm1) =        fem_int(2,sqeps_perp_inv)
!!sqeps    NLC(2,LQm1) = LQe1
!!sqeps
!!sqeps    ELM(1:NEMAX,1:4,3,LQm1) = - PZ * fem_int(2,sqeps_perp_inv)
!!sqeps    NLC(3,LQm1) = LQi1
!!sqeps
!!sqeps    ELM(1:NEMAX,1:4,4,LQm1) = - PZ * fem_int(2,sqeps_perp_inv) * BeamSW
!!sqeps    NLC(4,LQm1) = LQb1
!!sqeps
!!sqeps    ELM(1:NEMAX,1:4,5,LQm1) = - PZ * fem_int(28,rip_rat,sqeps_perp_inv) * BeamSW * RpplSW
!!sqeps    NLC(5,LQm1) = LQr1

    ! phi(b) : 0

    NLCMAX(LQm1) = 5
    RETURN
  END SUBROUTINE LQm1CC

!***************************************************************
!
!   Ampere's Law: Atheta'
!
!***************************************************************

  SUBROUTINE LQm2CC

    ! (r*Atheta)'(0) : 0

    ELM(1:NEMAX,1:4,0,LQm2) = fem_int(1) * EPS0 * invDT
    NLC(0,LQm2) = LQm2

    ! rot Bphi

    ELM(1:NEMAX,1:4,1,LQm2) = (- 4.D0 * fem_int(18,UNITY) - 4.D0 * fem_int(4)) * AMPm5
    NLC(1,LQm2) = LQm5

    ! Electron current

    ELM(1:NEMAX,1:4,2,LQm2) = -      AEE * 1.D20 * fem_int(1)
    NLC(2,LQm2) = LQe3

    ! Ion current

    ELM(1:NEMAX,1:4,3,LQm2) =   PZ * AEE * 1.D20 * fem_int(1)
    NLC(3,LQm2) = LQi3

    ! Beam ion current

    ELM(1:NEMAX,1:4,4,LQm2) =   PZ * AEE * 1.D20 * fem_int(1) * BeamSW
    NLC(4,LQm2) = LQb3

    ! (1/r)*(r*Atheta)'(NRMAX) : BB

    NLCMAX(LQm2) = 4
    RETURN
  END SUBROUTINE LQm2CC

!***************************************************************
!
!   Ampere's Law: Aphi'
!
!***************************************************************

  SUBROUTINE LQm3CC

    ! Aphi'(0) : 0

    ELM(1:NEMAX,1:4,0,LQm3) = fem_int(1) * EPS0 * rMUb1 * invDT
    NLC(0,LQm3) = LQm3

    ! rot Btheta

    ELM(1:NEMAX,1:4,1,LQm3) = - 4.D0 * fem_int(18,UNITY)
    NLC(1,LQm3) = LQm4

    ! Electron current

    ELM(1:NEMAX,1:4,2,LQm3) = - rMUb1      * AEE * 1.D20 * fem_int(1) * AMPe4
    NLC(2,LQm3) = LQe4

    ! Ion current

    ELM(1:NEMAX,1:4,3,LQm3) =   rMUb1 * PZ * AEE * 1.D20 * fem_int(1)
    NLC(3,LQm3) = LQi4

    ! Beam ion current

    ELM(1:NEMAX,1:4,4,LQm3) =   rMUb1 * PZ * AEE * 1.D20 * fem_int(1) * BeamSW
    NLC(4,LQm3) = LQb4

!   ! Virtual current for helical system
!
    PELM(1:NEMAX,1:4,5,LQm3) =  rMUb1 * fem_int(-1,AJV)
    NLC(5,LQm3) = 0

    ! Aphi'(NRMAX) : -Bthb

    NLCMAX(LQm3) = 5
    RETURN
  END SUBROUTINE LQm3CC

!**************************************************************
!
!   Faraday's Law : Aphi
!
!***************************************************************

  SUBROUTINE LQm4CC

    ELM(1:NEMAX,1:4,0,LQm4) = fem_int(1) * invDT
    NLC(0,LQm4) = LQm4

    ! Aphi'

    ELM(1:NEMAX,1:4,1,LQm4) = fem_int(1) / rMUb2
    NLC(1,LQm4) = LQm3

    NLCMAX(LQm4) = 1
    RETURN
  END SUBROUTINE LQm4CC

!***************************************************************
!
!   Faraday's Law : Atheta
!
!***************************************************************

  SUBROUTINE LQm5CC

    ELM(1:NEMAX,1:4,0,LQm5) = fem_int(1) * invDT * AMPm5
    NLC(0,LQm5) = LQm5

    ! r * Atheta'

    ELM(1:NEMAX,1:4,1,LQm5) = fem_int(1) / rMU0
    NLC(1,LQm5) = LQm2

    NLCMAX(LQm5) = 1
    RETURN
  END SUBROUTINE LQm5CC

!***************************************************************
!
!   Electron Density Equation
!
!***************************************************************

  SUBROUTINE LQe1CC

    ELM(1:NEMAX,1:4,0,LQe1) = fem_int(1) * invDT
    NLC(0,LQe1) = LQe1

    ! Divergence

    ELM(1:NEMAX,1:4,1,LQe1) = - 2.D0 * fem_int(4)
    NLC(1,LQe1) = LQe2

    ! Ionization of n01 and n02

    ELM(1:NEMAX,1:4,2,LQe1) =   fem_int(2,rNuIN0)
    NLC(2,LQe1) = LQn1

    ELM(1:NEMAX,1:4,3,LQe1) =   fem_int(2,rNuIN0)
    NLC(3,LQe1) = LQn2

    ! Loss to divertor

    ELM(1:NEMAX,1:4,4,LQe1) = - fem_int(2,rNuL)
    NLC(4,LQe1) = LQe1

    PELM(1:NEMAX,1:4,5,LQe1) =  PNeDIV * fem_int(-1,rNuL)
    NLC(5,LQe1) = 0

    ! Generated by NBI (Ionization)

    PELM(1:NEMAX,1:4,6,LQe1) =   (1.D0 - RatCX) * fem_int(-1,SNBe)
    NLC(6,LQe1) = 0

    !  Diffusion of electrons (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,7,LQe1) = - 4.D0 * fem_int(18,DMAGe)
    NLC(7,LQe1) = LQe1

    NLCMAX(LQe1) = 7
    RETURN
  END SUBROUTINE LQe1CC

!***************************************************************
!
!   Electron Radial Flow (SUPG)
!
!***************************************************************

  SUBROUTINE LQe2CC

    ! Ns*Usr(0) : fixed

    ELM(1:NEMAX,1:4,0,LQe2) = fem_int(1) * invDT &
         &                  + fem_int(8) * fem_int(0) * invDT
    NLC(0,LQe2) = LQe2

    ! Advection
    
    ELM(1:NEMAX,1:4,1,LQe2) = - 2.D0 * fem_int( 3,RUerV) + fem_int(2,UerVR) &
         &                  +(- 2.D0 * fem_int(10,RUerV) + fem_int(9,UerVR)) * fem_int(0)
    ELM(1      ,1:4,1,LQe2) = - 2.D0 * fem_int_point( 3,0,RUerV) + fem_int_point(2,0,UerVR) &
         &                  +(- 2.D0 * fem_int_point(10,0,RUerV) + fem_int_point(9,0,UerVR)) &
         &                  * fem_int_point(0,1)
    NLC(1,LQe2) = LQe2

    ! Centrifugal force

    ELM(1:NEMAX,1:4,2,LQe2) = fem_int(2,UethVR) &
         &                  + fem_int(9,UethVR) * fem_int(0)
    ELM(1      ,1:4,2,LQe2) = fem_int_point(2,0,UethVR) &
         &                  + fem_int_point(9,0,UethVR) * fem_int_point(0,1)
    NLC(2,LQe2) = LQe3

    ! Pressure gradient force

    ELM(1:NEMAX,1:4,3,LQe2) = - 2.D0 * rKeV / AME * fem_int(17,UNITY) &
         &                    - 2.D0 * rKeV / AME * fem_int(18,UNITY) * fem_int(0)
    NLC(3,LQe2) = LQe5

    ! Radial E force

    ELM(1:NEMAX,1:4,4,LQe2) =   2.D0 * (AEE / AME) * fem_int(17,PNeV) &
         &                  +   2.D0 * (AEE / AME) * fem_int(18,PNeV) * fem_int(0)
    NLC(4,LQe2) = LQm1

    ! v x B force

    ELM(1:NEMAX,1:4,5,LQe2) = -        (AEE / AME) * fem_int( 2,BphV) &
         &                    -        (AEE / AME) * fem_int( 9,BphV) * fem_int(0)
    NLC(5,LQe2) = LQe3

    ELM(1:NEMAX,1:4,6,LQe2) =(- 2.D0 * (AEE / AME) * fem_int(16,AphV) &
         &                    - 2.D0 * (AEE / AME) * fem_int(20,UNITY,AphV) * fem_int(0)) * AMPe4
    NLC(6,LQe2) = LQe4

    NLCMAX(LQe2) = 6

    IF(MDFIXT /= 0) THEN
       ELM(1:NEMAX,1:4,3,LQe2) = - 2.D0 * rKeV / AME &
            &                  * (  fem_int(17,PTeV) + fem_int(18,PTeV) * fem_int(0) &
            &                     + fem_int(16,PTeV) + fem_int(34,PSI,PTeV) * fem_int(0))
       NLC(3,LQe2) = LQe1
    END IF

    RETURN
  END SUBROUTINE LQe2CC

!***************************************************************
!
!   Electron Poloidal Flow
!
!***************************************************************

  SUBROUTINE LQe3CC

    integer(4) :: N

    ! Ns*UsTheta(0) : 0

    ELM(1:NEMAX,1:4, 0,LQe3) = fem_int(1) * invDT
    NLC( 0,LQe3) = LQe3

    ! Advection

    ELM(1:NEMAX,1:4, 1,LQe3) = - 2.D0 * fem_int(3,RUerV)
    NLC( 1,LQe3) = LQe3

    ! Viscosity force

!    ELM(1:NEMAX,1:4, 2,LQe3) = - 4.D0 * fem_int(18,rMue) &
!         &                     - 4.D0 * fem_int(39,rMueNe,dPNeV) &
!         &                     - 4.D0 * fem_int( 3,rMue)
    ELM(1:NEMAX,1:4, 2,LQe3) = - 4.D0 * fem_int(18,rMue) &
         &                     - 4.D0 * fem_int( 3,rMue)
    NLC( 2,LQe3) = LQe3

    ELM(1:NEMAX,1:4, 3,LQe3) =   4.D0 * fem_int(41,rMue,RUethV)
    NLC( 3,LQe3) = LQe1

    ! Poloidal E force

    ELM(1:NEMAX,1:4, 4,LQe3) =   (AEE / AME) * fem_int(2,PNeV)
    NLC( 4,LQe3) = LQm2

    ! v x B force

    ELM(1:NEMAX,1:4, 5,LQe3) =   (AEE / AME) * fem_int(2,BphV)
    NLC( 5,LQe3) = LQe2

    ! Neoclassical viscosity force

    ELM(1:NEMAX,1:4, 6,LQe3) = - fem_int(2,rNueNC)
    NLC( 6,LQe3) = LQe3

    ! Collisional friction force with ions

    ELM(1:NEMAX,1:4, 7,LQe3) = - fem_int(2,rNuei1)
    NLC( 7,LQe3) = LQe3

    ELM(1:NEMAX,1:4, 8,LQe3) =   fem_int(2,rNuei1EI)
    NLC( 8,LQe3) = LQi3

    ELM(1:NEMAX,1:4, 9,LQe3) =   2.D0 * fem_int(37,rNuei2Bth,AphV) * AMPe4
    NLC( 9,LQe3) = LQe4

    ELM(1:NEMAX,1:4,10,LQe3) = - 2.D0 * fem_int(37,rNuei2BthEI,AphV)
    NLC(10,LQe3) = LQi4

    ! Collisional friction with beam ions

    ELM(1:NEMAX,1:4,11,LQe3) = - (AMB / AME) * fem_int(2,rNube1BE)
    NLC(11,LQe3) = LQe3

    ELM(1:NEMAX,1:4,12,LQe3) =   (AMB / AME) * fem_int(2,rNube1)
    NLC(12,LQe3) = LQb3

    ELM(1:NEMAX,1:4,13,LQe3) =   2.D0 * (AMB / AME) * fem_int(37,rNube2BthBE,AphV) * AMPe4
    NLC(13,LQe3) = LQe4

    ELM(1:NEMAX,1:4,14,LQe3) = - 2.D0 * (AMB / AME) * fem_int(37,rNube2Bth,AphV)
    NLC(14,LQe3) = LQb4

    IF(MDLWTB == 0) THEN

       ! Wave interaction force (electron driven)

       ELM(1:NEMAX,1:4,15,LQe3) = - 1.D0 / AME * fem_int(2,FWthe)
       NLC(15,LQe3) = LQe3

       ELM(1:NEMAX,1:4,16,LQe3) =   1.D0 / AME * fem_int(2,FWthe)
       NLC(16,LQe3) = LQi3

       ELM(1:NEMAX,1:4,17,LQe3) = - 2.D0 / AME * fem_int(36,AphV,FWthphe) * AMPe4
       NLC(17,LQe3) = LQe4

       ELM(1:NEMAX,1:4,18,LQe3) =   2.D0 / AME * fem_int(36,AphV,FWthphe)
       NLC(18,LQe3) = LQi4

       ELM(1:NEMAX,1:4,19,LQe3) =   1.D0 / AME * fem_int(44,FWthe,WPM)
       NLC(19,LQe3) = LQe1

!!ion       ! Wave interaction force (ion driven)
!!ion
!!ion       ELM(1:NEMAX,1:4,15,LQe3) =   1.D0 / AME * fem_int(2,FWthi)
!!ion       NLC(15,LQe3) = LQi3
!!ion
!!ion       ELM(1:NEMAX,1:4,16,LQe3) = - 1.D0 / AME * fem_int(2,FWthi)
!!ion       NLC(16,LQe3) = LQe3
!!ion
!!ion       ELM(1:NEMAX,1:4,17,LQe3) =   2.D0 / AME * fem_int(36,AphV,FWthphi)
!!ion       NLC(17,LQe3) = LQi4
!!ion
!!ion       ELM(1:NEMAX,1:4,18,LQe3) = - 2.D0 / AME * fem_int(36,AphV,FWthphi) * AMPe4
!!ion       NLC(18,LQe3) = LQe4
!!ion
!!ion       ELM(1:NEMAX,1:4,19,LQe3) = - 1.D0 / AME * fem_int(44,FWthi,WPM)
!!ion       NLC(19,LQe3) = LQi1

       N = 0
    ELSEIF(MDLWTB == 1) THEN

       ! Wave interaction force (electron driven)

       ELM(1:NEMAX,1:4,15,LQe3) = - 2.D0 / AME * fem_int(17,WNthe)
       NLC(15,LQe3) = LQe1

!      ELM(1:NEMAX,1:4,16,LQe3) = - 2.D0 / AME * fem_int(36,PhiV,WEMthe)
!      NLC(16,LQe3) = LQe1
!
!      ELM(1:NEMAX,1:4,17,LQe3) = - 2.D0 / AME * fem_int(15,WWthe)
!      NLC(17,LQe3) = LQe1
!
!      ELM(1:NEMAX,1:4,18,LQe3) = - 2.D0 / AME * rKeV * fem_int(17,WT1the)
!      NLC(18,LQe3) = LQe5
!
!      ELM(1:NEMAX,1:4,19,LQe3) =   2.D0 / AME * fem_int(17,WT2the)
!      NLC(19,LQe3) = LQe1

       ! Wave interaction force (ion driven)

       ELM(1:NEMAX,1:4,20,LQe3) =   2.D0 / AME * fem_int(17,WNthi)
       NLC(20,LQe3) = LQi1

!      ELM(1:NEMAX,1:4,21,LQe3) =   2.D0 / AME * fem_int(36,PhiV,WEMthi)
!      NLC(21,LQe3) = LQi1
!
!      ELM(1:NEMAX,1:4,22,LQe3) =   1.D0 / AME * fem_int(15,WWthi)
!      NLC(22,LQe3) = LQi1
!
!      ELM(1:NEMAX,1:4,23,LQe3) =   2.D0 / AME * rKeV * fem_int(17,WT1thi)
!      NLC(23,LQe3) = LQi5
!
!      ELM(1:NEMAX,1:4,24,LQe3) = - 2.D0 / AME * fem_int(17,WT2thi)
!      NLC(24,LQe3) = LQi1

       N = 4
    END IF

    ! Loss to divertor

    ELM(1:NEMAX,1:4,20+N,LQe3) = - 2.D0 * fem_int(2,rNuL) * fact
    NLC(20+N,LQe3) = LQe3

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,21+N,LQe3) = - fem_int(2,rNu0e)
    NLC(21+N,LQe3) = LQe3

    ! Helical neoclassical viscosity force (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,22+N,LQe3) = - UHth * UHth * fem_int(2,rNueHL)
    NLC(22+N,LQe3) = LQe3

    ELM(1:NEMAX,1:4,23+N,LQe3) = UHph * UHth * fem_int(2,rNueHL) * AMPe4
    NLC(23+N,LQe3) = LQe4

    !  Diffusion of electrons (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,24+N,LQe3) = - 4.D0 * fem_int(18,DMAGe)
    NLC(24+N,LQe3) = LQe3

    ! Ns*UsTheta(NRMAX) : 0

    NLCMAX(LQe3) = 24+N
    RETURN
  END SUBROUTINE LQe3CC

!***************************************************************
!
!   Electron Toroidal Flow
!
!***************************************************************

  SUBROUTINE LQe4CC

    ! Uephi(0)' : 0

    ELM(1:NEMAX,1:4, 0,LQe4) = fem_int(1) * invDT * AMPe4
    NLC( 0,LQe4) = LQe4

    ! Advection

    ELM(1:NEMAX,1:4, 1,LQe4) = - 2.D0 * fem_int(3,RUerV) * AMPe4
    NLC( 1,LQe4) = LQe4
    
    ! Viscosity force
    
!    ELM(1:NEMAX,1:4, 2,LQe4) = - 4.D0 *(fem_int(18,rMue) + fem_int(39,rMueNe,dPNeV)) * AMPe4
    ELM(1:NEMAX,1:4, 2,LQe4) = - 4.D0 * fem_int(18,rMue) * AMPe4
    NLC( 2,LQe4) = LQe4

    ELM(1:NEMAX,1:4, 3,LQe4) =   4.D0 * fem_int(41,rMue,UephV)
    NLC( 3,LQe4) = LQe1

    ! Toroidal E force

    ELM(1:NEMAX,1:4, 4,LQe4) =   (AEE / AME) * fem_int(2,PNeV)
    NLC( 4,LQe4) = LQm3

    ! v x B force

    ELM(1:NEMAX,1:4, 5,LQe4) =   2.D0 * (AEE / AME) * fem_int(6,AphV)
    NLC( 5,LQe4) = LQe2

    ! Collisional friction with bulk ions

    ELM(1:NEMAX,1:4, 6,LQe4) = - fem_int(2,rNuei3) * AMPe4
    NLC( 6,LQe4) = LQe4

    ELM(1:NEMAX,1:4, 7,LQe4) =   fem_int(2,rNuei3EI)
    NLC( 7,LQe4) = LQi4

    ELM(1:NEMAX,1:4, 8,LQe4) =   2.D0 * fem_int(29,rNuei2Bth,AphV)
    NLC( 8,LQe4) = LQe3

    ELM(1:NEMAX,1:4, 9,LQe4) = - 2.D0 * fem_int(29,rNuei2BthEI,AphV)
    NLC( 9,LQe4) = LQi3

    ! Collisional friction with beam ions

    ELM(1:NEMAX,1:4,10,LQe4) = - (AMB / AME) * fem_int(2,rNube3BE) * AMPe4
    NLC(10,LQe4) = LQe4

    ELM(1:NEMAX,1:4,11,LQe4) =   (AMB / AME) * fem_int(2,rNube3)
    NLC(11,LQe4) = LQb4

    ELM(1:NEMAX,1:4,12,LQe4) =   2.D0 * (AMB / AME) * fem_int(29,rNube2BthBE,AphV)
    NLC(12,LQe4) = LQe3

    ELM(1:NEMAX,1:4,13,LQe4) = - 2.D0 * (AMB / AME) * fem_int(29,rNube2Bth,AphV)
    NLC(13,LQe4) = LQb3

    IF(MDLWTB == 0) THEN

       ! Wave interaction force (electron driven)

       ELM(1:NEMAX,1:4,14,LQe4) =   1.D0 / AME * fem_int(2,FWpheBB)
       NLC(14,LQe4) = LQe3

       ELM(1:NEMAX,1:4,15,LQe4) = - 1.D0 / AME * fem_int(2,FWpheBB)
       NLC(15,LQe4) = LQi3

       ELM(1:NEMAX,1:4,16,LQe4) = - 1.D0 / AME * fem_int(2,FWpheBB2) * AMPe4
       NLC(16,LQe4) = LQe4

       ELM(1:NEMAX,1:4,17,LQe4) =   1.D0 / AME * fem_int(2,FWpheBB2)
       NLC(17,LQe4) = LQi4

       ELM(1:NEMAX,1:4,18,LQe4) = - 1.D0 / AME * fem_int(44,FWpheBB,WPM)
       NLC(18,LQe4) = LQe1

!!ion       ! Wave interaction force (ion driven)
!!ion
!!ion       ELM(1:NEMAX,1:4,14,LQe4) = - 1.D0 / AME * fem_int(2,FWphiBB)
!!ion       NLC(14,LQe4) = LQi3
!!ion
!!ion       ELM(1:NEMAX,1:4,15,LQe4) =   1.D0 / AME * fem_int(2,FWphiBB)
!!ion       NLC(15,LQe4) = LQe3
!!ion
!!ion       ELM(1:NEMAX,1:4,16,LQe4) =   1.D0 / AME * fem_int(2,FWphiBB2)
!!ion       NLC(16,LQe4) = LQi4
!!ion
!!ion       ELM(1:NEMAX,1:4,17,LQe4) = - 1.D0 / AME * fem_int(2,FWphiBB2) * AMPe4
!!ion       NLC(17,LQe4) = LQe4
!!ion
!!ion       ELM(1:NEMAX,1:4,18,LQe4) =   1.D0 / AME * fem_int(44,FWphiBB,WPM)
!!ion       NLC(18,LQe4) = LQi1

    END IF

    ! Loss to divertor

    ELM(1:NEMAX,1:4,19,LQe4) = - 2.D0 * fem_int(2,rNuL) * AMPe4 * fact
    NLC(19,LQe4) = LQe4

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,20,LQe4) = - fem_int(2,rNu0e) * AMPe4
    NLC(20,LQe4) = LQe4

    ! Helical neoclassical viscosity force (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,21,LQe4) =  UHth * UHph * fem_int(2,rNueHL)
    NLC(21,LQe4) = LQe3

    ELM(1:NEMAX,1:4,22,LQe4) = - UHph * UHph * fem_int(2,rNueHL) * AMPe4
    NLC(22,LQe4) = LQe4

    !  Diffusion of electrons (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,23,LQe4) = - 4.D0 * fem_int(18,DMAGe) * AMPe4
    NLC(23,LQe4) = LQe4

    NLCMAX(LQe4) = 23
    RETURN
  END SUBROUTINE LQe4CC

!***************************************************************
!
!   Electron Energy Transport: Te
!
!***************************************************************

  SUBROUTINE LQe5CC

    ! Temperature evolution
    
    IF(MDFIXT == 0) THEN
       ELM(1:NEMAX,1:4, 0,LQe5) =   1.5D0 * fem_int(1) * invDT
       NLC( 0,LQe5) = LQe5

       ! Advection

       ELM(1:NEMAX,1:4, 1,LQe5) = - 5.D0 * fem_int(3,RUerV)
       NLC( 1,LQe5) = LQe5

       ! Conduction transport

       ELM(1:NEMAX,1:4, 2,LQe5) = - 4.D0 * fem_int(18,Chie1)
       NLC( 2,LQe5) = LQe5

       ELM(1:NEMAX,1:4, 3,LQe5) =   4.D0 * fem_int(41,Chie2,PTeV)
       NLC( 3,LQe5) = LQe1

       ! Redundant heat advection

       ELM(1:NEMAX,1:4, 4,LQe5) = 2.D0 * fem_int(5,RUerV)
       NLC( 4,LQe5) = LQe5

       ! Joule heating

       ELM(1:NEMAX,1:4, 5,LQe5) = - AEE / rKeV * fem_int(2,EthVR)
       ELM(1      ,1:4, 5,LQe5) = - AEE / rKeV * fem_int_point(2,0,EthVR)
       NLC( 5,LQe5) = LQe3

       ELM(1:NEMAX,1:4, 6,LQe5) = - AEE / rKeV * fem_int(2,EphV) * AMPe4
       NLC( 6,LQe5) = LQe4

       ! Collisional transfer with ions (Energy equilibration)

       ELM(1:NEMAX,1:4, 7,LQe5) = - 1.5d0 * fem_int(2,rNuTei)
       NLC( 7,LQe5) = LQe5

       ELM(1:NEMAX,1:4, 8,LQe5) =   1.5d0 * fem_int(2,rNuTeiEI)
       NLC( 8,LQe5) = LQi5

       ! Loss to diverter

       ELM(1:NEMAX,1:4, 9,LQe5) = -                  fem_int( 2,rNuL)
       NLC( 9,LQe5) = LQe5

       PELM(1:NEMAX,1:4,10,LQe5) =          PNeDIV * fem_int(-2,rNuL,PTeV)
       NLC(10,LQe5) = 0

       ELM(1:NEMAX,1:4,11,LQe5) = - 1.5D0          * fem_int( 2,rNuLTe)
       NLC(11,LQe5) = LQe5

       PELM(1:NEMAX,1:4,12,LQe5) =  1.5D0 * PTeDIV * fem_int(-2,rNuLTe,PNeV)
       NLC(12,LQe5) = 0

       ! Direct heating (RF)

       PELM(1:NEMAX,1:4,13,LQe5) =   1.D0 / (1.D20 * rKeV) * fem_int(-1,PRFe)
       NLC(13,LQe5) = 0

       ! Radiation loss  (Bremsstrahlung)

       PELM(1:NEMAX,1:4,14,LQe5) = - 1.D0 / (1.D20 * rKeV) * fem_int(-1,PBr)
       NLC(14,LQe5) = 0

       ! Ionization loss of n01 and n02

       ELM(1:NEMAX,1:4,15,LQe5) = - (EION * 1.D-3) * fem_int(2,rNuIN0)
       NLC(15,LQe5) = LQn1

       ELM(1:NEMAX,1:4,16,LQe5) = - (EION * 1.D-3) * fem_int(2,rNuIN0)
       NLC(16,LQe5) = LQn2

       ! Collisional NBI heating (Perp + Tan)

       PELM(1:NEMAX,1:4,17,LQe5) = Eb * fem_int(-2,SNB,PNBcol_e)
       NLC(17,LQe5) = 0

       ! Collisional heating with beam

!!oldNBI      ELM(1:NEMAX,1:4,18,LQe5) = - 0.5D0 * AMB * (PNBCD * Vb) / rKeV &
!!oldNBI           &                * fem_int(2,rNubeBE) * AMPe4
!!oldNBI      NLC(18,LQe5) = LQe4
!!oldNBI
!!oldNBI      ELM(1:NEMAX,1:4,19,LQe5) =   0.5D0 * AMB * (PNBCD * Vb) / rKeV &
!!oldNBI           &                * fem_int(2,rNube)
!!oldNBI      NLC(19,LQe5) = LQb4

       ! Simpified Alpha heating

       PELM(1:NEMAX,1:4,18,LQe5) = 1.D0 / (1.D20 * rKeV) * fem_int(-1,PALFe)
       NLC(18,LQe5) = 0
       
       !  Diffusion of electrons (***AF 2008-06-08)

       ELM(1:NEMAX,1:4,19,LQe5) = - 4.D0 * fem_int(18,DMAGe)
       NLC(19,LQe5) = LQe5

       NLCMAX(LQe5) = 19

    ELSE

       !  Fixed temperature profile

       ELM(1:NEMAX,1:4,0,LQe5) = fem_int(1) * invDT
       NLC(0,LQe5) = LQe5

       NLCMAX(LQe5) = 0
    END IF

    RETURN
  END SUBROUTINE LQe5CC

!***************************************************************
!
!   Ion Density Equation
!
!***************************************************************

  SUBROUTINE LQi1CC

    ELM(1:NEMAX,1:4,0,LQi1) = fem_int(1) * invDT
    NLC(0,LQi1) = LQi1

    ! Divergence

    ELM(1:NEMAX,1:4,1,LQi1) = - 2.D0 * fem_int(4)
    NLC(1,LQi1) = LQi2

    ! Ionization of n01 and n02

    ELM(1:NEMAX,1:4,2,LQi1) =     1.D0 / PZ * fem_int(2,rNuIN0)
    NLC(2,LQi1) = LQn1

    ELM(1:NEMAX,1:4,3,LQi1) =     1.D0 / PZ * fem_int(2,rNuIN0)
    NLC(3,LQi1) = LQn2

    ! Loss to divertor

    ELM(1:NEMAX,1:4,4,LQi1) = -   1.D0 / PZ * fem_int(2,rNuL)
    NLC(4,LQi1) = LQe1

    PELM(1:NEMAX,1:4,5,LQi1) =  PNeDIV / PZ * fem_int(-1,rNuL)
    NLC(5,LQi1) = 0

    ! Particle source from beam ion

    ELM(1:NEMAX,1:4,6,LQi1) =   fem_int(2,rNuB)
    NLC(6,LQi1) = LQb1

    ! Particle source from ripple trapped beam ions

    ELM(1:NEMAX,1:4,7,LQi1) =   fem_int(28,rNuB,rip_rat) * RpplSW
    NLC(7,LQi1) = LQr1

    ! NBI kick up ions (Charge exchange)

    PELM(1:NEMAX,1:4,8,LQi1) = - RatCX * fem_int(-1,SNBi)
    NLC(8,LQi1) = 0

    ! Parallel loss reduction due to the potential
    ! induced by the parallel loss of the beam ions

    ELM(1:NEMAX,1:4,9,LQi1) =   fem_int(2,rNuLB) * BeamSW
    NLC(9,LQi1) = LQb1

    ! Loss cone loss

    PELM(1:NEMAX,1:4,10,LQi1) =   fem_int(-1,SiLC)
    NLC(10,LQi1) = 0
 
    ! Ion orbit loss

    ELM(1:NEMAX,1:4,11,LQi1) = - fem_int(2,rNuOL)
    NLC(11,LQi1) = LQi1

    !  Diffusion of ions (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,12,LQi1) = - 4.D0 * fem_int(18,DMAGi)
    NLC(12,LQi1) = LQi1

    NLCMAX(LQi1) = 12
    RETURN
  END SUBROUTINE LQi1CC

!***************************************************************
!
!   Ion Radial Flow (SUPG)
!
!***************************************************************
  
  SUBROUTINE LQi2CC

    ! Ns*Usr(0) : fixed

    ELM(1:NEMAX,1:4,0,LQi2) = fem_int(1) * invDT &
         &                  + fem_int(8) * fem_int(0) * invDT
    NLC(0,LQi2) = LQi2

    ! Advection

    ELM(1:NEMAX,1:4,1,LQi2) = - 2.D0 * fem_int( 3,RUirV) + fem_int(2,UirVR) &
         &                  +(- 2.D0 * fem_int(10,RUirV) + fem_int(9,UirVR)) * fem_int(0)
    ELM(1      ,1:4,1,LQi2) = - 2.D0 * fem_int_point( 3,0,RUirV) + fem_int_point(2,0,UirVR) &
         &                  +(- 2.D0 * fem_int_point(10,0,RUirV) + fem_int_point(9,0,UirVR)) &
         &                  * fem_int_point(0,1)
    NLC(1,LQi2) = LQi2

    ! Centrifugal force

    ELM(1:NEMAX,1:4,2,LQi2) = fem_int(2,UithVR) &
         &                  + fem_int(9,UithVR) * fem_int(0)
    ELM(1      ,1:4,2,LQi2) = fem_int_point(2,0,UithVR) &
         &                  + fem_int_point(9,0,UithVR) * fem_int_point(0,1)
    NLC(2,LQi2) = LQi3

    ! Pressure gradient force

    ELM(1:NEMAX,1:4,3,LQi2) = - 2.D0 * rKeV / AMI * fem_int(17,UNITY) &
         &                    - 2.D0 * rKeV / AMI * fem_int(18,UNITY) * fem_int(0)
    NLC(3,LQi2) = LQi5

    ! Radial E force

    ELM(1:NEMAX,1:4,4,LQi2) = - 2.D0 * (PZ * AEE / AMI) * fem_int(17,PNiV) &
         &                    - 2.D0 * (PZ * AEE / AMI) * fem_int(18,PNiV) * fem_int(0)
    NLC(4,LQi2) = LQm1

    ! v x B force

    ELM(1:NEMAX,1:4,5,LQi2) =          (PZ * AEE / AMI) * fem_int( 2,BphV) &
         &                  +          (PZ * AEE / AMI) * fem_int( 9,BphV) * fem_int(0)
    NLC(5,LQi2) = LQi3

    ELM(1:NEMAX,1:4,6,LQi2) =   2.D0 * (PZ * AEE / AMI) * fem_int(16,AphV) &
         &                  +   2.D0 * (PZ * AEE / AMI) * fem_int(20,UNITY,AphV) * fem_int(0)
    NLC(6,LQi2) = LQi4

    NLCMAX(LQi2) = 6

    IF(MDFIXT /= 0) THEN
       ELM(1:NEMAX,1:4,3,LQi2) = - 2.D0 * rKeV / AMI &
            &                  *(  fem_int(17,PTiV) + fem_int(18,PTiV) * fem_int(0) &
            &                    + fem_int(16,PTiV) + fem_int(34,PSI,PTiV) * fem_int(0))
       NLC(3,LQi2) = LQi1
    END IF

    RETURN
  END SUBROUTINE LQi2CC

!***************************************************************
!
!   Ion Poloidal Flow
!
!***************************************************************

  SUBROUTINE LQi3CC

    integer(4) :: N

    ! Ni*UiTheta(0) : 0

    ELM(1:NEMAX,1:4, 0,LQi3) = fem_int(1) * invDT
    NLC( 0,LQi3) = LQi3

    ! Advection

    ELM(1:NEMAX,1:4, 1,LQi3) = - 2.D0 * fem_int(3,RUirV)
    NLC( 1,LQi3) = LQi3

    ! Viscosity force

!    ELM(1:NEMAX,1:4, 2,LQi3) = - 4.D0 * fem_int(18,rMui) &
!         &                     - 4.D0 * fem_int(39,rMuiNi,dPNiV) &
!         &                     - 4.D0 * fem_int( 3,rMui)
    ELM(1:NEMAX,1:4, 2,LQi3) = - 4.D0 * fem_int(18,rMui) &
         &                     - 4.D0 * fem_int( 3,rMui)
    NLC( 2,LQi3) = LQi3

    ELM(1:NEMAX,1:4, 3,LQi3) =   4.D0 * fem_int(41,rMui,RUithV)
    NLC( 3,LQi3) = LQi1

    ! Poroidal E force

    ELM(1:NEMAX,1:4, 4,LQi3) = - (PZ * AEE / AMI) * fem_int(2,PNiV)
    NLC( 4,LQi3) = LQm2

    ! v x B force

    ELM(1:NEMAX,1:4, 5,LQi3) = - (PZ * AEE / AMI) * fem_int(2,BphV)
    NLC( 5,LQi3) = LQi2

    ! Neoclassical viscosity force

    ELM(1:NEMAX,1:4, 6,LQi3) = - fem_int(2,rNuiNC)
    NLC( 6,LQi3) = LQi3

    ! Collisional friction force

    ELM(1:NEMAX,1:4, 7,LQi3) = - (AME / AMI) * fem_int(2,rNuei1EI)
    NLC( 7,LQi3) = LQi3

    ELM(1:NEMAX,1:4, 8,LQi3) =   (AME / AMI) * fem_int(2,rNuei1)
    NLC( 8,LQi3) = LQe3

    ELM(1:NEMAX,1:4, 9,LQi3) =   2.D0 * (AME / AMI) * fem_int(37,rNuei2BthEI,AphV)
    NLC( 9,LQi3) = LQi4

    ELM(1:NEMAX,1:4,10,LQi3) = - 2.D0 * (AME / AMI) * fem_int(37,rNuei2Bth,AphV) * AMPe4
    NLC(10,LQi3) = LQe4

    ! Collisional friction with beam ions

    ELM(1:NEMAX,1:4,11,LQi3) = - (AMB / AMI) * fem_int(2,rNubiBI)
    NLC(11,LQi3) = LQi3

    ELM(1:NEMAX,1:4,12,LQi3) =   (AMB / AMI) * fem_int(2,rNubi)
    NLC(12,LQi3) = LQb3

    IF(MDLWTB == 0) THEN

       ! Wave interaction force (electron driven)

       ELM(1:NEMAX,1:4,13,LQi3) =   1.D0 / AMI * fem_int(2,FWthe)
       NLC(13,LQi3) = LQe3

       ELM(1:NEMAX,1:4,14,LQi3) = - 1.D0 / AMI * fem_int(2,FWthe)
       NLC(14,LQi3) = LQi3

       ELM(1:NEMAX,1:4,15,LQi3) =   2.D0 / AMI * fem_int(36,AphV,FWthphe) * AMPe4
       NLC(15,LQi3) = LQe4

       ELM(1:NEMAX,1:4,16,LQi3) = - 2.D0 / AMI * fem_int(36,AphV,FWthphe)
       NLC(16,LQi3) = LQi4

       ELM(1:NEMAX,1:4,17,LQi3) = - 1.D0 / AMI * fem_int(44,FWthe,WPM)
       NLC(17,LQi3) = LQe1

!!ion       ! Wave interaction force (ion driven)
!!ion
!!ion       ELM(1:NEMAX,1:4,13,LQi3) = - 1.D0 / AMI * fem_int(2,FWthi)
!!ion       NLC(13,LQi3) = LQi3
!!ion
!!ion       ELM(1:NEMAX,1:4,14,LQi3) =   1.D0 / AMI * fem_int(2,FWthi)
!!ion       NLC(14,LQi3) = LQe3
!!ion
!!ion       ELM(1:NEMAX,1:4,15,LQi3) = - 2.D0 / AMI * fem_int(36,AphV,FWthphi)
!!ion       NLC(15,LQi3) = LQi4
!!ion
!!ion       ELM(1:NEMAX,1:4,16,LQi3) =   2.D0 / AMI * fem_int(36,AphV,FWthphi) * AMPe4
!!ion       NLC(16,LQi3) = LQe4
!!ion
!!ion       ELM(1:NEMAX,1:4,17,LQi3) =   1.D0 / AMI * fem_int(44,FWthi,WPM)
!!ion       NLC(17,LQi3) = LQi1

       N = 0
    ELSEIF(MDLWTB == 1) THEN

       ! Wave interaction force (electron driven)

       ELM(1:NEMAX,1:4,13,LQi3) =   2.D0 / AMI * fem_int(17,WNthe)
       NLC(13,LQi3) = LQe1

!       ELM(1:NEMAX,1:4,14,LQi3) =   2.D0 / AMI * fem_int(36,PhiV,WEMthe)
!       NLC(14,LQi3) = LQe1
       
!       ELM(1:NEMAX,1:4,15,LQi3) =   2.D0 / AMI * fem_int(15,WWthe)
!       NLC(15,LQi3) = LQe1
       
!       ELM(1:NEMAX,1:4,16,LQi3) =   2.D0 / AMI * rKeV * fem_int(17,WT1the)
!       NLC(16,LQi3) = LQe5
       
!       ELM(1:NEMAX,1:4,17,LQi3) = - 2.D0 / AMI * fem_int(17,WT2the)
!       NLC(17,LQi3) = LQe1

       ! Wave interaction force (ion driven)

       ELM(1:NEMAX,1:4,18,LQi3) = - 2.D0 / AMI * fem_int(17,WNthi)
       NLC(18,LQi3) = LQi1

!       ELM(1:NEMAX,1:4,19,LQi3) = - 2.D0 / AMI * fem_int(36,PhiV,WEMthi)
!       NLC(19,LQi3) = LQi1
       
!       ELM(1:NEMAX,1:4,20,LQi3) = - 2.D0 / AMI * fem_int(15,WWthi)
!       NLC(20,LQi3) = LQi1
       
!       ELM(1:NEMAX,1:4,21,LQi3) = - 2.D0 / AMI * rKeV * fem_int(17,WT1thi)
!       NLC(21,LQi3) = LQi5
       
!       ELM(1:NEMAX,1:4,22,LQi3) =   2.D0 / AMI * fem_int(17,WT2thi)
!       NLC(22,LQi3) = LQi1

       N = 4
    END IF

    ! Loss to divertor

    ELM(1:NEMAX,1:4,18+N,LQi3) = - 2.D0 * fem_int(2,rNuL) * fact
    NLC(18+N,LQi3) = LQi3

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,19+N,LQi3) = - fem_int(2,rNu0i)
    NLC(19+N,LQi3) = LQi3

    ! Charge exchange force

    ELM(1:NEMAX,1:4,20+N,LQi3) = - fem_int(2,rNuiCX)
    NLC(20+N,LQi3) = LQi3

    ! Loss cone loss

    PELM(1:NEMAX,1:4,21+N,LQi3) = fem_int(-1,SiLCth)
    NLC(21+N,LQi3) = 0

    ! Ion orbit loss

    ELM(1:NEMAX,1:4,22+N,LQi3) = - fem_int(2,rNuOL)
    NLC(22+N,LQi3) = LQi3

    ! Helical Neoclassical viscosity force

    ELM(1:NEMAX,1:4,23+N,LQi3) = - UHth * UHth * fem_int(2,rNuiHL)
    NLC(23+N,LQi3) = LQi3

    ELM(1:NEMAX,1:4,24+N,LQi3) = UHph * UHth * fem_int(2,rNuiHL)
    NLC(24+N,LQi3) = LQi4

    !  Diffusion of ions (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,25+N,LQi3) = - 4.D0 * fem_int(18,DMAGi)
    NLC(25+N,LQi3) = LQi3

    ! Ns*UsTheta(NRMAX) : 0

    NLCMAX(LQi3) = 25+N
    RETURN
  END SUBROUTINE LQi3CC

!***************************************************************
!
!   Ion Toroidal Flow
!
!***************************************************************

  SUBROUTINE LQi4CC

    ! Uiphi'(0) : 0

    ELM(1:NEMAX,1:4, 0,LQi4) = fem_int(1) * invDT
    NLC( 0,LQi4) = LQi4

    ! Advection

    ELM(1:NEMAX,1:4, 1,LQi4) = - 2.D0 * fem_int(3,RUirV)
    NLC( 1,LQi4) = LQi4

    ! Viscosity force

!    ELM(1:NEMAX,1:4, 2,LQi4) = - 4.D0 *(fem_int(18,rMui) + fem_int(39,rMuiNi,dPNiV))
    ELM(1:NEMAX,1:4, 2,LQi4) = - 4.D0 * fem_int(18,rMui)
    NLC( 2,LQi4) = LQi4

    ELM(1:NEMAX,1:4, 3,LQi4) =   4.D0 * fem_int(41,rMui,UiphV)
    NLC( 3,LQi4) = LQi1

    ! Toroidal E force

    ELM(1:NEMAX,1:4, 4,LQi4) = - (PZ * AEE / AMI) * fem_int(2,PNiV)
    NLC( 4,LQi4) = LQm3

    ! v x B force

    ELM(1:NEMAX,1:4, 5,LQi4) = - 2.D0 * (PZ * AEE / AMI) * fem_int(6,AphV)
    NLC( 5,LQi4) = LQi2

    ! Collisional friction with bulk ions

    ELM(1:NEMAX,1:4, 6,LQi4) = - (AME / AMI) * fem_int(2,rNuei3EI)
    NLC( 6,LQi4) = LQi4

    ELM(1:NEMAX,1:4, 7,LQi4) =   (AME / AMI) * fem_int(2,rNuei3) * AMPe4
    NLC( 7,LQi4) = LQe4

    ELM(1:NEMAX,1:4, 8,LQi4) =   2.D0 * (AME / AMI) * fem_int(29,rNuei2BthEI,AphV)
    NLC( 8,LQi4) = LQi3

    ELM(1:NEMAX,1:4, 9,LQi4) = - 2.D0 * (AME / AMI) * fem_int(29,rNuei2Bth,AphV)
    NLC( 9,LQi4) = LQe3

    ! Collisional friction with beam ions

    ELM(1:NEMAX,1:4,10,LQi4) = - (AMB / AMI) * fem_int(2,rNubiBI)
    NLC(10,LQi4) = LQi4

    ELM(1:NEMAX,1:4,11,LQi4) =   (AMB / AMI) * fem_int(2,rNubi)
    NLC(11,LQi4) = LQb4

    IF(MDLWTB == 0) THEN

       ! Wave interaction force (electron driven)

       ELM(1:NEMAX,1:4,12,LQi4) = - 1.D0 / AMI * fem_int(2,FWpheBB)
       NLC(12,LQi4) = LQe3

       ELM(1:NEMAX,1:4,13,LQi4) =   1.D0 / AMI * fem_int(2,FWpheBB)
       NLC(13,LQi4) = LQi3

       ELM(1:NEMAX,1:4,14,LQi4) =   1.D0 / AMI * fem_int(2,FWpheBB2) * AMPe4
       NLC(14,LQi4) = LQe4

       ELM(1:NEMAX,1:4,15,LQi4) = - 1.D0 / AMI * fem_int(2,FWpheBB2) 
       NLC(15,LQi4) = LQi4

       ELM(1:NEMAX,1:4,16,LQi4) =   1.D0 / AMI * fem_int(44,FWpheBB,WPM)
       NLC(16,LQi4) = LQe1

!!ion       ! Wave interaction force (ion driven)
!!ion
!!ion       ELM(1:NEMAX,1:4,12,LQi4) =   1.D0 / AMI * fem_int(2,FWphiBB)
!!ion       NLC(12,LQi4) = LQi3
!!ion
!!ion       ELM(1:NEMAX,1:4,13,LQi4) = - 1.D0 / AMI * fem_int(2,FWphiBB)
!!ion       NLC(13,LQi4) = LQe3
!!ion
!!ion       ELM(1:NEMAX,1:4,14,LQi4) = - 1.D0 / AMI * fem_int(2,FWphiBB2)
!!ion       NLC(14,LQi4) = LQi4
!!ion
!!ion       ELM(1:NEMAX,1:4,15,LQi4) =   1.D0 / AMI * fem_int(2,FWphiBB2) * AMPe4
!!ion       NLC(15,LQi4) = LQe4
!!ion
!!ion       ELM(1:NEMAX,1:4,16,LQi4) = - 1.D0 / AMI * fem_int(44,FWphiBB,WPM)
!!ion       NLC(16,LQi4) = LQi1

    END IF

    ! Loss to divertor

    ELM(1:NEMAX,1:4,17,LQi4) = - 2.D0 * fem_int(2,rNuL) * fact
    NLC(17,LQi4) = LQi4

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,18,LQi4) = - fem_int(2,rNu0i)
    NLC(18,LQi4) = LQi4

    ! Charge exchange force

    ELM(1:NEMAX,1:4,19,LQi4) = - fem_int(2,rNuiCX)
    NLC(19,LQi4) = LQi4

    ! Loss cone loss

    PELM(1:NEMAX,1:4,20,LQi4) =   fem_int(-1,SiLCph)
    NLC(20,LQi4) = 0

    ! Ion orbit loss

    ELM(1:NEMAX,1:4,21,LQi4) = - fem_int(2,rNuOL)
    NLC(21,LQi4) = LQi4

   ! Helical Neoclassical viscosity force

    ELM(1:NEMAX,1:4,22,LQi4) = UHth * UHph * fem_int(2,rNuiHL)
    NLC(22,LQi4) = LQi3

    ELM(1:NEMAX,1:4,23,LQi4) = - UHph * UHph * fem_int(2,rNuiHL)
    NLC(23,LQi4) = LQi4

    !  Diffusion of ions (***AF 2008-06-08)

    ELM(1:NEMAX,1:4,24,LQi4) = - 4.D0 * fem_int(18,DMAGi)
    NLC(24,LQi4) = LQi4

    !  Virtual torque input

    PELM(1:NEMAX,1:4,25,LQi4) =   1.D0 / (AMI * RR * 1.D20) * fem_int(-1,Tqi)
    NLC(25,LQi4) = 0

    NLCMAX(LQi4) = 25
    RETURN
  END SUBROUTINE LQi4CC

!***************************************************************
!
!  Ion Energy Transport: Ti
!
!***************************************************************

  SUBROUTINE LQi5CC

    ! Temperature evolution

    IF(MDFIXT == 0) THEN
       ELM(1:NEMAX,1:4, 0,LQi5) =   1.5D0 * fem_int(1) * invDT
       NLC( 0,LQi5) = LQi5

       ! Advection

       ELM(1:NEMAX,1:4, 1,LQi5) = - 5.D0 * fem_int(3,RUirV)
       NLC( 1,LQi5) = LQi5

       ! Conduction transport

       ELM(1:NEMAX,1:4, 2,LQi5) = - 4.D0 * fem_int(18,Chii1)
       NLC( 2,LQi5) = LQi5

       ELM(1:NEMAX,1:4, 3,LQi5) =   4.D0 * fem_int(41,Chii2,PTiV)
       NLC( 3,LQi5) = LQi1

       ! Redundant heat convection term

       ELM(1:NEMAX,1:4,4,LQi5) = 2.D0 * fem_int(5,RUirV)
       NLC(4,LQi5) = LQi5

       ! Joule heating

       ELM(1:NEMAX,1:4, 5,LQi5) =   PZ * AEE / rKeV * fem_int(2,EthVR)
       ELM(1      ,1:4, 5,LQi5) =   PZ * AEE / rKeV * fem_int_point(2,0,EthVR)
       NLC( 5,LQi5) = LQi3

       ELM(1:NEMAX,1:4, 6,LQi5) =   PZ * AEE / rKeV * fem_int(2,EphV)
       NLC( 6,LQi5) = LQi4

       ! Collisional transfer with electrons (Energy equilibration)

       ELM(1:NEMAX,1:4, 7,LQi5) = - 1.5d0 * fem_int(2,rNuTeiEI)
       NLC( 7,LQi5) = LQi5

       ELM(1:NEMAX,1:4, 8,LQi5) =   1.5d0 * fem_int(2,rNuTei)
       NLC( 8,LQi5) = LQe5

       ! Heating due to beam momentum deposition

!!oldNBI      ELM(1:NEMAX,1:4, 9,LQi5) = - 0.5D0 * AMB * (PNBCD * Vb) / rKeV &
!!oldNBI           &                * fem_int(2,rNubiBI)
!!oldNBI      NLC( 9,LQi5) = LQi4
!!oldNBI
!!oldNBI      ELM(1:NEMAX,1:4,10,LQi5) =   0.5D0 * AMB * (PNBCD * Vb) / rKeV &
!!oldNBI           &                * fem_int(2,rNubi)
!!oldNBI      NLC(10,LQi5) = LQb4

       ELM(1:NEMAX,1:4, 9,LQi5) = AMB * Vb / rKeV * fem_int(28,BthBNi,MNB)
       ELM(1      ,1:4, 9,LQi5) = AMB * Vb / rKeV * fem_int_point(28,0,BthBNi,MNB)
       NLC( 9,LQi5) = LQi3

       ELM(1:NEMAX,1:4,10,LQi5) = AMB * Vb / rKeV * fem_int(28,BphBNi,MNB)
       NLC(10,LQi5) = LQi4

       ! Loss to diverter

       ELM(1:NEMAX,1:4,11,LQi5) = -         1.D0   / PZ * fem_int( 2,rNuL)
       NLC(11,LQi5) = LQi5

       PELM(1:NEMAX,1:4,12,LQi5) =          PNeDIV / PZ * fem_int(-2,rNuL,PTiV)
       NLC(12,LQi5) = 0

       ELM(1:NEMAX,1:4,13,LQi5) = - 1.5D0               * fem_int( 2,rNuLTi)
       NLC(13,LQi5) = LQi5

       PELM(1:NEMAX,1:4,14,LQi5) =  1.5D0 * PTiDIV      * fem_int(-2,rNuLTi,PNiV)
       NLC(14,LQi5) = 0

       ! Direct heating (RF)

       PELM(1:NEMAX,1:4,15,LQi5) = 1.D0 / (1.D20 * rKeV) * fem_int(-1,PRFi)
       NLC(15,LQi5) = 0

       ! Ionization heating of n01 and n02

       ELM(1:NEMAX,1:4,16,LQi5) = 1.5D0 / PZ * fem_int(28,rNuIN0,PT01V)
       NLC(16,LQi5) = LQn1

       ELM(1:NEMAX,1:4,17,LQi5) = 1.5D0 / PZ * fem_int(28,rNuIN0,PT02V)
       NLC(17,LQi5) = LQn2

       ! Charge exchange force

       ELM(1:NEMAX,1:4,18,LQi5)  = - 1.5D0 * fem_int( 2,rNuiCX)
       NLC(18,LQi5) = LQi5

       PELM(1:NEMAX,1:4,19,LQi5) =   1.5D0 * fem_int(-1,rNuCXN1)
       NLC(19,LQi5) = 0

       ! Collisional NBI heating (Perp + Tan)

       PELM(1:NEMAX,1:4,20,LQi5) = Eb * fem_int(-2,SNB,PNBcol_i)
       NLC(20,LQi5) = 0

       ! Simplified Alpha heating

       PELM(1:NEMAX,1:4,21,LQi5) = 1.D0 / (1.D20 * rKeV) * fem_int(-1,PALFi)
       NLC(21,LQi5) = 0

       !  Diffusion of ions (***AF 2008-06-08)

       ELM(1:NEMAX,1:4,22,LQi5) = - 4.D0 * fem_int(18,DMAGi)
       NLC(22,LQi5) = LQi5

       NLCMAX(LQi5) = 22
    ELSE

       !  Fixed temperature profile

       ELM(1:NEMAX,1:4,0,LQi5) = fem_int(1) * invDT
       NLC(0,LQi5) = LQi5

       NLCMAX(LQi5) = 0
    END IF

    RETURN
  END SUBROUTINE LQi5CC

!***************************************************************
!
!   Beam Ion Density
!
!***************************************************************

  SUBROUTINE LQb1CC

    ELM(1:NEMAX,1:4,0,LQb1) = fem_int(1) * invDT
    NLC(0,LQb1) = LQb1

    ! NBI particle source (Both charge exchange and ionization)

    PELM(1:NEMAX,1:4,1,LQb1) =   fem_int(-1,SNBb)
    NLC(1,LQb1) = 0

    ! Extracted NBI perpendicular component

    PELM(1:NEMAX,1:4,2,LQb1) = - fem_int(-2,SNBPDi,rip_rat)
    NLC(2,LQb1) = 0

    ! Relaxation to thermal ions

    ELM(1:NEMAX,1:4,3,LQb1) = - fem_int(2,rNuB)
    NLC(3,LQb1) = LQb1

    ! Loss to divertor

    ELM(1:NEMAX,1:4,4,LQb1) = - fem_int(2,rNuLB) * BeamSW
    NLC(4,LQb1) = LQb1

    ! Ripple trapped beam ions collision with otherwise beam ions

    ELM(1:NEMAX,1:4,5,LQb1) = - fem_int(28,rNubrp2,rip_rat) * RpplSW
    NLC(5,LQb1) = LQb1

    ELM(1:NEMAX,1:4,6,LQb1) =   fem_int(28,rNubrp1,rip_rat) * RpplSW
    NLC(6,LQb1) = LQr1

!!rp_conv    PELM(1:NEMAX,1:4,6,LQb1) =  fem_int(-2,PNbrpLV,rNubLL)
!!rp_conv    NLC(6,LQb1) = 0

    ! Ripple diffusion

    ELM(1:NEMAX,1:4,7,LQb1) = - 4.d0 * fem_int(18,Dbrpft)
    NLC(7,LQb1) = LQb1

    NLCMAX(LQb1) = 7
    RETURN
  END SUBROUTINE LQb1CC

!***************************************************************
!
!   Beam Ion Poloidal Flow
!
!***************************************************************

  SUBROUTINE LQb3CC

    ! Ubth(0) : 0

    ELM(1:NEMAX,1:4,0,LQb3) = fem_int(1) * invDT
    NLC(0,LQb3) = LQb3

    ! Poroidal E force

    ELM(1:NEMAX,1:4,1,LQb3) = - (PZ * AEE / AMB) * fem_int(2,PNbV)
    NLC(1,LQb3) = LQm2

    ! Collisional friction force with electrons

    ELM(1:NEMAX,1:4,2,LQb3) = - fem_int(2,rNube1)
    NLC(2,LQb3) = LQb3

    ELM(1:NEMAX,1:4,3,LQb3) =   fem_int(2,rNube1BE)
    NLC(3,LQb3) = LQe3

    ELM(1:NEMAX,1:4,4,LQb3) =   2.D0 * fem_int(37,rNube2Bth,AphV)
    NLC(4,LQb3) = LQb4

    ELM(1:NEMAX,1:4,5,LQb3) = - 2.D0 * fem_int(37,rNube2BthBE,AphV) * AMPe4
    NLC(5,LQb3) = LQe4

    ! Collisional friction force with ions

    ELM(1:NEMAX,1:4,6,LQb3) = - fem_int(2,rNubi)
    NLC(6,LQb3) = LQb3

    ELM(1:NEMAX,1:4,7,LQb3) =   fem_int(2,rNubiBI)
    NLC(7,LQb3) = LQi3

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,8,LQb3) = - fem_int(2,rNu0b) * BeamSW
    NLC(8,LQb3) = LQb3

    ! Charge exchange force

    ELM(1:NEMAX,1:4,9,LQb3) = - fem_int(2,rNubCX) * BeamSW
    NLC(9,LQb3) = LQb3

    ! NBI momentum source

    PELM(1:NEMAX,1:4,10,LQb3) =  fem_int(-2,RVbparath,MNB)
    NLC(10,LQb3) = 0

    ! Loss to divertor

    ELM(1:NEMAX,1:4,11,LQb3) = - fem_int(2,rNuLB) * BeamSW * fact
    NLC(11,LQb3) = LQb3

    ! Momentum loss due to collisional ripple trapping

    ELM(1:NEMAX,1:4,12,LQb3) = - fem_int(28,rNubrp2,rip_rat) * RpplSW
    NLC(12,LQb3) = LQb3

    ! Momentum diffusion arising from beam ion convective flux due to ripple

    ELM(1:NEMAX,1:4,13,LQb3) = - 4.D0 * fem_int(18,Dbrpft)
    NLC(13,LQb3) = LQb3

    ELM(1:NEMAX,1:4,14,LQb3) = - 4.D0 * fem_int(45,Dbrpft,RUbthV)
    NLC(14,LQb3) = LQb1

!!neo    ! Neoclassical viscosity force
!!neo
!!neo    ELM(1:NEMAX,1:4,15,LQb3) = - fem_int(2,rNuiNC)
!!neo    NLC(15,LQb3) = LQb3

    ! Ubth(NRMAX) : 0

    NLCMAX(LQb3) = 14
    RETURN
  END SUBROUTINE LQb3CC

!***************************************************************
!
!   Beam Ion Toroidal Flow
!
!***************************************************************
 
 SUBROUTINE LQb4CC

    ! - UbPhi(0)' : 0

    ELM(1:NEMAX,1:4,0,LQb4) =   fem_int(1) * invDT
    NLC(0,LQb4) = LQb4

    ! Toroidal E force

    ELM(1:NEMAX,1:4,1,LQb4) = - (PZ * AEE / AMB) * fem_int(2,PNbV)
    NLC(1,LQb4) = LQm3

    ! Collisional friction with electrons

    ELM(1:NEMAX,1:4,2,LQb4) = - fem_int(2,rNube3)
    NLC(2,LQb4) = LQb4

    ELM(1:NEMAX,1:4,3,LQb4) =   fem_int(2,rNube3BE) * AMPe4
    NLC(3,LQb4) = LQe4

    ELM(1:NEMAX,1:4,4,LQb4) =   2.D0 * fem_int(29,rNube2Bth,AphV)
    NLC(4,LQb4) = LQb3

    ELM(1:NEMAX,1:4,5,LQb4) = - 2.D0 * fem_int(29,rNube2BthBE,AphV)
    NLC(5,LQb4) = LQe3

    ! Collisional friction with ions

    ELM(1:NEMAX,1:4,6,LQb4) = - fem_int(2,rNubi)
    NLC(6,LQb4) = LQb4

    ELM(1:NEMAX,1:4,7,LQb4) =   fem_int(2,rNubiBI)
    NLC(7,LQb4) = LQi4

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,8,LQb4) = - fem_int(2,rNu0b) * BeamSW
    NLC(8,LQb4) = LQb4

    ! Charge exchange force

    ELM(1:NEMAX,1:4,9,LQb4) = - fem_int(2,rNubCX) * BeamSW
    NLC(9,LQb4) = LQb4

    ! NBI momentum source

    PELM(1:NEMAX,1:4,10,LQb4) =  fem_int(-2,Vbparaph,MNB)
    NLC(10,LQb4) = 0

    ! Loss to divertor

    ELM(1:NEMAX,1:4,11,LQb4) = - fem_int(2,rNuLB) * BeamSW * fact
    NLC(11,LQb4) = LQb4

    ! Momentum loss due to collisional ripple trapping

    ELM(1:NEMAX,1:4,12,LQb4) = - fem_int(28,rNubrp2,rip_rat) * RpplSW
    NLC(12,LQb4) = LQb4

    ! Momentum diffusion arising from beam ion convective flux due to ripple

    ELM(1:NEMAX,1:4,13,LQb4) = - 4.D0 * fem_int(18,Dbrpft)
    NLC(13,LQb4) = LQb4

    ELM(1:NEMAX,1:4,14,LQb4) = - 4.D0 * fem_int(45,Dbrpft,UbphV)
    NLC(14,LQb4) = LQb1

    ! Ubphi(NRMAX) : 0

    NLCMAX(LQb4) = 14
    RETURN
  END SUBROUTINE LQb4CC

!***************************************************************
!
!   Slow Neutral Transport: n01
!
!***************************************************************

  SUBROUTINE LQn1CC

    ELM(1:NEMAX,1:4,0,LQn1) = fem_int(1) * invDT
    NLC(0,LQn1) = LQn1

    !  Diffusion of neutrals

    ELM(1:NEMAX,1:4,1,LQn1) = - 4.D0 * fem_int(18,D01)
    NLC(1,LQn1) = LQn1

    ! Ionization

    ELM(1:NEMAX,1:4,2,LQn1) = - 1.D0 / PZ * fem_int(2,rNuIN0)
    NLC(2,LQn1) = LQn1

    ! Generation of fast neutrals by charge exchange

    ELM(1:NEMAX,1:4,3,LQn1) = - fem_int(2,rNuCXN0)
    NLC(3,LQn1) = LQn1

    ! Recycling from divertor

    ELM(1:NEMAX,1:4,4,LQn1) =   rGamm0 / PZ * fem_int(2,rNuL)
    NLC(4,LQn1) = LQe1

    PELM(1:NEMAX,1:4,5,LQn1) = - rGamm0 * PNeDIV / PZ * fem_int(-1,rNuL)
    NLC(5,LQn1) = 0

    NLCMAX(LQn1) = 5
    RETURN
  END SUBROUTINE LQn1CC

!***************************************************************
!
!   Fast Neutral Transport: n02
!
!***************************************************************

  SUBROUTINE LQn2CC

    ELM(1:NEMAX,1:4,0,LQn2) = fem_int(1) * invDT
    NLC(0,LQn2) = LQn2

    !  Diffusion of neutrals

    ELM(1:NEMAX,1:4,1,LQn2) = - 4.D0 * fem_int(18,D02)
    NLC(1,LQn2) = LQn2

    ! Ionization

    ELM(1:NEMAX,1:4,2,LQn2) = - 1.D0 / PZ * fem_int(2,rNuIN0)
    NLC(2,LQn2) = LQn2

    ! Generation of fast neutrals by charge exchange

    ELM(1:NEMAX,1:4,3,LQn2) = fem_int(2,rNuCXN0)
    NLC(3,LQn2) = LQn1

    ! NBI particle source (Charge exchange)

    PELM(1:NEMAX,1:4,4,LQn2) = RatCX * fem_int(-1,SNBi)
    NLC(4,LQn2) = 0

    NLCMAX(LQn2) = 4
    RETURN
  END SUBROUTINE LQn2CC

!***************************************************************
!
!   Ripple Trapped Beam Ion Density (SUPG)
!
!***************************************************************

  SUBROUTINE LQr1CC

    integer(4) :: ne
    real(8) :: RUbrpl, SDbrpl, peclet, coef

    do ne = 1, nemax
       RUbrpl = RUbrp(ne-1) + RUbrp(ne)
       SDbrpl = Dbrp(ne-1)*S(ne-1) + Dbrp(ne)*S(ne)
       if (RUbrpl == 0.d0) then
          coef = 0.d0
       elseif (SDbrpl == 0.d0) then
          coef = 1.d0 / sqrt(15.d0) * hpsi(ne)
       else
          peclet = 0.5d0 * RUbrpl * hpsi(ne) / SDbrpl
          coef = 0.5d0 * langevin(peclet) * hpsi(ne)
       end if
       coef=0.d0 !!! no SUPG

       ELM(NE,1:4,0,LQr1) =   fem_int_point(1,NE) * invDT &
            &               + fem_int_point(8,NE) * coef * invDT
       NLC(0,LQr1) = LQr1

       ! NBI perpendicular particle source (Both charge exchange and ionization)
       ! (Beam ions by perpendicular NBI have few parallel velocity, hence
       !  they can be easily trapped by ripple wells if they go into the ripple well region.)
       !  (M.H. Redi, et al., NF 35 (1995) 1191, p.1201 sixth line from the bottom
       !   at right-hand-side column)
       
       PELM(NE,1:4,1,LQr1) =(  fem_int_point(-1,NE,SNBPDi) &
            &                + fem_int_point(-8,NE,SNBPDi) * coef) * RpplSW
       NLC(1,LQr1) = 0

       ! Ripple trapped beam ions collision with otherwise beam ions

       ELM(NE,1:4,2,LQr1) =(  fem_int_point(2,NE,rNubrp2) &
            &               + fem_int_point(9,NE,rNubrp2) * coef) * RpplSW
       NLC(2,LQr1) = LQb1

       ELM(NE,1:4,3,LQr1) =(- fem_int_point(2,NE,rNubrp1) &
            &               - fem_int_point(9,NE,rNubrp1) * coef) * RpplSW
       NLC(3,LQr1) = LQr1

       ! Relaxation to thermal ions

       ELM(NE,1:4,4,LQr1) =(- fem_int_point(2,NE,rNuB) &
            &               - fem_int_point(9,NE,rNuB) * coef) * RpplSW
       NLC(4,LQr1) = LQr1

       ! Ripple loss transport (convective)

       ELM(NE,1:4,5,LQr1) =(- 2.d0 * fem_int_point( 3,NE,RUbrp) &
            &               - 2.d0 * fem_int_point(10,NE,RUbrp) * coef) * RpplSW
       NLC(5,LQr1) = LQr1

!!rp_conv       ELM(NE,1:4,5,LQr1) = - fem_int_point(2,NE,rNubL)
!!rp_conv       NLC(5,LQr1) = LQr1
!!rp_conv       PELM(NE,1:4,5,LQr1) = - fem_int_point(-2,NE,rNubL,PNbrpV)
!!rp_conv       NLC(5,LQr1) = 0

       ! Ripple loss transport (diffusive)

       ELM(NE,1:4,6,LQr1) = - 4.d0 * fem_int_point(18,NE,Dbrp) * RpplSW
       NLC(6,LQr1) = LQr1
    end do

    NLCMAX(LQr1) = 6
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
    type(list), save :: IDX(1:NQM,0:1,0:50)

    IF(ID == 0) THEN
       ! Initialize ALC, BLC and CLC
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

    if (x < -3) then
       y = - 1.d0 - 1.d0 / x
    elseif (x > 3) then
       y =   1.d0 - 1.d0 / x
    else
       y = x / 3.d0 * (1.d0 - abs(x) / 9.d0)
    end if

  end function langevin

end module tx_coefficients
