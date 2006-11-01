!     $Id$
module coefficients
  use commons
  use core_module
  implicit none
  private
  real(8), dimension(:,:,:,:), allocatable :: ELM, PELM
  real(8), dimension(0:NRM) :: rNuIN0, rNuCXN0, rNueiEI, rNubeBE, rNubiBI, rNuTeiEI,&
       &                       ChieNe, ChiiNi, rMueNe, rMuiNi, dPNeV, dPNiV, &
       &                       UethVR, UithVR, EthVR, RUerV, RUirV, UerVR, UirVR, &
       &                       FWpheBB, FWphiBB, dAphV, FWpheBB2, FWphiBB2
  real(8), dimension(0:NRM) :: UNITY = 1.D0
  public :: TXCALA

contains

!***************************************************************
!
!   Calculate ALC, BLC, CLC, PLC
!
!***************************************************************

  SUBROUTINE TXCALA

    INTEGER :: NE, NR, NC, NQ

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

    allocate(ELM(1:4,0:NCM,1:NQMAX,1:NEMAX))
    allocate(PELM(1:4,1:NCM,1:NQMAX,1:NEMAX))

    ELM(1:4,0:NCM,1:NQMAX,1:NEMAX) = 0.D0
    PELM(1:4,1:NCM,1:NQMAX,1:NEMAX) = 0.D0

    ALC(0:NCM,1:NQMAX,0:NRMAX) = 0.D0
    BLC(0:NCM,1:NQMAX,0:NRMAX) = 0.D0
    CLC(0:NCM,1:NQMAX,0:NRMAX) = 0.D0
    NLCR(0:NCM,1:NQMAX,0:NRMAX) = 1
    NLC(0:NCM,1:NQMAX) = 1
    PLC(1:NCM,1:NQMAX,0:NRMAX) = 0.D0
    NLCMAX(1:NQMAX) = 0

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

    !     Elemental equations -> Nodal equations

    DO NE = 1, NEMAX
       NR = NE - 1
       DO NQ = 1, NQMAX
          NC = 0
             BLC(NC,NQ,NR) = BLC(NC,NQ,NR) + ELM(1,NC,NQ,NE)
             ALC(NC,NQ,NR) = ALC(NC,NQ,NR) + ELM(2,NC,NQ,NE)
          DO NC = 1, NLCMAX(NQ)
             BLC(NC,NQ,NR) = BLC(NC,NQ,NR) + ELM(1,NC,NQ,NE)
             ALC(NC,NQ,NR) = ALC(NC,NQ,NR) + ELM(2,NC,NQ,NE)
             PLC(NC,NQ,NR) = PLC(NC,NQ,NR) + PELM(1,NC,NQ,NE) + PELM(2,NC,NQ,NE)
          END DO
       END DO
       NR = NE
       DO NQ = 1, NQMAX
          NC = 0
             CLC(NC,NQ,NR) = CLC(NC,NQ,NR) + ELM(3,NC,NQ,NE)
             BLC(NC,NQ,NR) = BLC(NC,NQ,NR) + ELM(4,NC,NQ,NE)
          DO NC = 1, NLCMAX(NQ)
             CLC(NC,NQ,NR) = CLC(NC,NQ,NR) + ELM(3,NC,NQ,NE)
             BLC(NC,NQ,NR) = BLC(NC,NQ,NR) + ELM(4,NC,NQ,NE)
             PLC(NC,NQ,NR) = PLC(NC,NQ,NR) + PELM(3,NC,NQ,NE) + PELM(4,NC,NQ,NE)
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
    CALL BOUNDARY(NRMAX,LQe2,0)
    CALL BOUNDARY(0    ,LQe3,0)
    CALL BOUNDARY(NRMAX,LQe3,0)
    CALL BOUNDARY(NRMAX,LQe4,0)
    CALL BOUNDARY(0    ,LQi2,0)
    CALL BOUNDARY(NRMAX,LQi2,0)
    CALL BOUNDARY(0    ,LQi3,0)
    CALL BOUNDARY(NRMAX,LQi3,0)
    CALL BOUNDARY(NRMAX,LQi4,0)
    CALL BOUNDARY(NRMAX,LQn2,0)

    !     Integral term

    CALL BOUNDARY(NRMAX,LQm2,1, 2.D0*PSI(NRMAX)*BB)
    CALL BOUNDARY(NRMAX,LQm3,1,-2.D0*R(NRMAX)*Bthb)
    CALL BOUNDARY(NRMAX,LQn1,1, 2.D0*R(NRMAX)*rGASPF)

    deallocate(ELM,PELM)

    RETURN
  END SUBROUTINE TXCALA

!***************************************************************
!
!   Coefficients for Equations
!
!**************************************************************

  SUBROUTINE LQCOEF

    use libraries, only : DERIVS
    INTEGER :: NR
    REAL(8) :: AITKEN4P

    rNuIN0(0:NRMAX)   = rNuION(0:NRMAX) * PNeV(0:NRMAX) &
         &            / (PN01V(0:NRMAX) + PN02V(0:NRMAX))
    rNuCXN0(0:NRMAX)  = rNuiCX(0:NRMAX) * PNiV(0:NRMAX) &
         &            / (PN01V(0:NRMAX) + PN02V(0:NRMAX))
    rNueiEI(0:NRMAX)  = rNuei(0:NRMAX)  * PNeV(0:NRMAX) / PNiV(0:NRMAX)
    rNuTeiEI(0:NRMAX) = rNuTei(0:NRMAX) * PNeV(0:NRMAX) / PNiV(0:NRMAX)
    rNubeBE(0:NRMAX)  = rNube(0:NRMAX)  * PNbV(0:NRMAX) / PNeV(0:NRMAX)
    rNubiBI(0:NRMAX)  = rNubi(0:NRMAX)  * PNbV(0:NRMAX) / PNiV(0:NRMAX)
    ChieNe(0:NRMAX)   = Chie(0:NRMAX)   / PNeV(0:NRMAX)
    ChiiNi(0:NRMAX)   = Chii(0:NRMAX)   / PNiV(0:NRMAX)
    rMueNe(0:NRMAX)   = rMue(0:NRMAX)   / PNeV(0:NRMAX)
    rMuiNi(0:NRMAX)   = rMui(0:NRMAX)   / PNiV(0:NRMAX)
    RUerV(0:NRMAX)    = R(0:NRMAX)      * UerV(0:NRMAX)
    RUirV(0:NRMAX)    = R(0:NRMAX)      * UirV(0:NRMAX)
    UerVR(1:NRMAX)    = UerV(1:NRMAX) / R(1:NRMAX)
    UerVR(0)          = AITKEN4P(PSI(0), &
         &                       UerVR(1),UerVR(2),UerVR(3),UerVR(4),UerVR(5), &
         &                       PSI(1),PSI(2),PSI(3),PSI(4),PSI(5))
    UirVR(1:NRMAX)    = UirV(1:NRMAX) / R(1:NRMAX)
    UirVR(0)          = AITKEN4P(PSI(0), &
         &                       UirVR(1),UirVR(2),UirVR(3),UirVR(4),UirVR(5), &
         &                       PSI(1),PSI(2),PSI(3),PSI(4),PSI(5))
    UethVR(1:NRMAX)   = UethV(1:NRMAX) / R(1:NRMAX)
    UethVR(0)         = AITKEN4P(PSI(0), &
         &                       UethVR(1),UethVR(2),UethVR(3),UethVR(4),UethVR(5), &
         &                       PSI(1),PSI(2),PSI(3),PSI(4),PSI(5))
    UithVR(1:NRMAX)   = UithV(1:NRMAX) / R(1:NRMAX)
    UithVR(0)         = AITKEN4P(PSI(0), &
         &                       UithVR(1),UithVR(2),UithVR(3),UithVR(4),UithVR(5), &
         &                       PSI(1),PSI(2),PSI(3),PSI(4),PSI(5))
    EthVR(1:NRMAX)    = EthV(1:NRMAX) / R(1:NRMAX)
    EthVR(0)          = AITKEN4P(PSI(0), &
         &                       EthVR(1),EthVR(2),EthVR(3),EthVR(4),EthVR(5), &
         &                       PSI(1),PSI(2),PSI(3),PSI(4),PSI(5))
    CALL DERIVS(PSI,X(LQe1,0:NRMAX),dPNeV,NRMAX)
    CALL DERIVS(PSI,X(LQi1,0:NRMAX),dPNiV,NRMAX)
    CALL DERIVS(PSI,X(LQm4,0:NRMAX),dAphV,NRMAX)
    FWpheBB(0:NRMAX)  =- 2.D0 * dAphV(0:NRMAX) * FWthphe(0:NRMAX)
    FWphiBB(0:NRMAX)  =- 2.D0 * dAphV(0:NRMAX) * FWthphi(0:NRMAX)
    FWpheBB2(0:NRMAX) =(BthV(0:NRMAX) / BphV(0:NRMAX))**2 * FWthe(0:NRMAX)
    FWphiBB2(0:NRMAX) =(BthV(0:NRMAX) / BphV(0:NRMAX))**2 * FWthi(0:NRMAX)

  END SUBROUTINE LQCOEF

!***************************************************************
!
!   Poisson Equation: phi
!
!**************************************************************

  SUBROUTINE LQm1CC

    use physical_constants, only : AEE, EPS0

    INTEGER :: NE

    ! phi'(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4,1,LQm1,NE) =   4.D0 * EPS0 / (AEE * 1.D20) * fem_int(18,NE,UNITY)
       NLC(1,LQm1) = LQm1

       ELM(1:4,2,LQm1,NE) =        fem_int(1,NE)
       NLC(2,LQm1) = LQe1

       ELM(1:4,3,LQm1,NE) = - PZ * fem_int(1,NE)
       NLC(3,LQm1) = LQi1

       ELM(1:4,4,LQm1,NE) = - PZ * fem_int(1,NE)
       NLC(4,LQm1) = LQb1
    END DO

    ! phi(b) : 0

    NLCMAX(LQm1) = 4
    RETURN
  END SUBROUTINE LQm1CC

!***************************************************************
!
!   Ampere's Law: Atheta'
!
!***************************************************************

  SUBROUTINE LQm2CC

    use physical_constants, only : AEE, VC, rMU0

    INTEGER :: NE

    ! (r*Atheta)'(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4,0,LQm2,NE) = 1.D0 / (VC**2 * DT) * lump_int(1,NE)
       NLC(0,LQm2) = LQm2

       ! rot Bphi

       ELM(1:4,1,LQm2,NE) = - 4.D0 * fem_int(18,NE,UNITY) &
       &                    - 4.D0 * fem_int( 4,NE)
       NLC(1,LQm2) = LQm5

       ! Electron current

       ELM(1:4,2,LQm2,NE) = - rMU0      * AEE * 1.D20 * fem_int(1,NE)
       NLC(2,LQm2) = LQe3

       ! Ion current

       ELM(1:4,3,LQm2,NE) =   rMU0 * PZ * AEE * 1.D20 * fem_int(1,NE)
       NLC(3,LQm2) = LQi3

       ! Beam ion current

       ELM(1:4,4,LQm2,NE) =   rMU0 * PZ * AEE * 1.D20 * fem_int(1,NE)
       NLC(4,LQm2) = LQb3
    END DO

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

    use physical_constants, only : AEE, VC, rMU0

    INTEGER :: NE

    ! Aphi'(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4,0,LQm3,NE) = 1.D0 / (VC**2 * DT) * lump_int(1,NE)
       NLC(0,LQm3) = LQm3

       ! rot Btheta

       ELM(1:4,1,LQm3,NE) = - 4.D0 * fem_int(18,NE,UNITY)
       NLC(1,LQm3) = LQm4

       ! Electron current

       ELM(1:4,2,LQm3,NE) = - rMU0 *      AEE * 1.D20 * fem_int(1,NE)
       NLC(2,LQm3) = LQe4

       ! Ion current

       ELM(1:4,3,LQm3,NE) =   rMU0 * PZ * AEE * 1.D20 * fem_int(1,NE)
       NLC(3,LQm3) = LQi4

       ! Beam ion current

       ELM(1:4,4,LQm3,NE) =   rMU0 * PZ * AEE * 1.D20 * fem_int(1,NE)
       NLC(4,LQm3) = LQb4

!       ! Virtual current for helical system
!
!       PELM(1:4,5,LQm3,NE) =   rMU0 * fem_int(-1,NE,AJV)
!       NLC(5,LQm3) = 0
    END DO

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

    INTEGER :: NE

    DO NE = 1, NEMAX
       ELM(1:4,0,LQm4,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC(0,LQm4) = LQm4

       ! Aphi'

       ELM(1:4,1,LQm4,NE) = fem_int(1,NE)
       NLC(1,LQm4) = LQm3
    END DO

    NLCMAX(LQm4) = 1
    RETURN
  END SUBROUTINE LQm4CC

!***************************************************************
!
!   Faraday's Law : Atheta
!
!***************************************************************

  SUBROUTINE LQm5CC

    INTEGER :: NE

    DO NE = 1, NEMAX
       ELM(1:4,0,LQm5,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC(0,LQm5) = LQm5

       ! r * Atheta'

       ELM(1:4,1,LQm5,NE) = fem_int(1,NE)
       NLC(1,LQm5) = LQm2
    END DO

    NLCMAX(LQm5) = 1
    RETURN
  END SUBROUTINE LQm5CC

!***************************************************************
!
!   Electron Density Equation
!
!***************************************************************

  SUBROUTINE LQe1CC

    INTEGER :: NE

    DO NE = 1, NEMAX
       ELM(1:4,0,LQe1,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC(0,LQe1) = LQe1

       ! Convection

       ELM(1:4,1,LQe1,NE) = - 2.D0 * fem_int(4,NE)
       NLC(1,LQe1) = LQe2

       ! Ionization of n01 and n02

       ELM(1:4,2,LQe1,NE) =   fem_int(2,NE,rNuIN0)
       NLC(2,LQe1) = LQn1

       ELM(1:4,3,LQe1,NE) =   fem_int(2,NE,rNuIN0)
       NLC(3,LQe1) = LQn2

       ! Loss to divertor

       ELM(1:4,4,LQe1,NE) = - fem_int(2,NE,rNuL)
       NLC(4,LQe1) = LQe1

       PELM(1:4,5,LQe1,NE) =  PNeDIV * fem_int(-1,NE,rNuL)
       NLC(5,LQe1) = 0
    END DO

    NLCMAX(LQe1) = 5
    RETURN
  END SUBROUTINE LQe1CC

!***************************************************************
!
!   Electron Radial Flow
!
!***************************************************************

  SUBROUTINE LQe2CC

    use physical_constants, only : AEE, AME, rKeV

    INTEGER :: NE
    REAL(8) :: P

    ! Ns*Usr(0) : fixed

    DO NE = 1, NEMAX
       ELM(1:4,0,LQe2,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC(0,LQe2) = LQe2

       ! Nonlinear term

       ELM(1:4,1,LQe2,NE) = - 2.D0 * fem_int(3,NE,RUerV) &
            &               +        fem_int(2,NE,UerVR)
       NLC(1,LQe2) = LQe2

       ! Nonlinear centrifugal force

       ELM(1:4,2,LQe2,NE) =   fem_int(2,NE,UethVR)
       NLC(2,LQe2) = LQe3

       ! Pressure gradient force

       ELM(1:4,3,LQe2,NE) = - 2.D0 * rKeV / AME * fem_int(17,NE,UNITY)
       NLC(3,LQe2) = LQe5

       ! Radial E force

       ELM(1:4,4,LQe2,NE) =   2.D0 * (AEE / AME) * fem_int(17,NE,PNeV)
       NLC(4,LQe2) = LQm1

       ! v x B force

       ELM(1:4,5,LQe2,NE) = -        (AEE / AME) * fem_int( 2,NE,BphV)
       NLC(5,LQe2) = LQe3

       ELM(1:4,6,LQe2,NE) = - 2.D0 * (AEE / AME) * fem_int(16,NE,AphV)
       NLC(6,LQe2) = LQe4

       ! *** Streamline Upwind Petrov Galerkin Method ***

       P = HPSI(NE) / SQRT(15.D0)

       ELM(1:4,0,LQe2,NE) = ELM(1:4,0,LQe2,NE) + 1.D0 / DT * lump_int(8,NE) * P

       ELM(1:4,1,LQe2,NE) = ELM(1:4,1,LQe2,NE) &
            &             + (- 2.D0 * fem_int(10,NE,RUerV) + fem_int(9,NE,UerVR)) * P

       ELM(1:4,2,LQe2,NE) = ELM(1:4,2,LQe2,NE) + fem_int( 9,NE,UethVR) * P

       ELM(1:4,3,LQe2,NE) = ELM(1:4,3,LQe2,NE) &
            &             - 2.D0 *  rKeV / AME * fem_int(18,NE,UNITY) * P

       ELM(1:4,4,LQe2,NE) = ELM(1:4,4,LQe2,NE) &
            &             + 2.D0 * (AEE / AME) * fem_int(18,NE,PNeV) * P

       ELM(1:4,5,LQe2,NE) = ELM(1:4,5,LQe2,NE) &
            &             -        (AEE / AME) * fem_int( 9,NE,BphV) * P

       ELM(1:4,6,LQe2,NE) = ELM(1:4,6,LQe2,NE) &
            &             - 2.D0 * (AEE / AME) * fem_int(20,NE,UNITY,AphV) * P
    END DO

    ! Ns*Usr(NRMAX) : fixed or finite gradient

    NLCMAX(LQe2) = 6
    RETURN
  END SUBROUTINE LQe2CC

!***************************************************************
!
!   Electron Poloidal Flow
!
!***************************************************************

  SUBROUTINE LQe3CC

    use physical_constants, only : AEE, AME, rKeV

    INTEGER :: NE, N
    REAL(8) :: P

    ! Ns*UsTheta(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4, 0,LQe3,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC( 0,LQe3) = LQe3

       ! Nonlinear term

       ELM(1:4, 1,LQe3,NE) = - 2.D0 * fem_int(3,NE,RUerV)
       NLC( 1,LQe3) = LQe3

       ! Viscosity force

       ELM(1:4, 2,LQe3,NE) = - 4.D0 * fem_int(18,NE,rMue) &
            &                - 4.D0 * fem_int(39,NE,rMueNe,dPNeV) &
            &                - 4.D0 * fem_int( 3,NE,rMue)
       NLC( 2,LQe3) = LQe3

       ! Poloidal E force

       ELM(1:4, 3,LQe3,NE) =   (AEE / AME) * fem_int(2,NE,PNeV)
       NLC( 3,LQe3) = LQm2

       ! v x B force

       ELM(1:4, 4,LQe3,NE) =   (AEE / AME) * fem_int(2,NE,BphV)
       NLC( 4,LQe3) = LQe2

       ! Neoclassical viscosity force

       ELM(1:4, 5,LQe3,NE) = - fem_int(2,NE,rNueNC)
       NLC( 5,LQe3) = LQe3

       ! Collisional friction force with ions

       ELM(1:4, 6,LQe3,NE) = - fem_int(2,NE,rNuei)
       NLC( 6,LQe3) = LQe3

       ELM(1:4, 7,LQe3,NE) =   fem_int(2,NE,rNueiEI)
       NLC( 7,LQe3) = LQi3

       ! Collisional friction with beam ions

       ELM(1:4, 8,LQe3,NE) = - (AMB / AME) * fem_int(2,NE,rNubeBE)
       NLC( 8,LQe3) = LQe3

       ELM(1:4, 9,LQe3,NE) =   (AMB / AME) * fem_int(2,NE,rNube)
       NLC( 9,LQe3) = LQb3

       IF(MDLWTB == 0) THEN

          ! Wave interaction force (electron driven)

          ELM(1:4,10,LQe3,NE) = - 1.D0 / AME * fem_int(2,NE,FWthe)
          NLC(10,LQe3) = LQe3

          ELM(1:4,11,LQe3,NE) = - 2.D0 / AME * fem_int(36,NE,AphV,FWthphe)
          NLC(11,LQe3) = LQe4

          ELM(1:4,12,LQe3,NE) =   1.D0 / AME * fem_int(44,NE,FWthe,WPM)
          NLC(12,LQe3) = LQe1

          ! Wave interaction force (ion driven)

          ELM(1:4,13,LQe3,NE) =   1.D0 / AME * fem_int(2,NE,FWthi)
          NLC(13,LQe3) = LQi3

          ELM(1:4,14,LQe3,NE) =   2.D0 / AME * fem_int(36,NE,AphV,FWthphi)
          NLC(14,LQe3) = LQi4

          ELM(1:4,15,LQe3,NE) = - 1.D0 / AME * fem_int(44,NE,FWthi,WPM)
          NLC(15,LQe3) = LQi1

          N = 0
       ELSEIF(MDLWTB == 1) THEN

          ! Wave interaction force (electron driven)

          ELM(1:4,10,LQe3,NE) = - 2.D0 / AME * fem_int(17,NE,WNthe)
          NLC(10,LQe3) = LQe1

!!$          ELM(1:4,11,LQe3,NE) = - 2.D0 / AME * fem_int(36,NE,Phi,WEMthe)
!!$          NLC(11,LQe3) = LQe1
!!$
!!$          ELM(1:4,12,LQe3,NE) = - 2.D0 / AME * fem_int(15,NE,WWthe)
!!$          NLC(12,LQe3) = LQe1
!!$
!!$          ELM(1:4,13,LQe3,NE) = - 2.D0 / AME * rKeV * fem_int(17,NE,WT1the)
!!$          NLC(13,LQe3) = LQe5
!!$
!!$          ELM(1:4,14,LQe3,NE) =   2.D0 / AME * fem_int(17,NE,WT2the)
!!$          NLC(14,LQe3) = LQe1

          ! Wave interaction force (ion driven)

          ELM(1:4,15,LQe3,NE) =   2.D0 / AME * fem_int(17,NE,WNthi)
          NLC(15,LQe3) = LQi1

!!$          ELM(1:4,16,LQe3,NE) =   2.D0 / AME * fem_int(36,NE,Phi,WEMthi)
!!$          NLC(16,LQe3) = LQi1
!!$
!!$          ELM(1:4,17,LQe3,NE) =   1.D0 / AME * fem_int(15,NE,WWthi)
!!$          NLC(17,LQe3) = LQi1
!!$
!!$          ELM(1:4,18,LQe3,NE) =   2.D0 / AME * rKeV * fem_int(17,NE,WT1thi)
!!$          NLC(18,LQe3) = LQi5
!!$
!!$          ELM(1:4,19,LQe3,NE) = - 2.D0 / AME * fem_int(17,NE,WT2thi)
!!$          NLC(19,LQe3) = LQi1

          N = 4
       END IF

       ! Loss to divertor

       ELM(1:4,16+N,LQe3,NE) = - fem_int(2,NE,rNuL)
       NLC(16+N,LQe3) = LQe3

       ! Collisional friction force with neutrals

       ELM(1:4,17+N,LQe3,NE) = - fem_int(2,NE,rNu0e)
       NLC(17+N,LQe3) = LQe3

!       ! Helical neoclassical viscosity force
!
!       ELM(1:4,18+N,LQe3,NE) = - (1.D0 - UHth * UHth) * fem_int(2,NE,rNueHL)
!       NLC(18+N,LQe3) = LQe3
!
!       ELM(1:4,19+N,LQe3,NE) = UHph * UHth / 2.D0 * fem_int(22,NE,rNueHL)
!       NLC(19+N,LQe3) = LQe4
    END DO

    ! Ns*UsTheta(NRMAX) : 0

    NLCMAX(LQe3) = 17+N
    RETURN
  END SUBROUTINE LQe3CC

!***************************************************************
!
!   Electron Toroidal Flow
!
!***************************************************************

  SUBROUTINE LQe4CC

    use physical_constants, only : AEE, AME

    INTEGER :: NE

    ! Uephi(0)' : 0

    DO NE = 1, NEMAX
       ELM(1:4, 0,LQe4,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC( 0,LQe4) = LQe4

       ! Nonlinear term

       ELM(1:4, 1,LQe4,NE) = - 2.D0 * fem_int(3,NE,RUerV)
       NLC( 1,LQe4) = LQe4

       ! Viscosity force

       ELM(1:4, 2,LQe4,NE) = - 4.D0 * fem_int(18,NE,rMue) &
            &                - 4.D0 * fem_int(39,NE,rMueNe,dPNeV)
       NLC( 2,LQe4) = LQe4

       ! Toroidal E force

       ELM(1:4, 3,LQe4,NE) =   (AEE / AME) * fem_int(2,NE,PNeV)
       NLC( 3,LQe4) = LQm3

       ! v x B force

       ELM(1:4, 4,LQe4,NE) =   2.D0 * (AEE / AME) * fem_int(6,NE,AphV)
       NLC( 4,LQe4) = LQe2

       ! Collisional friction with bulk ions

       ELM(1:4, 5,LQe4,NE) = - fem_int(2,NE,rNuei)
       NLC( 5,LQe4) = LQe4

       ELM(1:4, 6,LQe4,NE) =   fem_int(2,NE,rNueiEI)
       NLC( 6,LQe4) = LQi4

       ! Collisional friction with beam ions

       ELM(1:4, 7,LQe4,NE) = - (AMB / AME) * fem_int(2,NE,rNubeBE)
       NLC( 7,LQe4) = LQe4

       ELM(1:4, 8,LQe4,NE) =   (AMB / AME) * fem_int(2,NE,rNube)
       NLC( 8,LQe4) = LQb4

       IF(MDLWTB == 0) THEN

          ! Wave interaction force (electron driven)

          ELM(1:4, 9,LQe4,NE) =   1.D0 / AME * fem_int(2,NE,FWpheBB)
          NLC( 9,LQe4) = LQe3

          ELM(1:4,10,LQe4,NE) = - 1.D0 / AME * fem_int(2,NE,FWpheBB2)
          NLC(10,LQe4) = LQe4

          ELM(1:4,11,LQe4,NE) = - 1.D0 / AME * fem_int(44,NE,FWpheBB,WPM)
          NLC(11,LQe4) = LQe1

          ! Wave interaction force (ion driven)

          ELM(1:4,12,LQe4,NE) = - 1.D0 / AME * fem_int(2,NE,FWphiBB)
          NLC(12,LQe4) = LQi3

          ELM(1:4,13,LQe4,NE) =   1.D0 / AME * fem_int(2,NE,FWphiBB2)
          NLC(13,LQe4) = LQi4

          ELM(1:4,14,LQe4,NE) =   1.D0 / AME * fem_int(44,NE,FWphiBB,WPM)
          NLC(14,LQe4) = LQi1

       END IF

       ! Loss to divertor

       ELM(1:4,15,LQe4,NE) = - fem_int(2,NE,rNuL)
       NLC(15,LQe4) = LQe4

       ! Collisional friction force with neutrals

       ELM(1:4,16,LQe4,NE) = - fem_int(2,NE,rNu0e)
       NLC(16,LQe4) = LQe4

!       ! Helical neoclassical viscosity force
!
!       ELM(1:4,17,LQe4,NE) =  UHth * UHph / 2.D0 * fem_int(2,NE,rNueHL)
!       NLC(17,LQe4) = LQe3
!
!       ELM(1:4,18,LQe4,NE) = - (1.D0 - UHph * UHph) * fem_int(2,NE,rNueHL)
!       NLC(18,LQe4) = LQe4
    END DO

    NLCMAX(LQe4) = 16
    RETURN
  END SUBROUTINE LQe4CC

!***************************************************************
!
!   Electron Energy Transport: Te
!
!***************************************************************

  SUBROUTINE LQe5CC

    use physical_constants, only : AEE, rKeV

    INTEGER :: NE
    REAL(8) :: Veff, Deff, Pec, P

    ! Temperature evolution
    
    IF(MDFIXT == 0) THEN
       DO NE = 1, NEMAX
          ELM(1:4, 0,LQe5,NE) =   1.5D0 / DT * lump_int(1,NE)
          NLC( 0,LQe5) = LQe5

          ! Convection transport

          ELM(1:4, 1,LQe5,NE) = - 5.D0 * fem_int(3,NE,RUerV)
          NLC( 1,LQe5) = LQe5

          ! Conduction transport

          ELM(1:4, 2,LQe5,NE) = - 6.D0 * fem_int(18,NE,Chie) &
               &                - 6.D0 * fem_int(39,NE,ChieNe,dPNeV)
          NLC( 2,LQe5) = LQe5

          ! Joule heating

          ELM(1:4, 3,LQe5,NE) = - AEE / rKeV * fem_int(2,NE,EthVR)
          NLC( 3,LQe5) = LQe3

          ELM(1:4, 4,LQe5,NE) = - AEE / rKeV * fem_int(2,NE,EphV)
          NLC( 4,LQe5) = LQe4

          ! Collisional transfer with ions

          ELM(1:4, 5,LQe5,NE) = - fem_int(2,NE,rNuTei)
          NLC( 5,LQe5) = LQe5

          ELM(1:4, 6,LQe5,NE) =   fem_int(2,NE,rNuTeiEI)
          NLC( 6,LQe5) = LQi5

          ! Collisional heating with beam

          ELM(1:4, 7,LQe5,NE) = - 0.5D0 * AMb * (PNBCD * Vb) / rKeV &
               &                * fem_int(2,NE,rNubeBE)
          NLC( 7,LQe5) = LQe4

          ELM(1:4, 8,LQe5,NE) =   0.5D0 * AMb * (PNBCD * Vb) / rKeV &
               &                * fem_int(2,NE,rNube)
          NLC( 8,LQe5) = LQb4

          ! Loss to diverter

          ELM(1:4, 9,LQe5,NE) = -                  fem_int( 2,NE,rNuL)
          NLC( 9,LQe5) = LQe5

          PELM(1:4,10,LQe5,NE) =          PNeDIV * fem_int(-2,NE,rNuLTe,PTeV)
          NLC(10,LQe5) = 0

          ELM(1:4,11,LQe5,NE) = - 1.5D0          * fem_int( 2,NE,rNuLTe)
          NLC(11,LQe5) = LQe5

          PELM(1:4,12,LQe5,NE) =  1.5D0 * PTeDIV * fem_int(-2,NE,rNuLTe,PNeV)
          NLC(12,LQe5) = 0

          ! Direct heating (RF)

          PELM(1:4,13,LQe5,NE) =   1.D0 / (1.D20 * rKeV) * fem_int(-1,NE,PRFe)
          NLC(13,LQe5) = 0

          ! Radiation loss

          PELM(1:4,14,LQe5,NE) = - 1.D0 / (1.D20 * rKeV) * fem_int(-1,NE,PBr)
          NLC(14,LQe5) = 0

!!$          ! *** Streamline Upwind Petrov Galerkin Method ***
!!$          
!!$          Veff = 5.D0 * (0.5D0 * (UerVR (NE-1) + UerVR (NE))) &
!!$               & + 6.D0 * (0.5D0 * (ChieNe(NE-1) + ChieNe(NE)))  &
!!$!               &         * 0.5D0 * (dPNeV (NE-1) + dPNeV (NE)))
!!$               &         *(PNeV (NE-1) - PNeV (NE))/HPSI(NE)
!!$          Deff = 6.D0 * (0.5D0 * (Chie  (NE-1) + Chie  (NE)))
!!$          Pec  = 0.5D0 * Veff * HPSI(NE) / Deff
!!$          P    = 0.5D0 * HPSI(NE) * falpha(Pec)
!!$
!!$          ELM(1:4, 0,LQe5,NE) = ELM(1:4, 0,LQe5,NE) + 1.5D0 / DT * lump_int(8,NE) * P
!!$
!!$          ELM(1:4, 1,LQe5,NE) = ELM(1:4, 1,LQe5,NE) - 5.D0 * fem_int(10,NE,RUerV) * P
!!$
!!$!          ELM(1:4, 2,LQe5,NE) = ELM(1:4, 2,LQe5,NE) - 6.D0 * fem_int(42,NE,ChieNe,dPNeV) * P
!!$          ELM(1:4, 2,LQe5,NE) = ELM(1:4, 2,LQe5,NE) - 6.D0 * fem_int(46,NE,ChieNe,PNeV) * P
!!$
!!$          ELM(1:4, 3,LQe5,NE) = ELM(1:4, 3,LQe5,NE) - AEE / rKeV * fem_int(9,NE,EthVR) * P
!!$
!!$          ELM(1:4, 4,LQe5,NE) = ELM(1:4, 4,LQe5,NE) - AEE / rKeV * fem_int(9,NE,EphV) * P
!!$
!!$          ELM(1:4, 5,LQe5,NE) = ELM(1:4, 5,LQe5,NE) - fem_int(9,NE,rNuTei) * P
!!$
!!$          ELM(1:4, 6,LQe5,NE) = ELM(1:4, 6,LQe5,NE) + fem_int(9,NE,rNuTeiEI) * P
!!$
!!$          ELM(1:4, 7,LQe5,NE) = ELM(1:4, 7,LQe5,NE) - 0.5D0 * AMb * (PNBCD * Vb) / rKeV &
!!$               &                * fem_int(9,NE,rNubeBE) * P
!!$
!!$          ELM(1:4, 8,LQe5,NE) = ELM(1:4, 8,LQe5,NE) + 0.5D0 * AMb * (PNBCD * Vb) / rKeV &
!!$               &                * fem_int(9,NE,rNube) * P
!!$
!!$          ELM(1:4, 9,LQe5,NE) = ELM(1:4, 9,LQe5,NE) - fem_int(9,NE,rNuL) * P
!!$
!!$          PELM(1:4,10,LQe5,NE) = PELM(1:4,10,LQe5,NE) + PNeDIV * fem_int(-9,NE,rNuLTe,PTeV) * P
!!$
!!$          ELM(1:4,11,LQe5,NE) = ELM(1:4,11,LQe5,NE) - 1.5D0 * fem_int(9,NE,rNuLTe) * P
!!$
!!$          PELM(1:4,12,LQe5,NE) = PELM(1:4,12,LQe5,NE) &
!!$               &               + 1.5D0 * PTeDIV * fem_int(-9,NE,rNuLTe,PNeV) * P
!!$
!!$          PELM(1:4,13,LQe5,NE) = PELM(1:4,13,LQe5,NE) &
!!$               &               + 1.D0 / (1.D20 * rKeV) * fem_int(-8,NE,PRFe) * P
!!$
!!$          PELM(1:4,14,LQe5,NE) = PELM(1:4,14,LQe5,NE) &
!!$               &               - 1.D0 / (1.D20 * rKeV) * fem_int(-8,NE,PBr) * P
       END DO

       NLCMAX(LQe5) = 14
    ELSE

       !  Fixed temperature profile

       DO NE = 1, NEMAX
          ELM(1:4,0,LQe5,NE) = 1.D0 / DT * lump_int(1,NE)
          NLC(0,LQe5) = LQe5

          ! Convection

          ELM(1:4,1,LQe5,NE) = - 2.D0 * fem_int(5,NE,PTeV)
          NLC(1,LQe5) = LQe2

          ! Ionization of n01 and n02

          ELM(1:4,2,LQe5,NE) =   fem_int(28,NE,rNuIN0,PTeV)
          NLC(2,LQe5) = LQn1

          ELM(1:4,3,LQe5,NE) =   fem_int(28,NE,rNuIN0,PTeV)
          NLC(3,LQe5) = LQn2

          ! Loss to divertor

          ELM(1:4,4,LQe5,NE) = - fem_int(2,NE,rNuL)
          NLC(4,LQe5) = LQe5

          PELM(1:4,5,LQe5,NE) =  PNeDIV * fem_int(-2,NE,rNuL,PTeV)
          NLC(5,LQe5) = 0
       END DO

       NLCMAX(LQe5) = 5
    END IF

    RETURN
  END SUBROUTINE LQe5CC

!***************************************************************
!
!   Ion Density Equation
!
!***************************************************************

  SUBROUTINE LQi1CC

    INTEGER :: NE

    DO NE = 1, NEMAX
       ELM(1:4,0,LQi1,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC(0,LQi1) = LQi1

       ! Convection

       ELM(1:4,1,LQi1,NE) = - 2.D0 * fem_int(4,NE)
       NLC(1,LQi1) = LQi2

       ! Ionization of n01 and n02

       ELM(1:4,2,LQi1,NE) =     1.D0 / PZ * fem_int(2,NE,rNuIN0)
       NLC(2,LQi1) = LQn1

       ELM(1:4,3,LQi1,NE) =     1.D0 / PZ * fem_int(2,NE,rNuIN0)
       NLC(3,LQi1) = LQn2

       ! Loss to divertor

       ELM(1:4,4,LQi1,NE) = -   1.D0 / PZ * fem_int(2,NE,rNuL)
       NLC(4,LQi1) = LQe1

       PELM(1:4,5,LQi1,NE) =  PNeDIV / PZ * fem_int(-1,NE,rNuL)
       NLC(5,LQi1) = 0

       ! Particle source from beam ion

       ELM(1:4,6,LQi1,NE) =   fem_int(2,NE,rNuB)
       NLC(6,LQi1) = LQb1

       ! NBI kick up ions

       PELM(1:4,7,LQi1,NE) = - fem_int(-1,NE,SNB)
       NLC(7,LQi1) = 0

       ! Loss cone loss

       PELM(1:4,8,LQi1,NE) =   fem_int(-1,NE,SiLC)
       NLC(8,LQi1) = 0
    END DO

    NLCMAX(LQi1) = 8
    RETURN
  END SUBROUTINE LQi1CC

!***************************************************************
!
!   Ion Radial Flow
!
!***************************************************************
  
  SUBROUTINE LQi2CC

    use physical_constants, only : AEE, rKeV

    INTEGER :: NE
    REAL(8) :: P

    ! Ns*Usr(0) : fixed

    DO NE = 1, NEMAX
       ELM(1:4,0,LQi2,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC(0,LQi2) = LQi2

       ! Nonlinear term

       ELM(1:4,1,LQi2,NE) = - 2.D0 * fem_int(3,NE,RUirV) &
            &               +        fem_int(2,NE,UirVR)
       NLC(1,LQi2) = LQi2

       ! Nonlinear centrifugal force

       ELM(1:4,2,LQi2,NE) =   fem_int(2,NE,UithVR)
       NLC(2,LQi2) = LQi3

       ! Pressure gradient force

       ELM(1:4,3,LQi2,NE) = - 2.D0 * rKeV / AMI * fem_int(17,NE,UNITY)
       NLC(3,LQi2) = LQi5

       ! Radial E force

       ELM(1:4,4,LQi2,NE) = - 2.D0 * (PZ * AEE / AMI) * fem_int(17,NE,PNiV)
       NLC(4,LQi2) = LQm1

       ! v x B force

       ELM(1:4,5,LQi2,NE) =          (PZ * AEE / AMI) * fem_int( 2,NE,BphV)
       NLC(5,LQi2) = LQi3

       ELM(1:4,6,LQi2,NE) =   2.D0 * (PZ * AEE / AMI) * fem_int(16,NE,AphV)
       NLC(6,LQi2) = LQi4

       ! *** Streamline Upwind Petrov Galerkin Method ***

       P = HPSI(NE) / SQRT(15.D0)

       ELM(1:4,0,LQi2,NE) = ELM(1:4,0,LQi2,NE) + 1.D0 / DT * lump_int(8,NE) * P

       ELM(1:4,1,LQi2,NE) = ELM(1:4,1,LQi2,NE) &
            &             + (- 2.D0 * fem_int(10,NE,RUirV) + fem_int(9,NE,UirVR)) * P

       ELM(1:4,2,LQi2,NE) = ELM(1:4,2,LQi2,NE) + fem_int( 9,NE,UithVR) * P

       ELM(1:4,3,LQi2,NE) = ELM(1:4,3,LQi2,NE) &
            &             - 2.D0 * rKeV / AMI * fem_int(18,NE,UNITY) * P

       ELM(1:4,4,LQi2,NE) = ELM(1:4,4,LQi2,NE) &
            &             - 2.D0 * (PZ * AEE / AMI) * fem_int(18,NE,PNiV) * P

       ELM(1:4,5,LQi2,NE) = ELM(1:4,5,LQi2,NE) &
            &             +        (PZ * AEE / AMI) * fem_int( 9,NE,BphV) * P

       ELM(1:4,6,LQi2,NE) = ELM(1:4,6,LQi2,NE) &
            &             + 2.D0 * (PZ * AEE / AMI) * fem_int(20,NE,UNITY,AphV) * P
    END DO

    ! Ns*Usr(NRMAX) : fixed or finite gradient

    NLCMAX(LQi2) = 6
    RETURN
  END SUBROUTINE LQi2CC

!***************************************************************
!
!   Ion Poloidal Flow
!
!***************************************************************

  SUBROUTINE LQi3CC

    use physical_constants, only : AEE, AME, rKeV

    INTEGER :: NE, N
    REAL(8) :: P

    ! Ni*UiTheta(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4, 0,LQi3,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC( 0,LQi3) = LQi3

       ! Nonlinear term

       ELM(1:4, 1,LQi3,NE) = - 2.D0 * fem_int(3,NE,RUirV)
       NLC( 1,LQi3) = LQi3

       ! Viscosity force

       ELM(1:4, 2,LQi3,NE) = - 4.D0 * fem_int(18,NE,rMui) &
            &                - 4.D0 * fem_int(39,NE,rMuiNi,dPNiV) &
            &                - 4.D0 * fem_int( 3,NE,rMui)
       NLC( 2,LQi3) = LQi3

       ! Poroidal E force

       ELM(1:4, 3,LQi3,NE) = - (PZ * AEE / AMI) * fem_int(2,NE,PNiV)
       NLC( 3,LQi3) = LQm2

       ! v x B force

       ELM(1:4, 4,LQi3,NE) = - (PZ * AEE / AMI) * fem_int(2,NE,BphV)
       NLC( 4,LQi3) = LQi2

       ! Neoclassical viscosity force

       ELM(1:4, 5,LQi3,NE) = - fem_int(2,NE,rNuiNC)
       NLC( 5,LQi3) = LQi3

       ! Collisional friction force

       ELM(1:4, 6,LQi3,NE) = - (AME / AMI) * fem_int(2,NE,rNueiEI)
       NLC( 6,LQi3) = LQi3

       ELM(1:4, 7,LQi3,NE) =   (AME / AMI) * fem_int(2,NE,rNuei)
       NLC( 7,LQi3) = LQe3

       ! Collisional friction with beam ions

       ELM(1:4, 8,LQi3,NE) = - (AMB / AMI) * fem_int(2,NE,rNubiBI)
       NLC( 8,LQi3) = LQi3

       ELM(1:4, 9,LQi3,NE) =   (AMB / AMI) * fem_int(2,NE,rNubi)
       NLC( 9,LQi3) = LQb3

       IF(MDLWTB == 0) THEN

          ! Wave interaction force (electron driven)

          ELM(1:4,10,LQi3,NE) =   1.D0 / AMI * fem_int(2,NE,FWthe)
          NLC(10,LQi3) = LQe3

          ELM(1:4,11,LQi3,NE) =   2.D0 / AMI * fem_int(36,NE,AphV,FWthphe)
          NLC(11,LQi3) = LQe4

          ELM(1:4,12,LQi3,NE) = - 1.D0 / AMI * fem_int(44,NE,FWthe,WPM)
          NLC(12,LQi3) = LQe1

          ! Wave interaction force (ion driven)

          ELM(1:4,13,LQi3,NE) = - 1.D0 / AMI * fem_int(2,NE,FWthi)
          NLC(13,LQi3) = LQi3

          ELM(1:4,14,LQi3,NE) = - 2.D0 / AMI * fem_int(36,NE,AphV,FWthphi)
          NLC(14,LQi3) = LQi4

          ELM(1:4,15,LQi3,NE) =   1.D0 / AMI * fem_int(44,NE,FWthi,WPM)
          NLC(15,LQi3) = LQi1

          N = 0
       ELSEIF(MDLWTB == 1) THEN

          ! Wave interaction force (electron driven)

          ELM(1:4,10,LQi3,NE) =   2.D0 / AMI * fem_int(17,NE,WNthe)
          NLC(10,LQi3) = LQe1

!!$          ELM(1:4,11,LQi3,NE) =   2.D0 / AMI * fem_int(36,NE,Phi,WEMthe)
!!$          NLC(11,LQi3) = LQe1
!!$
!!$          ELM(1:4,12,LQi3,NE) =   2.D0 / AMI * fem_int(15,NE,WWthe)
!!$          NLC(12,LQi3) = LQe1
!!$
!!$          ELM(1:4,13,LQi3,NE) =   2.D0 / AMI * rKeV * fem_int(17,NE,WT1the)
!!$          NLC(13,LQi3) = LQe5
!!$
!!$          ELM(1:4,14,LQi3,NE) = - 2.D0 / AMI * fem_int(17,NE,WT2the)
!!$          NLC(14,LQi3) = LQe1

          ! Wave interaction force (ion driven)

          ELM(1:4,15,LQi3,NE) = - 2.D0 / AMI * fem_int(17,NE,WNthi)
          NLC(15,LQi3) = LQi1

!!$          ELM(1:4,16,LQi3,NE) = - 2.D0 / AMI * fem_int(36,NE,Phi,WEMthi)
!!$          NLC(16,LQi3) = LQi1
!!$
!!$          ELM(1:4,17,LQi3,NE) = - 2.D0 / AMI * fem_int(15,NE,WWthi)
!!$          NLC(17,LQi3) = LQi1
!!$
!!$          ELM(1:4,18,LQi3,NE) = - 2.D0 / AMI * rKeV * fem_int(17,NE,WT1thi)
!!$          NLC(18,LQi3) = LQi5
!!$
!!$          ELM(1:4,19,LQi3,NE) =   2.D0 / AMI * fem_int(17,NE,WT2thi)
!!$          NLC(19,LQi3) = LQi1

          N = 4
       END IF

       ! Loss to divertor

       ELM(1:4,16+N,LQi3,NE) = - fem_int(2,NE,rNuL)
       NLC(16+N,LQi3) = LQi3

       ! Collisional friction force with neutrals

       ELM(1:4,17+N,LQi3,NE) = - fem_int(2,NE,rNu0i)
       NLC(17+N,LQi3) = LQi3

       ! Charge exchange force

       ELM(1:4,18+N,LQi3,NE) = - fem_int(2,NE,rNuiCX)
       NLC(18+N,LQi3) = LQi3

       ! Loss cone loss

       PELM(1:4,19+N,LQi3,NE) = fem_int(21,NE,SiLCth)
       NLC(19+N,LQi3) = 0

!       ! Helical Neoclassical viscosity force
!
!       ELM(1:4,20+N,LQi3,NE) = - (1.D0 - UHth * UHth) * fem_int(2,NE,rNuiHL)
!       NLC(20+N,LQi3) = LQi3
!
!       ELM(1:4,21+N,LQi3,NE) = UHph * UHth * fem_int(22,NE,rNuiHL)
!       NLC(21+N,LQi3) = LQi4
    END DO

    ! Ns*UsTheta(NRMAX) : 0

    NLCMAX(LQi3) = 19+N
    RETURN
  END SUBROUTINE LQi3CC

!***************************************************************
!
!   Ion Toroidal Flow
!
!***************************************************************

  SUBROUTINE LQi4CC

    use physical_constants, only : AEE, AME

    INTEGER :: NE, N

    ! Uiphi'(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4, 0,LQi4,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC( 0,LQi4) = LQi4

       ! Nonlinear term
       
       ELM(1:4, 1,LQi4,NE) = - 2.D0 * fem_int(3,NE,RUirV)
       NLC( 1,LQi4) = LQi4

       ! Viscosity force

       ELM(1:4, 2,LQi4,NE) = - 4.D0 * fem_int(18,NE,rMui) &
            &                - 4.D0 * fem_int(39,NE,rMuiNi,dPNiV)
       NLC( 2,LQi4) = LQi4

       ! Toroidal E force

       ELM(1:4, 3,LQi4,NE) = - (PZ * AEE / AMI) * fem_int(2,NE,PNiV)
       NLC( 3,LQi4) = LQm3

       ! v x B force

       ELM(1:4, 4,LQi4,NE) = - 2.D0 * (PZ * AEE / AMI) * fem_int(6,NE,AphV)
       NLC( 4,LQi4) = LQi2

       ! Collisional friction with bulk ions

       ELM(1:4, 5,LQi4,NE) = - (AME / AMI) * fem_int(2,NE,rNueiEI)
       NLC( 5,LQi4) = LQi4

       ELM(1:4, 6,LQi4,NE) =   (AME / AMI) * fem_int(2,NE,rNuei)
       NLC( 6,LQi4) = LQe4

       ! Collisional friction with beam ions

       ELM(1:4, 7,LQi4,NE) = - (AMB / AMI) * fem_int(2,NE,rNubiBI)
       NLC( 7,LQi4) = LQi4

       ELM(1:4, 8,LQi4,NE) =   (AMB / AMI) * fem_int(2,NE,rNubi)
       NLC( 8,LQi4) = LQb4

       IF(MDLWTB == 0) THEN

          ! Wave interaction force (electron driven)

          ELM(1:4, 9,LQi4,NE) = - 1.D0 / AMI * fem_int(2,NE,FWpheBB)
          NLC( 9,LQi4) = LQe3

          ELM(1:4,10,LQi4,NE) =   1.D0 / AMI * fem_int(2,NE,FWpheBB2)
          NLC(10,LQi4) = LQe4

          ELM(1:4,11,LQi4,NE) =   1.D0 / AMI * fem_int(44,NE,FWpheBB,WPM)
          NLC(11,LQi4) = LQe1

          ! Wave interaction force (ion driven)

          ELM(1:4,12,LQi4,NE) =   1.D0 / AMI * fem_int(2,NE,FWphiBB)
          NLC(12,LQi4) = LQi3

          ELM(1:4,13,LQi4,NE) = - 1.D0 / AMI * fem_int(2,NE,FWphiBB2)
          NLC(13,LQi4) = LQi4

          ELM(1:4,14,LQi4,NE) = - 1.D0 / AMI * fem_int(44,NE,FWphiBB,WPM)
          NLC(14,LQi4) = LQi1

       END IF

       ! Loss to divertor

       ELM(1:4,15,LQi4,NE) = - fem_int(2,NE,rNuL)
       NLC(15,LQi4) = LQi4

       ! Collisional friction force with neutrals

       ELM(1:4,16,LQi4,NE) = - fem_int(2,NE,rNu0i)
       NLC(16,LQi4) = LQi4

       ! Charge exchange force

       ELM(1:4,17,LQi4,NE) = - fem_int(2,NE,rNuiCX)
       NLC(17,LQi4) = LQi4

       ! Loss cone loss

       PELM(1:4,18,LQi4,NE) =   fem_int(-1,NE,SiLCph)
       NLC(18,LQi4) = 0

!       ! Helical Neoclassical viscosity force
!
!       ELM(1:4,19,LQi4,NE) = UHth * UHph / 2.D0 * fem_int(2,NE,rNuiHL)
!       NLC(19,LQi4) = LQi3
!
!       ELM(1:4,20,LQi4,NE) = - (1.D0 - UHph * UHph) * fem_int(2,NE,rNuiHL)
!       NLC(20,LQi4) = LQi4
    END DO

    NLCMAX(LQi4) = 18
    RETURN
  END SUBROUTINE LQi4CC

!***************************************************************
!
!  Ion Energy Transport: Ti
!
!***************************************************************

  SUBROUTINE LQi5CC

    use physical_constants, only : AEE, rKeV

    INTEGER :: NE
    REAL(8) :: Veff, Deff, Pec, P

    ! Temperature evolution

    IF(MDFIXT == 0) THEN
       DO NE = 1, NEMAX
          ELM(1:4, 0,LQi5,NE) =   1.5D0 / DT * lump_int(1,NE)
          NLC( 0,LQi5) = LQi5

          ! Convection transport

          ELM(1:4, 1,LQi5,NE) = - 5.D0 * fem_int(3,NE,RUirV)
          NLC( 1,LQi5) = LQi5

          ! Conduction transport

          ELM(1:4, 2,LQi5,NE) = - 6.D0 * fem_int(18,NE,Chii) &
               &                - 6.D0 * fem_int(39,NE,ChiiNi,dPNiV)
          NLC( 2,LQi5) = LQi5

          ! Joule heating

          ELM(1:4, 3,LQi5,NE) =   PZ * AEE / rKeV * fem_int(2,NE,EthVR)
          NLC( 3,LQi5) = LQi3

          ELM(1:4, 4,LQi5,NE) =   PZ * AEE / rKeV * fem_int(2,NE,EphV)
          NLC( 4,LQi5) = LQi4

          ! Collisional transfer with electrons

          ELM(1:4, 5,LQi5,NE) = - fem_int(2,NE,rNuTeiEI)
          NLC( 5,LQi5) = LQi5

          ELM(1:4, 6,LQi5,NE) =   fem_int(2,NE,rNuTei)
          NLC( 6,LQi5) = LQe5

          ! Collisional heating with beam

          ELM(1:4, 7,LQi5,NE) = - 0.5D0 * AMb * (PNBCD * Vb) / rKeV &
               &                * fem_int(2,NE,rNubiBI)
          NLC( 7,LQi5) = LQi4

          ELM(1:4, 8,LQi5,NE) =   0.5D0 * AMb * (PNBCD * Vb) / rKeV &
               &                * fem_int(2,NE,rNubi)
          NLC( 8,LQi5) = LQb4

          ! Loss to diverter

          ELM(1:4, 9,LQi5,NE) = -                       fem_int( 2,NE,rNuLTi)
          NLC( 9,LQi5) = LQi5

          PELM(1:4,10,LQi5,NE) =          PNeDIV / PZ * fem_int(-2,NE,rNuLTi,PTiV)
          NLC(10,LQi5) = 0

          ELM(1:4,11,LQi5,NE) = - 1.5D0               * fem_int( 2,NE,rNuLTi)
          NLC(11,LQi5) = LQi5

          PELM(1:4,12,LQi5,NE) =  1.5D0 * PTiDIV      * fem_int(-2,NE,rNuLTi,PNiV)
          NLC(12,LQi5) = 0

          ! Direct heating (RF)

          PELM(1:4,13,LQi5,NE) = 1.D0 / (1.D20 * rKeV) * fem_int(-1,NE,PRFi)
          NLC(13,LQi5) = 0

          ! *** Streamline Upwind Petrov Galerkin Method ***
          
!!$          Veff = 5.D0 * (0.5D0 * (UirVR (NE-1) + UirVR (NE))) &
!!$               & + 6.D0 * (0.5D0 * (ChiiNi(NE-1) + ChiiNi(NE)))  &
!!$!               &         * 0.5D0 * (dPNiV (NE-1) + dPNiV (NE)))
!!$               &         * (PNiV(NE) - PNiV(NE-1)) / HPSI(NE)
!!$          Deff = 6.D0 * (0.5D0 * (Chii  (NE-1) + Chii  (NE)))
!!$          Pec  = 0.5D0 * Veff * HPSI(NE) / Deff
!!$          P    = 0.5D0 * HPSI(NE) * falpha(Pec)
!!$
!!$          ELM(1:4, 0,LQi5,NE) = ELM(1:4, 0,LQi5,NE) + 1.5D0 / DT * fem_int(8,NE) * P
!!$
!!$          ELM(1:4, 1,LQi5,NE) = ELM(1:4, 1,LQi5,NE) - 5.D0 * fem_int(10,NE,RUirV) * P
!!$
!!$!          ELM(1:4, 2,LQi5,NE) = ELM(1:4, 2,LQi5,NE) - 6.D0 * fem_int(42,NE,ChiiNi,dPNiV) * P
!!$          ELM(1:4, 2,LQi5,NE) = ELM(1:4, 2,LQi5,NE) - 6.D0 * fem_int(46,NE,ChiiNi,PNiV) * P
!!$
!!$          ELM(1:4, 3,LQi5,NE) = ELM(1:4, 3,LQi5,NE) + PZ * AEE / rKeV * fem_int(9,NE,EthVR) * P
!!$
!!$          ELM(1:4, 4,LQi5,NE) = ELM(1:4, 4,LQi5,NE) + PZ * AEE / rKeV * fem_int(9,NE,EphV) * P
!!$
!!$          ELM(1:4, 5,LQi5,NE) = ELM(1:4, 5,LQi5,NE) - fem_int(9,NE,rNuTeiEI) * P
!!$
!!$          ELM(1:4, 6,LQi5,NE) = ELM(1:4, 6,LQi5,NE) + fem_int(9,NE,rNuTei) * P
!!$
!!$          ELM(1:4, 7,LQi5,NE) = ELM(1:4, 7,LQi5,NE) - 0.5D0 * AMb * (PNBCD * Vb) / rKeV &
!!$               &                * fem_int(9,NE,rNubiBI) * P
!!$
!!$          ELM(1:4, 8,LQi5,NE) = ELM(1:4, 8,LQi5,NE) + 0.5D0 * AMb * (PNBCD * Vb) / rKeV &
!!$               &                * fem_int(9,NE,rNubi) * P
!!$
!!$          ELM(1:4, 9,LQi5,NE) = ELM(1:4, 9,LQi5,NE) - fem_int(9,NE,rNuLTi) * P
!!$
!!$          PELM(1:4,10,LQi5,NE) = PELM(1:4,10,LQi5,NE) &
!!$               &               + PNeDIV / PZ * fem_int(-9,NE,rNuLTi,PTiV) * P
!!$
!!$          ELM(1:4,11,LQi5,NE) = ELM(1:4,11,LQi5,NE) &
!!$               &              - 1.5D0 * fem_int(9,NE,rNuLTi) * P
!!$
!!$          PELM(1:4,12,LQi5,NE) = PELM(1:4,12,LQi5,NE) &
!!$               &               + 1.5D0 * PTiDIV * fem_int(-9,NE,rNuLTi,PNiV) * P
!!$
!!$          PELM(1:4,13,LQi5,NE) = PELM(1:4,11,LQi5,NE) &
!!$               &               + 1.D0 / (1.D20 * rKeV) * fem_int(-8,NE,PRFi) * P
       END DO

       NLCMAX(LQi5) = 13
    ELSE

       !  Fixed temperature profile

       DO NE = 1, NEMAX
          ELM(1:4,0,LQi5,NE) = 1.D0 / DT * lump_int(1,NE)
          NLC(0,LQi5) = LQi5

          ! Convection

          ELM(1:4,1,LQi5,NE) = - 2.D0 * fem_int(5,NE,PTiV)
          NLC(1,LQi5) = LQi2

          ! Ionization of n01 and n02

          ELM(1:4,2,LQi5,NE) =     1.D0 / PZ * fem_int(28,NE,rNuIN0,PTiV)
          NLC(2,LQi5) = LQn1

          ELM(1:4,3,LQi5,NE) =     1.D0 / PZ * fem_int(28,NE,rNuIN0,PTiV)
          NLC(3,LQi5) = LQn2

          ! Loss to divertor

          ELM(1:4,4,LQi5,NE) = - fem_int(2,NE,rNuL)
          NLC(4,LQi5) = LQi5

          PELM(1:4,5,LQi5,NE) =  PNeDIV / PZ * fem_int(-2,NE,rNuL,PTiV)
          NLC(5,LQi5) = 0

          ! Particle source from beam ion

          ELM(1:4,6,LQi5,NE) =   fem_int(28,NE,rNuB,PTiV)
          NLC(6,LQi5) = LQb1

          ! NBI kick up ions

          PELM(1:4,7,LQi5,NE) = - fem_int(-2,NE,SNB,PTiV)
          NLC(7,LQi5) = 0

          ! Loss cone loss

          PELM(1:4,8,LQi5,NE) =   fem_int(-2,NE,SiLC,PTiV)
          NLC(8,LQi5) = 0
       END DO

       NLCMAX(LQi5) = 8
    END IF

    RETURN
  END SUBROUTINE LQi5CC

!***************************************************************
!
!   Beam Ion Density
!
!***************************************************************

  SUBROUTINE LQb1CC

    INTEGER :: NE

    DO NE = 1, NEMAX
       ELM(1:4,0,LQb1,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC(0,LQb1) = LQb1

       ! NBI particle source

       PELM(1:4,1,LQb1,NE) =   fem_int(-1,NE,SNB)
       NLC(1,LQb1) = 0

       ! Relaxation to thermal ions

       ELM(1:4,2,LQb1,NE) = - fem_int(2,NE,rNuB)
       NLC(2,LQb1) = LQb1
    END DO

    NLCMAX(LQb1) = 2
    RETURN
  END SUBROUTINE LQb1CC

!***************************************************************
!
!   Beam Ion Poloidal Flow
!
!***************************************************************

  SUBROUTINE LQb3CC

    use physical_constants, only : AEE, AME

    INTEGER :: NE

    ! Ubth(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4,0,LQb3,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC(0,LQb3) = LQb3

       ! Poroidal E force

       ELM(1:4,1,LQb3,NE) = - (PZ * AEE / AMB) * fem_int(2,NE,PNbV)
       NLC(1,LQb3) = LQm2

       ! Neoclassical viscosity force

       ELM(1:4,2,LQb3,NE) = - fem_int(2,NE,rNuiNC)
       NLC(2,LQb3) = LQb3

       ! Collisional friction force with electrons

       ELM(1:4,3,LQb3,NE) = - (AME / AMB) * fem_int(2,NE,rNube)
       NLC(3,LQb3) = LQb3

       ELM(1:4,4,LQb3,NE) =   (AME / AMB) * fem_int(2,NE,rNubeBE)
       NLC(4,LQb3) = LQe3

       ! Collisional friction force with ions

       ELM(1:4,5,LQb3,NE) = - (AMI / AMB) * fem_int(2,NE,rNubi)
       NLC(5,LQb3) = LQb3

       ELM(1:4,6,LQb3,NE) =   (AMI / AMB) * fem_int(2,NE,rNubiBI)
       NLC(6,LQb3) = LQi3

       ! Collisional friction force with neutrals

       ELM(1:4,7,LQb3,NE) = - fem_int(2,NE,rNu0b)
       NLC(7,LQb3) = LQb3

       ! Charge exchange force

       ELM(1:4,8,LQb3,NE) = - fem_int(2,NE,rNuiCX)
       NLC(8,LQb3) = LQb3
    END DO

    ! Ubth(NRMAX) : 0

    NLCMAX(LQb3) = 8
    RETURN
  END SUBROUTINE LQb3CC

!***************************************************************
!
!   Beam Ion Toroidal Flow
!
!***************************************************************
 
 SUBROUTINE LQb4CC

    INTEGER :: NE

    ! - UbPhi(0)' : 0

    DO NE = 1, NEMAX
       ELM(1:4,0,LQb4,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC(0,LQb4) = LQb4

       ! Collisional friction with electrons

       ELM(1:4,1,LQb4,NE) = - fem_int(2,NE,rNube)
       NLC(1,LQb4) = LQb4

       ELM(1:4,2,LQb4,NE) =   fem_int(2,NE,rNubeBE)
       NLC(2,LQb4) = LQe4

       ! Collisional friction with ions

       ELM(1:4,3,LQb4,NE) = - fem_int(2,NE,rNubi)
       NLC(3,LQb4) = LQb4

       ELM(1:4,4,LQb4,NE) =   fem_int(2,NE,rNubiBI)
       NLC(4,LQb4) = LQi4

       ! Collisional friction force with neutrals

       ELM(1:4,5,LQb4,NE) = - fem_int(2,NE,rNu0b)
       NLC(5,LQb4) = LQb4

       ! Charge exchange force

       ELM(1:4,6,LQb4,NE) = - fem_int(2,NE,rNuiCX)
       NLC(6,LQb4) = LQb4

       ! NBI momentum source

       PELM(1:4,7,LQb4,NE) = (PNBCD * Vb) * fem_int(-1,NE,SNB)
       NLC(7,LQb4) = 0
    END DO

    ! Ubphi(NRMAX) : 0

    NLCMAX(LQb4) = 7
    RETURN
  END SUBROUTINE LQb4CC

!***************************************************************
!
!   Slow Neutral Transport: n01
!
!***************************************************************

  SUBROUTINE LQn1CC

    INTEGER :: NE

    DO NE = 1, NEMAX
       ELM(1:4,0,LQn1,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC(0,LQn1) = LQn1

       !  Diffusion of neutrals

       ELM(1:4,1,LQn1,NE) = - 4.D0 * fem_int(18,NE,D01)
       NLC(1,LQn1) = LQn1

       ! Ionization

       ELM(1:4,2,LQn1,NE) = - 1.D0 / PZ * fem_int(2,NE,rNuIN0)
       NLC(2,LQn1) = LQn1

       ! Generation of fast neutrals by charge exchange

       ELM(1:4,3,LQn1,NE) = - fem_int(2,NE,rNuCXN0)
       NLC(3,LQn1) = LQn1

       ! Recycling from divertor

       ELM(1:4,4,LQn1,NE) =   rGamm0 / PZ * fem_int(2,NE,rNuL)
       NLC(4,LQn1) = LQe1

       PELM(1:4,5,LQn1,NE) = - rGamm0 * PNeDIV / PZ * fem_int(-1,NE,rNuL)
       NLC(5,LQn1) = 0
    END DO

    NLCMAX(LQn1) = 5
    RETURN
  END SUBROUTINE LQn1CC

!***************************************************************
!
!   Fast Neutral Transport: n02
!
!***************************************************************

  SUBROUTINE LQn2CC

    INTEGER :: NE

    DO NE = 1, NEMAX
       ELM(1:4,0,LQn2,NE) = 1.D0 / DT * lump_int(1,NE)
       NLC(0,LQn2) = LQn2

       !  Diffusion of neutrals

       ELM(1:4,1,LQn2,NE) = - 4.D0 * fem_int(18,NE,D02)
       NLC(1,LQn2) = LQn2

       ! Ionization

       ELM(1:4,2,LQn2,NE) = - 1.D0 / PZ * fem_int(2,NE,rNuIN0)
       NLC(2,LQn2) = LQn2

       ! Generation of fast neutrals by charge exchange

       ELM(1:4,3,LQn2,NE) = fem_int(2,NE,rNuCXN0)
       NLC(3,LQn2) = LQn1

       ! NBI particle source

       PELM(1:4,4,LQn2,NE) = fem_int(-1,NE,SNB)
       NLC(4,LQn2) = 0
    END DO

    NLCMAX(LQn2) = 4
    RETURN
  END SUBROUTINE LQn2CC

!***************************************************************
!
!   Dirichlet condition
!
!***************************************************************

  SUBROUTINE BOUNDARY(NR,LQ,ID,VAL)

    integer, intent(in) :: NR, LQ, ID
    real(8), intent(in), optional :: VAL
    integer :: NQ

    IF(ID == 0) THEN
       ! Initialize ALC, BLC and CLC
       ALC(0:NCM,LQ,NR) = 0.D0
       BLC(0:NCM,LQ,NR) = 0.D0
       CLC(0:NCM,LQ,NR) = 0.D0
       PLC(1:NCM,LQ,NR) = 0.D0
       ! Diagonal term on handled variable
       !    Although I have no precise proof, I need to set BLC to the value 
       !    comparable to ALC, BLC and CLC on the other grid if X=0 on axis is imposed.
       BLC(1,LQ,NR) = 1.D10
       NLCR(1,LQ,NR) = LQ

       IF(PRESENT(VAL)) THEN
          BLC(1,LQ,NR) = 1.D0
          PLC(2,LQ,NR) = - VAL
          NLCR(2,LQ,NR) = 0
          IF(NLCMAX(LQ) <= 1) NLCMAX(LQ) = 2
       END IF
    ELSE
       NLCMAX(LQ) = NLCMAX(LQ) + 1
       PLC(NLCMAX(LQ),LQ,NR) = VAL
       NLCR(NLCMAX(LQ),LQ,NR) = 0
    END IF

  END SUBROUTINE BOUNDARY

!***************************************************************
!
!   Langevin function for SUPG method
!
!***************************************************************

  real(8) function falpha(x)
    real(8), intent(in) :: x

    if (x < -3.d0) then
       falpha = - 1.d0 - 1.d0 / x
    elseif (x > 3.d0) then
       falpha =   1.d0 - 1.d0 / x
    else
       falpha = x / 3.d0 * ( 1.d0 - abs(x) / 9.d0)
    end if

  end function falpha

end module coefficients
