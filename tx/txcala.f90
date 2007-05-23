!     $Id$
module coefficients
  use commons
  use core_module
  implicit none
  private
  real(8), dimension(:,:,:,:), allocatable :: ELM, PELM, test
  real(8), dimension(0:NRM) :: rNuIN0, rNuCXN0, rNueiEI, rNubeBE, rNubiBI, rNuTeiEI,&
       &                       ChieNe, ChiiNi, rMueNe, rMuiNi, dPNeV, dPNiV, &
       &                       UethVR, UithVR, EthVR, RUerV, RUirV, UerVR, UirVR, &
       &                       FWpheBB, FWphiBB, dAphV, FWpheBB2, FWphiBB2
  real(8), dimension(0:NRM) :: UNITY = 1.D0
  real(8) :: DTt, DTf(1:NQM), invDT
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

    allocate(ELM(1:NEMAX,1:4,0:NCM,1:NQMAX))
    allocate(PELM(1:NEMAX,1:4,1:NCM,1:NQMAX))

    ELM(1:NEMAX,1:4,0:NCM,1:NQMAX) = 0.D0
    PELM(1:NEMAX,1:4,1:NCM,1:NQMAX) = 0.D0

    ALC(0:NCM,1:NQMAX,0:NRMAX) = 0.D0
    BLC(0:NCM,1:NQMAX,0:NRMAX) = 0.D0
    CLC(0:NCM,1:NQMAX,0:NRMAX) = 0.D0
    NLCR(0:NCM,1:NQMAX,0:NRMAX) = 1
    NLC(0:NCM,1:NQMAX) = 1
    PLC(1:NCM,1:NQMAX,0:NRMAX) = 0.D0
    NLCMAX(1:NQMAX) = 0

    invDT = 1.d0 / DT
!    IF(DT <= 3.D-6) THEN
    IF(DT <= 1.D-5) THEN
       DTt          = 1.d2
       DTf(1)       = 1.d0
       DTf(2:NQMAX) = 1.d2
    ELSE IF(DT <= 1.D-4) THEN
       DTt          = 1.d3
       DTf(1)       = 1.d0
       DTf(2:NQMAX) = 1.d3
    ELSE
       DTt          = 1.d0
       DTf(1:NQMAX) = 1.d0
    END IF

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

    CALL BOUNDARY(NRMAX,LQm2,1, 2.D0*PSI(NRMAX)*BB/sqeps0)
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
    CALL DERIVS(PSI,X,LQe1,NQMAX,NRMAX,dPNeV)
    CALL DERIVS(PSI,X,LQi1,NQMAX,NRMAX,dPNiV)
    CALL DERIVS(PSI,X,LQm4,NQMAX,NRMAX,dAphV)
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

    use physical_constants, only : AEE

    ! phi'(0) : 0

    ELM(1:NEMAX,1:4,1,LQm1) =   4.D0 * sqeps0 / (AEE * 1.D20) * fem_int(18,UNITY)
    NLC(1,LQm1) = LQm1

    ELM(1:NEMAX,1:4,2,LQm1) =        fem_int(1) / sqeps0
    NLC(2,LQm1) = LQe1

    ELM(1:NEMAX,1:4,3,LQm1) = - PZ * fem_int(1) / sqeps0
    NLC(3,LQm1) = LQi1

    ELM(1:NEMAX,1:4,4,LQm1) = - PZ * fem_int(1) / sqeps0
    NLC(4,LQm1) = LQb1

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

    use physical_constants, only : AEE, rMU0

    ! (r*Atheta)'(0) : 0

    ELM(1:NEMAX,1:4,0,LQm2) = sqeps0 * fem_int(1) * invDT
    NLC(0,LQm2) = LQm2

    ! rot Bphi

    ELM(1:NEMAX,1:4,1,LQm2) =(- 4.D0 * fem_int(18,UNITY) - 4.D0 * fem_int(4))
    NLC(1,LQm2) = LQm5

    ! Electron current

    ELM(1:NEMAX,1:4,2,LQm2) = - rMU0      * AEE * 1.D20 * fem_int(1) / sqeps0
    NLC(2,LQm2) = LQe3

    ! Ion current

    ELM(1:NEMAX,1:4,3,LQm2) =   rMU0 * PZ * AEE * 1.D20 * fem_int(1) / sqeps0
    NLC(3,LQm2) = LQi3

    ! Beam ion current

    ELM(1:NEMAX,1:4,4,LQm2) =   rMU0 * PZ * AEE * 1.D20 * fem_int(1) / sqeps0
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

    use physical_constants, only : AEE, VC, rMU0

    ! Aphi'(0) : 0

    ELM(1:NEMAX,1:4,0,LQm3) = fem_int(1) / (VC**2 * DT)
    NLC(0,LQm3) = LQm3

    ! rot Btheta

    ELM(1:NEMAX,1:4,1,LQm3) = - 4.D0 * fem_int(18,UNITY)
    NLC(1,LQm3) = LQm4

    ! Electron current

    ELM(1:NEMAX,1:4,2,LQm3) = - rMU0 *      AEE * 1.D20 * fem_int(1)
    NLC(2,LQm3) = LQe4

    ! Ion current

    ELM(1:NEMAX,1:4,3,LQm3) =   rMU0 * PZ * AEE * 1.D20 * fem_int(1)
    NLC(3,LQm3) = LQi4

    ! Beam ion current

    ELM(1:NEMAX,1:4,4,LQm3) =   rMU0 * PZ * AEE * 1.D20 * fem_int(1)
    NLC(4,LQm3) = LQb4

!   ! Virtual current for helical system
!
!   PELM(1:NEMAX,1:4,5,LQm3) =   rMU0 * fem_int(-1,AJV)
!   NLC(5,LQm3) = 0

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

    ELM(1:NEMAX,1:4,1,LQm4) = fem_int(1)
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

    use physical_constants, only : rMU0

    ELM(1:NEMAX,1:4,0,LQm5) = fem_int(1) * invDT 
    NLC(0,LQm5) = LQm5

    ! r * Atheta'

    ELM(1:NEMAX,1:4,1,LQm5) = fem_int(1) / (rMU0 * sqeps0)
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

    ! Convection

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

    ! Ns*Usr(0) : fixed

    ELM(1:NEMAX,1:4,0,LQe2) = fem_int(1) * invDT &
         &                  + fem_int(8) * fem_int(0) * invDT
    NLC(0,LQe2) = LQe2

    ! Nonlinear term
    
    ELM(1:NEMAX,1:4,1,LQe2) = - 2.D0 * fem_int(3,RUerV) + fem_int(2,UerVR) &
         &                  +(- 2.D0 * fem_int(10,RUerV) + fem_int(9,UerVR)) * fem_int(0)
    NLC(1,LQe2) = LQe2

    ! Nonlinear centrifugal force

    ELM(1:NEMAX,1:4,2,LQe2) = fem_int(2,UethVR) &
         &                  + fem_int( 9,UethVR) * fem_int(0)
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

    ELM(1:NEMAX,1:4,6,LQe2) = - 2.D0 * (AEE / AME) * fem_int(16,AphV) &
         &                    - 2.D0 * (AEE / AME) * fem_int(20,UNITY,AphV) * fem_int(0)
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

    use physical_constants, only : AEE, AME, rKeV, rMU0

    integer :: N

    ! Ns*UsTheta(0) : 0

    ELM(1:NEMAX,1:4, 0,LQe3) = fem_int(1) * invDT
    NLC( 0,LQe3) = LQe3

    ! Nonlinear term

    ELM(1:NEMAX,1:4, 1,LQe3) = - 2.D0 * fem_int(3,RUerV)
    NLC( 1,LQe3) = LQe3

    ! Viscosity force

    ELM(1:NEMAX,1:4, 2,LQe3) = - 4.D0 * fem_int(18,rMue) &
         &                     - 4.D0 * fem_int(39,rMueNe,dPNeV) &
         &                     - 4.D0 * fem_int( 3,rMue)
    NLC( 2,LQe3) = LQe3

    ! Poloidal E force

    ELM(1:NEMAX,1:4, 3,LQe3) =   (AEE / AME) * fem_int(2,PNeV) / rMU0
    NLC( 3,LQe3) = LQm2

    ! v x B force

    ELM(1:NEMAX,1:4, 4,LQe3) =   (AEE / AME) * fem_int(2,BphV)
    NLC( 4,LQe3) = LQe2

    ! Neoclassical viscosity force

    ELM(1:NEMAX,1:4, 5,LQe3) = - fem_int(2,rNueNC)
    NLC( 5,LQe3) = LQe3

    ! Collisional friction force with ions

    ELM(1:NEMAX,1:4, 6,LQe3) = - fem_int(2,rNuei)
    NLC( 6,LQe3) = LQe3

    ELM(1:NEMAX,1:4, 7,LQe3) =   fem_int(2,rNueiEI)
    NLC( 7,LQe3) = LQi3

    ! Collisional friction with beam ions

    ELM(1:NEMAX,1:4, 8,LQe3) = - (AMB / AME) * fem_int(2,rNubeBE)
    NLC( 8,LQe3) = LQe3

    ELM(1:NEMAX,1:4, 9,LQe3) =   (AMB / AME) * fem_int(2,rNube)
    NLC( 9,LQe3) = LQb3

    IF(MDLWTB == 0) THEN

       ! Wave int_arrayeraction force (electron driven)

       ELM(1:NEMAX,1:4,10,LQe3) = - 1.D0 / AME * fem_int(2,FWthe)
       NLC(10,LQe3) = LQe3

       ELM(1:NEMAX,1:4,11,LQe3) = - 2.D0 / AME * fem_int(36,AphV,FWthphe)
       NLC(11,LQe3) = LQe4

       ELM(1:NEMAX,1:4,12,LQe3) =   1.D0 / AME * fem_int(44,FWthe,WPM)
       NLC(12,LQe3) = LQe1

       ! Wave int_arrayeraction force (ion driven)

       ELM(1:NEMAX,1:4,13,LQe3) =   1.D0 / AME * fem_int(2,FWthi)
       NLC(13,LQe3) = LQi3

       ELM(1:NEMAX,1:4,14,LQe3) =   2.D0 / AME * fem_int(36,AphV,FWthphi)
       NLC(14,LQe3) = LQi4

       ELM(1:NEMAX,1:4,15,LQe3) = - 1.D0 / AME * fem_int(44,FWthi,WPM)
       NLC(15,LQe3) = LQi1

       N = 0
    ELSEIF(MDLWTB == 1) THEN

       ! Wave interaction force (electron driven)

       ELM(1:NEMAX,1:4,10,LQe3) = - 2.D0 / AME * fem_int(17,WNthe)
       NLC(10,LQe3) = LQe1

!      ELM(1:NEMAX,1:4,11,LQe3) = - 2.D0 / AME * fem_int(36,Phi,WEMthe)
!      NLC(1:NEMAX,1:4,11,LQe3) = LQe1
!
!      ELM(1:NEMAX,1:4,12,LQe3) = - 2.D0 / AME * fem_int(15,WWthe)
!      NLC(12,LQe3) = LQe1
!
!      ELM(1:NEMAX,1:4,13,LQe3) = - 2.D0 / AME * rKeV * fem_int(17,WT1the)
!      NLC(13,LQe3) = LQe5
!
!      ELM(1:NEMAX,1:4,14,LQe3) =   2.D0 / AME * fem_int(17,WT2the)
!      NLC(14,LQe3) = LQe1

       ! Wave interaction force (ion driven)

       ELM(1:NEMAX,1:4,15,LQe3) =   2.D0 / AME * fem_int(17,WNthi)
       NLC(15,LQe3) = LQi1

!      ELM(1:NEMAX,1:4,16,LQe3) =   2.D0 / AME * fem_int(36,Phi,WEMthi)
!      NLC(16,LQe3) = LQi1
!
!      ELM(1:NEMAX,1:4,17,LQe3) =   1.D0 / AME * fem_int(15,WWthi)
!      NLC(17,LQe3) = LQi1
!
!      ELM(1:NEMAX,1:4,18,LQe3) =   2.D0 / AME * rKeV * fem_int(17,WT1thi)
!      NLC(18,LQe3) = LQi5
!
!      ELM(1:NEMAX,1:4,19,LQe3) = - 2.D0 / AME * fem_int(17,WT2thi)
!      NLC(19,LQe3) = LQi1

       N = 4
    END IF

    ! Loss to divertor

    ELM(1:NEMAX,1:4,16+N,LQe3) = - 2.D0 * fem_int(2,rNuL)
    NLC(16+N,LQe3) = LQe3

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,17+N,LQe3) = - fem_int(2,rNu0e)
    NLC(17+N,LQe3) = LQe3

!   ! Helical neoclassical viscosity force
!
!   ELM(1:NEMAX,1:4,18+N,LQe3) = - (1.D0 - UHth * UHth) * fem_int(2,rNueHL)
!   NLC(18+N,LQe3) = LQe3
!
!   ELM(1:NEMAX,1:4,19+N,LQe3) = UHph * UHth / 2.D0 * fem_int(22,rNueHL)
!   NLC(19+N,LQe3) = LQe4

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

    ! Uephi(0)' : 0

    ELM(1:NEMAX,1:4, 0,LQe4) = fem_int(1) * invDT
    NLC( 0,LQe4) = LQe4

    ! Nonlinear term

    ELM(1:NEMAX,1:4, 1,LQe4) = - 2.D0 * fem_int(3,RUerV)
    NLC( 1,LQe4) = LQe4
    
    ! Viscosity force
    
    ELM(1:NEMAX,1:4, 2,LQe4) = - 4.D0 *(fem_int(18,rMue) + fem_int(39,rMueNe,dPNeV))
    NLC( 2,LQe4) = LQe4

    ! Toroidal E force

    ELM(1:NEMAX,1:4, 3,LQe4) =   (AEE / AME) * fem_int(2,PNeV)
    NLC( 3,LQe4) = LQm3

    ! v x B force

    ELM(1:NEMAX,1:4, 4,LQe4) =   2.D0 * (AEE / AME) * fem_int(6,AphV)
    NLC( 4,LQe4) = LQe2

    ! Collisional friction with bulk ions

    ELM(1:NEMAX,1:4, 5,LQe4) = - fem_int(2,rNuei)
    NLC( 5,LQe4) = LQe4

    ELM(1:NEMAX,1:4, 6,LQe4) =   fem_int(2,rNueiEI)
    NLC( 6,LQe4) = LQi4

    ! Collisional friction with beam ions

    ELM(1:NEMAX,1:4, 7,LQe4) = - (AMB / AME) * fem_int(2,rNubeBE)
    NLC( 7,LQe4) = LQe4

    ELM(1:NEMAX,1:4, 8,LQe4) =   (AMB / AME) * fem_int(2,rNube)
    NLC( 8,LQe4) = LQb4

    IF(MDLWTB == 0) THEN

       ! Wave int_arrayeraction force (electron driven)

       ELM(1:NEMAX,1:4, 9,LQe4) =   1.D0 / AME * fem_int(2,FWpheBB)
       NLC( 9,LQe4) = LQe3

       ELM(1:NEMAX,1:4,10,LQe4) = - 1.D0 / AME * fem_int(2,FWpheBB2)
       NLC(10,LQe4) = LQe4

       ELM(1:NEMAX,1:4,11,LQe4) = - 1.D0 / AME * fem_int(44,FWpheBB,WPM)
       NLC(11,LQe4) = LQe1

       ! Wave int_arrayeraction force (ion driven)

       ELM(1:NEMAX,1:4,12,LQe4) = - 1.D0 / AME * fem_int(2,FWphiBB)
       NLC(12,LQe4) = LQi3

       ELM(1:NEMAX,1:4,13,LQe4) =   1.D0 / AME * fem_int(2,FWphiBB2)
       NLC(13,LQe4) = LQi4

       ELM(1:NEMAX,1:4,14,LQe4) =   1.D0 / AME * fem_int(44,FWphiBB,WPM)
       NLC(14,LQe4) = LQi1

    END IF

    ! Loss to divertor

    ELM(1:NEMAX,1:4,15,LQe4) = - 2.D0 * fem_int(2,rNuL)
    NLC(15,LQe4) = LQe4

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,16,LQe4) = - fem_int(2,rNu0e)
    NLC(16,LQe4) = LQe4

!    ! Helical neoclassical viscosity force

!    ELM(1:NEMAX,1:4,17,LQe4) =  UHth * UHph / 2.D0 * fem_int(2,rNueHL)
!    NLC(17,LQe4) = LQe3

!    ELM(1:NEMAX,1:4,18,LQe4) = - (1.D0 - UHph * UHph) * fem_int(2,rNueHL)
!    NLC(18,LQe4) = LQe4

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

    ! Temperature evolution
    
    IF(MDFIXT == 0) THEN
       ELM(1:NEMAX,1:4, 0,LQe5) =   1.5D0 * fem_int(1) * invDT
       NLC( 0,LQe5) = LQe5

       ! Convection transport

       ELM(1:NEMAX,1:4, 1,LQe5) = - 5.D0 * fem_int(3,RUerV)
       NLC( 1,LQe5) = LQe5

       ! Conduction transport

!      ELM(1:NEMAX,1:4, 2,LQe5) = - 6.D0 *(fem_int(18,Chie) + 6.D0 * fem_int(39,ChieNe,dPNeV))
!      NLC( 2,LQe5) = LQe5
       ELM(1:NEMAX,1:4, 2,LQe5) = - 6.D0 * fem_int(18,Chie)
       NLC( 2,LQe5) = LQe5

       ELM(1:NEMAX,1:4,15,LQe5) =   6.D0 * fem_int(41,Chie,PTeV)
       NLC(15,LQe5) = LQe1

       ! Joule heating

       ELM(1:NEMAX,1:4, 3,LQe5) = - AEE / rKeV * fem_int(2,EthVR)
       NLC( 3,LQe5) = LQe3

       ELM(1:NEMAX,1:4, 4,LQe5) = - AEE / rKeV * fem_int(2,EphV)
       NLC( 4,LQe5) = LQe4

       ! Collisional transfer with ions

       ELM(1:NEMAX,1:4, 5,LQe5) = - 1.5d0 * fem_int(2,rNuTei)
       NLC( 5,LQe5) = LQe5

       ELM(1:NEMAX,1:4, 6,LQe5) =   1.5d0 * fem_int(2,rNuTeiEI)
       NLC( 6,LQe5) = LQi5

       ! Collisional heating with beam

!      ELM(1:NEMAX,1:4, 7,LQe5) = - 0.5D0 * AMb * (PNBCD * Vb) / rKeV &
!           &                * fem_int(2,rNubeBE)
!      NLC( 7,LQe5) = LQe4
!
!      ELM(1:NEMAX,1:4, 8,LQe5) =   0.5D0 * AMb * (PNBCD * Vb) / rKeV &
!           &                * fem_int(2,rNube)
!      NLC( 8,LQe5) = LQb4

       ELM(1:NEMAX,1:4, 7,LQe5) = - 0.5D0 * AMb / rKeV * fem_int(28,Vbedir,rNubeBE)
       NLC( 7,LQe5) = LQe4

       ELM(1:NEMAX,1:4, 8,LQe5) =   0.5D0 * AMb / rKeV * fem_int(28,Vbedir,rNube)
       NLC( 8,LQe5) = LQb4

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

       ! Radiation loss

       PELM(1:NEMAX,1:4,14,LQe5) = - 1.D0 / (1.D20 * rKeV) * fem_int(-1,PBr)
       NLC(14,LQe5) = 0

!       NLCMAX(LQe5) = 14
       NLCMAX(LQe5) = 15
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

    ! Convection

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

    ! NBI kick up ions

    PELM(1:NEMAX,1:4,7,LQi1) = - fem_int(-1,SNB)
    NLC(7,LQi1) = 0

    ! Loss cone loss

    PELM(1:NEMAX,1:4,8,LQi1) =   fem_int(-1,SiLC)
    NLC(8,LQi1) = 0

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

    ! Ns*Usr(0) : fixed

    ELM(1:NEMAX,1:4,0,LQi2) = fem_int(1) * invDT &
         &                  + fem_int(8) * fem_int(0) * invDT
    NLC(0,LQi2) = LQi2

    ! Nonlinear term

    ELM(1:NEMAX,1:4,1,LQi2) = - 2.D0 * fem_int( 3,RUirV) + fem_int(2,UirVR) &
         &                  +(- 2.D0 * fem_int(10,RUirV) + fem_int(9,UirVR)) * fem_int(0)
    NLC(1,LQi2) = LQi2

    ! Nonlinear centrifugal force

    ELM(1:NEMAX,1:4,2,LQi2) = fem_int(2,UithVR) &
         &                  + fem_int(9,UithVR) * fem_int(0)
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

    use physical_constants, only : AEE, AME, rKeV, rMU0

    integer :: N

    ! Ni*UiTheta(0) : 0

    ELM(1:NEMAX,1:4, 0,LQi3) = fem_int(1) * invDT
    NLC( 0,LQi3) = LQi3

    ! Nonlinear term

    ELM(1:NEMAX,1:4, 1,LQi3) = - 2.D0 * fem_int(3,RUirV)
    NLC( 1,LQi3) = LQi3

    ! Viscosity force

    ELM(1:NEMAX,1:4, 2,LQi3) = - 4.D0 * fem_int(18,rMui) &
         &                     - 4.D0 * fem_int(39,rMuiNi,dPNiV) &
         &                     - 4.D0 * fem_int( 3,rMui)
    NLC( 2,LQi3) = LQi3

    ! Poroidal E force

    ELM(1:NEMAX,1:4, 3,LQi3) = - (PZ * AEE / AMI) * fem_int(2,PNiV) / rMU0
    NLC( 3,LQi3) = LQm2

    ! v x B force

    ELM(1:NEMAX,1:4, 4,LQi3) = - (PZ * AEE / AMI) * fem_int(2,BphV)
    NLC( 4,LQi3) = LQi2

    ! Neoclassical viscosity force

    ELM(1:NEMAX,1:4, 5,LQi3) = - fem_int(2,rNuiNC)
    NLC( 5,LQi3) = LQi3

    ! Collisional friction force

    ELM(1:NEMAX,1:4, 6,LQi3) = - (AME / AMI) * fem_int(2,rNueiEI)
    NLC( 6,LQi3) = LQi3

    ELM(1:NEMAX,1:4, 7,LQi3) =   (AME / AMI) * fem_int(2,rNuei)
    NLC( 7,LQi3) = LQe3

    ! Collisional friction with beam ions

    ELM(1:NEMAX,1:4, 8,LQi3) = - (AMB / AMI) * fem_int(2,rNubiBI)
    NLC( 8,LQi3) = LQi3

    ELM(1:NEMAX,1:4, 9,LQi3) =   (AMB / AMI) * fem_int(2,rNubi)
    NLC( 9,LQi3) = LQb3

    IF(MDLWTB == 0) THEN

       ! Wave interaction force (electron driven)

       ELM(1:NEMAX,1:4,10,LQi3) =   1.D0 / AMI * fem_int(2,FWthe)
       NLC(10,LQi3) = LQe3

       ELM(1:NEMAX,1:4,11,LQi3) =   2.D0 / AMI * fem_int(36,AphV,FWthphe)
       NLC(11,LQi3) = LQe4

       ELM(1:NEMAX,1:4,12,LQi3) = - 1.D0 / AMI * fem_int(44,FWthe,WPM)
       NLC(12,LQi3) = LQe1

       ! Wave interaction force (ion driven)

       ELM(1:NEMAX,1:4,13,LQi3) = - 1.D0 / AMI * fem_int(2,FWthi)
       NLC(13,LQi3) = LQi3

       ELM(1:NEMAX,1:4,14,LQi3) = - 2.D0 / AMI * fem_int(36,AphV,FWthphi)
       NLC(14,LQi3) = LQi4

       ELM(1:NEMAX,1:4,15,LQi3) =   1.D0 / AMI * fem_int(44,FWthi,WPM)
       NLC(15,LQi3) = LQi1

       N = 0
    ELSEIF(MDLWTB == 1) THEN

       ! Wave interaction force (electron driven)

       ELM(1:NEMAX,1:4,10,LQi3) =   2.D0 / AMI * fem_int(17,WNthe)
       NLC(10,LQi3) = LQe1

!       ELM(1:NEMAX,1:4,11,LQi3) =   2.D0 / AMI * fem_int(36,Phi,WEMthe)
!       NLC(11,LQi3) = LQe1
       
!       ELM(1:NEMAX,1:4,12,LQi3) =   2.D0 / AMI * fem_int(15,WWthe)
!       NLC(12,LQi3) = LQe1
       
!       ELM(1:NEMAX,1:4,13,LQi3) =   2.D0 / AMI * rKeV * fem_int(17,WT1the)
!       NLC(13,LQi3) = LQe5
       
!       ELM(1:NEMAX,1:4,14,LQi3) = - 2.D0 / AMI * fem_int(17,WT2the)
!       NLC(14,LQi3) = LQe1

       ! Wave interaction force (ion driven)

       ELM(1:NEMAX,1:4,15,LQi3) = - 2.D0 / AMI * fem_int(17,WNthi)
       NLC(15,LQi3) = LQi1

!       ELM(1:NEMAX,1:4,16,LQi3) = - 2.D0 / AMI * fem_int(36,Phi,WEMthi)
!       NLC(16,LQi3) = LQi1
       
!       ELM(1:NEMAX,1:4,17,LQi3) = - 2.D0 / AMI * fem_int(15,WWthi)
!       NLC(17,LQi3) = LQi1
       
!       ELM(1:NEMAX,1:4,18,LQi3) = - 2.D0 / AMI * rKeV * fem_int(17,WT1thi)
!       NLC(18,LQi3) = LQi5
       
!       ELM(1:NEMAX,1:4,19,LQi3) =   2.D0 / AMI * fem_int(17,WT2thi)
!       NLC(19,LQi3) = LQi1

       N = 4
    END IF

    ! Loss to divertor

    ELM(1:NEMAX,1:4,16+N,LQi3) = - 2.D0 * fem_int(2,rNuL)
    NLC(16+N,LQi3) = LQi3

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,17+N,LQi3) = - fem_int(2,rNu0i)
    NLC(17+N,LQi3) = LQi3

    ! Charge exchange force

    ELM(1:NEMAX,1:4,18+N,LQi3) = - fem_int(2,rNuiCX)
    NLC(18+N,LQi3) = LQi3

!   ! Loss cone loss

!   PELM(1:NEMAX,1:4,19+N,LQi3) = fem_int(21,SiLCth)
!   NLC(19+N,LQi3) = 0

!   ! Helical Neoclassical viscosity force
!
!   ELM(1:NEMAX,1:4,20+N,LQi3) = - (1.D0 - UHth * UHth) * fem_int(2,rNuiHL)
!   NLC(20+N,LQi3) = LQi3
!
!   ELM(1:NEMAX,1:4,21+N,LQi3) = UHph * UHth * fem_int(22,rNuiHL)
!   NLC(21+N,LQi3) = LQi4

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

    ! Uiphi'(0) : 0

    ELM(1:NEMAX,1:4, 0,LQi4) = fem_int(1) * invDT
    NLC( 0,LQi4) = LQi4

    ! Nonlinear term

    ELM(1:NEMAX,1:4, 1,LQi4) = - 2.D0 * fem_int(3,RUirV)
    NLC( 1,LQi4) = LQi4

    ! Viscosity force

    ELM(1:NEMAX,1:4, 2,LQi4) = - 4.D0 *(fem_int(18,rMui) + fem_int(39,rMuiNi,dPNiV))
    NLC( 2,LQi4) = LQi4

    ! Toroidal E force

    ELM(1:NEMAX,1:4, 3,LQi4) = - (PZ * AEE / AMI) * fem_int(2,PNiV)
    NLC( 3,LQi4) = LQm3

    ! v x B force

    ELM(1:NEMAX,1:4, 4,LQi4) = - 2.D0 * (PZ * AEE / AMI) * fem_int(6,AphV)
    NLC( 4,LQi4) = LQi2

    ! Collisional friction with bulk ions

    ELM(1:NEMAX,1:4, 5,LQi4) = - (AME / AMI) * fem_int(2,rNueiEI)
    NLC( 5,LQi4) = LQi4

    ELM(1:NEMAX,1:4, 6,LQi4) =   (AME / AMI) * fem_int(2,rNuei)
    NLC( 6,LQi4) = LQe4

    ! Collisional friction with beam ions

    ELM(1:NEMAX,1:4, 7,LQi4) = - (AMB / AMI) * fem_int(2,rNubiBI)
    NLC( 7,LQi4) = LQi4

    ELM(1:NEMAX,1:4, 8,LQi4) =   (AMB / AMI) * fem_int(2,rNubi)
    NLC( 8,LQi4) = LQb4

    IF(MDLWTB == 0) THEN

       ! Wave int_arrayeraction force (electron driven)

       ELM(1:NEMAX,1:4, 9,LQi4) = - 1.D0 / AMI * fem_int(2,FWpheBB)
       NLC( 9,LQi4) = LQe3

       ELM(1:NEMAX,1:4,10,LQi4) =   1.D0 / AMI * fem_int(2,FWpheBB2)
       NLC(10,LQi4) = LQe4

       ELM(1:NEMAX,1:4,11,LQi4) =   1.D0 / AMI * fem_int(44,FWpheBB,WPM)
       NLC(11,LQi4) = LQe1

       ! Wave int_arrayeraction force (ion driven)

       ELM(1:NEMAX,1:4,12,LQi4) =   1.D0 / AMI * fem_int(2,FWphiBB)
       NLC(12,LQi4) = LQi3

       ELM(1:NEMAX,1:4,13,LQi4) = - 1.D0 / AMI * fem_int(2,FWphiBB2)
       NLC(13,LQi4) = LQi4

       ELM(1:NEMAX,1:4,14,LQi4) = - 1.D0 / AMI * fem_int(44,FWphiBB,WPM)
       NLC(14,LQi4) = LQi1

    END IF

    ! Loss to divertor

    ELM(1:NEMAX,1:4,15,LQi4) = - 2.D0 * fem_int(2,rNuL)
    NLC(15,LQi4) = LQi4

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,16,LQi4) = - fem_int(2,rNu0i)
    NLC(16,LQi4) = LQi4

    ! Charge exchange force

    ELM(1:NEMAX,1:4,17,LQi4) = - fem_int(2,rNuiCX)
    NLC(17,LQi4) = LQi4

    ! Loss cone loss

    PELM(1:NEMAX,1:4,18,LQi4) =   fem_int(-1,SiLCph)
    NLC(18,LQi4) = 0

!   ! Helical Neoclassical viscosity force

!   ELM(1:NEMAX,1:4,19,LQi4) = UHth * UHph / 2.D0 * fem_int(2,rNuiHL)
!   NLC(19,LQi4) = LQi3

!   ELM(1:NEMAX,1:4,20,LQi4) = - (1.D0 - UHph * UHph) * fem_int(2,rNuiHL)
!   NLC(20,LQi4) = LQi4

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

    ! Temperature evolution

    IF(MDFIXT == 0) THEN
       ELM(1:NEMAX,1:4, 0,LQi5) =   1.5D0 * fem_int(1) * invDT
       NLC( 0,LQi5) = LQi5

       ! Convection transport

       ELM(1:NEMAX,1:4, 1,LQi5) = - 5.D0 * fem_int(3,RUirV)
       NLC( 1,LQi5) = LQi5

       ! Conduction transport

!      ELM(1:NEMAX,1:4, 2,LQi5) = - 6.D0 *(fem_int(18,Chii) + fem_int(39,ChiiNi,dPNiV))
!      NLC( 2,LQi5) = LQi5
       ELM(1:NEMAX,1:4, 2,LQi5) = - 6.D0 * fem_int(18,Chii)
       NLC( 2,LQi5) = LQi5

       ELM(1:NEMAX,1:4,14,LQi5) =   6.D0 * fem_int(41,Chii,PTiV)
       NLC(14,LQi5) = LQi1

       ! Joule heating

       ELM(1:NEMAX,1:4, 3,LQi5) =   PZ * AEE / rKeV * fem_int(2,EthVR)
       NLC( 3,LQi5) = LQi3

       ELM(1:NEMAX,1:4, 4,LQi5) =   PZ * AEE / rKeV * fem_int(2,EphV)
       NLC( 4,LQi5) = LQi4

       ! Collisional transfer with electrons

       ELM(1:NEMAX,1:4, 5,LQi5) = - 1.5d0 * fem_int(2,rNuTeiEI)
       NLC( 5,LQi5) = LQi5

       ELM(1:NEMAX,1:4, 6,LQi5) =   1.5d0 * fem_int(2,rNuTei)
       NLC( 6,LQi5) = LQe5

       ! Collisional heating with beam

!      ELM(1:NEMAX,1:4, 7,LQi5) = - 0.5D0 * AMb * (PNBCD * Vb) / rKeV &
!           &                * fem_int(2,rNubiBI)
!      NLC( 7,LQi5) = LQi4
!
!      ELM(1:NEMAX,1:4, 8,LQi5) =   0.5D0 * AMb * (PNBCD * Vb) / rKeV &
!           &                * fem_int(2,rNubi)
!      NLC( 8,LQi5) = LQb4

       ELM(1:NEMAX,1:4, 7,LQi5) = - 0.5D0 * AMb / rKeV * fem_int(28,Vbidir,rNubiBI)
       NLC( 7,LQi5) = LQi4

       ELM(1:NEMAX,1:4, 8,LQi5) =   0.5D0 * AMb / rKeV * fem_int(28,Vbidir,rNubi)
       NLC( 8,LQi5) = LQb4

       ! Loss to diverter

       ELM(1:NEMAX,1:4, 9,LQi5) = -         1.D0   / PZ * fem_int( 2,rNuL)
       NLC( 9,LQi5) = LQi5

       PELM(1:NEMAX,1:4,10,LQi5) =          PNeDIV / PZ * fem_int(-2,rNuL,PTiV)
       NLC(10,LQi5) = 0

       ELM(1:NEMAX,1:4,11,LQi5) = - 1.5D0               * fem_int( 2,rNuLTi)
       NLC(11,LQi5) = LQi5

       PELM(1:NEMAX,1:4,12,LQi5) =  1.5D0 * PTiDIV      * fem_int(-2,rNuLTi,PNiV)
       NLC(12,LQi5) = 0

       ! Direct heating (RF)

       PELM(1:NEMAX,1:4,13,LQi5) = 1.D0 / (1.D20 * rKeV) * fem_int(-1,PRFi)
       NLC(13,LQi5) = 0

!       NLCMAX(LQi5) = 13
       NLCMAX(LQi5) = 14
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

    ! NBI particle source

    PELM(1:NEMAX,1:4,1,LQb1) =   fem_int(-1,SNB)
    NLC(1,LQb1) = 0

    ! Relaxation to thermal ions

    ELM(1:NEMAX,1:4,2,LQb1) = - fem_int(2,rNuB)
    NLC(2,LQb1) = LQb1

    NLCMAX(LQb1) = 2
    RETURN
  END SUBROUTINE LQb1CC

!***************************************************************
!
!   Beam Ion Poloidal Flow
!
!***************************************************************

  SUBROUTINE LQb3CC

    use physical_constants, only : AEE, rMU0

    ! Ubth(0) : 0

    ELM(1:NEMAX,1:4,0,LQb3) = fem_int(1) * invDT
    NLC(0,LQb3) = LQb3

    ! Poroidal E force

    ELM(1:NEMAX,1:4,1,LQb3) = - (PZ * AEE / AMB) * fem_int(2,PNbV) / rMU0
    NLC(1,LQb3) = LQm2

!    ! Neoclassical viscosity force

!    ELM(1:NEMAX,1:4,2,LQb3) = - fem_int(2,rNuiNC)
!    NLC(2,LQb3) = LQb3

    ! Collisional friction force with electrons

    ELM(1:NEMAX,1:4,3,LQb3) = - fem_int(2,rNube)
    NLC(3,LQb3) = LQb3

    ELM(1:NEMAX,1:4,4,LQb3) =   fem_int(2,rNubeBE)
    NLC(4,LQb3) = LQe3

    ! Collisional friction force with ions

    ELM(1:NEMAX,1:4,5,LQb3) = - fem_int(2,rNubi)
    NLC(5,LQb3) = LQb3

    ELM(1:NEMAX,1:4,6,LQb3) =   fem_int(2,rNubiBI)
    NLC(6,LQb3) = LQi3

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,7,LQb3) = - fem_int(2,rNu0b)
    NLC(7,LQb3) = LQb3

    ! Charge exchange force

    ELM(1:NEMAX,1:4,8,LQb3) = - fem_int(2,rNuiCX)
    NLC(8,LQb3) = LQb3

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

    use physical_constants, only : AEE

    ! - UbPhi(0)' : 0

    ELM(1:NEMAX,1:4,0,LQb4) = fem_int(1) * invDT
    NLC(0,LQb4) = LQb4

    ! Toroidal E force

    ELM(1:NEMAX,1:4,1,LQb4) = - (PZ * AEE / AMB) * fem_int(2,PNbV)
    NLC(1,LQb4) = LQm3

    ! Collisional friction with electrons

    ELM(1:NEMAX,1:4,2,LQb4) = - fem_int(2,rNube)
    NLC(2,LQb4) = LQb4

    ELM(1:NEMAX,1:4,3,LQb4) =   fem_int(2,rNubeBE)
    NLC(3,LQb4) = LQe4

    ! Collisional friction with ions

    ELM(1:NEMAX,1:4,4,LQb4) = - fem_int(2,rNubi)
    NLC(4,LQb4) = LQb4

    ELM(1:NEMAX,1:4,5,LQb4) =   fem_int(2,rNubiBI)
    NLC(5,LQb4) = LQi4

    ! Collisional friction force with neutrals

    ELM(1:NEMAX,1:4,6,LQb4) = - fem_int(2,rNu0b)
    NLC(6,LQb4) = LQb4

    ! Charge exchange force

    ELM(1:NEMAX,1:4,7,LQb4) = - fem_int(2,rNuiCX)
    NLC(7,LQb4) = LQb4

    ! NBI momentum source

    PELM(1:NEMAX,1:4,8,LQb4) = (PNBCD * Vb) * fem_int(-1,SNB)
    NLC(8,LQb4) = 0

    ! Ubphi(NRMAX) : 0

    NLCMAX(LQb4) = 8
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

    ! NBI particle source

    PELM(1:NEMAX,1:4,4,LQn2) = fem_int(-1,SNB)
    NLC(4,LQn2) = 0

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
       PLC(NLCMAX(LQ),LQ,NR) = VAL * DTf(LQ)
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
