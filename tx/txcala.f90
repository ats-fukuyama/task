!     $Id$
module coefficients
  use core_module
  implicit none
  private
  real(8), dimension(:,:,:,:), allocatable :: ELM, PELM
  public :: TXCALA
  
contains

!***************************************************************
!
!   Calculate ALC, BLC, CLC, PLC
!
!***************************************************************

  SUBROUTINE TXCALA

    use commons
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

    CALL BOUNDARY(0    ,LQm1,0)
    CALL BOUNDARY(0    ,LQm2,0)
    CALL BOUNDARY(0    ,LQm5,0)
    CALL BOUNDARY(0    ,LQe2,0)
!    CALL BOUNDARY(NRMAX,LQe2,0)
!    CALL BOUNDARY(0    ,LQe3,0)
!    CALL BOUNDARY(NRMAX,LQe3,0)
    CALL BOUNDARY(0    ,LQi2,0)
!    CALL BOUNDARY(NRMAX,LQi2,0)
!    CALL BOUNDARY(0    ,LQi3,0)
!    CALL BOUNDARY(NRMAX,LQi3,0)
!    CALL BOUNDARY(0    ,LQb3,0)
!    CALL BOUNDARY(NRMAX,LQb3,0)

    !  Gas puffing on NRMAX only

    CALL BOUNDARY(NRMAX,LQn1,1,R(NRMAX)*rGASPF)

    deallocate(ELM,PELM)

    RETURN
  END SUBROUTINE TXCALA

!***************************************************************
!
!   Poisson Equation
!
!**************************************************************

  SUBROUTINE LQm1CC

    use commons
    use physical_constants, only : AEE, EPS0

    INTEGER :: NE

    ! Er(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4,1,LQm1,NE) = EPS0 / (AEE * 1.D20) * fem_integral(3,NE,R)
       NLC(1,LQm1) = LQm1

       ELM(1:4,2,LQm1,NE) =        fem_integral(2,NE,R)
       NLC(2,LQm1) = LQe1

       ELM(1:4,3,LQm1,NE) = - PZ * fem_integral(2,NE,R)
       NLC(3,LQm1) = LQi1

       ELM(1:4,4,LQm1,NE) = - PZ * fem_integral(2,NE,R)
       NLC(4,LQm1) = LQb1
    END DO

    NLCMAX(LQm1) = 4
    RETURN
  END SUBROUTINE LQm1CC

!***************************************************************
!
!   Ampere's Law: Eth
!
!***************************************************************

  SUBROUTINE LQm2CC

    use commons
    use physical_constants, only : AEE, VC, rMU0

    INTEGER :: NE

    ! Etheta(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4,0,LQm2,NE) = 1.D0 / (VC**2 * DT) * fem_integral(2,NE,R)
       NLC(0,LQm2) = LQm2

       ! rot Bphi

       ELM(1:4,1,LQm2,NE) = - fem_integral(2,NE,R)
       NLC(1,LQm2) = LQm5

       ! Electron current

       ELM(1:4,2,LQm2,NE) =   rMU0      * AEE * 1.D20 * fem_integral(21,NE)
       NLC(2,LQm2) = LQe3

       ! Ion current

       ELM(1:4,3,LQm2,NE) = - rMU0 * PZ * AEE * 1.D20 * fem_integral(21,NE)
       NLC(3,LQm2) = LQi3

       ! Beam ion current

       ELM(1:4,4,LQm2,NE) = - rMU0 * PZ * AEE * 1.D20 * fem_integral(21,NE)
       NLC(4,LQm2) = LQb3
    END DO

    ! Etheta(NRMAX) : r * Etheta = const.

    NLCMAX(LQm2) = 4
    RETURN
  END SUBROUTINE LQm2CC

!***************************************************************
!
!   Ampere's Law: Ephi
!
!***************************************************************

  SUBROUTINE LQm3CC

    use commons
    use physical_constants, only : AEE, VC, rMU0

    INTEGER :: NE

    DO NE = 1, NEMAX
       ELM(1:4,0,LQm3,NE) = 1.D0 / (VC**2 * DT) * fem_integral(2,NE,R)
       NLC(0,LQm3) = LQm3

       ! rot Btheta

       ELM(1:4,1,LQm3,NE) = fem_integral(2,NE,R)
       NLC(1,LQm3) = LQm4

       ! Electron current

       ELM(1:4,2,LQm3,NE) =   rMU0 *      AEE * 1.D20 * fem_integral(2,NE,R)
       NLC(2,LQm3) = LQe4

       ! Ion current

       ELM(1:4,3,LQm3,NE) = - rMU0 * PZ * AEE * 1.D20 * fem_integral(2,NE,R)
       NLC(3,LQm3) = LQi4

       ! Beam ion current

       ELM(1:4,4,LQm3,NE) = - rMU0 * PZ * AEE * 1.D20 * fem_integral(2,NE,R)
       NLC(4,LQm3) = LQb4

       ! Virtual current for helical system

       PELM(1:4,5,LQm3,NE) = - rMU0 * fem_integral(15,NE,AJV)
       NLC(5,LQm3) = 0
    END DO

    NLCMAX(LQm3) = 5
    RETURN
  END SUBROUTINE LQm3CC

!**************************************************************
!
!   Faraday's Law : Btheta
!
!***************************************************************

  SUBROUTINE LQm4CC

    use commons
    INTEGER :: NE

    ! Btheta(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4,0,LQm4,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC(0,LQm4) = LQm4

       ! rot Ephi

       ELM(1:4,1,LQm4,NE) = - fem_integral(12,NE,R)
       NLC(1,LQm4) = LQm3
    END DO

    ! Btheta(NRMAX) : fixed

    NLCMAX(LQm4) = 1
    RETURN
  END SUBROUTINE LQm4CC

!***************************************************************
!
!   Faraday's Law : Bphi
!
!***************************************************************

  SUBROUTINE LQm5CC

    use commons
    INTEGER :: NE

    DO NE = 1, NEMAX
       ELM(1:4,0,LQm5,NE) = 1.D0 / DT * fem_integral(21,NE)
       NLC(0,LQm5) = LQm5

       ! rot Etheta

       ELM(1:4,1,LQm5,NE) = - fem_integral(5,NE,R)
       NLC(1,LQm5) = LQm2

       ELM(1:4,2,LQm5,NE) =   fem_integral(1,NE)
       NLC(2,LQm5) = LQm2

       ELM(1:4,3,LQm5,NE) =   fem_integral(11,NE)
       NLC(3,LQm5) = LQm2
    END DO

    NLCMAX(LQm5) = 3
    RETURN
  END SUBROUTINE LQm5CC

!***************************************************************
!
!   Electron Density Equation
!
!***************************************************************

  SUBROUTINE LQe1CC

    use commons
    INTEGER :: NE
    REAL(8), DIMENSION(0:NRMAX) :: TMP

    TMP(0:NRMAX) = rNuION(0:NRMAX) * PNeV(0:NRMAX) / (PN01V(0:NRMAX) + PN02V(0:NRMAX))

    DO NE = 1, NEMAX
       ELM(1:4,0,LQe1,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC(0,LQe1) = LQe1

       ! Convection

       ELM(1:4,1,LQe1,NE) = - fem_integral(3,NE,R)
       NLC(1,LQe1) = LQe2

       ! Ionization of n01 and n02

       ELM(1:4,2,LQe1,NE) =   fem_integral(16,NE,TMP)
       NLC(2,LQe1) = LQn1

       ELM(1:4,3,LQe1,NE) =   fem_integral(16,NE,TMP)
       NLC(3,LQe1) = LQn2

       ! Loss to divertor

       ELM(1:4,4,LQe1,NE) = - fem_integral(16,NE,rNuL)
       NLC(4,LQe1) = LQe1

       PELM(1:4,5,LQe1,NE) =  PNeDIV * fem_integral(15,NE,rNuL)
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

    use commons
    use physical_constants, only : AEE, AME, rKEV

    INTEGER :: NE

    ! Ns*Usr(0) : fixed

    DO NE = 1, NEMAX
       ELM(1:4,0,LQe2,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC(0,LQe2) = LQe2

       ! Nonlinear term

       ELM(1:4,1,LQe2,NE) = - fem_integral(2,NE,UerV) - fem_integral(17,NE,UerV) &
            &               - fem_integral(18,NE,UerV)
       NLC(1,LQe2) = LQe2

       ! Nonlinear centrifugal force

       ELM(1:4,2,LQe2,NE) =   fem_integral(16,NE,UethV)
       NLC(2,LQe2) = LQe3

       ! Pressure gradient force

       ELM(1:4,3,LQe2,NE) = - rKeV / AME * fem_integral(5,NE,R)
       NLC(3,LQe2) = LQe5

       ! Radial E force

       ELM(1:4,4,LQe2,NE) = - (AEE / AME) * fem_integral(16,NE,PNeV)
       NLC(4,LQe2) = LQm1

       ! v x B force

       ELM(1:4,5,LQe2,NE) = - (AEE / AME) * fem_integral(22,NE,BphV)
       NLC(5,LQe2) = LQe3

       ELM(1:4,6,LQe2,NE) =   (AEE / AME) * fem_integral(16,NE,BthV)
       NLC(6,LQe2) = LQe4
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

    use commons
    use physical_constants, only : AEE, AME, rKEV

    INTEGER :: NE
    REAL(8), DIMENSION(0:NRMAX) :: TMP1, TMP2, TMP3, TMP4, TMP5, TMP6, TMP7, TMP8, TMP9, TMP

    TMP1(0:NRMAX) = rMue (0:NRMAX) * UethV(0:NRMAX)
    TMP2(0:NRMAX) = rNuei(0:NRMAX) * PNeV(0:NRMAX) / PNiV(0:NRMAX)
    TMP3(0:NRMAX) = rNube(0:NRMAX) * PNbV(0:NRMAX) / PNeV(0:NRMAX)
    TMP4(0:NRMAX) = FWthe(0:NRMAX) * WPM(0:NRMAX)
    TMP5(0:NRMAX) = FWthi(0:NRMAX) * WPM(0:NRMAX)
!!$    TMP6(0:NRMAX) = WWthe(0:NRMAX) * PNeV(0:NRMAX) /(PTeV(0:NRMAX) * rKEV)
!!$    TMP7(0:NRMAX) = WDthe(0:NRMAX)                 /(PTeV(0:NRMAX) * rKEV)
!!$    TMP8(0:NRMAX) = WWthi(0:NRMAX) * PNiV(0:NRMAX) /(PTiV(0:NRMAX) * rKEV)
!!$    TMP9(0:NRMAX) = WDthi(0:NRMAX)                 /(PTiV(0:NRMAX) * rKEV)
    TMP(0:NRMAX) = FWthi(0:NRMAX) * PNiV(0:NRMAX) / PNeV(0:NRMAX)

    ! Ns*UsTheta(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4,0,LQe3,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC(0,LQe3) = LQe3

       ! Nonlinear term

       ELM(1:4, 1,LQe3,NE) = - 3.d0 * fem_integral(16,NE,UerV) - fem_integral(23,NE,UerV) &
            &               -        fem_integral(24,NE,UerV)
       NLC( 1,LQe3) = LQe3

       ! Viscosity force

       ELM(1:4, 2,LQe3,NE) =   fem_integral(18,NE,rMue) - fem_integral(25,NE,rMue)
       NLC( 2,LQe3) = LQe3

       ELM(1:4, 3,LQe3,NE) = - fem_integral(5,NE,TMP1)  + fem_integral(19,NE,TMP1)
       NLC( 3,LQe3) = LQe1

       ! Poloidal E force

       ELM(1:4, 4,LQe3,NE) = - (AEE / AME) * fem_integral(16,NE,PNeV)
       NLC( 4,LQe3) = LQm2

       ! v x B force

       ELM(1:4, 5,LQe3,NE) =   (AEE / AME) * fem_integral(16,NE,BphV)
       NLC( 5,LQe3) = LQe2

       ! Neoclassical viscosity force

       ELM(1:4, 6,LQe3,NE) = - fem_integral(22,NE,rNueNC)
       NLC( 6,LQe3) = LQe3

       ! Collisional friction force with ions

       ELM(1:4, 7,LQe3,NE) = - fem_integral(22,NE,rNuei)
       NLC( 7,LQe3) = LQe3

       ELM(1:4, 8,LQe3,NE) =   fem_integral(22,NE,TMP2)
       NLC( 8,LQe3) = LQi3

       ! Collisional friction with beam ions

       ELM(1:4, 9,LQe3,NE) = - (AMB / AME) * fem_integral(22,NE,TMP3)
       NLC( 9,LQe3) = LQe3

       ELM(1:4,10,LQe3,NE) =   (AMB / AME) * fem_integral(22,NE,rNube)
       NLC(10,LQe3) = LQb3

       ! Wave interaction force (electron driven)

       ELM(1:4,11,LQe3,NE) = - 1.D0 / AME * fem_integral(22,NE,FWthe)
       NLC(11,LQe3) = LQe3

       ELM(1:4,12,LQe3,NE) =   1.D0 / AME * fem_integral(22,NE,TMP4)
       NLC(12,LQe3) = LQe1

       ! Wave interaction force (NRon driven)

       ELM(1:4,13,LQe3,NE) =   1.D0 / AME * fem_integral(22,NE,FWthi)
       NLC(13,LQe3) = LQi3

       ELM(1:4,14,LQe3,NE) = - 1.D0 / AME * fem_integral(22,NE,TMP5)
       NLC(14,LQe3) = LQi1
!!$
!!$       ELM(1:4,13,LQe3,NE) =   1.D0 / AME * fem_integral(22,NE,FWthi)
!!$       NLC(13,LQe3) = LQi3
!!$
!!$       ELM(1:4,14,LQe3,NE) = - 1.D0 / AME * fem_integral(22,NE,TMP)
!!$       NLC(14,LQe3) = LQe3

       ! Loss to divertor

       ELM(1:4,15,LQe3,NE) = - 2.D0 * fem_integral(22,NE,rNuL)
       NLC(15,LQe3) = LQe3

       ! Collisional friction force with neutrals

       ELM(1:4,16,LQe3,NE) = - fem_integral(22,NE,rNu0e)
       NLC(16,LQe3) = LQe3

       ! Helical neoclassical viscosity force

       ELM(1:4,17,LQe3,NE) = - (1.D0 - UHth * UHth) * fem_integral(22,NE,rNueHL)
       NLC(17,LQe3) = LQe3

       ELM(1:4,18,LQe3,NE) = UHph * UHth / 2.D0 * fem_integral(16,NE,rNueHL)
       NLC(18,LQe3) = LQe4

!!$       ! Wave interaction force
!!$
!!$       ELM(1:4,19,LQe3,NE) =       1.D0 / AME * fem_integral(18,NE,WWthe)
!!$       NLC(19,LQe3) = LQe1
!!$
!!$       ELM(1:4,20,LQe3,NE) = -      AEE / AME * fem_integral(16,NE,TMP6)
!!$       NLC(20,LQe3) = LQm1
!!$
!!$       ELM(1:4,21,LQe3,NE) =       rKEV / AME * fem_integral(18,NE,TMP7)
!!$       NLC(21,LQe3) = LQe5
!!$
!!$       ELM(1:4,22,LQe3,NE) = -     1.D0 / AME * fem_integral(18,NE,WDthe)
!!$       NLC(22,LQe3) = LQe1
!!$
!!$       ELM(1:4,23,LQe3,NE) = -     1.D0 / AME * fem_integral(18,NE,WWthi)
!!$       NLC(23,LQe3) = LQi1
!!$
!!$       ELM(1:4,24,LQe3,NE) =   PZ * AEE / AME * fem_integral(16,NE,TMP8)
!!$       NLC(24,LQe3) = LQm1
!!$
!!$       ELM(1:4,25,LQe3,NE) = -     rKEV / AME * fem_integral(18,NE,TMP9)
!!$       NLC(25,LQe3) = LQi5
!!$
!!$       ELM(1:4,26,LQe3,NE) =       1.D0 / AME * fem_integral(18,NE,WDthi)
!!$       NLC(26,LQe3) = LQi1
    END DO

    ! Ns*UsTheta(NRMAX) : 0

!    NLCMAX(LQe3) = 26
    NLCMAX(LQe3) = 18
    RETURN
  END SUBROUTINE LQe3CC

!***************************************************************
!
!   Electron Toroidal Flow
!
!***************************************************************

  SUBROUTINE LQe4CC

    use commons
    use physical_constants, only : AEE, AME

    INTEGER :: NE
    REAL(8), DIMENSION(0:NRMAX) :: TMP1, TMP2, TMP3

    TMP1(0:NRMAX) = rMue (0:NRMAX) * UephV(0:NRMAX)
    TMP2(0:NRMAX) = rNuei(0:NRMAX) * PNeV(0:NRMAX) / PNiV(0:NRMAX)
    TMP3(0:NRMAX) = rNube(0:NRMAX) * PNbV(0:NRMAX) / PNeV(0:NRMAX)

    ! Uephi(0)' : 0

    DO NE = 1, NEMAX
       ELM(1:4, 0,LQe4,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC( 0,LQe4) = LQe4

       ! Nonlinear term

       ELM(1:4, 1,LQe4,NE) = - fem_integral(2,NE,UerV) - fem_integral(17,NE,UerV) &
            &                - fem_integral(18,NE,UerV)
       NLC( 1,LQe4) = LQe4

       ! Viscosity force

       ELM(1:4, 2,LQe4,NE) = - fem_integral(19,NE,rMue)
       NLC( 2,LQe4) = LQe4

       ELM(1:4, 3,LQe4,NE) =   fem_integral(19,NE,TMP1)
       NLC( 3,LQe4) = LQe1

       ! Toroidal E force

       ELM(1:4, 4,LQe4,NE) = - (AEE / AME) * fem_integral(16,NE,PNeV)
       NLC( 4,LQe4) = LQm3

       ! v x B force

       ELM(1:4, 5,LQe4,NE) = - (AEE / AME) * fem_integral(16,NE,BthV)
       NLC( 5,LQe4) = LQe2

       ! Collisional friction with bulk ions

       ELM(1:4, 6,LQe4,NE) = - fem_integral(16,NE,rNuei)
       NLC( 6,LQe4) = LQe4

       ELM(1:4, 7,LQe4,NE) =   fem_integral(16,NE,TMP2)
       NLC( 7,LQe4) = LQi4

       ! Collisional friction with beam ions

       ELM(1:4, 8,LQe4,NE) = - (AMB / AME) * fem_integral(16,NE,TMP3)
       NLC( 8,LQe4) = LQe4

       ELM(1:4, 9,LQe4,NE) =   (AMB / AME) * fem_integral(16,NE,rNube)
       NLC( 9,LQe4) = LQb4

       ! Loss to divertor

       ELM(1:4,10,LQe4,NE) = - 2.D0 * fem_integral(16,NE,rNuL)
       NLC(10,LQe4) = LQe4

       ! Collisional friction force with neutrals

       ELM(1:4,11,LQe4,NE) = - fem_integral(16,NE,rNu0e)
       NLC(11,LQe4) = LQe4

       ! Helical neoclassical viscosity force

       ELM(1:4,12,LQe4,NE) =  UHth * UHph / 2.D0 * fem_integral(22,NE,rNueHL)
       NLC(12,LQe4) = LQe3

       ELM(1:4,13,LQe4,NE) = - (1.D0 - UHph * UHph) * fem_integral(16,NE,rNueHL)
       NLC(13,LQe4) = LQe4
    END DO

    NLCMAX(LQe4) = 13
    RETURN
  END SUBROUTINE LQe4CC

!***************************************************************
!
!   Electron Energy Transport: Te
!
!***************************************************************

  SUBROUTINE LQe5CC

    use commons
    use physical_constants, only : AEE, rKEV

    INTEGER :: NE
    REAL(8), DIMENSION(0:NRMAX) :: TMP1, TMP2, TMP3, TMP4

    ! Fixed Temperature

    IF (ABS(Chie0) == 0.D0) THEN
       DO NE = 1, NEMAX
          ELM(1:4,0,LQe5,NE) = fem_integral(2,NE,R) / DT
       END DO
       NLC(0,LQe5) = LQe5
       NLCMAX(LQe5) = 0
    ELSE
       TMP1(0:NRMAX) = Chie  (0:NRMAX) * PTeV(0:NRMAX)
       TMP2(0:NRMAX) = rNuTei(0:NRMAX) * PNeV(0:NRMAX) / PNiV(0:NRMAX)
       TMP3(0:NRMAX) = rNube (0:NRMAX) * PNbV(0:NRMAX) / PNeV(0:NRMAX)
       TMP4(0:NRMAX) = rNuL  (0:NRMAX) * PNeV(0:NRMAX)

       ! Temperature evolution

       DO NE = 1, NEMAX
          ELM(1:4,0,LQe5,NE) =   1.5D0 / DT * fem_integral(2,NE,R)
          NLC(0,LQe5) = LQe5

          ! Convection transport

          ELM(1:4,1,LQe5,NE) = - 2.5D0 * (fem_integral(2,NE,UerV) &
               &               + fem_integral(17,NE,UerV) + fem_integral(18,NE,UerV))
          NLC(1,LQe5) = LQe5

          ! Conduction transport

          ELM(1:4,2,LQe5,NE) = - 1.5D0 * fem_integral(19,NE,Chie)
          NLC(2,LQe5) = LQe5

          ELM(1:4,3,LQe5,NE) =   1.5D0 * fem_integral(19,NE,TMP1)
          NLC(3,LQe5) = LQe1

          ! Joule heating

          ELM(1:4,4,LQe5,NE) = - AEE / rKeV * fem_integral(22,NE,EthV)
          NLC(4,LQe5) = LQe3

          ELM(1:4,5,LQe5,NE) = - AEE / rKeV * fem_integral(16,NE,EphV)
          NLC(5,LQe5) = LQe4

          ! Collisional transfer with ions

          ELM(1:4,6,LQe5,NE) = - fem_integral(16,NE,rNuTei)
          NLC(6,LQe5) = LQe5

          ELM(1:4,7,LQe5,NE) =   fem_integral(16,NE,TMP2)
          NLC(7,LQe5) = LQi5

          ! Collisional heating with beam

          ELM(1:4,8,LQe5,NE) = - 0.5D0 * AMb * (PNBCD * Vb) / rKeV &
               &             * fem_integral(16,NE,TMP3)
          NLC(8,LQe5) = LQe4

          ELM(1:4,9,LQe5,NE) =   0.5D0 * AMb * (PNBCD * Vb) / rKeV &
               &             * fem_integral(16,NE,rNube)
          NLC(9,LQe5) = LQb4

          ! Loss to diverter

          ELM(1:4,10,LQe5,NE) = - 2.5D0 * fem_integral(16,NE,rNuL)
          NLC(10,LQe5) = LQe5

          PELM(1:4,11,LQe5,NE) =  2.5D0 * PTeDIV * fem_integral(15,NE,TMP4)
          NLC(11,LQe5) = 0

          ! Direct heating (RF)

          PELM(1:4,12,LQe5,NE) =   1.D0 / (1.D20 * rKeV) * fem_integral(15,NE,PRFe)
          NLC(12,LQe5) = 0
       END DO

       NLCMAX(LQe5) = 12
    END IF
    RETURN
  END SUBROUTINE LQe5CC

!***************************************************************
!
!   Ion Density Equation
!
!***************************************************************

  SUBROUTINE LQi1CC

    use commons
    INTEGER :: NE
    REAL(8), DIMENSION(0:NRMAX) :: TMP

    TMP(0:NRMAX) = rNuION(0:NRMAX) * PNeV(0:NRMAX) / (PN01V(0:NRMAX) + PN02V(0:NRMAX))

    DO NE = 1, NEMAX
       ELM(1:4,0,LQi1,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC(0,LQi1) = LQi1

       ! Convection

       ELM(1:4,1,LQi1,NE) = - fem_integral(3,NE,R)
       NLC(1,LQi1) = LQi2

       ! Ionization of n01 and n02

       ELM(1:4,2,LQi1,NE) =     1.D0 / PZ * fem_integral(16,NE,TMP)
       NLC(2,LQi1) = LQn1

       ELM(1:4,3,LQi1,NE) =     1.D0 / PZ * fem_integral(16,NE,TMP)
       NLC(3,LQi1) = LQn2

       ! Loss to divertor

       ELM(1:4,4,LQi1,NE) = -   1.D0 / PZ * fem_integral(16,NE,rNuL)
       NLC(4,LQi1) = LQe1

       PELM(1:4,5,LQi1,NE) =  PNeDIV / PZ * fem_integral(15,NE,rNuL)
       NLC(5,LQi1) = 0

       ! Particle source from beam ion

       ELM(1:4,6,LQi1,NE) =   fem_integral(16,NE,rNuB)
       NLC(6,LQi1) = LQb1

       ! NBI kick up ions

       PELM(1:4,7,LQi1,NE) = - fem_integral(15,NE,SNB)
       NLC(7,LQi1) = 0

       ! Loss cone loss

       PELM(1:4,8,LQi1,NE) =   fem_integral(15,NE,SiLC)
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

    use commons
    use physical_constants, only : AEE, rKEV

    INTEGER :: NE

    ! Ns*Usr(0) : fixed

    DO NE = 1, NEMAX
       ELM(1:4,0,LQi2,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC(0,LQi2) = LQi2

       ! Nonlinear term

       ELM(1:4,1,LQi2,NE) = - fem_integral(2,NE,UirV) - fem_integral(17,NE,UirV) &
            &               - fem_integral(18,NE,UirV)
       NLC(1,LQi2) = LQi2

       ! Nonlinear centrifugal force

       ELM(1:4,2,LQi2,NE) =   fem_integral(16,NE,UithV)
       NLC(2,LQi2) = LQi3

       ! Pressure gradient force

       ELM(1:4,3,LQi2,NE) = - rKEV / AMI * fem_integral(5,NE,R)
       NLC(3,LQi2) = LQi5

       ! Radial E force

       ELM(1:4,4,LQi2,NE) =   (PZ * AEE / AMI) * fem_integral(16,NE,PNiV)
       NLC(4,LQi2) = LQm1

       ! v x B force

       ELM(1:4,5,LQi2,NE) =   (PZ * AEE / AMI) * fem_integral(22,NE,BphV)
       NLC(5,LQi2) = LQi3

       ELM(1:4,6,LQi2,NE) = - (PZ * AEE / AMI) * fem_integral(16,NE,BthV)
       NLC(6,LQi2) = LQi4
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

    use commons
    use physical_constants, only : AEE, AME, rKEV

    INTEGER :: NE
    REAL(8), DIMENSION(0:NRMAX) :: TMP1, TMP2, TMP3, TMP4, TMP5, TMP6, TMP7, TMP8, TMP9, TMP

    TMP1(0:NRMAX) = rMui (0:NRMAX) * UithV(0:NRMAX)
    TMP2(0:NRMAX) = rNuei(0:NRMAX) * PNeV(0:NRMAX) / PNiV(0:NRMAX)
    TMP3(0:NRMAX) = rNubi(0:NRMAX) * PNbV(0:NRMAX) / PNiV(0:NRMAX)
    TMP4(0:NRMAX) = FWthe(0:NRMAX) * WPM(0:NRMAX) 
    TMP5(0:NRMAX) = FWthi(0:NRMAX) * WPM(0:NRMAX)
!!$    TMP6(0:NRMAX) = WWthe(0:NRMAX) * PNeV(0:NRMAX) /(PTeV(0:NRMAX) * rKEV)
!!$    TMP7(0:NRMAX) = WDthe(0:NRMAX)                 /(PTeV(0:NRMAX) * rKEV)
!!$    TMP8(0:NRMAX) = WWthi(0:NRMAX) * PNiV(0:NRMAX) /(PTiV(0:NRMAX) * rKEV)
!!$    TMP9(0:NRMAX) = WDthi(0:NRMAX)                 /(PTiV(0:NRMAX) * rKEV)
    TMP (0:NRMAX) = FWthi(0:NRMAX) * PNiV(0:NRMAX) / PNeV(0:NRMAX)

    ! Ni*UiTheta(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4,0,LQi3,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC(0,LQi3) = LQi3

       ! Nonlinear term

       ELM(1:4, 1,LQi3,NE) = - 3.d0 * fem_integral(16,NE,UirV) - fem_integral(23,NE,UirV) &
            &                -        fem_integral(24,NE,UirV)
       NLC( 1,LQi3) = LQi3
 
      ! Viscosity force

       ELM(1:4, 2,LQi3,NE) =   fem_integral(18,NE,rMui) - fem_integral(25,NE,rMui)
       NLC( 2,LQi3) = LQi3

       ELM(1:4, 3,LQi3,NE) = - fem_integral(5,NE,TMP1)  + fem_integral(19,NE,TMP1)
       NLC( 3,LQi3) = LQi1

       ! Poroidal E force

       ELM(1:4, 4,LQi3,NE) =   (PZ * AEE / AMI) * fem_integral(16,NE,PNiV)
       NLC( 4,LQi3) = LQm2

       ! v x B force

       ELM(1:4, 5,LQi3,NE) = - (PZ * AEE / AMI) * fem_integral(16,NE,BphV)
       NLC( 5,LQi3) = LQi2

       ! Neoclassical viscosity force

       ELM(1:4, 6,LQi3,NE) = - fem_integral(22,NE,rNuiNC)
       NLC( 6,LQi3) = LQi3

       ! Collisional friction force

       ELM(1:4, 7,LQi3,NE) = - (AME / AMI) * fem_integral(22,NE,TMP2)
       NLC( 7,LQi3) = LQi3

       ELM(1:4, 8,LQi3,NE) =   (AME / AMI) * fem_integral(22,NE,rNuei)
       NLC( 8,LQi3) = LQe3

       ! Collisional friction with beam ions

       ELM(1:4, 9,LQi3,NE) = - (AMB / AMI) * fem_integral(22,NE,TMP3)
       NLC( 9,LQi3) = LQi3

       ELM(1:4,10,LQi3,NE) =   (AMB / AMI) * fem_integral(22,NE,rNubi)
       NLC(10,LQi3) = LQb3

       ! Wave interaction force (electron driven)

       ELM(1:4,11,LQi3,NE) =   1.D0 / AMI * fem_integral(22,NE,FWthe)
       NLC(11,LQi3) = LQe3

       ELM(1:4,12,LQi3,NE) = - 1.D0 / AMI * fem_integral(22,NE,TMP4)
       NLC(12,LQi3) = LQe1

       ! Wave interaction force (NRon driven)

       ELM(1:4,13,LQi3,NE) = - 1.D0 / AMI * fem_integral(22,NE,FWthi)
       NLC(13,LQi3) = LQi3

       ELM(1:4,14,LQi3,NE) =   1.D0 / AMI * fem_integral(22,NE,TMP5)
       NLC(14,LQi3) = LQi1
!!$
!!$       ELM(1:4,13,LQi3,NE) = - 1.D0 / AMI * fem_integral(22,NE,FWthi)
!!$       NLC(13,LQi3) = LQi3
!!$
!!$       ELM(1:4,14,LQi3,NE) =   1.D0 / AMI * fem_integral(22,NE,TMP)
!!$       NLC(14,LQi3) = LQe3

       ! Loss to divertor

       ELM(1:4,15,LQi3,NE) = - 2.D0 * fem_integral(22,NE,rNuL)
       NLC(15,LQi3) = LQi3

       ! Collisional friction force with neutrals

       ELM(1:4,16,LQi3,NE) = - fem_integral(22,NE,rNu0i)
       NLC(16,LQi3) = LQi3

       ! Charge exchange force

       ELM(1:4,17,LQi3,NE) = - fem_integral(22,NE,rNuiCX)
       NLC(17,LQi3) = LQi3

       ! Loss cone loss

       PELM(1:4,18,LQi3,NE) = fem_integral(15,NE,SiLCth)
       NLC(18,LQi3) = 0

       ! Helical Neoclassical viscosity force

       ELM(1:4,19,LQi3,NE) = - (1.D0 - UHth * UHth) * fem_integral(22,NE,rNuiHL)
       NLC(19,LQi3) = LQi3

       ELM(1:4,20,LQi3,NE) = UHph * UHth * fem_integral(16,NE,rNuiHL)
       NLC(20,LQi3) = LQi4

       ! Wave interaction force

!!$       ELM(1:4,21,LQi3,NE) = -     1.D0 / AMI * fem_integral(18,NE,WWthe)
!!$       NLC(21,LQi3) = LQe1
!!$
!!$       ELM(1:4,22,LQi3,NE) =        AEE / AMI * fem_integral(16,NE,TMP6)
!!$       NLC(22,LQi3) = LQm1
!!$
!!$       ELM(1:4,23,LQi3,NE) = -     rKEV / AMI * fem_integral(18,NE,TMP7)
!!$       NLC(23,LQi3) = LQe5
!!$
!!$       ELM(1:4,24,LQi3,NE) =       1.D0 / AMI * fem_integral(18,NE,WDthe)
!!$       NLC(24,LQi3) = LQe1
!!$
!!$       ELM(1:4,25,LQi3,NE) =       1.D0 / AMI * fem_integral(18,NE,WWthi)
!!$       NLC(25,LQi3) = LQi1
!!$
!!$       ELM(1:4,26,LQi3,NE) = - PZ * AEE / AMI * fem_integral(16,NE,TMP8)
!!$       NLC(26,LQi3) = LQm1
!!$
!!$       ELM(1:4,27,LQi3,NE) =       rKEV / AMI * fem_integral(18,NE,TMP9)
!!$       NLC(27,LQi3) = LQi5
!!$
!!$       ELM(1:4,28,LQi3,NE) = -     1.D0 / AMI * fem_integral(18,NE,WDthi)
!!$       NLC(28,LQi3) = LQi1
    END DO

    ! Ns*UsTheta(NRMAX) : 0

!    NLCMAX(LQi3) = 28
    NLCMAX(LQi3) = 20
    RETURN
  END SUBROUTINE LQi3CC

!***************************************************************
!
!   Ion Toroidal Flow
!
!***************************************************************

  SUBROUTINE LQi4CC

    use commons
    use physical_constants, only : AEE, AME

    INTEGER :: NE
    REAL(8), DIMENSION(0:NRMAX) :: TMP1, TMP2, TMP3

    TMP1(0:NRMAX) = rMui (0:NRMAX) * UiphV(0:NRMAX)
    TMP2(0:NRMAX) = rNuei(0:NRMAX) * PNeV(0:NRMAX) / PNiV(0:NRMAX)
    TMP3(0:NRMAX) = rNubi(0:NRMAX) * PNbV(0:NRMAX) / PNiV(0:NRMAX)

    ! Uiphi'(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4, 0,LQi4,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC( 0,LQi4) = LQi4

       ! Nonlinear term
       
       ELM(1:4, 1,LQi4,NE) = - fem_integral(2,NE,UirV) - fem_integral(17,NE,UirV) &
            &                - fem_integral(18,NE,UirV)
       NLC( 1,LQi4) = LQi4

       ! Viscosity force

       ELM(1:4, 2,LQi4,NE) = - fem_integral(19,NE,rMui)
       NLC( 2,LQi4) = LQi4

       ELM(1:4, 3,LQi4,NE) =   fem_integral(19,NE,TMP1)
       NLC( 3,LQi4) = LQi1

       ! Toroidal E force

       ELM(1:4, 4,LQi4,NE) = (PZ * AEE / AMI) * fem_integral(16,NE,PNiV)
       NLC( 4,LQi4) = LQm3

       ! v x B force

       ELM(1:4, 5,LQi4,NE) = (PZ * AEE / AMI) * fem_integral(16,NE,BthV)
       NLC( 5,LQi4) = LQi2

       ! Collisional friction with bulk ions

       ELM(1:4, 6,LQi4,NE) = - (AME / AMI) * fem_integral(16,NE,TMP2)
       NLC( 6,LQi4) = LQi4

       ELM(1:4, 7,LQi4,NE) =   (AME / AMI) * fem_integral(16,NE,rNuei)
       NLC( 7,LQi4) = LQe4

       ! Collisional friction with beam ions

       ELM(1:4, 8,LQi4,NE) = - (AMB / AMI) * fem_integral(16,NE,TMP3)
       NLC( 8,LQi4) = LQi4

       ELM(1:4, 9,LQi4,NE) =   (AMB / AMI) * fem_integral(16,NE,rNubi)
       NLC( 9,LQi4) = LQb4

       ! Loss to divertor

       ELM(1:4,10,LQi4,NE) = - 2.D0 * fem_integral(16,NE,rNuL)
       NLC(10,LQi4) = LQi4

       ! Collisional friction force with neutrals

       ELM(1:4,11,LQi4,NE) = - fem_integral(16,NE,rNu0i)
       NLC(11,LQi4) = LQi4

       ! Charge exchange force

       ELM(1:4,12,LQi4,NE) = - fem_integral(16,NE,rNuiCX)
       NLC(12,LQi4) = LQi4

       ! Loss cone loss

       PELM(1:4,13,LQi4,NE) =   fem_integral(15,NE,SiLCph)
       NLC(13,LQi4) = 0

       ! Helical Neoclassical viscosity force

       ELM(1:4,14,LQi4,NE) = UHth * UHph / 2.D0 * fem_integral(22,NE,rNuiHL)
       NLC(14,LQi4) = LQi3

       ELM(1:4,15,LQi4,NE) = - (1.D0 - UHph * UHph) * fem_integral(16,NE,rNuiHL)
       NLC(15,LQi4) = LQi4
    END DO

    NLCMAX(LQi4) = 15
    RETURN
  END SUBROUTINE LQi4CC

!***************************************************************
!
!  Ion Energy Transport: Ti
!
!***************************************************************

  SUBROUTINE LQi5CC

    use commons
    use physical_constants, only : AEE, rKEV

    INTEGER :: NE
    REAL(8), DIMENSION(0:NRMAX) :: TMP1, TMP2, TMP3, TMP4

    ! Fixed temperature

    IF (ABS(Chii0) == 0.D0) THEN
       DO NE = 1, NEMAX
          ELM(1:4,0,LQi5,NE) = fem_integral(2,NE,R) / DT
       END DO
       NLC(0,LQi5) = LQi5
       NLCMAX(LQi5) = 0

    ELSE
       ! Temperature evolution

       TMP1(0:NRMAX) = Chii  (0:NRMAX) * PTiV(0:NRMAX)
       TMP2(0:NRMAX) = rNuTei(0:NRMAX) * PNeV(0:NRMAX) / PNiV(0:NRMAX)
       TMP3(0:NRMAX) = rNubi (0:NRMAX) * PNbV(0:NRMAX) / PNiV(0:NRMAX)
       TMP4(0:NRMAX) = rNuL  (0:NRMAX) * PNiV(0:NRMAX)

       DO NE = 1, NEMAX
          ELM(1:4,0,LQi5,NE) = 1.5D0 / DT * fem_integral(2,NE,R)
          NLC(0,LQi5) = LQi5

          ! Convection transport

          ELM(1:4,1,LQi5,NE) = - 2.5D0 *(fem_integral(2,NE,UirV) &
               &               + fem_integral(17,NE,UirV) + fem_integral(18,NE,UirV))
          NLC(1,LQi5) = LQi5

          ! Conduction transport
          
          ELM(1:4,2,LQi5,NE) = - 1.5D0 * fem_integral(19,NE,Chii)
          NLC(2,LQi5) = LQi5

          ELM(1:4,3,LQi5,NE) =   1.5D0 * fem_integral(19,NE,TMP1)
          NLC(3,LQi5) = LQi1

          ! Joule heating

          ELM(1:4,4,LQi5,NE) =   PZ * AEE / rKeV * fem_integral(22,NE,EthV)
          NLC(4,LQi5) = LQi3

          ELM(1:4,5,LQi5,NE) =   PZ * AEE / rKeV * fem_integral(16,NE,EphV)
          NLC(5,LQi5) = LQi4

          ! Collisional transfer with electrons

          ELM(1:4,6,LQi5,NE) = - fem_integral(16,NE,TMP2)
          NLC(6,LQi5) = LQi5

          ELM(1:4,7,LQi5,NE) =   fem_integral(16,NE,rNuTei)
          NLC(7,LQi5) = LQe5

          ! Collisional heating with beam

          ELM(1:4,8,LQi5,NE) = - 0.5D0 * AMb * (PNBCD * Vb) / rKeV &
               &             * fem_integral(16,NE,TMP3)
          NLC(8,LQi5) = LQi4

          ELM(1:4,9,LQi5,NE) =   0.5D0 * AMb * (PNBCD * Vb) / rKeV &
               &             * fem_integral(16,NE,rNubi)
          NLC(9,LQi5) = LQb4

          ! Loss to diverter

          ELM(1:4,10,LQi5,NE) = - 2.5D0 * fem_integral(16,NE,rNuL)
          NLC(10,LQi5) = LQi5

          PELM(1:4,11,LQi5,NE) =  2.5D0 * PTiDIV * fem_integral(15,NE,TMP4)
          NLC(11,LQi5) = 0

          ! Direct heating (RF)

          PELM(1:4,12,LQi5,NE) = 1.D0 / (1.D20 * rKeV) &
               &              * fem_integral(15,NE,PRFi)
          NLC(12,LQi5) = 0
       END DO

       NLCMAX(LQi5) = 12
    END IF
    RETURN
  END SUBROUTINE LQi5CC

!***************************************************************
!
!   Beam Ion Density
!
!***************************************************************

  SUBROUTINE LQb1CC

    use commons
    INTEGER :: NE

    DO NE = 1, NEMAX
       ELM(1:4,0,LQb1,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC(0,LQb1) = LQb1

       ! NBI particle source

       PELM(1:4,1,LQb1,NE) =   fem_integral(15,NE,SNB)
       NLC(1,LQb1) = 0

       ! Relaxation to thermal ions

       ELM(1:4,2,LQb1,NE) = - fem_integral(16,NE,rNuB)
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

    use commons
    use physical_constants, only : AEE

    INTEGER :: NE

    ! Ubth(0) : 0

    DO NE = 1, NEMAX
       ELM(1:4,0,LQb3,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC(0,LQb3) = LQb3

       ! Nonlinear centrifugal force

       ELM(1:4,1,LQb3,NE) = fem_integral(16,NE,UbthV)
       NLC(1,LQb3) = LQb3

       ! Radial E force

       ELM(1:4,2,LQb3,NE) =   PZ * (AEE / AMB) * fem_integral(16,NE,PNbV)
       NLC(2,LQb3) = LQm1

       ! v x B force

       ELM(1:4,3,LQb3,NE) =   PZ * (AEE / AMB) * fem_integral(22,NE,BphV)
       NLC(3,LQb3) = LQb3

       ELM(1:4,4,LQb3,NE) = - PZ * (AEE / AMB) * fem_integral(16,NE,BthV)
       NLC(4,LQb3) = LQb4
    END DO

    ! Ubth(NRMAX) : 0

    NLCMAX(LQb3) = 4
    RETURN
  END SUBROUTINE LQb3CC

!***************************************************************
!
!   Beam Ion Toroidal Flow
!
!***************************************************************

  SUBROUTINE LQb4CC

    use commons
    INTEGER :: NE
    REAL(8), DIMENSION(0:NRMAX) :: TMP1, TMP2

    ! - UbPhi(0)' : 0

    TMP1(0:NRMAX) = rNube(0:NRMAX) * PNbV(0:NRMAX) / PNeV(0:NRMAX)
    TMP2(0:NRMAX) = rNubi(0:NRMAX) * PNbV(0:NRMAX) / PNiV(0:NRMAX)

    DO NE = 1, NEMAX
       ELM(1:4,0,LQb4,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC(0,LQb4) = LQb4

       ! Collisional friction with electrons

       ELM(1:4,1,LQb4,NE) = - fem_integral(16,NE,rNube)
       NLC(1,LQb4) = LQb4

       ELM(1:4,2,LQb4,NE) =   fem_integral(16,NE,TMP1)
       NLC(2,LQb4) = LQe4

       ! Collisional friction with ions

       ELM(1:4,3,LQb4,NE) = - fem_integral(16,NE,rNubi)
       NLC(3,LQb4) = LQb4

       ELM(1:4,4,LQb4,NE) =   fem_integral(16,NE,TMP2)
       NLC(4,LQb4) = LQi4

       ! NBI momentum source

       PELM(1:4,5,LQb4,NE) = (PNBCD * Vb) * fem_integral(15,NE,SNB)
       NLC(5,LQb4) = 0
    END DO

    ! Ubphi(NRMAX) : 0

    NLCMAX(LQb4) = 5
    RETURN
  END SUBROUTINE LQb4CC

!***************************************************************
!
!   Slow Neutral Transport: n01
!
!***************************************************************

  SUBROUTINE LQn1CC

    use commons
    INTEGER :: NE
    REAL(8), DIMENSION(0:NRMAX) :: TMP1, TMP2

    TMP1(0:NRMAX) = rNuION(0:NRMAX) * PNeV(0:NRMAX) / (PN01V(0:NRMAX) + PN02V(0:NRMAX))
    TMP2(0:NRMAX) = rNuiCX(0:NRMAX) * PNiV(0:NRMAX) / (PN01V(0:NRMAX) + PN02V(0:NRMAX))

    DO NE = 1, NEMAX
       ELM(1:4,0,LQn1,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC(0,LQn1) = LQn1

       !  Diffusion of neutrals

       ELM(1:4,1,LQn1,NE) = - fem_integral(19,NE,D01)
       NLC(1,LQn1) = LQn1

       ! Ionization

       ELM(1:4,2,LQn1,NE) = - 1.D0 / PZ * fem_integral(16,NE,TMP1)
       NLC(2,LQn1) = LQn1

       ! Generation of fast neutrals by charge exchange

       ELM(1:4,3,LQn1,NE) = - fem_integral(16,NE,TMP2)
       NLC(3,LQn1) = LQn1

       ! Recycling from divertor

       ELM(1:4,4,LQn1,NE) =   rGamm0 / PZ * fem_integral(16,NE,rNuL)
       NLC(4,LQn1) = LQe1

       PELM(1:4,5,LQn1,NE) = - rGamm0 * PNeDIV / PZ * fem_integral(15,NE,rNuL)
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

    use commons
    INTEGER :: NE
    REAL(8), DIMENSION(0:NRMAX) :: TMP1, TMP2

    TMP1(0:NRMAX) = rNuION(0:NRMAX) * PNeV(0:NRMAX) / (PN01V(0:NRMAX) + PN02V(0:NRMAX))
    TMP2(0:NRMAX) = rNuiCX(0:NRMAX) * PNiV(0:NRMAX) / (PN01V(0:NRMAX) + PN02V(0:NRMAX))

    DO NE = 1, NEMAX
       ELM(1:4,0,LQn2,NE) = 1.D0 / DT * fem_integral(2,NE,R)
       NLC(0,LQn2) = LQn2

       !  Diffusion of neutrals

       ELM(1:4,1,LQn2,NE) = - fem_integral(19,NE,D02)
       NLC(1,LQn2) = LQn2

       ! Ionization

       ELM(1:4,2,LQn2,NE) = - 1.D0 / PZ * fem_integral(16,NE,TMP1)
       NLC(2,LQn2) = LQn2

       ! Generation of fast neutrals by charge exchange

       ELM(1:4,3,LQn2,NE) = fem_integral(16,NE,TMP2)
       NLC(3,LQn2) = LQn1

       ! NBI particle source

       PELM(1:4,4,LQn2,NE) = fem_integral(15,NE,SNB)
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

    use commons, only : ALC, BLC, CLC, PLC, NLCR, NLCMAX, NCM
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

end module coefficients
