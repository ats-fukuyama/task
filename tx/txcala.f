!     $Id$
!
!     ***************************************************************
!
!        Calculate ALC, BLC, CLC, PLC
!
!     ***************************************************************
!
      SUBROUTINE TXCALA

      INCLUDE 'txcomm.inc'

      ALC(0:NCM,1:NQMAX,0:NRMAX) = 0.D0
      BLC(0:NCM,1:NQMAX,0:NRMAX) = 0.D0
      CLC(0:NCM,1:NQMAX,0:NRMAX) = 0.D0
      NLC(0:NCM,1:NQMAX,0:NRMAX) = 1
      PLC(1:NCM,1:NQMAX,0:NRMAX) = 0.D0
      NLCMAX(1:NQM) = 0.D0

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

      RETURN
      END
!
!     ***************************************************************
!
!        Poisson Equation  (HI)
!
!     ***************************************************************
!
      SUBROUTINE LQm1CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR
      REAL(8) :: FACTOR
! Er(0) = 0

      NR = 0
          BLC(1,LQm1,NR) = 1.D0
          NLC(1,LQm1,NR) = LQm1

      FACTOR = EPS0 / (AEE * 1.D20)

!!!      DO NR = 1, NRMAX-1
      DO NR = 1, NRMAX
          BLC(1,LQm1,NR) =  FACTOR * R(NR  ) / (DR * RHI(NR-1))
          CLC(1,LQm1,NR) = -FACTOR * R(NR-1) / (DR * RHI(NR-1))
          NLC(1,LQm1,NR) = LQm1
!
          CLC(2,LQm1,NR) = 1.D0
          NLC(2,LQm1,NR) = LQe1
!
          CLC(3,LQm1,NR) = - PZ
          NLC(3,LQm1,NR) = LQi1
!
          CLC(4,LQm1,NR) = - PZ
          NLC(4,LQm1,NR) = LQb1
      ENDDO
!!!
!!!      NR = NRMAX
!!!          BLC(1,LQm1,NR) = 1.D0
!!!          NLC(1,LQm1,NR) = LQm1
!
      NLCMAX(LQm1) = 4
      RETURN
      END
!
!     ***************************************************************
!
!        Ampere's Law: Eth (NR)
!
!     ***************************************************************
!
      SUBROUTINE LQm2CC

      INCLUDE 'txcomm.inc'
      
      INTEGER :: NR
! Etheta(0) : 0

      NR = 0
          BLC(1,LQm2,NR) = 1.D0
          NLC(1,LQm2,NR) = LQm2

      DO NR = 1, NRMAX - 1
          BLC(0,LQm2,NR) = 1.D0 / (VC**2 * DT)
          NLC(0,LQm2,NR) = LQm2

! rot Bphi

          BLC(1,LQm2,NR) = - 1.D0 / DR
          CLC(1,LQm2,NR) =   1.D0 / DR
          NLC(1,LQm2,NR) = LQm5

! Electron current

          BLC(2,LQm2,NR) =   rMU0 *      AEE * 1.D20 
          NLC(2,LQm2,NR) = LQe3

! Ion current

          BLC(3,LQm2,NR) = - rMU0 * PZ * AEE * 1.D20
          NLC(3,LQm2,NR) = LQi3

! Beam ion current

          BLC(4,LQm2,NR) = - rMU0 * PZ * AEE * 1.D20
          NLC(4,LQm2,NR) = LQb3
      ENDDO

! Etheta(NRMAX) : r * Etheta = const.

      NR = NRMAX
          BLC(1,LQm2,NR) = R(NR)
          CLC(1,LQm2,NR) =-R(NR-1)
          NLC(1,LQm2,NR) = LQm2

      NLCMAX(LQm2) = 4
      RETURN
      END
!
!     ***************************************************************
!
!        Ampere's Law: Ephi (HI)
!
!     ***************************************************************
!
      SUBROUTINE LQm3CC

      INCLUDE 'txcomm.inc'
      
      INTEGER :: NR

      DO NR = 0, NRMAX - 1
          BLC(0,LQm3,NR) = 1.D0 / (VC**2 * DT)
          NLC(0,LQm3,NR) = LQm3

! rot Btheta

          ALC(1,LQm3,NR) =   R(NR+1) / (DR * RHI(NR))
          BLC(1,LQm3,NR) = - R(NR  ) / (DR * RHI(NR))
          NLC(1,LQm3,NR) = LQm4

! Electron current

          BLC(2,LQm3,NR) =   rMU0 *      AEE * 1.D20
          NLC(2,LQm3,NR) = LQe4

! Ion current

          BLC(3,LQm3,NR) = - rMU0 * PZ * AEE * 1.D20
          NLC(3,LQm3,NR) = LQi4

! Beam ion current

          BLC(4,LQm3,NR) = - rMU0 * PZ * AEE * 1.D20
          NLC(4,LQm3,NR) = LQb4

! Virtual current for helical sysytem

          PLC(5,LQm3,NR) = - rMU0 * AJV(NR)
      ENDDO

! Out of region

      NR = NRMAX
          BLC(1,LQm3,NR) = 1.D0
          NLC(1,LQm3,NR) = LQm3

      NLCMAX(LQm3) = 5
      RETURN
      END
!
!     **************************************************************
!
!        Faraday's Law : Btheta (NR)
!
!     ***************************************************************
!
      SUBROUTINE LQm4CC

      INCLUDE 'txcomm.inc'
      
      INTEGER :: NR

! Btheta(0) : 0

      NR = 0
          BLC(1,LQm4,NR) = 1.D0
          NLC(1,LQm4,NR) = LQm4

      DO NR = 1, NRMAX - 1
          BLC(0,LQm4,NR) = 1.D0 / DT
          NLC(0,LQm4,NR) = LQm4

! rot Ephi

          BLC(1,LQm4,NR) =   1.D0 / DR
          CLC(1,LQm4,NR) = - 1.D0 / DR
          NLC(1,LQm4,NR) = LQm3
      ENDDO

! Btheta(NRMAX) : fixed

      NR = NRMAX
          BLC(1,LQm4,NR) = 1.D0
          NLC(1,LQm4,NR) = LQm4

          PLC(2,LQm4,NR) = - Bthb

      NLCMAX(LQm4) = 2
      RETURN
      END
!
!     ***************************************************************
!
!        Fraday's Law : Bphi (HI)
!
!     ***************************************************************
!
      SUBROUTINE LQm5CC

      INCLUDE 'txcomm.inc'
      
      INTEGER :: NR

      DO NR = 0, NRMAX - 1
          BLC(0,LQm5,NR) = 1.D0 / DT
          NLC(0,LQm5,NR) = LQm5

! rot Etheta

          ALC(1,LQm5,NR) = - R(NR+1) / (DR * RHI(NR))
          BLC(1,LQm5,NR) =   R(NR  ) / (DR * RHI(NR))
          NLC(1,LQm5,NR) = LQm2
      ENDDO

! Out of region

      NR = NRMAX
          BLC(1,LQm5,NR) = 1.D0
          NLC(1,LQm5,NR) = LQm5

      NLCMAX(LQm5) = 1
      RETURN
      END
!
!     ***************************************************************
!
!        Electron Density Equation (HI)
!
!     ***************************************************************
!
      SUBROUTINE LQe1CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR

      DO NR = 0, NRMAX - 1
          BLC(0,LQe1,NR) = 1.D0 / DT
          NLC(0,LQe1,NR) = LQe1

! Convection

          ALC(1,LQe1,NR) = - R(NR+1) / (RHI(NR) * DR)
          BLC(1,LQe1,NR) = + R(NR  ) / (RHI(NR) * DR)
          NLC(1,LQe1,NR) = LQe2

! Ionization of n01 and n02

          BLC(2,LQe1,NR) = rNuION(NR) 
     &                    * PNeHI(NR) / (PN01HI(NR) + PN02HI(NR))
          NLC(2,LQe1,NR) = LQn1

          BLC(3,LQe1,NR) = rNuION(NR) 
     &                    * PNeHI(NR) / (PN01HI(NR) + PN02HI(NR))
          NLC(3,LQe1,NR) = LQn2

! Loss to divertor

          BLC(4,LQe1,NR) = - rNuL(NR)
          NLC(4,LQe1,NR) = LQe1

          PLC(5,LQe1,NR) =   rNuL(NR) * PNeDIV

      ENDDO

!     OUT OF REGION

      NR = NRMAX
          BLC(1,LQe1,NR) = 1.D0
          NLC(1,LQe1,NR) = LQe1

      NLCMAX(LQe1) = 5
      RETURN
      END
!
!     ***************************************************************
!
!        Electron Radial Flow (NR)
!
!     ***************************************************************
!
      SUBROUTINE LQe2CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR

! Ns*Usr(0) : fixed

      NR = 0
          BLC(1,LQe2,NR) = 1.D0
          NLC(1,LQe2,NR) = LQe2

      DO NR = 1, NRMAX-1
          BLC(0,LQe2,NR) = 1.D0 / DT
          NLC(0,LQe2,NR) = LQe2

! Nonlinear term

          ALC(1,LQe2,NR) = - RHI(NR  ) * UerHI(NR  ) / (2 * R(NR) * DR)
          BLC(1,LQe2,NR) = - RHI(NR  ) * UerHI(NR  ) / (2 * R(NR) * DR)
     &                     + RHI(NR-1) * UerHI(NR-1) / (2 * R(NR) * DR)
          CLC(1,LQe2,NR) = + RHI(NR-1) * UerHI(NR-1) / (2 * R(NR) * DR)
          NLC(1,LQe2,NR) = LQe2

! Nonlinear centrifugal force

          BLC(2,LQe2,NR) = UethI(NR) / R(NR)
          NLC(2,LQe2,NR) = LQe3

! Pressure gradient force

          BLC(3,LQe2,NR) = - rKeV / (AME * DR)
          CLC(3,LQe2,NR) = + rKeV / (AME * DR)
          NLC(3,LQe2,NR) = LQe5

! Radial E force

          BLC(4,LQe2,NR) = - PNeI(NR) * (AEE/AME)
          NLC(4,LQe2,NR) = LQm1

! v x B force

          BLC(5,LQe2,NR) = - BphI(NR) * (AEE/AME)
          NLC(5,LQe2,NR) = LQe3

          BLC(6,LQe2,NR) =   BthI(NR) * (AEE/AME) / 2
          CLC(6,LQe2,NR) =   BthI(NR) * (AEE/AME) / 2
          NLC(6,LQe2,NR) = LQe4

      ENDDO

! Ns*Usr(NRMAX) : fixed or finite gradient

      NR = NRMAX
          BLC(1,LQe2,NR) = 1.D0
          NLC(1,LQe2,NR) = LQe2

      NLCMAX(LQe2) = 6
      RETURN
      END
!
!     ***************************************************************
!
!        Electron Poloidal Flow (NR)
!
!     ***************************************************************
!
      SUBROUTINE LQe3CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR
      REAL(8) :: rNueNCL, rNueiL, rNubeL, FWtheL, WPML, TMP, FWthiL,
     &           rNuLL, rNu0eL, rNueHLL

! Ns*UsTheta(0) : 0

      NR = 0
          BLC(1,LQe3,NR) = 1.D0
          NLC(1,LQe3,NR) = LQe3

      DO NR = 1, NRMAX - 1
          BLC(0,LQe3,NR) = 1.D0 / DT
          NLC(0,LQe3,NR) = LQe3

! Nonlinear term

          ALC(1,LQe3,NR) = - RHI(NR  )**2 * UerHI(NR  )
     &                       / (2 * R(NR)**2 * DR)
          BLC(1,LQe3,NR) = - RHI(NR  )**2 * UerHI(NR  )
     &                       / (2 * R(NR)**2 * DR)
     &                     + RHI(NR-1)**2 * UerHI(NR-1)
     &                       / (2 * R(NR)**2 * DR)
          CLC(1,LQe3,NR) = + RHI(NR-1)**2 * UerHI(NR-1)
     &                       / (2 * R(NR)**2 * DR)
          NLC(1,LQe3,NR) = LQe3

! Viscosity force

         IF (NR .EQ. 1) THEN
             ALC(2,LQe3,NR) =   RHI(NR  )**3 * PNeHI(NR  ) * rMue(NR  )
     &                       / (PNeI(NR+1) * R(NR+1) * R(NR)**2 * DR**2)
             BLC(2,LQe3,NR) = - RHI(NR  )**3 * PNeHI(NR  ) * rMue(NR  )
     &                       / (PNeI(NR  ) * R(NR  ) * R(NR)**2 * DR**2)
         ELSE
             ALC(2,LQe3,NR) =   RHI(NR  )**3 * PNeHI(NR  ) * rMue(NR  )
     &                       / (PNeI(NR+1) * R(NR+1) * R(NR)**2 * DR**2)
             BLC(2,LQe3,NR) = - RHI(NR  )**3 * PNeHI(NR  ) * rMue(NR  )
     &                       / (PNeI(NR  ) * R(NR  ) * R(NR)**2 * DR**2)
     &                       - RHI(NR-1)**3 * PNeHI(NR-1) * rMue(NR-1)
     &                       / (PNeI(NR  ) * R(NR  ) * R(NR)**2 * DR**2)
             CLC(2,LQe3,NR) =   RHI(NR-1)**3 * PNeHI(NR-1) * rMue(NR-1)
     &                       / (PNeI(NR-1) * R(NR-1) * R(NR)**2 * DR**2)
         ENDIF
          NLC(2,LQe3,NR) = LQe3

! Poloidal E force

          BLC(3,LQe3,NR) = - PNeI(NR) * (AEE/AME)
          NLC(3,LQe3,NR) = LQm2

! v x B force

          BLC(4,LQe3,NR) = BphI(NR) * (AEE/AME)
          NLC(4,LQe3,NR) = LQe2

! Neoclassical viscosity force

          rNueNCL = 0.5D0 * (rNueNC(NR) + rNueNC(NR-1))
          BLC(5,LQe3,NR) = - rNueNCL
          NLC(5,LQe3,NR) = LQe3

! Collisional friction force with ions

          rNueiL = 0.5D0 * (rNuei(NR) + rNuei(NR-1))
          BLC(6,LQe3,NR) = - rNueiL
          NLC(6,LQe3,NR) = LQe3

          BLC(7,LQe3,NR) = + rNueiL * PNeI(NR) / PNiI(NR)
          NLC(7,LQe3,NR) = LQi3

! Collisional friction with beam ions

          rNubeL = 0.5D0 * (rNube(NR) + rNube(NR-1))
          BLC(8,LQe3,NR) = - (AMB/AME) * rNubeL * PNbI(NR) / PNeI(NR)
          NLC(8,LQe3,NR) = LQe3

          BLC(9,LQe3,NR) =   (AMB/AME) * rNubeL
          NLC(9,LQe3,NR) = LQb3

! Wave interaction force (electron driven)

          FWtheL = 0.5D0 * (FWthe(NR) + FWthe(NR-1))
          WPML = 0.5D0 * (  WPM(NR) +   WPM(NR-1))
          BLC(10,LQe3,NR) = - FWtheL / AME
          NLC(10,LQe3,NR) = LQe3

          TMP =   FWtheL * WPML * R(NR) / AME
          BLC(11,LQe3,NR) = TMP / 2
          CLC(11,LQe3,NR) = TMP / 2
          NLC(11,LQe3,NR) = LQe1

! Wave interaction force (NRon driven)

          FWthiL = 0.5D0 * (FWthi(NR) + FWthi(NR-1))
          BLC(12,LQe3,NR) =   FWthiL / AME
          NLC(12,LQe3,NR) = LQi3

          TMP = - FWthiL * WPML * R(NR) / AME
          BLC(13,LQe3,NR) = TMP / 2
          CLC(13,LQe3,NR) = TMP / 2
          NLC(13,LQe3,NR) = LQi1

! Loss to divertor

          rNuLL = 0.5D0 * (rNuL(NR) + rNuL(NR-1))
          BLC(14,LQe3,NR) = - 2.D0 * rNuLL
          NLC(14,LQe3,NR) = LQe3

! Collisional friction force with neutrals

          rNu0eL = 0.5D0 * (rNu0e(NR) + rNu0e(NR-1))
          BLC(15,LQe3,NR) = - rNu0eL
          NLC(15,LQe3,NR) = LQe3

! Helical Neoclassical viscosity force

          rNueHLL =rNueHL(NR)
          BLC(16,LQe3,NR) = - rNueHLL*(1.D0-UHth*UHth)
          NLC(16,LQe3,NR) = LQe3
          BLC(17,LQe3,NR) =   rNueHLL*UHph*UHth/2
          CLC(17,LQe3,NR) =   rNueHLL*UHph*UHth/2
          NLC(17,LQe3,NR) = LQe4

      ENDDO

! Ns*UsTheta(NRMAX) : 0

      NR = NRMAX
          BLC(1,LQe3,NR) = 1
          NLC(1,LQe3,NR) = LQe3

      NLCMAX(LQe3) = 17
      RETURN
      END
!
!     ***************************************************************
!
!        Electron Toroidal Flow (NR)
!
!     ***************************************************************
!
      SUBROUTINE LQe4CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR
      REAL(8) :: rMueP, rMueM, rNueHLL

! Uephi(0)' : 0

      DO NR = 0, NRMAX - 1
          BLC(0,LQe4,NR) = 1.D0 / DT
          NLC(0,LQe4,NR) = LQe4

! Nonlinear term

          ALC(1,LQe4,NR) = - R(NR+1) * UerI(NR+1) / (2 * RHI(NR) * DR)
          BLC(1,LQe4,NR) = - R(NR+1) * UerI(NR+1) / (2 * RHI(NR) * DR)
     &                    + R(NR  ) * UerI(NR  ) / (2 * RHI(NR) * DR)
          CLC(1,LQe4,NR) = + R(NR  ) * UerI(NR  ) / (2 * RHI(NR) * DR)
          NLC(1,LQe4,NR) = LQe4

! Viscosity force

         IF(NR.EQ.0) THEN
             rMueP = 0.5D0*(rMue(NR) + rMue(NR+1))
             rMueM =        rMue(NR)
             ALC(2,LQe4,NR) =   R(NR+1) * PNeI(NR+1) * rMueP
     &                    / (PNeHI(NR+1) * RHI(NR) * DR**2)
             BLC(2,LQe4,NR) = - R(NR+1) * PNeI(NR+1) * rMueP
     &                    / (PNeHI(NR  ) * RHI(NR) * DR**2)
     &                       - R(NR  ) * PNeI(NR  ) * rMueM
     &                    / (PNeHI(NR  ) * RHI(NR) * DR**2)
     &                       + R(NR  ) * PNeI(NR  ) * rMueM
     &                    / (PNeHI(NR  ) * RHI(NR) * DR**2)
          ELSEIF(NR.EQ.NRMAX-1) THEN
             rMueP =        rMue(NR)
             rMueM = 0.5D0*(rMue(NR-1) + rMue(NR))
             BLC(2,LQe4,NR) =   R(NR+1) * PNeI(NR+1) * rMueP
     &                    / (PNeHI(NR  ) * RHI(NR) * DR**2)
     &                       - R(NR+1) * PNeI(NR+1) * rMueP
     &                    / (PNeHI(NR  ) * RHI(NR) * DR**2)
     &                       - R(NR  ) * PNeI(NR  ) * rMueM
     &                    / (PNeHI(NR  ) * RHI(NR) * DR**2)
             CLC(2,LQe4,NR) =   R(NR  ) * PNeI(NR  ) * rMueM
     &                    / (PNeHI(NR-1) * RHI(NR) * DR**2)
          ELSE
             rMueP = 0.5D0*(rMue(NR) + rMue(NR+1))
             rMueM = 0.5D0*(rMue(NR-1) + rMue(NR))
             ALC(2,LQe4,NR) =   R(NR+1) * PNeI(NR+1) * rMueP
     &                    / (PNeHI(NR+1) * RHI(NR) * DR**2)
             BLC(2,LQe4,NR) = - R(NR+1) * PNeI(NR+1) * rMueP
     &                    / (PNeHI(NR  ) * RHI(NR) * DR**2)
     &                       - R(NR  ) * PNeI(NR  ) * rMueM
     &                    / (PNeHI(NR  ) * RHI(NR) * DR**2)
             CLC(2,LQe4,NR) =   R(NR  ) * PNeI(NR  ) * rMueM
     &                    / (PNeHI(NR-1) * RHI(NR) * DR**2)
          ENDIF
          NLC(2,LQe4,NR) = LQe4

! Toroidal E force

          BLC(3,LQe4,NR) = - PNeHI(NR  ) * (AEE/AME)
          NLC(3,LQe4,NR) = LQm3

! v x B force

          ALC(4,LQe4,NR) = - BthI(NR+1) * (AEE/AME) / 2
          BLC(4,LQe4,NR) = - BthI(NR  ) * (AEE/AME) / 2
          NLC(4,LQe4,NR) = LQe2

! Collisional friction with bulk ions

          BLC(5,LQe4,NR) = - rNuei(NR)
          NLC(5,LQe4,NR) = LQe4

          BLC(6,LQe4,NR) = + rNuei(NR) * PNeI(NR) / PNiI(NR)
          NLC(6,LQe4,NR) = LQi4

! Collisional friction with beam ions

          BLC(7,LQe4,NR) = - (AMB/AME) * rNube(NR) * PNbI(NR) / PNeI(NR)
          NLC(7,LQe4,NR) = LQe4

          BLC(8,LQe4,NR) =   (AMB/AME) * rNube(NR)
          NLC(8,LQe4,NR) = LQb4

! Loss to divertor

          BLC(9,LQe4,NR) = - 2.D0 * rNuL(NR)
          NLC(9,LQe4,NR) = LQe4

! Collisional friction force with neutrals

          BLC(10,LQe4,NR) = - rNu0e(NR)
          NLC(10,LQe4,NR) = LQe4

! Helical Neoclassical viscosity force

          rNueHLL = 0.5D0 * (rNueHL(NR) + rNueHL(NR+1))
          ALC(11,LQe4,NR) =   rNueHL(NR+1)*UHth*UHph/2
          BLC(11,LQe4,NR) =   rNueHL(NR  )*UHth*UHph/2
          NLC(11,LQe4,NR) = LQe3
          BLC(12,LQe4,NR) = - rNueHLL*(1.D0-UHph*UHph)
          NLC(12,LQe4,NR) = LQe4

      ENDDO

! Out of region

      NR = NRMAX
          BLC(1,LQe4,NR) = 1.D0
          NLC(1,LQe4,NR) = LQe4

      NLCMAX(LQe4) = 12
      RETURN
      END
!
!     ***************************************************************
!
!        Electron Energy Transport: Te (HI)
!
!     ***************************************************************
!
      SUBROUTINE LQe5CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR
      REAL(8) :: ChieLP, ChieLM

! Fixed Temperature

      IF (ABS(Chie0) .EQ. 0.D0) THEN
         DO NR = 0, NRMAX
             BLC(0,LQe5,NR) = 1.D0
             NLC(0,LQe5,NR) = LQe5
         ENDDO

         NLCMAX(LQe5) = 0
      ELSE

! Temperature evolution

         DO NR = 0, NRMAX - 1
             BLC(0,LQe5,NR) = 1.5D0 / DT
             NLC(0,LQe5,NR) = LQe5

! Convection transport

            IF (NR.EQ.0) THEN
               ALC(1,LQe5,NR) = - 2.5D0 * R(NR+1) * UerI(NR+1)
     &                            / (2 * RHI(NR) * DR)
               BLC(1,LQe5,NR) = - 2.5D0 * R(NR+1) * UerI(NR+1)
     &                            / (2 * RHI(NR) * DR)
            ELSEIF (NR.EQ.NRMAX-1) THEN
               BLC(1,LQe5,NR) = - 2.5D0 * R(NR+1) * UerI(NR+1)
     &                            / (2 * RHI(NR) * DR)
     &                         - 2.5D0 * R(NR+1) * UerI(NR+1)
     &                            / (2 * RHI(NR) * DR)
     &                         + 2.5D0 * R(NR  ) * UerI(NR  )
     &                            / (2 * RHI(NR) * DR)
               CLC(1,LQe5,NR) = + 2.5D0 * R(NR  ) * UerI(NR  )
     &                            / (2 * RHI(NR) * DR)
            ELSE
               ALC(1,LQe5,NR) = - 2.5D0 * R(NR+1) * UerI(NR+1)
     &                            / (2 * RHI(NR) * DR)
               BLC(1,LQe5,NR) = - 2.5D0 * R(NR+1) * UerI(NR+1)
     &                            / (2 * RHI(NR) * DR)
     &                         + 2.5D0 * R(NR  ) * UerI(NR  )
     &                            / (2 * RHI(NR) * DR)
               CLC(1,LQe5,NR) = + 2.5D0 * R(NR  ) * UerI(NR  )
     &                            / (2 * RHI(NR) * DR)
            ENDIF
             NLC(1,LQe5,NR) = LQe5

! Conduction transport

            IF (NR.EQ.0) THEN
               ChieLP = 0.5D0 * (Chie(NR) + Chie(NR+1))
                ALC(2,LQe5,NR) = + 1.5D0 * R(NR+1) * PNeI(NR+1) * ChieLP
     &                          / (RHI(NR) * PNeHI(NR+1))
                BLC(2,LQe5,NR) = - 1.5D0 * R(NR+1) * PNeI(NR+1) * ChieLP
     &                          / (RHI(NR) * PNeHI(NR))
            ELSEIF (NR.EQ.NRMAX-1) THEN
               ChieLP =          Chie(NR)
               ChieLM = 0.5D0 * (Chie(NR-1) + Chie(NR))
                BLC(2,LQe5,NR) = + 1.5D0 * R(NR+1) * PNeI(NR+1) * ChieLP
     &                          / (RHI(NR) * PNeHI(NR))
     &                          - 1.5D0 * R(NR+1) * PNeI(NR+1) * ChieLP
     &                          / (RHI(NR) * PNeHI(NR))
     &                          - 1.5D0 * R(NR  ) * PNeI(NR  ) * ChieLM
     &                          / (RHI(NR) * PNeHI(NR))
                CLC(2,LQe5,NR) = + 1.5D0 * R(NR  ) * PNeI(NR  ) * ChieLM
     &                          / (RHI(NR) * PNeHI(NR-1))
            ELSE
               ChieLP = 0.5D0 * (Chie(NR) + Chie(NR+1))
               ChieLM = 0.5D0 * (Chie(NR-1) + Chie(NR))
                ALC(2,LQe5,NR) = + 1.5D0 * R(NR+1) * PNeI(NR+1) * ChieLP
     &                          / (RHI(NR) * PNeHI(NR+1))
                BLC(2,LQe5,NR) = - 1.5D0 * R(NR+1) * PNeI(NR+1) * ChieLP
     &                          / (RHI(NR) * PNeHI(NR))
     &                          - 1.5D0 * R(NR  ) * PNeI(NR  ) * ChieLM
     &                          / (RHI(NR) * PNeHI(NR))
                CLC(2,LQe5,NR) = + 1.5D0 * R(NR  ) * PNeI(NR  ) * ChieLM
     &                          / (RHI(NR) * PNeHI(NR-1))
            ENDIF
             NLC(2,LQe5,NR) = LQe5

! Joule heating

             ALC(3,LQe5,NR) = - AEE * EthI(NR+1) / (2 * rKeV)
             BLC(3,LQe5,NR) = - AEE * EthI(NR  ) / (2 * rKeV)
             NLC(3,LQe5,NR) = LQe3

             BLC(4,LQe5,NR) = - AEE * EphHI(NR) / rKeV
             NLC(4,LQe5,NR) = LQe4

! Collisional transfer with ions

             BLC(5,LQe5,NR) = - rNuTei(NR)
             NLC(5,LQe5,NR) = LQe5

             BLC(6,LQe5,NR) =   rNuTei(NR) * (PNeHI(NR)/PNiHI(NR))
             NLC(6,LQe5,NR) = LQi5

! Collisional heating with beam

             BLC(7,LQe5,NR) = - AMb * Vb * rNube(NR) / (2 * rKeV)
     &                       * (PNbHI(NR) / PNeHI(NR))
             NLC(7,LQe5,NR) = LQe4

             BLC(8,LQe5,NR) =   AMb * Vb * rNube(NR) / (2 * rKeV)
             NLC(8,LQe5,NR) = LQb4

! Loss to diverter

             BLC(9,LQe5,NR) = - 2.5D0 * rNuL(NR)
             NLC(9,LQe5,NR) = LQe5

             PLC(10,LQe5,NR) =  2.5D0 * rNuL(NR) * PNeHI(NR) * PTeDIV

! Direct heating (RF)

             PLC(11,LQe5,NR) =  PRFe(NR) / (1.D20 * rKeV)
         ENDDO

! Out of region

         NR = NRMAX
             BLC(1,LQe5,NR) = 1.D0
             NLC(1,LQe5,NR) = LQe5

         NLCMAX(LQe5) = 11
      ENDIF
      RETURN
      END
!
!     ***************************************************************
!
!        Ion Density Equation (HI)
!
!     ***************************************************************
!
      SUBROUTINE LQi1CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR

      DO NR = 0, NRMAX - 1
          BLC(0,LQi1,NR) = 1.D0 / DT
          NLC(0,LQi1,NR) = LQi1

! Convection

          ALC(1,LQi1,NR) = - R(NR+1) / (RHI(NR) * DR)
          BLC(1,LQi1,NR) = + R(NR  ) / (RHI(NR) * DR)
          NLC(1,LQi1,NR) = LQi2

! Ionization of n01 and n02

          BLC(2,LQi1,NR) = rNuION(NR) / PZ
     &                   * PNeHI(NR) / (PN01HI(NR) + PN02HI(NR))
          NLC(2,LQi1,NR) = LQn1

          BLC(3,LQi1,NR) = rNuION(NR) / PZ
     &                   * PNeHI(NR) / (PN01HI(NR) + PN02HI(NR))
          NLC(3,LQi1,NR) = LQn2

! Loss to divertor

          BLC(4,LQi1,NR) = - rNuL(NR) / PZ
          NLC(4,LQi1,NR) = LQe1

          PLC(5,LQi1,NR) =   rNuL(NR) * PNeDIV / PZ

! Particle source from beam ion

          BLC(6,LQi1,NR) = rNuB(NR)
          NLC(6,LQi1,NR) = LQb1

! NBI kick up ions

         PLC(7,LQi1,NR) = - SNB(NR)

! Loss cone loss

         PLC(8,LQi1,NR) = SiLC(NR)

      ENDDO

!     OUT OF REGION

      NR = NRMAX
          BLC(1,LQi1,NR) = 1.D0
          NLC(1,LQi1,NR) = LQi1

      NLCMAX(LQi1) = 8
      RETURN
      END
!
!     ***************************************************************
!
!        Ion Radial Flow (NR)
!
!     ***************************************************************
!
      SUBROUTINE LQi2CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR

! Ns*Usr(0) : fixed

      NR = 0
          BLC(1,LQi2,NR) = 1.D0
          NLC(1,LQi2,NR) = LQi2

      DO NR = 1, NRMAX-1

          BLC(0,LQi2,NR) = 1.D0 / DT
          NLC(0,LQi2,NR) = LQi2

! Nonlinear term

          ALC(1,LQi2,NR) = - RHI(NR  ) * UirHI(NR  ) / (2 * R(NR) * DR)
          BLC(1,LQi2,NR) = - RHI(NR  ) * UirHI(NR  ) / (2 * R(NR) * DR)
     &                    + RHI(NR-1) * UirHI(NR-1) / (2 * R(NR) * DR)
          CLC(1,LQi2,NR) = + RHI(NR-1) * UirHI(NR-1) / (2 * R(NR) * DR)
          NLC(1,LQi2,NR) = LQi2

! Nonlinear centrifugal force

          BLC(2,LQi2,NR) = UithI(NR) / R(NR)
          NLC(2,LQi2,NR) = LQi3

! Pressure gradient force

          BLC(3,LQi2,NR) = - rKEV / (AMI * DR)
          CLC(3,LQi2,NR) = + rKEV / (AMI * DR)
          NLC(3,LQi2,NR) = LQi5

! Radial E force

          BLC(4,LQi2,NR) =   PZ * PNiI(NR) * (AEE/AMI)
          NLC(4,LQi2,NR) = LQm1

! v x B force

          BLC(5,LQi2,NR) =   PZ * BphI(NR) * (AEE/AMI)
          NLC(5,LQi2,NR) = LQi3

          BLC(6,LQi2,NR) = - PZ * BthI(NR) * (AEE/AMI) / 2
          CLC(6,LQi2,NR) = - PZ * BthI(NR) * (AEE/AMI) / 2
          NLC(6,LQi2,NR) = LQi4
      ENDDO

! Ns*Usr(NRMAX) : fixed or finite gradient

      NR = NRMAX
          BLC(1,LQi2,NR) = 1.D0
          NLC(1,LQi2,NR) = LQi2

      NLCMAX(LQi2) = 6
      RETURN
      END
!
!     ***************************************************************
!
!        Ion Poloidal Flow (NR)
!
!     ***************************************************************
!
      SUBROUTINE LQi3CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR
      REAL(8) :: rNuiNCL, rNueiL, rNubiL, FWtheL, WPML, TMP, FWthiL,
     &           rNuLL, rNu0iL, rNuiCXL, SiLCthL, rNuiHLL

! Ni*UiTheta(0) : 0

      NR = 0
          BLC(1,LQi3,NR) = 1.D0
          NLC(1,LQi3,NR) = LQi3

      DO NR = 1, NRMAX - 1

          BLC(0,LQi3,NR) = 1.D0 / DT
          NLC(0,LQi3,NR) = LQi3

! Nonlinear term

          ALC(1,LQi3,NR) = - RHI(NR  )**2 * UirHI(NR  )
     &                    / (2 * R(NR)**2 * DR)
          BLC(1,LQi3,NR) = - RHI(NR  )**2 * UirHI(NR  )
     &                    / (2 * R(NR)**2 * DR)
     &                    + RHI(NR-1)**2 * UirHI(NR-1)
     &                    / (2 * R(NR)**2 * DR)
          CLC(1,LQi3,NR) = + RHI(NR-1)**2 * UirHI(NR-1)
     &                    / (2 * R(NR)**2 * DR)
          NLC(1,LQi3,NR) = LQi3

! Viscosity force

         IF (NR .EQ. 1) THEN
             ALC(2,LQi3,NR) =   RHI(NR  )**3 * PNiHI(NR  ) * rMui(NR  )
     &                      / (PNiI(NR+1) * R(NR+1) * R(NR)**2 * DR**2)
             BLC(2,LQi3,NR) = - RHI(NR  )**3 * PNiHI(NR  ) * rMui(NR  )
     &                      / (PNiI(NR  ) * R(NR  ) * R(NR)**2 * DR**2)
         ELSE
             ALC(2,LQi3,NR) =   RHI(NR  )**3 * PNiHI(NR  ) * rMui(NR  )
     &                      / (PNiI(NR+1) * R(NR+1) * R(NR)**2 * DR**2)
             BLC(2,LQi3,NR) = - RHI(NR  )**3 * PNiHI(NR  ) * rMui(NR  )
     &                      / (PNiI(NR  ) * R(NR  ) * R(NR)**2 * DR**2)
     &                       - RHI(NR-1)**3 * PNiHI(NR-1) * rMui(NR-1)
     &                      / (PNiI(NR  ) * R(NR  ) * R(NR)**2 * DR**2)
             CLC(2,LQi3,NR) =   RHI(NR-1)**3 * PNiHI(NR-1) * rMui(NR-1)
     &                      / (PNiI(NR-1) * R(NR-1) * R(NR)**2 * DR**2)
         ENDIF
          NLC(2,LQi3,NR) = LQi3

! Poroidal E force

          BLC(3,LQi3,NR) = PZ * PNiI(NR  ) * (AEE/AMI)
          NLC(3,LQi3,NR) = LQm2

! v x B force

          BLC(4,LQi3,NR) = - PZ * BphI(NR) * (AEE/AMI)
          NLC(4,LQi3,NR) = LQi2

! Neoclassical viscosity force

          rNuiNCL = 0.5D0 * (rNuiNC(NR) + rNuiNC(NR-1))
          BLC(5,LQi3,NR) = - rNuiNCL
          NLC(5,LQi3,NR) = LQi3

! Collisional friction force

          rNueiL = 0.5D0 * (rNuei(NR) + rNuei(NR-1))
          BLC(6,LQi3,NR) = - (AME/AMI) * rNueiL * PNeI(NR) / PNiI(NR)
          NLC(6,LQi3,NR) = LQi3

          BLC(7,LQi3,NR) =   (AME/AMI) * rNueiL
          NLC(7,LQi3,NR) = LQe3

! Collisional friction with beam ions

          rNubiL = 0.5D0 * (rNubi(NR) + rNubi(NR-1))
          BLC(8,LQi3,NR) = - (AMB/AMI) * rNubiL * PNbI(NR) / PNiI(NR)
          NLC(8,LQi3,NR) = LQi3

          BLC(9,LQi3,NR) =   (AMB/AMI) * rNubiL
          NLC(9,LQi3,NR) = LQb3

! Wave interaction force (electron driven)

          FWtheL = 0.5D0 * (FWthe(NR) + FWthe(NR-1))
          WPML = 0.5D0 * (  WPM(NR) +   WPM(NR-1))
          BLC(10,LQi3,NR) = FWtheL / AMI
          NLC(10,LQi3,NR) = LQe3

          TMP = - FWtheL * WPML * R(NR) / AMI
          BLC(11,LQi3,NR) = TMP / 2
          CLC(11,LQi3,NR) = TMP / 2
          NLC(11,LQi3,NR) = LQe1

! Wave interaction force (NRon driven)

          FWthiL = 0.5D0 * (FWthi(NR) + FWthi(NR-1))
          BLC(12,LQi3,NR) = - FWthiL / AMI
          NLC(12,LQi3,NR) = LQi3

          TMP = FWthiL * WPML * R(NR) / AMI
          BLC(13,LQi3,NR) = TMP / 2
          CLC(13,LQi3,NR) = TMP / 2
          NLC(13,LQi3,NR) = LQi1

! Loss to divertor

          rNuLL = 0.5D0 * (rNuL(NR) + rNuL(NR-1))
          BLC(14,LQi3,NR) = - 2.D0 * rNuLL
          NLC(14,LQi3,NR) = LQi3

! Collisional friction force with neutrals

          rNu0iL = 0.5D0 * (rNu0i(NR) + rNu0i(NR-1))
          BLC(15,LQi3,NR) = - rNu0iL
          NLC(15,LQi3,NR) = LQi3

! Charge exchange force

          rNuiCXL = 0.5D0 * (rNuiCX(NR) + rNuiCX(NR-1))
          BLC(16,LQi3,NR) = - rNuiCXL
          NLC(16,LQi3,NR) = LQi3

! Loss cone loss

          SiLCthL = 0.5D0 * (SiLCth(NR) + SiLCth(NR-1))
          PLC(17,LQi3,NR) = + SiLCthL

! Helical Neoclassical viscosity force

          rNuiHLL = rNuiHL(NR)
          BLC(18,LQi3,NR) = - rNuiHLL*(1.D0-UHth*UHth)
          NLC(18,LQi3,NR) = LQi3
          BLC(19,LQi3,NR) =   0.5D0*rNuiHLL*UHph*UHth
          CLC(19,LQi3,NR) =   0.5D0*rNuiHLL*UHph*UHth
          NLC(19,LQi3,NR) = LQi4
      ENDDO

! Ns*UsTheta(NRMAX) : 0

      NR = NRMAX
          BLC(1,LQi3,NR) = 1.D0
          NLC(1,LQi3,NR) = LQi3

      NLCMAX(LQi3) = 19
      RETURN
      END
!
!     ***************************************************************
!
!        Ion Toroidal Flow (NR)
!
!     ***************************************************************
!
      SUBROUTINE LQi4CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR
      REAL(8) :: rMuiP, rMuiM, rNuiHLL

! Uiphi'(0) : 0

      DO NR = 0, NRMAX - 1
          BLC(0,LQi4,NR) = 1.D0 / DT
          NLC(0,LQi4,NR) = LQi4

! Nonlinear term

          ALC(1,LQi4,NR) = - R(NR+1) * UirI(NR+1) / (2 * RHI(NR) * DR)
          BLC(1,LQi4,NR) = - R(NR+1) * UirI(NR+1) / (2 * RHI(NR) * DR)
     &                     + R(NR  ) * UirI(NR  ) / (2 * RHI(NR) * DR)
          CLC(1,LQi4,NR) = + R(NR  ) * UirI(NR  ) / (2 * RHI(NR) * DR)
          NLC(1,LQi4,NR) = LQi4

! Viscosity force

          IF(NR.EQ.0) THEN
             rMuiP = 0.5D0*(rMui(NR) + rMui(NR+1))
             rMuiM =        rMui(NR)
             ALC(2,LQi4,NR) =   R(NR+1) * PNiI(NR+1) * rMuiP
     &                    / (PNiHI(NR+1) * RHI(NR) * DR**2)
             BLC(2,LQi4,NR) = - R(NR+1) * PNiI(NR+1) * rMuiP
     &                    / (PNiHI(NR  ) * RHI(NR) * DR**2)
     &                       - R(NR  ) * PNiI(NR  ) * rMuiM
     &                    / (PNiHI(NR  ) * RHI(NR) * DR**2)
     &                       + R(NR  ) * PNiI(NR  ) * rMuiM
     &                    / (PNiHI(NR  ) * RHI(NR) * DR**2)
          ELSEIF(NR.EQ.NRMAX-1) THEN
             rMuiP =        rMui(NR)
             rMuiM = 0.5D0*(rMui(NR-1) + rMui(NR))
             BLC(2,LQi4,NR) =   R(NR+1) * PNiI(NR+1) * rMuiP
     &                    / (PNiHI(NR  ) * RHI(NR) * DR**2)
     &                       - R(NR+1) * PNiI(NR+1) * rMuiP
     &                    / (PNiHI(NR  ) * RHI(NR) * DR**2)
     &                       - R(NR  ) * PNiI(NR  ) * rMuiM
     &                    / (PNiHI(NR  ) * RHI(NR) * DR**2)
             CLC(2,LQi4,NR) =   R(NR  ) * PNiI(NR  ) * rMuiM
     &                    / (PNiHI(NR-1) * RHI(NR) * DR**2)
          ELSE
             rMuiP = 0.5D0*(rMui(NR) + rMui(NR+1))
             rMuiM = 0.5D0*(rMui(NR-1) + rMui(NR))
             ALC(2,LQi4,NR) =   R(NR+1) * PNiI(NR+1) * rMuiP
     &                    / (PNiHI(NR+1) * RHI(NR) * DR**2)
             BLC(2,LQi4,NR) = - R(NR+1) * PNiI(NR  ) * rMuiP
     &                    / (PNiHI(NR  ) * RHI(NR) * DR**2)
     &                       - R(NR  ) * PNiI(NR-1) * rMuiM
     &                    / (PNiHI(NR  ) * RHI(NR) * DR**2)
             CLC(2,LQi4,NR) =   R(NR  ) * PNiI(NR-1) * rMuiM
     &                    / (PNiHI(NR-1) * RHI(NR) * DR**2)
          ENDIF
          NLC(2,LQi4,NR) = LQi4

! Toroidal E force

          BLC(3,LQi4,NR) = PZ * PNiHI(NR  ) * (AEE/AMI)
          NLC(3,LQi4,NR) = LQm3

! v x B force

          ALC(4,LQi4,NR) = PZ * BthI(NR+1) * (AEE/AMI) / 2
          BLC(4,LQi4,NR) = PZ * BthI(NR  ) * (AEE/AMI) / 2 
          NLC(4,LQi4,NR) = LQi2

! Collisional friction with bulk ions

          BLC(5,LQi4,NR) = - (AME/AMI) * rNuei(NR) * PNeI(NR) / PNiI(NR)
          NLC(5,LQi4,NR) = LQi4

          BLC(6,LQi4,NR) =   (AME/AMI) * rNuei(NR)
          NLC(6,LQi4,NR) = LQe4

! Collisional friction with beam ions

          BLC(7,LQi4,NR) = - (AMB/AMI) * rNubi(NR) * PNbI(NR) / PNiI(NR)
          NLC(7,LQi4,NR) = LQi4

          BLC(8,LQi4,NR) =   (AMB/AMI) * rNubi(NR)
          NLC(8,LQi4,NR) = LQb4

! Loss to divertor

          BLC(9,LQi4,NR) = - 2.D0 * rNuL(NR)
          NLC(9,LQi4,NR) = LQi4

! Collisional friction force with neutrals

          BLC(10,LQi4,NR) = - rNu0i(NR)
          NLC(10,LQi4,NR) = LQi4

! Charge exchange force

          BLC(11,LQi4,NR) = - rNuiCX(NR)
          NLC(11,LQi4,NR) = LQi4

! Loss conde loss

          PLC(12,LQi4,NR) = SiLCph(NR)

! Helical Neoclassical viscosity force

          rNuiHLL = 0.5D0 * (rNuiHL(NR) + rNuiHL(NR+1))
          ALC(13,LQi4,NR) =   rNuiHL(NR+1)*UHth*UHph/2
          BLC(13,LQi4,NR) =   rNuiHL(NR  )*UHth*UHph/2
          NLC(13,LQi4,NR) = LQi3
          BLC(14,LQi4,NR) = - rNuiHLL*(1.D0-UHph*UHph)
          NLC(14,LQi4,NR) = LQi4

      ENDDO

! Out of region

      NR = NRMAX
          BLC(1,LQi4,NR) = 1.D0
          NLC(1,LQi4,NR) = LQi4

      NLCMAX(LQi4) = 14
      RETURN
      END
!
!     ***************************************************************
!
!        Ion Energy Transport: Ti (HI)
!
!     ***************************************************************
!
      SUBROUTINE LQi5CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR
      REAL(8) :: ChiiLP, ChiiLM

! Fixed temperature

      IF (ABS(Chii0) .EQ. 0.D0) THEN
         DO NR = 0, NRMAX
             BLC(0,LQi5,NR) = 1.D0
             NLC(0,LQi5,NR) = LQi5
         ENDDO
         NLCMAX(LQi5) = 0

! Temperature evolution

      ELSE
         DO NR = 0, NRMAX-1
             BLC(0,LQi5,NR) = 1.5D0 / DT
             NLC(0,LQi5,NR) = LQi5

! Convection transport

            IF (NR.EQ.0) THEN
               ALC(1,LQi5,NR) = - 2.5D0 * R(NR+1) * UirI(NR+1)
     &                            / (2 * RHI(NR) * DR)
               BLC(1,LQi5,NR) = - 2.5D0 * R(NR+1) * UirI(NR+1)
     &                            / (2 * RHI(NR) * DR)
            ELSEIF (NR.EQ.NRMAX-1) THEN
               BLC(1,LQi5,NR) = - 2.5D0 * R(NR+1) * UirI(NR+1)
     &                            / (2 * RHI(NR) * DR)
     &                         - 2.5D0 * R(NR+1) * UirI(NR+1)
     &                            / (2 * RHI(NR) * DR)
     &                         + 2.5D0 * R(NR  ) * UirI(NR  )
     &                            / (2 * RHI(NR) * DR)
               CLC(1,LQi5,NR) = + 2.5D0 * R(NR  ) * UirI(NR  )
     &                            / (2 * RHI(NR) * DR)
            ELSE
               ALC(1,LQi5,NR) = - 2.5D0 * R(NR+1) * UirI(NR+1)
     &                            / (2 * RHI(NR) * DR)
               BLC(1,LQi5,NR) = - 2.5D0 * R(NR+1) * UirI(NR+1)
     &                            / (2 * RHI(NR) * DR)
     &                         + 2.5D0 * R(NR  ) * UirI(NR  )
     &                            / (2 * RHI(NR) * DR)
               CLC(1,LQi5,NR) = + 2.5D0 * R(NR  ) * UirI(NR  )
     &                            / (2 * RHI(NR) * DR)
            ENDIF
             NLC(1,LQi5,NR) = LQi5

! Conduction transport

            IF (NR.EQ.0) THEN
               ChiiLP = 0.5D0 * (Chii(NR) + Chii(NR+1))
                ALC(2,LQi5,NR) = + 1.5D0 * R(NR+1) * PNiI(NR+1) * ChiiLP
     &                          / (RHI(NR) * PNiHI(NR+1))
                BLC(2,LQi5,NR) = - 1.5D0 * R(NR+1) * PNiI(NR+1) * ChiiLP
     &                          / (RHI(NR) * PNiHI(NR))
            ELSEIF (NR.EQ.NRMAX-1) THEN
               ChiiLP =          Chii(NR) 
               ChiiLM = 0.5D0 * (Chii(NR-1) + Chii(NR))
                BLC(2,LQi5,NR) = + 1.5D0 * R(NR+1) * PNiI(NR+1) * ChiiLP
     &                          / (RHI(NR) * PNiHI(NR))
     &                          - 1.5D0 * R(NR+1) * PNiI(NR+1) * ChiiLP
     &                          / (RHI(NR) * PNiHI(NR))
     &                          - 1.5D0 * R(NR  ) * PNiI(NR  ) * ChiiLM
     &                          / (RHI(NR) * PNiHI(NR))
                CLC(2,LQi5,NR) = + 1.5D0 * R(NR  ) * PNiI(NR  ) * ChiiLM
     &                          / (RHI(NR) * PNiHI(NR-1))
            ELSE
               ChiiLP = 0.5D0 * (Chii(NR) + Chii(NR+1))
               ChiiLM = 0.5D0 * (Chii(NR-1) + Chii(NR))
                ALC(2,LQi5,NR) = + 1.5D0 * R(NR+1) * PNiI(NR+1) * ChiiLP
     &                          / (RHI(NR) * PNiHI(NR+1))
                BLC(2,LQi5,NR) = - 1.5D0 * R(NR+1) * PNiI(NR+1) * ChiiLP
     &                          / (RHI(NR) * PNiHI(NR))
     &                          - 1.5D0 * R(NR  ) * PNiI(NR  ) * ChiiLM
     &                          / (RHI(NR) * PNiHI(NR))
                CLC(2,LQi5,NR) = + 1.5D0 * R(NR  ) * PNiI(NR  ) * ChiiLM
     &                          / (RHI(NR) * PNiHI(NR-1))
            ENDIF
             NLC(2,LQi5,NR) = LQi5

! Joule heating

             ALC(3,LQi5,NR) =   PZ * AEE * EthI(NR+1) / (2 * rKeV)
             BLC(3,LQi5,NR) =   PZ * AEE * EthI(NR  ) / (2 * rKeV)
             NLC(3,LQi5,NR) = LQi3

             BLC(4,LQi5,NR) =   PZ * AEE * EphHI(NR) / rKeV
             NLC(4,LQi5,NR) = LQi4

! Collisional transfer with electrons

             BLC(5,LQi5,NR) = - rNuTei(NR) * PNeHI(NR) / PNiHI(NR)
             NLC(5,LQi5,NR) = LQi5

             BLC(6,LQi5,NR) =   rNuTei(NR)
             NLC(6,LQi5,NR) = LQe5

! Collisional heating with beam

             BLC(7,LQi5,NR) = - AMb * Vb * rNubi(NR) / (2 * rKeV)
     &                       * (PNbHI(NR) / PNiHI(NR))
             NLC(7,LQi5,NR) = LQi4

             BLC(8,LQi5,NR) =   AMb * Vb * rNubi(NR) / (2 * rKeV)
             NLC(8,LQi5,NR) = LQb4

! Loss to diverter

             BLC(9,LQi5,NR) = - 2.5D0 * rNuL(NR)
             NLC(9,LQi5,NR) = LQi5

             PLC(10,LQi5,NR) =  2.5D0 * rNuL(NR) * PNiHI(NR) * PTiDIV

! Direct heating (RF)

             PLC(11,LQi5,NR) =  PRFi(NR) / (1.D20 * rKeV)
         ENDDO

! Out of region

         NR = NRMAX
             BLC(1,LQi5,NR) = 1.D0
             NLC(1,LQi5,NR) = LQi5

         NLCMAX(LQi5) = 11
      ENDIF
      RETURN
      END
!
!     ***************************************************************
!
!        Beam Ion Density (HI)
!
!     ***************************************************************
!
      SUBROUTINE LQb1CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR

      DO NR = 0, NRMAX - 1
          BLC(0,LQb1,NR) = 1.D0 / DT
          NLC(0,LQb1,NR) = LQb1

! NBI particle source

          PLC(1,LQb1,NR) = SNB(NR)

! Relaxation to thermal ions

          BLC(2,LQb1,NR) = - rNuB(NR)
          NLC(2,LQb1,NR) = LQb1
      ENDDO

! Out of region

      NR = NRMAX
          BLC(1,LQb1,NR) = 1.D0
          NLC(1,LQb1,NR) = LQb1

      NLCMAX(LQb1) = 2
      RETURN
      END
!
!     ***************************************************************
!
!        Beam Ion Poloidal Flow (NR)
!
!     ***************************************************************
!
      SUBROUTINE LQb3CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR

! Ubth(0) : 0

      NR = 0
          BLC(1,LQb3,NR) = 1.D0
          NLC(1,LQb3,NR) = LQb3

      DO NR = 1, NRMAX-1

! Nonlinear centrifugal force

          BLC(1,LQb3,NR) = UbthI(NR) / R(NR)
          NLC(1,LQb3,NR) = LQb3

! Radial E force

          BLC(2,LQb3,NR) = PZ * PNbI(NR) * (AEE/AMB)
          NLC(2,LQb3,NR) = LQm1

! v x B force

          BLC(3,LQb3,NR) = PZ * BphI(NR) * (AEE/AMB)
          NLC(3,LQb3,NR) = LQb3

          BLC(4,LQb3,NR) = - PZ * BthI(NR) * (AEE/AMB)
          NLC(4,LQb3,NR) = LQb4
      ENDDO

! Ubth(NRMAX) : 0

      NR = NRMAX
          BLC(1,LQb3,NR) = 1.D0
          NLC(1,LQb3,NR) = LQb3

      NLCMAX(LQb3) = 4
      RETURN
      END
!
!     ***************************************************************
!
!        Beam Ion Toroildal Flow (NR)
!
!     ***************************************************************
!
      SUBROUTINE LQb4CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR
      REAL(8) :: SNBL

! Ubphi'(0) : 0

      NR = 0
          IF (PNbI(0)*PNbI(1) .EQ. 0.D0) THEN
             ALC(1,LQb4,NR) = -1.D0
             BLC(1,LQb4,NR) =  1.D0
          ELSE
             ALC(1,LQb4,NR) = -1.D0 / PNbI(NR+1)
             BLC(1,LQb4,NR) =  1.D0 / PNbI(NR  )
          ENDIF
          NLC(1,LQb4,NR) = LQb4

      DO NR = 1, NRMAX - 1
          BLC(0,LQb4,NR) = 1.D0 / DT
          NLC(0,LQb4,NR) = LQb4

! Collisional friction with electrons

          BLC(1,LQb4,NR) = - rNube(NR)
          NLC(1,LQb4,NR) = LQb4

          BLC(2,LQb4,NR) =   rNube(NR) * PNbI(NR) / PNeI(NR)
          NLC(2,LQb4,NR) = LQe4

! Collisional friction with ions

          BLC(3,LQb4,NR) = - rNubi(NR)
          NLC(3,LQb4,NR) = LQb4

          BLC(4,LQb4,NR) =   rNubi(NR)
          NLC(4,LQb4,NR) = LQi4

! NBI momentum source

         SNBL = 0.5D0 * (SNB(NR) + SNB(NR-1))
         PLC(5,LQb4,NR) = PNBCD * Vb * SNBL
      ENDDO

! Ubphi(NRMAX) : 0

      NR = NRMAX
          BLC(1,LQb4,NR) = 1.D0
          NLC(1,LQb4,NR) = LQb4

      NLCMAX(LQb4) = 5
      RETURN
      END
!
!     ***************************************************************
!
!        Slow Neutral Transport: n01 (HI)
!
!     ***************************************************************
!
      SUBROUTINE LQn1CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR
      REAL(8) :: D0LP, D0LM, rNuiCXL

      DO NR = 0, NRMAX-1
          BLC(0,LQn1,NR) = 1.D0 / DT
          NLC(0,LQn1,NR) = LQn1

!  Diffusion of neutrals

         IF (NR.EQ.0) THEN
             D0LP = 0.5D0 * (D01(NR) + D01(NR+1))
             ALC(1,LQn1,NR) =   D0LP * R(NR+1) / (RHI(NR)*DR*DR)
             BLC(1,LQn1,NR) = - D0LP * R(NR+1) / (RHI(NR)*DR*DR)
         ELSEIF (NR.EQ.NRMAX-1) THEN
             D0LM = 0.5D0 * (D01(NR-1) + D01(NR))
             BLC(1,LQn1,NR) = - D0LM * R(NR  ) / (RHI(NR)*DR*DR)
             CLC(1,LQn1,NR) =   D0LM * R(NR  ) / (RHI(NR)*DR*DR)
             PLC(1,LQn1,NR) =          R(NR+1) / (RHI(NR)*DR)
     &                     * rGASPF
         ELSE
             D0LP = 0.5D0 * (D01(NR) + D01(NR+1))
             D0LM = 0.5D0 * (D01(NR-1) + D01(NR))
             ALC(1,LQn1,NR) =   D0LP * R(NR+1) / (RHI(NR)*DR*DR)
             BLC(1,LQn1,NR) = - D0LP * R(NR+1) / (RHI(NR)*DR*DR)
     &                       - D0LM * R(NR  ) / (RHI(NR)*DR*DR)
             CLC(1,LQn1,NR) =   D0LM * R(NR  ) / (RHI(NR)*DR*DR)
         ENDIF
          NLC(1,LQn1,NR) = LQn1

! Ionization

          BLC(2,LQn1,NR) = - rNuION(NR) / PZ
     &                    * PNeHI(NR) / (PN01HI(NR) + PN02HI(NR))
          NLC(2,LQn1,NR) = LQn1

! Generation of fast neutrals by charge exchange

         IF(NR.NE.0) THEN
             rNuiCXL = ( rNuiCX(NR  )/(PN01HI(NR  ) + PN02HI(NR  ))
     &                 + rNuiCX(NR-1)/(PN01HI(NR-1) + PN02HI(NR-1)))/2
         ELSE
             rNuiCXL =   rNuiCX(NR  )/(PN01HI(NR  ) + PN02HI(NR  ))
         ENDIF
          BLC(3,LQn1,NR) = - rNuiCXL*PNiHI(NR)
          NLC(3,LQn1,NR) = LQn1

! Recycling from divertor

          BLC(4,LQn1,NR) =   rGamm0 * rNuL(NR) / PZ
          NLC(4,LQn1,NR) = LQe1

          PLC(5,LQn1,NR) = - rGamm0 * rNuL(NR) * PneDIV / PZ

      ENDDO

! Out of region

      NR = NRMAX
          BLC(1,LQn1,NR) = 1.D0
          NLC(1,LQn1,NR) = LQn1

      NLCMAX(LQn1) = 5

      RETURN
      END
!
!     ***************************************************************
!
!        Fast Neutral Transport: n02 (HI)
!
!     ***************************************************************
!
      SUBROUTINE LQn2CC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR
      REAL(8) :: D0LP, D0LM, rNuiCXL

      DO NR = 0, NRMAX-1
          BLC(0,LQn2,NR) = 1.D0 / DT
          NLC(0,LQn2,NR) = LQn2

!  Diffusion of neutrals

         IF (NR.EQ.0) THEN
             D0LP = 0.5D0 * (D02(NR) + D02(NR+1))
             ALC(1,LQn2,NR) =   D0LP * R(NR+1) / (RHI(NR)*DR*DR)
             BLC(1,LQn2,NR) = - D0LP * R(NR+1) / (RHI(NR)*DR*DR)
         ELSEIF (NR.EQ.NRMAX-1) THEN
             D0LM = 0.5D0 * (D02(NR-1) + D02(NR))
             BLC(1,LQn2,NR) = - D0LM * R(NR  ) / (RHI(NR)*DR*DR)
             CLC(1,LQn2,NR) =   D0LM * R(NR  ) / (RHI(NR)*DR*DR)
         ELSE
             D0LP = 0.5D0 * (D02(NR) + D02(NR+1))
             D0LM = 0.5D0 * (D02(NR-1) + D02(NR))
             ALC(1,LQn2,NR) =   D0LP * R(NR+1) / (RHI(NR)*DR*DR)
             BLC(1,LQn2,NR) = - D0LP * R(NR+1) / (RHI(NR)*DR*DR)
     &                       - D0LM * R(NR  ) / (RHI(NR)*DR*DR)
             CLC(1,LQn2,NR) =   D0LM * R(NR  ) / (RHI(NR)*DR*DR)
         ENDIF
          NLC(1,LQn2,NR) = LQn2

! Ionization

          BLC(2,LQn2,NR) = - rNuION(NR) / PZ
     &                    * PNeHI(NR) / (PN01HI(NR) + PN02HI(NR))
          NLC(2,LQn2,NR) = LQn2

! Generation of fast neutrals by charge exchange

         IF(NR.NE.0) THEN
             rNuiCXL = ( rNuiCX(NR  )/(PN01HI(NR  ) + PN02HI(NR  ))
     &                 + rNuiCX(NR-1)/(PN01HI(NR-1) + PN02HI(NR-1)))/2
         ELSE
             rNuiCXL =   rNuiCX(NR  )/(PN01HI(NR  ) + PN02HI(NR  ))
         ENDIF
          BLC(3,LQn2,NR) = rNuiCXL * PNiHI(NR)
          NLC(3,LQn2,NR) = LQn1

! NBI particle source

          PLC(4,LQn2,NR) = SNB(NR)

      ENDDO

! Out of region

      NR = NRMAX
          BLC(1,LQn2,NR) = 1.D0
          NLC(1,LQn2,NR) = LQn2

      NLCMAX(LQn2) = 4

      RETURN
      END
