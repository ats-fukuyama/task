C     $Id$
C
C     ***************************************************************
C
C        Calculate ALC, BLC, CLC, PLC
C
C     ***************************************************************
C
      SUBROUTINE TXCALA
C
      INCLUDE 'txcomm.inc'
C
      DO NR = 0, NRMAX
         DO NQ = 1, NQMAX
            DO NC = 0, NCM
                ALC(NC,NQ,NR) = 0.D0
                BLC(NC,NQ,NR) = 0.D0
                CLC(NC,NQ,NR) = 0.D0
                NLC(NC,NQ,NR) = 1
            ENDDO
         ENDDO
      ENDDO
C
      DO NR = 0, NRMAX
         DO NQ = 1, NQMAX
            DO NC = 1, NCM
               PLC(NC,NQ,NR) = 0.D0
            ENDDO
         ENDDO
      ENDDO
C
      DO NQ = 1, NQMAX
         NLCMAX(NQ) = 0
      ENDDO
C
C     Maxwell
C
      CALL LQm1CC
      CALL LQm2CC
      CALL LQm3CC
      CALL LQm4CC
      CALL LQm5CC
C
C     Electron
C
      CALL LQe1CC
      CALL LQe2CC
      CALL LQe3CC
      CALL LQe4CC
      CALL LQe5CC
C
C     Ion
C
      CALL LQi1CC
      CALL LQi2CC
      CALL LQi3CC
      CALL LQi4CC
      CALL LQi5CC
C
C     Beam
C
      CALL LQb1CC
      CALL LQb3CC
      CALL LQb4CC
C
C     Neutral
C
      CALL LQn1CC
      CALL LQn2CC
C
      RETURN
      END
C
C     ***************************************************************
C
C        Poisson Equation  (HI)
C
C     ***************************************************************
C
      SUBROUTINE LQm1CC
C
      INCLUDE 'txcomm.inc'
C
C Er(0) = 0
C
      NR = 0
          BLC(1,LQm1,NR) = 1.D0
          NLC(1,LQm1,NR) = LQm1
C
      FACTOR = EPS0 / (AEE * 1.D20)
C
CCC      DO NR = 1, NRMAX-1
      DO NR = 1, NRMAX
          BLC(1,LQm1,NR) =  FACTOR * R(NR  ) / (DR * RHI(NR-1))
          CLC(1,LQm1,NR) = -FACTOR * R(NR-1) / (DR * RHI(NR-1))
          NLC(1,LQm1,NR) = LQm1
C
          CLC(2,LQm1,NR) = 1.D0
          NLC(2,LQm1,NR) = LQe1
C
          CLC(3,LQm1,NR) = - PZ
          NLC(3,LQm1,NR) = LQi1
C
          CLC(4,LQm1,NR) = - PZ
          NLC(4,LQm1,NR) = LQb1
      ENDDO
CCC
CCC      NR = NRMAX
CCC          BLC(1,LQm1,NR) = 1.D0
CCC          NLC(1,LQm1,NR) = LQm1
C
      NLCMAX(LQm1) = 4
      RETURN
      END
C
C     ***************************************************************
C
C        Ampere's Law: Eth (NR)
C
C     ***************************************************************
C
      SUBROUTINE LQm2CC
C
      INCLUDE 'txcomm.inc'
C
C Etheta(0) : 0
C
      NR = 0
          BLC(1,LQm2,NR) = 1.D0
          NLC(1,LQm2,NR) = LQm2
C
      DO NR = 1, NRMAX - 1
          BLC(0,LQm2,NR) = 1.D0 / (VC**2 * DT)
          NLC(0,LQm2,NR) = LQm2
C
C rot Bphi
C
          BLC(1,LQm2,NR) = - 1.D0 / DR
          CLC(1,LQm2,NR) =   1.D0 / DR
          NLC(1,LQm2,NR) = LQm5
C
C Electron current
C
          BLC(2,LQm2,NR) =   rMU0 *      AEE * 1.D20 
          NLC(2,LQm2,NR) = LQe3
C
C Ion current
C
          BLC(3,LQm2,NR) = - rMU0 * PZ * AEE * 1.D20
          NLC(3,LQm2,NR) = LQi3
C
C Beam ion current
C
          BLC(4,LQm2,NR) = - rMU0 * PZ * AEE * 1.D20
          NLC(4,LQm2,NR) = LQb3
      ENDDO
C
C Etheta(NRMAX) : r * Etheta = const.
C
      NR = NRMAX
          BLC(1,LQm2,NR) = R(NR)
          CLC(1,LQm2,NR) =-R(NR-1)
          NLC(1,LQm2,NR) = LQm2
C
      NLCMAX(LQm2) = 4
      RETURN
      END
C
C     ***************************************************************
C
C        Ampere's Law: Ephi (HI)
C
C     ***************************************************************
C
      SUBROUTINE LQm3CC
C
      INCLUDE 'txcomm.inc'
C
      DO NR = 0, NRMAX - 1
          BLC(0,LQm3,NR) = 1.D0 / (VC**2 * DT)
          NLC(0,LQm3,NR) = LQm3
C
C rot Btheta
C
          ALC(1,LQm3,NR) =   R(NR+1) / (DR * RHI(NR))
          BLC(1,LQm3,NR) = - R(NR  ) / (DR * RHI(NR))
          NLC(1,LQm3,NR) = LQm4
C
C Electron current
C
          BLC(2,LQm3,NR) =   rMU0 *      AEE * 1.D20
          NLC(2,LQm3,NR) = LQe4
C
C Ion current
C
          BLC(3,LQm3,NR) = - rMU0 * PZ * AEE * 1.D20
          NLC(3,LQm3,NR) = LQi4
C
C Beam ion current
C
          BLC(4,LQm3,NR) = - rMU0 * PZ * AEE * 1.D20
          NLC(4,LQm3,NR) = LQb4
      ENDDO
C
C Out of region
C
      NR = NRMAX
          BLC(1,LQm3,NR) = 1.D0
          NLC(1,LQm3,NR) = LQm3
C
      NLCMAX(LQm3) = 4
      RETURN
      END
C
C     **************************************************************
C
C        Faraday's Law : Btheta (NR)
C
C     ***************************************************************
C
      SUBROUTINE LQm4CC
C
      INCLUDE 'txcomm.inc'
C
C Btheta(0) : 0
C
      NR = 0
          BLC(1,LQm4,NR) = 1.D0
          NLC(1,LQm4,NR) = LQm4
C
      DO NR = 1, NRMAX - 1
          BLC(0,LQm4,NR) = 1.D0 / DT
          NLC(0,LQm4,NR) = LQm4
C
C rot Ephi
C
          BLC(1,LQm4,NR) =   1.D0 / DR
          CLC(1,LQm4,NR) = - 1.D0 / DR
          NLC(1,LQm4,NR) = LQm3
      ENDDO
C
C Btheta(NRMAX) : fixed
C
      NR = NRMAX
          BLC(1,LQm4,NR) = 1.D0
          NLC(1,LQm4,NR) = LQm4
C
          PLC(2,LQm4,NR) = - Bthb
C
      NLCMAX(LQm4) = 2
      RETURN
      END
C
C     ***************************************************************
C
C        Fraday's Law : Bphi (HI)
C
C     ***************************************************************
C
      SUBROUTINE LQm5CC
C
      INCLUDE 'txcomm.inc'
C
      DO NR = 0, NRMAX - 1
          BLC(0,LQm5,NR) = 1.D0 / DT
          NLC(0,LQm5,NR) = LQm5
C
C rot Etheta
C
          ALC(1,LQm5,NR) = - R(NR+1) / (DR * RHI(NR))
          BLC(1,LQm5,NR) =   R(NR  ) / (DR * RHI(NR))
          NLC(1,LQm5,NR) = LQm2
      ENDDO
C
C Out of region
C
      NR = NRMAX
          BLC(1,LQm5,NR) = 1.D0
          NLC(1,LQm5,NR) = LQm5
C
      NLCMAX(LQm5) = 1
      RETURN
      END
C
C     ***************************************************************
C
C        Electron Density Equation (HI)
C
C     ***************************************************************
C
      SUBROUTINE LQe1CC
C
      INCLUDE 'txcomm.inc'
C
      DO NR = 0, NRMAX - 1
          BLC(0,LQe1,NR) = 1.D0 / DT
          NLC(0,LQe1,NR) = LQe1
C
C Convection
C
          ALC(1,LQe1,NR) = - R(NR+1) / (RHI(NR) * DR)
          BLC(1,LQe1,NR) = + R(NR  ) / (RHI(NR) * DR)
          NLC(1,LQe1,NR) = LQe2
C
C Ionization of n01 and n02
C
          BLC(2,LQe1,NR) = rNuION(NR) 
     &                    * PNeHI(NR) / (PN01HI(NR) + PN02HI(NR))
          NLC(2,LQe1,NR) = LQn1
C
          BLC(3,LQe1,NR) = rNuION(NR) 
     &                    * PNeHI(NR) / (PN01HI(NR) + PN02HI(NR))
          NLC(3,LQe1,NR) = LQn2
C
C Loss to divertor
C
          BLC(4,LQe1,NR) = - rNuL(NR)
          NLC(4,LQe1,NR) = LQe1
C
          PLC(5,LQe1,NR) =   rNuL(NR) * PNeDIV
C
      ENDDO
C
C     OUT OF REGION
C
      NR = NRMAX
          BLC(1,LQe1,NR) = 1.D0
          NLC(1,LQe1,NR) = LQe1
C
      NLCMAX(LQe1) = 5
      RETURN
      END
C
C     ***************************************************************
C
C        Electron Radial Flow (NR)
C
C     ***************************************************************
C
      SUBROUTINE LQe2CC
C
      INCLUDE 'txcomm.inc'
C
C Ns*Usr(0) : fixed
C
      NR = 0
          BLC(1,LQe2,NR) = 1.D0
          NLC(1,LQe2,NR) = LQe2
C
      DO NR = 1, NRMAX-1
          BLC(0,LQe2,NR) = 1.D0 / DT
          NLC(0,LQe2,NR) = LQe2
C
C Nonlinear term
C
          ALC(1,LQe2,NR) = - RHI(NR  ) * UerHI(NR  ) / (2 * R(NR) * DR)
          BLC(1,LQe2,NR) = - RHI(NR  ) * UerHI(NR  ) / (2 * R(NR) * DR)
     &                     + RHI(NR-1) * UerHI(NR-1) / (2 * R(NR) * DR)
          CLC(1,LQe2,NR) = + RHI(NR-1) * UerHI(NR-1) / (2 * R(NR) * DR)
          NLC(1,LQe2,NR) = LQe2
C
C Nonlinear centrifugal force
C
          BLC(2,LQe2,NR) = UethI(NR) / R(NR)
          NLC(2,LQe2,NR) = LQe3
C
C Pressure gradient force
C
          BLC(3,LQe2,NR) = - rKeV / (AME * DR)
          CLC(3,LQe2,NR) = + rKeV / (AME * DR)
          NLC(3,LQe2,NR) = LQe5
C
C Radial E force
C
          BLC(4,LQe2,NR) = - PNeI(NR) * (AEE/AME)
          NLC(4,LQe2,NR) = LQm1
C
C v x B force
C
          BLC(5,LQe2,NR) = - BphI(NR) * (AEE/AME)
          NLC(5,LQe2,NR) = LQe3
C
          BLC(6,LQe2,NR) =   BthI(NR) * (AEE/AME) / 2
          CLC(6,LQe2,NR) =   BthI(NR) * (AEE/AME) / 2
          NLC(6,LQe2,NR) = LQe4
C
      ENDDO
C
C Ns*Usr(NRMAX) : fixed or finite gradient
C
      NR = NRMAX
          BLC(1,LQe2,NR) = 1.D0
          NLC(1,LQe2,NR) = LQe2
C
      NLCMAX(LQe2) = 6
      RETURN
      END
C
C     ***************************************************************
C
C        Electron Poloidal Flow (NR)
C
C     ***************************************************************
C
      SUBROUTINE LQe3CC
C
      INCLUDE 'txcomm.inc'
C
C Ns*UsTheta(0) : 0
C
      NR = 0
          BLC(1,LQe3,NR) = 1.D0
          NLC(1,LQe3,NR) = LQe3
C
      DO NR = 1, NRMAX - 1
          BLC(0,LQe3,NR) = 1.D0 / DT
          NLC(0,LQe3,NR) = LQe3
C
C Nonlinear term
C
          ALC(1,LQe3,NR) = - RHI(NR  )**2 * UerHI(NR  )
     &                       / (2 * R(NR)**2 * DR)
          BLC(1,LQe3,NR) = - RHI(NR  )**2 * UerHI(NR  )
     &                       / (2 * R(NR)**2 * DR)
     &                     + RHI(NR-1)**2 * UerHI(NR-1)
     &                       / (2 * R(NR)**2 * DR)
          CLC(1,LQe3,NR) = + RHI(NR-1)**2 * UerHI(NR-1)
     &                       / (2 * R(NR)**2 * DR)
          NLC(1,LQe3,NR) = LQe3
C
C Viscosity force
C
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
C
C Poloidal E force
C
          BLC(3,LQe3,NR) = - PNeI(NR) * (AEE/AME)
          NLC(3,LQe3,NR) = LQm2
C
C v x B force
C
          BLC(4,LQe3,NR) = BphI(NR) * (AEE/AME)
          NLC(4,LQe3,NR) = LQe2
C
C Neoclassical viscosity force
C
          rNueNCL = 0.5D0 * (rNueNC(NR) + rNueNC(NR-1))
          BLC(5,LQe3,NR) = - rNueNCL
          NLC(5,LQe3,NR) = LQe3
C
C Collisional friction force with ions
C
          rNueiL = 0.5D0 * (rNuei(NR) + rNuei(NR-1))
          BLC(6,LQe3,NR) = - rNueiL
          NLC(6,LQe3,NR) = LQe3
C
          BLC(7,LQe3,NR) = + rNueiL * PNeI(NR) / PNiI(NR)
          NLC(7,LQe3,NR) = LQi3
C
C Collisional friction with beam ions
C
          rNubeL = 0.5D0 * (rNube(NR) + rNube(NR-1))
          BLC(8,LQe3,NR) = - (AMB/AME) * rNubeL * PNbI(NR) / PNeI(NR)
          NLC(8,LQe3,NR) = LQe3
C
          BLC(9,LQe3,NR) =   (AMB/AME) * rNubeL
          NLC(9,LQe3,NR) = LQb3
C
C Wave interaction force (electron driven)
C
         FWtheL = 0.5D0 * (FWthe(NR) + FWthe(NR-1))
           WPML = 0.5D0 * (  WPM(NR) +   WPM(NR-1))
          BLC(10,LQe3,NR) = - FWtheL / AME
          NLC(10,LQe3,NR) = LQe3
C
         TMP =   FWtheL * WPML * R(NR) / AME
          BLC(11,LQe3,NR) = TMP / 2
          CLC(11,LQe3,NR) = TMP / 2
          NLC(11,LQe3,NR) = LQe1
C
C Wave interaction force (NRon driven)
C
         FWthiL = 0.5D0 * (FWthi(NR) + FWthi(NR-1))
          BLC(12,LQe3,NR) =   FWthiL / AME
          NLC(12,LQe3,NR) = LQi3
C
         TMP = - FWthiL * WPML * R(NR) / AME
          BLC(13,LQe3,NR) = TMP / 2
          CLC(13,LQe3,NR) = TMP / 2
          NLC(13,LQe3,NR) = LQi1
C
C Loss to divertor
C
          rNuLL = 0.5D0 * (rNuL(NR) + rNuL(NR-1))
          BLC(14,LQe3,NR) = - 2.D0 * rNuLL
          NLC(14,LQe3,NR) = LQe3
C
C Collisional friction force with neutrals
C
          rNu0eL = 0.5D0 * (rNu0e(NR) + rNu0e(NR-1))
          BLC(15,LQe3,NR) = - rNu0eL
          NLC(15,LQe3,NR) = LQe3
C
      ENDDO
C
C Ns*UsTheta(NRMAX) : 0
C
      NR = NRMAX
          BLC(1,LQe3,NR) = 1
          NLC(1,LQe3,NR) = LQe3
C
      NLCMAX(LQe3) = 15
      RETURN
      END
C
C     ***************************************************************
C
C        Electron Toroidal Flow (NR)
C
C     ***************************************************************
C
      SUBROUTINE LQe4CC
C
      INCLUDE 'txcomm.inc'
C
C Uephi(0)' : 0
C
      DO NR = 0, NRMAX - 1
          BLC(0,LQe4,NR) = 1.D0 / DT
          NLC(0,LQe4,NR) = LQe4
C
C Nonlinear term
C
          ALC(1,LQe4,NR) = - R(NR+1) * UerI(NR+1) / (2 * RHI(NR) * DR)
          BLC(1,LQe4,NR) = - R(NR+1) * UerI(NR+1) / (2 * RHI(NR) * DR)
     &                    + R(NR  ) * UerI(NR  ) / (2 * RHI(NR) * DR)
          CLC(1,LQe4,NR) = + R(NR  ) * UerI(NR  ) / (2 * RHI(NR) * DR)
          NLC(1,LQe4,NR) = LQe4
C
C Viscosity force
C
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
C
C Toroidal E force
C
          BLC(3,LQe4,NR) = - PNeHI(NR  ) * (AEE/AME)
          NLC(3,LQe4,NR) = LQm3
C
C v x B force
C
          ALC(4,LQe4,NR) = - BthI(NR+1) * (AEE/AME) / 2
          BLC(4,LQe4,NR) = - BthI(NR  ) * (AEE/AME) / 2
          NLC(4,LQe4,NR) = LQe2
C
C Collisional friction with bulk ions
C
          BLC(5,LQe4,NR) = - rNuei(NR)
          NLC(5,LQe4,NR) = LQe4
C
          BLC(6,LQe4,NR) = + rNuei(NR) * PNeI(NR) / PNiI(NR)
          NLC(6,LQe4,NR) = LQi4
C
C Collisional friction with beam ions
C
          BLC(7,LQe4,NR) = - (AMB/AME) * rNube(NR) * PNbI(NR) / PNeI(NR)
          NLC(7,LQe4,NR) = LQe4
C
          BLC(8,LQe4,NR) =   (AMB/AME) * rNube(NR)
          NLC(8,LQe4,NR) = LQb4
C
C Loss to divertor
C
          BLC(9,LQe4,NR) = - 2.D0 * rNuL(NR)
          NLC(9,LQe4,NR) = LQe4
C
C Collisional friction force with neutrals
C
          BLC(10,LQe4,NR) = - rNu0e(NR)
          NLC(10,LQe4,NR) = LQe4
C
      ENDDO
C
C Out of region
C
      NR = NRMAX
          BLC(1,LQe4,NR) = 1.D0
          NLC(1,LQe4,NR) = LQe4
C
      NLCMAX(LQe4) = 10
      RETURN
      END
C
C     ***************************************************************
C
C        Electron Energy Transport: Te (HI)
C
C     ***************************************************************
C
      SUBROUTINE LQe5CC
C
      INCLUDE 'txcomm.inc'
C
C Fixed Temperature
C
      IF (ABS(Chie0) .EQ. 0.D0) THEN
         DO NR = 0, NRMAX
             BLC(0,LQe5,NR) = 1.D0
             NLC(0,LQe5,NR) = LQe5
         ENDDO
C
         NLCMAX(LQe5) = 0
      ELSE
C
C Temperature evolution
C
         DO NR = 0, NRMAX - 1
             BLC(0,LQe5,NR) = 1.5D0 / DT
             NLC(0,LQe5,NR) = LQe5
C
C Convection transport
C
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
C
C Conduction transport
C
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
C
C Joule heating
C
             ALC(3,LQe5,NR) = - AEE * EthI(NR+1) / (2 * rKeV)
             BLC(3,LQe5,NR) = - AEE * EthI(NR  ) / (2 * rKeV)
             NLC(3,LQe5,NR) = LQe3
C
             BLC(4,LQe5,NR) = - AEE * EphHI(NR) / rKeV
             NLC(4,LQe5,NR) = LQe4
C
C Collisional transfer with ions
C
             BLC(5,LQe5,NR) = - rNuTei(NR)
             NLC(5,LQe5,NR) = LQe5
C
             BLC(6,LQe5,NR) =   rNuTei(NR) * (PNeHI(NR)/PNiHI(NR))
             NLC(6,LQe5,NR) = LQi5
C
C Collisional heating with beam
C
             BLC(7,LQe5,NR) = - AMb * Vb * rNube(NR) / (2 * rKeV)
     &                       * (PNbHI(NR) / PNeHI(NR))
             NLC(7,LQe5,NR) = LQe4
C
             BLC(8,LQe5,NR) =   AMb * Vb * rNube(NR) / (2 * rKeV)
             NLC(8,LQe5,NR) = LQb4
C
C Loss to diverter
C
             BLC(9,LQe5,NR) = - 2.5D0 * rNuL(NR)
             NLC(9,LQe5,NR) = LQe5
C
             PLC(10,LQe5,NR) =  2.5D0 * rNuL(NR) * PNeHI(NR) * PTeDIV
C
C Direct heating (RF)
C
             PLC(11,LQe5,NR) =  PRFe(NR) / (1.D20 * rKeV)
         ENDDO
C
C Out of region
C
         NR = NRMAX
             BLC(1,LQe5,NR) = 1.D0
             NLC(1,LQe5,NR) = LQe5
C
         NLCMAX(LQe5) = 11
      ENDIF
      RETURN
      END
C
C     ***************************************************************
C
C        Ion Density Equation (HI)
C
C     ***************************************************************
C
      SUBROUTINE LQi1CC
C
      INCLUDE 'txcomm.inc'
C
      DO NR = 0, NRMAX - 1
          BLC(0,LQi1,NR) = 1.D0 / DT
          NLC(0,LQi1,NR) = LQi1
C
C Convection
C
          ALC(1,LQi1,NR) = - R(NR+1) / (RHI(NR) * DR)
          BLC(1,LQi1,NR) = + R(NR  ) / (RHI(NR) * DR)
          NLC(1,LQi1,NR) = LQi2
C
C Ionization of n01 and n02
C
          BLC(2,LQi1,NR) = rNuION(NR) / PZ
     &                   * PNeHI(NR) / (PN01HI(NR) + PN02HI(NR))
          NLC(2,LQi1,NR) = LQn1
C
          BLC(3,LQi1,NR) = rNuION(NR) / PZ
     &                   * PNeHI(NR) / (PN01HI(NR) + PN02HI(NR))
          NLC(3,LQi1,NR) = LQn2
C
C Loss to divertor
C
          BLC(4,LQi1,NR) = - rNuL(NR) / PZ
          NLC(4,LQi1,NR) = LQe1
C
          PLC(5,LQi1,NR) =   rNuL(NR) * PNeDIV / PZ
C
C Particle source from beam ion
C
          BLC(6,LQi1,NR) = rNuB(NR)
          NLC(6,LQi1,NR) = LQb1
C
C NBI kick up ions
C
         PLC(7,LQi1,NR) = - SNB(NR)
C
C Loss cone loss
C
         PLC(8,LQi1,NR) = SiLC(NR)
C
      ENDDO
C
C     OUT OF REGION
C
      NR = NRMAX
          BLC(1,LQi1,NR) = 1.D0
          NLC(1,LQi1,NR) = LQi1
C
      NLCMAX(LQi1) = 8
      RETURN
      END
C
C     ***************************************************************
C
C        Ion Radial Flow (NR)
C
C     ***************************************************************
C
      SUBROUTINE LQi2CC
C
      INCLUDE 'txcomm.inc'
C
C Ns*Usr(0) : fixed
C
      NR = 0
          BLC(1,LQi2,NR) = 1.D0
          NLC(1,LQi2,NR) = LQi2
C
      DO NR = 1, NRMAX-1
C
          BLC(0,LQi2,NR) = 1.D0 / DT
          NLC(0,LQi2,NR) = LQi2
C
C Nonlinear term
C
          ALC(1,LQi2,NR) = - RHI(NR  ) * UirHI(NR  ) / (2 * R(NR) * DR)
          BLC(1,LQi2,NR) = - RHI(NR  ) * UirHI(NR  ) / (2 * R(NR) * DR)
     &                    + RHI(NR-1) * UirHI(NR-1) / (2 * R(NR) * DR)
          CLC(1,LQi2,NR) = + RHI(NR-1) * UirHI(NR-1) / (2 * R(NR) * DR)
          NLC(1,LQi2,NR) = LQi2
C
C Nonlinear centrifugal force
C
          BLC(2,LQi2,NR) = UithI(NR) / R(NR)
          NLC(2,LQi2,NR) = LQi3
C
C Pressure gradient force
C
          BLC(3,LQi2,NR) = - rKEV / (AMI * DR)
          CLC(3,LQi2,NR) = + rKEV / (AMI * DR)
          NLC(3,LQi2,NR) = LQi5
C
C Radial E force
C
          BLC(4,LQi2,NR) =   PZ * PNiI(NR) * (AEE/AMI)
          NLC(4,LQi2,NR) = LQm1
C
C v x B force
C
          BLC(5,LQi2,NR) =   PZ * BphI(NR) * (AEE/AMI)
          NLC(5,LQi2,NR) = LQi3
C
          BLC(6,LQi2,NR) = - PZ * BthI(NR) * (AEE/AMI) / 2
          CLC(6,LQi2,NR) = - PZ * BthI(NR) * (AEE/AMI) / 2
          NLC(6,LQi2,NR) = LQi4
      ENDDO
C
C Ns*Usr(NRMAX) : fixed or finite gradient
C
      NR = NRMAX
          BLC(1,LQi2,NR) = 1.D0
          NLC(1,LQi2,NR) = LQi2
C
      NLCMAX(LQi2) = 6
      RETURN
      END
C
C     ***************************************************************
C
C        Ion Poloidal Flow (NR)
C
C     ***************************************************************
C
      SUBROUTINE LQi3CC
C
      INCLUDE 'txcomm.inc'
C
C Ni*UiTheta(0) : 0
C
      NR = 0
          BLC(1,LQi3,NR) = 1.D0
          NLC(1,LQi3,NR) = LQi3
C
      DO NR = 1, NRMAX - 1
C
          BLC(0,LQi3,NR) = 1.D0 / DT
          NLC(0,LQi3,NR) = LQi3
C
C Nonlinear term
C
          ALC(1,LQi3,NR) = - RHI(NR  )**2 * UirHI(NR  )
     &                    / (2 * R(NR)**2 * DR)
          BLC(1,LQi3,NR) = - RHI(NR  )**2 * UirHI(NR  )
     &                    / (2 * R(NR)**2 * DR)
     &                    + RHI(NR-1)**2 * UirHI(NR-1)
     &                    / (2 * R(NR)**2 * DR)
          CLC(1,LQi3,NR) = + RHI(NR-1)**2 * UirHI(NR-1)
     &                    / (2 * R(NR)**2 * DR)
          NLC(1,LQi3,NR) = LQi3
C
C Viscosity force
C
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
C
C Poroidal E force
C
          BLC(3,LQi3,NR) = PZ * PNiI(NR  ) * (AEE/AMI)
          NLC(3,LQi3,NR) = LQm2
C
C v x B force
C
          BLC(4,LQi3,NR) = - PZ * BphI(NR) * (AEE/AMI)
          NLC(4,LQi3,NR) = LQi2
C
C Neoclassical viscosity force
C
          rNuiNCL = 0.5D0 * (rNuiNC(NR) + rNuiNC(NR-1))
          BLC(5,LQi3,NR) = - rNuiNCL
          NLC(5,LQi3,NR) = LQi3
C
C Collisional friction force
C
          rNueiL = 0.5D0 * (rNuei(NR) + rNuei(NR-1))
          BLC(6,LQi3,NR) = - (AME/AMI) * rNueiL * PNeI(NR) / PNiI(NR)
          NLC(6,LQi3,NR) = LQi3
C
          BLC(7,LQi3,NR) =   (AME/AMI) * rNueiL
          NLC(7,LQi3,NR) = LQe3
C
C Collisional friction with beam ions
C
          rNubiL = 0.5D0 * (rNubi(NR) + rNubi(NR-1))
          BLC(8,LQi3,NR) = - (AMB/AMI) * rNubiL * PNbI(NR) / PNiI(NR)
          NLC(8,LQi3,NR) = LQi3
C
          BLC(9,LQi3,NR) =   (AMB/AMI) * rNubiL
          NLC(9,LQi3,NR) = LQb3
C
C Wave interaction force (electron driven)
C
         FWtheL = 0.5D0 * (FWthe(NR) + FWthe(NR-1))
           WPML = 0.5D0 * (  WPM(NR) +   WPM(NR-1))
          BLC(10,LQi3,NR) = FWtheL / AMI
          NLC(10,LQi3,NR) = LQe3
C
         TMP = - FWtheL * WPML * R(NR) / AMI
          BLC(11,LQi3,NR) = TMP / 2
          CLC(11,LQi3,NR) = TMP / 2
          NLC(11,LQi3,NR) = LQe1
C
C Wave interaction force (NRon driven)
C
         FWthiL = 0.5D0 * (FWthi(NR) + FWthi(NR-1))
          BLC(12,LQi3,NR) = - FWthiL / AMI
          NLC(12,LQi3,NR) = LQi3
C
         TMP = FWthiL * WPML * R(NR) / AMI
          BLC(13,LQi3,NR) = TMP / 2
          CLC(13,LQi3,NR) = TMP / 2
          NLC(13,LQi3,NR) = LQi1
C
C Loss to divertor
C
          rNuLL = 0.5D0 * (rNuL(NR) + rNuL(NR-1))
          BLC(14,LQi3,NR) = - 2.D0 * rNuLL
          NLC(14,LQi3,NR) = LQi3
C
C Collisional friction force with neutrals
C
          rNu0iL = 0.5D0 * (rNu0i(NR) + rNu0i(NR-1))
          BLC(15,LQi3,NR) = - rNu0iL
          NLC(15,LQi3,NR) = LQi3
C
C Charge exchange force
C
          rNuiCXL = 0.5D0 * (rNuiCX(NR) + rNuiCX(NR-1))
          BLC(16,LQi3,NR) = - rNuiCXL
          NLC(16,LQi3,NR) = LQi3
C
C Loss cone loss
C
          SiLCthL = 0.5D0 * (SiLCth(NR) + SiLCth(NR-1))
          PLC(17,LQi3,NR) = + SiLCthL
      ENDDO
C
C Ns*UsTheta(NRMAX) : 0
C
      NR = NRMAX
          BLC(1,LQi3,NR) = 1.D0
          NLC(1,LQi3,NR) = LQi3
C
      NLCMAX(LQi3) = 17
      RETURN
      END
C
C     ***************************************************************
C
C        Ion Toroidal Flow (NR)
C
C     ***************************************************************
C
      SUBROUTINE LQi4CC
C
      INCLUDE 'txcomm.inc'
C
C Uiphi'(0) : 0
C
      DO NR = 0, NRMAX - 1
          BLC(0,LQi4,NR) = 1.D0 / DT
          NLC(0,LQi4,NR) = LQi4
C
C Nonlinear term
C
          ALC(1,LQi4,NR) = - R(NR+1) * UirI(NR+1) / (2 * RHI(NR) * DR)
          BLC(1,LQi4,NR) = - R(NR+1) * UirI(NR+1) / (2 * RHI(NR) * DR)
     &                     + R(NR  ) * UirI(NR  ) / (2 * RHI(NR) * DR)
          CLC(1,LQi4,NR) = + R(NR  ) * UirI(NR  ) / (2 * RHI(NR) * DR)
          NLC(1,LQi4,NR) = LQi4
C
C Viscosity force
C
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
C
C Toroidal E force
C
          BLC(3,LQi4,NR) = PZ * PNiHI(NR  ) * (AEE/AMI)
          NLC(3,LQi4,NR) = LQm3
C
C v x B force
C
          ALC(4,LQi4,NR) = PZ * BthI(NR+1) * (AEE/AMI) / 2
          BLC(4,LQi4,NR) = PZ * BthI(NR  ) * (AEE/AMI) / 2 
          NLC(4,LQi4,NR) = LQi2
C
C Collisional friction with bulk ions
C
          BLC(5,LQi4,NR) = - (AME/AMI) * rNuei(NR) * PNeI(NR) / PNiI(NR)
          NLC(5,LQi4,NR) = LQi4
C
          BLC(6,LQi4,NR) =   (AME/AMI) * rNuei(NR)
          NLC(6,LQi4,NR) = LQe4
C
C Collisional friction with beam ions
C
          BLC(7,LQi4,NR) = - (AMB/AMI) * rNubi(NR) * PNbI(NR) / PNiI(NR)
          NLC(7,LQi4,NR) = LQi4
C
          BLC(8,LQi4,NR) =   (AMB/AMI) * rNubi(NR)
          NLC(8,LQi4,NR) = LQb4
C
C Loss to divertor
C
          BLC(9,LQi4,NR) = - 2.D0 * rNuL(NR)
          NLC(9,LQi4,NR) = LQi4
C
C Collisional friction force with neutrals
C
          BLC(10,LQi4,NR) = - rNu0i(NR)
          NLC(10,LQi4,NR) = LQi4
C
C Charge exchange force
C
          BLC(11,LQi4,NR) = - rNuiCX(NR)
          NLC(11,LQi4,NR) = LQi4
C
C Loss conde loss
C
          PLC(12,LQi4,NR) = SiLCph(NR)
C
      ENDDO
C
C Out of region
C
      NR = NRMAX
          BLC(1,LQi4,NR) = 1.D0
          NLC(1,LQi4,NR) = LQi4
C
      NLCMAX(LQi4) = 12
      RETURN
      END
C
C     ***************************************************************
C
C        Ion Energy Transport: Ti (HI)
C
C     ***************************************************************
C
      SUBROUTINE LQi5CC
C
      INCLUDE 'txcomm.inc'
C
C Fixed temperature
C
      IF (ABS(Chii0) .EQ. 0.D0) THEN
         DO NR = 0, NRMAX
             BLC(0,LQi5,NR) = 1.D0
             NLC(0,LQi5,NR) = LQi5
         ENDDO
         NLCMAX(LQi5) = 0
C
C Temperature evolution
C
      ELSE
         DO NR = 0, NRMAX-1
             BLC(0,LQi5,NR) = 1.5D0 / DT
             NLC(0,LQi5,NR) = LQi5
C
C Convection transport
C
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
C
C Conduction transport
C
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
C
C Joule heating
C
             ALC(3,LQi5,NR) =   PZ * AEE * EthI(NR+1) / (2 * rKeV)
             BLC(3,LQi5,NR) =   PZ * AEE * EthI(NR  ) / (2 * rKeV)
             NLC(3,LQi5,NR) = LQi3
C
             BLC(4,LQi5,NR) =   PZ * AEE * EphHI(NR) / rKeV
             NLC(4,LQi5,NR) = LQi4
C
C Collisional transfer with electrons
C
             BLC(5,LQi5,NR) = - rNuTei(NR) * PNeHI(NR) / PNiHI(NR)
             NLC(5,LQi5,NR) = LQi5
C
             BLC(6,LQi5,NR) =   rNuTei(NR)
             NLC(6,LQi5,NR) = LQe5
C
C Collisional heating with beam
C
             BLC(7,LQi5,NR) = - AMb * Vb * rNubi(NR) / (2 * rKeV)
     &                       * (PNbHI(NR) / PNiHI(NR))
             NLC(7,LQi5,NR) = LQi4
C
             BLC(8,LQi5,NR) =   AMb * Vb * rNubi(NR) / (2 * rKeV)
             NLC(8,LQi5,NR) = LQb4
C
C Loss to diverter
C
             BLC(9,LQi5,NR) = - 2.5D0 * rNuL(NR)
             NLC(9,LQi5,NR) = LQi5
C
             PLC(10,LQi5,NR) =  2.5D0 * rNuL(NR) * PNiHI(NR) * PTiDIV
C
C Direct heating (RF)
C
             PLC(11,LQi5,NR) =  PRFi(NR) / (1.D20 * rKeV)
         ENDDO
C
C Out of region
C
         NR = NRMAX
             BLC(1,LQi5,NR) = 1.D0
             NLC(1,LQi5,NR) = LQi5
C
         NLCMAX(LQi5) = 11
      ENDIF
      RETURN
      END
C
C     ***************************************************************
C
C        Beam Ion Density (HI)
C
C     ***************************************************************
C
      SUBROUTINE LQb1CC
C
      INCLUDE 'txcomm.inc'
C
      DO NR = 0, NRMAX - 1
          BLC(0,LQb1,NR) = 1.D0 / DT
          NLC(0,LQb1,NR) = LQb1
C
C NBI particle source
C
          PLC(1,LQb1,NR) = SNB(NR)
C
C Relaxation to thermal ions
C
          BLC(2,LQb1,NR) = - rNuB(NR)
          NLC(2,LQb1,NR) = LQb1
      ENDDO
C
C Out of region
C
      NR = NRMAX
          BLC(1,LQb1,NR) = 1.D0
          NLC(1,LQb1,NR) = LQb1
C
      NLCMAX(LQb1) = 2
      RETURN
      END
C
C     ***************************************************************
C
C        Beam Ion Poloidal Flow (NR)
C
C     ***************************************************************
C
      SUBROUTINE LQb3CC
C
      INCLUDE 'txcomm.inc'
C
C Ubth(0) : 0
C
      NR = 0
          BLC(1,LQb3,NR) = 1.D0
          NLC(1,LQb3,NR) = LQb3
C
      DO NR = 1, NRMAX-1
C
C Nonlinear centrifugal force
C
          BLC(1,LQb3,NR) = UbthI(NR) / R(NR)
          NLC(1,LQb3,NR) = LQb3
C
C Radial E force
C
          BLC(2,LQb3,NR) = PZ * PNbI(NR) * (AEE/AMB)
          NLC(2,LQb3,NR) = LQm1
C
C v x B force
C
          BLC(3,LQb3,NR) = PZ * BphI(NR) * (AEE/AMB)
          NLC(3,LQb3,NR) = LQb3
C
          BLC(4,LQb3,NR) = - PZ * BthI(NR) * (AEE/AMB)
          NLC(4,LQb3,NR) = LQb4
      ENDDO
C
C Ubth(NRMAX) : 0
C
      NR = NRMAX
          BLC(1,LQb3,NR) = 1.D0
          NLC(1,LQb3,NR) = LQb3
C
      NLCMAX(LQb3) = 4
      RETURN
      END
C
C     ***************************************************************
C
C        Beam Ion Toroildal Flow (NR)
C
C     ***************************************************************
C
      SUBROUTINE LQb4CC
C
      INCLUDE 'txcomm.inc'
C
C Ubphi'(0) : 0
C
      NR = 0
          IF (PNbI(0)*PNbI(1) .EQ. 0.D0) THEN
             ALC(1,LQb4,NR) = -1.D0
             BLC(1,LQb4,NR) =  1.D0
          ELSE
             ALC(1,LQb4,NR) = -1.D0 / PNbI(NR+1)
             BLC(1,LQb4,NR) =  1.D0 / PNbI(NR  )
          ENDIF
          NLC(1,LQb4,NR) = LQb4
C
      DO NR = 1, NRMAX - 1
          BLC(0,LQb4,NR) = 1.D0 / DT
          NLC(0,LQb4,NR) = LQb4
C
C Collisional friction with electrons
C
          BLC(1,LQb4,NR) = - rNube(NR)
          NLC(1,LQb4,NR) = LQb4
C
          BLC(2,LQb4,NR) =   rNube(NR) * PNbI(NR) / PNeI(NR)
          NLC(2,LQb4,NR) = LQe4
C
C Collisional friction with ions
C
          BLC(3,LQb4,NR) = - rNubi(NR)
          NLC(3,LQb4,NR) = LQb4
C
          BLC(4,LQb4,NR) =   rNubi(NR)
          NLC(4,LQb4,NR) = LQi4
C
C NBI momentum source
C
         SNBL = 0.5D0 * (SNB(NR) + SNB(NR-1))
         PLC(5,LQb4,NR) = PNBCD * Vb * SNBL
      ENDDO
C
C Ubphi(NRMAX) : 0
C
      NR = NRMAX
          BLC(1,LQb4,NR) = 1.D0
          NLC(1,LQb4,NR) = LQb4
C
      NLCMAX(LQb4) = 5
      RETURN
      END
C
C     ***************************************************************
C
C        Slow Neutral Transport: n01 (HI)
C
C     ***************************************************************
C
      SUBROUTINE LQn1CC
C
      INCLUDE 'txcomm.inc'
C
      DO NR = 0, NRMAX-1
          BLC(0,LQn1,NR) = 1.D0 / DT
          NLC(0,LQn1,NR) = LQn1
C
C  Diffusion of neutrals
C
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
C
C Ionization
C
          BLC(2,LQn1,NR) = - rNuION(NR) / PZ
     &                    * PNeHI(NR) / (PN01HI(NR) + PN02HI(NR))
          NLC(2,LQn1,NR) = LQn1
C
C Generation of fast neutrals by charge exchange
C
         IF(NR.NE.0) THEN
             rNuiCXL = ( rNuiCX(NR  )/(PN01HI(NR  ) + PN02HI(NR  ))
     &                 + rNuiCX(NR-1)/(PN01HI(NR-1) + PN02HI(NR-1)))/2
         ELSE
             rNuiCXL =   rNuiCX(NR  )/(PN01HI(NR  ) + PN02HI(NR  ))
         ENDIF
          BLC(3,LQn1,NR) = - rNuiCXL*PNiHI(NR)
          NLC(3,LQn1,NR) = LQn1
C
C Recycling from divertor
C
          BLC(4,LQn1,NR) =   rGamm0 * rNuL(NR) / PZ
          NLC(4,LQn1,NR) = LQe1
C
          PLC(5,LQn1,NR) = - rGamm0 * rNuL(NR) * PneDIV / PZ
C
      ENDDO
C
C Out of region
C
      NR = NRMAX
          BLC(1,LQn1,NR) = 1.D0
          NLC(1,LQn1,NR) = LQn1
C
      NLCMAX(LQn1) = 5
C
      RETURN
      END
C
C     ***************************************************************
C
C        Fast Neutral Transport: n02 (HI)
C
C     ***************************************************************
C
      SUBROUTINE LQn2CC
C
      INCLUDE 'txcomm.inc'
C
      DO NR = 0, NRMAX-1
          BLC(0,LQn2,NR) = 1.D0 / DT
          NLC(0,LQn2,NR) = LQn2
C
C  Diffusion of neutrals
C
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
C
C Ionization
C
          BLC(2,LQn2,NR) = - rNuION(NR) / PZ
     &                    * PNeHI(NR) / (PN01HI(NR) + PN02HI(NR))
          NLC(2,LQn2,NR) = LQn2
C
C Generation of fast neutrals by charge exchange
C
         IF(NR.NE.0) THEN
             rNuiCXL = ( rNuiCX(NR  )/(PN01HI(NR  ) + PN02HI(NR  ))
     &                 + rNuiCX(NR-1)/(PN01HI(NR-1) + PN02HI(NR-1)))/2
         ELSE
             rNuiCXL =   rNuiCX(NR  )/(PN01HI(NR  ) + PN02HI(NR  ))
         ENDIF
          BLC(3,LQn2,NR) = rNuiCXL * PNiHI(NR)
          NLC(3,LQn2,NR) = LQn1
C
C NBI particle source
C
          PLC(4,LQn2,NR) = SNB(NR)
C
      ENDDO
C
C Out of region
C
      NR = NRMAX
          BLC(1,LQn2,NR) = 1.D0
          NLC(1,LQn2,NR) = LQn2
C
      NLCMAX(LQn2) = 4
C
      RETURN
      END
