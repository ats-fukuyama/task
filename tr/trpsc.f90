!     ***********************************************************

!           Continuous particle source

!     ***********************************************************

      SUBROUTINE TRPSC

      USE TRCOMM, ONLY : MDLPSC
      IMPLICIT NONE

      IF(MDLPSC.EQ.0) RETURN

      IF(MDLPSC.EQ.1) THEN
         CALL TRPSCA
      ENDIF

      RETURN
      END  SUBROUTINE TRPSC

!     ***********************************************************

!           Particle source (GAUSSIAN PROFILE)

!     ***********************************************************

      SUBROUTINE TRPSCA

      USE TRCOMM, ONLY : DR, RA, RM, DVRHO, NRMAX, PZ, &
           MDLPSC, NPSCMAX, PSCTOT, PSCR0, PSCRW, NSPSC, SPSC
      IMPLICIT NONE
      REAL(8)    :: SSUM, S0, SPSCL
      INTEGER(4) :: NR, NS, NPSC

      DO NPSC=1,NPSCMAX
         SSUM = 0.D0
         DO NR=1,NRMAX
            SSUM=SSUM+DEXP(-((RA*RM(NR)-PSCR0(NPSC))/PSCRW(NPSC))**2) &
                     *DVRHO(NR)*DR
         ENDDO

         S0 = PSCTOT(NPSC)/SSUM
         NS = NSPSC(NPSC)

         DO NR=1,NRMAX
            SPSCL = S0*DEXP(-((RA*RM(NR)-PSCR0(NPSC))/PSCRW(NPSC))**2)
            SPSC(NR, 1)=SPSC(NR, 1) + PZ(NS)*SPSCL  ! Electron
            SPSC(NR,NS)=SPSC(NR,NS) +        SPSCL  ! Ion
         ENDDO
      ENDDO

      RETURN
    END  SUBROUTINE TRPSCA

