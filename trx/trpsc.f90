!     ***********************************************************

!           Continuous particle source

!     ***********************************************************

      SUBROUTINE TRPSC

      USE TRCOMM, ONLY : MDLPSC,NPSCMAX
      IMPLICIT NONE
      INTEGER:: NPSC

      DO NPSC=1,NPSCMAX
         IF(MDLPSC(NPSC).NE.0) GO TO 100
      END DO
      RETURN

100   CONTINUE
      CALL TRPSCA

      RETURN
      END  SUBROUTINE TRPSC

!     ***********************************************************

!           Particle source (GAUSSIAN PROFILE)

!     ***********************************************************

      SUBROUTINE TRPSCA

      USE TRCOMM, ONLY : DR, RA, RM, DVRHO, NRMAX, PZ, &
           NPSCMAX, PSCIN, PSCR0, PSCRW, NSPSC, SPSC_NSNR, rkind, NSMAX
      IMPLICIT NONE
      REAL(rkind)    :: SSUM, S0, SPSCL
      INTEGER :: NR, NS, NPSC
      
      SPSC_NSNR(1:NSMAX,1:NRMAX)=0.D0
      
      DO NPSC=1,NPSCMAX
         SSUM = 0.D0
         DO NR=1,NRMAX
            SSUM=SSUM+DEXP(-((RA*RM(NR)-PSCR0(NPSC))/PSCRW(NPSC))**2) &
                     *DVRHO(NR)*DR
         ENDDO

         S0 = PSCIN(NPSC)/SSUM
         NS = NSPSC(NPSC)

         DO NR=1,NRMAX
            SPSCL = S0*DEXP(-((RA*RM(NR)-PSCR0(NPSC))/PSCRW(NPSC))**2)
            SPSC_NSNR( 1,NR)=SPSC_NSNR( 1,NR) + PZ(NS)*SPSCL  ! Electron
            SPSC_NSNR(NS,NR)=SPSC_NSNR(NS,NR) +        SPSCL  ! Ion
         ENDDO
      ENDDO

      RETURN
    END  SUBROUTINE TRPSCA

