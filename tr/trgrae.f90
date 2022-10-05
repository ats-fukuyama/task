!     ***********************************************************

!           GRAPHIC : CONTROL ROUTINE

!     ***********************************************************

      SUBROUTINE TRGRE0(K2,INQ)

      USE TRCOMM, ONLY : NRAMAX, NRMAX, NROMAX, RHOA, MODELG
      IMPLICIT NONE
      CHARACTER(LEN=1),INTENT(IN):: K2
      INTEGER,      INTENT(IN):: INQ

      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      IF(K2.EQ.'1') CALL TRGRE1(INQ)
      IF(K2.EQ.'2') CALL TRGRE2(INQ)
      IF(K2.EQ.'3') CALL TRGRE3(INQ)
      IF(K2.EQ.'4') CALL TRGRE4(INQ)
      IF(K2.EQ.'5') CALL TRGRE5(INQ)
      IF(K2.EQ.'6') CALL TRGRE6(INQ)
      IF(K2.EQ.'Q'.AND.MODELG.GE.3) CALL EQGOUT(0)

      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
      RETURN
      END SUBROUTINE TRGRE0

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : VGR

!     ***********************************************************

      SUBROUTINE TRGRE1(INQ)

      USE TRCOMM,ONLY : NRMAX, NRMP, GRM, GYR, DVRHO, ABRHO, ARRHO, AR1RHO
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: INQ
      INTEGER:: NR
      REAL   :: GUCLIP

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(DVRHO(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,1, &
     &            '@DVRHO vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(ABRHO(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,1, &
     &            '@ABRHO vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(ARRHO(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1, &
     &            '@ARRHO vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(AR1RHO(NR))
      ENDDO

      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1, &
     &            '@AR1RHO vs RHO@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRE1

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : VGR

!     ***********************************************************

      SUBROUTINE TRGRE2(INQ)

      USE TRCOMM,ONLY : NRMAX, NRMP, GRM, GYR, AR2RHO,RKPRHO
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: INQ
      INTEGER:: NR
      REAL   :: GUCLIP

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(AR2RHO(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,1, &
     &            '@AR2RHO vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(RKPRHO(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,1, &
     &            '@RKPRHO vs RHO@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRE2

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : VGR

!     ***********************************************************

      SUBROUTINE TRGRE3(INQ)

      USE TRCOMM,ONLY : NRMAX, NRMP, GRG, GYR, DVRHOG, ABRHOG, ARRHOG, AR1RHOG
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: INQ
      INTEGER:: NR
      REAL   :: GUCLIP

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(DVRHOG(NR))
      ENDDO
      GYR(1,1) = 0.0
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1, &
     &            '@DVRHOG vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(ABRHOG(NR))
      ENDDO
      GYR(1,1) = (4*GYR(2,1)-GYR(3,1))/3.0
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1, &
     &            '@ABRHOG vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(ARRHOG(NR))
      ENDDO
      GYR(1,1) = (4*GYR(2,1)-GYR(3,1))/3.0
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,1, &
     &            '@ARRHOG vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(AR1RHOG(NR))
      ENDDO
      GYR(1,1) = (4*GYR(2,1)-GYR(3,1))/3.0
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,1, &
     &            '@AR1RHOG vs RHO@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRE3

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : VGR

!     ***********************************************************

      SUBROUTINE TRGRE4(INQ)

      USE TRCOMM,ONLY : NRMAX, NRMP,GRG, GYR, AR2RHOG, ABB2RHOG, AIB2RHOG, &
           & ARHBRHOG
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: INQ
      INTEGER:: NR
      REAL   :: GUCLIP

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(AR2RHOG(NR))
      ENDDO
      GYR(1,1) = (4*GYR(2,1)-GYR(3,1))/3.0
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1, &
     &            '@AR2RHOG vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(ABB2RHOG(NR))
      ENDDO
      GYR(1,1) = (4*GYR(2,1)-GYR(3,1))/3.0
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1, &
     &            '@ABB2RHOG vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(AIB2RHOG(NR))
      ENDDO
      GYR(1,1) = (4*GYR(2,1)-GYR(3,1))/3.0
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,1, &
     &            '@AIB2RHOG vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(ARHBRHOG(NR))
      ENDDO
      GYR(1,1) = (4*GYR(2,1)-GYR(3,1))/3.0
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,1, &
     &            '@ARHBRHOG vs RHO@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRE4

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : VGR

!     ***********************************************************

      SUBROUTINE TRGRE5(INQ)

      USE TRCOMM,ONLY : NRMAX, NRMP, GRG, GYR, EPSRHO, RMJRHO, RMNRHO, RKPRHOG
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: INQ
      INTEGER:: NR
      REAL   :: GUCLIP

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(EPSRHO(NR))
      ENDDO
      GYR(1,1) = 0.0
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1, &
     &            '@EPSRHO vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(RMJRHO(NR))
      ENDDO
      GYR(1,1) = (4*GYR(2,1)-GYR(3,1))/3.0
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1, &
     &            '@RMJRHO  vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(RMNRHO(NR))
      ENDDO
      GYR(1,1) = 0.D0
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,1, &
     &            '@RMNRHO vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(RKPRHOG(NR))
      ENDDO
      GYR(1,1) = (4*GYR(2,1)-GYR(3,1))/3.0
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,1, &
     &            '@RKPRHOG vs RHO@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRE5

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : VGR

!     ***********************************************************

      SUBROUTINE TRGRE6(INQ)

      USE TRCOMM,ONLY : NRMAX, NRMP, GRM, GYR, PPPRHO, PIQRHO, TTRHO, PIRHO
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: INQ
      INTEGER:: NR
      REAL   :: GUCLIP

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(PPPRHO(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,1, &
     &            '@PPPRHO vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(PIQRHO(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,1, &
     &            '@PIQRHO  vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(TTRHO(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1, &
     &            '@TTRHO vs RHO@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(PIRHO(NR))
      ENDDO

      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1, &
     &            '@PIRHO vs RHO@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRE6
