!     ***********************************************************

!           GRAPHIC : CONTROL ROUTINE

!     ***********************************************************

      SUBROUTINE TRGRG0(K2,INQ)

      USE TRCOMM, ONLY : NRAMAX, NRMAX, NROMAX, RHOA
      IMPLICIT NONE
      CHARACTER(LEN=1),INTENT(IN):: K2
      INTEGER,      INTENT(IN):: INQ


      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      IF(K2.EQ.'1') CALL TRGRG1(INQ)
      IF(K2.EQ.'2') CALL TRGRG2(INQ)
      IF(K2.EQ.'3') CALL TRGRG3(INQ)
      IF(K2.EQ.'4') CALL TRGRG4(INQ)
      IF(K2.EQ.'5') CALL TRGRG5(INQ)
      IF(K2.EQ.'6') CALL TRGRG6(INQ)
      IF(K2.EQ.'7') CALL TRGRG7(INQ)
      IF(K2.EQ.'A') CALL TRGRGA(INQ)
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
      RETURN
      END SUBROUTINE TRGRG0
!
!     ***********************************************************

!           GRAPHIC : CONTROL ROUTINE

!     ***********************************************************

      SUBROUTINE TRGRP0(K2,INQ)

      IMPLICIT NONE
      CHARACTER(LEN=1),INTENT(IN):: K2
      INTEGER,      INTENT(IN):: INQ


      IF(K2.EQ.'1') CALL TRGRP1(INQ)
      IF(K2.EQ.'2') CALL TRGRP2(INQ)
      IF(K2.EQ.'3') CALL TRGRP3(INQ)
      IF(K2.EQ.'4') CALL TRGRP4(INQ)
      IF(K2.EQ.'5') CALL TRGRP5(INQ)
      RETURN
      END SUBROUTINE TRGRP0

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : NE,ND,NT,NA

!     ***********************************************************

      SUBROUTINE TRGRG1(INQ)

      USE TRCOMM, ONLY : GRM, GVR, NGR, NRMAX, NRMP
      IMPLICIT NONE
      INTEGER,INTENT(IN):: INQ


      CALL PAGES

      CALL TRGR1D( 3.0,12.0,11.0,17.0, &
     &            GRM,GVR(1:NRMAX,1,1),NRMP,NRMAX,NGR, &
     &            '@NE [10^20/m^3]  vs r@',0+INQ)

      CALL TRGR1D(15.5,24.5,11.0,17.0, &
     &            GRM,GVR(1:NRMAX,1,2),NRMP,NRMAX,NGR, &
     &            '@ND [10^20/m^3]  vs r@',0+INQ)

      CALL TRGR1D( 3.0,12.0, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX,1,5),NRMP,NRMAX,NGR, &
     &            '@TE [keV]  vs r@',0+INQ)

      CALL TRGR1D(15.5,24.5, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX,1,6),NRMP,NRMAX,NGR, &
     &            '@TD [keV]  vs r@',0+INQ)

      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRG1

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : TE,TD,TT,TA

!     ***********************************************************

      SUBROUTINE TRGRG2(INQ)

      USE TRCOMM, ONLY : GRM, GVR, NGR, NRMAX, NRMP
      IMPLICIT NONE
      INTEGER,INTENT(IN):: INQ


      CALL PAGES

      CALL TRGR1D( 3.0,12.0,11.0,17.0, &
     &            GRM,GVR(1:NRMAX,1,3),NRMP,NRMAX,NGR, &
     &            '@NT [10^20/m^3]  vs r@',0+INQ)

      CALL TRGR1D(15.5,24.5,11.0,17.0, &
     &            GRM,GVR(1:NRMAX,1,4),NRMP,NRMAX,NGR, &
     &            '@NA [10^20/m^3]  vs r@',0+INQ)

      CALL TRGR1D( 3.0,12.0, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX,1,7),NRMP,NRMAX,NGR, &
     &            '@TT [keV]  vs r@',0+INQ)

      CALL TRGR1D(15.5,24.5, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX,1,8),NRMP,NRMAX,NGR, &
     &            '@TA [keV]  vs r@',0+INQ)

      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRG2

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : Q,J,EZOH,JOH

!     ***********************************************************

      SUBROUTINE TRGRG3(INQ)

      USE TRCOMM, ONLY : GRG, GRM, GVR, NGR, NRMAX, NRMP
      IMPLICIT NONE
      INTEGER,INTENT(IN):: INQ


      CALL PAGES

      CALL TRGR1D( 3.0,12.0,11.0,17.0, &
     &            GRG,GVR(1:NRMAX+1,1,9),NRMP,NRMAX+1,NGR, &
     &           '@QP  vs r@',0+INQ)

      CALL TRGR1D(15.5,24.5,11.0,17.0, &
     &            GRM,GVR(1:NRMAX,1,10),NRMP,NRMAX,NGR, &
     &           '@AJ [MA/m^2]  vs r@',0+INQ)

      CALL TRGR1D( 3.0,12.0, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX,1,11),NRMP,NRMAX,NGR, &
     &           '@EZ [V/m]  vs r@',0+INQ)

      CALL TRGR1D(15.5,24.5, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX,1,12),NRMP,NRMAX,NGR, &
     &           '@AJOH [MA/m^2]  vs r@',0+INQ)

      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRG3

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : JNB,JBS,WB,WF

!     ***********************************************************

      SUBROUTINE TRGRG4(INQ)

      USE TRCOMM, ONLY : GRM, GVR, NGR, NRMAX, NRMP
      IMPLICIT NONE
      INTEGER,INTENT(IN):: INQ


      CALL PAGES

      CALL TRGR1D( 3.0,12.0,11.0,17.0, &
     &            GRM,GVR(1:NRMAX,1,13),NRMP,NRMAX,NGR, &
     &            '@AJNB+AJRF [MA/m^2]  vs r@',0+INQ)

      CALL TRGR1D(15.5,24.5,11.0,17.0, &
     &            GRM,GVR(1:NRMAX,1,14),NRMP,NRMAX,NGR, &
     &            '@AJBS [MA/m^2]  vs r@',0+INQ)

      CALL TRGR1D( 3.0,12.0, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX,1,15),NRMP,NRMAX,NGR, &
     &            '@PIN [MW/m^3]  vs r@',0+INQ)

      CALL TRGR1D(15.5,24.5, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX,1,16),NRMP,NRMAX,NGR, &
     &            '@POH [MW/m^3]  vs r@',0+INQ)

      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRG4

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : G,s,alpha,chiD

!     ***********************************************************

      SUBROUTINE TRGRG5(INQ)

      USE TRCOMM, ONLY : GRM, GVR, NGR, NRMAX, NRMP
      IMPLICIT NONE
      INTEGER,INTENT(IN):: INQ


      CALL PAGES

      CALL TRGR1D( 3.0,12.0,11.0,17.0, &
     &            GRM,GVR(1:NRMAX,1,17),NRMP,NRMAX,NGR, &
     &            '@s  vs r@',0+INQ)

      CALL TRGR1D(15.5,24.5,11.0,17.0, &
     &            GRM,GVR(1:NRMAX,1,18),NRMP,NRMAX,NGR, &
     &            '@G  vs r@',0+INQ)

      CALL TRGR1D( 3.0,12.0, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX,1,19),NRMP,NRMAX,NGR, &
     &            '@alpha  vs r@',0+INQ)

      CALL TRGR1D(15.5,24.5, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX-1,1,20),NRMP,NRMAX-1,NGR, &
     &            '@AKD [m^2/s]  vs r@',0+INQ)

      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRG5

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : s-alpha,s,alpha,chiD

!     ***********************************************************

      SUBROUTINE TRGRG6(INQ)

      USE TRCOMM, ONLY : GRM, GVR, NGR, NRMAX, NRMP
      IMPLICIT NONE
      INTEGER,INTENT(IN):: INQ


      CALL PAGES

      CALL TRGR1DX( 3.0,12.0,11.0,17.0, &
     &            GVR(1:NRMAX,1,19),GVR(1:NRMAX,1,18),NRMP,NRMAX,NGR, &
     &            '@s vs alpha@',0+INQ)

      CALL TRGR1D(15.5,24.5,11.0,17.0, &
     &            GRM,GVR(1:NRMAX,1,18),NRMP,NRMAX,NGR, &
     &            '@s  vs r@',0+INQ)

      CALL TRGR1D( 3.0,12.0, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX,1,19),NRMP,NRMAX,NGR, &
     &            '@alpha  vs r@',0+INQ)

      CALL TRGR1D(15.5,24.5, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX-1,1,20),NRMP,NRMAX-1,NGR, &
     &            '@AKD [m^2/s]  vs r@',0+INQ)

      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRG6


!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : BP,PSI

!     ***********************************************************

      SUBROUTINE TRGRG7(INQ)

      USE TRCOMM, ONLY : GRG, GRM, GVR, NGR, NRMAX, NRMP
      IMPLICIT NONE
      INTEGER,INTENT(IN):: INQ


      CALL PAGES

      CALL TRGR1D( 3.0,12.0,11.0,17.0, &
     &            GRG,GVR(1:NRMAX,1,21),NRMP,NRMAX,NGR, &
     &            '@BP vs r@',0+INQ)

      CALL TRGR1D(15.5,24.5,11.0,17.0, &
     &            GRM,GVR(1:NRMAX,1,22),NRMP,NRMAX,NGR, &
     &            '@PSI vs r@',0+INQ)

!$$$      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,
!$$$     &            GRM,GVR(1,1,23),NRMP,NRMAX,NGR,
!$$$     &            '@alpha  vs r@',0+INQ)
!$$$C
!$$$      CALL TRGR1D(15.5,24.5, 2.0, 8.0,
!$$$     &            GRM,GVR(1,1,24),NRMP,NRMAX-1,NGR,
!$$$     &            '@AKD [m^2/s]  vs r@',0+INQ)

      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRG7

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : Pressure,J(r),s(r),chi(r)

!     ***********************************************************

      SUBROUTINE TRGRGA(INQ)

      USE TRCOMM, ONLY : GRM, GVR, NGM, NGR, NGR, NRMAX, NRMP, RKEV
      IMPLICIT NONE
      INTEGER,      INTENT(IN):: INQ
      INTEGER                 :: NG, NR
      REAL,DIMENSION(NRMP,NGM):: GWR
      REAL :: GUCLIP


      CALL PAGES

      DO NG=1,NGR
      DO NR=1,NRMAX
         GWR(NR,NG) =  GVR(NR,NG,1)*GVR(NR,NG,5)*GUCLIP(RKEV*1.D14) &
     &              +  GVR(NR,NG,2)*GVR(NR,NG,6)*GUCLIP(RKEV*1.D14)
      ENDDO
      ENDDO

      CALL TRGR1D( 3.0,12.0,11.0,17.0, &
     &            GRM,GWR,NRMP,NRMAX,NGR, &
     &            '@P [MP]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GWR(NR,1) =  GVR(NR,1,10)
         GWR(NR,2) =  GVR(NR,1,12)
         GWR(NR,3) =  0.0
         GWR(NR,4) =  GVR(NR,1,14)
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0, &
     &            GRM,GWR,NRMP,NRMAX,4, &
     &            '@J,JOH,JBS [MA/m^2]  vs r@',2+INQ)

      CALL TRGR1D( 3.0,12.0, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX,1,18),NRMP,NRMAX,NGR, &
     &            '@s  vs r@',2+INQ)

      CALL TRGR1D(15.5,24.5, 2.0, 8.0, &
     &            GRM,GVR(1:NRMAX-1,1,20),NRMP,NRMAX-1,NGR, &
     &            '@AKD [m^2/s]  vs r@',2+INQ)

      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRGA

!     ***********************************************************

!           THREE DIMENSIONAL : DENSITY PROFILE

!     ***********************************************************

      SUBROUTINE TRGRP1(INQ)

      USE TRCOMM, ONLY : GVR, NGR, NRMAX, NRMP
      IMPLICIT NONE
      INTEGER,      INTENT(IN):: INQ


      CALL PAGES

      CALL TRGR2D( 1.0,12.0,10.0,17.0,GVR(1,1,1),NRMP,NRMAX,NGR, &
     &            '@NE [10^20/m^3]@',INQ)

      CALL TRGR2D(13.5,24.5,10.0,17.0,GVR(1,1,2),NRMP,NRMAX,NGR, &
     &            '@ND [10^20/m^3]@',INQ)

      CALL TRGR2D( 1.0,12.0, 1.0, 8.0,GVR(1,1,3),NRMP,NRMAX,NGR, &
     &            '@NT [10^20/m^3]@',INQ)

      CALL TRGR2D(13.5,24.5, 1.0, 8.0,GVR(1,1,4),NRMP,NRMAX,NGR, &
     &            '@NA [10^20/m^3]@',INQ)

      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRP1

!     ***********************************************************

!           THREE DIMENSIONAL : TEMPERATURE PROFILE

!     ***********************************************************

      SUBROUTINE TRGRP2(INQ)

      USE TRCOMM, ONLY : GVR, NGR, NRMAX, NRMP
      IMPLICIT NONE
      INTEGER,      INTENT(IN):: INQ


      CALL PAGES

      CALL TRGR2D( 1.0,12.0,10.0,17.0,GVR(1,1,5),NRMP,NRMAX,NGR, &
     &            '@TE [keV]@',INQ)

      CALL TRGR2D(13.5,24.5,10.0,17.0,GVR(1,1,6),NRMP,NRMAX,NGR, &
     &            '@TD [keV]@',INQ)

      CALL TRGR2D( 1.0,12.0, 1.0, 8.0,GVR(1,1,7),NRMP,NRMAX,NGR, &
     &            '@TT [keV]@',INQ)

      CALL TRGR2D(13.5,24.5, 1.0, 8.0,GVR(1,1,8),NRMP,NRMAX,NGR, &
     &            '@TA [keV]@',INQ)

      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRP2

!     ***********************************************************

!           THREE DIMENSIONAL : QP,AJ,EZOH,AJOH

!     ***********************************************************

      SUBROUTINE TRGRP3(INQ)

      USE TRCOMM, ONLY : GVR, NGR, NRMAX, NRMP
      IMPLICIT NONE
      INTEGER,      INTENT(IN):: INQ

      CALL PAGES

      CALL TRGR2D( 1.0,12.0,10.0,17.0,GVR(1,1, 9),NRMP,NRMAX+1,NGR, &
     &            '@QP@',INQ)

      CALL TRGR2D(13.5,24.5,10.0,17.0,GVR(1,1,10),NRMP,NRMAX,NGR, &
     &            '@AJ [MA/m^2]@',INQ)

      CALL TRGR2D( 1.0,12.0, 1.0, 8.0,GVR(1,1,11),NRMP,NRMAX,NGR, &
     &            '@EZOH [V/m]@',INQ)

      CALL TRGR2D(13.5,24.5, 1.0, 8.0,GVR(1,1,12),NRMP,NRMAX,NGR, &
     &            '@AJOH [MA/m^2]@',INQ)

      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRP3

!     ***********************************************************

!           THREE DIMENSIONAL : AJNB,AJBS,WB,WF

!     ***********************************************************

      SUBROUTINE TRGRP4(INQ)

      USE TRCOMM, ONLY : GVR, NGR, NRMAX, NRMP
      IMPLICIT NONE
      INTEGER,      INTENT(IN):: INQ

      CALL PAGES

      CALL TRGR2D( 1.0,12.0,10.0,17.0,GVR(1,1,13),NRMP,NRMAX,NGR, &
     &            '@AJNB+AJRF [MA/m^2]@',INQ)

      CALL TRGR2D(13.5,24.5,10.0,17.0,GVR(1,1,14),NRMP,NRMAX,NGR, &
     &            '@AJBS [MA/m^2]@',INQ)

      CALL TRGR2D( 1.0,12.0, 1.0, 8.0,GVR(1,1,15),NRMP,NRMAX,NGR, &
     &            '@PIN [MW/m^3]@',INQ)

      CALL TRGR2D(13.5,24.5, 1.0, 8.0,GVR(1,1,16),NRMP,NRMAX,NGR, &
     &            '@POH [MW/m^3]@',INQ)

      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRP4

!     ***********************************************************

!           THREE DIMENSIONAL : PIN,POH,PNB,PNF

!     ***********************************************************

      SUBROUTINE TRGRP5(INQ)

      USE TRCOMM, ONLY : GVR, NGR, NRMAX, NRMP
      IMPLICIT NONE
      INTEGER,      INTENT(IN):: INQ


      CALL PAGES

      CALL TRGR2D( 1.0,12.0,10.0,17.0,GVR(1,1,17),NRMP,NRMAX,NGR, &
     &            '@WB [MJ/m^3]@',INQ)

      CALL TRGR2D(13.5,24.5,10.0,17.0,GVR(1,1,18),NRMP,NRMAX,NGR, &
     &            '@WF [MJ/m^3]@',INQ)

      CALL TRGR2D( 1.0,12.0, 1.0, 8.0,GVR(1,1,19),NRMP,NRMAX,NGR, &
     &            '@PNB [MW/m^3]@',INQ)

      CALL TRGR2D(13.5,24.5, 1.0, 8.0,GVR(1,1,20),NRMP,NRMAX,NGR, &
     &            '@PNF [MW/m^3]@',INQ)

      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRP5
