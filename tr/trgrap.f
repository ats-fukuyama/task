C     $Id$
C  
C     ***********************************************************
C
C           GRAPHIC : CONTROL ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TRGRG0(K2,INQ)
C
      INCLUDE 'trcomm.inc'
      CHARACTER K2*1
C
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
      END
C  
C     ***********************************************************
C
C           GRAPHIC : CONTROL ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TRGRP0(K2,INQ)
C
      CHARACTER K2*1
C
      IF(K2.EQ.'1') CALL TRGRP1(INQ)
      IF(K2.EQ.'2') CALL TRGRP2(INQ)
      IF(K2.EQ.'3') CALL TRGRP3(INQ)
      IF(K2.EQ.'4') CALL TRGRP4(INQ)
      IF(K2.EQ.'5') CALL TRGRP5(INQ)
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : NE,ND,NT,NA
C
C     ***********************************************************
C
      SUBROUTINE TRGRG1(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      CALL TRGR1D( 3.0,12.0,11.0,17.0,
     &            GRM,GVR(1,1,1),NRMP,NRMAX,NGR,
     &            '@NE [10^20/m^3]  vs r@',0+INQ)
C
      CALL TRGR1D(15.5,24.5,11.0,17.0,
     &            GRM,GVR(1,1,2),NRMP,NRMAX,NGR,
     &            '@ND [10^20/m^3]  vs r@',0+INQ)
C
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,
     &            GRM,GVR(1,1,5),NRMP,NRMAX,NGR,
     &            '@TE [keV]  vs r@',0+INQ)
C
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,
     &            GRM,GVR(1,1,6),NRMP,NRMAX,NGR,
     &            '@TD [keV]  vs r@',0+INQ)
C
      CALL PAGEE
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : TE,TD,TT,TA
C
C     ***********************************************************
C
      SUBROUTINE TRGRG2(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      CALL TRGR1D( 3.0,12.0,11.0,17.0,
     &            GRM,GVR(1,1,3),NRMP,NRMAX,NGR,
     &            '@NT [10^20/m^3]  vs r@',0+INQ)
C
      CALL TRGR1D(15.5,24.5,11.0,17.0,
     &            GRM,GVR(1,1,4),NRMP,NRMAX,NGR,
     &            '@NA [10^20/m^3]  vs r@',0+INQ)
C
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,
     &            GRM,GVR(1,1,7),NRMP,NRMAX,NGR,
     &            '@TT [keV]  vs r@',0+INQ)
C
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,
     &            GRM,GVR(1,1,8),NRMP,NRMAX,NGR,
     &            '@TA [keV]  vs r@',0+INQ)
C
      CALL PAGEE
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : Q,J,EZOH,JOH
C
C     ***********************************************************
C
      SUBROUTINE TRGRG3(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      CALL TRGR1D( 3.0,12.0,11.0,17.0,
     &            GRG,GVR(1,1,9),NRMP,NRMAX+1,NGR,
     &           '@QP  vs r@',0+INQ)
C
      CALL TRGR1D(15.5,24.5,11.0,17.0,
     &            GRM,GVR(1,1,10),NRMP,NRMAX,NGR,
     &           '@AJ [MA/m^2]  vs r@',0+INQ)
C
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,
     &            GRM,GVR(1,1,11),NRMP,NRMAX,NGR,
     &           '@EZ [V/m]  vs r@',0+INQ)
C
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,
     &            GRM,GVR(1,1,12),NRMP,NRMAX,NGR,
     &           '@AJOH [MA/m^2]  vs r@',0+INQ)
C
      CALL PAGEE
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : JNB,JBS,WB,WF
C
C     ***********************************************************
C
      SUBROUTINE TRGRG4(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      CALL TRGR1D( 3.0,12.0,11.0,17.0,
     &            GRM,GVR(1,1,13),NRMP,NRMAX,NGR,
     &            '@AJNB+AJRF [MA/m^2]  vs r@',0+INQ)
C
      CALL TRGR1D(15.5,24.5,11.0,17.0,
     &            GRM,GVR(1,1,14),NRMP,NRMAX,NGR,
     &            '@AJBS [MA/m^2]  vs r@',0+INQ)
C
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,
     &            GRM,GVR(1,1,15),NRMP,NRMAX,NGR,
     &            '@PIN [MW/m^3]  vs r@',0+INQ)
C
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,
     &            GRM,GVR(1,1,16),NRMP,NRMAX,NGR,
     &            '@POH [MW/m^3]  vs r@',0+INQ)
C
      CALL PAGEE
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : G,s,alpha,chiD
C
C     ***********************************************************
C
      SUBROUTINE TRGRG5(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      CALL TRGR1D( 3.0,12.0,11.0,17.0,
     &            GRM,GVR(1,1,17),NRMP,NRMAX,NGR,
     &            '@s  vs r@',0+INQ)
C
      CALL TRGR1D(15.5,24.5,11.0,17.0,
     &            GRM,GVR(1,1,18),NRMP,NRMAX,NGR,
     &            '@G  vs r@',0+INQ)
C
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,
     &            GRM,GVR(1,1,19),NRMP,NRMAX,NGR,
     &            '@alpha  vs r@',0+INQ)
C
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,
     &            GRM,GVR(1,1,20),NRMP,NRMAX-1,NGR,
     &            '@AKD [m^2/s]  vs r@',0+INQ)
C
      CALL PAGEE
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : s-alpha,s,alpha,chiD
C
C     ***********************************************************
C
      SUBROUTINE TRGRG6(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      CALL TRGR1DX( 3.0,12.0,11.0,17.0,
     &            GVR(1,1,19),GVR(1,1,18),NRMP,NRMAX,NGR,
     &            '@s vs alpha@',0+INQ)
C
      CALL TRGR1D(15.5,24.5,11.0,17.0,
     &            GRM,GVR(1,1,18),NRMP,NRMAX,NGR,
     &            '@s  vs r@',0+INQ)
C
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,
     &            GRM,GVR(1,1,19),NRMP,NRMAX,NGR,
     &            '@alpha  vs r@',0+INQ)
C
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,
     &            GRM,GVR(1,1,20),NRMP,NRMAX-1,NGR,
     &            '@AKD [m^2/s]  vs r@',0+INQ)
C
      CALL PAGEE
      RETURN
      END
C
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : BP,PSI
C
C     ***********************************************************
C
      SUBROUTINE TRGRG7(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      CALL TRGR1D( 3.0,12.0,11.0,17.0,
     &            GRG,GVR(1,1,21),NRMP,NRMAX,NGR,
     &            '@BP vs r@',0+INQ)
C
      CALL TRGR1D(15.5,24.5,11.0,17.0,
     &            GRM,GVR(1,1,22),NRMP,NRMAX,NGR,
     &            '@PSI vs r@',0+INQ)
C
c$$$      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,
c$$$     &            GRM,GVR(1,1,23),NRMP,NRMAX,NGR,
c$$$     &            '@alpha  vs r@',0+INQ)
c$$$C
c$$$      CALL TRGR1D(15.5,24.5, 2.0, 8.0,
c$$$     &            GRM,GVR(1,1,24),NRMP,NRMAX-1,NGR,
c$$$     &            '@AKD [m^2/s]  vs r@',0+INQ)
C
      CALL PAGEE
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : Pressure,J(r),s(r),chi(r)
C
C     ***********************************************************
C
      SUBROUTINE TRGRGA(INQ)
C
      INCLUDE 'trcomm.inc'
      DIMENSION GWR(NRMP,NGM)
C
      CALL PAGES
C
      DO NG=1,NGR
      DO NR=1,NRMAX
         GWR(NR,NG) =  GVR(NR,NG,1)*GVR(NR,NG,5)*GUCLIP(RKEV*1.D14)
     &              +  GVR(NR,NG,2)*GVR(NR,NG,6)*GUCLIP(RKEV*1.D14)
      ENDDO
      ENDDO
C
      CALL TRGR1D( 3.0,12.0,11.0,17.0,
     &            GRM,GWR,NRMP,NRMAX,NGR,
     &            '@P [MP]  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GWR(NR,1) =  GVR(NR,1,10)
         GWR(NR,2) =  GVR(NR,1,12)
         GWR(NR,3) =  0.0
         GWR(NR,4) =  GVR(NR,1,14)
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,
     &            GRM,GWR,NRMP,NRMAX,4,
     &            '@J,JOH,JBS [MA/m^2]  vs r@',2+INQ)
C
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,
     &            GRM,GVR(1,1,18),NRMP,NRMAX,NGR,
     &            '@s  vs r@',2+INQ)
C
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,
     &            GRM,GVR(1,1,20),NRMP,NRMAX-1,NGR,
     &            '@AKD [m^2/s]  vs r@',2+INQ)
C
      CALL PAGEE
      RETURN
      END
C
C     ***********************************************************
C
C           THREE DIMENSIONAL : DENSITY PROFILE
C
C     ***********************************************************
C
      SUBROUTINE TRGRP1(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      CALL TRGR2D( 1.0,12.0,10.0,17.0,GVR(1,1,1),NRMP,NRMAX,NGR,
     &            '@NE [10^20/m^3]@',INQ)
C
      CALL TRGR2D(13.5,24.5,10.0,17.0,GVR(1,1,2),NRMP,NRMAX,NGR,
     &            '@ND [10^20/m^3]@',INQ)
C
      CALL TRGR2D( 1.0,12.0, 1.0, 8.0,GVR(1,1,3),NRMP,NRMAX,NGR,
     &            '@NT [10^20/m^3]@',INQ)
C
      CALL TRGR2D(13.5,24.5, 1.0, 8.0,GVR(1,1,4),NRMP,NRMAX,NGR,
     &            '@NA [10^20/m^3]@',INQ)
C
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           THREE DIMENSIONAL : TEMPERATURE PROFILE
C
C     ***********************************************************
C
      SUBROUTINE TRGRP2(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      CALL TRGR2D( 1.0,12.0,10.0,17.0,GVR(1,1,5),NRMP,NRMAX,NGR,
     &            '@TE [keV]@',INQ)
C
      CALL TRGR2D(13.5,24.5,10.0,17.0,GVR(1,1,6),NRMP,NRMAX,NGR,
     &            '@TD [keV]@',INQ)
C
      CALL TRGR2D( 1.0,12.0, 1.0, 8.0,GVR(1,1,7),NRMP,NRMAX,NGR,
     &            '@TT [keV]@',INQ)
C
      CALL TRGR2D(13.5,24.5, 1.0, 8.0,GVR(1,1,8),NRMP,NRMAX,NGR,
     &            '@TA [keV]@',INQ)
C
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           THREE DIMENSIONAL : QP,AJ,EZOH,AJOH
C
C     ***********************************************************
C
      SUBROUTINE TRGRP3(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      CALL TRGR2D( 1.0,12.0,10.0,17.0,GVR(1,1, 9),NRMP,NRMAX+1,NGR,
     &            '@QP@',INQ)
C
      CALL TRGR2D(13.5,24.5,10.0,17.0,GVR(1,1,10),NRMP,NRMAX,NGR,
     &            '@AJ [MA/m^2]@',INQ)
C
      CALL TRGR2D( 1.0,12.0, 1.0, 8.0,GVR(1,1,11),NRMP,NRMAX,NGR,
     &            '@EZOH [V/m]@',INQ)
C
      CALL TRGR2D(13.5,24.5, 1.0, 8.0,GVR(1,1,12),NRMP,NRMAX,NGR,
     &            '@AJOH [MA/m^2]@',INQ)
C
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           THREE DIMENSIONAL : AJNB,AJBS,WB,WF
C
C     ***********************************************************
C
      SUBROUTINE TRGRP4(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      CALL TRGR2D( 1.0,12.0,10.0,17.0,GVR(1,1,13),NRMP,NRMAX,NGR,
     &            '@AJNB+AJRF [MA/m^2]@',INQ)
C
      CALL TRGR2D(13.5,24.5,10.0,17.0,GVR(1,1,14),NRMP,NRMAX,NGR,
     &            '@AJBS [MA/m^2]@',INQ)
C
      CALL TRGR2D( 1.0,12.0, 1.0, 8.0,GVR(1,1,15),NRMP,NRMAX,NGR,
     &            '@PIN [MW/m^3]@',INQ)
C
      CALL TRGR2D(13.5,24.5, 1.0, 8.0,GVR(1,1,16),NRMP,NRMAX,NGR,
     &            '@POH [MW/m^3]@',INQ)
C
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           THREE DIMENSIONAL : PIN,POH,PNB,PNF
C
C     ***********************************************************
C
      SUBROUTINE TRGRP5(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      CALL TRGR2D( 1.0,12.0,10.0,17.0,GVR(1,1,17),NRMP,NRMAX,NGR,
     &            '@WB [MJ/m^3]@',INQ)
C
      CALL TRGR2D(13.5,24.5,10.0,17.0,GVR(1,1,18),NRMP,NRMAX,NGR,
     &            '@WF [MJ/m^3]@',INQ)
C
      CALL TRGR2D( 1.0,12.0, 1.0, 8.0,GVR(1,1,19),NRMP,NRMAX,NGR,
     &            '@PNB [MW/m^3]@',INQ)
C
      CALL TRGR2D(13.5,24.5, 1.0, 8.0,GVR(1,1,20),NRMP,NRMAX,NGR,
     &            '@PNF [MW/m^3]@',INQ)
C
      CALL PAGEE
C
      RETURN
      END
