C     $Id$
C  
C     ***********************************************************
C
C           GRAPHIC : CONTROL ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TRGRE0(K2,INQ)
C
      INCLUDE 'trcomm.inc'
      CHARACTER K2*1
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      IF(K2.EQ.'1') CALL EQGC2D
      IF(K2.EQ.'2') CALL EQGC1D
      IF(K2.EQ.'3') CALL EQGS2D
      IF(K2.EQ.'4') CALL EQGS1D(0)
      IF(K2.EQ.'5') CALL EQGS1D(1)
      IF(K2.EQ.'6') CALL EQGSDD
      IF(K2.EQ.'7') CALL EQGC1M
      IF(K2.EQ.'8') CALL TRGRE8(INQ)
      IF(K2.EQ.'9') CALL TRGRE9(INQ)
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : VGR
C
C     ***********************************************************
C
      SUBROUTINE TRGRE8(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(QRHO(NR))
         GYR(NR,2) = GUCLIP(QP(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &            '@QRHO,QP vs RHO@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(BPRHO(NR))
         GYR(NR+1,2) = GUCLIP(BP(NR))
      ENDDO
      GYR(1,1) = 0.0
      GYR(1,2) = 0.0
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,2,
     &            '@BPRHO,BP [T] vs RHO@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(DVRHO(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@DVRHO vs RHO@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(DVRHO(NR)/(2.D0*PI*RR))
      ENDDO
C
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@DSRHO vs RHO@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : AD,AV,AVK,TAUB,TAUF
C
C     ***********************************************************
C
      SUBROUTINE TRGRE9(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(ABRHO(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@ABRHO vs RHO@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(ARRHO(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@ARRHO vs RHO@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(AR1RHO(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@AR1RHO vs RHO@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(AR2RHO(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@AR2RHO vs RHO@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
