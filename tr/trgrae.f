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
      CHARACTER K2*1
C
      IF(K2.EQ.'1') CALL EQGC2D
      IF(K2.EQ.'2') CALL EQGC1D
      IF(K2.EQ.'3') CALL EQGS2D
      IF(K2.EQ.'4') CALL EQGS1D(0)
      IF(K2.EQ.'5') CALL EQGS1D(1)
      IF(K2.EQ.'6') CALL EQGSDD
      IF(K2.EQ.'7') CALL EQGC1M
      IF(K2.EQ.'8') CALL TRGRE8(INQ)
      IF(K2.EQ.'9') CALL TRGRE9(INQ)
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
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO NR=1,NRMAX
         GYR(NR,1) = GCLIP(QRHO(NR))
         GYR(NR,2) = GCLIP(QP(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &            '@QRHO,QP vs RHO@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GCLIP(BPRHO(NR))
         GYR(NR,2) = GCLIP(BP(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &            '@BPRHO,BP [T] vs RHO@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GCLIP(DVRHO(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@DVRHO vs RHO@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GCLIP(DSRHO(NR))
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
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO NR=1,NRMAX
         GYR(NR,1) = GCLIP(ABRHO(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@ABRHO vs RHO@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GCLIP(ARRHO(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@ARRHO vs RHO@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GCLIP(AR1RHO(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@AR1RHO vs RHO@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GCLIP(AR2RHO(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@AR2RHO vs RHO@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
