C     $Id$
C
      PARAMETER (NRGM=101,NZGM=101)
      PARAMETER (NPSM=101)
C
      COMMON /EQCNS1/ PI,RMU0,AEE,AMP,BLTZ,AN0,RGAS
C
      COMMON /EQPRM1/ RR,BB,RIP
      COMMON /EQPRM2/ RA,RKAP,RDLT,RB
      COMMON /EQPRM3/ PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2
      COMMON /EQPRM4/ PP0,PP1,PP2,PROFP0,PROFP1,PROFP2
      COMMON /EQPRM6/ PT0,PT1,PT2,PROFT0,PROFT1,PROFT2,PTS,PN0
      COMMON /EQPRM7/ PV0,PV1,PV2,PROFV0,PROFV1,PROFV2
      COMMON /EQPRM8/ PROFR0,PROFR1,PROFR2,RHOITB
      COMMON /EQPRM9/ EPSEQ
      COMMON /EQPRN1/ NSGMAX,NTGMAX
      COMMON /EQPRN2/ NRGMAX,NZGMAX,NPSMAX
      COMMON /EQPRN3/ NRMAX,NTHMAX,NSUMAX
C
      COMMON /EQGLB1/ RAXIS,ZAXIS,SAXIS
      COMMON /EQGLB2/ PVOL,RAAVE,BETAT,BETAP,QAXIS,QSURF
C
      COMMON /EQOUT1/ RG(NRGM),ZG(NZGM),PSIRZ(NRGM,NZGM)
      COMMON /EQOUT2/ PSIPS(NPSM),PPPS(NPSM),TTPS(NPSM)
      COMMON /EQOUT3/ TEPS(NPSM),OMPS(NPSM)
      COMMON /EQRSP1/ URZ(4,4,NRGM,NZGM)
C
      COMMON /EQNAM1/ KNAMEQ
      CHARACTER KNAMEQ*32
