C     $Id$
C
      PARAMETER (NRGM=101,NZGM=101)
      PARAMETER (NPSM=101)
C
      COMMON /EQCNS1/ PI,RMU0
C
      COMMON /EQPRM1/ RR,BB,RIP
      COMMON /EQPRM2/ RA,RKAP,RDLT,RB
      COMMON /EQPRM3/ PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2
      COMMON /EQPRM4/ PP0,PP1,PP2,PROFP0,PROFP1,PROFP2
      COMMON /EQPRM5/ PROFR0,PROFR1,PROFR2,RHOITB
      COMMON /EQPRM6/ NSGMAX,NTGMAX
      COMMON /EQPRM7/ NRGMAX,NZGMAX,NPSMAX
      COMMON /EQPRM8/ EPSEQ
      COMMON /EQPRM9/ NRMAX,NTHMAX,NSUMAX
C
      COMMON /EQGLB1/ RAXIS,ZAXIS,SAXIS
      COMMON /EQGLB2/ PVOL,RAAVE,BETAT,BETAP,QAXIS,QSURF
C
      COMMON /EQOUT1/ RG(NRGM),ZG(NZGM),PSIRZ(NRGM,NZGM)
      COMMON /EQOUT2/ PSIPS(NPSM),PPPS(NPSM),TTPS(NPSM)
      COMMON /EQRSP1/ URZ(4,4,NRGM,NZGM)
C
      COMMON /EQNAM1/ KNAMEQ
      CHARACTER KNAMEQ*32
