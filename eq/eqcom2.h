C     $Id$
C
      PARAMETER (NSGM=32,NTGM=32)
C
      PARAMETER (MLM=NSGM*NTGM,MWM=4*NTGM-1)
      PARAMETER (NSGMP=NSGM+1,NTGMP=NTGM+1)
      PARAMETER (NSGPM=NSGM+2,NTGPM=NTGM+2)
      PARAMETER (NXM=2*NSGM)
C
      COMMON /EQMSH1/ DSG,DTG
      COMMON /EQMSH2/ SIGM(NSGM),THGM(NTGM),RHOM(NTGM)
      COMMON /EQMSH3/ RMG(NSGM,NTGMP),RGM(NSGMP,NTGM)
      COMMON /EQMSH4/ SIGG(NSGMP),THGG(NTGMP),RHOG(NTGMP)
      COMMON /EQMSH5/ RMM(NTGM,NSGM)
C
      COMMON /EQBND1/ Q(MWM,MLM)
      COMMON /EQCEF1/ AA(NSGMP,NTGM),AB(NSGMP,NTGM)
      COMMON /EQCEF2/ AC(NSGM,NTGMP),AD(NSGM,NTGMP)
C
      COMMON /EQPSI1/ PSI(NTGM,NSGM),PSI0
      COMMON /EQPSI2/ DELPSI(NTGM,NSGM)
C
      COMMON /EQRHS1/ HJT(NTGM,NSGM),PP(NTGM,NSGM),TT(NTGM,NSGM),TJ,RRC
      COMMON /EQRHS2/ HJP1(NTGM,NSGM),HJP2(NTGM,NSGM)
      COMMON /EQRHS3/ HJT1(NTGM,NSGM),HJT2(NTGM,NSGM)
      COMMON /EQRHS3/ RHO(NTGM,NSGM)
C
      COMMON /EQSPL1/ NSGPMAX,NTGPMAX
      COMMON /EQSPL2/ SIGMX(NSGPM),THGMX(NTGPM)
      COMMON /EQSPL3/ PSIX(NTGPM,NSGPM)
      COMMON /EQSPL5/ U(4,4,NTGPM,NSGPM)
C
      COMMON /EQZBR1/ ZBRF
C
      COMMON /EQPRF1/ XAX(NXM),RPSI(NXM),RPP(NXM),RTT(NXM)
      COMMON /EQPRF2/ RHJP(NXM),RHJT(NXM)
