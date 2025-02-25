C     $Id$
C
C vmecF modelg=4
C      PARAMETER (NSD=31,MPOL=6,NMAX=3,NTHETA=18,NZETA=12)
C vmecS modelg=5
C      PARAMETER (NSD=31,MPOL=8,NMAX=8,NTHETA=24,NZETA=32)
C vmecT modelg=6
C      PARAMETER (NSD=61,MPOL=8,NMAX=6,NTHETA=32,NZETA=32)
C max
C      PARAMETER (NSD=61,MPOL=8,NMAX=8,NTHETA=32,NZETA=32)
C
      PARAMETER (NSD=61,MPOL=8,NMAX=6,NTHETA=22,NZETA=16)
C
      PARAMETER (NMIN=-NMAX,NTOR=NMAX-NMIN+1)
      PARAMETER (MND=MPOL*NTOR,NTZ=NTHETA*NZETA)
      PARAMETER (MNSD=MND*NSD)
C      PARAMETER (NTHPTS=48,NRT=NSD*NTHPTS,NPL=NRT-NTHPTS)
      PARAMETER (NTHPTS=49,NRT=NSD*NTHPTS,NPL=NRT-NTHPTS)
C      PARAMETER (NSRM=62,NMNM=256,NSMNM=NSRM*NMNM)
C     PARAMETER (NSRM=nsd,NMNM=mpol*(2*ntor+1)-ntor,NSMNM=NSRM*NMNM) 
C                       ! nmnm = 2**p  ??
C
      PARAMETER (NSRM=60,NMNM=276)
      PARAMETER (NSRMP=NSRM+3,NSMNM=NSRMP*NMNM)
C
      COMMON /WMHVR1/ NSRMAX,MNMAX
      COMMON /WMHVR2/ XM(NSMNM),XN(NSMNM),RMNC(NSMNM),ZMNS(NSMNM)
      COMMON /WMHVR3/ RLMNS(NSMNM),BMOD(NSMNM),RGMOD(NSMNM)
      COMMON /WMHVR4/ BSU(NSMNM),BSV(NSMNM)
      COMMON /WMHVR5/ RIOTAS(NSRMP),RMASS(NSRMP),PRES(NSRMP)
      COMMON /WMHVR6/ PHIPS(NSRMP),BPCO(NSRMP),BACO(NSRMP)
      COMMON /WMHVR7/ PHI(NSRMP),RCHI(NSRMP),VP(NSRMP)
      COMMON /WMHVR8/ RJTHETA(NSRMP),RJZETA(NSRMP),SPECW(NSRMP)
      COMMON /WMHVR9/ XS(NSRMP),XSH(NSRMP),XSHRHO(NSRMP)
      COMMON /WMHBZ1/ BBOZH(NMNM,NSRMP)
      COMMON /WMHBZ2/ RBOZH(NMNM,NSRMP)
      COMMON /WMHBZ3/ ZBOZH(NMNM,NSRMP)
      COMMON /WMHBZ4/ PBOZH(NMNM,NSRMP)
      COMMON /WMHBZ5/ SHALF(NSRMP),WJS(NSRMP),WIS(NSRMP)
C
      COMMON /WMCEF1/ BSQ1(NTZ,NSD),PHIMOD(MNSD),PHIMOD2(MNSD)
      COMMON /WMCEF2/ LMNS(MNSD),BETAP(NSD),BPOL(NSD),BTOR1(NSD)
      COMMON /WMCEF3/ BZCO(NSD)
      COMMON /WMCEF4/ UB(NSD),IXM(MND),IXN(MND)
C
      COMMON /WMVOL1/ RJTOR(NSD)
      COMMON /WMVOL2/ RJROL(NSD)
      COMMON /WMVOL3/ RGAM,VOLI
C
      COMMON /WMPLT1/ FSQT(100),WDOT(100),AMVM(0:10)
      COMMON /WMHRD1/ MN0,NTOR0
      COMMON /WMHRD2/ HS,OHS,DNORM,RC
C
      COMMON /WMHFR1/ YRBS(NSRMP),YZBS(NSRMP)
      COMMON /WMHFR2/ R1DATA(12),R2DATA(11,NSRMP),R3DATA(2,100)
C
      COMMON /WMHSP1/ SRMNC(NMNM,NRM),DRMNC(NMNM,NRM)
      COMMON /WMHSP2/ SZMNS(NMNM,NRM),DZMNS(NMNM,NRM)
      COMMON /WMHSP3/ SPMNS(NMNM,NRM),DPMNS(NMNM,NRM)
      COMMON /WMHSP4/ SBMNC(NMNM,NRM),DBMNC(NMNM,NRM)
      COMMON /WMHSP5/ RMNCC(NMNM,NSRMP),ZMNSS(NMNM,NSRMP)
C
      COMMON /WMHFT1/ BSTH(NMNM,NRM),BSPH(NMNM,NRM)
      COMMON /WMHFT2/ U1(4,NSRMP,NMNM),FX1(NSRMP)
      COMMON /WMHFT3/ U2(4,NSRMP,NMNM),FX2(NSRMP)
      COMMON /WMHFT4/ U3(4,NSRMP,NMNM),FX3(NSRMP)
      COMMON /WMHFT5/ U4(4,NSRMP,NMNM),FX4(NSRMP)
      COMMON /WMHFT6/ U5(4,NSRMP),FX5(NSRMP)
      COMMON /WMHFT7/ U6(4,NSRMP),FX6(NSRMP)
      COMMON /WMHFT8/ BSTHSV(NMNM),BSTHSD(NMNM)
      COMMON /WMHFT9/ BSPHSV(NMNM),BSPHSD(NMNM)
