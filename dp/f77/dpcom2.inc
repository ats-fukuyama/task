

      PARAMETER(NHM=200)
      PARAMETER(NXM=1001)
      PARAMETER(NPM=201,NTHM=101,NRM=201)
      PARAMETER(NGXM=201,NGYM=201)
      PARAMETER(NSAM=4)
      PARAMETER(NCHM=200)

      REAL*8 CHI,CCHI

      COMMON /DPPRM3/ RF0,RFI0,RKX0,RKY0,RKZ0,RX0,RY0,RZ0
      COMMON /DPPRM4/ RF1,RFI1,RKX1,RKY1,RKZ1,RX1
      COMMON /DPPRM5/ RF2,RFI2,RKX2,RKY2,RKZ2,RX2
      COMMON /DPPRM6/ NXMAX
      COMMON /DPPRM7/ EPSRT,LMAXRT

      COMMON /DPDSP5/ CRF0,CKX0,CKY0,CKZ0,XPOS0,YPOS0,ZPOS0,ILIST
      COMMON /DPDSP6/ CALAM(0:NHM)

      COMMON /DPFPT1/ NPMAX,NTHMAX,NRMAX
      COMMON /DPFPT2/ DELTH,DELR,RMIN,RMAX
      COMMON /DPFPT3/ FP(NTHM,NPM,NRM)
      COMMON /DPFPT4/ PN0,PT0,PTH0
      COMMON /DPDPT5/ NS_NSA(NSAM),NSAMAX
      COMMON /DPDPT6/ AEFP(NSAM),AMFP(NSAM),RNFP0(NSAM),RTFP0(NSAM)
      COMMON /DPDPT7/ DELP(NSAM),PMAX(NSAM)

      COMMON /DPFPM1/ FM(NPM,NTHM),DFP(NPM,NTHM),DFT(NPM,NTHM)
!      COMMON /DPFPM2/ DFR(NPM,NTHM)

      COMMON /DPFPM2/ THM(NTHM),THG(NTHM),RM(NRM)
      COMMON /DPFPM3/ TSNM(NTHM),TSNG(NTHM)
      COMMON /DPFPM4/ TCSM(NTHM),TCSG(NTHM)
      COMMON /DPFPM5/ TTNM(NTHM),TTNG(NTHM)
      COMMON /DPFPM6/ RGMM(NPM),RGMG(NPM)
      COMMON /DPFPM7/ DGP1(NPM,NTHM),DGP2(NPM,NTHM)
      COMMON /DPFPM8/ DGT1(NPM,NTHM),DGT2(NPM,NTHM)
      COMMON /DPFPM9/ ADJ(0:NHM),ADJD(0:NHM)
      COMMON /DPFPM10/ CHI(NCHM)
      COMMON /DPFPM11/ CCHI(NCHM),SCHI(NCHM)
      COMMON /DPFPM12/ PM(NPM,NSAM),PG(NPM,NSAM)

      COMMON /DPFP01/ RHON_MIN,RHON_MAX
      COMMON /DPFP02/ RFP(NRM+2),URFP(4,NRM+2,NPM,NTHM)
      COMMON /DPFP03/ FPDATA(3),FPDATAS(NSAM),NFPDAT(3)
      COMMON /DPFP04/ FNS(NTHM,NPM,NRM,NSAM)

      COMMON /DPFA01/ US(4,4,NPM,NRM,NTHM),PMa0(NPM),RHOa0(NRM)
      COMMON /DPFA02/ fa0(NPM,NRM,NTHM),dfpra0(NPM,NRM,NTHM)
      COMMON /DPFA03/ dfpa0(NPM,NRM,NTHM),dfra0(NPM,NRM,NTHM)
      COMMON /DPFA04/ MODEFA,NRMAXFP,NTHMAXFP