C
      COMMON /WMPRL1/ PNA,PNAL,PTA,ZEFF
      COMMON /WMPRL2/ NRMAX,NTHMAX,NHHMAX
      COMMON /WMPRL3/ NTH0,NPH0,NHC
      COMMON /WMPRL4/ NPRINT,NGRAPH
      COMMON /WMPRL5/ MODELJ,MODELA,MODELK,MODELM,MODELW
      COMMON /WMPRL6/ PRFIN,NSUMAX,NSWMAX
      COMMON /WMPRL7/ MODELPR(NRM,NSM),MODELVR(NRM,NSM)
      COMMON /WMPRL8/ MDLWMF,MDLWMX,MDLWMD,MDLWM_NPHS,NPHSMAX
      COMMON /WMPRL9/ PFRACS(NPHSM),NPH0S(NPHSM)
C
      COMMON /WMANT1/ CRF,RD,BETAJ,ANTANG,MWGMAX
      COMMON /WMANT2/ AJ(NAM),APH(NAM),NAMAX
      COMMON /WMANT3/ THJ1(NAM),THJ2(NAM),PHJ1(NAM),PHJ2(NAM)
C
      COMMON /WMNOD1/ XR(NRM),XRHO(NRM)
      COMMON /WMNOD2/ MDSIZ,MDMIN,MDMAX,LDSIZ,LDMIN,LDMAX,MDSIZX
      COMMON /WMNOD3/ NDSIZ,NDMIN,NDMAX,KDSIZ,KDMIN,KDMAX,NDSIZX
      COMMON /WMNOD4/ CJANT(3,MDM,NDM),MODEWG
C
C      COMMON /WMMTX1/ CEM(MBNDM,MSIZP),CFV(MSIZP)
C      COMMON /WMMTX2/ MSIZ,MBND
C      COMMON /WMMTX3/ CFVG(MSIZM)
C
C      COMMON /WMTNS1/ CTNSR(3,3,-MDMX:MDMX,-NDMX:NDMX,MDM,NDM)
C      COMMON /WMTNS2/ CGD(3,3,MDM,MDM,NDM,NDM,3)
C
      COMMON /WMEQG1/ RG11(MDM,NDM,NRM),RG12(MDM,NDM,NRM)
      COMMON /WMEQG2/ RG13(MDM,NDM,NRM),RG22(MDM,NDM,NRM)
      COMMON /WMEQG3/ RG23(MDM,NDM,NRM),RG33(MDM,NDM,NRM)
      COMMON /WMEQG4/ RJ(MDM,NDM,NRM)
C      COMMON /WMEQG5/ CGF11(MDM,NDM,3),CGF12(MDM,NDM,3)
C      COMMON /WMEQG6/ CGF13(MDM,NDM,3),CGF22(MDM,NDM,3)
C      COMMON /WMEQG7/ CGF23(MDM,NDM,3),CGF33(MDM,NDM,3)
C
      COMMON /WMEQV1/ DVOL(MDM,NDM,NRM),DVOLS(NRM)
C
      COMMON /WMFLD1/ CEFLD(3,MDM,NDM,NRM),CEFLDK(3,MDM,NDM,NRM)
      COMMON /WMFLD2/ CBFLD(3,MDM,NDM,NRM),CBFLDK(3,MDM,NDM,NRM)
      COMMON /WMFLD3/ CEN(3,MDM,NDM,NRM),CEP(3,MDM,NDM,NRM)
      COMMON /WMFLD4/ CBN(3,MDM,NDM,NRM),CBP(3,MDM,NDM,NRM)
C
      COMMON /WMPWR1/ PABS(MDM,NDM,NRM,NSM),PABSK(MDM,NDM,NRM,NSM)
      COMMON /WMPWR2/ PABSR(NRM,NSM),PABSKT(MDM,NDM,NSM)
      COMMON /WMPWR3/ CFLX(MDM,NDM,NRM),CFLXK(MDM,NDM,NRM),CFLXR(NRM)
      COMMON /WMPWR4/ CRADK(MDM,NDM,NAM),CRADKT(MDM,NDM)
      COMMON /WMPWR5/ CRADT(NAM),CRADTT,PABST(NSM),PABSTT
      COMMON /WMPWR6/ PCUR(MDM,NDM,NRM),PCURR(NRM),PCURT
      COMMON /WMPWR7/ CPABS(MDM,MDM,NDM,NDM,NSM,NRM)
      COMMON /WMPWR8/ PTABS(MDM,NDM,NRM,NSM)
C
      COMMON /WMEQX1/ RPST(NTHM,NHHM,NRM),ZPST(NTHM,NHHM,NRM)
      COMMON /WMEQX2/ PPST(NTHM,NHHM,NRM),BPST(NTHM,NHHM,NRM)
      COMMON /WMEQX3/ BFLD(3,NTHM,NHHM,NRM)
      COMMON /WMEQX4/ RGMIN,RGMAX,ZGMIN,ZGMAX
      COMMON /WMEQX5/ RSU(NSUM,NDM),ZSU(NSUM,NDM)
      COMMON /WMEQX6/ RSW(NSUM,NDM),ZSW(NSUM,NDM)
C
      COMMON /WMEQU1/ RHOT(NRM),PSIP(NRM),PSIPA,PSITA
      COMMON /WMEQU2/ RPS(NTHM,NRM),ZPS(NTHM,NRM)
      COMMON /WMEQU4/ DRPSI(NTHM,NRM),DZPSI(NTHM,NRM)
      COMMON /WMEQU5/ DRCHI(NTHM,NRM),DZCHI(NTHM,NRM)
      COMMON /WMEQU6/ PPS(NRM),QPS(NRM),RBPS(NRM),VPS(NRM),RLEN(NRM)
      COMMON /WMEQU7/ BPR(NTHM,NRM),BPZ(NTHM,NRM)
      COMMON /WMEQU8/ BPT(NTHM,NRM),BTP(NTHM,NRM)
      COMMON /WMEQU9/ SIOTA(NRM)
C
C      COMMON /EQPRM4/ BTV,TTCU,TTPR,BETS,BETA,BETJ
C      COMMON /EQPRM5/ SAXIS,RAXIS,ZAXIS,ELL,TRG
      COMMON /EQPRM/ RAXIS,ZAXIS
C
      COMMON /WMFFT1/ CT(MDM*2),RFFT(MDM*2),LFFT(MDM*2)
C
      COMMON /WMEGN1/ FRMIN,FRMAX,FIMIN,FIMAX,FI0
      COMMON /WMEGN2/ FRINI,FIINI
      COMMON /WMEGN3/ NGFMAX,NGXMAX,NGYMAX
      COMMON /WMEGN4/ SCMIN,SCMAX,NSCMAX
      COMMON /WMEGN5/ LISTEG
      COMMON /WMEGN6/ DLTNW,EPSNW,LMAXNW,LISTNW,MODENW
C
      COMMON /WMTAE1/ WAEMIN,WAEMAX
C
      COMMON /WMPRF1/ RNPRF(NRM,NSM),RTPRF(NRM,NSM)
C
      COMMON /WMEGP1/ AMPEIGEN,MODEEG
C
      COMMON /WMGRA1/ KACONT(8,NRGM,NTHGM)
      COMMON /WMGRA2/ RPSG(NTHGM,NRM),ZPSG(NTHGM,NRM)
C
C      COMMON /WMKNM1/ KNAMEQ,KNAMPF
C      CHARACTER KNAMEQ*80,KNAMPF*80
C
      COMMON /WMNAM1/ KNAMEQ_SAVE
      CHARACTER KNAMEQ_SAVE*80
      COMMON /WMNAM2/ KNAMTR_SAVE
      CHARACTER KNAMTR_SAVE*80
