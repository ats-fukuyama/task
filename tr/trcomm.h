C     $Id$
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      PARAMETER (NRM=50,NTM=10001,NGM=101)
      PARAMETER (NSM=4,NSZM=2,N0M=2,NFM=2)
      PARAMETER (NSTM=NSM+NSZM+N0M,NEQM=3*NSTM+1)
      PARAMETER (NCGM=20,NCTM=300)
      PARAMETER (NVM=3*NSTM+1,MWM=4*NVM-1,MLM=NVM*NRM)
      PARAMETER (NRMP=NRM+1,NGLF=NRM+2)
      PARAMETER (NRMU=52)
      PARAMETER (NTURM=900,NTUM=10001)
      CHARACTER KGR1*40,KGR2*40,KGR3*40,KGR4*60
      CHARACTER KUFDEV*80,KUFDCG*80
      CHARACTER KFNLOG*40
C
C     ****** CONSTANTS ******
C
      COMMON /TRCNS1/ PI,AME,AMM,AEE,VC,AMYU0,AEPS0,RKEV
C
C     ****** PARAMETERS ******
C
      COMMON /TRPRM1/ RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RIPSS
      COMMON /TRPRM2/ PA(NSTM),PZ(NSTM),PN(NSTM),PNS(NSTM),
     &                PT(NSTM),PTS(NSTM)
      COMMON /TRPRM3/ PNC,PNFE,PNNU,PNNUS
      COMMON /TRPRM4/ PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2
      COMMON /TRPRM5/ PROFJ1,PROFJ2,ALP(3)
      COMMON /TRPRM6/ AD0,AV0,CNC,CDW(8),MDLKAI,MDLETA,MDLAD,MDLAVK
      COMMON /TRPRM7/ DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST
      COMMON /TRPRM8/ EPSLTR,LMAXTR
      COMMON /TRPRM9/ CHP,CK0,CKALFA,CKBETA,CKGUMA
      COMMON /TRPRF1/ MODELG
C
C     ****** MODEL PARAMETERS ******
C
      COMMON /TRMDL1/ TPRST,MDLST,MDLNF,IZERO
      COMMON /TRMDL2/ PNBTOT,PNBR0,PNBRW,PNBENG,PNBRTG,MDLNB
      COMMON /TRMDL3/ PECTOT,PECR0,PECRW,PECTOE,PECNPR,MDLEC
      COMMON /TRMDL4/ PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,MDLLH
      COMMON /TRMDL5/ PICTOT,PICR0,PICRW,PICTOE,PICNPR,MDLIC
      COMMON /TRMDL6/ PNBCD,PECCD,PLHCD,PICCD,PBSCD,MDLCD
      COMMON /TRMDL7/ PELTOT,PELR0,PELRW,PELRAD,PELVEL
      COMMON /TRMDL8/ PELTIM,PELPAT(NSTM),MDLPEL,MDLJBS
C
C     ****** CONTROL VARIABLES ******
C
      COMMON /TRCTL1/ T,TST,TPRE,WPPRE,RIP,NT
      COMMON /TRCTL2/ DR,FKAP,PNSS(NSTM),RKAPS,RHOA,NRAMAX,NROMAX
      COMMON /TRCTL3/ VSEC
      COMMON /TRCTL4/ GTCPU1
C
C     ****** MATRIX VARIBALES ******
C
      COMMON /TRMTX1/ XV(NVM,NRM),YV(NFM,NRM)
      COMMON /TRMTX2/ AX(MWM,MLM),X(MLM),AY(NFM,NRM),Y(NFM,NRM)
C     
C     ****** FUNDAMENTAL VARIABLES ******
C
      COMMON /TRVAR1/ RG(NRM),RM(NRM)
      COMMON /TRVAR2/ RN(NRM,NSTM),RT(NRM,NSTM),RU(NRM,NSTM)
      COMMON /TRVAR3/ RW(NRM,NFM),BP(NRM)
C
C     ****** PROFILE VARIABLES ******
C
      COMMON /TRPFL1/ RNF(NRM,NFM),RTF(NRM,NFM)
      COMMON /TRPFL2/ ANC(NRM),ANFE(NRM),ANNU(NRM)
      COMMON /TRPFL3/ ZEFF(NRM),PZC(NRM),PZFE(NRM)
      COMMON /TRPFL4/ BETA(NRM),BETAP(NRM)
      COMMON /TRPFL5/ BETAL(NRM),BETAPL(NRM),BETAQ(NRM)
C
C     ****** SOURCE VARIABLES ******
C
      COMMON /TRSRC1/ PIN(NRM,NSTM),SSIN(NRM,NSTM)
      COMMON /TRSRC2/ AJ(NRM),AJOH(NRM),EZOH(NRM),QP(NRM)
      COMMON /TRSRC3/ AJNB(NRM),AJRF(NRM),AJBS(NRM)
      COMMON /TRSRC4/ PNB(NRM),SNB(NRM),PBIN(NRM),PBCL(NRM,NSTM)
      COMMON /TRSRC5/ PNF(NRM),SNF(NRM),PFIN(NRM),PFCL(NRM,NSTM)
      COMMON /TRSRC6/ POH(NRM),PRF(NRM,NSTM),SPE(NRM,NSTM)
      COMMON /TRSRC7/ PRL(NRM),PCX(NRM),PIE(NRM),SIE(NRM),SCX(NRM),
     &                TSIE(NRM),TSCX(NRM)
      COMMON /TRSRC8/ PRFV(NRM,NSTM,3),AJRFV(NRM,3)
C
C     ****** COEFFICIENT VARIABLES ******
C
      COMMON /TRCEF1/ ETA(NRM)
      COMMON /TRCEF2/ AK(NRM,NSTM),AVK(NRM,NSTM),AD(NRM,NSTM),
     &                AV(NRM,NSTM)
      COMMON /TRCEF3/ AKNC(NRM,NSTM),AKDW(NRM,NSTM),MDLKNC
      COMMON /TRCEF4/ ADNC(NRM,NSTM),ADDW(NRM,NSTM)
      COMMON /TRCEF5/ AVNC(NRM,NSTM),AVDW(NRM,NSTM)
      COMMON /TRCEF6/ TAUB(NRM),TAUF(NRM)
      COMMON /TRCEF7/ VGR1(NRM,4),VGR2(NRM,4),VGR3(NRM,4),VGR4(NRM,4)
      COMMON /TRCEF8/ KGR1,KGR2,KGR3,KGR4
C
C     ****** GLOBAL VARIABLES ******
C
      COMMON /TRGLB1/ ANS0(NSTM),TS0(NSTM),ANSAV(NSTM),
     &                TSAV(NSTM),WST(NSTM)
      COMMON /TRGLB2/ ANF0(NFM),TF0(NFM),ANFAV(NFM),TFAV(NFM),WFT(NFM)
      COMMON /TRGLB3/ WBULKT,WTAILT,WPT
      COMMON /TRGLB4/ AJT,AJOHT,AJNBT,AJRFT,AJBST,AJRFVT(3)
      COMMON /TRGLB5/ PINT,POHT,PNBT,PRFT(NSTM),PNFT,PRFVT(NSTM,3)
      COMMON /TRGLB6/ PBINT,PBCLT(NSTM),PFINT,PFCLT(NSTM)
      COMMON /TRGLB7/ POUT,PCXT,PIET,PRLT,PLT(NSTM)
      COMMON /TRGLB8/ SINT,SIET,SNBT,SNFT,SPET(NSTM),SOUT,SLT(NSTM)
      COMMON /TRGLB9/ VLOOP,ALI,RQ1,RPE,Q0,ZEFF0,QF
      COMMON /TRGLBA/ WPDOT,TAUE1,TAUE2,TAUEP
      COMMON /TRGLBB/ BETAP0,BETAPA,BETA0,BETAA,BETAQ0,BETAN
C
C     ****** GRAPHIC DATA VARIABLES ******
C
      COMMON /TRPLR1/ GVR(NRMP,NGM,NCGM),GRM(NRM),GRG(NRMP),GTR(NGM),NGR
      COMMON /TRPLT1/ GVT(NTM,NCTM),GT(NTM),GTS(NTM),NGT,NGST
      COMMON /TRPLG1/ GYR(NRMP,8),GYT(NTM,8)
      COMMON /TRPLD1/ G3D(NRM,NTM,NCTM)
      COMMON /TRPLC1/ GJB(NRMP,6),GET(NRMP,4),GAD(NRMP,4),GAK(NRMP,4)
C
C     ****** EQUILIBRIUM INTERFACE VARIABLES ******
C
      COMMON /TREQV1/ RHOTR(NRM)
      COMMON /TREQV2/ PRHO(NRM),HJRHO(NRM),VTRHO(NRM),TRHO(NRM)
      COMMON /TREQV3/ QRHO(NRM),TTRHO(NRM),DVRHO(NRM),DSRHO(NRM)
      COMMON /TREQV4/ ABRHO(NRM),ARRHO(NRM)
      COMMON /TREQV5/ AR1RHO(NRM),AR2RHO(NRM)
      COMMON /TREQV6/ EPSRHO(NRM),BPRHO(NRM)
      COMMON /TREQV7/ EKAPPA(NRM),RMJRHO(NRM),RMNRHO(NRM)
      COMMON /TREQV8/ BPSEQ,RIPEQ,FACTJ,FACTQ(NRM)
C
C     ****** LOG FILE NAME ******
C
      COMMON /TRLOG1/ KFNLOG
C
C     ****** ADDITIONAL PART BY HONDA MITSURU ******
C
      COMMON /TRADD1/ VV(NVM,NVM,4,3)
      COMMON /TRADD2/ DD(NVM,NVM,4,3)
      COMMON /TRADD3/ VI(NVM,NVM,2,3),DI(NVM,NVM,2,3)
      COMMON /TRADD4/ FA(3,3),FB(3,3),FC(3,3)
      COMMON /TRADD5/ RTM(NSTM),AMZ(NSTM)
      COMMON /TRADD6/ PEX(NRM,NSTM),PEXT(NSTM),SEX(NRM,NSTM)
      COMMON /TRADD7/ NTEQIT
      COMMON /TRADD8/ PNSSA(NSTM),PNSA(NSTM),PTSA(NSTM)
      COMMON /TRADD9/ VOID
C
C     ****** MODEL SELECTION VARIABLES ******
C
      COMMON /TRMDS1/ NSS(NEQM),NSV(NEQM),NNS(NEQM),NST(NEQM)
      COMMON /TRMDS2/ NEQMAX
      COMMON /TRMDS3/ MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE
      COMMON /TRMDS4/ NSMAX,NSZMAX,NSNMAX,NSCMAX,NSTMAX
C
C     ****** FOR NCLASS ******
C
      COMMON /TRNCL1/ MDNCLS
      COMMON /TRNCL2/ AJBSNC(NRM),        ETANC(NRM),        AJEXNC(NRM)
      COMMON /TRNCL3/ CJBSP(NRM,NSM),     CJBST(NRM,NSM)
      COMMON /TRNCL4/ AKNCP(NRM,NSM,NSM), AKNCT(NRM,NSM,NSM)
      COMMON /TRNCL5/ ADNCP(NRM,NSM,NSM), ADNCT(NRM,NSM,NSM)
      COMMON /TRNCL6/ RGFLS(NRM,5,NSM),   RQFLS(NRM,5,NSM)
      COMMON /TRNCL7/ AVKNCS(NRM,NSM),    AVNCS(NRM,NSM)
      COMMON /TRNCL8/ AKNCLA(NRM,NSM,NSM),AKNCLB(NRM,NSM,NSM)
      COMMON /TRNCL9/ ADNCLA(NRM,NSM,NSM,2),ADNCLB(NRM,NSM,NSM,2)
      COMMON /TRNCLA/ AR1RHOG(NRM),AR2RHOG(NRM)
      COMMON /TRNCLB/ RMJRHOG(NRM),RMNRHOG(NRM)
C
C     ****** STORED VARIABLES FOR UFILE ******
C
      COMMON /TRSVU1/ QPU(NRMP,NTUM),     RNU(NRMP,NSM,NTUM),
     &                RTU(NRMP,NSM,NTUM), AJU(NRMP,NTUM)
      COMMON /TRSVU2/ PNBU(NRMP,NSM,NTUM),PICU(NRMP,NSM,NTUM),
     &                PRLU(NRMP,NTUM),    PECU(NRMP,NTUM)
      COMMON /TRSVU3/ DVRHOU(NRMP,NTUM),  DSRHOU(NRMP,NTUM)
      COMMON /TRSVU4/ RMJRHOU(NRMP,NTUM), RMNRHOU(NRMP,NTUM)
      COMMON /TRSVU5/ ARRHOU(NRMP,NTUM),  AR1RHOU(NRMP,NTUM),
     &                                    AR2RHOU(NRMP,NTUM)
      COMMON /TRSVU6/ ABRHOU(NRMP,NTUM),  TTRHOU(NRMP,NTUM)
      COMMON /TRSVU7/ PTSU(NSTM,NTUM),    PNSU(NSTM,NTUM)
      COMMON /TRSVU7/ PTSUA(NSTM,NTUM),   PNSUA(NSTM,NTUM)
      COMMON /TRSVU8/ SNBU(NRMP,NSM,NTUM)
      COMMON /TRSVU9/ RRU(NTUM),RAU(NTUM)
      COMMON /TRSVUA/ RIPU(NTUM),BBU(NTUM),RKAPU(NTUM)
      COMMON /TRSVUB/ AJBSU(NRMP,NTUM),   ZEFFU(NRMP,NTUM)
C
C     ****** UFILE CONTROL ******
C
      COMMON /TRUFL1/ KUFDEV,KUFDCG
      COMMON /TRUFL2/ MDLUF
      COMMON /TRUFL3/ MODEP,MDNI,MDCURT,MDNM1
      COMMON /TRUFL4/ TIME_INT
