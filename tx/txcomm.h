C     $Id$
C
      IMPLICIT REAL*8(A-F, H, O-Z)
      PARAMETER (NRM = 100, NQM = 20, NCM = 17, 
     &           NGRM = 20, NGTM = 1000,  NGVM = 1000, 
     &           NGYRM = 36, NGYTM = 39,  NGYVM = 43, 
     &           NGPRM = 13, NGPTM = 13,  NGPVM = 7)
      PARAMETER (NSM=2, NFM=2)
C
      PARAMETER (LQm1=1, LQm2=2, LQm3=3, LQm4=4, LQm5=5)
      PARAMETER (LQe1=6, LQe2=7, LQe3=8, LQe4=9, LQe5=10)
      PARAMETER (LQi1=11,LQi2=12,LQi3=13,LQi4=14,LQi5=15)
      PARAMETER (LQb1=16,        LQb3=17,LQb4=18)
      PARAMETER (LQn1=19,LQn2=20)
C
      COMMON /TXABC1/ ALC(0:NCM,NQM,0:NRM),BLC(0:NCM,NQM,0:NRM)
      COMMON /TXABC2/ CLC(0:NCM,NQM,0:NRM),PLC(NCM,NQM,0:NRM)
      COMMON /TXABC3/ NLC(0:NCM,NQM,0:NRM),NLCMAX(NQM)
      COMMON /TXCAL1/ DR,TIME,TMAX,NT,NQMAX,IERR
      COMMON /TXCAL2/ AMi,AMb,Vb,Bthb
      COMMON /TXCAL3/ rIP
      COMMON /TXCEF1/ rNuION(0:NRM),D01(0:NRM),D02(0:NRM)
      COMMON /TXCEF2/ rNu0e(0:NRM),rNu0i(0:NRM)
      COMMON /TXCEF3/ rNuL(0:NRM),rNuiCX(0:NRM)
      COMMON /TXCEF4/ rNuei(0:NRM),rNuii(0:NRM),rNuTei(0:NRM)
      COMMON /TXCEF5/ rNube(0:NRM),rNubi(0:NRM)
      COMMON /TXCEF6/ rNueNC(0:NRM),rNuiNC(0:NRM)
      COMMON /TXCEF7/ FWthe(0:NRM),FWthi(0:NRM),WPM(0:NRM)
      COMMON /TXCEF9/ rMue(0:NRM),rMui(0:NRM)
      COMMON /TXCEFA/ Chie(0:NRM),Chii(0:NRM)
      COMMON /TXCEFB/ PNB(0:NRM),SNB(0:NRM),rNuB(0:NRM)
      COMMON /TXCEFC/ PRFe(0:NRM),PRFi(0:NRM),PBIN(0:NRM),PBCL(0:NRM)
      COMMON /TXCEFD/ De(0:NRM),Di(0:NRM),rG1h2(0:NRM),FCDBM(0:NRM)
      COMMON /TXCEFE/ SiLC(0:NRM),SiLCth(0:NRM),SiLCph(0:NRM)
      COMMON /TXCEFF/ S(0:NRM), Alpha(0:NRM), rKappa(0:NRM)
      COMMON /TXCEFG/ AJ(0:NRM),AJOH(0:NRM),POH(0:NRM)
      COMMON /TXCEFH/ AJRF(0:NRM),AJNB(0:NRM),AJBS(0:NRM)
      COMMON /TXCNS1/ AEE,AME,AMP,VC,PI,rMU0,EPS0,rKEV
      COMMON /TXGLB1/ ANS0(NSM),TS0(NSM),ANSAV(NSM),TSAV(NSM),WST(NSM)
      COMMON /TXGLB2/ ANF0(NFM),TF0(NFM),ANFAV(NFM),TFAV(NFM),WFT(NFM)
      COMMON /TXGLB3/ WBULKT,WTAILT,WPT
      COMMON /TXGLB4/ AJT,AJOHT,AJNBT,AJRFT,AJBST
      COMMON /TXGLB5/ PINT,POHT,PNBT,PRFT,PNFT
      COMMON /TXGLB6/ PBINT,PBCLT(NSM),PFINT,PFCLT(NSM)
      COMMON /TXGLB7/ POUT,PCXT,PIET,PRLT,PLT(NSM)
      COMMON /TXGLB8/ SINT,SIET,SNBT,SNFT,SPET(NSM),SOUT,SLT(NSM)
      COMMON /TXGLB9/ VLOOP,ALI,RQ1,RPE,Q0,ZEFF0,QF
      COMMON /TXGLBA/ WPDOT,TAUE1,TAUE2,TAUEP
      COMMON /TXGLBB/ BETAP0,BETAPA,BETA0,BETAA,BETAQ0
      COMMON /TXGLBB/ TPRE,WPPRE
      COMMON /TXMTX1/ BA(4*NQM-1,NQM*(NRM+1)),BX(NQM*(NRM+1))
      COMMON /TXPRM1/ RA,RB,RR,BB
      COMMON /TXPRM2/ PA,PZ,Zeff
      COMMON /TXPRM3/ PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ
      COMMON /TXPRM4/ De0,Di0,rMue0,rMui0,WPM0
      COMMON /TXPRM5/ Chie0,Chii0
      COMMON /TXPRM6/ FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
      COMMON /TXPRM7/ FSCX,FSLC,FSNC,FSLP,FSION,FSD0
      COMMON /TXPRM8/ rLn,rLT
      COMMON /TXPRM9/ Eb,RNB,PNBH,PNBCD,rNRF,RRF,PRFH
      COMMON /TXPRMA/ PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
      COMMON /TXPRMB/ DLT,DT,EPS,ICMAX
      COMMON /TXPRMC/ NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
      COMMON /TXPRMD/ DelR,DelN
      COMMON /TXPRME/ rG1
      COMMON /TXPRMF/ rIPs,rIPe
      COMMON /TXSTP1/ X(NQM,0:NRM),XN(NQM,0:NRM)
      COMMON /TXVAG1/ NGR,NGT,NGVV
      COMMON /TXVAG4/ GX(0:NRM,2)
      COMMON /TXVAG5/ GY(0:NRM,0:NGRM,NGYRM),GT(0:NGRM)
      COMMON /TXVAG6/ GTX(0:NGTM),GTY(0:NGTM,NGYTM)
      COMMON /TXVAG7/ MODEG,gDIV(NGYRM),MODEAV,MODEL
      COMMON /TXVAG8/ GVX(0:NGVM),GVY(0:NGVM,NGYVM)
      COMMON /TXVAR0/ R(0:NRM),RHI(0:NRM)
      COMMON /TXVAR1/ ErI(0:NRM),ErHI(0:NRM)
      COMMON /TXVAR2/ EthI(0:NRM),EthHI(0:NRM)
      COMMON /TXVAR3/ EphHI(0:NRM),EphI(0:NRM)
      COMMON /TXVAR4/ BthI(0:NRM),BthHI(0:NRM)
      COMMON /TXVAR5/ BphHI(0:NRM),BphI(0:NRM)
      COMMON /TXVAR6/ PNeHI(0:NRM),PNeI(0:NRM)
      COMMON /TXVAR7/ UerI(0:NRM),UerHI(0:NRM)
      COMMON /TXVAR8/ UethI(0:NRM),UethHI(0:NRM)
      COMMON /TXVAR9/ UephI(0:NRM),UephHI(0:NRM)
      COMMON /TXVARA/ PTeHI(0:NRM),PTeI(0:NRM)
      COMMON /TXVARB/ PNiHI(0:NRM),PNiI(0:NRM)
      COMMON /TXVARC/ UirI(0:NRM),UirHI(0:NRM)
      COMMON /TXVARD/ UithI(0:NRM),UithHI(0:NRM)
      COMMON /TXVARE/ UiphI(0:NRM),UiphHI(0:NRM)
      COMMON /TXVARF/ PTiHI(0:NRM),PTiI(0:NRM)
      COMMON /TXVAS0/ PNbHI(0:NRM),PNbI(0:NRM)
      COMMON /TXVAS1/ UbthI(0:NRM),UbthHI(0:NRM)
      COMMON /TXVAS2/ UbphI(0:NRM),UbphHI(0:NRM)
      COMMON /TXVAS3/ PN01HI(0:NRM),PN02HI(0:NRM)
      COMMON /TXVAS4/ Q(0:NRM),QHI(0:NRM)
      COMMON /TXVER1/ SLID
      CHARACTER SLID*20




