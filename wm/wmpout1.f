C     $Id$
C
C     ****** CALCULATE ABSORBED POWER ******
C
      SUBROUTINE WMPABS1
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RU(NSM)
      DIMENSION DS(NRM),DSS(NTHM,NPHM,NRM)
      DIMENSION CPF1(MDM,NDM),CPF2(MDM,NDM)
      DIMENSION CDV(3,3,3),CDW(3,3,3)
      DIMENSION CFA(NRM*NSM*MDM*MDM*NDM*NDM)
      DIMENSION CFB(NRM*NSM*MDM*MDM*NDM*NDM)
C
      NM=NRM*NSM*MDM*MDM*NDM*NDM
      CW=2*PI*CRF*1.D6
      DTH=2.D0*PI/DBLE(NTHMAX)
      DPH=2.D0*PI/DBLE(NPHMAX)
C
      IF(MYRANK.EQ.0) THEN
         NRS=NBST
      ELSE
         NRS=NBST-1
      ENDIF
      IF(MYRANK.EQ.NPROCS-1) THEN
         NRE=NBED+1
      ELSE
         NRE=NBED
      ENDIF
C         
      DO NR=1,NRMAX+1
      DO NS=1,NSMAX
      DO NKX=1,NDSIZ
      DO KDX=1,KDSIZ
      DO MLX=1,MDSIZ
      DO LDX=1,LDSIZ
         CPABS(LDX,MLX,KDX,NKX,NS,NR)=(0.D0,0.D0)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      DO NS=1,NSMAX
C
         CALL WMSETF(NRS,NS)
C
         DO NDX=1,NDSIZ
         DO KDX=1,KDSIZ
         DO MDX=1,MDSIZ
         DO LDX=1,LDSIZ
         DO J=1,3
         DO I=1,3
            CGD(I,J,LDX,MDX,KDX,NDX,2)=CGD(I,J,LDX,MDX,KDX,NDX,3)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
C
         DO NR=NRS,NBED
C
            CALL WMSETF(NR+1,NS)
C
            DO NDX=1,NDSIZ
            DO KDX=1,KDSIZ
            DO MDX=1,MDSIZ
            DO LDX=1,LDSIZ
            DO J=1,3
            DO I=1,3
               CGD(I,J,LDX,MDX,KDX,NDX,1)=CGD(I,J,LDX,MDX,KDX,NDX,2)
               CGD(I,J,LDX,MDX,KDX,NDX,2)=CGD(I,J,LDX,MDX,KDX,NDX,3)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
C
            IF(NR.EQ.1) THEN
               XRHOM=XRHO(2)/1.D6
            ELSE
               XRHOM =XRHO(NR)
            ENDIF
            XRHOC =       XRHO(NR+1)
            IF(NR.EQ.NRMAX) THEN
               XRHOP=XRHO(NR+1)
            ELSE
               XRHOP=XRHO(NR+2)
            ENDIF
            XRHOMH=0.5D0*(XRHOC+XRHOM)
            XRHOPH=0.5D0*(XRHOC+XRHOP)
C
            DRHOM =XRHO(NR+1)-XRHO(NR)
            IF(NR.EQ.NRMAX) THEN
               NRPP=NR+1
               DRHOP=XRHO(NR+1)-XRHO(NR)
            ELSE
               NRPP=NR+2
               DRHOP=XRHO(NR+2)-XRHO(NR+1)
            ENDIF
C
            IF(MODELG.EQ.3) THEN
               QPMH=0.5D0*(QPS(NR)+QPS(NR+1))
               QPC=QPS(NR+1)
               DPSIPDRHOMH=2.D0*PSITA*XRHOMH/QPMH
               DPSIPDRHOC =2.D0*PSITA*XRHOC /QPC
            ELSE
               DPSIPDRHOMH=2.D0*PSIPA*XRHOMH
               DPSIPDRHOC =2.D0*PSIPA*XRHOC
            ENDIF
C
            DRHOPM= 0.5D0*(DRHOM+DRHOP)
C
            FMHM=0.5D0
            FMHC=0.5D0
C
            FCMH  = DRHOM/(2.D0*DRHOPM)
            FCPH  = DRHOP/(2.D0*DRHOPM)
C
C        ND : (n - n0) / Np                NDMIN -> NDMAX
C        NC : (n' - n0) / Np               NDMIN -> NDMAX
C        NK : (n + k - n0) / Np            NDMIN -> NDMAX
C
C        KD : (NK - ND) / Np : k           -KDISZ -> KDSIZ
C        KK : (NK - NC) / Np : k + n - n'  -KDISZ -> KDSIZ
C
C        MD : m - m0                  MDMIN -> MDMAX
C        MC : m' - m0                 MDMIN -> MDMAX
C        ML : m + l - m0              MDMIN -> MDMAX
C        
C        LD : ML - MD : l             -LDISZ -> LDSIZ
C        LL : ML - MC : l + m - m'    -LDISZ -> LDSIZ
C
            DO ND=NDMIN,NDMAX
               NDX=ND-NDMIN+1
            DO NC=NDMIN,NDMAX
               NCX=NC-NDMIN+1
            DO NK=NDMIN,NDMAX
               KDX=MOD(NK-ND-KDMIN+2*KDSIZ,KDSIZ)+1
               KD=KDX+KDMIN-1
               KKX=MOD(NK-NC-KDMIN+2*KDSIZ,KDSIZ)+1
               KK=KKX+KDMIN-1
C
            DO MD=MDMIN,MDMAX
               MDX=MD-MDMIN+1
            DO MC=MDMIN,MDMAX
               MCX=MC-MDMIN+1
            DO ML=MDMIN,MDMAX
               LDX=MOD(ML-MD-LDMIN+2*LDSIZ,LDSIZ)+1
               LD=LDX+LDMIN-1
               LLX=MOD(ML-MC-LDMIN+2*LDSIZ,LDSIZ)+1
               LL=LLX+LDMIN-1
C
               DO K=1,2
               DO J=1,3
               DO I=1,3
                  CDV(I,J,K)=CGD(I,J,LDX,MDX,KDX,NDX,K)
                  CDW(I,J,K)=CGD(I,J,LDX,MDX,KDX,NDX,K)
               ENDDO
               ENDDO
               ENDDO
C
               FACT1M=XRHOMH/XRHOM
               FACT1C=XRHOMH/XRHOC
C               FACT1P=XRHOPH/XRHOC
C               FACT2M=XRHOMH/XRHOM
               FACT2C=XRHOMH/XRHOC
               FACT2P=XRHOPH/XRHOC
C               FACT3M=XRHOMH/XRHOM
               FACT3C=XRHOMH/XRHOC
               FACT3P=XRHOPH/XRHOC
C
               CDV11M=CHERMIT(CDV(1,1,1),CDW(1,1,1))
               CDV11C=CHERMIT(CDV(1,1,2),CDW(1,1,2))
               CDV12M=CHERMIT(CDV(1,2,1),CDW(2,1,1))
               CDV12C=CHERMIT(CDV(1,2,2),CDW(2,1,2))
               CDV13M=CHERMIT(CDV(1,3,1),CDW(3,1,1))
               CDV13C=CHERMIT(CDV(1,3,2),CDW(3,1,2))
               CDV21C=CHERMIT(CDV(2,1,2),CDW(1,2,2))
               CDV22C=CHERMIT(CDV(2,2,2),CDW(2,2,2))
               CDV23C=CHERMIT(CDV(2,3,2),CDW(3,2,2))
               CDV31C=CHERMIT(CDV(3,1,2),CDW(1,3,2))
               CDV32C=CHERMIT(CDV(3,2,2),CDW(2,3,2))
               CDV33C=CHERMIT(CDV(3,3,2),CDW(3,3,2))
C
C               CDV11M=CDV(1,1,1)
C               CDV11C=CDV(1,1,2)
C               CDV12M=CDV(1,2,1)
C               CDV12C=CDV(1,2,2)
C               CDV13M=CDV(1,3,1)
C               CDV13C=CDV(1,3,2)
C               CDV21C=CDV(2,1,2)
C               CDV22C=CDV(2,2,2)
C               CDV23C=CDV(2,3,2)
C               CDV31C=CDV(3,1,2)
C               CDV32C=CDV(3,2,2)
C               CDV33C=CDV(3,3,2)
C
C     --- R COMPONENT OF MAXWELL EQUATION ---
C
               CEMM11=0.5D0*CDV11M*FACT1M /XRHOMH/XRHOMH
               CEMC11=0.5D0*CDV11C*FACT1C /XRHOMH/XRHOMH
               CEMM12= FMHM*CDV12M        *XRHOM /XRHOMH
               CEMC12= FMHC*CDV12C        *XRHOC /XRHOMH
               CEMM13= FMHM*CDV13M               /XRHOMH
               CEMC13= FMHC*CDV13C               /XRHOMH
C
C     --- THETA COMPONENT OF MAXWELL EQUATION ---
C
               CEMC21= FCMH*CDV21C*FACT2C /XRHOMH*XRHOC
               CEMP21= FCPH*CDV21C*FACT2P /XRHOPH*XRHOC
               CEMC22=      CDV22C        *XRHOC *XRHOC
               CEMC23=      CDV23C               *XRHOC
C
C     --- PHI COMPONENT OF MAXWELL EQUATION ---
C
               CEMC31= FCMH*CDV31C*FACT3C /XRHOMH
               CEMP31= FCPH*CDV31C*FACT3P /XRHOPH
               CEMC32=      CDV32C        *XRHOC
               CEMC33=      CDV33C
C
               CCE1=CEFLDK(1,MCX,NCX,NR+1)
               CCE2=CEFLDK(2,MCX,NCX,NR+1)
               CCE3=CEFLDK(3,MCX,NCX,NR+1)
C
               CJM1=CEMM11*CEFLDK(1,MDX,NDX,NR+1)
     &             +CEMM12*CEFLDK(2,MDX,NDX,NR)
     &             +CEMM13*CEFLDK(3,MDX,NDX,NR)
               CJC1=CEMC11*CEFLDK(1,MDX,NDX,NR+1)
     &             +CEMC12*CEFLDK(2,MDX,NDX,NR+1)
     &             +CEMC13*CEFLDK(3,MDX,NDX,NR+1)
               CJC2=CEMC21*CEFLDK(1,MDX,NDX,NR+1)
     &             +CEMC22*CEFLDK(2,MDX,NDX,NR+1)
     &             +CEMC23*CEFLDK(3,MDX,NDX,NR+1)
               CJC3=CEMC31*CEFLDK(1,MDX,NDX,NR+1)
     &             +CEMC32*CEFLDK(2,MDX,NDX,NR+1)
     &             +CEMC33*CEFLDK(3,MDX,NDX,NR+1)
               CJP2=CEMP21*CEFLDK(1,MDX,NDX,NRPP)
               CJP3=CEMP31*CEFLDK(1,MDX,NDX,NRPP)
C
               CPM1=DCONJG(CCE1)*CJM1
               CPC1=DCONJG(CCE1)*CJC1
               CPC2=DCONJG(CCE2)*CJC2
               CPC3=DCONJG(CCE3)*CJC3
               CPP2=DCONJG(CCE2)*CJP2
               CPP3=DCONJG(CCE3)*CJP3
C     
               CPABSM=CI*EPS0*(VC*VC/CW)*CPM1
     &               *DPSIPDRHOMH*DRHOM
C
               CPABSC=CI*EPS0*(VC*VC/CW)*CPC1
     &               *DPSIPDRHOMH*DRHOM
     &               +CI*EPS0*(VC*VC/CW)*(CPC2+CPC3+CPP2+CPP3)
     &               *DPSIPDRHOC*DRHOPM
C
               CPABS(LLX,MDX,KKX,NDX,NS,NR)
     &        =CPABS(LLX,MDX,KKX,NDX,NS,NR)
     &        +0.5D0*CPABSM
C
               CPABS(LLX,MDX,KKX,NDX,NS,NR+1)
     &        =CPABS(LLX,MDX,KKX,NDX,NS,NR+1)
     &        +0.5D0*CPABSC
C
C               KKY=MOD(NC-NK-KDMIN+2*KDSIZ,KDSIZ)+1
C               LLY=MOD(MC-ML-LDMIN+2*LDSIZ,LDSIZ)+1
C
C               CPM1=CCE1*DCONJG(CJM1)
C               CPC1=CCE1*DCONJG(CJC1)
C               CPC2=CCE2*DCONJG(CJC2)
C               CPC3=CCE3*DCONJG(CJC3)
C               CPP2=CCE2*DCONJG(CJP2)
C               CPP3=CCE3*DCONJG(CJP3)
C     
C               CPABSM=CI*EPS0*(VC*VC/CW)*CPM1
C     &                *DPSIPDRHOMH*DRHOM
C
C               CPABSC=CI*EPS0*(VC*VC/CW)*CPC1
C     &                *DPSIPDRHOMH*DRHOM
C     &               +CI*EPS0*(VC*VC/CW)*(CPC2+CPC3+CPP2+CPP3)
C     &                *DPSIPDRHOC*DRHOPM
C
C               CPABS(LLY,MDX,KKY,NDX,NS,NR)
C     &        =CPABS(LLY,MDX,KKY,NDX,NS,NR)
C     &        -0.5D0*CPABSM
C
C               CPABS(LLY,MDX,KKY,NDX,NS,NR+1)
C     &        =CPABS(LLY,MDX,KKY,NDX,NS,NR+1)
C     &        -0.5D0*CPABSC
C
            ENDDO
            ENDDO
            ENDDO
C
            ENDDO
            ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
      NRS=NBST
      IF(MYRANK.EQ.NPROCS-1) THEN
         NRE=NBED+1
      ELSE
         NRE=NBED
      ENDIF
C
      MN=0
      DO NR=NRS,NRE
      DO NS=1,NSMAX
      DO NDX=1,NDSIZ
      DO KKX=1,KDSIZ
      DO MDX=1,MDSIZ
      DO LLX=1,LDSIZ
         MN=MN+1
         CFB(MN)=CPABS(LLX,MDX,KKX,NDX,NS,NR)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      CALL MPGTCV(CFB,MN,CFA,NVTOT,NM)
C
      IF(MYRANK.EQ.0) THEN
         MN=0
         DO NR=1,NRMAX
         DO NS=1,NSMAX
         DO NDX=1,NDSIZ
         DO KKX=1,KDSIZ
         DO MDX=1,MDSIZ
         DO LLX=1,LDSIZ
            MN=MN+1
            CPABS(LLX,MDX,KKX,NDX,NS,NR)=CFA(MN)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
C
C     +++++ CALCULATE PABS IN MODE NUMBER SPACE +++++
C
         DO NR=1,NRMAX+1
         DO NS=1,NSMAX
         DO NDX=1,NDSIZ
            KK=0
            KKX=KK-KDMIN+1
         DO MDX=1,MDSIZ
            LL=0
            LLX=LL-LDMIN+1
            PABSK(MDX,NDX,NR,NS)=DBLE(CPABS(LLX,MDX,KKX,NDX,NS,NR))
         ENDDO
         ENDDO
         ENDDO
         ENDDO
C
C     +++++ CALCULATE PABS IN REAL SPACE +++++
C
         DO NS=1,NSMAX
         DO NR=1,NRMAX+1
            DO NPH=1,NPHMAX
            DO NTH=1,NTHMAX
               PABS(NTH,NPH,NR,NS)=0.D0
            ENDDO
            ENDDO
            DO NDX=1,NDSIZ
            DO MDX=1,MDSIZ
               DO KK=KDMIN,KDMAX
                  KKX=KK-KDMIN+1
               DO LL=LDMIN,LDMAX
                  LLX=LL-LDMIN+1
                  CPF1(LLX,KKX)=CPABS(LLX,MDX,KKX,NDX,NS,NR)
               ENDDO
               ENDDO
               CALL WMSUBE(CPF1,CPF2)
               DO NPH=1,NPHMAX
               DO NTH=1,NTHMAX
                  PABS(NTH,NPH,NR,NS)=PABS(NTH,NPH,NR,NS)
     &                               +DBLE(CPF2(NTH,NPH))
               ENDDO
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
C
C     +++++ CALCULATE DRIVEN CURRENT IN REAL SPACE +++++
C
      NS=1
      DO NR=1,NRMAX+1
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            PCUR(NTH,NPH,NR)=0.D0
         ENDDO
         ENDDO
         CALL WMCDEN(NR,RN,RTPR,RTPP,RU)
         VTE=SQRT(RTPR(1)*AEE*1.D3/(PA(1)*AMP))
         WW=DBLE(CW)
         IF(RN(1).LE.0.D0) THEN
            RLNLMD=15.D0
         ELSE
            RT=(RTPR(1)+2*RTPP(1))/3.D0
            RLNLMD=16.1D0 - 1.15D0*LOG10(RN(1))
     &           + 2.30D0*LOG10(RT)
         ENDIF
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
            NN=NPH0+ND
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            MM=NTH0+MD
            DO KKX=1,KDSIZ
            DO LLX=1,LDSIZ
               CPF1(LLX,KKX)=CPABS(LLX,MDX,KKX,NDX,NS,NR)
            ENDDO
            ENDDO
            CALL WMSUBE(CPF1,CPF2)
            DO NPH=1,NPHMAX
            DO NTH=1,NTHMAX
               CALL WMCMAG(NR,NTH,NPH,BABS,BSUPTH,BSUPPH)
               RKPR=MM*BSUPTH/BABS+NN*BSUPPH/BABS
              IF(ABS(RKPR).LT.1.D-8) RKPR=1.D-8
              IF(ABS(WW/RKPR).LT.VC) THEN
                 W=WW/(RKPR*VTE)
                 XL=(RPST(NTH,NPH,NR)-RR  )/RR
                 YL=(ZPST(NTH,NPH,NR)-0.D0)/RR
                 EFCD=W1CDEF(ABS(W),ZEFF,XL,YL,1)
                 IF(W.LT.0.D0) EFCD=-EFCD
                 IF (RN(1).GT.0.D0) THEN
                    PCUR(NTH,NPH,NR)=PCUR(NTH,NPH,NR)
     &                   +0.384D0*RTPR(1)/(AEE*1.D3)*EFCD
     &                   /((RN(1)/1.D20)*RLNLMD)*DBLE(CPF2(NTH,NPH))
     &                   /(2.D0*PI*RPST(NTH,NPH,NR))
                 END IF
              ENDIF
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         PCURR(NR)=0.D0
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            PCURR(NR)=PCURR(NR)+PCUR(NTH,NPH,NR)*DTH
         ENDDO
         ENDDO
      ENDDO
C
      DO NS=1,NSMAX
      DO NR=1,NRMAX
         PABSR(NR,NS)=0.D0
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            PABSR(NR,NS)=PABSR(NR,NS)+PABS(NTH,NPH,NR,NS)*DTH*DPH
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      PCURT=0.D0
      DO NR=1,NRMAX
         PCURT=PCURT+PCURR(NR)
      ENDDO
C
      DO NS=1,NSMAX
         PABST(NS)=0.D0
         DO NR=1,NRMAX
            PABST(NS)=PABST(NS)+PABSR(NR,NS)
         ENDDO
      ENDDO
C
      PABSTT=0.D0
      DO NS=1,NSMAX
         PABSTT=PABSTT+PABST(NS)
      ENDDO
C
      FACT=1.D0
C
      IF(PRFIN.GT.0.D0.AND.PABSTT.GT.0.D0) THEN
         FACT=PRFIN/PABSTT
      ENDIF
      FACTSQ=SQRT(FACT)
C
      NR=1
         DS(NR)=0.D0
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            DSS(NTH,NPH,NR)=0.D0
         ENDDO
         ENDDO
      DO NR=2,NRMAX
         DS(NR)=0.D0
         DRHO=0.5D0*(XRHO(NR+1)-XRHO(NR-1))
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            IF(MODELG.EQ.3) THEN
               DPSIPDRHO=2.D0*PSITA*XRHO(NR)/QPS(NR)
            ELSE
               DPSIPDRHO=2.D0*PSIPA*XRHO(NR)
            ENDIF
            DSSS=DPSIPDRHO*DRHO*RJ(NTH,NPH,NR)
            DSS(NTH,NPH,NR)=1.D0/DSSS
            DS(NR)=DS(NR)+DSSS*DTH*DPH
         ENDDO
         ENDDO
         DS(NR)=1.D0/DS(NR)
      ENDDO
C
      PABSTT=FACT*PABSTT
      DO NS=1,NSMAX
         PABST(NS)=FACT*PABST(NS)
         DO NR=1,NRMAX
            PABSR(NR,NS)=FACT*PABSR(NR,NS)*DS(NR)
            DO NPH=1,NPHMAX
            DO NTH=1,NTHMAX
               PABS(NTH,NPH,NR,NS)=FACT*PABS(NTH,NPH,NR,NS)
     &                            *DSS(NTH,NPH,NR)
            ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO NS=1,NSMAX
      DO NR=1,NRMAX
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            PABSK(MDX,NDX,NR,NS)=FACT*PABSK(MDX,NDX,NR,NS)*DS(NR)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      PCURT=FACT*PCURT
      DO NR=1,NRMAX
         PCURR(NR)=FACT*PCURR(NR)*DS(NR)
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            PCUR(NTH,NPH,NR)=FACT*PCUR(NTH,NPH,NR)*DSS(NTH,NPH,NR)
         ENDDO
         ENDDO
      ENDDO
      ENDIF
C
      DO NR=1,NRMAX+1
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            CEFLDK(1,MDX,NDX,NR)=FACTSQ*CEFLDK(1,MDX,NDX,NR)
            CEFLDK(2,MDX,NDX,NR)=FACTSQ*CEFLDK(2,MDX,NDX,NR)
            CEFLDK(3,MDX,NDX,NR)=FACTSQ*CEFLDK(3,MDX,NDX,NR)
            CBFLDK(1,MDX,NDX,NR)=FACTSQ*CBFLDK(1,MDX,NDX,NR)
            CBFLDK(2,MDX,NDX,NR)=FACTSQ*CBFLDK(2,MDX,NDX,NR)
            CBFLDK(3,MDX,NDX,NR)=FACTSQ*CBFLDK(3,MDX,NDX,NR)
         ENDDO
         ENDDO
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            CEFLD(1,NTH,NPH,NR)=FACTSQ*CEFLD(1,NTH,NPH,NR)
            CEFLD(2,NTH,NPH,NR)=FACTSQ*CEFLD(2,NTH,NPH,NR)
            CEFLD(3,NTH,NPH,NR)=FACTSQ*CEFLD(3,NTH,NPH,NR)
            CBFLD(1,NTH,NPH,NR)=FACTSQ*CBFLD(1,NTH,NPH,NR)
            CBFLD(2,NTH,NPH,NR)=FACTSQ*CBFLD(2,NTH,NPH,NR)
            CBFLD(3,NTH,NPH,NR)=FACTSQ*CBFLD(3,NTH,NPH,NR)
            CEN(1,NTH,NPH,NR)  =FACTSQ*CEN(1,NTH,NPH,NR)
            CEN(2,NTH,NPH,NR)  =FACTSQ*CEN(2,NTH,NPH,NR)
            CEN(3,NTH,NPH,NR)  =FACTSQ*CEN(3,NTH,NPH,NR)
            CEP(1,NTH,NPH,NR)  =FACTSQ*CEP(1,NTH,NPH,NR)
            CEP(2,NTH,NPH,NR)  =FACTSQ*CEP(2,NTH,NPH,NR)
            CEP(3,NTH,NPH,NR)  =FACTSQ*CEP(3,NTH,NPH,NR)
         ENDDO
         ENDDO
      ENDDO
      CALL MPSYNC
C
      RETURN
      END
