C  
C     ***********************************************************
C
C            Weiland Model
C
C     ***********************************************************
C
      SUBROUTINE TR_WEILAND
C
C***********************************************************************
C  <INPUT>
C     ENL    : 2 Lne/Lb (Lb=R)
C     EIL    : Lni/Lti
C     EEL    : Lne/Lte
C     TAUL   : Te/Ti
C     FLL    : Finite Larmor Radius parameter: (k_perp*rhos)**2=0.1
C     FTL    : trapped particle fraction
C     PDEVL  : maybe not used in this code
C     PSTL   : maybe not used in this code
C     BQL    : fraction of an impurity
C     EQL    : Lnq/Ltq (q represents an impurity)
C     ENQL   : 2 Lnq/Lb
C     ZL     : Z value for an impurity
C     BETAEL : electron beta
C     AZL    : impurity mass number
C     COLL   : factor multiplying collisionality (0 or 1)
C     ELL    : factor multiplying electromagnetic effects (0 or 1)
C     TEL    : electron temperature
C     TAUZL  : Te/Tq
C     RA     : minor radius (meter)
C     QL     : safety factor
C     SL     : magnetic shear
C     EPS    : inverse aspect ratio (a/R)
C     RNL    : electron density
C     RLIST  : controling printout in disp9t (0: off)
C     RNEQL  : number of equations (NDISP)
C     ALAL   : normalized pressure gradient
C     RKAP   : ellipticity
C     TRL    : not used in this code
C     RIWL   : controling printout in DIFFTD (0: off)
C     RISBL  : if ISB=1, we use the strong ballooning approximation
C              (GAV=1)
C     BB     : toroidal magnetic field
C     SEARCH : The way the code chooses between ion and electron
C              eigenvalues for the use in the eigenfunction (WZ).
C              1 : Only eigenvalues with negative real part are used.
C                 (Real(WZ)<0)
C              2 : Eigenvalues with positive real part are used unless
C                  there are eigenvalues with negative real part.
C                 (Real(WZ)>0)
C              3 : The fastest growing mode with positive real part is
C                  used for eigenvalues with positive real part 
C                  if there is such a root.
C     PMA    : mass number for an main ion
C     RGKL   : GRHO1/(a*GRHO2)
C     WEXBL  : ExB shearing rate
C     ROTL   : factor multiplying WEXBL (0 or 1)
C     NR     : radial grid number
C     IST    : 1   : iterations starting from an analytical
C                    approximation
C                    You should use only the first time step.
C              else: iterations starting from the value stored
C                    in WZJ(IK) which is the eigenvalue from the
C                    previous time step
C
C  <OUTPUT>
C     CHIL(5) : ion thermal transport coefficients for Ti, Te, Ne, Tq
C               and Nq equations
C     CHEL(5) : electron thermal transport coefficients
C     DL(5)   : ion particle transport coefficients
C     CHQL(5) : impurity thermal transport coefficients
C     DQL(5)  : impurity particle transport coefficients
C     SCHI    : effective ion thermal transport coefficient
C     SCHE    : effective electron thermal transport coefficient
C     SD      : effective ion particle transport coefficient
C     SCHQ    : effective impurity thermal transport coefficient
C     SDQ     : effective impurity particle transport coefficient
C
C***********************************************************************
C
      INCLUDE 'trcomm.inc'
      DIMENSION CHIL(5),CHEL(5),DL(5),CHQL(5),DQL(5)
C
      IF(NT.EQ.0) THEN
         IST=1
      ELSE
         IST=0
      ENDIF
      DO NR=1,NRMAX-1
         DRL   = RJCB(NR)/DR
         EPS   = EPSRHO(NR)
         SLNEL = 0.5D0*(RN(NR+1,1)+RN(NR,1))/((RN(NR+1,1)-RN(NR,1))*DRL)
         SLNIL = 0.5D0*(RN(NR+1,2)+RN(NR,2))/((RN(NR+1,2)-RN(NR,2))*DRL)
         SLNQL = 0.5D0*(RN(NR+1,3)+RN(NR,3))/((RN(NR+1,3)-RN(NR,3))*DRL)
         SLTEL = 0.5D0*(RT(NR+1,1)+RT(NR,1))/((RT(NR+1,1)-RT(NR,1))*DRL)
         SLTIL = 0.5D0*(RT(NR+1,2)+RT(NR,2))/((RT(NR+1,2)-RT(NR,2))*DRL)
         SLTQL = 0.5D0*(RT(NR+1,3)+RT(NR,3))/((RT(NR+1,3)-RT(NR,3))*DRL)
         SLBL  = RR
         ENL   =-2.D0*(SLNEL/SLBL )
         EIL   =       SLNIL/SLTIL
         EEL   =       SLNEL/SLTEL
         TAUL  = (RT(NR+1,1)+RT(NR,1))/(RT(NR+1,2)+RT(NR,2))
         FLL   = 1.D-1
         FTL   = 1.D0-(1.D0-EPS)**2.D0
     &          /(DSQRT(1.D0-EPS**2)*(1.D0+1.46D0*DSQRT(EPS)))
         PDEVL = 8.D0
         PSTL  = 1.D0
         BQL   = (RN(NR+1,3)+RN(NR,3))/(RN(NR+1,1)+RN(NR,1))
         EQL   =-      SLNQL/SLTQL
         ENQL  =-2.D0*(SLNQL/SLBL )
         ZL    = PZ(3)
         BETAEL= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
     &          *RKEV*1.D20/(BB**2/(2.D0*AMYU0))
         AZL   = PA(3)
         COLL  = 1.D0
         ELL   = 1.D0
         TEL   = 0.5D0*(RT(NR+1,1)+RT(NR,1))
         TAUZL = TEL/(0.5D0*(RT(NR+1,3)+RT(NR,3)))
         QL    = QP(NR)
         IF(NR.EQ.1) THEN
            DQ = (QP(NR+1)-QP(NR))/1.5D0*DRL
         ELSE
            DQ = (QP(NR+1)-QP(NR-1))/2.D0*DRL
         ENDIF
         SL    = RR*EPS*DQ/QL
         RNL   = 0.5D0*(RN(NR+1,1)+RN(NR,1))*1.D1
         RLIST = 1.D0
         RNEQL = 9.D0 
         RNTP=0.D0
         RNTM=0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+RN(NR+1,NS)*RT(NR+1,NS)
            RNTM=RNTM+RN(NR  ,NS)*RT(NR  ,NS)
         ENDDO
         RPP   = RNTP+RN(NR+1,1)*RT(NR+1,1)
         RPM   = RNTM+RN(NR  ,1)*RT(NR  ,1)
         DPP   = (RPP-RPM)*DRL
         DBDR  = DPP*1.D20*RKEV*RA/(BB**2/(2*AMYU0))
         ALAL  =-QL*QL*DBDR*RR/RA
         TRL   = 1.D0
         RIWL  = 2.D0
         RISBL = 2.D0
         SEARCH= 1.D0
         PMA   = PA(2)
         RGKL  = AR1RHOG(NR)/(RA*AR2RHOG(NR))
         WEXBL = 0.D0
         ROTL  = 0.D0
C
         CALL TR_WEILAND_BRIDGE
     &     (ENL,EIL,EEL,TAUL,FLL,FTL,PDEVL,PSTL,
     &      BQL,EQL,ENQL,ZL,BETAEL,AZL,COLL,ELL,TEL,TAUZL,RA,
     &      QL,SL,EPS,RNL,RLIST,RNEQL,ALAL,RKAP,TRL,RIWL,RISBL,BB,
     &      SEARCH,PMA,RGKL,WEXBL,ROTL,NR,IST,
     &      CHIL,CHEL,DL,CHQL,DQL,SCHI,SCHE,SD,SCHQ,SDQ)
C
C     MODE : mode selector
C            0    : using effective transport coefficients
C            else : using transport coefficients' vectors
         MODE=0
         IF(MODE.EQ.0) THEN
            AKDW(NR,1)=SCHE
            AKDW(NR,2)=SCHI
            AKDW(NR,3)=SCHQ
            AKDW(NR,4)=SCHQ
            ADDW(NR,1)=SD
            ADDW(NR,2)=SD
            ADDW(NR,3)=SDQ
            ADDW(NR,4)=SDQ
         ELSE
            DO NS=1,NSM
               DO NS1=1,NSM
                  DO NA=1,2
                     AKWLDW(NR,NS,NS1,NA)=0.D0
                     ADWLDW(NR,NS,NS1,NA)=0.D0
                  ENDDO
               ENDDO
            ENDDO
            DO NEQ=1,NEQMAX
               NSSN=NSS(NEQ)
               NSVN=NSV(NEQ)
               IF(NSVN.EQ.1) THEN
                  IF(NSSN.EQ.1) THEN
                     AKWLDW(NR,1,NSSN,NSVN)=CHEL(3)
                     ADWLDW(NR,1,NSSN,NSVN)=DL  (3)
                     AKWLDW(NR,2,NSSN,NSVN)=CHIL(3)
                     AKWLDW(NR,3,NSSN,NSVN)=CHQL(3)
                     ADWLDW(NR,3,NSSN,NSVN)=DQL (3)
                  ELSEIF(NSSN.EQ.3) THEN
                     AKWLDW(NR,1,NSSN,NSVN)=CHEL(5)
                     ADWLDW(NR,1,NSSN,NSVN)=DL  (5)
                     AKWLDW(NR,2,NSSN,NSVN)=CHIL(5)
                     AKWLDW(NR,3,NSSN,NSVN)=CHQL(5)
                     ADWLDW(NR,3,NSSN,NSVN)=DQL (5)
                  ENDIF
               ELSEIF(NSVN.EQ.2) THEN
                  IF(NSSN.EQ.1) THEN
                     AKWLDW(NR,1,NSSN,NSVN)=CHEL(2)
                     ADWLDW(NR,1,NSSN,NSVN)=DL  (2)
                     AKWLDW(NR,2,NSSN,NSVN)=CHIL(2)
                     AKWLDW(NR,3,NSSN,NSVN)=CHQL(2)
                     ADWLDW(NR,3,NSSN,NSVN)=DQL (2)
                  ELSEIF(NSSN.EQ.2) THEN
                     AKWLDW(NR,1,NSSN,NSVN)=CHEL(1)
                     ADWLDW(NR,1,NSSN,NSVN)=DL  (1)
                     AKWLDW(NR,2,NSSN,NSVN)=CHIL(1)
                     AKWLDW(NR,3,NSSN,NSVN)=CHQL(1)
                     ADWLDW(NR,3,NSSN,NSVN)=DQL (1)
                  ELSEIF(NSSN.EQ.3) THEN
                     AKWLDW(NR,1,NSSN,NSVN)=CHEL(4)
                     ADWLDW(NR,1,NSSN,NSVN)=DL  (4)
                     AKWLDW(NR,2,NSSN,NSVN)=CHIL(4)
                     AKWLDW(NR,3,NSSN,NSVN)=CHQL(4)
                     ADWLDW(NR,3,NSSN,NSVN)=DQL (4)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO
C     
      NR=NRMAX
         DRL   = RJCB(NR)/DR
         EPS   = EPSRHO(NR)
         SLNEL = PNSS(1)/(2.D0*(PNSS(1)-RN(NR,1))*DRL)
         SLNIL = PNSS(2)/(2.D0*(PNSS(2)-RN(NR,2))*DRL)
         SLNQL = PNSS(3)/(2.D0*(PNSS(3)-RN(NR,3))*DRL)
         SLTEL = PTS (1)/(2.D0*(PTS (1)-RT(NR,1))*DRL)
         SLTIL = PTS (2)/(2.D0*(PTS (2)-RT(NR,2))*DRL)
         SLTQL = PTS (3)/(2.D0*(PTS (3)-RT(NR,3))*DRL)
         SLBL  = RR
         ENL   =-2.D0*(SLNEL/SLBL )
         EIL   =       SLNIL/SLTIL
         EEL   =       SLNEL/SLTEL
         TAUL  = PTS(1)/PTS(2)
         FLL   = 1.D-1
         FTL   = 1.D0-(1.D0-EPS)**2.D0
     &          /(DSQRT(1.D0-EPS**2)*(1.D0+1.46D0*DSQRT(EPS)))
         PDEVL = 8.D0
         PSTL  = 1.D0
         BQL   = PNSS(3)/PNSS(1)
         EQL   =-      SLNQL/SLTQL
         ENQL  =-2.D0*(SLNQL/SLBL )
         ZL    = PZ(3)
         BETAEL= PNSS(1)*PTS(1)*RKEV*1.D20/(BB**2/(2.D0*AMYU0))
         AZL   = PA(3)
         COLL  = 1.D0
         ELL   = 1.D0
         TEL   = PTS(1)
         TAUZL = TEL/PTS(3)
         QL    = QP(NR)
         DQ    = (QP(NR)-QP(NR-1))*DRL
         SL    = RR*EPS*DQ/QL
         RNL   = PNSS(1)*1.D1
         RLIST = 1.D0
         RNEQL = 9.D0
         RNTP=0.D0
         RNTM=0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+PNSS(NS)*PTS(NS)
            RNTM=RNTM+RN(NR-1,NS)*RT(NR-1,NS)
     &               +RN(NR  ,NS)*RT(NR  ,NS)
         ENDDO
         RPP   = RNTP+PNSS(1)*PTS(1)
         RPM   = RNTM+RN(NR-1,1)*RT(NR-1,1)
     &               +RN(NR  ,1)*RT(NR  ,1)
         DPP   = (RPP-RPM)*DRL
         DBDR  = DPP*1.D20*RKEV*RA/(BB**2/(2*AMYU0))
         ALAL  =-QL*QL*DBDR*RR/RA
         TRL   = 1.D0
         RIWL  = 2.D0
         RISBL = 2.D0
         SEARCH= 2.D0
         PMA   = PA(2)
         RGKL  = AR1RHOG(NR)/(RA*AR2RHOG(NR))
         WEXBL = 0.D0
         ROTL  = 0.D0
C
         CALL TR_WEILAND_BRIDGE
     &     (ENL,EIL,EEL,TAUL,FLL,FTL,PDEVL,PSTL,
     &      BQL,EQL,ENQL,ZL,BETAEL,AZL,COLL,ELL,TEL,TAUZL,RA,
     &      QL,SL,EPS,RNL,RLIST,RNEQL,ALAL,RKAP,TRL,RIWL,RISBL,BB,
     &      SEARCH,PMA,RGKL,WEXBL,ROTL,NR,IST,
     &      CHIL,CHEL,DL,CHQL,DQL,SCHI,SCHE,SD,SCHQ,SDQ)
C
         IF(MODE.EQ.0) THEN
            AKDW(NR,1)=SCHE
            AKDW(NR,2)=SCHI
            AKDW(NR,3)=SCHQ
            AKDW(NR,4)=SCHQ
            ADDW(NR,1)=SD
            ADDW(NR,2)=SD
            ADDW(NR,3)=SDQ
            ADDW(NR,4)=SDQ
         ELSE
            DO NS=1,NSM
               DO NS1=1,NSM
                  DO NA=1,2
                     AKWLDW(NR,NS,NS1,NA)=0.D0
                     ADWLDW(NR,NS,NS1,NA)=0.D0
                  ENDDO
               ENDDO
            ENDDO
            DO NEQ=1,NEQMAX
               NSSN=NSS(NEQ)
               NSVN=NSV(NEQ)
               IF(NSVN.EQ.1) THEN
                  IF(NSSN.EQ.1) THEN
                     AKWLDW(NR,1,NSSN,NSVN)=CHEL(3)
                     ADWLDW(NR,1,NSSN,NSVN)=DL  (3)
                     AKWLDW(NR,2,NSSN,NSVN)=CHIL(3)
                     AKWLDW(NR,3,NSSN,NSVN)=CHQL(3)
                     ADWLDW(NR,3,NSSN,NSVN)=DQL (3)
                  ELSEIF(NSSN.EQ.3) THEN
                     AKWLDW(NR,1,NSSN,NSVN)=CHEL(5)
                     ADWLDW(NR,1,NSSN,NSVN)=DL  (5)
                     AKWLDW(NR,2,NSSN,NSVN)=CHIL(5)
                     AKWLDW(NR,3,NSSN,NSVN)=CHQL(5)
                     ADWLDW(NR,3,NSSN,NSVN)=DQL (5)
                  ENDIF
               ELSEIF(NSVN.EQ.2) THEN
                  IF(NSSN.EQ.1) THEN
                     AKWLDW(NR,1,NSSN,NSVN)=CHEL(2)
                     ADWLDW(NR,1,NSSN,NSVN)=DL  (2)
                     AKWLDW(NR,2,NSSN,NSVN)=CHIL(2)
                     AKWLDW(NR,3,NSSN,NSVN)=CHQL(2)
                     ADWLDW(NR,3,NSSN,NSVN)=DQL (2)
                  ELSEIF(NSSN.EQ.2) THEN
                     AKWLDW(NR,1,NSSN,NSVN)=CHEL(1)
                     ADWLDW(NR,1,NSSN,NSVN)=DL  (1)
                     AKWLDW(NR,2,NSSN,NSVN)=CHIL(1)
                     AKWLDW(NR,3,NSSN,NSVN)=CHQL(1)
                     ADWLDW(NR,3,NSSN,NSVN)=DQL (1)
                  ELSEIF(NSSN.EQ.3) THEN
                     AKWLDW(NR,1,NSSN,NSVN)=CHEL(4)
                     ADWLDW(NR,1,NSSN,NSVN)=DL  (4)
                     AKWLDW(NR,2,NSSN,NSVN)=CHIL(4)
                     AKWLDW(NR,3,NSSN,NSVN)=CHQL(4)
                     ADWLDW(NR,3,NSSN,NSVN)=DQL (4)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
C
      RETURN
      END
C
C     ****************************************************************
C
      SUBROUTINE TR_WEILAND_BRIDGE
     &     (EN,EI,EE,TAU,FL,FT,PDEV,PST,
     &      BQ,EQ,ENQ,Z,BETAE,AZ,COL,EL,TE,TAUZ,PR,
     &      Q,S,E,N,LIST,RNEQ,ALA,KAPPA,TR,RIW,RISB,BTOR,
     &      SEARCH,MA,GKL,WEXBL,ROTL,IKL,ISTL,
     &      CHI,CHE,D,CHQ,DQ,SCHIL,SCHEL,SDL,SCHQL,SDQL)
C
      IMPLICIT NONE
      INTEGER FNBI,I,IW,IX
      REAL*8 CT(0:4),U(5,100)
      COMPLEX*16 ZZ(10),RP(10),W
      PARAMETER (FNBI=20)
      INTEGER IR,PDEV,IK,IST,ITS,ITL,ITERA,ITC,ISB
      REAL*8 RAT,RT,ENI,EEI,TE,BTOR
      REAL*8 A,B,C,DS,ZX,TR,RIW,RISB
      REAL*8 EI,TAU,FL,THRD,TVR,STR,XIH
      REAL*8 EN,ENH,ENN,PST,RFL,H,H1
      REAL*8 LBT,TAUI,FTR,EIH,EEH
      REAL*8 FT,EE
      REAL*8 BQ,EQ,ENQ,Z,GQ,BF,ZE,TAUZ,ZEFF,NQ,NI,G
      REAL*8 CETAIN(32),BETAE,MA
      REAL*8 E,Q,S,CS,KAPPA,KPPA,RAV
      REAL*8 ALP,ALF,KPC,PR,FTRT,WST,D1,SI,KIQ,KXQ
      REAL*8 N,WR,WI,RNEQ,WDE,EPS,WRS,WIS
      REAL*8 SCHI,SCHE,SD,SCHQ,SDQ,EA,HPT(5),GK,DTOT
      REAL*8 ETE,ETI,ETQ,AZ,AZL,ALA,ALAF,GAV
      REAL*8 VEI,VEF,BTA,COL,EL,EM,LAMB
      REAL*8 CHI(5),CHE(5),D(5),CHQ(5),DQ(5)
      REAL*8 GNH,GNE,GNQ,GTH,GTE,GTQ
      INTEGER LPRINTIN,NDIM,NEQ,NDISP,IRET
      REAL*8 EIC,EEC,ENC,TAUC,FLC,FTC,EQC,ENQC,BETAEC,TAUZC,QC,SC
      REAL*8 ENHC,LIST,ZFS,KPS,CHIC,R
      REAL*8 WEXB,ROT
      REAL*8 WZIMAX,TOL,SHPE,SCHEF,DEF
      REAL*8 LTH,LTE,LN,LNH,LTQ,LNQ
      REAL*8 ZVR(10,10),ZVI(10,10)
C      REAL*8 SEARCHMODE,SEARCH
      INTEGER SEARCHMODE
      REAL*8 SEARCH
      INTEGER ISHW,IKL,ISTL
      REAL*8 GKL,WEXBL,ROTL,SCHIL,SCHEL,SDL,SCHQL,SDQL
      COMPLEX*16 HQ,WZ,WZP
      COMMON/PARAM/ CT
      COMMON/IRET/ IRET
      COMMON/PAR/ THRD,FTR,STR,RFL,D1,XIH,IX,IW
      COMMON/EFF DIFF/ SCHI,SCHE,SD,SCHQ,SDQ
      COMMON/GRAD/ ENC,ENHC,EIC,EEC
      COMMON/COLL/ BTA,VEF
      COMMON/IMP/ BF,GQ,SI,ZE,TAUZC,ENQC,EQC,KIQ,KXQ,ZFS
      COMMON/ETAWN6/ CETAIN,ETE,ETI,ETQ,TAUC,AZL,
     &FTC,FLC,LPRINTIN,NDIM
      COMMON/W/ WR,WI
      COMMON/BETAE/ BETAEC,QC,SC,CS,EM
      COMMON/TEST/ ALP,ALF,KPC,KPS
      COMMON/HQ/ H1,HQ,GAV,WZ,WZP,WZIMAX
      COMMON/SHAFRS/ ALAF
      COMMON/KAPPA/ KPPA,RAV
      COMMON/IK/ IK,IST,ITC,ITL,ITS,ITERA,TOL
      COMMON/TP/ EA,HPT
      COMMON/WROT/ WEXB,ROT
      COMMON/GRKVOT/ GK
      COMMON/ZV/ ZVR,ZVI
      COMMON/ZZ/ ZZ
      COMMON/NEQ/ NEQ
      COMMON/EM/ SHPE,CHIC,SCHEF,DEF
      COMMON/LT/ LTH,LTE,LN,LTQ,LNQ
      COMMON/ISB/ ISB
      COMMON/SEARCHMODE/ SEARCHMODE
      COMMON/WDE/ WDE
C
C     ISHW : controling printout in TR_WEILAND_BRIDGE (0: off)
      ISHW=0
C
      GK   = GKL
      WEXB = WEXBL
      ROT  = ROTL
      IK   = IKL
      IST  = ISTL
C
      ZFS = 0.D0
C
      THRD=1.D0/3.D0
      TVR=2.D0*THRD
      FTR=5.D0/3.D0
      STR=7.D0/3.D0
      BTA=1.5D0
      EPS=E
      SEARCHMODE=INT(SEARCH)
      AZL=AZ
      EM=EL
      ZE=Z
      LBT=EI/EN
      ENI=EN
      EEI=EE
      BF=BQ
      G=1.D0-Z*BQ
      NI=G*N
      ENN=1.D0-Z*BQ*EN/ENQ
      IF(ABS(ENN).GE.0.001D0) GO TO 89
      IF(ENN.LT.0.D0) ENN=-0.001D0
      IF(ENN.GT.0.D0) ENN= 0.001D0
   89 CONTINUE
      ENH=G*EN/ENN
      EIC=EI
      EEC=EE
      TAUC=TAU
      FLC=DSQRT(FL)
      FTC=FT
      EQC=EQ
      ENQC=EN
      ENHC=ENH
      BETAEC=BETAE
      TAUZC=TAUZ
      QC=Q
      SC=S
      NEQ=INT(RNEQ)
      ISB=INT(RISB)
      ALAF=ALA
      NDISP=NEQ
      NDIM=5
      KPPA=KAPPA
      EA=E
      IW=INT(RIW)
      R=PR/E
      D1=6.462D0*DSQRT(MA)/(R*BTOR**2)
C
      IX=IK
C
      DO I=1,31
         CETAIN(I)=0.D0
      ENDDO
      CETAIN(32)=1.D-15
      LPRINTIN=INT(LIST)
C
      FTRT=FTR/TAU
      TAUI=1.D0/TAU
      RFL=SQRT(FL)
      RAT=EI/EN
      GQ=1.D0-Z*BQ
      NQ=BQ*N
      NI=GQ*N
      ZEFF=(NI+Z*Z*NQ)/N
C
      EIC=EI
      EEC=EE
      ENC=EN
      FTC=FT
      FLC=DSQRT(FL)
      BETAEC=BETAE
      KPPA=KAPPA
      LN=0.5D0*R*EN/PR
      LTE=LN/EE
      ENN=1.D0-Z*BQ*EN/ENQ
      IF(ABS(ENN).GE.0.001D0) GO TO 90
      IF(ENN.LT.0.) ENN=-0.001D0
      IF(ENN.GT.0.D0) ENN=0.001D0
   90 CONTINUE
      ENH=GQ*EN/ENN
      ENHC=ENH
      LNH=0.5D0*R*ENH/PR
      LNQ=0.5D0*R*ENQ/PR
      LTH=LNH/EI
      LTQ=LNQ/EQ
C
C
      TAUI=1.D0/TAU
      RT=EN/ENI
C      EI=RAT*EN
C      EE=RT*RT*EEI
      EIH=EI-7.D0/3.D0+FTR*EN
      EEH=EE-7.D0/3.D0+FTR*EN
      GQ=1.D0-Z*BQ
      GNE=2.D0/EN
      GNH=2.D0/ENH
      GNQ=2.D0/ENQ
      GTE=EE*GNE
      GTH=EI*GNH
      GTQ=EQ*GNQ
      ETI=ENH/EI
      ETE=EN/EE
      ETQ=ENQ/EQ
C
      KIQ=ENQ/ENH
      KXQ=EN/ENQ
C
      NI=GQ*N
      CS=3.09501D5*DSQRT(TE/MA)
      IF(ISHW.NE.0) THEN
         WRITE(*,00126) EN,EI,EE,FL,TAU
00126    FORMAT(2X,'EN=',F8.3,' EI=',F8.3,' EE=',F8.3,' FL=',F8.3,
     &        ' TAU=',G12.4)
         WRITE(*,00299) ENQ,ENH
00299    FORMAT(2X,'ENQ=',G11.3,' ENH=',G11.3)
C         WRITE(*,00127) ENI,EEI,RAT,RT
C00127    FORMAT(2X,'ENI=',G10.4,' EEI=',G10.4,' RAT=',G10.4
C     &        ,' RT=',G10.4)
         WRITE(*,00129) BETAE
00129    FORMAT(2X,'BETAE=',G11.3)
         WRITE(*,00131) Q,S,CS
00131    FORMAT(2X,'q=',G12.4,' S=',G12.4,' CS=',G12.4)
         WRITE(*,00132) FT
00132    FORMAT(2X,'FT=',G12.4)
         WRITE(*,00134) ALAF
00134    FORMAT(2X,' Ballooning alpha =',G11.3)
      ENDIF
C
      WST=DSQRT(FL)*CS/(PR*ABS(LN))
      WDE=ABS(EN)*WST
      LAMB=15.95D0-DLOG(DSQRT(N)/TE)
      VEI=9.19D2*NI*LAMB/TE**1.5D0
      VEF=VEI/(EPS*WDE)
      VEF=COL*VEF
C
      IF(ISHW.NE.0) THEN
         WRITE(*,00130) WST,VEF,VEI,ZEFF,COL
00130    FORMAT(2X,'WST=',G11.3,' VEF=',G11.3,' VEI=',G11.3,
     &        ' ZEFF=',G11.3,' COL=',G11.3)
      ENDIF
C
      H=0.5D0*ABS(S)/q
      IF(ISHW.NE.0) THEN
         WRITE(*,00133) H
00133    FORMAT(2X,'H=',G12.4)
      ENDIF
C   -----------------------------------------------
C
      A=1.D0-EN*(1.D0+10.D0/(3.D0*TAU))-FL*(1.D0+EI+5.D0*EN/3.D0)/TAU
      B=EI-7.D0/3.D0+5.D0*EN*(1.D0+1.D0/TAU)/3.D0
     &              +5.D0*FL*(1.D0+EI)/(3.D0*TAU)
      B=B*EN/TAU
      C=A/(2.D0*(1.D0+FL))
      DS=C*C-B/(1.D0+FL)
      IF(DS.LT.0.D0) GOTO 140
      WR=C+SQRT(DS)
      WI=0.D0
      GO TO 160
  140 WR=C
      WI=SQRT(-DS)
      IF(ISHW.NE.0) WRITE(*,170) WR,WI
  160 CONTINUE
  170 FORMAT(2X,'WR=',F7.3,' WI=',F7.3)
C
      ITC=1
      ITL=80
      TOL=0.01D0
C      IST=1
c
      CALL disp9t(NDISP,ZZ)
c
      IF(ISHW.NE.0) THEN
         WRITE(*,174) ISB
 174     FORMAT(' ISB=',I5)
         WRITE(*,175) ITC,ITS,ITERA
 175     FORMAT('  ITC=',I5,' ITS=',I5,' ITER=',I5)
         WRITE(*,177) WZ,WZP
 177     FORMAT('  WZ=',2G11.3,' WZP=',2G11.3)
      ENDIF
C
      IR=0
C
      DO 00199 I=1,NDISP
      ZX=DIMAG(ZZ(I))
      IF(ZX.LE. 0.001) GOTO 00199
      IR=IR+1
      RP(IR)=ZZ(I)
00199 CONTINUE
C
      IF(ISHW.NE.0) THEN
         WRITE(*,310) IR
 310     FORMAT(2X,' IR=',I4)
         WRITE(*,00134) ALAF
      ENDIF
C
      DO 0200 I=1,IR
      W=RP(I)
      WR=EN*DREAL(W)
      WI=EN*DIMAG(W)
      IF(ISHW.NE.0) WRITE(*,311) WR,WI,I
      WRS=WDE*DREAL(W)
      WIS=WDE*DIMAG(W)
      IF(ISHW.NE.0) WRITE(*,321) WRS,WIS
 0200 CONTINUE
  311 FORMAT(//,2X,'WR=',G11.3,' WI=',G11.3,' I=',I5)
  321 FORMAT(' WRS=',G11.3,' WIS=',G11.3)
C
      WZ=EN*WZ
      HQ=EN*HQ
      IF(ISHW.NE.0) THEN
         WRITE(*,00128) ALF,ALP,WZ,KAPPA
00128    FORMAT(/,2X,' ALF=',G11.3,' ALP=',G11.3,' WZ=',2G11.3,
     &        /,' KAPPA=',G11.3)
         WRITE(*,312) H1,HQ,GAV,RAV
 312     FORMAT(//,2X,'H1=',G12.4,' HQ=',2G12.4,' GAV=',G12.4,
     &        /,' RAV=',G12.4)
      ENDIF
      U(1,IX)=TAUI*TE
      U(2,IX)=TE
      U(3,IX)=N
      U(4,IX)=TE/TAUZ
      U(5,IX)=BQ*N
C
      CALL DIFF(RP,IR,TAUI,FT,U,CHI,CHE,D,CHQ,DQ)
C
      IF(ISHW.NE.0) THEN
         WRITE(*,330) SCHI,SCHE,SD,SCHQ,SDQ
 330     FORMAT(/,2X,'CHIEFF=',G11.3,' CHEEFF=',G11.3,' DEFF=',G11.3,
     &        ' CHQEFF=',G11.3,' DQEFF=',G11.3)
C     
         WRITE(*,331) CHE(2),SCHEF
 331     FORMAT(' CHE(2)=',G11.3,' SCHEF=',G11.3)
         DTOT=SD+DEF
         WRITE(*,332) D(3),DEF,DTOT
 332     FORMAT(' D(3)=',G11.3,' DEF=',G11.3,' DTOT=',G11.3)
         WRITE(*,335) XIH,SHPE,CHIC,SCHEF
 335     FORMAT('  XIH=',G11.3,' SHPE=',G11.3,' CHIC=',G11.3,
     &        ' SCHEF=',G11.3)
      ENDIF
C00150 CONTINUE
      SCHIL = SCHI
      SCHEL = SCHEL
      SDL   = SD
      SCHQL = SCHQ
      SDQL  = SDQ
C
      RETURN
      END
