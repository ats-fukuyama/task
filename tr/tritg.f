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
C     ENL    : 2 Lne/Lb (Lb=R) (epsilon_ne)
C     EIL    : Lni/Lti (eta_i)
C     EEL    : Lne/Lte (eta_e)
C     TAUL   : Te/Ti (tau_i)
C     FLL    : Finite Larmor Radius parameter: (k_perp*rhos)**2=0.1
C     FTL    : trapped particle fraction
C     BQL    : fraction of an impurity
C     EQL    : Lnq/Ltq (q represents an impurity) (eta_q)
C     ENQL   : 2 Lnq/Lb (epsilon_nq)
C     ZL     : Z value for an impurity
C     BETAEL : electron beta
C     AZL    : impurity mass number
C     COLL   : factor multiplying collisionality (0 or 1)
C     ELL    : factor multiplying electromagnetic effects (0 or 1)
C     TEL    : electron temperature [keV]
C     TAUZL  : Te/Tq (tau_q)
C     RA     : minor radius (a)[m]
C     QL     : safety factor
C     SL     : magnetic shear
C     EPS    : local inverse aspect ratio (r/R)
C     EPSA   : inverse aspect ratio (a/R)
C     RNL    : electron density [10^19(20?) m^-3]
C     RLIST  : controlling printout in disp9t (0: off)
C     RNEQL  : number of equations (NDISP)
C     RKAP   : ellipticity
C     RIWL   : controlling printout in DIFFTD (0: off)
C     RISBL  : If ISB=1, we use the strong ballooning approximation
C              (GAV=1).
C     BB     : toroidal magnetic field [T]
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
C               and Nq equations [m^2/s (same as above)]
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
      MDDW=1
      IF(NT.EQ.0) THEN
         IST=1
      ELSE
         IST=0
      ENDIF
      ZL    = PZ(3)
      AZL   = PA(3)
      COLL  = 1.D0
      ELL   = 1.D0
      RLIST = 1.D0
      RNEQL = 9.D0 
      RIWL  = 2.D0
      RISBL = 2.D0
      SEARCH= 2.D0
      PMA   = PA(2)
      ROTL  = 1.D0
      EPSA  = RA/RR
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
C         ENL   = ABS(2.D0*(SLNEL/SLBL ))
C         EIL   = ABS(      SLNIL/SLTIL)
C         EEL   = ABS(      SLNEL/SLTEL)
         TAUL  = (RT(NR+1,1)+RT(NR,1))/(RT(NR+1,2)+RT(NR,2))
         FLL   = 1.D-1
         FTL   = FTPF(MDLTPF,EPS)
         BQL   = (RN(NR+1,3)+RN(NR,3))/(RN(NR+1,1)+RN(NR,1))
         EQL   =       SLNQL/SLTQL
         ENQL  =-2.D0*(SLNQL/SLBL )
C         EQL   = ABS(      SLNQL/SLTQL)
C         ENQL  = ABS(2.D0*(SLNQL/SLBL ))
         BETAEL= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
     &          *RKEV*1.D20/(BB**2/(2.D0*RMU0))
         TEL   = 0.5D0*(RT(NR+1,1)+RT(NR,1))
         TAUZL = (RT(NR+1,1)+RT(NR,1))/(RT(NR+1,3)+RT(NR,3))
         QL    = QP(NR)
         IF(NR.EQ.1) THEN
            DQ = (QP(NR+1)-QP(NR))/1.5D0*DRL
         ELSE
            DQ = (QP(NR+1)-QP(NR-1))/2.D0*DRL
         ENDIF
         SL    = RR*EPS*DQ/QL
         RNL   = 0.5D0*(RN(NR+1,1)+RN(NR,1))*1.D1
C
         RGKL  = AR1RHOG(NR)/(RA*AR2RHOG(NR))
         WEXBL = WEXB(NR)
C
C         COEF = PZ(2)**2*AEE**4*1.D20
C     &         /(6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)*RKEV**1.5D0)
C         RLAMB =15.2D0-0.5D0*DLOG(RNL)+DLOG(TEL)
C         VEI  = COEF*RNL*RLAMB/TEL**1.5D0
C         write(6,*) NR,COEF,VEI
C
         CALL TR_WEILAND_BRIDGE
     &     (ENL,EIL,EEL,TAUL,FLL,FTL,
     &      BQL,EQL,ENQL,ZL,BETAEL,AZL,COLL,ELL,TEL,TAUZL,RA,
     &      QL,SL,EPS,EPSA,RNL,RLIST,RNEQL,RKAP,RIWL,RISBL,BB,
     &      SEARCH,PMA,RGKL,WEXBL,ROTL,NR,IST,
     &      CHIL,CHEL,DL,CHQL,DQL,SCHI,SCHE,SD,SCHQ,SDQ)
C
         CALL WEILAND_COEF(NR,CHIL,CHEL,DL,CHQL,DQL,
     &                        SCHI,SCHE,SD,SCHQ,SDQ)
C       write(6,'(I3,5F15.7)') NR,CHEL(1),CHEL(2),CHEL(3),CHEL(4),CHEL(5)
C       write(6,'(I3,5F15.7)') NR,CHIL(1),CHIL(2),CHIL(3),CHIL(4),CHIL(5)
C       write(6,'(I3,5F15.7)') NR,CHQL(1),CHQL(2),CHQL(3),CHQL(4),CHQL(5)
C       write(6,'(I3,5F15.7)') NR,DL(1),DL(2),DL(3),DL(4),DL(5)
C         if(nr.eq.4) write(6,'(I4,5F15.7)') NT,DL(3),CHQL(4)
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
C         ENL   = ABS(2.D0*(SLNEL/SLBL ))
C         EIL   = ABS(      SLNIL/SLTIL)
C         EEL   = ABS(      SLNEL/SLTEL)
         TAUL  = PTS(1)/PTS(2)
         FLL   = 1.D-1
         FTL   = FTPF(MDLTPF,EPS)
         BQL   = PNSS(3)/PNSS(1)
         EQL   =       SLNQL/SLTQL
         ENQL  =-2.D0*(SLNQL/SLBL )
C         EQL   = ABS(      SLNQL/SLTQL)
C         ENQL  = ABS(2.D0*(SLNQL/SLBL ))
         BETAEL= PNSS(1)*PTS(1)*RKEV*1.D20/(BB**2/(2.D0*RMU0))
         TEL   = PTS(1)
         TAUZL = PTS(1)/PTS(3)
         QL    = QP(NR)
         DQ    = (QP(NR)-QP(NR-1))*DRL
         SL    = RR*EPS*DQ/QL
         RNL   = PNSS(1)*1.D1
C
         RGKL  = AR1RHOG(NR)/(RA*AR2RHOG(NR))
         WEXBL = WEXB(NR)
C
         CALL TR_WEILAND_BRIDGE
     &     (ENL,EIL,EEL,TAUL,FLL,FTL,
     &      BQL,EQL,ENQL,ZL,BETAEL,AZL,COLL,ELL,TEL,TAUZL,RA,
     &      QL,SL,EPS,EPSA,RNL,RLIST,RNEQL,RKAP,RIWL,RISBL,BB,
     &      SEARCH,PMA,RGKL,WEXBL,ROTL,NR,IST,
     &      CHIL,CHEL,DL,CHQL,DQL,SCHI,SCHE,SD,SCHQ,SDQ)
C
         CALL WEILAND_COEF(NR,CHIL,CHEL,DL,CHQL,DQL,
     &                        SCHI,SCHE,SD,SCHQ,SDQ)
C
      RETURN
      END
C
C     ****************************************************************
C
      SUBROUTINE WEILAND_COEF(NR,CHIL,CHEL,DL,CHQL,DQL,
     &                           SCHI,SCHE,SD,SCHQ,SDQ)
C
      INCLUDE 'trcomm.inc'
      DIMENSION CHIL(5),CHEL(5),DL(5),CHQL(5),DQL(5)
C
C     The diagonal value of the transport coefficient matrix 
C     is set to be zero if it becomes negative.
C
      IF(DL  (3).LT.0.D0) DL  (3)=0.D0
      IF(CHEL(2).LT.0.D0) CHEL(2)=0.D0
      IF(CHIL(1).LT.0.D0) CHIL(1)=0.D0
      IF(DQL (5).LT.0.D0) DQL (5)=0.D0
      IF(CHQL(4).LT.0.D0) CHQL(4)=0.D0
C
C     It is assumed that De=Di in the followings.
C
      IF(PA(3).EQ.3.D0) THEN
C         AKDW(NR,1)=SCHE
         AKDW(NR,1)=CHEL(2)
         DO NS=2,NSM
C            AKDW(NR,NS)=SCHI
            AKDW(NR,NS)=CHIL(1)
         ENDDO
         DO NS=1,NSM
C            ADDW(NR,NS)=SD
            ADDW(NR,NS)=DL(3)
         ENDDO
         IF(MDLWLD.NE.0) THEN
            ADDWD(NR,1,1)=DL(3)
            ADDWP(NR,1,1)=CHEL(3)
            DO NS=2,NSM
               ADDWD(NR,NS,1)=DL(3)
               ADDWP(NR,NS,1)=CHIL(3)
            ENDDO
            AKDWD(NR,1,1)=DL(2)
            AKDWP(NR,1,1)=CHEL(2)
            DO NS=2,NSM
               AKDWD(NR,NS,1)=DL(2)
               AKDWP(NR,NS,1)=CHIL(2)
            ENDDO
            DO NS1=2,NSM
               ADDWD(NR,1,NS1)=DL(3)
               ADDWP(NR,1,NS1)=CHEL(3)
               DO NS=2,NSM
                  ADDWD(NR,NS,NS1)=DL(3)
                  ADDWP(NR,NS,NS1)=CHIL(3)
               ENDDO
               AKDWD(NR,1,NS1)=DL(1)
               AKDWP(NR,1,NS1)=CHEL(1)
               DO NS=2,NSM
                  AKDWD(NR,NS,NS1)=DL(1)
                  AKDWP(NR,NS,NS1)=CHIL(1)
               ENDDO
            ENDDO
         ENDIF
      ELSE
         AKDW(NR,1)=CHEL(2)
         AKDW(NR,2)=CHIL(1)
         AKDW(NR,3)=CHQL(4)
         AKDW(NR,4)=CHQL(4)
         ADDW(NR,1)=DL(3)
         ADDW(NR,2)=DL(3)
C         ADDW(NR,2)=0.D0
         ADDW(NR,3)=DQL(5)
         ADDW(NR,4)=DQL(5)
         IF(MDLWLD.NE.0) THEN
            ADDWD(NR,1,1)=DL(3)
            ADDWP(NR,1,1)=CHEL(3)
            ADDWD(NR,2,1)=DL(3)
            ADDWP(NR,2,1)=CHIL(3)
            DO NS=3,NSM
               ADDWD(NR,NS,1)=DQL(3)
               ADDWP(NR,NS,1)=CHQL(3)
            ENDDO
            AKDWD(NR,1,1)=DL(2)
            AKDWP(NR,1,1)=CHEL(2)
            AKDWD(NR,2,1)=DL(2)
            AKDWP(NR,2,1)=CHIL(2)
            DO NS=3,NSM
               AKDWD(NR,NS,1)=DQL(2)
               AKDWP(NR,NS,1)=CHQL(2)
            ENDDO
            ADDWD(NR,1,2)=DL(3)
            ADDWP(NR,1,2)=CHEL(3)
            ADDWD(NR,2,2)=DL(3)
            ADDWP(NR,2,2)=CHIL(3)
            DO NS=3,NSM
               ADDWD(NR,NS,2)=DQL(3)
               ADDWP(NR,NS,2)=CHQL(3)
            ENDDO
            AKDWD(NR,1,2)=DL(1)
            AKDWP(NR,1,2)=CHEL(1)
            AKDWD(NR,2,2)=DL(1)
            AKDWP(NR,2,2)=CHIL(1)
            DO NS=3,NSM
               AKDWD(NR,NS,2)=DQL(1)
               AKDWP(NR,NS,2)=CHQL(1)
            ENDDO
            DO NS1=3,NSM
               ADDWD(NR,1,NS1)=DL(5)
               ADDWP(NR,1,NS1)=CHEL(5)
               ADDWD(NR,2,NS1)=DL(5)
               ADDWP(NR,2,NS1)=CHIL(5)
               DO NS=3,NSM
                  ADDWD(NR,NS,NS1)=DQL(5)
                  ADDWP(NR,NS,NS1)=CHQL(5)
               ENDDO
               AKDWD(NR,1,NS1)=DL(4)
               AKDWP(NR,1,NS1)=CHEL(4)
               AKDWD(NR,2,NS1)=DL(4)
               AKDWP(NR,2,NS1)=CHIL(4)
               DO NS=3,NSM
                  AKDWD(NR,NS,NS1)=DQL(4)
                  AKDWP(NR,NS,NS1)=CHQL(4)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
C
      RETURN
      END
C
C     ****************************************************************
C
      SUBROUTINE TR_WEILAND_BRIDGE
     &     (EN,EI,EE,TAU,FL,FT,
     &      BQ,EQ,ENQ,Z,BETAE,AZ,COL,EL,TE,TAUZ,PR,
     &      Q,S,EPSR,E,N,LIST,RNEQ,KAPPA,RIW,RISB,BTOR,
     &      SEARCH,MA,GKL,WEXBL,ROTL,IKL,ISTL,
     &      CHI,CHE,D,CHQ,DQ,SCHIL,SCHEL,SDL,SCHQL,SDQL)
C
      IMPLICIT NONE
      INTEGER I,IW,IX
      REAL*8 U(5,100)
      COMPLEX*16 ZZ(10),RP(10),W
      INTEGER IR,IK,IST,ITS,ITL,ITERA,ITC,ISB
      REAL*8 TE,BTOR
      REAL*8 A,B,C,DS,ZX,RIW,RISB
      REAL*8 EI,TAU,FL,THRD,TVR,STR,XIH
      REAL*8 EN,ENN,RFL,H,H1
      REAL*8 LBT,TAUI,FTR,EIH,EEH
      REAL*8 FT,EE
      REAL*8 BQ,EQ,ENQ,Z,GQ,BF,ZE,TAUZ,ZEFF,NQ,NI,G
      REAL*8 CETAIN(32),BETAE,MA
      REAL*8 E,Q,S,CS,KAPPA,KPPA,RAV
      REAL*8 ALP,ALF,KPC,PR,FTRT,WST,D1,SI,KIQ,KXQ
      REAL*8 N,WR,WI,RNEQ,WDE,EPS,EPSR,WRS,WIS
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
      INTEGER ICP,IKL,ISTL
      REAL*8 GKL,WEXBL,ROTL,SCHIL,SCHEL,SDL,SCHQL,SDQL
      COMPLEX*16 HQ,WZ,WZP
      COMMON/IRET/ IRET
      COMMON/PAR/ THRD,FTR,STR,RFL,D1,XIH,IX,IW
      COMMON/EFF DIFF/ SCHI,SCHE,SD,SCHQ,SDQ
      COMMON/GRAD/ ENC,ENHC,EIC,EEC
      COMMON/COLL/ BTA,VEF
      COMMON/IMP/ BF,GQ,SI,ZE,TAUZC,ENQC,EQC,KIQ,KXQ,ZFS
      COMMON/ETAWN6/ CETAIN,ETE,ETI,ETQ,TAUC,AZL,
     &FTC,FLC,LPRINTIN,NDIM
C      COMMON/W/ WR,WI ! not used
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
      COMMON/WDE/ WDE ! added
C
C     ICP : controlling printout in TR_WEILAND_BRIDGE (0: off)
      ICP=0
      IF(ICP.NE.0) THEN
         WRITE(6,*) '//////////////////////////////////////////////////'
      ENDIF
C
      GK   = GKL     ! GRHO1/(a*GRHO2)
      WEXB = WEXBL   ! EXB shearing rate according to Hahm and Burrell
      ROT  = ROTL    ! factor that multiplies the EXB shearing rate
      IK   = IKL     ! the space profile index for WZL(IK)
      IST  = ISTL    ! 1:analytical form used, other:previous value used
C
      ZFS = 0.D0     ! the product of charge and fraction (to Ne) of fast particles
C
      THRD=1.D0/3.D0 ! fixed coefficient (THiRD)
      TVR=2.D0*THRD  ! fixed coefficient (Two thiRd)
      FTR=5.D0/3.D0  ! fixed coefficient (Five ThiRd)
      STR=7.D0/3.D0  ! fixed coefficient (Seven ThiRd)
      BTA=1.5D0      ! fixed, parameter used in the collision model
C
      EPS=EPSR       ! local inverse aspect ratio
      SEARCHMODE=INT(SEARCH) ! the way to search eigenvalues
      AZL=AZ         ! Aq (atomic number of impurity)
      EM=EL          ! factor that multiplies electromagnetic effects
      ZE=Z           ! Zq
      BF=BQ          ! nq/ne
C
      G=1.D0-Z*BQ    ! transforming factor from ne to ni
      ENN=1.D0-Z*BQ*EN/ENQ ! (Lne/Lni)*(ni/ne)
      IF(ABS(ENN).GE.0.001D0) GO TO 89
      IF(ENN.LT.0.D0) ENN=-0.001D0
      IF(ENN.GT.0.D0) ENN= 0.001D0
   89 CONTINUE
      ENC=EN         ! 2*Lne/R
      ENHC=G*EN/ENN  ! 2*Lni/R
      ENQC=ENQ       ! 2*Lnq/R
C      ENQC=EN        ! 2*Lne/R (original definition)
C
      EIC=EI         ! Lni/Lti
      EEC=EE         ! Lne/Lte
      EQC=EQ         ! Lnq/Ltq
      TAUC=TAU       ! Te/Ti
      TAUZC=TAUZ     ! Te/Tq
C
      FLC=DSQRT(FL)  ! FLR parameter (sqrt((Kperp*rhos)**2=0.1))
      FTC=FT         ! trapped particle fraction
      BETAEC=BETAE   ! electron beta
      QC=Q           ! safety factor
      SC=S           ! magnetic shear
      NEQ=INT(RNEQ)  ! number of equations (NDISP)
      ISB=INT(RISB)  ! ballooning parameter
      NDISP=NEQ      ! number of equations (NDISP)
      NDIM=5         ! dimension of transport matrix
      KPPA=KAPPA     ! elongation
      EA=E           ! inverse aspect ratio
      IW=INT(RIW)    ! controls printout in DIFFTD
      R=PR/E         ! major radius
      D1=6.462D0*DSQRT(MA)/(R*BTOR**2) ! machine dependent parameter
C
      IX=IK          ! the space profile index for U(J,IX)
C
      DO I=1,31
         CETAIN(I)=0.D0 ! control vector
      ENDDO
C      CETAIN(32)=1.D-15
      CETAIN(32)=0.001  ! accuracy parameter in the NAG routine
      LPRINTIN=INT(LIST)
C
      TAUI=1.D0/TAU  ! Ti/Te
      RFL=SQRT(FL)   ! FLR parameter (same with FLC)
      GQ=1.D0-Z*BQ   ! transforming factor from ne to ni (same with G)
      NQ=BQ*N        ! nq
      NI=GQ*N        ! ni
      ZEFF=(NI+Z*Z*NQ)/N ! Zeff
C
      LN =0.5D0*R*ENC/PR  ! Lne/a
      LNH=0.5D0*R*ENHC/PR ! Lni/a
      LNQ=0.5D0*R*ENQC/PR ! Lnq/a
      LTE=LN/EEC          ! Lte/a
      LTH=LNH/EIC         ! Lti/a
      LTQ=LNQ/EQC         ! Ltq/a
C
      ETI=ENHC/EIC   ! 2*Lti/R
      ETE=ENC/EEC    ! 2*Lte/R
      ETQ=ENQC/EQC   ! 2*Ltq/R
C
      KIQ=ENQC/ENHC  ! Lnq/Lni
      KXQ=ENC/ENQC   ! Lne/Lnq
C
      CS=3.095D5*DSQRT(TE/MA) ! sound speed
      IF(ICP.NE.0) THEN
         WRITE(*,00126) EN,EI,EE,FL,TAU
00126    FORMAT(2X,'EN=',F8.3,' EI=',F8.3,' EE=',F8.3,' FL=',F8.3,
     &        ' TAU=',G12.4)
         WRITE(*,00299) ENQ,ENHC
00299    FORMAT(2X,'ENQ=',G11.3,' ENH=',G11.3)
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
      WST=DSQRT(FL)*CS/(PR*ABS(LN))  ! diamagnetic frequency (omega_star)
      WDE=ABS(EN)*WST           ! curvature drift frequency (omega_drift_electron)
      LAMB=15.95D0-DLOG(DSQRT(N)/TE) ! coulomb logarithm (lambda)
      VEI=9.19D2*NI*LAMB/TE**1.5D0   ! collisionality (nu_electron_ion)
      VEF=COL*VEI/(EPS*WDE)     ! effective electron ion collision frequency for trapped electrons, normalized by omega_de
      WEXB=WEXB/WDE ! ExB shearing rate should be normalized with WDE
C
      IF(ICP.NE.0) THEN
         WRITE(*,00130) WST,VEF,VEI,ZEFF,COL
00130    FORMAT(2X,'WST=',G11.3,' VEF=',G11.3,' VEI=',G11.3,
     &        ' ZEFF=',G11.3,' COL=',G11.3)
C
         H=0.5D0*ABS(S)/q       ! not used
         WRITE(*,00133) H
00133    FORMAT(2X,'H=',G12.4)
      ENDIF
C   -----------------------------------------------
C
      IF(ICP.NE.0) THEN
         A=1.D0-EN*(1.D0+10.D0/(3.D0*TAU))-FL*(1.D0+EI+5.D0*EN/3.D0)/TAU
         B=EI-7.D0/3.D0+5.D0*EN*(1.D0+1.D0/TAU)/3.D0
     &                 +5.D0*FL*(1.D0+EI)/(3.D0*TAU)
         B=B*EN/TAU
         C=A/(2.D0*(1.D0+FL))
         DS=C*C-B/(1.D0+FL)
         IF(DS.LT.0.D0) GOTO 140
         WR=C+SQRT(DS)          ! not used
         WI=0.D0                ! not used
         GO TO 160
 140     WR=C                   ! not used
         WI=SQRT(-DS)           ! not used
         WRITE(*,170) WR,WI
 160     CONTINUE
 170     FORMAT(2X,'WR=',F7.3,' WI=',F7.3)
      ENDIF
C
      ITC=1      ! 1:iteration, other:previous eigenvalues used
      ITL=80     ! Maximum number of iterations
      TOL=0.01D0 ! relative error for convergence
C      IST=1
c
C      IF(IX.GE.48) write(6,'(I3,2F15.7)') IX,HQ
      CALL disp9t(NDISP,ZZ)
c
      IF(ICP.NE.0) THEN
         WRITE(*,174) ISB
 174     FORMAT(' ISB=',I5)
         WRITE(*,175) ITC,ITS,ITERA
 175     FORMAT('  ITC=',I5,' ITS=',I5,' ITER=',I5)
         WRITE(*,177) WZ,WZP
 177     FORMAT('  WZ=',2G11.3,' WZP=',2G11.3)
      ENDIF
C
      IR=0  ! the number of unstable roots (given by the following)
C
      DO 00199 I=1,NDISP
      ZX=DIMAG(ZZ(I))
      IF(ZX.LE. 0.001) GOTO 00199
      IR=IR+1
      RP(IR)=ZZ(I) ! unstable roots found by disp9t
00199 CONTINUE
C
      IF(ICP.NE.0) THEN
         WRITE(*,310) IR
 310     FORMAT(2X,' IR=',I4)
         WRITE(*,00134) ALAF
      ENDIF
C
      DO 0200 I=1,IR
      W=RP(I)
      WR=EN*DREAL(W) ! real part of unstable roots
      WI=EN*DIMAG(W) ! imaginary part of unstable roots
      IF(ICP.NE.0) WRITE(*,311) WR,WI,I
      WRS=WDE*DREAL(W)
      WIS=WDE*DIMAG(W)
      IF(ICP.NE.0) WRITE(*,321) WRS,WIS
 0200 CONTINUE
  311 FORMAT(//,2X,'WR=',G11.3,' WI=',G11.3,' I=',I5)
  321 FORMAT(' WRS=',G11.3,' WIS=',G11.3)
C
      WZ=EN*WZ
      HQ=EN*HQ
      IF(ICP.NE.0) THEN
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
      IF(ICP.NE.0) THEN
         WRITE(*,330) SCHI,SCHE,SD,SCHQ,SDQ
 330     FORMAT(/,2X,'CHIEFF=',G11.3,' CHEEFF=',G11.3,' DEFF=',G11.3,
     &        ' CHQEFF=',G11.3,' DQEFF=',G11.3)
C     
C         WRITE(*,331) CHE(2),SCHEF
C 331     FORMAT(' CHE(2)=',G11.3,' SCHEF=',G11.3)
         WRITE(*,331) CHQ(5),SCHQ
 331     FORMAT(' CHQ(5)=',G11.3,' SCHQ=',G11.3)
         DTOT=SD+DEF
         WRITE(*,332) D(3),DEF,DTOT
 332     FORMAT(' D(3)=',G11.3,' DEF=',G11.3,' DTOT=',G11.3)
         WRITE(*,335) XIH,SHPE,CHIC,SCHEF
 335     FORMAT('  XIH=',G11.3,' SHPE=',G11.3,' CHIC=',G11.3,
     &        ' SCHEF=',G11.3)
      ENDIF
C00150 CONTINUE
      SCHIL = SCHI
      SCHEL = SCHE
      SDL   = SD
      SCHQL = SCHQ
      SDQL  = SDQ
C
      RETURN
      END
C
C  
C     ***********************************************************
C
C            IFS/PPPL Model
C
C     ***********************************************************
C
      SUBROUTINE IFSPPPL_DRIVER(NRM,NSM,NSTM,NRMAX,RN,RR,DR,RJCB,QP,
     &                          S_AR,EPSRHO,RKPRHOG,RT,BB,AMM,AME,
     &                          PNSS,PTS,RNFL,RNFEDG,MDLUF,NSMAX,
     &                          AR1RHOG,AR2RHOG,AKDW)
C
      IMPLICIT NONE
C
      INTEGER NRM,NSM,NSTM,NRMAX,NR,MDLUF,NSMAX
      REAL*8 RN(NRM,NSM),RR,DR,RJCB(NRM),QP(NRM),S_AR(NRM),
     &       EPSRHO(NRM),RKPRHOG(NRM),RT(NRM,NSM),BB,AMM,AME,
     &       PNSS(NSM),PTS(NSM),RNFL(NRM),RNFEDG,
     &       AR1RHOG(NRM),AR2RHOG(NRM),AKDW(NRM,NSTM)
      integer switches(32), ipin, ipout, iptmp, screen, ii, ierr
      parameter (ipin=7,iptmp=8,ipout=9,screen=6)
      real znine, zncne, znbne, zrlt, zrln, zq, zshat, zeps,
     &       ne19, tekev, tikev, rmajor, grhoi, gvti, gnu,
     &       chii, chie, zkappa, btesla, gtau, omegaexb,
     &       zchiicyc, zchii1, zchii2, zchie1, zchie2,
     &       zrlt1, zrlt2
C
      EXTERNAL FEDG
C
      ierr=0
C
      do ii = 1, 32
         switches(ii) = 0
      end do
      switches(1)  = 0 ! 0: it produces no diagnostic output
      switches(2)  = 0 ! 0: it uses inputs ne19, tekev, tilev and btesla
                       ! 1: it uses inputs grhoi, gvti, gnu and gtau
      switches(3)  = 0 ! 0: the 1995 model, 1: the 1994 model
      switches(4)  = 1 ! 0: use gnu as given
                       ! 1: the definition in Dorland's IFS/PPPL routine
      switches(5)  = 1 ! 0: won't relax the restrictions on znu
                       ! 1: allows znu to be larger than 10.0
      switches(30) = 0
      switches(31) = 0
      switches(32) = 0
C
      DO NR=1,NRMAX-1
         znine  = SNGL( (RN(NR+1,2)+RN(NR,2))
     &                 /(RN(NR+1,1)+RN(NR,1)))
         IF(MDLUF.NE.0.AND.NSMAX.EQ.3) THEN
            zncne  = SNGL( (RN(NR+1,3)+RN(NR,3))
     &                    /(RN(NR+1,1)+RN(NR,1)))
         ELSE
            zncne  = 0.0
         ENDIF
         znbne  = SNGL( (RNFL(NR+1 )+RNFL(NR ))
     &                 /(RN (NR+1,1)+RN (NR,1)))
         zrlt   =-SNGL(RR/(0.5D0*(RT(NR+1,2)+RT(NR,2)))*
     &                           (RT(NR+1,2)-RT(NR,2))/DR*RJCB(NR))
         zrln   =-SNGL(RR/(0.5D0*(RN(NR+1,1)+RN(NR,1)))*
     &                           (RN(NR+1,1)-RN(NR,1))/DR*RJCB(NR))
         zq     = SNGL(QP(NR))
         zshat  = SNGL(S_AR(NR))
         zeps   = SNGL(EPSRHO(NR))
         zkappa = SNGL(RKPRHOG(NR))
         gnu    = SNGL((AME/AMM)*1.5625D-15*RN(NR,2)*1D20
     &                 /RT(NR,1)**1.5D0)
C         gnu    = 2.1*rmajor*ne19/(tekev**1.5 * tikev**0.5)
         gtau   = SNGL( (RT(NR+1,2)+RT(NR,2))
     &                 /(RT(NR+1,1)+RT(NR,1)))
C
         ne19   = SNGL(0.5D0*(RN(NR+1,1)+RN(NR,1))*1.D1)
         tekev  = SNGL(0.5D0*(RT(NR+1,1)+RT(NR,1)))
         tikev  = SNGL(0.5D0*(RT(NR+1,2)+RT(NR,2)))
         rmajor = SNGL(RR)
         btesla = SNGL(BB)
         grhoi  = SNGL(6.46D-3*SQRT(0.5D0*(RT(NR+1,2)+RT(NR,2)))/BB)
         gvti   = SNGL(2.19D5*SQRT(0.5D0*(RT(NR+1,2)+RT(NR,2))))
C
         CALL IFSPPPL( znine, zncne, znbne, zrlt, zrln,
     &                 zq, zshat, zeps, zkappa, omegaexb,
     &                 ne19, tekev, tikev, rmajor, btesla,
     &                 switches, grhoi, gvti, gnu, gtau,
     &                 chii, chie,
     &                 zchiicyc, zchii1, zchii2, zchie1, zchie2,
     &                 zrlt1, zrlt2, ierr )
C         IF(IERR.NE.0) THEN
C            WRITE(6,*) 'XX IFS/PPPL : ERROR IERR=',IERR
C            STOP
C         ENDIF
C
         AKDW(NR,1) = DBLE(chie)*AR1RHOG(NR)/AR2RHOG(NR)
         AKDW(NR,2) = DBLE(chii)*AR1RHOG(NR)/AR2RHOG(NR)
      ENDDO
C
      NR=NRMAX
         znine  = SNGL(PNSS(2)/PNSS(1))
         IF(MDLUF.NE.0.AND.NSMAX.EQ.3) THEN
            zncne  = SNGL(PNSS(3)/PNSS(1))
         ELSE
            zncne  = 0.0
         ENDIF
         znbne  = SNGL(RNFEDG)
         zrlt   =-SNGL(RR/PTS(2)*2.D0*(PTS (2)-RT(NR,2))/DR*RJCB(NR))
         zrln   =-SNGL(RR/PTS(1)*2.D0*(PNSS(1)-RN(NR,1))/DR*RJCB(NR))
         zq     = SNGL(QP(NR))
         zshat  = SNGL(S_AR(NR))
         zeps   = SNGL(EPSRHO(NR))
         zkappa = SNGL(RKPRHOG(NR))
         gnu    = SNGL((AME/AMM)*1.5625D-15*RN(NR,2)*1D20
     &                 /RT(NR,1)**1.5D0)
         gtau   = SNGL(PTS(2)/PTS(1))
C
         ne19   = SNGL(PNSS(1)*1.D1)
         tekev  = SNGL(PTS(1))
         tikev  = SNGL(PTS(2))
         rmajor = SNGL(RR)
         btesla = SNGL(BB)
         grhoi  = SNGL(6.46D-3*SQRT(PTS(2))/BB)
         gvti   = SNGL(2.19D5*SQRT(PTS(2)))
C
         CALL IFSPPPL( znine, zncne, znbne, zrlt, zrln,
     &                 zq, zshat, zeps, zkappa, omegaexb,
     &                 ne19, tekev, tikev, rmajor, btesla,
     &                 switches, grhoi, gvti, gnu, gtau,
     &                 chii, chie,
     &                 zchiicyc, zchii1, zchii2, zchie1, zchie2,
     &                 zrlt1, zrlt2, ierr )
C         IF(IERR.NE.0) THEN
C            WRITE(6,*) 'XX IFS/PPPL : ERROR IERR=',IERR
C            STOP
C         ENDIF
C
         AKDW(NR,1) = DBLE(chie)*AR1RHOG(NR)/AR2RHOG(NR)
         AKDW(NR,2) = DBLE(chii)*AR1RHOG(NR)/AR2RHOG(NR)
C
      RETURN
      END
