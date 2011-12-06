MODULE tritg

  PRIVATE
  PUBLIC tr_weiland,tr_ifspppl

CONTAINS

!     ***********************************************************

!            Weiland Model

!     ***********************************************************

  SUBROUTINE tr_weiland

!***********************************************************************
!  <INPUT>
!     ENL    : 2 Lne/Lb (Lb=R) (epsilon_ne)
!     EIL    : Lni/Lti (eta_i)
!     EEL    : Lne/Lte (eta_e)
!     TAUL   : Te/Ti (tau_i)
!     FLL    : Finite Larmor Radius parameter: (k_perp*rhos)**2=0.1
!     FTL    : trapped particle fraction
!     BQL    : fraction of an impurity
!     EQL    : Lnq/Ltq (q represents an impurity) (eta_q)
!     ENQL   : 2 Lnq/Lb (epsilon_nq)
!     ZL     : Z value for an impurity
!     BETAEL : electron beta
!     AZL    : impurity mass number
!     COLL   : factor multiplying collisionality (0 or 1)
!     ELL    : factor multiplying electromagnetic effects (0 or 1)
!     TEL    : electron temperature [keV]
!     TAUZL  : Te/Tq (tau_q)
!     RA     : minor radius (a)[m]
!     QL     : safety factor
!     SL     : magnetic shear
!     EPS    : local inverse aspect ratio (r/R)
!     EPSA   : inverse aspect ratio (a/R)
!     RNL    : electron density [10^19(20?) m^-3]
!     RLIST  : controlling printout in disp9t (0: off)
!     RNEQL  : number of equations (NDISP)
!     RKAP   : ellipticity
!     RIWL   : controlling printout in DIFFTD (0: off)
!     RISBL  : If ISB=1, we use the strong ballooning approximation
!              (GAV=1).
!     BB     : toroidal magnetic field [T]
!     SEARCH : The way the code chooses between ion and electron
!              eigenvalues for the use in the eigenfunction (WZ).
!              1 : Only eigenvalues with negative real part are used.
!                 (Real(WZ)<0)
!              2 : Eigenvalues with positive real part are used unless
!                  there are eigenvalues with negative real part.
!                 (Real(WZ)>0)
!              3 : The fastest growing mode with positive real part is
!                  used for eigenvalues with positive real part
!                  if there is such a root.
!     PMA    : mass number for an main ion
!     RGKL   : GRHO1/(a*GRHO2)
!     WEXBL  : ExB shearing rate
!     ROTL   : factor multiplying WEXBL (0 or 1)
!     NR     : radial grid number
!     IST    : 1   : iterations starting from an analytical
!                    approximation
!                    You should use only the first time step.
!              else: iterations starting from the value stored
!                    in WZJ(IK) which is the eigenvalue from the
!                    previous time step

!  <OUTPUT>
!     CHIL(5) : ion thermal transport coefficients for Ti, Te, Ne, Tq
!               and Nq equations [m^2/s (same as above)]
!     CHEL(5) : electron thermal transport coefficients
!     DL(5)   : ion particle transport coefficients
!     CHQL(5) : impurity thermal transport coefficients
!     DQL(5)  : impurity particle transport coefficients
!     SCHI    : effective ion thermal transport coefficient
!     SCHE    : effective electron thermal transport coefficient
!     SD      : effective ion particle transport coefficient
!     SCHQ    : effective impurity thermal transport coefficient
!     SDQ     : effective impurity particle transport coefficient

!***********************************************************************

    USE TRCOMM, ONLY : &
         AR1RHOG, AR2RHOG, BB, DR, EPSRHO, MDDW, MDLKAI, MDLTPF, NGLF, &
         NRMAX, NT, PA, PNSS, PTS, PZ, QP, RA, RHOG, RHOM, RJCB, RKAP, &
         RKEV, RMU0, RN, RR, RT, WEXB, S
    IMPLICIT NONE
    INTEGER(4):: ist, nr
    REAL(8)   :: &
           azl, betael, bql, coll, drl, eel, eil, ell, enl, enql, eps, &
           epsa, eql, fll, fls, ftl, ftpf, pma, ql, rgkl, risbl, riwl, &
           rlist, rneql, rnl, rotl, sche, schi, schq, sd, sdq, search, &
           shat, sl, slbl, slnel, slnil, slnql, sltel, sltil, sltql, taul, &
           tauzl, tel, wexbl, zl
      REAL(8)   :: deriv3p
      REAL(8),DIMENSION(5):: CHEL, CHIL, CHQL, DL, DQL

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
         SLNEL =-0.5D0*(RN(NR+1,1)+RN(NR,1))/((RN(NR+1,1)-RN(NR,1))*DRL)
         SLNIL =-0.5D0*(RN(NR+1,2)+RN(NR,2))/((RN(NR+1,2)-RN(NR,2))*DRL)
         SLNQL =-0.5D0*(RN(NR+1,3)+RN(NR,3))/((RN(NR+1,3)-RN(NR,3))*DRL)
         SLTEL =-0.5D0*(RT(NR+1,1)+RT(NR,1))/((RT(NR+1,1)-RT(NR,1))*DRL)
         SLTIL =-0.5D0*(RT(NR+1,2)+RT(NR,2))/((RT(NR+1,2)-RT(NR,2))*DRL)
         SLTQL =-0.5D0*(RT(NR+1,3)+RT(NR,3))/((RT(NR+1,3)-RT(NR,3))*DRL)
         SLBL  = RR
         ENL   = 2.D0*(SLNEL/SLBL )
         EIL   =       SLNIL/SLTIL
         EEL   =       SLNEL/SLTEL
         TAUL  = (RT(NR+1,1)+RT(NR,1))/(RT(NR+1,2)+RT(NR,2))
         FLL   = 1.D-1
         FTL   = FTPF(MDLTPF,EPS)
         BQL   = (RN(NR+1,3)+RN(NR,3))/(RN(NR+1,1)+RN(NR,1))
         EQL   =       SLNQL/SLTQL
         ENQL  = 2.D0*(SLNQL/SLBL )
         BETAEL= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))*RKEV*1.D20 &
                 /(BB**2/(2.D0*RMU0))
         TEL   = 0.5D0*(RT(NR+1,1)+RT(NR,1))
         TAUZL = (RT(NR+1,1)+RT(NR,1))/(RT(NR+1,3)+RT(NR,3))
         QL    = QP(NR)
         SL    = S(NR)
         RNL   = 0.5D0*(RN(NR+1,1)+RN(NR,1))*1.D1

         RGKL  = AR1RHOG(NR)/(RA*AR2RHOG(NR))
         WEXBL = WEXB(NR)
         IF(MDLKAI.EQ.64) THEN
            SHAT  = SQRT(2.D0*SL-1.D0+RKAP**2*(SL-1.D0)**2)
            FLS   = (0.7D0+2.4D0/(7.14D0*QL*SHAT+0.1D0))*FLL
            FLL   = 2.D0*FLS/(1.D0+1.D0/TAUL)
         ENDIF

!         COEF = PZ(2)**2*AEE**4*1.D20
!     &         /(6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)*RKEV**1.5D0)
!         RLAMB =15.2D0-0.5D0*DLOG(RNL)+DLOG(TEL)
!         VEI  = COEF*RNL*RLAMB/TEL**1.5D0
!         write(6,*) NR,COEF,VEI

         CALL TR_WEILAND_BRIDGE (ENL,EIL,EEL,TAUL,FLL,FTL,BQL,EQL,ENQL,ZL, &
              BETAEL,AZL,COLL,ELL,TEL,TAUZL,RA,QL,SL,EPS,EPSA,RNL,RLIST, &
              RNEQL,RKAP,RIWL,RISBL,BB,SEARCH,PMA,RGKL,WEXBL,ROTL,NR,IST, &
              CHIL,CHEL,DL,CHQL,DQL,SCHI,SCHE,SD,SCHQ,SDQ)

!         write(6,'(I3,5F15.7)') NR,CHIL(2),CHEL(2),DL(2),CHQL(2),DQL(2)

         CALL WEILAND_COEF(NR,CHIL,CHEL,DL,CHQL,DQL) !,SCHI,SCHE,SD,SCHQ,SDQ)
!       write(6,'(I3,5F15.7)') NR,CHEL(1),CHEL(2),CHEL(3),CHEL(4),CHEL(5)
!       write(6,'(I3,5F15.7)') NR,CHIL(1),CHIL(2),CHIL(3),CHIL(4),CHIL(5)
!       write(6,'(I3,5F15.7)') NR,CHQL(1),CHQL(2),CHQL(3),CHQL(4),CHQL(5)
!       write(6,'(I3,5F15.7)') NR,DL(1),DL(2),DL(3),DL(4),DL(5)
!         if(nr.eq.4) write(6,'(I4,5F15.7)') NT,DL(3),CHQL(4)
      ENDDO

      NR=NRMAX
         EPS   = EPSRHO(NR)
         SLNEL =-PNSS(1)/DERIV3P(PNSS(1),RN(NR,1),RN(NR-1,1), &
                                 RHOG(NR),RHOM(NR),RHOM(NR-1))
         SLNIL =-PNSS(2)/DERIV3P(PNSS(2),RN(NR,2),RN(NR-1,2), &
                                 RHOG(NR),RHOM(NR),RHOM(NR-1))
         SLNQL =-PNSS(3)/DERIV3P(PNSS(3),RN(NR,3),RN(NR-1,3), &
                                 RHOG(NR),RHOM(NR),RHOM(NR-1))
         SLTEL =-PTS (1)/DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1), &
                                 RHOG(NR),RHOM(NR),RHOM(NR-1))
         SLTIL =-PTS (2)/DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2), &
                                 RHOG(NR),RHOM(NR),RHOM(NR-1))
         SLTQL =-PTS (3)/DERIV3P(PTS(3),RT(NR,3),RT(NR-1,3), &
                                 RHOG(NR),RHOM(NR),RHOM(NR-1))
         SLBL  = RR
         ENL   = 2.D0*(SLNEL/SLBL )
         EIL   =       SLNIL/SLTIL
         EEL   =       SLNEL/SLTEL
         TAUL  = PTS(1)/PTS(2)
         FLL   = 1.D-1
         FTL   = FTPF(MDLTPF,EPS)
         BQL   = PNSS(3)/PNSS(1)
         EQL   =       SLNQL/SLTQL
         ENQL  = 2.D0*(SLNQL/SLBL )
         BETAEL= PNSS(1)*PTS(1)*RKEV*1.D20/(BB**2/(2.D0*RMU0))
         TEL   = PTS(1)
         TAUZL = PTS(1)/PTS(3)
         QL    = QP(NR)
         SL    = S(NR)
         RNL   = PNSS(1)*1.D1

         RGKL  = AR1RHOG(NR)/(RA*AR2RHOG(NR))
         WEXBL = WEXB(NR)
         IF(MDLKAI.EQ.64) THEN
            SHAT  = SQRT(2.D0*SL-1.D0+RKAP**2*(SL-1.D0)**2)
            FLS   = (0.7D0+2.4D0/(7.14D0*QL*SHAT+0.1D0))*FLL
            FLL   = 2.D0*FLS/(1.D0+1.D0/TAUL)
         ENDIF

         CALL TR_WEILAND_BRIDGE (ENL,EIL,EEL,TAUL,FLL,FTL,BQL,EQL,ENQL,ZL, &
              BETAEL,AZL,COLL,ELL,TEL,TAUZL,RA,QL,SL,EPS,EPSA,RNL,RLIST, &
              RNEQL,RKAP,RIWL,RISBL,BB,SEARCH,PMA,RGKL,WEXBL,ROTL,NR,IST, &
              CHIL,CHEL,DL,CHQL,DQL,SCHI,SCHE,SD,SCHQ,SDQ)

         CALL WEILAND_COEF(NR,CHIL,CHEL,DL,CHQL,DQL) !,SCHI,SCHE,SD,SCHQ,SDQ)

      RETURN
    END SUBROUTINE tr_weiland

!     ****************************************************************

    SUBROUTINE WEILAND_COEF(NR,CHIL,CHEL,DL,CHQL,DQL) !,SCHI,SCHE,SD,SCHQ,SDQ)

      USE TRCOMM, ONLY : &
           ADDW, ADDWD, ADDWP, AKDW, AKDWD, AKDWP, MDLWLD, NGLF, NSM, PA
      IMPLICIT NONE
      INTEGER(4)  :: NR
!      REAL(8)     :: SCHE, SCHI, SCHQ, SD, SDQ
      REAL(8),DIMENSION(5):: CHEL, CHIL, CHQL, DL, DQL
      INTEGER(4):: ns, ns1


!     The diagonal value of the transport coefficient matrix
!     is set to be zero if it becomes negative.

      IF(DL  (3).LT.0.D0) DL  (3)=0.D0
      IF(CHEL(2).LT.0.D0) CHEL(2)=0.D0
      IF(CHIL(1).LT.0.D0) CHIL(1)=0.D0
      IF(DQL (5).LT.0.D0) DQL (5)=0.D0
      IF(CHQL(4).LT.0.D0) CHQL(4)=0.D0

!     It is assumed that De=Di in the followings.

      IF(PA(3).EQ.3.D0) THEN
!         AKDW(NR,1)=SCHE
         AKDW(NR,1)=CHEL(2)
         DO NS=2,NSM
!            AKDW(NR,NS)=SCHI
            AKDW(NR,NS)=CHIL(1)
         ENDDO
         DO NS=1,NSM
!            ADDW(NR,NS)=SD
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
!         ADDW(NR,2)=0.D0
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
      RETURN
    END SUBROUTINE WEILAND_COEF

!     ****************************************************************

    SUBROUTINE TR_WEILAND_BRIDGE (EN,EI,EE,TAU,FL,FT,BQ,EQ,ENQ,Z,BETAE,AZ, &
           COL,EL,TE,TAUZ,PR,Q,S,EPSR,E,N,LIST,RNEQ,KAPPA,RIW,RISB,BTOR, &
           SEARCH,MA,GKL,WEXBL,ROTL,IKL,ISTL,CHI,CHE,D,CHQ,DQ,SCHIL,SCHEL, &
           SDL,SCHQL,SDQL)

      IMPLICIT NONE
      REAL(8),   INTENT(IN):: &
           TE,BTOR,RIW,RISB,EI,TAU,FL,EN,FT,EE,BQ,EQ,ENQ,Z,TAUZ,BETAE,MA,E, &
           Q,S,KAPPA,PR,N,RNEQ,EPSR,AZ,COL,EL,LIST,SEARCH,GKL,WEXBL,ROTL
      INTEGER(4),INTENT(IN):: IKL,ISTL
      REAL(8),DIMENSION(5),INTENT(INOUT):: CHI,CHE,D,CHQ,DQ
      REAL(8),INTENT(OUT)  :: SCHIL,SCHEL,SDL,SCHQL,SDQL
      INTEGER(4):: ICP, IR, NDISP
      REAL(8)   :: &
           A, B, C, DS, DTOT, ENN, G, H, LAMB, LNH, NI, NQ, R, TAUI, TVR, &
           VEI, WI, WIS, WR, WRS, WST, ZEFF, ZX
      REAL(8),DIMENSION(5,100):: U
      COMPLEX(8):: W
      COMPLEX(8),DIMENSION(10):: RP

      INTEGER(4):: &
           I,IW,IX,IK,IST,ITS,ITL,ITERA,ITC,ISB,LPRINTIN,NDIM,NEQ,IRET, &
           SEARCHMODE
      REAL(8)   :: &
           THRD,STR,XIH,RFL,H1,FTR,GQ,BF,ZE,CS,KPPA,RAV,ALP,ALF,KPC,D1,SI, &
           KIQ,KXQ,WDE,EPS,SCHI,SCHE,SD,SCHQ,SDQ,EA,GK,ETE,ETI,ETQ,AZL,ALAF, &
           GAV,VEF,BTA,EM,EIC,EEC,ENC,TAUC,FLC,FTC,EQC,ENQC,BETAEC,TAUZC,QC, &
           SC,ENHC,ZFS,KPS,CHIC,WEXB,ROT,WZIMAX,TOL,SHPE,SCHEF,DEF,LTH,LTE, &
           LN,LTQ,LNQ
      REAL(8),DIMENSION(5):: HPT
      REAL(8),DIMENSION(32):: CETAIN
      REAL(8),DIMENSION(10,10):: ZVR,ZVI
      COMPLEX(8):: HQ,WZ,WZP
      COMPLEX(8),DIMENSION(10):: ZZ
      COMMON/IRET/ IRET
      COMMON/PAR/ THRD,FTR,STR,RFL,D1,XIH,IX,IW
      COMMON/EFFDIFF/ SCHI,SCHE,SD,SCHQ,SDQ
!      COMMON/EFF DIFF/ SCHI,SCHE,SD,SCHQ,SDQ
      COMMON/GRAD/ ENC,ENHC,EIC,EEC
      COMMON/COLL/ BTA,VEF
      COMMON/IMP/ BF,GQ,SI,ZE,TAUZC,ENQC,EQC,KIQ,KXQ,ZFS
      COMMON/ETAWN6/ CETAIN,ETE,ETI,ETQ,TAUC,AZL,FTC,FLC,LPRINTIN,NDIM
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


!     ICP : controlling printout in TR_WEILAND_BRIDGE (0: off)
      ICP=0
      IF(ICP.NE.0) THEN
         WRITE(6,*) '//////////////////////////////////////////////////'
      ENDIF

      GK   = GKL     ! GRHO1/(a*GRHO2)
      WEXB = WEXBL   ! EXB shearing rate according to Hahm and Burrell
      ROT  = ROTL    ! factor that multiplies the EXB shearing rate
      IK   = IKL     ! the space profile index for WZL(IK)
      IST  = ISTL    ! 1:analytical form used, other:previous value used

      ZFS = 0.D0     ! the product of charge and fraction (to Ne) of fast particles

      THRD=1.D0/3.D0 ! fixed coefficient (THiRD)
      TVR=2.D0*THRD  ! fixed coefficient (Two thiRd)
      FTR=5.D0/3.D0  ! fixed coefficient (Five ThiRd)
      STR=7.D0/3.D0  ! fixed coefficient (Seven ThiRd)
      BTA=1.5D0      ! fixed, parameter used in the collision model

      EPS=EPSR       ! local inverse aspect ratio
      SEARCHMODE=INT(SEARCH) ! the way to search eigenvalues
      AZL=AZ         ! Aq (atomic number of impurity)
      EM=EL          ! factor that multiplies electromagnetic effects
      ZE=Z           ! Zq
      BF=BQ          ! nq/ne

      G=1.D0-Z*BQ    ! transforming factor from ne to ni
      ENN=1.D0-Z*BQ*EN/ENQ ! (Lne/Lni)*(ni/ne)
      IF(ABS(ENN).GE.0.001D0) GO TO 89
      IF(ENN.LT.0.D0) ENN=-0.001D0
      IF(ENN.GT.0.D0) ENN= 0.001D0
   89 CONTINUE
      ENC=EN         ! 2*Lne/R
      ENHC=G*EN/ENN  ! 2*Lni/R
      ENQC=ENQ       ! 2*Lnq/R
!      ENQC=EN        ! 2*Lne/R (original definition)

      EIC=EI         ! Lni/Lti
      EEC=EE         ! Lne/Lte
      EQC=EQ         ! Lnq/Ltq
      TAUC=TAU       ! Te/Ti
      TAUZC=TAUZ     ! Te/Tq

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

      IX=IK          ! the space profile index for U(J,IX)

      DO I=1,31
         CETAIN(I)=0.D0 ! control vector
      ENDDO
!      CETAIN(32)=1.D-15
      CETAIN(32)=0.001  ! accuracy parameter in the NAG routine
      LPRINTIN=INT(LIST)

      TAUI=1.D0/TAU  ! Ti/Te
      RFL=SQRT(FL)   ! FLR parameter (same with FLC)
      GQ=1.D0-Z*BQ   ! transforming factor from ne to ni (same with G)
      NQ=BQ*N        ! nq
      NI=GQ*N        ! ni
      ZEFF=(NI+Z*Z*NQ)/N ! Zeff

      LN =0.5D0*R*ENC/PR  ! Lne/a
      LNH=0.5D0*R*ENHC/PR ! Lni/a
      LNQ=0.5D0*R*ENQC/PR ! Lnq/a
      LTE=LN/EEC          ! Lte/a
      LTH=LNH/EIC         ! Lti/a
      LTQ=LNQ/EQC         ! Ltq/a

      ETI=ENHC/EIC   ! 2*Lti/R
      ETE=ENC/EEC    ! 2*Lte/R
      ETQ=ENQC/EQC   ! 2*Ltq/R

      KIQ=ENQC/ENHC  ! Lnq/Lni
      KXQ=ENC/ENQC   ! Lne/Lnq

      CS=3.095D5*DSQRT(TE/MA) ! sound speed
      IF(ICP.NE.0) THEN
         WRITE(*,00126) EN,EI,EE,FL,TAU
00126    FORMAT(2X,'EN=',F8.3,' EI=',F8.3,' EE=',F8.3,' FL=',F8.3, ' TAU=',G12.4)
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

      WST=DSQRT(FL)*CS/(PR*ABS(LN))  ! diamagnetic frequency (omega_star)
      WDE=ABS(EN)*WST           ! curvature drift frequency (omega_drift_electron)
      LAMB=15.95D0-DLOG(DSQRT(N)/TE) ! coulomb logarithm (lambda)
      VEI=9.19D2*NI*LAMB/TE**1.5D0   ! collisionality (nu_electron_ion)
      VEF=COL*VEI/(EPS*WDE)     ! effective electron ion collision frequency for trapped electrons, normalized by omega_de
      WEXB=WEXB/WDE ! ExB shearing rate should be normalized with WDE

      IF(ICP.NE.0) THEN
         WRITE(*,00130) WST,VEF,VEI,ZEFF,COL
00130    FORMAT(2X,'WST=',G11.3,' VEF=',G11.3,' VEI=',G11.3, ' ZEFF=',G11.3,' COL=',G11.3)

         H=0.5D0*ABS(S)/q       ! not used
         WRITE(*,00133) H
00133    FORMAT(2X,'H=',G12.4)
      ENDIF
!   -----------------------------------------------

      IF(ICP.NE.0) THEN
         A=1.D0-EN*(1.D0+10.D0/(3.D0*TAU))-FL*(1.D0+EI+5.D0*EN/3.D0)/TAU
         B=EI-7.D0/3.D0+5.D0*EN*(1.D0+1.D0/TAU)/3.D0 +5.D0*FL*(1.D0+EI)/(3.D0*TAU)
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

      ITC=1      ! 1:iteration, other:previous eigenvalues used
      ITL=80     ! Maximum number of iterations
      TOL=0.01D0 ! relative error for convergence
!      IST=1

!      IF(IX.GE.48) write(6,'(I3,2F15.7)') IX,HQ
      CALL disp9t(NDISP,ZZ)

      IF(ICP.NE.0) THEN
         WRITE(*,174) ISB
 174     FORMAT(' ISB=',I5)
         WRITE(*,175) ITC,ITS,ITERA
 175     FORMAT('  ITC=',I5,' ITS=',I5,' ITER=',I5)
         WRITE(*,177) WZ,WZP
 177     FORMAT('  WZ=',2G11.3,' WZP=',2G11.3)
      ENDIF

      IR=0  ! the number of unstable roots (given by the following)

      DO I=1,NDISP
         ZX=DIMAG(ZZ(I))
         IF(ZX.GT. 0.001) THEN
            IR=IR+1
            RP(IR)=ZZ(I) ! unstable roots found by disp9t
         ENDIF
      END DO
      IF(ICP.NE.0) THEN
         WRITE(*,310) IR
 310     FORMAT(2X,' IR=',I4)
         WRITE(*,00134) ALAF
      ENDIF

      DO I=1,IR
         W=RP(I)
         WR=EN*DREAL(W) ! real part of unstable roots
         WI=EN*DIMAG(W) ! imaginary part of unstable roots
         IF(ICP.NE.0) WRITE(*,311) WR,WI,I
         WRS=WDE*DREAL(W)
         WIS=WDE*DIMAG(W)
         IF(ICP.NE.0) WRITE(*,321) WRS,WIS
      END DO
  311 FORMAT(//,2X,'WR=',G11.3,' WI=',G11.3,' I=',I5)
  321 FORMAT(' WRS=',G11.3,' WIS=',G11.3)

      WZ=EN*WZ
      HQ=EN*HQ
      IF(ICP.NE.0) THEN
         WRITE(*,00128) ALF,ALP,WZ,KAPPA
00128    FORMAT(/,2X,' ALF=',G11.3,' ALP=',G11.3,' WZ=',2G11.3,/,' KAPPA=',G11.3)
         WRITE(*,312) H1,HQ,GAV,RAV
 312     FORMAT(//,2X,'H1=',G12.4,' HQ=',2G12.4,' GAV=',G12.4,/,' RAV=',G12.4)
      ENDIF
      U(1,IX)=TAUI*TE
      U(2,IX)=TE
      U(3,IX)=N
      U(4,IX)=TE/TAUZ
      U(5,IX)=BQ*N

      CALL DIFF(RP,IR,TAUI,FT,U,CHI,CHE,D,CHQ,DQ)

      IF(ICP.NE.0) THEN
         WRITE(*,330) SCHI,SCHE,SD,SCHQ,SDQ
 330     FORMAT(/,2X,'CHIEFF=',G11.3,' CHEEFF=',G11.3,' DEFF=',G11.3, ' CHQEFF=',G11.3,' DQEFF=',G11.3)

!         WRITE(*,331) CHE(2),SCHEF
! 331     FORMAT(' CHE(2)=',G11.3,' SCHEF=',G11.3)
         WRITE(*,331) CHQ(5),SCHQ
 331     FORMAT(' CHQ(5)=',G11.3,' SCHQ=',G11.3)
         DTOT=SD+DEF
         WRITE(*,332) D(3),DEF,DTOT
 332     FORMAT(' D(3)=',G11.3,' DEF=',G11.3,' DTOT=',G11.3)
         WRITE(*,335) XIH,SHPE,CHIC,SCHEF
 335     FORMAT('  XIH=',G11.3,' SHPE=',G11.3,' CHIC=',G11.3,' SCHEF=',G11.3)
      ENDIF
!00150 CONTINUE
      SCHIL = SCHI
      SCHEL = SCHE
      SDL   = SD
      SCHQL = SCHQ
      SDQL  = SDQ

      RETURN
    END SUBROUTINE TR_WEILAND_BRIDGE
!     ***********************************************************

!            IFS/PPPL Model

!     ***********************************************************

    SUBROUTINE tr_ifspppls(NSTM,NRMAX,RN,RR,DR,RJCB,RHOG,RHOM,QP, &
           S,EPSRHO,RKPRHOG,RT,BB,AMM,AME,PNSS,PTS,RNFL,RBEEDG,MDLUF,NSMAX, &
           AR1RHOG,AR2RHOG,AKDW)

      IMPLICIT NONE

      INTEGER(4),INTENT(IN):: NSTM,NRMAX,MDLUF,NSMAX
      REAL(8)   ,INTENT(IN):: RR,DR,BB,AMM,AME,RBEEDG
      REAL(8),DIMENSION(NRMAX,NSMAX),INTENT(IN):: RN, RT
      REAL(8),DIMENSION(NRMAX),INTENT(IN):: &
           RJCB,RHOG,RHOM,QP,S,EPSRHO,RKPRHOG,RNFL,AR1RHOG,AR2RHOG
      REAL(8),DIMENSION(NSMAX)    ,INTENT(IN):: PNSS,PTS
      REAL(8),DIMENSION(NRMAX,NSTM),INTENT(OUT)::AKDW
      integer(4) :: ii, ierr, nr
      integer(4),parameter :: ipin=7,iptmp=8,ipout=9,screen=6
      integer(4),dimension(32):: switches
      real(4) :: &
           btesla, chie, chii, gnu, grhoi, gtau, gvti, ne19, omegaexb, &
           rmajor, tekev, tikev, zchie1, zchie2, zchii1, zchii2, zchiicyc, &
           zeps, zkappa, znbne, zncne, znine, zq, zrln, zrlt, zrlt1, zrlt2, &
           zshat
      REAL(8)   :: DERIV3P


      ierr=0

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

      DO NR=1,NRMAX-1
         znine  = SNGL( (RN(NR+1,2)+RN(NR,2)) /(RN(NR+1,1)+RN(NR,1)))
         IF(MDLUF.NE.0.AND.NSMAX.EQ.3) THEN
            zncne  = SNGL( (RN(NR+1,3)+RN(NR,3))/(RN(NR+1,1)+RN(NR,1)))
         ELSE
            zncne  = 0.0
         ENDIF
         znbne  = SNGL( (RNFL(NR+1 )+RNFL(NR ))/(RN (NR+1,1)+RN (NR,1)))
         zrlt   =-SNGL(RR/(0.5D0*(RT(NR+1,2)+RT(NR,2))) &
                           *(RT(NR+1,2)-RT(NR,2))/DR*RJCB(NR))
         zrln   =-SNGL(RR/(0.5D0*(RN(NR+1,1)+RN(NR,1))) &
                           *(RN(NR+1,1)-RN(NR,1))/DR*RJCB(NR))
         zq     = SNGL(QP(NR))
         zshat  = SNGL(S(NR))
         zeps   = SNGL(EPSRHO(NR))
         zkappa = SNGL(RKPRHOG(NR))
         gnu    = SNGL((AME/AMM)*1.5625D-15*RN(NR,2)*1D20/RT(NR,1)**1.5D0)
!         gnu    = 2.1*rmajor*ne19/(tekev**1.5 * tikev**0.5)
         gtau   = SNGL( (RT(NR+1,2)+RT(NR,2))/(RT(NR+1,1)+RT(NR,1)))

         ne19   = SNGL(0.5D0*(RN(NR+1,1)+RN(NR,1))*1.D1)
         tekev  = SNGL(0.5D0*(RT(NR+1,1)+RT(NR,1)))
         tikev  = SNGL(0.5D0*(RT(NR+1,2)+RT(NR,2)))
         rmajor = SNGL(RR)
         btesla = SNGL(BB)
         grhoi  = SNGL(6.46D-3*SQRT(0.5D0*(RT(NR+1,2)+RT(NR,2)))/BB)
         gvti   = SNGL(2.19D5*SQRT(0.5D0*(RT(NR+1,2)+RT(NR,2))))

         CALL IFSPPPL( znine, zncne, znbne, zrlt, zrln, zq, zshat, zeps, &
              zkappa, omegaexb, ne19, tekev, tikev, rmajor, btesla, &
              switches, grhoi, gvti, gnu, gtau, chii, chie, zchiicyc, &
              zchii1, zchii2, zchie1, zchie2, zrlt1, zrlt2, ierr )
!         IF(IERR.NE.0) THEN
!            WRITE(6,*) 'XX IFS/PPPL : ERROR IERR=',IERR
!            STOP
!         ENDIF

         AKDW(NR,1) = DBLE(chie)*AR1RHOG(NR)/AR2RHOG(NR)
         AKDW(NR,2) = DBLE(chii)*AR1RHOG(NR)/AR2RHOG(NR)
      ENDDO

      NR=NRMAX
         znine  = SNGL(PNSS(2)/PNSS(1))
         IF(MDLUF.NE.0.AND.NSMAX.EQ.3) THEN
            zncne  = SNGL(PNSS(3)/PNSS(1))
         ELSE
            zncne  = 0.0
         ENDIF
         znbne  = SNGL(RBEEDG)
         zrlt   =-SNGL(RR/PTS(2)*DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2), &
                                         RHOG(NR),RHOM(NR),RHOM(NR-1)))
         zrln   =-SNGL(RR/PTS(1)*DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1), &
                                         RHOG(NR),RHOM(NR),RHOM(NR-1)))
         zq     = SNGL(QP(NR))
         zshat  = SNGL(S(NR))
         zeps   = SNGL(EPSRHO(NR))
         zkappa = SNGL(RKPRHOG(NR))
         gnu    = SNGL((AME/AMM)*1.5625D-15*RN(NR,2)*1D20/RT(NR,1)**1.5D0)
         gtau   = SNGL(PTS(2)/PTS(1))

         ne19   = SNGL(PNSS(1)*1.D1)
         tekev  = SNGL(PTS(1))
         tikev  = SNGL(PTS(2))
         rmajor = SNGL(RR)
         btesla = SNGL(BB)
         grhoi  = SNGL(6.46D-3*SQRT(PTS(2))/BB)
         gvti   = SNGL(2.19D5*SQRT(PTS(2)))

         CALL IFSPPPL( znine, zncne, znbne, zrlt, zrln,zq, zshat, zeps, &
              zkappa, omegaexb, ne19, tekev, tikev, rmajor, btesla, &
              switches, grhoi, gvti, gnu, gtau,chii, chie, zchiicyc, &
              zchii1, zchii2, zchie1, zchie2,zrlt1, zrlt2, ierr )
!         IF(IERR.NE.0) THEN
!            WRITE(6,*) 'XX IFS/PPPL : ERROR IERR=',IERR
!            STOP
!         ENDIF

         AKDW(NR,1) = DBLE(chie)*AR1RHOG(NR)/AR2RHOG(NR)
         AKDW(NR,2) = DBLE(chii)*AR1RHOG(NR)/AR2RHOG(NR)
    END SUBROUTINE tr_ifspppls
  END MODULE tritg
