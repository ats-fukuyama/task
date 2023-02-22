!
!     ***********************************************************

!           CALCULATE TRANSPORT COEFFICIENTS AND SOURCE

!     ***********************************************************

      SUBROUTINE TRCALC(IERR)

      USE TRCOMM, ONLY : AJBS, AJRF, AJRFV, AR1RHOG, ARRHOG, &
           BP, DT, MDLEQ0, &
           MDLJBS, MDLNF, MDLUF, MDLPR, MDNCLS, NRAMAX, NRMAX, &
           NROMAX, NSM, NSMAX, PBCL, PBIN, PCX, PELTIM, PEX, PFCL, PFIN, &
           PI, PIE, PIN, PN, PNB, PNF, POH, PRB, PRC, PRF, PRFV, PRL, PRSUM, &
           Q0, QP, RDP, RG, RHOA, RR, SCX, SEX, SIE, SNB, SNF, SPE, SSIN, &
           T, TAUF, TTRHOG, RDPVRHOG, SPSC, &
           pellet_time_start,pellet_time_interval, &
           number_of_pellet_repeat,icount_of_pellet
      USE tr_cytran_mod
      USE libitp
      IMPLICIT NONE
      INTEGER,INTENT(OUT)    :: IERR
      INTEGER                :: NR,NS
      REAL(rkind),SAVE:: pellet_time_start_save=-1.D0
      REAL(rkind):: t_pellet

      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      IERR=0

      SIE(1:NRMAX)=0.D0
      SNF(1:NRMAX)=0.D0
      SNB(1:NRMAX)=0.D0
      POH(1:NRMAX)=0.D0
      PIE(1:NRMAX)=0.D0
      PCX(1:NRMAX)=0.D0
      PRB(1:NRMAX)=0.D0
      PRC(1:NRMAX)=0.D0
      PRL(1:NRMAX)=0.D0
      PRSUM(1:NRMAX)=0.D0
      PNB(1:NRMAX)=0.D0
      PNF(1:NRMAX)=0.D0
      PBIN(1:NRMAX)=0.D0
      PFIN(1:NRMAX)=0.D0
!      AJNB(1:NRMAX)=0.D0
      AJRFV(1:NRMAX,1)=0.D0
      AJRFV(1:NRMAX,2)=0.D0
      AJRFV(1:NRMAX,3)=0.D0
      AJRF(1:NRMAX)=0.D0
      AJBS(1:NRMAX)=0.D0
      SPSC(1:NRMAX,1:NSMAX)=0.D0
      SPE(1:NRMAX,1:NSMAX)=0.D0
      PBCL(1:NRMAX,1:NSMAX)=0.D0
      PFCL(1:NRMAX,1:NSMAX)=0.D0
      IF(MDLUF.NE.0) THEN
         PRFV(1:NRMAX,3:NSM,1)=0.D0
         PRFV(1:NRMAX,3:NSM,2)=0.D0
         PRFV(1:NRMAX,3:NSM,3)=0.D0
         PRF(1:NRMAX,3:NSM)=0.D0
      ELSE
         PRFV(1:NRMAX,1:NSM,1)=0.D0
         PRFV(1:NRMAX,1:NSM,2)=0.D0
         PRFV(1:NRMAX,1:NSM,3)=0.D0
         PRF(1:NRMAX,1:NSM)=0.D0
      ENDIF

      BP(1:NRMAX)=AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR
      QP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*RDPVRHOG(1:NRMAX))
      Q0=FCTR(RG(1),RG(2),QP(1),QP(2))

!     *** RADIAL ELECTRIC FIELD ***

      CALL TRERAD

      IF(pellet_time_start.NE.pellet_time_start_save) THEN
         icount_of_pellet=0
         pellet_time_start_save=pellet_time_start
!         WRITE(6,'(A,I6,2ES12.4)') &
!              '# pellet setup:',icount_of_pellet,pellet_time_start,t
      END IF
      IF(T.GT.pellet_time_start-0.5D0*DT) THEN
!         WRITE(6,*) '@@@',t,pellet_time_start+0.5D0*DT
         IF(icount_of_pellet.LT.number_of_pellet_repeat) THEN
            t_pellet=pellet_time_start+icount_of_pellet*pellet_time_interval
!            WRITE(6,'(A,3ES12.4)') '# t_pellet:', &
!                 T,t_pellet-0.5D0*DT,t_pellet+0.5D0*DT
            IF(T.GT.t_pellet-0.5D0*DT.AND.T.LE.t_pellet+0.5D0*DT) THEN
               CALL TRPELT
               icount_of_pellet=icount_of_pellet+1
               WRITE(6,'(A,I6,2ES12.4)') &
                    '# pellet injected:',icount_of_pellet,T,t_pellet
            ENDIF
         ENDIF
      END IF
            
      IF(T.LT.PELTIM+0.5D0*DT.AND. T.GE.PELTIM-0.5D0*DT) CALL TRPELT

      CALL TRPSC

      CALL TRZEFF

      IF(MDLPR.GT.0) CALL TR_CYTRAN

      IF(MDNCLS.NE.0) THEN
         CALL TR_NCLASS(IERR)
         IF(IERR.NE.0) RETURN
      ENDIF

      CALL TRCOEF
      CALL TRLOSS
      IF(MDLUF.NE.1.AND.MDLUF.NE.3) CALL TRPWRF
      CALL TRPWNB

      IF(MDNCLS.NE.0) THEN
         CALL TRAJBS_NCLASS
      ELSE
         select case(MDLJBS)
         case(1)
            CALL TRAJBS
         case(2)
            CALL TRAJBS
         case(3)
            CALL TRAJBS
         case(4)
            CALL TRAJBSNEW
         case(5)
            CALL TRAJBSSAUTER
         case default
            CALL TRAJBS
         end select
      ENDIF

      SELECT CASE(MDLNF)
      CASE(0)
         TAUF(1:NRMAX)=1.D0
      CASE(1:4)
         CALL TRNFDT
      CASE(5:6)
         CALL TRNFDHE3
      END SELECT
      
      CALL TRAJOH

      DO NR=1,NRMAX
         IF(MDLEQ0.EQ.0) THEN
            DO NS=1,NSMAX
               SELECT CASE(NS)
               CASE(1)
                  SSIN(NR,1)= SIE(NR) &
                                     +SNB(NR)+SEX(NR,1)+SPSC(NR,1)
               CASE(2)
                  SSIN(NR,2)= PN(2)*SIE(NR)/(PN(2)+PN(3)) &
                             -SNF(NR)+SNB(NR)+SEX(NR,2)+SPSC(NR,2)
               CASE(3)
                  SSIN(NR,3)= PN(3)*SIE(NR)/(PN(2)+PN(3)) &
                             -SNF(NR)        +SEX(NR,3)+SPSC(NR,3)
               CASE(4)
                  SSIN(NR,4)= SNF(NR)        +SEX(NR,4)+SPSC(NR,4)
               CASE(7)
                  SSIN(NR,7)=-SIE(NR)        -SCX(NR)
               CASE(8)
                  SSIN(NR,8)=         SNB(NR)+SCX(NR)
               END SELECT
            END DO
         ELSEIF(MDLEQ0.EQ.1) THEN
            DO NS=1,NSMAX
               SELECT CASE(NS)
               CASE(1)
                  SSIN(NR,1)=         SNB(NR)+SEX(NR,1)+SPSC(NR,1)
               CASE(2)
                  SSIN(NR,2)=-SNF(NR)+SNB(NR)+SEX(NR,2)+SPSC(NR,2)
               CASE(3)
                  SSIN(NR,3)=-SNF(NR)        +SEX(NR,3)+SPSC(NR,3)
               CASE(4)
                  SSIN(NR,4)= SNF(NR)        +SEX(NR,4)+SPSC(NR,4)
               CASE(7)
                  SSIN(NR,7)=0.D0
               CASE(8)
                  SSIN(NR,8)=         SNB(NR)
               END SELECT
            END DO
         ENDIF
         PIN(NR,1)=PBCL(NR,1)+PFCL(NR,1)+PRF(NR,1) &
              &   +POH(NR)-PRSUM(NR)-PIE(NR)+PEX(NR,1)
         PIN(NR,2)=PBCL(NR,2)+PFCL(NR,2)+PRF(NR,2) &
              &   -PN(2)*PCX(NR)/(PN(2)+PN(3))+PEX(NR,2)
         PIN(NR,3)=PBCL(NR,3)+PFCL(NR,3)+PRF(NR,3) &
              &   -PN(3)*PCX(NR)/(PN(2)+PN(3))+PEX(NR,3)
         PIN(NR,4)=PBCL(NR,4)+PFCL(NR,4)+PRF(NR,4)+PEX(NR,4)
      ENDDO

      IF(RHOA.NE.1.D0) NRMAX=NRAMAX

      RETURN
      END SUBROUTINE TRCALC

!     ***********************************************************

!           RADIAL ELECTRIC FIELD

!     **********************************************************

      SUBROUTINE TRERAD

      USE TRCOMM, ONLY : AEE, AMM, BB, BP, DR, EPSRHO, ER, MDLER, NRMAX, &
           & PA, PADD, PBM, PNSS, PTS, PZ, QP, RHOG, RHOM, &
           & RJCB, RKEV, RN, RNF, RT, SUMPBM, VPOL, VTOR, rkind
      USE libitp
      IMPLICIT NONE
      INTEGER:: NR
      REAL(rkind)   :: ALPHA_NEO, CS, DPD, DRL, EPS, F_UNTRAP, RHO_S, RLNI, &
           & RLTI, TEL, TERM_DP, TIL

      IF(SUMPBM.EQ.0.D0) THEN
         PADD(1:NRMAX)=0.D0
      ELSE
         PADD(1:NRMAX)=PBM(1:NRMAX)*1.D-20/RKEV-RNF(1:NRMAX,1)*RT(1:NRMAX,2)
      ENDIF
      DO NR=1,NRMAX
         DRL=RJCB(NR)/DR
         IF(NR.EQ.NRMAX) THEN
            DPD = DERIV3P(PNSS(2)*PTS(2),RN(NR  ,2)*RT(NR  ,2)-PADD(NR  ), &
     &                    RN(NR-1,2)*RT(NR-1,2)-PADD(NR-1),RHOG(NR),RHOM(NR),RHOM(NR-1))
            TERM_DP = DPD*RKEV/(PZ(2)*AEE*PNSS(2))
         ELSE
            DPD =(  RN(NR+1,2)*RT(NR+1,2)-PADD(NR+1)-(RN(NR  ,2)*RT(NR  ,2)-PADD(NR  )))*DRL
            TERM_DP = DPD*RKEV/(PZ(2)*AEE*0.5D0*(RN(NR+1,2)+RN(NR,2)))
         ENDIF
         IF(MDLER.EQ.0) THEN
!     pressure gradient only
            ER(NR) = TERM_DP
         ELSEIF(MDLER.EQ.1) THEN
!     nabla p + toroidal rotation
            ER(NR) = TERM_DP+VTOR(NR)*BP(NR)
         ELSEIF(MDLER.EQ.2) THEN
!     nabla p + V_tor + poloidal rotation
            ER(NR) = TERM_DP+VTOR(NR)*BP(NR)-VPOL(NR)*BB
         ELSEIF(MDLER.EQ.3) THEN
!     Waltz definition
            EPS = EPSRHO(NR)
            F_UNTRAP = 1.D0-1.46D0*SQRT(EPS)+0.46D0*EPS**1.5D0
            ALPHA_NEO = 1.D0-0.8839D0*F_UNTRAP/(0.3477D0+0.4058D0*F_UNTRAP)
            IF(NR.EQ.NRMAX) THEN
               TEL = PTS(1)
               TIL = PTS(2)
               RLNI = -DERIV3P(PNSS(2),RN(NR,2),RN(NR-1,2),RHOG(NR),RHOM(NR),RHOM(NR-1))/PNSS(2)
               RLTI = -DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2),RHOG(NR),RHOM(NR),RHOM(NR-1))/PTS(2)
            ELSE
               TEL = 0.5D0*(RT(NR,1)+RT(NR+1,1))
               TIL = 0.5D0*(RT(NR,2)+RT(NR+1,2))
               RLNI = -(LOG(RN(NR+1,2))-LOG(RN(NR,2)))*DRL
               RLTI = -(LOG(RT(NR+1,2))-LOG(RT(NR,2)))*DRL
            ENDIF
            CS = SQRT(TEL*RKEV/(PA(2)*AMM))
            RHO_S = CS*PA(2)*AMM/(PZ(2)*AEE*BB)
            ER(NR) =-BB*( (TIL/TEL)*RHO_S*CS*(RLNI+ALPHA_NEO*RLTI)-EPS/QP(NR)*VTOR(NR))
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE TRERAD

!     ***********************************************

!         BOOTSTRAP CURRENT (NCLASS)

!     ***********************************************

      SUBROUTINE TRAJBS_NCLASS

      USE TRCOMM, ONLY : AJBS, AJBSNC, BB, CJBSP, CJBST, DR, NRMAX, NSMAX, PBSCD, PNSS, PTS, RG, RM, RN, RT, rkind
      USE libitp  
      IMPLICIT NONE
      INTEGER:: NR, NS, NSW
      REAL(rkind)   :: DRPNW, DRTNW, RPNW, RTNW, SUML
      REAL(rkind),DIMENSION(NRMAX)::  AJBSL


      IF(PBSCD.LE.0.D0) RETURN

      NSW=1
      IF(NSW.EQ.0) THEN
         AJBSL(1:NRMAX)=PBSCD*AJBSNC(1:NRMAX)
      ELSE
         DO NR=1,NRMAX-1
            SUML=0.D0
            DO NS=1,NSMAX
               RTNW=0.5D0*(RT(NR+1,NS)+RT(NR,NS))
               RPNW=0.5D0*(RN(NR+1,NS)*RT(NR+1,NS)+RN(NR  ,NS)*RT(NR  ,NS))
               DRTNW=(RT(NR+1,NS)-RT(NR,NS))/DR
               DRPNW=(RN(NR+1,NS)*RT(NR+1,NS)-RN(NR,NS)*RT(NR,NS))/DR
               SUML=SUML+CJBST(NR,NS)*DRTNW/RTNW +CJBSP(NR,NS)*DRPNW/RPNW
            ENDDO
            AJBSL(NR)=-PBSCD*SUML/BB
         ENDDO

         NR=NRMAX
         SUML=0.D0
         DO NS=1,NSMAX
            RTNW=PTS(NS)
            RPNW=PNSS(NS)*PTS(NS)
            DRTNW=DERIV3P(PTS(NS),RT(NR,NS),RT(NR-1,NS),RG(NR),RM(NR),RM(NR-1))
            DRPNW=DERIV3P(PNSS(NS)*PTS(NS),RN(NR  ,NS)*RT(NR  ,NS),RN(NR-1,NS)*RT(NR-1,NS),RG(NR),RM(NR),RM(NR-1))
            SUML=SUML+CJBST(NR,NS)*DRTNW/RTNW +CJBSP(NR,NS)*DRPNW/RPNW
!CC            if(ns.eq.3) write(6,*) NR,PNSS(NS),PTS(NS)
         ENDDO
         AJBSL(NR)=-PBSCD*SUML/BB
      ENDIF
!
      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO

      RETURN
      END SUBROUTINE TRAJBS_NCLASS

!     ***********************************************

!         BOOTSTRAP CURRENT (O. SAUTER)

!     ***********************************************

      SUBROUTINE TRAJBSSAUTER

      USE TRCOMM, ONLY : AJBS, BB, DR, EPSRHO, MDLTPF, NRMAX, NSMAX,&
           & PADD, PBSCD, PNSS, PTS, PZ, QP, RDP, RHOG, RHOM, RKEV,&
           & RN, RPE, RR, RT, RW, TTRHOG, ZEFF, rkind
      USE libitp
      IMPLICIT NONE
      INTEGER:: NR, NS
      REAL(rkind)   :: ANE, DPE, DPI, DRL, DTE, DTI, EPS, EPSS,&
           & F31TEFF, F32EETEFF, F32EITEFF, F34TEFF, FT, FTPF, PE,&
           & PPI, QL, RL31, RL32, RL34, RLNLAME, RLNLAMII, RNM, RNP,&
           & RNTM, RNTP, RNUE, RNUI, RPIM, RPIP, SALFA, SALFA0, TE,&
           & TI, ZEFFL
      REAL(rkind),DIMENSION(NRMAX):: AJBSL, ANI
      REAL(rkind)   :: F31, F32EE, F32EI


      IF(PBSCD.LE.0.D0) RETURN

      DO NR=1,NRMAX-1

         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         QL=ABS(QP(NR))
         ZEFFL=0.5D0*(ZEFF(NR)+ZEFF(NR+1))
         DRL=1.D0/DR

         RNTP=0.D0
         RNP =0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+RN(NR+1,NS)*RT(NR+1,NS)
            RNP =RNP +RN(NR+1,NS)
            RNTM=RNTM+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR  ,NS)
         ENDDO
         RNTP=RNTP+RW(NR+1,1)+RW(NR+1,2)
         RNTM=RNTM+RW(NR  ,1)+RW(NR  ,2)
         RPIP=RNTP+PADD(NR+1)
         RPIM=RNTM+PADD(NR  )

!     ****** ION PARAMETER ******

!     *** ANI  is the the ion density (ni) ***
!     *** TI   is the ion temperature (Ti) ***
!     *** DTI  is the derivative of ion temperature (dTi/dr) ***
!     *** PPI  is the ion pressure (Pi) ***
!     *** DPI  is the derivative of ion pressure (dPi/dr) ***

         ANI(NR)=0.D0
         DO NS=2,NSMAX
            ANI(NR)=ANI(NR)+0.5D0*PZ(NS)*(RN(NR+1,NS)+RN(NR,NS))
         ENDDO
         PPI=0.5D0*(RPIP+RPIM)
         DPI=(RPIP-RPIM)*DRL
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         DTI=(RNTP/RNP-RNTM/RNM)*DRL

         rLnLamii=30.D0-LOG(PZ(2)**3*SQRT(ANI(NR)*1.D20)/(ABS(TI*1.D3)**1.5D0))
         RNUI=4.90D-18*QL*RR*ANI(NR)*1.D20*PZ(2)**4*rLnLamii /(ABS(TI*1.D3)**2*EPSS)
!
!     ****** ELECTORON PARAMETER ******

!     *** ANE  is the the electron density (ne) ***
!     *** TE   is the electron temperature (Te) ***
!     *** PE   is the electron pressure (Pe) ***
!     *** DTE  is the derivative of electron temperature (dTe/dr) ***
!     *** DPE  is the derivative of electron pressure (dPe/dr) ***

         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
         TE= 0.5D0*(RT(NR+1,1)+RT(NR,1))
         PE= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
         DTE=(RT(NR+1,1)-RT(NR,1))*DRL
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL

         rLnLame=31.3D0-LOG(SQRT(ANE*1.D20)/ABS(TE*1.D3))
         RNUE=6.921D-18*QL*RR*ANE*1.D20*ZEFFL*rLnLame /(ABS(TE*1.D3)**2*EPSS)

         RPE=PE/(PE+PPI)
         FT=FTPF(MDLTPF,EPS)

!         F33TEFF=FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE)+0.45D0*(1.D0-FT)*RNUE/ZEFFL**1.5)
         F31TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)   +0.5D0*(1.D0-FT)*RNUE/ZEFFL)
         F32EETEFF=FT/(1.D0+0.26D0*(1.D0-FT)*SQRT(RNUE)+0.18D0*(1.D0-0.37D0*FT)*RNUE/SQRT(ZEFFL))
         F32EITEFF=FT/(1.D0+(1.D0+0.6D0*FT)*SQRT(RNUE) +0.85D0*(1.D0-0.37D0*FT)*RNUE*(1.D0+ZEFFL))
         F34TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)   +0.5D0*(1.D0-0.5D0*FT)*RNUE/ZEFFL)

         SALFA0=-1.17D0*(1.D0-FT)/(1.D0-0.22D0*FT-0.19D0*FT**2)
         SALFA=((SALFA0+0.25D0*(1.D0-FT**2)*SQRT(RNUI))/(1.D0+0.5D0*SQRT(RNUI))+0.315D0*RNUI**2*FT**6) &
     &        /(1.D0+0.15D0*RNUI**2*FT**6)

!         RNZ=0.58D0+0.74D0/(0.76D0+ZEFFL)
!         SGMSPTZ=1.9012D4*(TE*1.D3)**1.5/(ZEFFL*RNZ*rLnLame)
!         SGMNEO=SGMSPTZ*F33(F33TEFF,ZEFFL)

         RL31=F31(F31TEFF,ZEFFL)
         RL32=F32EE(F32EETEFF,ZEFFL)+F32EI(F32EITEFF,ZEFFL)
         RL34=F31(F34TEFF,ZEFFL)

         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV *( RL31*(DPE/PE+DPI/PE)+RL32*DTE/TE &
     &             +RL34*SALFA*(1.D0-RPE)/RPE*DTI/TI)/RDP(NR)/BB
      ENDDO

      NR=NRMAX
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         QL=ABS(QP(NR))
         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
         DRL=1.D0/DR

!     In the following, we assume that
!        1. pressures of beam and fusion at rho=1 are negligible,
!        2. calibration of Pbeam (from UFILE) is not necessary at rho=1.

         RNTP=0.D0
         RNP =0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+PNSS(NS)*PTS(NS)
            RNP =RNP +PNSS(NS)
            RNTM=RNTM+RN(NR-1,NS)*RT(NR-1,NS)+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR-1,NS)+RN(NR  ,NS)
         ENDDO
         RNTM=RNTM+RW(NR-1,1)+RW(NR-1,2)+RW(NR  ,1)+RW(NR  ,2)
         RPIP=RNTP
         RPIM=RNTM+PADD(NR-1)+PADD(NR  )
         RNTM=0.5D0*RNTM
         RNM =0.5D0*RNM
         RPIM=0.5D0*RPIM

!     ****** ION PARAMETER ******

!     *** ANI  is the the ion density (ni) ***
!     *** TI   is the ion temperature (Ti) ***
!     *** DTI  is the derivative of ion temperature (dTi/dr) ***
!     *** PPI  is the ion pressure (Pi) ***
!     *** DPI  is the derivative of ion pressure (dPi/dr) ***

         ANI(NR)=0.D0
         DO NS=2,NSMAX
            ANI(NR)=ANI(NR)+PZ(NS)*PNSS(NS)
         ENDDO
         PPI=RPIP
         DPI=(RPIP-RPIM)*DRL
         TI =RNTP/RNP
         DTI=(RNTP/RNP-RNTM/RNM)*DRL

!         WRITE(6,'(A,I6,4ES12.4)') 'rLnLamii:',NR,PZ(2),ANI(NR),TI,QL
!         WRITE(6,'(A,I6,4ES12.4)') 'RNUI:    ',NR,RR,ANI(NR),rLnLamii,EPSS
         rLnLamii=30.D0-LOG(PZ(2)**3*SQRT(ANI(NR)*1.D20)/(ABS(TI*1.D3)**1.5D0))
         RNUI=4.90D-18*QL*RR*ANI(NR)*1.D20*PZ(2)**4*rLnLamii /(ABS(TI*1.D3)**2*EPSS)
!
!     ****** ELECTORON PARAMETER ******

!     *** ANE  is the the electron density (ne) ***
!     *** TE   is the electron temperature (Te) ***
!     *** PE   is the electron pressure (Pe) ***
!     *** DTE  is the derivative of electron temperature (dTe/dr) ***
!     *** DPE  is the derivative of electron pressure (dPe/dr) ***

         ANE=PNSS(1)
         TE =PTS(1)
         PE =PNSS(1)*PTS(1)
         DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPE=DERIV3P(PNSS(1)*PTS(1),RN(NR,1)*RT(NR,1),RN(NR-1,1)*RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
!
         rLnLame=31.3D0-LOG(SQRT(ANE*1.D20)/ABS(TE*1.D3))
         RNUE=6.921D-18*QL*RR*ANE*1.D20*ZEFFL*rLnLame /(ABS(TE*1.D3)**2*EPSS)
!
         RPE=PE/(PE+PPI)
         FT=FTPF(MDLTPF,EPS)

!     F33TEFF=FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE)+0.45D0*(1.D0-FT)*RNUE/ZEFFL**1.5)
         F31TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)   +0.5D0*(1.D0-FT)*RNUE/ZEFFL)
         F32EETEFF=FT/(1.D0+0.26D0*(1.D0-FT)*SQRT(RNUE)+0.18D0*(1.D0-0.37D0*FT)*RNUE/SQRT(ZEFFL))
         F32EITEFF=FT/(1.D0+(1.D0+0.6D0*FT)*SQRT(RNUE) +0.85D0*(1.D0-0.37D0*FT)*RNUE*(1.D0+ZEFFL))
         F34TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)   +0.5D0*(1.D0-0.5D0*FT)*RNUE/ZEFFL)

         SALFA0=-1.17D0*(1.D0-FT)/(1.D0-0.22D0*FT-0.19D0*FT**2)

!         WRITE(6,'(A,I6,3ES12.4)') 'SALFA:',NR,SALFA0,FT,RNUI
         SALFA=((SALFA0+0.25D0*(1.D0-FT**2)*SQRT(RNUI)) /(1.D0+0.5D0*SQRT(RNUI))+0.315D0*RNUI**2*FT**6) &
     &        /(1.D0+0.15D0*RNUI**2*FT**6)

!     RNZ=0.58D0+0.74D0/(0.76D0+ZEFFL)
!     SGMSPTZ=1.9012D4*(TE*1.D3)**1.5/(ZEFFL*RNZ*rLnLame)
!     SGMNEO=SGMSPTZ*F33(F33TEFF,ZEFFL)
!
         RL31=F31(F31TEFF,ZEFFL)
         RL32=F32EE(F32EETEFF,ZEFFL)+F32EI(F32EITEFF,ZEFFL)
         RL34=F31(F34TEFF,ZEFFL)

         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV*( RL31*(DPE/PE+DPI/PE)+RL32*DTE/TE &
     &             +RL34*SALFA*(1.D0-RPE)/RPE*DTI/TI)/RDP(NR)/BB

      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO

      RETURN
      END SUBROUTINE TRAJBSSAUTER

!     *********************
!     *  Fitting Function *
!     *********************

      FUNCTION F33(X,Z)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind):: X,Z,F33

      F33 = 1.D0-(1.D0+0.36D0/Z)*X+0.59D0/Z*X**2-0.23D0/Z*X**3

      RETURN
      END FUNCTION F33

      FUNCTION F31(X,Z)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind) X,Z,F31

      F31 = (1.D0+1.4D0/(Z+1.D0))*X-1.9D0/(Z+1.D0)*X**2 +0.3D0/(Z+1.D0)*X**3+0.2D0/(Z+1.D0)*X**4

      RETURN
      END FUNCTION F31

      FUNCTION F32EE(X,Z)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind) X,Z,F32EE

      F32EE = (0.05D0+0.62D0*Z)/(Z*(1.D0+0.44D0*Z))*(X-X**4)&
     &       +1.D0/(1.D0+0.22D0*Z)*(X**2-X**4-1.2D0*(X**3-X**4))+1.2D0/(1.D0+0.5D0*Z)*X**4

      RETURN
      END FUNCTION F32EE

      FUNCTION F32EI(X,Z)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind) X,Z,F32EI

      F32EI =-(0.56D0+1.93D0*Z)/(Z*(1.D0+0.44D0*Z))*(X-X**4) &
     &       +4.95D0/(1.D0+2.48D0*Z)*(X**2-X**4-0.55D0*(X**3-X**4))-1.2D0/(1.D0+0.5D0*Z)*X**4

      RETURN
      END FUNCTION F32EI

!     ************************************************

!         BOOTSTRAP CURRENT (Hirshman (Wilson))

!     ************************************************

      SUBROUTINE TRAJBSNEW

      USE TRCOMM, ONLY : AJBS, BB, DR, EPSRHO, NRMAX, NSMAX, PADD, PBSCD, PNSS, PTS, PZ, RDP, RHOG, RHOM, RKEV, RN, RT, &
     &                   RW, TTRHOG, rkind
      USE libitp
      IMPLICIT NONE
      INTEGER:: NR, NS
      REAL(rkind) :: DDD, DDX, DPE, DPI, DRL, DTE, DTI, EPS, FT, PE, PPI, RL31, RL32, RNM, RNP, RNTM, RNTP, RPIM, RPIP, TE, TI
      REAL(rkind),DIMENSION(NRMAX)::  AJBSL, ANI


      IF(PBSCD.LE.0.D0) RETURN

      DO NR=1,NRMAX-1

         EPS=EPSRHO(NR)
!         EPSS=SQRT(EPS)**3
!         QL=ABS(QP(NR))
!         ZEFFL=0.5D0*(ZEFF(NR)+ZEFF(NR+1))
         DRL=1.D0/DR

         RNTP=0.D0
         RNP= 0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+RN(NR+1,NS)*RT(NR+1,NS)
            RNP =RNP +RN(NR+1,NS)
            RNTM=RNTM+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR  ,NS)
         ENDDO
         RNTP=RNTP+RW(NR+1,1)+RW(NR+1,2)
         RNTM=RNTM+RW(NR  ,1)+RW(NR  ,2)
         RPIP=RNTP+PADD(NR+1)
         RPIM=RNTM+PADD(NR  )

!     ****** ION PARAMETER ******

!     ***** ANI  is the the ion density (ni) *****
!     ***** TI   is the ion temperature (Ti) *****
!     ***** DTI  is the derivative of ion temperature (dTi/dr) *****
!     ***** PPI  is the ion pressure (Pi) *****
!     ***** DPI  is the derivative of ion pressure (dPi/dr) *****
!     ***** VTI  is the ion velocity (VTi) *****

         ANI(NR)=0.D0
         DO NS=2,NSMAX
            ANI(NR)=ANI(NR)+0.5D0*PZ(NS)*(RN(NR+1,NS)+RN(NR,NS))
         ENDDO
         PPI=0.5D0*(RPIP+RPIM)
         DPI=(RPIP-RPIM)*DRL
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         DTI=(RNTP/RNP-RNTM/RNM)*DRL
!         VTI=SQRT(ABS(TI)*RKEV/AMM)

!         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
!         rLnLam=17.3D0-DLOG(ANE)*0.5D0+DLOG(ABS(TI))*1.5D0
!         TAUI=12.D0*PI*SQRT(PI)*EPS0**2*SQRT(AMM)
!     &             *(ABS(TI)*RKEV)**1.5D0/(ANI(NR)*1.D20
!     &             *ZEFFL**4*AEE**4*rLnLam)

!         RNUI=QL*RR/(TAUI*VTI*EPSS)

!     ****** ELECTORON PARAMETER ******

!     ***** ANE  is the the electron density (ne) *****
!     ***** TE   is the electron temperature (Te) *****
!     ***** DTE  is the derivative of electron temperature (dTe/dr) ****
!     ***** PE   is the electron pressure (Pe) *****
!     ***** DPE  is the derivative of electron pressure (dPe/dr) *****
!     ***** VTE  is the electron velocity (VTe) *****

         TE= 0.5D0*(RT(NR+1,1)+RT(NR,1))
         PE= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
         DTE=(RT(NR+1,1)-RT(NR,1))*DRL
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL
!         VTE=SQRT(ABS(TE)*RKEV/AME)

!         rLnLam=15.2D0-DLOG(ANE)*0.5D0+DLOG(ABS(TE))
!         TAUE=6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)*(ABS(TE)*RKEV)**1.5D0/(ANI*1.D20*ZEFFL**2*AEE**4*rLnLam)

!         RNUE=QL*RR/(TAUE*VTE*EPSS)

!         FT=FTPF(MDLTPF,EPS)
         FT=(1.46D0*SQRT(EPS)+2.4D0*EPS)/(1.D0-EPS)**1.5D0

!         DDX=2.4D0+5.4D0*FT+2.6D0*FT**2
!         DDD=-1.17D0/(1.D0+0.46D0*FT)
!         C1=(4.D0+2.6D0*FT)/((1.D0+1.02D0*DSQRT(RNUE)+1.07D0*RNUE)*(1.D0+1.07D0*EPSS*RNUE))
!         C2=C1*TI/TE
!         C3=(7.D0+6.5D0*FT)/((1.D0+0.57D0*DSQRT(RNUE)+0.61D0*RNUE)*(1.D0+0.61D0*EPSS*RNUE))-2.5D0*C1
!         C4=((DDD+0.35D0*DSQRT(RNUI))/(1.D0+0.7D0*DSQRT(RNUI))
!     &      +2.1D0*EPS**3*RNUI**2)*C2/((1.D0-EPS**3*RNUI**2)*(1.D0+EPS**3*RNUE**2))

!         AJBSL(NR)=-PBSCD*FT*PE*1.D20*RKEV/(DDX*BP(NR))*(C1*(DPE/PE)
!     &         +C2*(DPI/PPI)+C3*(DTE/TE)+C4*(DTI/TI))

!     *** S. P. Hirshman, Phys Fluids 31, 1988 3150 ***
!     *** cited (H.R. WILSON, Nucl.Fusion 32,no.2,1992 259-263) ***

         DDX=1.414D0*PZ(2)+PZ(2)**2+FT*(0.754D0+2.657D0*PZ(2) &
     &        +2.D0*PZ(2)**2)+FT**2*(0.348D0+1.243D0*PZ(2)+PZ(2)**2)
         RL31= FT*(0.754D0+2.210D0*PZ(2)+PZ(2)**2 +FT*(0.348D0+1.243D0*PZ(2)+PZ(2)**2))/DDX
         RL32=-FT*(0.884D0+2.074D0*PZ(2))/DDX
         DDD=-1.172D0/(1.D0+0.462D0*FT)

         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV*(RL31*((DPE/PE)+(TI/(PZ(2)*TE)) &
     &        *((DPI/PPI)+DDD*(DTI/TI)))+RL32*(DTE/TE))/RDP(NR)/BB
      ENDDO

      NR=NRMAX
         EPS=EPSRHO(NR)
!         EPSS=SQRT(EPS)**3
!         QL=ABS(QP(NR))
!         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
         DRL=1.D0/DR

!     In the following, we assume that
!        1. pressures of beam and fusion at rho=1 are negligible,
!        2. calibration of Pbeam (from UFILE) is not necessary at rho=1.

         RNTP=0.D0
         RNP= 0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+PNSS(NS)*PTS(NS)
            RNP =RNP +PNSS(NS)
            RNTM=RNTM+RN(NR-1,NS)*RT(NR-1,NS)+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR-1,NS)+RN(NR  ,NS)
         ENDDO
         RNTM=RNTP+RW(NR-1,1)+RW(NR-1,2)+RW(NR  ,1)+RW(NR  ,2)
         RPIP=RNTP
         RPIM=RNTM+PADD(NR-1)+PADD(NR  )
         RNTM=0.5D0*RNTM
         RNM =0.5D0*RNM
         RPIM=0.5D0*RPIM

!     ****** ION PARAMETER ******

!     ***** ANI  is the the ion density (ni) *****
!     ***** TI   is the ion temperature (Ti) *****
!     ***** DTI  is the derivative of ion temperature (dTi/dr) *****
!     ***** PPI  is the ion pressure (Pi) *****
!     ***** DPI  is the derivative of ion pressure (dPi/dr) *****
!     ***** VTI  is the ion velocity (VTi) *****

         ANI(NR)=0.D0
         DO NS=2,NSMAX
            ANI(NR)=ANI(NR)+PZ(NS)*PNSS(NS)
         ENDDO
         PPI=0.5D0*(RPIP+RPIM)
         DPI=(RPIP-RPIM)*DRL
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         DTI=(RNTP/RNP-RNTM/RNM)*DRL
!         VTI=SQRT(ABS(TI)*RKEV/AMM)

!         ANE=PNSS(1)
!         rLnLam=17.3D0-DLOG(ANE)*0.5D0+DLOG(ABS(TI))*1.5D0
!         TAUI=12.D0*PI*SQRT(PI)*EPS0**2*SQRT(AMM)
!     &             *(ABS(TI)*RKEV)**1.5D0/(ANI(NR)*1.D20
!     &             *ZEFFL**4*AEE**4*rLnLam)

!         RNUI=QL*RR/(TAUI*VTI*EPSS)

!     ****** ELECTORON PARAMETER ******

!     ***** ANE  is the the electron density (ne) *****
!     ***** TE   is the electron temperature (Te) *****
!     ***** DTE  is the derivative of electron temperature (dTe/dr) ****
!     ***** PE   is the electron pressure (Pe) *****
!     ***** DPE  is the derivative of electron pressure (dPe/dr) *****
!     ***** VTE  is the electron velocity (VTe) *****

         TE =PTS(1)
         PE =PNSS(1)*PTS(1)
         DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPE=DERIV3P(PNSS(1)*PTS(1),RN(NR,1)*RT(NR,1),RN(NR-1,1)*RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
!         VTE=SQRT(ABS(TE)*RKEV/AME)

!         rLnLam=15.2D0-DLOG(ANE)*0.5D0+DLOG(ABS(TE))
!         TAUE=6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)
!     &             *(ABS(TE)*RKEV)**1.5D0/(ANI*1.D20*ZEFFL**2*AEE**4*rLnLam)

!         RNUE=QL*RR/(TAUE*VTE*EPSS)

!         FT=FTPF(MDLTPF,EPS)
         FT=(1.46D0*SQRT(EPS)+2.4D0*EPS)/(1.D0-EPS)**1.5D0

!         DDX=2.4D0+5.4D0*FT+2.6D0*FT**2
!         DDD=-1.17D0/(1.D0+0.46D0*FT)
!         C1=(4.D0+2.6D0*FT)/((1.D0+1.02D0*DSQRT(RNUE)+1.07D0*RNUE)*(1.D0+1.07D0*EPSS*RNUE))
!         C2=C1*TI/TE
!         C3=(7.D0+6.5D0*FT)/((1.D0+0.57D0*DSQRT(RNUE)+0.61D0*RNUE)*(1.D0+0.61D0*EPSS*RNUE))
!     &      -2.5D0*C1
!         C4=((DDD+0.35D0*DSQRT(RNUI))/(1.D0+0.7D0*DSQRT(RNUI))+2.1D0*EPS**3*RNUI**2)*C2
!     &      /((1.D0-EPS**3*RNUI**2)*(1.D0+EPS**3*RNUE**2))

!         AJBSL(NR)=-PBSCD*FT*PE*1.D20*RKEV/(DDX*BP(NR))*(C1*(DPE/PE)
!     &         +C2*(DPI/PPI)+C3*(DTE/TE)+C4*(DTI/TI))

!     *** S. P. Hirshman, Phys Fluids 31, 1988 3150 ***
!     *** cited (H.R. WILSON, Nucl.Fusion 32,no.2,1992 259-263) ***

         DDX=1.414D0*PZ(2)+PZ(2)**2+FT*(0.754D0+2.657D0*PZ(2) &
     &        +2.D0*PZ(2)**2)+FT**2*(0.348D0+1.243D0*PZ(2)+PZ(2)**2)
         RL31=FT*( 0.754D0+2.21D0*PZ(2)+PZ(2)**2+FT*(0.348D0+1.243D0 &
     &            *PZ(2)+PZ(2)**2))/DDX
         RL32=-FT*(0.884D0+2.074D0*PZ(2))/DDX
         DDD=-1.172D0/(1.D0+0.462D0*FT)

         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV &
     &        *(RL31*((DPE/PE)+(TI/(PZ(2)*TE))*((DPI/PPI)+DDD*(DTI/TI)))+RL32*(DTE/TE))/RDP(NR)/BB

      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO

      RETURN
      END SUBROUTINE TRAJBSNEW

!     ***********************************************************

!           BOOTSTRAP CURRENT (Hinton & Hazeltine)

!     ***********************************************************

      SUBROUTINE TRAJBS

      USE TRCOMM, ONLY : AJBS, AME, AMM, BB, BP, DR, EPSRHO, NRMAX, NSMAX, PA, PBSCD, PNSS, PTS, PZ, QP, RHOG, RHOM, &
     &                   RJCB, RKEV, RN, RR, RT, ZEFF, rkind
      USE libitp
      IMPLICIT NONE
      INTEGER:: NR
      REAL(rkind)   :: A, AMA, AMD, AMT, ANA, ANDX, ANE, ANT, BPL, DPA, DPD, DPE, DPT, DRL, DTA, DTD, DTE, DTT, EPS, EPSS, FACT, &
     &             FTAUE, FTAUI, H, PAL, PDL, PEL, PTL, RK13E, RK23E, RK3A, RK3D, RK3T, RNUA, RNUD, RNUE, RNUT, TAL, TAUA,   &
     &             TAUD, TAUE, TAUT, TDL, TEL, TTL, VTA, VTD, VTE, VTT, ZEFFL
      REAL(rkind),DIMENSION(NRMAX):: AJBSL

      REAL(rkind):: RK13=2.30D0, RA13=1.02D0, RB13=1.07D0, RC13=1.07D0, RK23=4.19D0, RA23=0.57D0, RB23=0.61D0, RC23=0.61D0

!     ZEFF=1

!      DATA RK11,RA11,RB11,RC11/1.04D0,2.01D0,1.53D0,0.89D0/
!      DATA RK12,RA12,RB12,RC12/1.20D0,0.76D0,0.67D0,0.56D0/
!      DATA RK22,RA22,RB22,RC22/2.55D0,0.45D0,0.43D0,0.43D0/
!!      DATA RK13,RA13,RB13,RC13/2.30D0,1.02D0,1.07D0,1.07D0/
!!      DATA RK23,RA23,RB23,RC23/4.19D0,0.57D0,0.61D0,0.61D0/
!      DATA RK33,RA33,RB33,RC33/1.83D0,0.68D0,0.32D0,0.66D0/
!      DATA RK2 ,RA2 ,RB2 ,RC2 /0.66D0,1.03D0,0.31D0,0.74D0/

      IF(PBSCD.LE.0.D0) RETURN

      AMD=PA(2)*AMM
      AMT=PA(3)*AMM
      AMA=PA(4)*AMM

      DO NR=1,NRMAX-1

         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
         ANDX=0.5D0*(RN(NR+1,2)+RN(NR,2))
         ANT =0.5D0*(RN(NR+1,3)+RN(NR,3))
         ANA =0.5D0*(RN(NR+1,4)+RN(NR,4))
         TEL=ABS(0.5D0*(RT(NR+1,1)+RT(NR,1)))
         TDL=ABS(0.5D0*(RT(NR+1,2)+RT(NR,2)))
         TTL=ABS(0.5D0*(RT(NR+1,3)+RT(NR,3)))
         TAL=ABS(0.5D0*(RT(NR+1,4)+RT(NR,4)))
         PEL=0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
         PDL=0.5D0*(RN(NR+1,2)*RT(NR+1,2)+RN(NR,2)*RT(NR,2))
         PTL=0.5D0*(RN(NR+1,3)*RT(NR+1,3)+RN(NR,3)*RT(NR,3))
         PAL=0.5D0*(RN(NR+1,4)*RT(NR+1,4)+RN(NR,4)*RT(NR,4))
         ZEFFL=0.5D0*(ZEFF(NR+1)+ZEFF(NR))

         TAUE = FTAUE(ANE,ANDX,TEL,ZEFFL)
         TAUD = FTAUI(ANE,ANDX,TDL,PZ(2),PA(2))
         TAUT = FTAUI(ANE,ANT ,TTL,PZ(3),PA(3))
         TAUA = FTAUI(ANE,ANA ,TAL,PZ(4),PA(4))

         VTE=SQRT(TEL*RKEV/AME)
         VTD=SQRT(TDL*RKEV/AMD)
         VTT=SQRT(TTL*RKEV/AMT)
         VTA=SQRT(TAL*RKEV/AMA)

         RNUE=ABS(QP(NR))*RR/(TAUE*VTE*EPSS)
         RNUD=ABS(QP(NR))*RR/(TAUD*VTD*EPSS)
         RNUT=ABS(QP(NR))*RR/(TAUT*VTT*EPSS)
         RNUA=ABS(QP(NR))*RR/(TAUA*VTA*EPSS)

!         RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)+(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
!         RK12E=RK12*(1.D0/(1.D0+RA12*SQRT(RNUE)+RB12*RNUE)+(EPSS*RC12)**2/RB12*RNUE/(1.D0+RC12*RNUE*EPSS))
!         RK22E=RK22*(1.D0/(1.D0+RA22*SQRT(RNUE)+RB22*RNUE)+(EPSS*RC22)**2/RB22*RNUE/(1.D0+RC22*RNUE*EPSS))
         RK13E=RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)/(1.D0+RC13*RNUE*EPSS)
         RK23E=RK23/(1.D0+RA23*SQRT(RNUE)+RB23*RNUE)/(1.D0+RC23*RNUE*EPSS)
!         RK33E=RK33/(1.D0+RA33*SQRT(RNUE)+RB33*RNUE)/(1.D0+RC33*RNUE*EPSS)

!         RK2D =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUD)+RB2 *RNUD)+(EPSS*RC2 )**2/RB2 *RNUD/(1.D0+RC2 *RNUD*EPSS))
!         RK2T =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUT)+RB2 *RNUT)+(EPSS*RC2 )**2/RB2 *RNUT/(1.D0+RC2 *RNUT*EPSS))
!         RK2A =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUA)+RB2 *RNUA)+(EPSS*RC2 )**2/RB2 *RNUA/(1.D0+RC2 *RNUA*EPSS))
         RK3D=((1.17D0-0.35D0*SQRT(RNUD))/(1.D0+0.7D0*SQRT(RNUD))-2.1D0*(RNUD*EPSS)**2)/(1.D0+(RNUD*EPSS)**2)
         RK3T=((1.17D0-0.35D0*SQRT(RNUT))/(1.D0+0.7D0*SQRT(RNUT))-2.1D0*(RNUT*EPSS)**2)/(1.D0+(RNUT*EPSS)**2)
         RK3A=((1.17D0-0.35D0*SQRT(RNUA))/(1.D0+0.7D0*SQRT(RNUA))-2.1D0*(RNUA*EPSS)**2)/(1.D0+(RNUA*EPSS)**2)

         DRL=RJCB(NR)/DR
         DTE=(RT(NR+1,1)-RT(NR,1))*DRL
         DTD=(RT(NR+1,2)-RT(NR,2))*DRL
         DTT=(RT(NR+1,3)-RT(NR,3))*DRL
         DTA=(RT(NR+1,4)-RT(NR,4))*DRL
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL
         DPD=(RN(NR+1,2)*RT(NR+1,2)-RN(NR,2)*RT(NR,2))*DRL
         DPT=(RN(NR+1,3)*RT(NR+1,3)-RN(NR,3)*RT(NR,3))*DRL
         DPA=(RN(NR+1,4)*RT(NR+1,4)-RN(NR,4)*RT(NR,4))*DRL
         BPL=BP(NR)

         FACT=1.D0/(1.D0+(RNUE*EPS)**2)
         IF(NSMAX.EQ.2) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
         ELSEIF(NSMAX.EQ.3) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL &
     &       +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)+(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)
         ELSEIF(NSMAX.EQ.4) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL) &
     &       +(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)+(PAL/PEL)*(DPA/PAL-RK3A*FACT*DTA/TAL)
         ENDIF
         H=BB/(BB+BP(NR))
         AJBSL(NR)=-PBSCD*H*SQRT(EPS)*ANE*1.D20*TEL*RKEV/BPL *(RK13E*A+RK23E*DTE/TEL)
!
!        WRITE(6,'(I3,1P6E12.4)') NR,RNUE,RK13E,RK23E,A,BPL,AJBSL(NR)
!         IF(NR.EQ.NRMAX) THEN
!            WRITE(6,'(1P6E12.4)') DN,DTE,ANE,PNSS(1),TE,PTS(1)
!            WRITE(6,'(1P6E12.4)') RK13E*(RN(NR,1)*RT(NR,1)
!     &                   +RN(NR,2)*RT(NR,2)
!     &                   +RN(NR,3)*RT(NR,3)
!     &                   +RN(NR,4)*RT(NR,4)
!     &                   +RW(NR,1)
!     &                   +RW(NR,2)         )*DN/ANE,
!     &            +(RK23E-1.5D0*RK13E)*RN(NR,1)*DTE,
!     &            +RK13D*(RK23D-1.5D0)*RN(NR,2)*DTD,
!     &            BPL,AJBS(NR-1),AJBS(NR)
!         ENDIF
      ENDDO

      NR=NRMAX
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         ANE=PNSS(1)
         ANDX=PNSS(2)
!         ANT =PNSS(3)
!         ANA =PNSS(4)
         TEL=ABS(PTS(1))
         TDL=ABS(PTS(2))
         TTL=ABS(PTS(3))
         TAL=ABS(PTS(4))
         PEL=PNSS(1)*PTS(1)
         PDL=PNSS(2)*PTS(2)
         PTL=PNSS(3)*PTS(3)
         PAL=PNSS(4)*PTS(4)
         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)

         TAUE = FTAUE(ANE,ANDX,TEL,ZEFFL)
         TAUD = FTAUI(ANE,ANDX,TDL,PZ(2),PA(2))
         TAUT = FTAUI(ANE,ANT ,TTL,PZ(3),PA(3))
         TAUA = FTAUI(ANE,ANA ,TAL,PZ(4),PA(4))

         VTE=SQRT(TEL*RKEV/AME)
         VTD=SQRT(TDL*RKEV/AMD)
         VTT=SQRT(TTL*RKEV/AMT)
         VTA=SQRT(TAL*RKEV/AMA)

         RNUE=ABS(QP(NR))*RR/(TAUE*VTE*EPSS)
         RNUD=ABS(QP(NR))*RR/(TAUD*VTD*EPSS)
         RNUT=ABS(QP(NR))*RR/(TAUT*VTT*EPSS)
         RNUA=ABS(QP(NR))*RR/(TAUA*VTA*EPSS)

!         RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)+(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
!         RK12E=RK12*(1.D0/(1.D0+RA12*SQRT(RNUE)+RB12*RNUE)+(EPSS*RC12)**2/RB12*RNUE/(1.D0+RC12*RNUE*EPSS))
!         RK22E=RK22*(1.D0/(1.D0+RA22*SQRT(RNUE)+RB22*RNUE)+(EPSS*RC22)**2/RB22*RNUE/(1.D0+RC22*RNUE*EPSS))
         RK13E=RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE) /(1.D0+RC13*RNUE*EPSS)
         RK23E=RK23/(1.D0+RA23*SQRT(RNUE)+RB23*RNUE) /(1.D0+RC23*RNUE*EPSS)
!         RK33E=RK33/(1.D0+RA33*SQRT(RNUE)+RB33*RNUE)/(1.D0+RC33*RNUE*EPSS)

!         RK2D =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUD)+RB2 *RNUD)+(EPSS*RC2 )**2/RB2 *RNUD/(1.D0+RC2 *RNUD*EPSS))
!         RK2T =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUT)+RB2 *RNUT)+(EPSS*RC2 )**2/RB2 *RNUT/(1.D0+RC2 *RNUT*EPSS))
!         RK2A =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUA)+RB2 *RNUA)+(EPSS*RC2 )**2/RB2 *RNUA/(1.D0+RC2 *RNUA*EPSS))
         RK3D=((1.17D0-0.35D0*SQRT(RNUD))/(1.D0+0.7D0*SQRT(RNUD))-2.1D0*(RNUD*EPSS)**2)/(1.D0+(RNUD*EPSS)**2)
         RK3T=((1.17D0-0.35D0*SQRT(RNUT))/(1.D0+0.7D0*SQRT(RNUT))-2.1D0*(RNUT*EPSS)**2)/(1.D0+(RNUT*EPSS)**2)
         RK3A=((1.17D0-0.35D0*SQRT(RNUA))/(1.D0+0.7D0*SQRT(RNUA))-2.1D0*(RNUA*EPSS)**2)/(1.D0+(RNUA*EPSS)**2)

         DRL=RJCB(NR)/DR
         DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DTD=DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DTT=DERIV3P(PTS(3),RT(NR,3),RT(NR-1,3),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DTA=DERIV3P(PTS(4),RT(NR,4),RT(NR-1,4),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPE=DERIV3P(PNSS(1)*PTS(1),RN(NR,1)*RT(NR,1),RN(NR-1,1)*RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPD=DERIV3P(PNSS(2)*PTS(2),RN(NR,2)*RT(NR,2),RN(NR-1,2)*RT(NR-1,2),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPT=DERIV3P(PNSS(3)*PTS(3),RN(NR,3)*RT(NR,3),RN(NR-1,3)*RT(NR-1,3),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPA=DERIV3P(PNSS(4)*PTS(4),RN(NR,4)*RT(NR,4),RN(NR-1,4)*RT(NR-1,4),RHOG(NR),RHOM(NR),RHOM(NR-1))
         BPL=BP(NR)

         FACT=1.D0/(1.D0+(RNUE*EPS)**2)
         IF(NSMAX.EQ.2) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
         ELSEIF(NSMAX.EQ.3) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL &
     &       +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)+(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)
         ELSEIF(NSMAX.EQ.4) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL) &
     &       +(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)+(PAL/PEL)*(DPA/PAL-RK3A*FACT*DTA/TAL)
         ENDIF
         H=BB/(BB+BP(NR))
         AJBSL(NR)=-PBSCD*H*SQRT(EPS)*ANE*1.D20*TEL*RKEV/BPL*(RK13E*A+RK23E*DTE/TEL)
!
!        WRITE(6,'(I3,1P6E12.4)') NR,RNUE,RK13E,RK23E,A,BPL,AJBSL(NR)
!         IF(NR.EQ.NRMAX) THEN
!            WRITE(6,'(1P6E12.4)') DN,DTE,ANE,PNSS(1),TE,PTS(1)
!            WRITE(6,'(1P6E12.4)') RK13E*(RN(NR,1)*RT(NR,1)
!     &                   +RN(NR,2)*RT(NR,2)
!     &                   +RN(NR,3)*RT(NR,3)
!     &                   +RN(NR,4)*RT(NR,4)
!     &                   +RW(NR,1)
!     &                   +RW(NR,2)         )*DN/ANE,
!     &            +(RK23E-1.5D0*RK13E)*RN(NR,1)*DTE,
!     &            +RK13D*(RK23D-1.5D0)*RN(NR,2)*DTD,
!     &            BPL,AJBS(NR-1),AJBS(NR)
!         ENDIF

      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO

      RETURN
      END SUBROUTINE TRAJBS

!     ***********************************************************

!           CALCULATE AJ, AJOH, POH, EZOH

!     ***********************************************************

      SUBROUTINE TRAJOH

      USE TRCOMM, ONLY : &
           AJ, AJBS, AJNB, AJOH, AJRF, AJTOR, BB, DR, DVRHO, ETA, EZOH, &
           MDLEQB, MDLJQ, MDLUF, NRMAX, POH, RMU0, RR, TTRHO, TTRHOG, &
           RDPVRHOG, PI,abvrhog, rkind
      IMPLICIT NONE
      INTEGER:: NR
      REAL(rkind)   :: FACTOR0, FACTORM, FACTORP

      IF(MDLEQB.EQ.1.OR.MDLJQ.EQ.1.OR.(MDLUF.EQ.0.OR.MDLUF.EQ.3)) THEN
      NR=1
         FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
         FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
         AJ(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
      DO NR=2,NRMAX
         FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
         FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
         FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
         AJ(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
      ENDDO
      NR=1
         FACTOR0=RR/(RMU0*DVRHO(NR))
         FACTORP=ABVRHOG(NR  )
         AJTOR(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
      DO NR=2,NRMAX
         FACTOR0=RR/(RMU0*DVRHO(NR))
         FACTORM=ABVRHOG(NR-1)
         FACTORP=ABVRHOG(NR  )
         AJTOR(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
      ENDDO
      ENDIF

      AJOH(1:NRMAX) = AJ(1:NRMAX)-(AJNB(1:NRMAX)+AJRF(1:NRMAX)+AJBS(1:NRMAX))
      EZOH(1:NRMAX) = ETA(1:NRMAX)*AJOH(1:NRMAX)
!!!      IF(KUFDEV.EQ.'lhd') THEN
!!!         POH(1:NRMAX)  = 0.D0
!!!      ELSE
         POH(1:NRMAX)  = EZOH(1:NRMAX)*AJOH(1:NRMAX)
!!!      ENDIF

      RETURN
      END SUBROUTINE TRAJOH

!     ***********************************************************

!           SAWTOOTH OSCILLATION

!     ***********************************************************

      SUBROUTINE TRSAWT

      USE TRCOMM, ONLY : AR1RHOG, ARRHOG, BP, DR, DVRHO, DVRHOG, MDLST, NRMAX, NSMAX, PI, QP, RDP, RG, RM, RN, RR, RT, &
     &                   T, TTRHOG, RDPVRHOG, rkind
      IMPLICIT NONE
      INTEGER:: IONE, IZEROX, LN, LQ, LT, NR, NS
      REAL(rkind)   :: RNN, RTN, SUML, SUML1, SUML2
      REAL(rkind),DIMENSION(NRMAX):: QONE


      IF(MDLST.EQ.0) RETURN

      IF(MOD(MDLST,2).EQ.1) THEN
         LT = 1
      ELSE
         LT = 0
      ENDIF
      IF(MOD(MDLST/2,2).EQ.1) THEN
         LN = 1
      ELSE
         LN = 0
      ENDIF
      IF(MOD(MDLST/4,2).EQ.1) THEN
         LQ = 1
      ELSE
         LQ = 0
      ENDIF

      WRITE(6,601) MDLST,T
  601 FORMAT(' ','# SAWTOOTH OSCILLATION -TYPE ',I1, ' AT ',F7.3,' SEC')

      SUML=0.D0
      DO NR=1,NRMAX
         SUML = SUML+(1.D0/ABS(QP(NR))-1.D0)*(DVRHOG(NR)/2.D0*PI*RR)*DR
         IF(SUML.LT.0.D0) GOTO 1000
      ENDDO
      NR=NRMAX

 1000 IZEROX=NR

      DO NR=1,NRMAX
         IF(QP(NR).GE.1.D0) GOTO 2000
      ENDDO
      NR=NRMAX

 2000 IONE=NR

      SUML = SUM((1.D0/ABS(QP(IONE:IZEROX))-1.D0)*DVRHOG(IONE:IZEROX))/(2.D0*PI*RR)*DR

      QONE(1:IZEROX) = 1.D0+SUML*4.D0*RG(1:IZEROX)**2/RG(IZEROX)**4

      IF(LT.EQ.1) THEN
         DO NS=1,NSMAX
            SUML1 = SUM(RN(1:IZEROX,NS)                *DVRHO(1:IZEROX))
            SUML2 = SUM(RN(1:IZEROX,NS)*RT(1:IZEROX,NS)*DVRHO(1:IZEROX))
            RTN = SUML2/SUML1
            RT(1:IZEROX,NS) = RTN
         ENDDO
      ENDIF

      IF(LN.EQ.1) THEN
         DO NS=1,NSMAX
            SUML1 = SUM(                DVRHO(1:IZEROX))
            SUML2 = SUM(RN(1:IZEROX,NS)*DVRHO(1:IZEROX))
            RNN = SUML2/SUML1
            RN(1:IZEROX,NS) = RNN
         ENDDO
      ENDIF

      IF(LQ.EQ.1) THEN
         QP(1:IZEROX) = 1.D0/QONE(1:IZEROX)
         RDPVRHOG(1:IZEROX)=TTRHOG(1:IZEROX)*ARRHOG(1:IZEROX) &
     &                /(4.D0*PI**2*QP(1:IZEROX))
         RDP(1:IZEROX)=RDPVRHOG(1:IZEROX)*DVRHOG(1:IZEROX)
         BP(1:IZEROX) =AR1RHOG(1:IZEROX)*RDP(1:IZEROX)/RR
      ENDIF

      WRITE(6,602) RM(IONE),RM(IZEROX),RTN,RNN
602   FORMAT(' ',' R-ONE,R-ZERO,RTN,RNN = ',4F8.3)

      DO NR=1,IZEROX+2
         WRITE(6,'(A,I5,1P5E12.4)') 'NR,R,Q,B,Q,S=', &
                       NR,RM(NR),QP(NR),BP(NR),QONE(NR),SUML
      ENDDO
      RETURN
      END SUBROUTINE TRSAWT

!     ***********************************************************

!           TRAPPED PARTICLE FRACTION

!     ***********************************************************

      FUNCTION FTPF(ID,EPS)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      INTEGER                  :: ID, IERR, N, I
      INTEGER,PARAMETER        :: IMAX=20
      REAL(rkind)                     :: EPS, EPSC, FTLL, FTUL, OMEGA, PI, S
      REAL(rkind),DIMENSION(IMAX,IMAX):: TABLE
      REAL(rkind)                     :: FTPF
      REAL(rkind)                     :: FTL, FTU
      EXTERNAL FTL, FTU

      IF(ID.EQ.1) THEN
!  Y. R. Lin-Liu and R. L. Miller, PoP 2 1666 (1995), eqs(7)(13)(18)(19)
         PI=3.14159265358979323846D0
         EPSC=1.D-9
         I=IMAX
         FTUL=1.D0-(1.D0-1.5D0*SQRT(EPS)+0.5D0*EPS**1.5D0)/SQRT(1-EPS**2)
         CALL RMBRG(0.D0,2.D0*PI,EPSC,S,I,N,IERR,TABLE,EPS,FTL)
         FTLL=1.D0-(1.D0-EPS)**1.5D0/SQRT(1.D0+EPS)*(S/(2.D0*PI))
         OMEGA=(3.D0*SQRT(2.D0)/2.D0*0.69D0-3.D0*SQRT(2.D0)/PI)/(1.5D0-3.D0*SQRT(2.D0)/PI)
         FTPF=OMEGA*FTUL+(1.D0-OMEGA)*FTLL
      ELSEIF(ID.EQ.2) THEN
!  S. P. Hirshman et al., NF 17 611 (1977)
         FTPF=1.D0-(1.D0-EPS)**2.D0/(DSQRT(1.D0-EPS**2)*(1.D0+1.46D0*DSQRT(EPS)))
      ELSEIF(ID.EQ.3) THEN
!  Y. R. Lin-Liu and R. L. Miller, PoP 2 1666 (1995), eqs(16)(17)(18)
         PI=3.14159265358979323846D0
         FTUL=1.5D0*SQRT(EPS)
         FTLL=3.D0*SQRT(2.D0)/PI*SQRT(EPS)
         FTPF=0.75D0*FTUL+0.25D0*FTLL
      ELSEIF(ID.EQ.4) THEN
!  M. N. Rosenbluth et al., PoF 15 116 (1972)
         FTPF=1.46D0*SQRT(EPS)
      ELSE
!  Y. B. Kim et al., PoF B 3 2050 (1991) eq(C18), default
         FTPF=1.46D0*SQRT(EPS)-0.46D0*(EPS)**1.5D0
      ENDIF

      RETURN
      END FUNCTION FTPF

!     *********************************************

!           FUNCTION FOR ROMBERG INTEGRATION

!     *********************************************

      FUNCTION FTU(X,EPS)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind):: EPS, X, FTU

      FTU = X/SQRT(1.D0-X*(1.D0-EPS))

      RETURN
      END FUNCTION FTU

      FUNCTION FTL(X,EPS)
!
      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind):: X, EPS,FTL
      REAL(rkind):: H

      H = (1.D0 - EPS) / (1.D0 + EPS * COS(X))
      FTL = (1.D0 - SQRT(1.D0 - H) * (1.D0 + 0.5D0 * H)) / H**2

      RETURN
      END FUNCTION FTL
!
!     *********************************************

!           ROMBERG INTEGRATION METHOD

!     *********************************************

      SUBROUTINE RMBRG(A,B,EPS,S,IMAX,N,IERR,T,ARG,F)

!     <input>
!        A     : lower bound
!        B     : upper bound
!        EPS   : stopping criterion
!        IMAX  : maximum division number
!        F     : formula of integrand
!     <output>
!        S     : integration value
!        N     : division number
!        IERR  : error indicator
!        T     : Romberg T table

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind),INTENT(IN)      ::  A, B, ARG, EPS
      REAL(rkind),INTENT(OUT)     ::  S
      INTEGER,INTENT(INOUT):: IMAX
      INTEGER,INTENT(OUT)  :: IERR, N
      REAL(rkind),DIMENSION(IMAX,IMAX),INTENT(OUT)::  T
      INTEGER::  I, J, K, N2
      REAL(rkind)   ::  F, X, H, S1, Y1, Y2

      EXTERNAL F

      DO K=1,IMAX
         N=2**(K-1)
         N2=N/2
         H=(B-A)/N
         Y1=0
         IF(N.EQ.1) THEN
            Y2=(F(A,ARG)+F(B,ARG))/2
         ELSE
            DO I=1,N2
               X=A+(2*I-1)*H
               Y1=Y1+F(X,ARG)
            ENDDO
            Y2=Y2+Y1
            Y1=0.D0
         ENDIF
         S=H*Y2
         T(K,1)=S
         IF(K.LE.1) GOTO 10
         DO J=2,K
            T(K,J)=T(K,J-1)+(T(K,J-1)-T(K-1,J-1))/(4**(J-1)-1)
         ENDDO
         S=T(K,K)
         S1=T(K,K-1)
         IF(ABS(S-S1).LT.EPS) THEN
            IERR=0
            IMAX=K
            RETURN
         ENDIF
 10      CONTINUE
      ENDDO
      IERR=1

      RETURN
      END SUBROUTINE RMBRG

!     ***********************************************************

!           COULOMB LOGARITHM

!     ***********************************************************

      FUNCTION COULOG(NS1,NS2,ANEL,TL)

!     ANEL : electron density [10^20 /m^3]
!     TL   : electron or ion temperature [keV]
!            in case of ion-ion collision, TL becomes ion temp.

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      INTEGER:: NS1,NS2
      REAL(rkind)   :: ANEL,TL,COULOG

      IF(NS1.EQ.1.AND.NS2.EQ.1) THEN
         COULOG=14.9D0-0.5D0*LOG(ANEL)+LOG(TL)
      ELSE
         IF(NS1.EQ.1.OR.NS2.EQ.1) THEN
            COULOG=15.2D0-0.5D0*LOG(ANEL)+LOG(TL)
         ELSE
            COULOG=17.3D0-0.5D0*LOG(ANEL)+1.5D0*LOG(TL)
         ENDIF
      ENDIF

      RETURN
      END FUNCTION COULOG

!     ***********************************************************

!           COLLISION TIME

!     ***********************************************************

!     between electrons and ions

      FUNCTION FTAUE(ANEL,ANIL,TEL,ZL)

!     ANEL : electron density [10^20 /m^3]
!     ANIL : ion density [10^20 /m^3]
!     TEL  : electron temperature [kev]
!     ZL   : ion charge number

      USE TRCOMM, ONLY : AEE, AME, EPS0, PI, PZ, RKEV, rkind
      IMPLICIT NONE
      REAL(rkind) :: ANEL, ANIL, TEL, ZL, FTAUE
      REAL(rkind) :: COEF, COULOG

      COEF = 6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)/(AEE**4*1.D20)
      IF(ZL-PZ(2).LE.1.D-7) THEN
         FTAUE = COEF*(TEL*RKEV)**1.5D0/(ANIL*ZL**2*COULOG(1,2,ANEL,TEL))
      ELSE
!     If the plasma contains impurities, we need to consider the
!     effective charge number instead of ion charge number.
!     From the definition of Zeff=sum(n_iZ_i^2)/n_e,
!     n_iZ_i^2 is replaced by n_eZ_eff at the denominator of tau_e.
         FTAUE = COEF*(TEL*RKEV)**1.5D0/(ANEL*ZL*COULOG(1,2,ANEL,TEL))
      ENDIF

      RETURN
      END FUNCTION FTAUE

!     between ions and ions

      FUNCTION FTAUI(ANEL,ANIL,TIL,ZL,PAL)

!     ANEL : electron density [10^20 /m^3]
!     ANIL : ion density [10^20 /m^3]
!     TIL  : ion temperature [kev]
!     ZL   : ion charge number
!     PAL  : ion atomic number

      USE TRCOMM, ONLY : AEE, AMM, EPS0, PI, RKEV, rkind
      IMPLICIT NONE
      REAL(rkind):: ANEL, ANIL, PAL, TIL, ZL, FTAUI
      REAL(rkind):: COEF, COULOG

      COEF = 12.D0*PI*SQRT(PI)*EPS0**2*SQRT(PAL*AMM)/(AEE**4*1.D20)
      FTAUI = COEF*(TIL*RKEV)**1.5D0/(ANIL*ZL**4*COULOG(2,2,ANEL,TIL))

      RETURN
      END FUNCTION FTAUI
