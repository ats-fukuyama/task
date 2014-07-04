!     $Id$
!
! ************************************************************
!
!            CALCULATION OF DW (BOUNCE AVERAGED)
!
! ************************************************************
!
      MODULE fpcalwm

      USE fpcomm
      USE plprof
      USE fpwmin
      USE libbes,ONLY: bessjn

      contains
!
!     -----Calculate PSIN, PCOS, PSI -----
!
      SUBROUTINE FPDWRP2(NR,ETAL,RSIN,RCOS,PSIN,PCOS,PSI)
!
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      REAL(8),INTENT(IN):: ETAL,RSIN,RCOS
      REAL(8),INTENT(OUT):: PSIN,PCOS,PSI

      IF(MODELA.EQ.0) THEN
         PSI=1.D0
         PSIN=RSIN
         PCOS=RCOS
      ELSE
         PSI=(1.D0+EPSRM2(NR))/(1.D0+EPSRM2(NR)*COS(ETAL))
         PSIN=SQRT(PSI)*RSIN
         IF (RCOS.GT.0.0D0) THEN
            PCOS= SQRT(1.D0-PSI*RSIN**2)
         ELSE
            PCOS=-SQRT(1.D0-PSI*RSIN**2)
         END IF
      ENDIF
      RETURN
      END SUBROUTINE FPDWRP2

!-------------------------------------------------------

      SUBROUTINE FP_CALWM(NSA)

      IMPLICIT NONE
      real(8),DIMENSION(NSBMAX):: sum11,sum12,sum13,sum14,sum15,sum16
      integer:: NSA, NR, NRDO, NTH, NP, NSBA
      real(8):: FACT, ETA, DWPPS, DWPTS, DWTPS, DWTTS, RHOL, P
      real(8):: DWPPL, DWPTL, DWTPL, DWTTL
      real(4):: gut, gut1, gut2
!
! =============  CALCULATION OF DWPP AND DWPT  ===============
!
      FACT=0.5D0
      NSBA=NSB_NSA(NSA)

      DO NRDO=NRSTART,NREND
         NR=NRDO
         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            IF(NTH.eq.ITL(NR).or.NTH.eq.ITU(NR))THEN
!
            ELSE
               CALL FPSUMV(NTH,NP,NR,NSA,0,DWPPS,DWPTS)                     
               DWPP(NTH,NP,NR,NSA)=DWPPS
               DWPT(NTH,NP,NR,NSA)=DWPTS
            END IF
         END DO
         END DO

         IF(MODELA.EQ.1) THEN
            DO NP=1,NPMAX+1
               DO NTH=ITL(NR)+1,NTHMAX/2
                  DWPP(NTH,NP,NR,NSA)  =(DWPP(NTH,NP,NR,NSA)           &
                       +DWPP(NTHMAX-NTH+1,NP,NR,NSA))*0.5D0
                  DWPT(NTH,NP,NR,NSA)  =(DWPT(NTH,NP,NR,NSA)           &
                       +DWPT(NTHMAX-NTH+1,NP,NR,NSA))*0.5D0
                  DWPP(NTHMAX-NTH+1,NP,NR,NSA)  =DWPP(NTH,NP,NR,NSA)
                  DWPT(NTHMAX-NTH+1,NP,NR,NSA)  =DWPT(NTH,NP,NR,NSA)
               END DO
               DWPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)*0.25D0        &
                    *( DWPP(ITL(NR)-1,NP,NR,NSA)            &
                                     /RLAMDA(ITL(NR)-1,NR)  &
                      +DWPP(ITL(NR)+1,NP,NR,NSA)            &
                                     /RLAMDA(ITL(NR)+1,NR)  &
                      +DWPP(ITU(NR)-1,NP,NR,NSA)            &
                                     /RLAMDA(ITU(NR)-1,NR)  &
                      +DWPP(ITU(NR)+1,NP,NR,NSA)            &
                                     /RLAMDA(ITU(NR)+1,NR))
               
               DWPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)*0.25D0        &
                    *( DWPT(ITL(NR)-1,NP,NR,NSA)            &
                                     /RLAMDA(ITL(NR)-1,NR)  &
                      +DWPT(ITL(NR)+1,NP,NR,NSA)            &
                                     /RLAMDA(ITL(NR)+1,NR)  &
                      +DWPT(ITU(NR)-1,NP,NR,NSA)            &
                                     /RLAMDA(ITU(NR)-1,NR)  &
                      +DWPT(ITU(NR)+1,NP,NR,NSA)            &
                                    /RLAMDA(ITU(NR)+1,NR))
               DWPP(ITU(NR),NP,NR,NSA)  =DWPP(ITL(NR),NP,NR,NSA)
               DWPT(ITU(NR),NP,NR,NSA)  =-DWPT(ITL(NR),NP,NR,NSA)
            END DO
         ENDIF
! Dw=0 at rho=0, 1 
         IF(RM(NR).le.2.D-2)THEN
            DO NP=1,NPMAX+1
               DO NTH=1,NTHMAX
                  DWPP(NTH,NP,NR,NSA)=DWPP(NTH,NP,NR,NSA) &
                       *0.D0
                  DWPT(NTH,NP,NR,NSA)=DWPT(NTH,NP,NR,NSA) &
                       *0.D0
               END DO
            END DO
         ELSEIF(RM(NR).le.1.D-1)THEN
            DO NP=1,NPMAX+1
               DO NTH=1,NTHMAX
                  DWPP(NTH,NP,NR,NSA)=DWPP(NTH,NP,NR,NSA) &
                       *(RM(NR)/1.D-1)**4
                  DWPT(NTH,NP,NR,NSA)=DWPT(NTH,NP,NR,NSA) &
                       *(RM(NR)/1.D-1)**4
               END DO
            END DO
         ELSEIF(RM(NR).ge.0.98D0)THEN
            DO NP=1,NPMAX+1
               DO NTH=1,NTHMAX
                  DWPP(NTH,NP,NR,NSA)=DWPP(NTH,NP,NR,NSA) &
                       *0.D0
                  DWPT(NTH,NP,NR,NSA)=DWPT(NTH,NP,NR,NSA) &
                       *0.D0
               END DO
            END DO
         ELSEIF(RM(NR).ge.9.D-1)THEN
            DO NP=1,NPMAX+1
               DO NTH=1,NTHMAX
                  DWPP(NTH,NP,NR,NSA)=DWPP(NTH,NP,NR,NSA) &
                       *( (1.D0-RM(NR))/1.D-1)**4
                  DWPT(NTH,NP,NR,NSA)=DWPT(NTH,NP,NR,NSA) &
                       *( (1.D0-RM(NR))/1.D-1)**4
               END DO
            END DO
         END IF
      END DO
!
! =============  CALCULATION OF DWTP AND DWTT  ===============
!
      DO NRDO=NRSTART,NREND
         NR=NRDO
         RHOL=RM(NR)
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX+1
               IF(NTH.NE.NTHMAX/2+1) THEN
                  CALL FPSUMV(NTH,NP,NR,NSA,1,DWTPS,DWTTS)
                  DWTP(NTH,NP,NR,NSA)=DWTPS
                  DWTT(NTH,NP,NR,NSA)=DWTTS
               ELSE
!                  P=PM(NP,NSA)
!                  ETA=ETAG(NTH,NR)
!                  CALL FPWAVV(RHOL,ETA,SING(NTH),COSG(NTH),P,SING(NTH),COSG(NTH),NSA,1, &
!                       DWTPL,DWTTL)
                  DWTP(NTH,NP,NR,NSA)=0.D0
                  DWTT(NTH,NP,NR,NSA)=0.D0
!                  DWTT(NTH,NP,NR,NSA)=DWTTL
               ENDIF
            END DO
         END DO

         IF(MODELA.EQ.1) THEN
            DO NTH=ITL(NR)+1,NTHMAX/2
            DO NP=1,NPMAX
               DWTP(NTH,NP,NR,NSA)=(DWTP(NTH,NP,NR,NSA)                &
                                   +DWTP(NTHMAX-NTH+2,NP,NR,NSA))*0.5D0
               DWTT(NTH,NP,NR,NSA)=(DWTT(NTH,NP,NR,NSA)                &
                                   +DWTT(NTHMAX-NTH+2,NP,NR,NSA))*0.5D0
               DWTP(NTHMAX-NTH+2,NP,NR,NSA)= DWTP(NTH,NP,NR,NSA)
               DWTT(NTHMAX-NTH+2,NP,NR,NSA)= DWTT(NTH,NP,NR,NSA)
            END DO
            END DO
         ENDIF

         IF(RM(NR).le.2.D-2)THEN
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX+1
                  DWTP(NTH,NP,NR,NSA)=DWTP(NTH,NP,NR,NSA) &
                       *0.D0
                  DWTT(NTH,NP,NR,NSA)=DWTT(NTH,NP,NR,NSA) &
                       *0.D0
               END DO
            END DO
         ELSEIF(RM(NR).le.1.D-1)THEN
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX+1
                  DWTP(NTH,NP,NR,NSA)=DWTP(NTH,NP,NR,NSA) &
                       *(RM(NR)/1.D-1)**4
                  DWTT(NTH,NP,NR,NSA)=DWTT(NTH,NP,NR,NSA) &
                       *(RM(NR)/1.D-1)**4
               END DO
            END DO
         ELSEIF(RM(NR).ge.0.98D0)THEN
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX+1
                  DWTP(NTH,NP,NR,NSA)=DWTP(NTH,NP,NR,NSA) &
                       *0.D0
                  DWTT(NTH,NP,NR,NSA)=DWTT(NTH,NP,NR,NSA) &
                       *0.D0
               END DO
            END DO
         ELSEIF(RM(NR).ge.9.D-1)THEN
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX+1
                  DWTP(NTH,NP,NR,NSA)=DWTP(NTH,NP,NR,NSA) &
                       *( (1.D0-RM(NR))/1.D-1)**4
                  DWTT(NTH,NP,NR,NSA)=DWTT(NTH,NP,NR,NSA) &
                       *( (1.D0-RM(NR))/1.D-1)**4
               END DO
            END DO
         END IF
      END DO

!      WRITE(*,*)"NCM",NCMIN(NSa),NCMAX(NSA)
      RETURN
      END SUBROUTINE FP_CALWM
!
! =======================================================
!
!      SUBROUTINE FPSUMV(ETA,RSIN,RCOS,P,NR,DWPPS,DWPTS,DWTPS,DWTTS,NSA)
      SUBROUTINE FPSUMV(NTH,NP,NR,NSA,NSW_PT,DWPS,DWTS)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSW_PT ! =0, dwpp, dwpt; =1, dwtp, dwtt
      INTEGER,INTENT(IN):: NTH, NP, NR, NSA
      DOUBLE PRECISION,INTENT(OUT):: DWPS, DWTS
      DOUBLE PRECISION:: DWPL, DWTL
      real(8):: ETA, RSIN, RCOS, P, DWPPS, DWPTS, DWTPS, DWTTS
      REAL(8):: DELH, RHOL, ETAL, PSIN, PCOS, BMIN, PSI
      REAL(8):: DWPPL, DWPTL, DWTPL, DWTTL
      integer:: N
      real(4):: gut2, gut1, gut

      IF(NSW_PT.eq.0)THEN
         DELH=2.D0*ETAM(NTH,NR)/NAVMAX
         RCOS=COSM(NTH)
         RSIN=SINM(NTH)
         P=PG(NP,NSA)
      ELSEIF(NSW_PT.eq.1)THEN
         DELH=2.D0*ETAG(NTH,NR)/NAVMAX
         RCOS=COSG(NTH)
         RSIN=SING(NTH)
         P=PM(NP,NSA)
      END IF
      DWPS=0.D0
      DWTS=0.D0
      RHOL=RM(NR)
      
      CALL PL_BMIN(RHOL,BMIN)

      DO N=1,NAVMAX
         ETAL=DELH*(N-0.5D0)
         CALL FPDWRP2(NR,ETAL,RSIN,RCOS,PSIN,PCOS,PSI)
         CALL FPWAVV(RHOL,ETAL,PSIN,PCOS,P,RSIN,RCOS,NSA,NSW_PT,DWPL,DWTL)

         IF(NSW_PT.eq.0)THEN
            DWPS=DWPS+DWPL*ABS(RCOS)/ABS(PCOS)
            DWTS=DWTS+DWTL          /SQRT(PSI)
         ELSEIF(NSW_PT.eq.1)THEN
            DWPS=DWPS+DWPL          /SQRT(PSI)
            DWTS=DWTS+DWTL*PCOS/(RCOS*PSI)
         END IF
      END DO
      DWPS =  DWPS * DELH
      DWTS =  DWPS * DELH

      RETURN
      END SUBROUTINE FPSUMV
!
! ************************************************************
!
      SUBROUTINE FPWAVV(RHOL,ETAL,PSIN,PCOS,P,RSIN,RCOS, &
                        NSA,NSW_PT,DWPL,DWTL)

      IMPLICIT NONE

      INTEGER:: NJMAX
      INTEGER,INTENT(IN):: NSA, NSW_PT
      DOUBLE PRECISION,INTENT(IN):: RHOL, ETAL, PSIN, PCOS, P
      DOUBLE PRECISION,INTENT(OUT):: DWPL, DWTL
      PARAMETER(NJMAX=100)
      REAL(8),DIMENSION(0:NJMAX):: RJ,DRJ
      DOUBLE PRECISION:: DWPPL, DWPTL, DWTPL, DWTTL
      INTEGER:: NHMAX, N, NMI, NPI, NS
      REAL(8):: RKR, RKTH, RKPH, B0PH, B0TH, RW, B2, RABSE, RGAMMA
      COMPLEX(8):: CER, CETH, CEPH, CEPARA, CEPERP, CEPLUS, CEMINUS
      REAL(8):: PPARA, PPERP, VPARA, VPERP, RKPARA, RKPERP, RWC, RKW
      REAL(8):: RGZAI, DWC11, DWC12, DWC21, DWC22, RJN, RJNM, RJNP
      COMPLEX(8):: CTHETA, CEX, CEY
      REAL(8):: RTHETA, A11, A12, A21, A22, DWC, EX
      real(8):: ratioCE, ratioCE2, ratioBE, RTHETA2
      REAL(8):: RFWM, RSIN,RCOS

      NS=NS_NSA(NSA)
!      THETA0(NSA)=RTFP0(NSA)*1.D3*AEE/(AMFP(NSA)*VC*VC)

      CALL FPSETV(RHOL,ETAL,RFWM,RKR,RKTH,RKPH,CER,CETH,CEPH,NSA)
      CALL FPSETB(RHOL,ETAL,B0PH,B0TH)
!     BOPH=B_\phi=B_toroidal, B0TH=B_\theta=B_poloidal

      RW     =2.D0*PI*RFWM*1.D6
      B2     =B0TH**2+B0PH**2
!      CEPARA =(B0TH*CETH+B0PH*CEPH)/SQRT(B2)
!      CEPERP =(B0PH*CETH-B0TH*CEPH)/SQRT(B2)
      CEPARA =CEPH
      CEPERP =CETH

      RABSE = SQRT(ABS(CEPARA)**2+ABS(CEPERP)**2)
      CEPLUS =(CER+CI*CEPERP)/SQRT(2.D0)
      CEMINUS=(CER-CI*CEPERP)/SQRT(2.D0)

      RGAMMA =SQRT(1.D0+P*P*THETA0(NSA))
      PPARA  =PTFP0(NSA)*P*PCOS
      PPERP  =PTFP0(NSA)*P*PSIN
      VPARA  =PTFP0(NSA)*P*PCOS/(AMFP(NSA)*RGAMMA)
      VPERP  =PTFP0(NSA)*P*PSIN/(AMFP(NSA)*RGAMMA)
      RKPARA =(B0TH*RKTH+B0PH*RKPH)/SQRT(B2)
      RKPERP =(B0PH*RKTH-B0TH*RKPH)/SQRT(B2)
      RWC    =AEFP(NSA)*SQRT(B2)/AMFP(NSA)
 
      RKW  =RKPARA/RW
      RGZAI=RKPERP*VPERP/ABS(RWC)
      NHMAX=MAX(ABS(NCMIN(NS)-1),ABS(NCMAX(NS)+1),2)

      CALL BESSJN(RGZAI,NHMAX,RJ,DRJ)

      DWC11=0.D0
      DWC12=0.D0
      DWC21=0.D0
      DWC22=0.D0

      DO N=NCMIN(NS),NCMAX(NS)
         IF (N.LT.0) THEN
             RJN=(-1)**(-N)*RJ(-N)
         ELSE
             RJN=RJ(N)
         ENDIF
         NMI=N-1
         IF (NMI.LT.0) THEN
             RJNM=(-1)**(-NMI)*RJ(-NMI)
         ELSE
             RJNM=RJ(NMI)
         ENDIF
         NPI=N+1
         IF (NPI.LT.0) THEN
             RJNP=(-1)**(-NPI)*RJ(-NPI)
         ELSE
             RJNP=RJ(NPI)
         ENDIF
         IF (N.EQ.0) THEN ! CTHETA*VPERP
            CTHETA=VPERP                                 &
                  *(CEPLUS*RJNM+CEMINUS*RJNP)/SQRT(2.D0) &
                  +CEPARA*VPARA*RJN
         ELSE
            CTHETA=(CEPLUS*RJNM+CEMINUS*RJNP)/SQRT(2.D0) &
                    +CEPARA*VPARA*(RJNM+RJNP)*RKPERP     &
                    /(2*N*RWC)
         ENDIF
         RTHETA2=ABS(CTHETA)**2

! coordinate transform (ppara, pperp) -> (p_0, theta_0) not (p, theta)
         IF(N.EQ.0)THEN
            IF(NSW_PT.eq.0)THEN
               A11= (RKW*RCOS)**2  *RTHETA2
               A12=-RKW**2*RCOS*RSIN*RTHETA2
            ELSEIF(NSW_PT.eq.1)THEN
               A21=-RKW**2*RCOS*RSIN*RTHETA2
               A22= (RKW*RSIN)**2  *RTHETA2
            END IF
         ELSE
            IF(NSW_PT.eq.0)THEN
               A11=(N*RWC/RW/RGAMMA*RSIN + RKW*VPERP*RCOS )**2*RTHETA2
               A12=(N*RWC/RW/RGAMMA*RSIN + RKW*VPERP*RCOS )          &
                    *(N*RWC/RW/RGAMMA*RCOS - RKW*VPERP*RSIN )*RTHETA2
            ELSEIF(NSW_PT.eq.1)THEN
               A21=(N*RWC/RW/RGAMMA*RSIN + RKW*VPERP*RCOS )          &
                    *(N*RWC/RW/RGAMMA*RCOS - RKW*VPERP*RSIN )*RTHETA2
               A22=(N*RWC/RW/RGAMMA*RCOS - RKW*VPERP*RSIN )**2*RTHETA2
            END IF
         ENDIF

         IF(VPARA.EQ.0.D0) THEN
            DWC=0.D0  
         ELSE
! DELNPR = sqrt(2)*sigma
            EX=-( (RW-RKPARA*VPARA-N*RWC/RGAMMA) &
                 /(DELNPR*RW) )**2
            IF (EX.LT.-100.D0) THEN 
                DWC=0.D0
            ELSE
                DWC=0.5D0*SQRT(PI)*AEFP(NSA)**2*EXP(EX)/PTFP0(NSA)**2 &
                    /ABS(RW)/DELNPR
            ENDIF
         ENDIF
         IF(NSW_PT.eq.0)THEN
            DWC11=DWC11+DWC*A11
            DWC12=DWC12+DWC*A12
         ELSEIF(NSW_PT.eq.1)THEN
            DWC21=DWC21+DWC*A21
            DWC22=DWC22+DWC*A22
         END IF
      END DO

      IF(NSW_PT.eq.0)THEN
         DWPL = DWC11
         DWTL = DWC12
      ELSEIF(NSW_PT.eq.1)THEN
         DWPL = DWC21
         DWTL = DWC22
      END IF

      RETURN
      END SUBROUTINE FPWAVV
!
!     ****** CALCULATE LOCAL MAGNETIC FIELD ******
!
      SUBROUTINE FPSETB(RHOL,ETAL,BTL,BPL)
!
!     BTL=B_toroidal, BPL=B_poloidal
!
      IMPLICIT NONE
      REAL(8):: RHOL, ETAL, BTL, BPL
      REAL(8):: RS, X, RS1, RS2, QL, RRMINL, RRMAXL

      IF(MODELG.EQ.2) THEN
         RS=RSRHON(RHOL)
         X=RS*COS(ETAL)
         BTL=BB/(1.D0+X/RR)
         CALL PL_QPRF(RHOL,QL)
         BPL=RS*BTL/((RR+X)*QL)
      ELSEIF(MODELG.EQ.3) THEN
         CALL PL_RRMX(RHOL,RRMINL,RRMAXL)
         RS1=RRMAXL-RR
         RS2=RR-RRMINL
         RS=0.5D0*(RS1+RS2)+0.5D0*(RS1-RS2)*COS(ETAL)
         X=RS*COS(ETAL)
         CALL PL_QPRF(RHOL,QL)
         BTL=BB/(1.D0+X/RR)
         BPL=RS*BTL/((RR+X)*QL)
      ENDIF
      RETURN
      END SUBROUTINE FPSETB
!
!     ****** CALCULATE LOCAL WAVE ELECTRIC FIELD ******
!
      SUBROUTINE FPSETV(RHOL,ETAL,RFWM,RKR,RKTH,RKPH,CER,CETH,CEPH,NSA)

      IMPLICIT NONE
      REAL(8):: RHOL, ETAL, RFWM, RKR, RKTH, RKPH
      COMPLEX(8):: CER, CETH, CEPH
      COMPLEX(8):: CEWR1, CEWTH1, CEWPH1, CKWR1, CKWTH1, CKWPH1
      INTEGER:: NSA, IERR
      REAL(8):: Y, ARG, FACT, PHL

      IF(MODELW(NS_NSA(NSA)).EQ.3) THEN
         Y=RHOL*RA*SIN(ETAL)

         DREWY=1.D1
         REWY=0.D0

         ARG=(Y-REWY)**2/DREWY**2
         IF(ARG.GT.100.D0) THEN
            FACT=0.D0
         ELSE
            FACT=EXP(-ARG)
         ENDIF
         RFWM=RFDW
         CER= CEWR*FACT
         CETH=CEWTH*FACT
         CEPH=CEWPH*FACT
         RKR= RKWR
         RKTH=RKWTH
         RKPH=RKWPH
      ELSE
         PHL=0.D0
         CALL FPWMGET(RHOL,ETAL,PHL, RFWM,CEWR1,CEWTH1,CEWPH1, &
                                     CKWR1,CKWTH1,CKWPH1,IERR)
         CER =SQRT(PWAVE)*CEWR1
         CETH=SQRT(PWAVE)*CEWTH1
         CEPH=SQRT(PWAVE)*CEWPH1
         RKR = DBLE(CKWR1)
         RKTH= DBLE(CKWTH1)
         RKPH= DBLE(CKWPH1)
      ENDIF

      RETURN
      END SUBROUTINE FPSETV
!------------------------------
    END MODULE fpcalwm
