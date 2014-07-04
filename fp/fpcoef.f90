!     $Id$
!
! ************************************************************
!
!                      CALCULATION OF D AND F
!
! ************************************************************
!
      MODULE fpcoef

      USE fpcomm
      USE fpcalc
      USE fpcalw
      USE fpcalwm
      USE fpcalwr
      USE libbes,ONLY: besekn
!      USE fpcaleind
      USE libmtx

      contains
!-------------------------------------------------------------
      SUBROUTINE FP_COEF(NSA)

      IMPLICIT NONE
      integer:: NSA, NR, NTH, NP, NS
      real(kind8):: FPMAX
      integer:: NCONST_RF
!      real(kind8),dimension(NTHMAX,NPMAX+1,NRSTART:NREND):: DWPPM, DWPTM
!      real(kind8),dimension(NTHMAX+1,NPMAX,NRSTART:NREND):: DWTPM, DWTTM
      real(kind8):: DWTTEC, DWTTIC, DWTPEC, DWTPIC

      ISAVE=0
      NS=NS_NSA(NSA)

      DO NR=NRSTART,NREND
         DO NP=NPSTART,NPENDWG
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=0.D0
            DPT(NTH,NP,NR,NSA)=0.D0
            FPP(NTH,NP,NR,NSA)=0.D0
            FEPP(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
         DO NP=NPSTARTW,NPENDWM
         DO NTH=1,NTHMAX+1
            DTP(NTH,NP,NR,NSA)=0.D0
            DTT(NTH,NP,NR,NSA)=0.D0
            FTH(NTH,NP,NR,NSA)=0.D0
            FETH(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
      ENDDO

!
!     ----- Parallel electric field accleration term -----
!
!      IF(E0.ne.0.D0)THEN
         CALL FP_CALE(NSA)
!      END IF
!
!     ----- Quasi-linear wave-particle interaction term -----

!!!   reduce the number of call of fp_calwm and fp_calw
!     include 2 IF
      IF(N_IMPL.eq.0)THEN ! N_IMPL=0

!     ----- Initialize ------------------------------------- 
      DO NR=NRSTART,NREND
         DO NP=NPSTART,NPENDWG
           DO NTH=1,NTHMAX
               DWECPP(NTH,NP,NR,NSA)=0.D0
               DWECPT(NTH,NP,NR,NSA)=0.D0
            END DO
         END DO
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               DWECTP(NTH,NP,NR,NSA)=0.D0
               DWECTT(NTH,NP,NR,NSA)=0.D0
            END DO
         END DO
      END DO

!     ECRF 
      IF(DEC.ne.0.and.NSA.eq.1) THEN
         CALL FP_CALW(NSA)
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPENDWG
               DO NTH=1,NTHMAX
                  DWECPP(NTH,NP,NR,NSA)=DWPP(NTH,NP,NR,NSA)
                  DWECPT(NTH,NP,NR,NSA)=DWPT(NTH,NP,NR,NSA)
                  DWPP(NTH,NP,NR,NSA)=0.D0
                  DWPT(NTH,NP,NR,NSA)=0.D0
               END DO
            END DO
            DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX+1
                  DWECTP(NTH,NP,NR,NSA)=DWTP(NTH,NP,NR,NSA)
                  DWECTT(NTH,NP,NR,NSA)=DWTT(NTH,NP,NR,NSA)
                  DWTP(NTH,NP,NR,NSA)=0.D0
                  DWTT(NTH,NP,NR,NSA)=0.D0
               END DO
            END DO
         END DO
      END IF

!     ICRF
      IF(MODELW(NS).EQ.1) THEN
         CALL FP_CALWR(NSA)
      ELSEIF(MODELW(NS).EQ.2) THEN
         CALL FP_CALWR(NSA)
      ELSEIF(MODELW(NS).EQ.3) THEN
         CALL FP_CALWM(NSA)
      ELSEIF(MODELW(NS).EQ.4) THEN
         CALL FP_CALWM(NSA)
      ELSEIF(MODELW(NS).ne.0) THEN
         IF(nrank.eq.0) WRITE(6,*) 'XX UNKNOWN MODELW =',MODELW(NS)
      ENDIF

      DO NR=NRSTART,NREND
         DO NP=NPSTART,NPENDWG
            DO NTH=1,NTHMAX
               DWICPP(NTH,NP,NR,NSA)=DWPP(NTH,NP,NR,NSA)
               DWICPT(NTH,NP,NR,NSA)=DWPT(NTH,NP,NR,NSA)
               DWPP(NTH,NP,NR,NSA)=DWPP(NTH,NP,NR,NSA) + DWECPP(NTH,NP,NR,NSA)
               DWPT(NTH,NP,NR,NSA)=DWPT(NTH,NP,NR,NSA) + DWECPT(NTH,NP,NR,NSA)
            END DO
         END DO
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               DWTP(NTH,NP,NR,NSA)=DWTP(NTH,NP,NR,NSA) + DWECTP(NTH,NP,NR,NSA)
               DWTT(NTH,NP,NR,NSA)=DWTT(NTH,NP,NR,NSA) + DWECTT(NTH,NP,NR,NSA)
            END DO
         END DO
      END DO

!     POOL coef DW in order to reduce the number of call fp_calwm
      DO NR=NRSTART, NREND
         DO NP=NPSTART,NPENDWG
            DO NTH=1,NTHMAX
               DWPP_P(NTH,NP,NR,NSA) = DWPP(NTH,NP,NR,NSA)
               DWPT_P(NTH,NP,NR,NSA) = DWPT(NTH,NP,NR,NSA)
               DWICPP_P(NTH,NP,NR,NSA) = DWICPP(NTH,NP,NR,NSA)
               DWICPT_P(NTH,NP,NR,NSA) = DWICPT(NTH,NP,NR,NSA)
               DWECPP_P(NTH,NP,NR,NSA) = DWECPP(NTH,NP,NR,NSA)
               DWECPT_P(NTH,NP,NR,NSA) = DWECPT(NTH,NP,NR,NSA)
            END DO
         END DO
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               DWTP_P(NTH,NP,NR,NSA) = DWTP(NTH,NP,NR,NSA)
               DWTT_P(NTH,NP,NR,NSA) = DWTT(NTH,NP,NR,NSA)
!              IC, EC is not updated
            END DO
         END DO
      END DO

      END IF ! N_IMPL=0

      IF(N_IMPL.ne.0)THEN ! N_IMPL!=0
      DO NR=NRSTART, NREND
         DO NP=NPSTART,NPENDWG
            DO NTH=1,NTHMAX
               DWPP(NTH,NP,NR,NSA) = DWPP_P(NTH,NP,NR,NSA)
               DWPT(NTH,NP,NR,NSA) = DWPT_P(NTH,NP,NR,NSA)
               DWICPP(NTH,NP,NR,NSA) = DWICPP_P(NTH,NP,NR,NSA)
               DWICPT(NTH,NP,NR,NSA) = DWICPT_P(NTH,NP,NR,NSA)
               DWECPP(NTH,NP,NR,NSA) = DWECPP_P(NTH,NP,NR,NSA)
               DWECPT(NTH,NP,NR,NSA) = DWECPT_P(NTH,NP,NR,NSA)
            END DO
         END DO
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               DWTP(NTH,NP,NR,NSA) = DWTP_P(NTH,NP,NR,NSA)
               DWTT(NTH,NP,NR,NSA) = DWTT_P(NTH,NP,NR,NSA)
!              IC, EC is not updated
            END DO
         END DO
      END DO
      END IF ! N_IMPL!=0
!!!   end of the reduction of the number of calling fp_calwm and fp_calw


!     N_IMPL = 0 means initial state in fpprep
      IF(N_IMPL.ne.0) CALL FPWAVE_CONST
!     ----- Constant Dw
      NCONST_RF=3
      IF(MODELW(NSA).eq.4.and.NCONST_RF.eq.2.and.N_IMPL.ne.0)THEN ! TOTAL Pabs(r) invariant
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPENDWG
               DO NTH=1,NTHMAX
                  IF(RPW_IMPL(NR,NSA,N_IMPL).gt.0.D0)THEN
                     DWPP(NTH,NP,NR,NSA)=DWPP(NTH,NP,NR,NSA)*RPW_INIT(NR,NSA)/RPW_IMPL(NR,NSA,N_IMPL)
                     DWPT(NTH,NP,NR,NSA)=DWPT(NTH,NP,NR,NSA)*RPW_INIT(NR,NSA)/RPW_IMPL(NR,NSA,N_IMPL)
                     DWECPP(NTH,NP,NR,NSA)=DWECPP(NTH,NP,NR,NSA)*RPW_INIT(NR,NSA)/RPW_IMPL(NR,NSA,N_IMPL)
                     DWECPT(NTH,NP,NR,NSA)=DWECPT(NTH,NP,NR,NSA)*RPW_INIT(NR,NSA)/RPW_IMPL(NR,NSA,N_IMPL)
                  END IF
               END DO
            END DO
            DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX+1
                  IF(RPW_IMPL(NR,NSA,N_IMPL).gt.0.D0)THEN
                     DWTP(NTH,NP,NR,NSA)=DWTP(NTH,NP,NR,NSA)*RPW_INIT(NR,NSA)/RPW_IMPL(NR,NSA,N_IMPL)
                     DWTT(NTH,NP,NR,NSA)=DWTT(NTH,NP,NR,NSA)*RPW_INIT(NR,NSA)/RPW_IMPL(NR,NSA,N_IMPL)
                  END IF
               END DO
            END DO
         END DO
      ELSEIF(MODELW(NSA).eq.4.and.NCONST_RF.eq.3.and.N_IMPL.ne.0)THEN ! Pabs_EC(r), Pabs_IC(r) invariant
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPENDWG
               DO NTH=1,NTHMAX
                  IF(RPWEC_IMPL(NR,NSA,N_IMPL).gt.0.D0)THEN
                     DWECPP(NTH,NP,NR,NSA)=DWECPP(NTH,NP,NR,NSA)*RPWEC_INIT(NR,NSA)/RPWEC_IMPL(NR,NSA,N_IMPL)
                     DWECPT(NTH,NP,NR,NSA)=DWECPT(NTH,NP,NR,NSA)*RPWEC_INIT(NR,NSA)/RPWEC_IMPL(NR,NSA,N_IMPL)
                  END IF
                  IF(RPWIC_IMPL(NR,NSA,N_IMPL).gt.0.D0)THEN
                     DWICPP(NTH,NP,NR,NSA)=DWICPP(NTH,NP,NR,NSA)*RPWIC_INIT(NR,NSA)/RPWIC_IMPL(NR,NSA,N_IMPL)
                     DWICPT(NTH,NP,NR,NSA)=DWICPT(NTH,NP,NR,NSA)*RPWIC_INIT(NR,NSA)/RPWIC_IMPL(NR,NSA,N_IMPL)
                  END IF
                  DWPP(NTH,NP,NR,NSA)=DWECPP(NTH,NP,NR,NSA)+DWICPP(NTH,NP,NR,NSA)
                  DWPT(NTH,NP,NR,NSA)=DWECPT(NTH,NP,NR,NSA)+DWICPT(NTH,NP,NR,NSA)
               END DO
            END DO
            DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX+1
                  DWTPEC=0.D0
                  DWTPIC=0.D0
                  DWTTEC=0.D0
                  DWTTIC=0.D0
                  IF(RPWEC_IMPL(NR,NSA,N_IMPL).gt.0.D0)THEN
                     DWTPEC = DWECTP(NTH,NP,NR,NSA)*RPWEC_INIT(NR,NSA)/RPWEC_IMPL(NR,NSA,N_IMPL) 
                     DWTTEC = DWECTT(NTH,NP,NR,NSA)*RPWEC_INIT(NR,NSA)/RPWEC_IMPL(NR,NSA,N_IMPL) 
                  END IF
                  IF(RPWIC_IMPL(NR,NSA,N_IMPL).gt.0.D0)THEN
                     DWTPIC =( DWTP(NTH,NP,NR,NSA)-DWECTP(NTH,NP,NR,NSA) ) &
                          *RPWIC_INIT(NR,NSA)/RPWIC_IMPL(NR,NSA,N_IMPL) 
                     DWTTIC =( DWTT(NTH,NP,NR,NSA)-DWECTT(NTH,NP,NR,NSA) ) &
                          *RPWIC_INIT(NR,NSA)/RPWIC_IMPL(NR,NSA,N_IMPL) 
                  END IF
                  DWTP(NTH,NP,NR,NSA)=DWTPEC+DWTPIC
                  DWTT(NTH,NP,NR,NSA)=DWTTEC+DWTTIC
               END DO
            END DO
         END DO
      END IF
!     ----- Collisional slowing down and diffusion term -----
      CALL FP_CALC(NSA)
!     ----- Sum up velocity diffusion terms -----

!      IF(NRANK.eq.0) open(9,file='FE.dat')
      DO NR=NRSTART,NREND
         DO NP=NPSTART,NPENDWG
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=DCPP(NTH,NP,NR,NSA)+DWPP(NTH,NP,NR,NSA)
            DPT(NTH,NP,NR,NSA)=DCPT(NTH,NP,NR,NSA)+DWPT(NTH,NP,NR,NSA)
            FPP(NTH,NP,NR,NSA)=FEPP(NTH,NP,NR,NSA)+FCPP(NTH,NP,NR,NSA)
!            IF(NRANK.eq.0)THEN
!               WRITE(9,'(4E16.8)') PG(NP,NSA)*COSM(NTH), PG(NP,NSA)*SINM(NTH), &
!                 DCPP(NTH,NP,NR,NSA),FCPP(NTH,NP,NR,NSA)
!            ENDIF
         ENDDO
!         IF(NRANK.eq.0) WRITE(9,*) " "
!         IF(NRANK.eq.0) WRITE(9,*) " "
         ENDDO
!
         DO NP=NPSTARTW,NPENDWM
         DO NTH=1,NTHMAX+1
            DTP(NTH,NP,NR,NSA)=DCTP(NTH,NP,NR,NSA)+DWTP(NTH,NP,NR,NSA)
            DTT(NTH,NP,NR,NSA)=DCTT(NTH,NP,NR,NSA)+DWTT(NTH,NP,NR,NSA)
            FTH(NTH,NP,NR,NSA)=FETH(NTH,NP,NR,NSA)+FCTH(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
      ENDDO
!      IF(NRANK.eq.0) close(9)

!     ----- Radial diffusion term -----

      IF(MODELD.ge.1) THEN
            CALL FP_CALR2(NSA)
!      ELSEIF(MODELD.GT.0) THEN
!         CALL FP_CALR(NSA)
      END IF
!      IF(MODELD.GT.0) CALL FP_CALR(NSA)

!     ----- Particle source term -----

      CALL FP_CALS(NSA)

!
!     ****************************
!     Boundary condition at p=pmax
!     ****************************
!
      IF(NPENDWG.eq.NPMAX+1)THEN
      DO NR=NRSTART,NREND
         DO NTH=1,NTHMAX
            DPP(NTH,NPMAX+1,NR,NSA)=0.D0
            DPT(NTH,NPMAX+1,NR,NSA)=0.D0
!            FPP(NTH,NPMAX+1,NR,NSA)=0.D0
            FPP(NTH,NPMAX+1,NR,NSA)=max(0.D0,FPP(NTH,NPMAX+1,NR,NSA))
         END DO
      END DO
      END IF

      RETURN
      END SUBROUTINE fp_coef

! ****************************************
!     Parallel electric field
! ****************************************

      SUBROUTINE FP_CALE(NSA)

      IMPLICIT NONE
      integer:: NSA, NSB, NR, NTH, NP
      real(kind8):: PSP, SUML, ANGSP, SPL, FPMAX
      integer:: NG
      real(kind8):: FACT, DELH, sum11, ETAL, X, PSIB, PCOS, sum15, ARG

      DO NR=NRSTART,NREND
         DO NP=NPSTART,NPENDWG
            DO NTH=1,NTHMAX
               FEPP(NTH,NP,NR,NSA)= AEFP(NSA)*EP(NR)/PTFP0(NSA)*COSM(NTH)
            ENDDO
         ENDDO
      ENDDO

      DO NR=NRSTART,NREND
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               FETH(NTH,NP,NR,NSA)=-AEFP(NSA)*EP(NR)/PTFP0(NSA)*SING(NTH)
            ENDDO
         ENDDO
      ENDDO
      
      IF(MODELA.eq.1)THEN
         DO NR=NRSTART,NREND
            CALL FP_CALE_LAV(NR,NSA)
         ENDDO
      END IF
      
      RETURN
      END SUBROUTINE FP_CALE
! ****************************************
!     BOUNCE AVERAGING FEPP, FETH
! ****************************************
      SUBROUTINE FP_CALE_LAV(NR, NSA)

      IMPLICIT NONE
      integer,intent(in):: NR, NSA
      integer:: NTH, NP, NG
      double precision:: DELH, SUM, ETAL, X, PSIB

!     BOUNCE AVERAGE FEPP
      DO NP=NPSTART,NPENDWG
         DO NTH=1,ITL(NR)-1
            DELH=2.D0*ETAM(NTH,NR)/NAVMAX
            SUM=0.D0
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               X=EPSRM2(NR)*COS(ETAL)
               PSIB=(1.D0+EPSRM2(NR))/(1.D0+X) 
               
               SUM = SUM + FEPP(NTH,NP,NR,NSA)*PSIB
            END DO
            FEPP(NTH,NP,NR,NSA) = SUM*DELH * RCOEFNG(NR)
         END DO ! END NTH
         DO NTH=ITL(NR)+1,ITU(NR)-1
            FEPP(NTH,NP,NR,NSA)=0.D0
         END DO
         DO NTH=ITU(NR)+1,NTHMAX
            DELH=2.D0*ETAM(NTH,NR)/NAVMAX
            SUM=0.D0
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               X=EPSRM2(NR)*COS(ETAL)
               PSIB=(1.D0+EPSRM2(NR))/(1.D0+X) 
               
               SUM = SUM + FEPP(NTH,NP,NR,NSA)*PSIB
            END DO
            FEPP(NTH,NP,NR,NSA) = SUM*DELH * RCOEFNG(NR)
         END DO
         FEPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0 &
              *( FEPP(ITL(NR)-1,NP,NR,NSA)/RLAMDA(ITL(NR)-1,NR) &
              +FEPP(ITL(NR)+1,NP,NR,NSA)/RLAMDA(ITL(NR)+1,NR) &
              +FEPP(ITU(NR)-1,NP,NR,NSA)/RLAMDA(ITU(NR)-1,NR) &
              +FEPP(ITU(NR)+1,NP,NR,NSA)/RLAMDA(ITU(NR)+1,NR))
         FEPP(ITU(NR),NP,NR,NSA)=FEPP(ITL(NR),NP,NR,NSA)
      END DO ! END NP
!     BOUNCE AVERAGE FETH
      DO NP=NPSTARTW,NPENDWM
         DO NTH=1,ITL(NR)
            DELH=2.D0*ETAG(NTH,NR)/NAVMAX
            SUM=0.D0
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               X=EPSRM2(NR)*COS(ETAL)
               PSIB=(1.D0+EPSRM2(NR))/(1.D0+X)
               
               SUM = SUM + FETH(NTH,NP,NR,NSA)*PSIB
            END DO
            FETH(NTH,NP,NR,NSA) = SUM * DELH * RCOEFNG(NR)
         END DO
         DO NTH=ITL(NR)+1,ITU(NR)
            FETH(NTH,NP,NR,NSA)=FETH(NTH,NP,NR,NSA)*2.D0*ETAG(NTH,NR)*0.D0
         END DO
         DO NTH=ITU(NR)+1,NTHMAX+1
            DELH=2.D0*ETAG(NTH,NR)/NAVMAX
            SUM=0.D0
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               X=EPSRM2(NR)*COS(ETAL)
               PSIB=(1.D0+EPSRM2(NR))/(1.D0+X)
               
               SUM = SUM + FETH(NTH,NP,NR,NSA)*PSIB
            END DO
            FETH(NTH,NP,NR,NSA) = SUM * DELH * RCOEFNG(NR)
         END DO
      END DO ! END NP

      END SUBROUTINE FP_CALE_LAV

! ****************************************
!     Radial transport
! ****************************************

      SUBROUTINE FP_CALR(NSA)

      IMPLICIT NONE
      integer:: NSA, NSBA, NS, NR, NTH, NP, NG
      real(kind8):: RHON, RTFPL, FACTR, FACTP, FACTRN, FACTRT, SV
      real(kind8):: PSIB, PCOS, X, ETAL, sumd, sumf, DELH
      real(kind8):: DNDR, NEDGE, FACT, DINT_D, DINT_F, DFDR_R1, F_R1, WRL

      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)
      DO NR=NRSTART,NREND+1
!         IF(MODELD.EQ.2.OR.MODELD.EQ.3.OR.MODELD.EQ.4.OR.MODELD.EQ.6.OR.MODELD.EQ.7) THEN
         IF(MODELD.ne.1.and.MODELD.ne.5) THEN
            RTFPL=RTFP(NR,NSA)/RTFP0(NSA)
         ENDIF
         RHON=RG(NR)
         IF(MODELD.EQ.2.OR.MODELD.EQ.4.OR.MODELD.eq.5) THEN ! analytical pinch
            SV=MAX(PNS(NS)/PN(NS),1.D-3)
            FACTRN=PROFN1*PROFN2*RHON**(PROFN1-1.D0)/((1-RHON**PROFN1)+SV)
            SV=MAX(PTS(NS)/PTPP(NS),1.D-3)
            FACTRT=PROFT1*PROFT2*RHON**(PROFT1-1.D0)/((1-RHON**PROFT1)+SV)

            IF(PROFN2.ge.1.D0)THEN
               NEDGE=(RNFP0(NSA)-RNFPS(NSA))*(1.D0-RG(NR)**PROFN1)**PROFN2+RNFPS(NSA)
               DNDR=-PROFN1*RHON**(PROFN1-1.D0)          &
                    *PROFN2*( RNFP0(NSA)-RNFPS(NSA) )    &
                    *( 1.D0-RHON**PROFN1 )**(PROFN2-1.D0)
            ELSE
               IF(NR.ne.NRMAX+1)THEN
                  NEDGE=(RNFP0(NSA)-RNFPS(NSA))*(1.D0-RG(NR)**PROFN1)**PROFN2+RNFPS(NSA)
                  DNDR=-PROFN1*RHON**(PROFN1-1.D0)          &
                       *PROFN2*( RNFP0(NSA)-RNFPS(NSA) )    &
                       *( 1.D0-RHON**PROFN1 )**(PROFN2-1.D0)
               ELSEIF(NR.eq.NRMAX+1)THEN
                  NEDGE=(RNFP0(NSA)-RNFPS(NSA))*(1.D0-RG(NRMAX)**PROFN1)**PROFN2+RNFPS(NSA)
                  DNDR=( -NEDGE+RNFPS(NSA))/DELR
               END IF
            END IF
         ELSEIF(MODELD.ge.6)THEN ! numarical pinch
            DINT_D=0.D0
            DINT_F=0.D0
!            DO NP=1,NPMAX
            DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX
                  IF(NTG1.ne.1)THEN
                     WRL=WEIGHR(NTH,NP,NR,NSA)
                  ELSEIF(NTG1.eq.1.and.NR.ne.1)THEN
!                     WRL=0.5D0
                     WRL=(4.D0*NR-3.D0)/(2.D0*NR-2.D0)/4.D0
                  END IF
                  IF(NR.eq.1)THEN
                     DFDR_R1 = ( FNSP(NTH,NP,NR,NSA)-FS1(NTH,NP,NSA) ) / DELR * 2.D0
                     F_R1 = FS1(NTH,NP,NSA) 
                  ELSEIF(NR.eq.NRMAX+1)THEN ! FS2 = F_INIT at rho=1+delR/2
                     DFDR_R1 = ( FS2(NTH,NP,NSA)-FNSP(NTH,NP,NR-1,NSA) ) / DELR
                     F_R1 = ( (1.D0-WRL)*FS2(NTH,NP,NSA) + WRL*FNSP(NTH,NP,NR-1,NSA) )
                  ELSE
                     DFDR_R1 = ( FNSP(NTH,NP,NR,NSA)-FNSP(NTH,NP,NR-1,NSA) ) / DELR
                     F_R1 = ( (1.D0-WRL)*FNSP(NTH,NP,NR,NSA) + WRL*FNSP(NTH,NP,NR-1,NSA) )
                  END IF
                  IF(MODELD.eq.6)THEN
                     DINT_D = DINT_D + VOLP(NTH,NP,NSBA)/SQRT(RTFPL+PM(NP,NSBA)**2)*DFDR_R1 * RLAMDA_GG(NTH,NR)!PM?
                     DINT_F = DINT_F + VOLP(NTH,NP,NSBA)/SQRT(RTFPL+PM(NP,NSBA)**2)*F_R1 * RLAMDA_GG(NTH,NR)
                  ELSEIF(MODELD.eq.7)THEN
                     DINT_D = DINT_D + VOLP(NTH,NP,NSBA)/SQRT(1.D0+PG(NP,NSBA)/RTFPL)*DFDR_R1
                     DINT_F = DINT_F + VOLP(NTH,NP,NSBA)/SQRT(1.D0+PG(NP,NSBA)/RTFPL)*F_R1
                  ELSEIF(MODELD.eq.8)THEN
                     DINT_D = DINT_D + VOLP(NTH,NP,NSBA)/(RTFPL+PG(NP,NSBA)**2)*DFDR_R1 * RLAMDA_GG(NTH,NR)
                     DINT_F = DINT_F + VOLP(NTH,NP,NSBA)/(RTFPL+PG(NP,NSBA)**2)*F_R1 * RLAMDA_GG(NTH,NR)
                  ELSEIF(MODELD.eq.9)THEN
                     DINT_D = DINT_D + VOLP(NTH,NP,NSBA)*DFDR_R1 * RLAMDA_GG(NTH,NR)
                     DINT_F = DINT_F + VOLP(NTH,NP,NSBA)*F_R1 * RLAMDA_GG(NTH,NR)
                  END IF
               END DO
            END DO
            FACTR = DINT_D/DINT_F
         ENDIF

!         DO NP=1,NPMAX
         DO NP=NPSTART,NPEND
            SELECT CASE(MODELD)
            CASE(1)
               FACTR=0.D0
               FACTP=1.D0
            CASE(2)
               FACTR=-FACTRN+(1.5D0-0.5D0*PM(NP,NSBA)**2/RTFPL)*FACTRT
               FACTP=1.D0
            CASE(3)
               FACTR=0.D0
               FACTP=1.D0/SQRT(RTFPL+PG(NP,NSBA)**2)
            CASE(4)
               FACTR=-FACTRN+(1.5D0-0.5D0*PM(NP,NSBA)**2/RTFPL)*FACTRT
               FACTP=1.D0/SQRT(1.D0+PM(NP,NSBA)**2/RTFPL)
            CASE(5) ! case(1) with pinch, independent on p
               FACTR=DNDR/NEDGE
               FACTP=1.D0
            CASE(6) ! case(3) with pinch, depend on 1/p
               FACTP=1.D0/SQRT(RTFPL+PG(NP,NSBA)**2)
            CASE(7) ! case(3) with pinch, depend on 1/sqrt{p}
               FACTP=1.D0/SQRT(RTFPL+PG(NP,NSBA))
            CASE(8) ! case(3) with pinch, depend on 1/p^2
               FACTP=1.D0/(RTFPL+PG(NP,NSBA)**2)
            CASE(9) ! case(3) with pinch, = MODELD=5
               FACTP=1.D0
            END SELECT
            IF(NR.ne.NRMAX+1)THEN
               DO NTH=1,NTHMAX
                  FACT= (DRR0-DRRS)*(1.D0-RHON**2)+DRRS 
                  DRR(NTH,NP,NR,NSA)= FACT &
                       *FACTP      /(RA*RA)*RLAMDA_GG(NTH,NR)
                  FRR(NTH,NP,NR,NSA)= FACT &
                       *FACTP*FACTR/(RA*RA)*RLAMDA_GG(NTH,NR)
               ENDDO
            ELSE
               DO NTH=1,NTHMAX
                  FACT= (DRR0-DRRS)*(1.D0-RHON**2)+DRRS 
                  DRR(NTH,NP,NR,NSA)=0.D0 
                  FRR(NTH,NP,NR,NSA)=0.D0
               ENDDO               
            END IF
         ENDDO

         IF(MODELA.eq.1)THEN! Bounce average for radial diffusion coef.
!            DO NP=1,NPMAX
            DO NP=NPSTART,NPEND
               DO NTH=ITLG(NR)+1,NTHMAX/2
                  DRR(NTH,NP,NR,NSA) &
                       =(DRR(NTH,NP,NR,NSA) &
                       +DRR(NTHMAX-NTH+1,NP,NR,NSA))/2.D0
                  FRR(NTH,NP,NR,NSA) &
                       =(FRR(NTH,NP,NR,NSA) &
                       +FRR(NTHMAX-NTH+1,NP,NR,NSA))/2.D0
                  DRR(NTHMAX-NTH+1,NP,NR,NSA) &
                       =DRR(NTH,NP,NR,NSA)
                  FRR(NTHMAX-NTH+1,NP,NR,NSA) &
                       =FRR(NTH,NP,NR,NSA)
               END DO
            END DO
            IF(NR.eq.1)THEN
!               DO NP=1,NPMAX
               DO NP=NPSTART,NPEND
                  DRR(ITLG(NR),NP,NR,NSA) = RLAMDA_GG(ITLG(NR),NR)/4.D0 &
                       *( DRR(ITLG(NR)-1,NP,NR,NSA)/RLAMDA_GG(ITLG(NR)-1,NR) &
                       +DRR(ITLG(NR)+1,NP,NR,NSA)/RLAMDA_GG(ITLG(NR)+1,NR)   &
                       +DRR(ITUG(NR)-1,NP,NR,NSA)/RLAMDA_GG(ITUG(NR)-1,NR)   &
                       +DRR(ITUG(NR)+1,NP,NR,NSA)/RLAMDA_GG(ITUG(NR)+1,NR)  )
                  FRR(ITLG(NR),NP,NR,NSA) = RLAMDA_GG(ITLG(NR),NR)/4.D0 &
                       *( FRR(ITLG(NR)-1,NP,NR,NSA)/RLAMDA_GG(ITLG(NR)-1,NR) &
                       +FRR(ITLG(NR)+1,NP,NR,NSA)/RLAMDA_GG(ITLG(NR)+1,NR)   &
                       +FRR(ITUG(NR)-1,NP,NR,NSA)/RLAMDA_GG(ITUG(NR)-1,NR)   &
                       +FRR(ITUG(NR)+1,NP,NR,NSA)/RLAMDA_GG(ITUG(NR)+1,NR)  )
                  DRR(ITUG(NR),NP,NR,NSA)=DRR(ITLG(NR),NP,NR,NSA)
                  FRR(ITUG(NR),NP,NR,NSA)=FRR(ITLG(NR),NP,NR,NSA)
               END DO
            ELSE
               DO NP=NPSTART,NPEND
                     DRR(ITLG(NR),NP,NR,NSA) = RLAMDA_GG(ITLG(NR),NR)/4.D0 &
                          *( DRR(ITLG(NR)-1,NP,NR,NSA)/RLAMDA_GG(ITLG(NR)-1,NR) &
                          +DRR(ITLG(NR)+1,NP,NR,NSA)/RLAMDA_GG(ITLG(NR)+1,NR)   &
                          +DRR(ITUG(NR)-1,NP,NR,NSA)/RLAMDA_GG(ITUG(NR)-1,NR)   &
                          +DRR(ITUG(NR)+1,NP,NR,NSA)/RLAMDA_GG(ITUG(NR)+1,NR)  )
                     FRR(ITLG(NR),NP,NR,NSA) = RLAMDA_GG(ITLG(NR),NR)/4.D0 &
                          *( FRR(ITLG(NR)-1,NP,NR,NSA)/RLAMDA_GG(ITLG(NR)-1,NR) &
                          +FRR(ITLG(NR)+1,NP,NR,NSA)/RLAMDA_GG(ITLG(NR)+1,NR)   &
                          +FRR(ITUG(NR)-1,NP,NR,NSA)/RLAMDA_GG(ITUG(NR)-1,NR)   &
                          +FRR(ITUG(NR)+1,NP,NR,NSA)/RLAMDA_GG(ITUG(NR)+1,NR)  )
                     DRR(ITUG(NR),NP,NR,NSA)=DRR(ITLG(NR),NP,NR,NSA)
                     FRR(ITUG(NR),NP,NR,NSA)=FRR(ITLG(NR),NP,NR,NSA)
                  END DO
            END IF
         END IF! end of bounce average
      ENDDO ! NR      


      RETURN
      END SUBROUTINE FP_CALR

! ****************************************
!     Radial transport
! ****************************************

      SUBROUTINE FP_CALR2(NSA)

      USE plprof
      IMPLICIT NONE
      integer:: NSA, NSBA, NS, NR, NTH, NP, NG
      real(kind8):: RHON, RTFPL, FACTR, FACTP, SV
      real(kind8):: PSIB, PCOS, X, ETAL, sumd, sumf, DELH
      real(kind8):: DNDR, NEDGE, FACT, DINT_D, DINT_F, WRL
      real(kind8):: sum1, temp1, SRHODP, SRHODM, SRHOFP, SRHOFM
      real(kind8):: WRH, DFDR_D, DFDR_F, F_R2, DFDR_R2, F_R1, DFDR_R1
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      double precision:: densm, densp
      INTEGER:: ISW_D

      ISW_D=MOD(MODELD,5)
      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)
      DO NR=NRSTART,NRENDWG
         RHON=RG(NR)
         RTFPL=RTFP(NR,NSA)/RTFP0(NSA)
         DO NP=NPSTART,NPEND
            DO NTH=1,NTHMAX
!------------- SET P DEPENDENCE
               SELECT CASE(ISW_D)
               CASE(0) ! no transport
                  FACTP=0.D0
               CASE(1) ! no p dependence
                  FACTP=1.D0
               CASE(2) ! depend on 1/p
                  FACTP=1.D0/SQRT(RTFPL+PM(NP,NSBA)**2*SINM(NTH)**2)
               CASE(3) ! depend on 1/sqrt{p}
                  FACTP=1.D0/SQRT(RTFPL+PM(NP,NSBA)*SINM(NTH))
               CASE(4) ! depend on 1/p^2
                  FACTP=1.D0/(RTFPL+PM(NP,NSBA)**2*SINM(NTH)**2)
                END SELECT

               FACTR= (DRR0-DRRS)*(1.D0-RHON**2)+DRRS 
               DRR(NTH,NP,NR,NSA)= FACTR &
                    *FACTP      /(RA*RA)
            ENDDO
         ENDDO

! ------ SET PINCH TERM
         DINT_D=0.D0
         DINT_F=1.D0
         IF(MODELD.ge.6)THEN
            DINT_F=0.D0
            DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX
                  IF(NTG1.eq.0.and.N_IMPL.eq.0)THEN
                     IF(NR.eq.1)THEN
                        WRL=0.25D0 ! not necessary
                     ELSE
                        WRL=(4.D0*RG(NR)+DELR)/(8.D0*RG(NR))                     
                     END IF
                  ELSE
                     WRL=WEIGHR(NTH,NP,NR,NSA)
                  END IF

                  IF(NR.eq.1)THEN
                     DFDR_R1 = ( FNSP(NTH,NP,NR,NSA)-FS1(NTH,NP,NSA) ) / DELR *2.D0*0
                     F_R1 = FS1(NTH,NP,NSA)
                     SRHODM=DFDR_R1 * DRR(NTH,NP,NR,NSA)
                     SRHOFM=F_R1    * DRR(NTH,NP,NR,NSA)
                  ELSEIF(NR.eq.NRMAX+1)THEN ! FS2 = F_INIT at rho=1+delR/2
                     DFDR_R1 = ( FS2(NTH,NP,NSA)-FNSP(NTH,NP,NR-1,NSA) ) / DELR
                     F_R1 = ( (1.D0-WRL)*FS2(NTH,NP,NSA) + WRL*FNSP(NTH,NP,NR-1,NSA) )
!                     F_R1 = FS3(NTH,NP,NSA)
                     SRHODM=DFDR_R1 * DRR(NTH,NP,NR,NSA)
                     SRHOFM=F_R1    * DRR(NTH,NP,NR,NSA)
                  ELSE
                     DFDR_R1 = ( FNSP(NTH,NP,NR,NSA)-FNSP(NTH,NP,NR-1,NSA) ) / DELR
                     F_R1 = ( (1.D0-WRL)*FNSP(NTH,NP,NR,NSA) + WRL*FNSP(NTH,NP,NR-1,NSA) ) 
                     SRHODM=DFDR_R1 * DRR(NTH,NP,NR,NSA)
                     SRHOFM=F_R1    * DRR(NTH,NP,NR,NSA)
                  END IF
                  DINT_D = DINT_D + VOLP(NTH,NP,NSBA)*SRHODM
                  DINT_F = DINT_F + VOLP(NTH,NP,NSBA)*SRHOFM
               END DO
            END DO
! integration
            CALL mtx_set_communicator(comm_np) 
            CALL p_theta_integration(DINT_D) 
            CALL p_theta_integration(DINT_F) 
            CALL mtx_reset_communicator 
!         WRITE(*,'(A,2I3,2E14.6)') "DINT=", NSA,NR,DINT_D, DINT_F
         ENDIF

         FACT=DINT_D/DINT_F
         DO NP=NPSTART,NPEND
            DO NTH=1,NTHMAX
               FRR(NTH,NP,NR,NSA) = FACT * DRR(NTH,NP,NR,NSA)
            END DO
         END DO

!         IF(NR.eq.NREND+1.and.NSA.eq.1) WRITE(*,'(I4,3E16.8)') NR, DRR(1,1,NR,1), FRR(1,1,NR,1),WRL
!         IF(NPSTART.eq.1.and.NSA.eq.1) WRITE(*,'(I4,3E16.8)') NR, DRR(1,1,NR,1), FRR(1,1,NR,1), WRL

      ENDDO ! NR


      RETURN
      END SUBROUTINE FP_CALR2


! ****************************************
!     Particle source and loss
! ****************************************

      SUBROUTINE FP_CALS(NSA)

      USE fpnfrr
      IMPLICIT NONE
      integer:: NSA, NSB, NSBA, NR, NTH, NP, NS, ID
      integer:: NBEAM, NSABEAM, NSAX, ISW_LOSS
      real(kind8):: PSP, SUML, ANGSP, SPL, FL

      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)

!     ----- Particle source term -----

      DO NR=NRSTART,NREND
!         DO NP=1,NPMAX
         DO NP=NPSTART,NPEND
            DO NTH=1,NTHMAX
               PPL(NTH,NP,NR,NSA)=0.D0
               SPPB(NTH,NP,NR,NSA)=0.D0
               SPPS(NTH,NP,NR,NSA)=0.D0
               SPPD(NTH,NP,NSA)=0.D0
            ENDDO
         ENDDO
      ENDDO
!     ----- NBI source term -----

      IF(MODELA.eq.0)THEN
         CALL NBI_SOURCE_A0(NSA)
      ELSE
         CALL NBI_SOURCE_A1(NSA)
      END IF

!     ----- Fixed fusion source term -----

      IF(MODELS.EQ.1) THEN
         IF(NSSPF.EQ.NS_NSA(NSA)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*SPFENG*AEE)/PTFP0(NSA)
            SUML=0.D0
!            DO NP=1,NPMAX
            DO NP=NPSTART,NPEND
               IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                  DO NR=1,NRMAX
                     SPL=EXP(-(RM(NR)-SPFR0)**2/SPFRW**2)
                     DO NTH=1,NTHMAX
                        SUML=SUML &
                            +SPL*VOLP(NTH,NP,NSBA)*VOLR(NR)!*RLAMDAG(NTH,NR)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
            SUML=SUML*RNFP0(NSA)
!            DO NP=1,NPMAX
            DO NP=NPSTART,NPEND
               IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                  DO NR=NRSTART,NREND
                     SPL=EXP(-(RM(NR)-SPFR0)**2/SPFRW**2)
                     DO NTH=1,NTHMAX
                        SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                             + SPFTOT*SPL/SUML!*RLAMDA(NTH,NR)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         END IF

!     ----- electron -----

         IF(NS_NSA(NSA).EQ.1) THEN
            NSABEAM=0
            DO NSAX=1,NSAMAX
               IF(NS_NSA(NSAX).EQ.NSSPF) NSABEAM=NSAX
            ENDDO
            IF(NSABEAM.NE.0) THEN
            PSP=SQRT(2.D0*AMFP(NSA)**2*SPFENG*AEE &
                    /(AMFP(NSABEAM)))*PTFP0(NSA)
            SUML=0.D0
!            DO NP=1,NPMAX
            DO NP=NPSTART,NPEND
               IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                  DO NR=1,NRMAX
                     SPL=EXP(-(RM(NR)-SPFR0)**2/SPFRW**2)
                     DO NTH=1,NTHMAX
                        SUML=SUML &
                            +SPL*VOLP(NTH,NP,NSBA)*VOLR(NR)!*RLAMDAG(NTH,NR)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
            SUML=SUML*RNFP0(NSA)
!            DO NP=1,NPMAX
            DO NP=NPSTART,NPEND
               IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                  DO NR=NRSTART,NREND
                     SPL=EXP(-(RM(NR)-SPFR0)**2/SPFRW**2)
                     DO NTH=1,NTHMAX
                        SPPB(NTH,NP,NR,NSA)=SPPB(NTH,NP,NR,NSA) &
                             +SPFTOT*SPL/SUML!*RLAMDA(NTH,NR)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
            ENDIF
         END IF
      END IF ! MODELS=1

!     ----- Calcluated fusion source term -----
      IF(MODELS.EQ.2) THEN
         IF(MODELA.eq.0)THEN
            CALL FUSION_SOURCE_S2A0(NSA)
         ELSE
            CALL FUSION_SOURCE_S2A1(NSA)            
         END IF

      ENDIF ! MODELS=2
!
!     ----- Particle loss and source terms -----
!

      ISW_LOSS=1
      IF(TLOSS(NS).EQ.0.D0) THEN
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX
                  PPL(NTH,NP,NR,NSA)=0.D0
               ENDDO
            ENDDO
         ENDDO
      ELSE
         IF(ISW_LOSS.eq.0)THEN
            DO NR=NRSTART,NREND
               DO NP=NPSTART,NPEND
                  FL=FPMXWL(PM(NP,NSBA),NR,NS) 
                  DO NTH=1,NTHMAX
                     PPL(NTH,NP,NR,NSA)=-1.D0/TLOSS(NS)!*RLAMDA(NTH,NR)
                     SPPS(NTH,NP,NR,NSA)= FL /TLOSS(NS)!*RLAMDA(NTH,NR)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO NR=NRSTART,NREND
               DO NP=NPSTART,NPEND
                  FL=FPMXWL_LT(PM(NP,NSBA),NR,NS)
                  DO NTH=1,NTHMAX
                     PPL(NTH,NP,NR,NSA)=-1.D0/TLOSS(NS)
                     SPPS(NTH,NP,NR,NSA)=FL/TLOSS(NS)
                  ENDDO
!                  IF(NRANK.eq.0.and.N_IMPL.eq.1) &
!                       WRITE(*,'(2I3,3E14.6)') NSA,NP,SPPS(1,NP,NR,NSA),FNSP(1,NP,NR,NSA),FNSM(1,np,nr,nsa)
               ENDDO
            ENDDO
         END IF
      ENDIF

      RETURN
      END SUBROUTINE FP_CALS

! ****************************************
!     MAXWELLIAN VELOCITY DISTRIBUTION
! ****************************************

      FUNCTION FPMXWL(PML,NR,NS)

      USE plprof
      implicit none
      integer :: NR, NS
      real(kind8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      real(kind8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      TYPE(pl_plf_type),DIMENSION(NSMAX):: plf
      real(kind8):: FPMXWL

      AMFDL=PA(NS)*AMP
      AEFDL=PZ(NS)*AEE
      RNFD0L=PN(NS)
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0L=SQRT(RTFD0L*1.D3*AEE*AMFDL)

      IF(NR.eq.0)THEN
         RL=0.D0
         RHON=ABS(RL)
      ELSEIF(NR.EQ.NRMAX+1) THEN
         RL=RM(NRMAX)+DELR
         RHON=MIN(RL,1.D0)
      ELSE
         RL=RM(NR)
         RHON=RL
      ENDIF
      CALL PL_PROF(RHON,PLF)
      RNFDL=PLF(NS)%RN/RNFD0L
      RTFDL=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0

      IF(MODELR.EQ.0) THEN
         FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0L)**3
         EX=PML**2/(2.D0*RTFDL/RTFD0L)
         FPMXWL=FACT*EXP(-EX)
      ELSE
         THETA0L=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
         THETAL=THETA0L*RTFDL/RTFD0L
         Z=1.D0/THETAL
         DKBSL=BESEKN(2,Z)
         FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL) &
              *RTFD0L
         EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
         FPMXWL=FACT*EXP(EX)
      END IF

      RETURN
      END FUNCTION FPMXWL
!-------------------------------------------------------------
      FUNCTION FPMXWL_S(PML,NR,NS)

      USE plprof
      implicit none
      integer :: NR, NS
      real(kind8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      real(kind8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      TYPE(pl_plf_type),DIMENSION(NSMAX):: plf
      real(kind8):: FPMXWL_S

      AMFDL=PA(NS)*AMP
      AEFDL=PZ(NS)*AEE
      RNFD0L=PN(NS)
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0L=SQRT(RTFD0L*1.D3*AEE*AMFDL)

!      IF(NR.eq.NRMAX)THEN
!         RL=RM(NR)
!         RHON=RL
!      ELSEIF(NR.EQ.NRMAX+1) THEN
!         RL=RM(NRMAX)+DELR
!         RHON=MIN(RL,1.D0)
!      ENDIF

!      CALL PL_PROF(RHON,PLF)
!      RNFDL=PLF(NS)%RN/RNFD0L
!      RTFDL=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
      RNFDL=PNS(NS)/RNFD0L*1.D-1
      RTFDL=PTS(NS)*1.D-2

      IF(MODELR.EQ.0) THEN
         FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0L)**3
         EX=PML**2/(2.D0*RTFDL/RTFD0L)
         IF(EX.GT.100.D0) THEN
            FPMXWL_S=0.D0
         ELSE
            FPMXWL_S=FACT*EXP(-EX)
         ENDIF
      ELSE
         THETA0L=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
         THETAL=THETA0L*RTFDL/RTFD0L
         Z=1.D0/THETAL
            DKBSL=BESEKN(2,Z)
            FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL) &
             *RTFD0L
            EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
            FPMXWL_S=FACT*EXP(EX)
      END IF

      RETURN
      END FUNCTION FPMXWL_S
!-------------------------------------------------------------

      SUBROUTINE FPMXWL_EDGE(NP,NSA,FL)

      implicit none
!      integer,intent(in):: NP, NSA
      integer:: NP, NSA
      integer:: NSBA, NS
      real(kind8):: FL1, FL2
!      real(kind8),intent(out):: FL
      real(kind8):: FL

      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)

!      FL1=FPMXWL_S(PM(NP,NSBA),NRMAX,NS) ! at RM(NRMAX)
!      FL2=FPMXWL_S(PM(NP,NSBA),NRMAX+1,NS) ! at RG(NRMAX+1)

!     F at R=1.0+DELR/2
!      FL=FL2*1.D-1

      FL=FPMXWL_S(PM(NP,NSBA),NRMAX,NS) 

      RETURN
      END SUBROUTINE FPMXWL_EDGE
!-------------------------------------------------------------
      FUNCTION FPMXWL_LT(PML,NR,NS)

      USE plprof
      implicit none
      integer :: NR, NS
      real(kind8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      real(kind8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      TYPE(pl_plf_type),DIMENSION(NSMAX):: plf
      real(kind8):: FPMXWL_LT

      AMFDL=PA(NS)*AMP
      AEFDL=PZ(NS)*AEE
      RNFD0L=PN(NS)
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0L=SQRT(RTFD0L*1.D3*AEE*AMFDL)

      IF(NR.eq.0)THEN
         RL=0.D0
         RHON=ABS(RL)
      ELSEIF(NR.EQ.NRMAX+1) THEN
         RL=RM(NRMAX)+DELR
         RHON=MIN(RL,1.D0)
      ELSE
         RL=RM(NR)
         RHON=RL
      ENDIF

      CALL PL_PROF(RHON,PLF)
      RNFDL=PLF(NS)%RN/RNFD0L
!      RTFDL=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
      RTFDL=PTS(NS)*1.D-1

      IF(MODELR.EQ.0) THEN
         FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0L)**3
         EX=PML**2/(2.D0*RTFDL/RTFD0L)
         IF(EX.GT.100.D0) THEN
            FPMXWL_LT=0.D0
         ELSE
            FPMXWL_LT=FACT*EXP(-EX)
         ENDIF
      ELSE
         THETA0L=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
         THETAL=THETA0L*RTFDL/RTFD0L
         Z=1.D0/THETAL
            DKBSL=BESEKN(2,Z)
            FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL) &
             *RTFD0L
            EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
            FPMXWL_LT=FACT*EXP(EX)
      END IF

      RETURN
      END FUNCTION FPMXWL_LT
!-------------------------------------------------------------

      SUBROUTINE FPWAVE_CONST
!
      USE fpmpi
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NSBA, NP, NTH, NS, NPS
      integer:: IERR
      real(kind8):: RSUM_W,RSUM_EC,RSUM_IC
      real(kind8):: PV, WPL, WPM, WPP
      real(kind8):: DFP, DFT, FFP, FACT

      CALL mtx_set_communicator(comm_np) 
      DO NR=NRSTART,NREND
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)
            NSBA=NSB_NSA(NSA)

            RSUM_W=0.D0
            RSUM_EC=0.D0
            RSUM_IC=0.D0

            IF(NPSTART.eq.1)THEN
               NPS=2
            ELSE
               NPS=NPSTART
            END IF
!            DO NP=2,NPMAX
            DO NP=NPS,NPEND
               PV=SQRT(1.D0+THETA0(NSA)*PG(NP,NSBA)**2)
               DO NTH=1,NTHMAX
                  WPL=WEIGHP(NTH  ,NP,NR,NSA)
                  IF(NTH.EQ.1) THEN
                     WPM=0.D0
                  ELSE
                     WPM=WEIGHP(NTH-1,NP,NR,NSA)
                  ENDIF
                  IF(NTH.EQ.NTHMAX) THEN
                     WPP=0.D0
                  ELSE
                     WPP=WEIGHP(NTH+1,NP,NR,NSA)
                  ENDIF
                  DFP=    PG(NP,NSBA) &
                       /DELP(NSBA)*(FNSP(NTH,NP,NR,NSBA)-FNSP(NTH,NP-1,NR,NSBA))
                  IF(NTH.EQ.1) THEN
                     DFT=1.D0/DELTH                             &
                         *(                                     &
                            ((1.D0-WPP)*FNSP(NTH+1,NP  ,NR,NSBA)   &
                                  +WPP *FNSP(NTH+1,NP-1,NR,NSBA))&
                           -                                    &
                            ((1.D0-WPM)*FNSP(NTH,NP  ,NR,NSBA)     &
                                  +WPM *FNSP(NTH,NP-1,NR,NSBA))&
                          )

                  ELSE IF(NTH.EQ.NTHMAX) THEN
                     DFT=    1.D0/DELTH                         & 
                         *(-                                    &
                            ((1.D0-WPM)*FNSP(NTH-1,NP  ,NR,NSBA)   &
                                  +WPM *FNSP(NTH-1,NP-1,NR,NSBA))&
                          +                                     &
                            ((1.D0-WPP)*FNSP(NTH,NP  ,NR,NSBA)     &
                                  +WPP *FNSP(NTH,NP-1,NR,NSBA))&
                          )
                  ELSE
                     DFT=    1.D0/(2.D0*DELTH)                  &
                         *(                                     &
                            ((1.D0-WPP)*FNSP(NTH+1,NP  ,NR,NSBA)   &
                                  +WPP *FNSP(NTH+1,NP-1,NR,NSBA))&
                           -                                    &
                            ((1.D0-WPM)*FNSP(NTH-1,NP  ,NR,NSBA)   &
                                  +WPM *FNSP(NTH-1,NP-1,NR,NSBA))&
                                  )
                  ENDIF

                  RSUM_W = RSUM_W+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                         *(DWPP(NTH,NP,NR,NSA)*DFP           &
                          +DWPT(NTH,NP,NR,NSA)*DFT)
                  RSUM_IC = RSUM_IC+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                         *(DWICPP(NTH,NP,NR,NSA)*DFP         &
                          +DWICPT(NTH,NP,NR,NSA)*DFT)
                  RSUM_EC = RSUM_EC+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                         *(DWECPP(NTH,NP,NR,NSA)*DFP         &
                          +DWECPT(NTH,NP,NR,NSA)*DFT)
               ENDDO
            ENDDO
            CALL p_theta_integration(RSUM_W)
            CALL p_theta_integration(RSUM_IC)
            CALL p_theta_integration(RSUM_EC)
               
            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)
            RPW_IMPL(NR,NSA,N_IMPL)=-RSUM_W*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6 
            RPWIC_IMPL(NR,NSA,N_IMPL)=-RSUM_IC*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RPWEC_IMPL(NR,NSA,N_IMPL)=-RSUM_EC*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            IF(N_IMPL.eq.0)THEN
!               WRITE(6,'("ALERT ", 3I4)') NR, NSA, N_IMPL
               RPW_INIT(NR,NSA)=-RSUM_W*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6 
               RPWIC_INIT(NR,NSA)=-RSUM_IC*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
               RPWEC_INIT(NR,NSA)=-RSUM_EC*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            END IF
         ENDDO
      ENDDO

      CALL mtx_reset_communicator
      RETURN
      END SUBROUTINE FPWAVE_CONST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE NBI_SOURCE_A1(NSA)

      USE fpmpi
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSA
      INTEGER:: NTH, NP, NR, NBEAM, NS, NSBA, NSABEAM, NSAX
      DOUBLE PRECISION:: PSP, ANGSP, SUML, SPL, SPFS, PSI, PCOS
      DOUBLE PRECISION:: TH0B, PANGSP

      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)

!     NBI distribute Gaussian in rho space
!     and distribute as delta function in p, theta, poloidal angle, respectively. 

      CALL mtx_set_communicator(comm_np)
      DO NBEAM = 1, NBEAMMAX
         IF(NSSPB(NBEAM).EQ.NS_NSA(NSA)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*SPBENG(NBEAM)*AEE)/PTFP0(NSA)
            PANGSP=PI*SPBPANG(NBEAM)/180.D0 ! poloidal angle where the beam deposit
            ANGSP=PI*SPBANG(NBEAM)/180.D0 ! pitch angle at phi = PANGSP
            SUML = 0.D0
            DO NR=1, NRMAX
               PSI = (1.D0+EPSRM2(NR))/(1.D0+EPSRM2(NR)*COS(PANGSP) )
               TH0B = ASIN( SQRT(SIN(ANGSP)**2/PSI) )
!               DO NP=1, NPMAX-1
               DO NP=NPSTART, NPEND
                  IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                     DO NTH=1, NTHMAX
                        IF(THG(NTH).LE.TH0B.AND.THG(NTH+1).GT.TH0B) THEN 
!                           PCOS = COS(ANGSP)
!                           SPFS = VOLP(NTH,NP,NSBA)*COSM(NTH)/(RFSADG(NR)*PCOS)
                           SPFS = VOLP(NTH,NP,NSBA)!*RLAMDAG(NTH,NR)/RFSADG(NR)*RCOEFNG(NR)
                           SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2) 
                           SUML = SUML &
                                + SPFS * SPL * VOLR(NR)
                        END IF
                     END DO
                  END IF
               END DO ! NP
            END DO ! NR

            CALL p_theta_integration(SUML)

            SUML=SUML*RNFP0(NSA)
            DO NR=NRSTART, NREND
               PSI = (1.D0+EPSRM2(NR))/(1.D0+EPSRM2(NR)*COS(PANGSP) )
               TH0B = ASIN( SQRT(SIN(ANGSP)**2/PSI) )
               SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
!               DO NP=1, NPMAX-1
               DO NP=NPSTART, NPEND
                  IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                     DO NTH=1, NTHMAX
                        IF(THG(NTH).LE.TH0B.AND.THG(NTH+1).GT.TH0B) THEN
!                           WRITE(*,'(A,2I4,2E14.6)') "TH0B", NR, NTH, TH0B, ANGSP
                           SPPB(NTH,NP,NR,NSA)=SPPB(NTH,NP,NR,NSA) &
                                + SPBTOT(NBEAM)*SPL/SUML!/RFSADG(NR)
                        END IF
                     END DO
                  END IF
               END DO
            END DO
         END IF
      END DO ! NBEAM
!  ----  FOR ELECTRON (NS=1)
!      IF(NS_NSA(NSA).EQ.1) THEN
      IF(NS_NSA(NSA).EQ.1.and.NSSPB(1).ne.1) THEN
         DO NBEAM=1,NBEAMMAX
            NSABEAM=0
            DO NSAX=1,NSAMAX
               IF(NS_NSA(NSAX).EQ.NSSPB(NBEAM)) NSABEAM=NSAX
            ENDDO
            IF(NSABEAM.NE.0) THEN
               PSP=SQRT(2.D0*AMFP(NSA)**2*SPBENG(NBEAM)*AEE &
                    /AMFP(NSABEAM))/PTFP0(NSA)
               PANGSP=PI*SPBPANG(NBEAM)/180.D0 ! poloidal angle where the beam deposit
               ANGSP=PI*SPBANG(NBEAM)/180.D0 ! pitch angle at phi = PANGSP
               SUML = 0.D0
               DO NR=1,NRMAX
                  PSI = (1.D0+EPSRM2(NR))/(1.D0+EPSRM2(NR)*COS(PANGSP) )
                  TH0B = ASIN( SQRT(SIN(ANGSP)**2/PSI) )
!                  DO NP=1,NPMAX-1
                  DO NP=NPSTART,NPEND
                     IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                        DO NTH=1,NTHMAX
                           IF(THG(NTH).LE.TH0B.AND.THG(NTH+1).GT.TH0B) THEN
!                              PCOS = COS(ANGSP)
!                              SPFS = VOLP(NTH,NP,NSBA)*COSM(NTH)/(RFSADG(NR)*PCOS)
                              SPFS = VOLP(NTH,NP,NSBA)!*RLAMDAG(NTH,NR)/RFSADG(NR)
                              SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2) 
                              SUML = SUML &
                                   + SPFS * SPL * VOLR(NR)
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
               CALL p_theta_integration(SUML)
               DO NR=NRSTART,NREND
                  PSI = (1.D0+EPSRM2(NR))/(1.D0+EPSRM2(NR)*COS(PANGSP) )
                  TH0B = ASIN( SQRT(SIN(ANGSP)**2/PSI) )
!                  DO NP=1,NPMAX-1
                  DO NP=NPSTART,NPEND
                     IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                        DO NTH=1,NTHMAX
                           IF(THG(NTH).LE.TH0B.AND.THG(NTH+1).GT.TH0B) THEN
                              SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                              SPPB(NTH,NP,NR,NSA)=SPPB(NTH,NP,NR,NSA) &
                                   +PZ(NSABEAM)*SPBTOT(NBEAM)*SPL/SUML!/RFSADG(NR)
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO

            ENDIF ! NBEAM eq 0
         END DO
      END IF
      CALL mtx_reset_communicator

      END SUBROUTINE NBI_SOURCE_A1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE NBI_SOURCE_A0(NSA)

      USE fpmpi
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSA
      integer:: NSB, NSBA, NR, NTH, NP, NS, ID
      integer:: NBEAM, NSABEAM, NSAX
      DOUBLE PRECISION:: PSP, SUML, ANGSP, SPL

      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)
      CALL mtx_set_communicator(comm_np)

      DO NBEAM=1,NBEAMMAX
         IF(NSSPB(NBEAM).EQ.NS_NSA(NSA)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*SPBENG(NBEAM)*AEE)/PTFP0(NSA)
            ANGSP=PI*SPBANG(NBEAM)/180.D0
            SUML=0.D0
!            DO NP=1,NPMAX-1
            DO NP=NPSTART,NPEND
               IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                  DO NTH=1,NTHMAX
                     IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                        DO NR=1,NRMAX
                           SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                           SUML=SUML &
                                +SPL*VOLP(NTH,NP,NSBA)*VOLR(NR)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            SUML=SUML*RNFP0(NSA)
            CALL p_theta_integration(SUML)
!            DO NP=1,NPMAX-1
            DO NP=NPSTART,NPEND
               IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                  DO NTH=1,NTHMAX
                     IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                        DO NR=NRSTART,NREND
                           SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                           SPPB(NTH,NP,NR,NSA)=SPPB(NTH,NP,NR,NSA) &
                                + SPBTOT(NBEAM)*SPL/SUML
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
!  ----  FOR ELECTRON (NS=1)
!      IF(NS_NSA(NSA).EQ.1) THEN
      IF(NS_NSA(NSA).EQ.1.and.NSSPB(1).ne.1) THEN
         DO NBEAM=1,NBEAMMAX
            NSABEAM=0
            DO NSAX=1,NSAMAX
               IF(NS_NSA(NSAX).EQ.NSSPB(NBEAM)) NSABEAM=NSAX
            ENDDO

            IF(NSABEAM.NE.0) THEN
               PSP=SQRT(2.D0*AMFP(NSA)**2*SPBENG(NBEAM)*AEE &
                    /AMFP(NSABEAM))/PTFP0(NSA)
               ANGSP=PI*SPBANG(NBEAM)/180.D0
               SUML=0.D0
!               DO NP=1,NPMAX-1
               DO NP=NPSTART,NPEND
                  IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                     DO NTH=1,NTHMAX
                        IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                           DO NR=1,NRMAX
                              SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                              SUML=SUML &
                                   +SPL*VOLP(NTH,NP,NSBA)*VOLR(NR)
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
               CALL p_theta_integration(SUML)
!               DO NP=1,NPMAX-1
               DO NP=NPSTART,NPEND
                  IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                     DO NTH=1,NTHMAX
                        IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                           DO NR=NRSTART,NREND
                              SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                              SPPB(NTH,NP,NR,NSA)=SPPB(NTH,NP,NR,NSA) &
                                   +PZ(NSABEAM)*SPBTOT(NBEAM)*SPL/SUML
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF

         END DO
      END IF

      CALL mtx_reset_communicator
      END SUBROUTINE NBI_SOURCE_A0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FUSION_SOURCE_S2A0(NSA)

      USE fpnfrr
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSA
      integer:: NSB, NSBA, NR, NTH, NP, NS, ID
      integer:: NBEAM, NSABEAM, NSAX
      DOUBLE PRECISION:: PSP, SUML, ANGSP, SPL

      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)

      DO ID=1,6
         IF(NSA.EQ.NSA1_NF(ID)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*ENG1_NF(ID)*AEE)/PTFP0(NSA)
            DO NR=NRSTART,NREND
               CALL NF_REACTION_RATE(NR,ID)
!               DO NP=1,NPMAX-1
               DO NP=NPSTART,NPEND
                  IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NSBA)
                     ENDDO
                     DO NTH=1,NTHMAX
                        SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                             +RATE_NF(NR,ID)/SUML/RNFP0(NSA)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
            IF(PSP.ge.PG(NPMAX,NSBA))THEN
               NP=NPMAX
               IF(N_IMPL.eq.0.or.N_IMPL.gt.LMAXFP)THEN
                  write(6,'(A,I5,1P3E12.4)') '  |-NP,PSP,PG=',&
                       NP,PSP,PMAX(NSBA)
                  WRITE(6,*) ' |-  OUT OF RANGE PMAX'
               END IF
               IF(NPEND.eq.NPMAX)THEN
                  DO NR=NRSTART,NREND
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NSBA)
                     ENDDO
                     DO NTH=1,NTHMAX
                        SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                             +RATE_NF(NR,ID)/SUML/RNFP0(NSA)
                     ENDDO
                  ENDDO
               END IF
            END IF
         ENDIF ! NSA=NSA1_NF(ID)
         IF(NSA.EQ.NSA2_NF(ID)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*ENG2_NF(ID)*AEE)/PTFP0(NSA)
            DO NR=NRSTART,NREND
               CALL NF_REACTION_RATE(NR,ID)
!               DO NP=1,NPMAX-1
               DO NP=NPSTART,NPEND
                  IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NSBA)
                     ENDDO
                     DO NTH=1,NTHMAX
                        SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                             +RATE_NF(NR,ID)/SUML/RNFP0(NSA)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         IF( NSA.EQ.NSA1_NF(ID).or.NSA.EQ.NSA2_NF(ID) ) THEN
            DO NR=NRSTART,NREND
!               DO NP=1,NPMAX
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     SPPF(NTH,NP,NR,NSB1_NF(ID))=                  &
                          SPPF(NTH,NP,NR,NSB1_NF(ID))              &
                          -RATE_NF_D1(NTH,NP,NR,ID)                &
                          /RNFP0(NSB1_NF(ID))
                     SPPF(NTH,NP,NR,NSB2_NF(ID))=                  &
                          SPPF(NTH,NP,NR,NSB2_NF(ID))              &
                          -RATE_NF_D2(NTH,NP,NR,ID)                &
                          /RNFP0(NSB2_NF(ID))
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO ! ID
      
      END SUBROUTINE FUSION_SOURCE_S2A0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FUSION_SOURCE_S2A1(NSA)

      USE fpnfrr
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSA
      integer:: NSB, NSBA, NR, NTH, NP, NS, ID
      integer:: NBEAM, NSABEAM, NSAX
      DOUBLE PRECISION:: PSP, SUML, ANGSP, SPL

      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)

      DO ID=1,6
         IF(NSA.EQ.NSA1_NF(ID)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*ENG1_NF(ID)*AEE)/PTFP0(NSA)
            DO NR=NRSTART,NREND
               CALL NF_REACTION_RATE(NR,ID)
!               DO NP=1,NPMAX-1
               DO NP=NPSTART,NPEND
                  IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NSBA)!*RLAMDA(NTH,NR)/RFSADG(NR)
                     ENDDO
                     DO NTH=1,NTHMAX
                        SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                             +RATE_NF(NR,ID)/SUML/RNFP0(NSA)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
            IF(PSP.ge.PG(NPMAX,NSBA))THEN
               NP=NPMAX
               IF(N_IMPL.eq.0.or.N_IMPL.gt.LMAXFP)THEN
                  write(6,'(A,I5,1P3E12.4)') '  |-NP,PSP,PG=',&
                       NP,PSP,PMAX(NSBA)
                  WRITE(6,*) ' |-  OUT OF RANGE PMAX'
               END IF
               IF(NPEND.eq.NPMAX)THEN
                  DO NR=NRSTART,NREND
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NSBA)!*RLAMDA(NTH,NR)/RFSADG(NR)
                     ENDDO
                     DO NTH=1,NTHMAX
                        SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                             +RATE_NF(NR,ID)/SUML/RNFP0(NSA)
                     ENDDO
                  ENDDO
               END IF
            END IF
         ENDIF ! NSA=NSA1_NF(ID)
         IF(NSA.EQ.NSA2_NF(ID)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*ENG2_NF(ID)*AEE)/PTFP0(NSA)
            DO NR=NRSTART,NREND
               CALL NF_REACTION_RATE(NR,ID)
!               DO NP=1,NPMAX-1
               DO NP=NPSTART,NPEND
                  IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NSBA)!*RLAMDA(NTH,NR)/RFSADG(NR)
                     ENDDO
                     DO NTH=1,NTHMAX
                        SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                             +RATE_NF(NR,ID)/SUML/RNFP0(NSA)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         IF( NSA.EQ.NSA1_NF(ID).or.NSA.EQ.NSA2_NF(ID) ) THEN
            DO NR=NRSTART,NREND
!               DO NP=1,NPMAX
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     SPPF(NTH,NP,NR,NSB1_NF(ID))=                  &
                          SPPF(NTH,NP,NR,NSB1_NF(ID))              &
                          -RATE_NF_D1(NTH,NP,NR,ID)                &
                          /RNFP0(NSB1_NF(ID))
                     SPPF(NTH,NP,NR,NSB2_NF(ID))=                  &
                          SPPF(NTH,NP,NR,NSB2_NF(ID))              &
                          -RATE_NF_D2(NTH,NP,NR,ID)                &
                          /RNFP0(NSB2_NF(ID))
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO ! ID

      END SUBROUTINE FUSION_SOURCE_S2A1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE fusion_source_init

      IMPLICIT NONE

      SPPF(:,:,:,:)=0.D0

      END SUBROUTINE fusion_source_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE fpcoef
