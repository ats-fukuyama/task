!     $Id: fpcoef.f90,v 1.38 2013/02/08 07:36:24 nuga Exp $
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
      USE fpcalr
      USE libbes,ONLY: besekn
      USE libmtx
      USE FP_READ_FIT

      INTEGER,parameter:: MODEL_T_IMP=2
!      double precision,parameter::deltaB_B=1.D-4

      contains
!-------------------------------------------------------------
      SUBROUTINE FP_COEF(NT)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NT
      integer:: NSA, NR, NTH, NP, NS
      real(kind8):: FPMAX
      integer:: NCONST_RF
      real(kind8):: DWTTEC, DWTTIC, DWTPEC, DWTPIC

      ISAVE=0

      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
      DO NR=NRSTART,NREND
         DO NP=NPSTART,NPENDWG
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=0.D0
            DPT(NTH,NP,NR,NSA)=0.D0
            FPP(NTH,NP,NR,NSA)=0.D0
            FEPP(NTH,NP,NR,NSA)=0.D0

            DWPP(NTH,NP,NR,NSA)=0.D0
            DWPT(NTH,NP,NR,NSA)=0.D0
            DLPP(NTH,NP,NR,NSA)=0.D0
            FLPP(NTH,NP,NR,NSA)=0.D0
            FSPP(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
         DO NP=NPSTARTW,NPENDWM
         DO NTH=1,NTHMAX+1
            DTP(NTH,NP,NR,NSA)=0.D0
            DTT(NTH,NP,NR,NSA)=0.D0
            FTH(NTH,NP,NR,NSA)=0.D0
            FETH(NTH,NP,NR,NSA)=0.D0

            FSTH(NTH,NP,NR,NSA)=0.D0
            DWTP(NTH,NP,NR,NSA)=0.D0
            DWTT(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
      ENDDO
      ENDDO

!     ----- Parallel electric field accleration term -----

      IF(E0.ne.0.D0) CALL FP_CALE

!     ----- Quasi-linear wave-particle interaction term -----

      IF(MODEL_WAVE.ne.0) CALL FP_CALW

!     ----- Collisional slowing down and diffusion term -----

      IF(MOD(NT,NTSTEP_COLL).EQ.0) CALL FP_CALC

!     ----- Synchrotoron radiation -----

      IF(MODEL_synch.ne.0) CALL synchrotron

!     ----- Loss term -----

      IF(MODEL_loss.ne.0) CALL loss_for_CNL

!     ----- Radial diffusion term -----

      IF(MODELD.NE.0) CALL FP_CALR

!     ----- Particle source term -----

      CALL FP_CALS

!     ----- Charge Exchange Loss Term ----
      IF(MODEL_CX_LOSS.ne.0)THEN
         DO NSA=NSASTART, NSAEND
            NS=NS_NSA(NSA)
            IF(PA(NS).eq.1.D0.or.PA(NS).eq.2.D0)THEN
               CALL CX_LOSS_TERM(NSA) 
            END IF
         END DO
      END IF
!     ----- Sum up velocity diffusion terms -----

      DO NSA=NSASTART,NSAEND
      DO NR=NRSTART,NREND
         DO NP=NPSTART,NPENDWG
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=DCPP(NTH,NP,NR,NSA)+DWPP(NTH,NP,NR,NSA) &
                 + DLPP(NTH,NP,NR,NSA)
            DPT(NTH,NP,NR,NSA)=DCPT(NTH,NP,NR,NSA)+DWPT(NTH,NP,NR,NSA)
            FPP(NTH,NP,NR,NSA)=FEPP(NTH,NP,NR,NSA)+FCPP(NTH,NP,NR,NSA) &
                 + FSPP(NTH,NP,NR,NSA) + FLPP(NTH,NP,NR,NSA)
!            WRITE(*,'(A,5I5,3E14.6)') "TEST1 ", NRANK, NSA, NR, NP, NTH, DCPP(NTH,NP,NR,NSA), FEPP(NTH,NP,NR,NSA), EP(NR)
         ENDDO
         ENDDO
!
         DO NP=NPSTARTW,NPENDWM
         DO NTH=1,NTHMAX+1
            DTP(NTH,NP,NR,NSA)=DCTP(NTH,NP,NR,NSA)+DWTP(NTH,NP,NR,NSA)
            DTT(NTH,NP,NR,NSA)=DCTT(NTH,NP,NR,NSA)+DWTT(NTH,NP,NR,NSA)
            FTH(NTH,NP,NR,NSA)=FETH(NTH,NP,NR,NSA)+FCTH(NTH,NP,NR,NSA) &
                              +FSTH(NTH,NP,NR,NSA)
!            WRITE(*,'(A,5I5,3E14.6)') "TEST2 ", NRANK, NSA, NR, NP, NTH, DCTT(NTH,NP,NR,NSA), FETH(NTH,NP,NR,NSA), EP(NR)
         ENDDO
         ENDDO
      ENDDO
      ENDDO

!     ****************************
!     Boundary condition at p=pmax
!     ****************************
!
      IF(NPENDWG.eq.NPMAX+1)THEN
         DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
         DO NTH=1,NTHMAX
            DPP(NTH,NPMAX+1,NR,NSA)=0.D0
            DPT(NTH,NPMAX+1,NR,NSA)=0.D0
            IF(MODEL_disrupt.eq.1)THEN
               IF(MODEL_RE_pmax.eq.0)THEN
                  FPP(NTH,NPMAX+1,NR,NSA)=max(0.D0,FPP(NTH,NPMAX+1,NR,NSA))
               ELSE
                  FPP(NTH,NPMAX+1,NR,NSA)=0.D0
               END IF
            ELSE
               FPP(NTH,NPMAX+1,NR,NSA)=0.D0
            END IF
         END DO
         END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE fp_coef

! ****************************************
!     Parallel electric field
! ****************************************

      SUBROUTINE FP_CALE

      IMPLICIT NONE
      integer:: NSA, NSB, NR, NTH, NP
!      real(kind8):: PSP, SUML, ANGSP, SPL, FPMAX
      integer:: NG

      IF(MODELA.eq.0)THEN
      DO NSA=NSASTART,NSAEND
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
      END DO
      ELSE
         DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            CALL FP_CALE_LAV(NR,NSA)
         ENDDO
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
      integer:: NTH, NP, NG, ITLB, ITUB
      double precision:: DELH, SUM, ETAL, X, PSIB, PSIN

!     BOUNCE AVERAGE FEPP
      DO NP=NPSTART,NPENDWG
         DO NTH=1,ITL(NR)-1
            FEPP(NTH,NP,NR,NSA) = AEFP(NSA)*EP(NR)/PTFP0(NSA)*COSM(NTH)*Line_Element(NR)*A_chi0(NR)*2.D0*PI

         END DO ! END NTH
         DO NTH=ITL(NR),ITU(NR)
            FEPP(NTH,NP,NR,NSA)=0.D0
         END DO
         DO NTH=ITU(NR)+1,NTHMAX
            FEPP(NTH,NP,NR,NSA) = AEFP(NSA)*EP(NR)/PTFP0(NSA)*COSM(NTH) &
                                *Line_Element(NR)*A_chi0(NR)*2.D0*PI
         END DO
         ITLB=ITL(NR)-1
         ITUB=NTHMAX-ITLB+1

!         FEPP(ITLB,NP,NR,NSA)=RLAMDA(ITLB,NR)/4.D0 &
!              *( FEPP(ITLB-1,NP,NR,NSA)/RLAMDA(ITLB-1,NR) &
!              +FEPP(ITLB+1,NP,NR,NSA)/RLAMDA(ITLB+1,NR) &
!              -FEPP(ITUB-1,NP,NR,NSA)/RLAMDA(ITUB-1,NR) &
!              -FEPP(ITUB+1,NP,NR,NSA)/RLAMDA(ITUB+1,NR))
!         FEPP(ITUB,NP,NR,NSA)=-FEPP(ITLB,NP,NR,NSA)
      END DO ! END NP

!     BOUNCE AVERAGE FETH
      ITLB=ITL_judge(NR)
      ITUB=NTHMAX+2-ITLB
      DO NP=NPSTARTW,NPENDWM
         DO NTH=1,ITLB ! FETH(NTH=1)=0
            FETH(NTH,NP,NR,NSA) = -AEFP(NSA)*EP(NR)/PTFP0(NSA)*SING(NTH) &
                                  *Line_Element(NR)*A_chi0(NR)*2.D0*PI
         END DO
         DO NTH=ITLB+1,ITUB-1
            FETH(NTH,NP,NR,NSA)=0.D0
         END DO
         DO NTH=ITUB,NTHMAX+1 ! FETH(NTHMAX+1)=0 
            FETH(NTH,NP,NR,NSA) = -AEFP(NSA)*EP(NR)/PTFP0(NSA)*SING(NTH) &
                                  *Line_Element(NR)*A_chi0(NR)*2.D0*PI
         END DO
      END DO ! END NP

      END SUBROUTINE FP_CALE_LAV

! ****************************************
!     Particle source and sink
! ****************************************

      SUBROUTINE FP_CALS

      USE fpnfrr
      IMPLICIT NONE
      integer:: NSA, NSB, NR, NTH, NP, NS, ID
      integer:: NBEAM, NSABEAM, NSAX, ISW_LOSS
      real(kind8):: PSP, SUML, ANGSP, SPL, FL, const_inv_tau

      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)

!     ----- Particle source term -----

         DO NR=NRSTART,NREND
         DO NP=NPSTART,NPEND
            DO NTH=1,NTHMAX
               PPL(NTH,NP,NR,NSA)=0.D0
               SPPL(NTH,NP,NR,NSA)=0.D0
               SPPB(NTH,NP,NR,NSA)=0.D0
               SPPS(NTH,NP,NR,NSA)=0.D0
               SPPI(NTH,NP,NR,NSA)=0.D0
               SPPD(NTH,NP,NSA)=0.D0
               IF(MODEL_CX_LOSS.ne.0) SPPL_CX(NTH,NP,NR,NSA)=0.D0
            ENDDO
         ENDDO
         ENDDO
!     ----- NBI source term -----

      IF(MODEL_NBI.eq.1)THEN
!         IF(MODELA.eq.0)THEN
         CALL NBI_SOURCE_A0(NSA)
!         ELSE
!            CALL NBI_SOURCE_A1(NSA)
!         END IF
      ELSEIF(MODEL_NBI.eq.2)THEN
         CALL NBI_SOURCE_FIT3D(NSA)
      END IF
!     ----- Fixed fusion source term -----

      IF(MODELS.EQ.1) THEN
         IF(NSSPF.EQ.NS_NSA(NSA)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*SPFENG*AEE)/PTFP0(NSA)
            SUML=0.D0
            DO NP=NPSTART,NPEND
               IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                  DO NR=1,NRMAX
                     SPL=EXP(-(RM(NR)-SPFR0)**2/SPFRW**2)
                     DO NTH=1,NTHMAX
                        SUML=SUML &
                            +SPL*VOLP(NTH,NP,NS)*VOLR(NR)!*RLAMDAG(NTH,NR)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
            SUML=SUML*RNFP0(NSA)
            DO NP=NPSTART,NPEND
               IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
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
      END IF ! MODELS=1

!     ----- non-Maxwell fusion source term -----
      IF(MODELS.EQ.2.OR.MODELS.EQ.3) THEN
!         IF(MODELA.eq.0)THEN
            CALL FUSION_SOURCE_S2A0(NSA)
!         ELSE
!            CALL FUSION_SOURCE_S2A1(NSA)            
!         END IF
      ENDIF ! MODELS=2 or MODELS=3
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
         IF(ISW_LOSS.eq.0)THEN ! loss term depends on the present FNSP
            DO NR=NRSTART,NREND
               DO NP=NPSTART,NPEND
                  IF(PM(NP,NS).le.3.D0)THEN ! for beam benchmark
                     FL=FPMXWL(PM(NP,NS),NR,NS) 
                     DO NTH=1,NTHMAX
                        PPL(NTH,NP,NR,NSA)=-1.D0/TLOSS(NS)*RLAMDA(NTH,NR)
!                        SPPS(NTH,NP,NR,NSA)= FL /TLOSS(NS)!*RLAMDA(NTH,NR)
                     ENDDO
                  END IF
               ENDDO
            ENDDO
         ELSEIF(ISW_LOSS.eq.1)THEN ! loss term depends on initial FNS, const total density
            SUML=0.D0
            DO NR=NRSTART,NREND
               DO NP=NPSTART,NPEND
                  IF(PM(NP,NS).le.NP_BULK(NR,NSA))THEN ! for beam benchmark
                     FL=FPMXWL(PM(NP,NS),NR,NS) 
                     DO NTH=1,NTHMAX
                        SUML = SUML &
                             + FL*VOLP(NTH,NP,NS) * VOLR(NR)*RLAMDAG(NTH,NR)*RFSADG(NR)
!                        SPPL(NTH,NP,NR,NSA)=-FL/TLOSS(NS)*RLAMDA(NTH,NR)
!                        SPPS(NTH,NP,NR,NSA)= FL /TLOSS(NS)!*RLAMDA(NTH,NR)
                     ENDDO
                  END IF
               ENDDO
            ENDDO
            CALL mtx_set_communicator(comm_np)
            CALL p_theta_integration(SUML)
            CALL mtx_reset_communicator
            SUML = SUML*RNFP0(NSA)
            DO NR=NRSTART, NREND
               DO NP=NPSTART, NPEND
                  FL=FPMXWL(PM(NP,NS),NR,NS) 
                  IF(PM(NP,NS).le.5.D0)THEN ! for beam benchmark
                     DO NTH=1, NTHMAX
                        SPPL(NTH,NP,NR,NSA)=-SPBTOT(1)*FL/SUML
                     END DO
                  END IF
               END DO
            END DO
         ELSE
            DO NR=NRSTART,NREND
               DO NP=NPSTART,NPEND
                  FL=FPMXWL_LT(PM(NP,NS),NR,NS)
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

      const_inv_tau=1.d3 ! tau=1ms
      NS=NS_NSA(NSA)
      IF(MODEL_DELTA_F(NS).eq.1.and.MODEL_BULK_CONST.eq.2)THEN ! sink in bulk region for delta f
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPEND
               IF(NP.le.NP_BULK(NR,NSA))THEN
                  DO NTH=1,NTHMAX
                     PPL(NTH,NP,NR,NSA)=PPL(NTH,NP,NR,NSA) -const_inv_tau*RLAMDA(NTH,NR)
                  END DO
               ELSEIF(NP_BULK(NR,NSA).le.NP.and.NP.le.NP_BULK(NR,NSA)+1)THEN
                  DO NTH=1,NTHMAX
                     PPL(NTH,NP,NR,NSA)=PPL(NTH,NP,NR,NSA) +const_inv_tau*(NP-NP_BULK(NR,NSA)-1)*RLAMDA(NTH,NR)
                  END DO
               END IF
            END DO
         END DO
      END IF

      IF(MODEL_IMPURITY.ne.0.and.TIMEFP.le.5.D0*tau_quench)THEN
         CALL IMPURITY_SOURCE(NSA)
      END IF

      IF(MODEL_SINK.eq.1)THEN
         CALL DELTA_B_LOSS_TERM(NSA)
      END IF
      END DO ! NSA

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
      IF(MODEL_EX_READ_Tn.eq.0)THEN
         CALL PL_PROF(RHON,PLF)
         RNFDL=PLF(NS)%RN/RNFD0L
         RTFDL=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
      ELSEIF(MODEL_EX_READ_Tn.eq.1)THEN
         RNFDL=RN_TEMP(NR,NS)/RNFD0L
         RTFDL=RT_TEMP(NR,NS)
         IF(NR.eq.0)THEN
            RNFDL=RN_TEMP(1,NS)/RNFD0L
            RTFDL=RT_TEMP(1,NS)
         ELSEIF(NR.eq.NRMAX+1)THEN
            RNFDL=RNE_EXP_EDGE/RNFD0L
            RTFDL=RTE_EXP_EDGE
         END IF
      END IF

      IF(MODELR.EQ.0) THEN
         FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0L)**3
         EX=PML**2/(2.D0*RTFDL/RTFD0L)
         FPMXWL=FACT*EXP(-EX)
      ELSE
         THETA0L=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
         THETAL=THETA0L*RTFDL/RTFD0L
         Z=1.D0/THETAL
         DKBSL=BESEKN(2,Z)
!         WRITE(*,'(A,2I5,4E14.6)') "TEST MXWL ", NR, NS, PML, Z, THETAL, DKBSL
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

      IF(MODEL_EX_READ_Tn.eq.0)THEN
         CALL PL_PROF(RHON,PLF)
         RNFDL=PNS(NS)/RNFD0L*1.D-1
         RTFDL=PTS(NS)*1.D-2
      ELSEIF(MODEL_EX_READ_Tn.eq.1)THEN
         RNFDL=RNE_EXP_EDGE/RNFD0L*1.D-1
         RTFDL=RTE_EXP_EDGE*1.D-1
      END IF


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
      integer:: NS
      real(kind8):: FL1, FL2
!      real(kind8),intent(out):: FL
      real(kind8):: FL

      NS=NS_NSA(NSA)

!      FL1=FPMXWL_S(PM(NP,NS),NRMAX,NS) ! at RM(NRMAX)
!      FL2=FPMXWL_S(PM(NP,NS),NRMAX+1,NS) ! at RG(NRMAX+1)

!     F at R=1.0+DELR/2
!      FL=FL2*1.D-1

      FL=FPMXWL_S(PM(NP,NS),NRMAX,NS) 

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE NBI_SOURCE_A1(NSA)

      USE fpmpi
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSA
      INTEGER:: NTH, NP, NR, NBEAM, NS, NSABEAM, NSAX
      DOUBLE PRECISION:: PSP, ANGSP, SUML, SPL, SPFS, PSI, PCOS
      DOUBLE PRECISION:: TH0B, PANGSP

      NS=NS_NSA(NSA)

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
               DO NP=NPSTART, NPEND
                  IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                     DO NTH=1, NTHMAX
                        IF(THG(NTH).LE.TH0B.AND.THG(NTH+1).GT.TH0B) THEN 
                           SPFS = VOLP(NTH,NP,NS)!*RLAMDAG(NTH,NR)
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
                  IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                     DO NTH=1, NTHMAX
                        IF(THG(NTH).LE.TH0B.AND.THG(NTH+1).GT.TH0B) THEN
!                           WRITE(*,'(A,2I4,2E14.6)') "TH0B", NR, NTH, TH0B, ANGSP
                           SPPB(NTH,NP,NR,NSA)=SPPB(NTH,NP,NR,NSA) &
                                + SPBTOT(NBEAM)*SPL/SUML
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
                     IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                        DO NTH=1,NTHMAX
                           IF(THG(NTH).LE.TH0B.AND.THG(NTH+1).GT.TH0B) THEN
                              SPFS = VOLP(NTH,NP,NS)!*RLAMDAG(NTH,NR)
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
                     IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                        DO NTH=1,NTHMAX
                           IF(THG(NTH).LE.TH0B.AND.THG(NTH+1).GT.TH0B) THEN
                              SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                              SPPB(NTH,NP,NR,NSA)=SPPB(NTH,NP,NR,NSA) &
                                   +PZ(NSABEAM)*SPBTOT(NBEAM)*SPL/SUML
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
      integer:: NSB, NR, NTH, NP, NS, ID
      integer:: NBEAM, NSABEAM, NSAX
      DOUBLE PRECISION:: PSP, SUML, ANGSP, SPL

      NS=NS_NSA(NSA)
      CALL mtx_set_communicator(comm_np)

      DO NBEAM=1,NBEAMMAX
         IF(NSSPB(NBEAM).EQ.NS_NSA(NSA)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*SPBENG(NBEAM)*AEE)/PTFP0(NSA)
            ANGSP=PI*SPBANG(NBEAM)/180.D0
            SUML=0.D0
            DO NP=NPSTART,NPEND
               IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                  DO NTH=1,NTHMAX
                     IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                        DO NR=1,NRMAX
                           SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                           SUML=SUML &
                                +SPL*VOLP(NTH,NP,NS)*VOLR(NR)*RLAMDAG(NTH,NR)*RFSADG(NR)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            SUML=SUML*RNFP0(NSA)
            CALL p_theta_integration(SUML)
            DO NP=NPSTART,NPEND
               IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                  DO NTH=1,NTHMAX
                     IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                        DO NR=NRSTART,NREND
                           SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                           SPPB(NTH,NP,NR,NSA)=SPPB(NTH,NP,NR,NSA) &
                                + SPBTOT(NBEAM)*SPL/SUML!*RLAMDAG(NTH,NR)!*RFSADG(NR)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
!  ----  FOR ELECTRON (NS=1)
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
                  IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                     DO NTH=1,NTHMAX
                        IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                           DO NR=1,NRMAX
                              SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                              SUML=SUML &
                                   +SPL*VOLP(NTH,NP,NS)*VOLR(NR)*RLAMDAG(NTH,NR)*RFSADG(NR)
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
               SUML=SUML*RNFP0(NSA)
               CALL p_theta_integration(SUML)
!               DO NP=1,NPMAX-1
               DO NP=NPSTART,NPEND
                  IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                     DO NTH=1,NTHMAX
                        IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                           DO NR=NRSTART,NREND
                              SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                              SPPB(NTH,NP,NR,NSA)=SPPB(NTH,NP,NR,NSA) &
                                   +PZ(NSABEAM)*SPBTOT(NBEAM)*SPL/SUML!*RLAMDAG(NTH,NR)!*RFSADG(NR)
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
      USE fpnflg
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSA
      integer:: NSB, NR, NTH, NP, NS, ID, NSA1, NSA2
      integer:: NBEAM, NSABEAM, NSAX
      DOUBLE PRECISION:: PSP, SUML, ANGSP, SPL

      NS=NS_NSA(NSA)

      RATE_NF(:,:)=0.D0
      RATE_NF_BB(:,:)=0.D0

      DO ID=1,6
         IF(NSA.EQ.NSA1_NF(ID)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*ENG1_NF(ID)*AEE)/PTFP0(NSA)
            DO NR=NRSTART,NREND
               IF(MODELS.NE.3) CALL NF_REACTION_RATE(NR,ID)
               IF(MODELS.EQ.3) CALL NF_REACTION_RATE_LG(NR,ID)
               DO NP=NPSTART,NPEND
                  IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NS)*RLAMDAG(NTH,NR)*RFSADG(NR)
                     ENDDO
                     DO NTH=1,NTHMAX
                        SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                             +RATE_NF(NR,ID)/SUML/RNFP0(NSA)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
            IF(PSP.ge.PG(NPMAX,NS))THEN
               NP=NPMAX
               IF(N_IMPL.eq.0.or.N_IMPL.gt.LMAXFP)THEN
                  write(6,'(A,I5,1P3E12.4)') '   |-NP,PSP,PG=',&
                       NP,PSP,PMAX(NS)
                  WRITE(6,*) '  |-  OUT OF RANGE PMAX'
               END IF
               IF(NPEND.eq.NPMAX)THEN
                  DO NR=NRSTART,NREND
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NS)*RLAMDAG(NTH,NR)*RFSADG(NR)
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
               IF(MODELS.NE.3) CALL NF_REACTION_RATE(NR,ID)
               IF(MODELS.EQ.3) CALL NF_REACTION_RATE_LG(NR,ID)
!               DO NP=1,NPMAX-1
               DO NP=NPSTART,NPEND
                  IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NS)*RLAMDAG(NTH,NR)*RFSADG(NR)
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
            NSA1=NSA_NSB(NSB1_NF(ID))
            NSA2=NSA_NSB(NSB2_NF(ID))
            DO NR=NRSTART,NREND
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     SPPF(NTH,NP,NR,NSA1)=                  &
                          SPPF(NTH,NP,NR,NSA1)              &
                          -RATE_NF_D1(NTH,NP,NR,ID)                &
                          /RNFP0(NSA1)
                     SPPF(NTH,NP,NR,NSA2)=                  &
                          SPPF(NTH,NP,NR,NSA2)              &
                          -RATE_NF_D2(NTH,NP,NR,ID)                &
                          /RNFP0(NSA2)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO ! ID
      
      END SUBROUTINE FUSION_SOURCE_S2A0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FUSION_SOURCE_S2A1(NSA)

      USE fpnfrr
      USE fpnflg
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSA
      integer:: NSB, NR, NTH, NP, NS, ID, NSA1, NSA2
      integer:: NBEAM, NSABEAM, NSAX
      DOUBLE PRECISION:: PSP, SUML, ANGSP, SPL

      NS=NS_NSA(NSA)

      DO ID=1,6
         IF(NSA.EQ.NSA1_NF(ID)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*ENG1_NF(ID)*AEE)/PTFP0(NSA)
            DO NR=NRSTART,NREND
               IF(MODELS.NE.3) CALL NF_REACTION_RATE(NR,ID)
               IF(MODELS.EQ.3) CALL NF_REACTION_RATE_LG(NR,ID)
!               DO NP=1,NPMAX-1
               DO NP=NPSTART,NPEND
                  IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NS)!*RLAMDA(NTH,NR)
                     ENDDO
                     DO NTH=1,NTHMAX
                        SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                             +RATE_NF(NR,ID)/SUML/RNFP0(NSA)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
            IF(PSP.ge.PG(NPMAX,NS))THEN
               NP=NPMAX
               IF(N_IMPL.eq.0.or.N_IMPL.gt.LMAXFP)THEN
                  write(6,'(A,I5,1P3E12.4)') '  |-NP,PSP,PG=',&
                       NP,PSP,PMAX(NS)
                  WRITE(6,*) ' |-  OUT OF RANGE PMAX'
               END IF
               IF(NPEND.eq.NPMAX)THEN
                  DO NR=NRSTART,NREND
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NS)!*RLAMDA(NTH,NR)
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
               IF(MODELS.NE.3) CALL NF_REACTION_RATE(NR,ID)
               IF(MODELS.EQ.3) CALL NF_REACTION_RATE_LG(NR,ID)
!               DO NP=1,NPMAX-1
               DO NP=NPSTART,NPEND
                  IF(PG(NP,NS).LE.PSP.AND.PG(NP+1,NS).GT.PSP) THEN
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NS)!*RLAMDA(NTH,NR)
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
            NSA1=NSA_NSB(NSB1_NF(ID))
            NSA2=NSA_NSB(NSB2_NF(ID))
            DO NR=NRSTART,NREND
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     SPPF(NTH,NP,NR,NSA1)=                  &
                          SPPF(NTH,NP,NR,NSA1)              &
                          -RATE_NF_D1(NTH,NP,NR,ID)                &
                          /RNFP0(NSA1)
                     SPPF(NTH,NP,NR,NSA2)=                  &
                          SPPF(NTH,NP,NR,NSA2)              &
                          -RATE_NF_D2(NTH,NP,NR,ID)                &
                          /RNFP0(NSA2)
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
!      RATE_NF(:,:)=0.D0
!      RATE_NF_BB(:,:)=0.D0

      END SUBROUTINE fusion_source_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE synchrotron

      IMPLICIT NONE
      real(8):: alpha, rgama_para, u, rgama
      integer:: NSA, NTH, NP, NR, NS

      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
      DO NR=NRSTART, NREND

         DO NP=NPSTART, NPENDWG
            rgama=SQRT(1.D0+THETA0(NS)*PG(NP,NSA)**2)
            alpha=(2.D0*AEE**4)/(3*AMFP(NSA)**3*VC**5*RGAMA)*1.e30
            u=PTFP0(NSA)*PG(NP,NSA)/AMFP(NSA)
            DO NTH=1, NTHMAX
               rgama_para=SQRT(1.D0+THETA0(NS)*PG(NP,NSA)**2*COSM(NTH)**2)
               
               FSPP(NTH,NP,NR,NSA)= -&
                    alpha*BB**2*rgama_para**2*U*SINM(NTH)**2* &
                    (1.D0+ (U/(rgama_para*VC))**2*COSM(NTH)**2 ) &
                    *AMFP(NSA)/PTFP0(NSA)
            END DO
         END DO
         
         DO NP=NPSTARTW, NPENDWM
            rgama=SQRT(1.D0+THETA0(NS)*PM(NP,NSA)**2)
            alpha=(2.D0*AEE**4)/(3*AMFP(NSA)**3*VC**5*RGAMA)*1.e30
            u=PTFP0(NSA)*PM(NP,NSA)/AMFP(NSA)
            DO NTH=1, NTHMAX+1
               rgama_para=SQRT(1.D0+THETA0(NS)*PM(NP,NSA)**2*COSG(NTH)**2)
               
               FSTH(NTH,NP,NR,NSA)= -&
                    alpha*BB**2*rgama_para**2*U*SING(NTH)*COSG(NTH)* &
                    (1.D0-U**2/(rgama_para*VC)**2*SING(NTH)**2 ) &
                    *AMFP(NSA)/PTFP0(NSA)
            END DO
         END DO
         
      END DO
      END DO
      
      END SUBROUTINE synchrotron
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE loss_for_CNL

      IMPLICIT NONE
      real(8):: rgama
      INTEGER:: NSA,NTH, NP, NR, NS

      DO NSA=NSASTART, NSAEND
         NS=NS_NSA(NSA)
      DO NR=NRSTART, NREND
      DO NP=NPSTART, NPENDWG
         rgama=SQRT(1.D0+THETA0(NS)*PG(NP,NS)**2)
         IF(NP.ne.1)THEN
            DO NTH=1,NTHMAX
               FLPP(NTH,NP,NR,NSA)=-1.D0/tau_quench!/PG(NP,NSA)**2
               DLPP(NTH,NP,NR,NSA)=-FLPP(NTH,NP,NR,NSA) * &
                    (RT_quench(NR)*rgama)/(RTFP0(NSA)*PG(NP,NS))
            END DO
         ELSE
            DO NTH=1,NTHMAX
               FLPP(NTH,NP,NR,NSA)=0.D0
               DLPP(NTH,NP,NR,NSA)=0.D0
            END DO
         END IF
      END DO
      END DO
      END DO

      END SUBROUTINE loss_for_CNL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE IMPURITY_SOURCE(NSA)
! assume NSAMAX=1, NSBMAX>=2

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSA
      INTEGER:: NS,NTH,NP,NR,NSB
      real(kind8):: tau_imp, FACTZ, imp_charge, SUM_MGI

!      tau_imp=5.D0*tau_quench
      tau_imp=tau_mgi
!      N_impu=3
!      target_zeff_mgi=3.D0
!      SPITOT=0.5D0 ! m^-3 on axis ! desired value of e dens.

      NS=NS_NSA(NSA)

!     electron source
      IF(NS_NSA(NSA).eq.1) THEN
         DO NR=NRSTART, NREND
            DO NP=NPSTART, NPEND
               DO NTH=1, NTHMAX
                  SPPI(NTH,NP,NR,NSA)= &!SPPI(NTH,NP,NR,NSA) &
                       FPMXWL_IMP(PM(NP,NS),NR,NS,n_impu)/tau_imp
               END DO
            END DO
         END DO
      END IF

      END SUBROUTINE IMPURITY_SOURCE
!-------------------------------------------------------------
      FUNCTION FPMXWL_IMP(PML,NR,NSA,N_ion)

      implicit none
      integer,intent(in):: NR,NSA,N_ion
      integer:: NS
      real(kind8),intent(in):: PML
      real(kind8):: rnfp0l,rtfp0l,ptfp0l,RTFD0L
      real(kind8):: amfpl,rnfpl,rtfpl,fact,ex,theta0l,thetal,z,dkbsl
      real(kind8):: FPMXWL_IMP, target_Z, target_ne, target_ni

      NS=NS_NSA(NSA)
!      AMFPL=PA(NSA)*AMP
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      RNFP0L=PN(NS)
      RTFP0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFP0L=SQRT(RTFD0L*1.D3*AEE*AMFPL)

      target_z=ABS(PZ(NS))
      target_ni=(SPITOT-RNFP0(1))*RNFP(NR,NSA)/RNFP0(NSA)/target_z

      RNFPL=target_ni
      
      IF(MODEL_T_IMP.eq.0)THEN
         RTFPL=RT_quench(NR)
      ELSEIF(MODEL_T_IMP.eq.1)THEN
         RTFPL=( 0.1D0*RT_quench(NR) + 0.9D0*RT_quench_f(NR) )
      ELSEIF(MODEL_T_IMP.eq.2)THEN
         RTFPL=( 0.01D0*RT_quench(NR) + 0.99D0*RT_quench_f(NR) )
      END IF

      IF(MODELR.EQ.0) THEN
         FACT=RNFPL/SQRT(2.D0*PI*RTFPL/RTFP0L)**3
         EX=PML**2/(2.D0*RTFPL/RTFP0L)
         IF(EX.GT.100.D0) THEN
            FPMXWL_IMP=0.D0
         ELSE
            FPMXWL_IMP=FACT*EXP(-EX)
         ENDIF
      ELSE
         THETA0L=THETA0(NS)
         THETAL=THETA0L*RTFPL/RTFP0L
         Z=1.D0/THETAL
            DKBSL=BESEKN(2,Z)
            FACT=RNFPL*SQRT(THETA0L)/(4.D0*PI*RTFPL*DKBSL) &
             *RTFP0L
            EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
            FPMXWL_IMP=FACT*EXP(EX)
      END IF

      RETURN
      END FUNCTION FPMXWL_IMP
!-------------------------------------------------------------
      SUBROUTINE DELTA_B_LOSS_TERM(NSA)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSA
      INTEGER:: NS,NTH,NP,NR
      double precision:: tau_dB, rgama, factp, factr, diff_c

      NS=NS_NSA(NSA)

      DO NR=NRSTART, NREND
         DO NP=NPSTART, NPEND
            DO NTH=1, NTHMAX
               RGAMA=SQRT(1.D0+THETA0(NS)*PM(NP,NS)**2)
               FACTP=PI*RR*PM(NP,NS)*ABS(COSM(NTH))/RGAMA *PTFP0(NSA)/AMFP(NSA)
               FACTR=QLM(NR)*deltaB_B**2
               diff_c = FACTR*FACTP/(RA**2)
               tau_dB = 3.D0*(1.D0-RM(NR)**2)/(4.D0*diff_c)

               SPPL(NTH,NP,NR,NSA) = -FNSP(NTH,NP,NR,NSA)/tau_dB
            END DO
         END DO
      END DO


      END SUBROUTINE DELTA_B_LOSS_TERM
!-------------------------------------------------------------
      SUBROUTINE CX_LOSS_TERM(NSA)

      IMPLICIT NONE
      integer,intent(in):: NSA
      double precision:: sigma_cx, k_energy, log_energy
      double precision:: log10_neu0, log10_neus, alpha, beta
      integer:: NP, NTH, NR, NS
      double precision,dimension(NRSTART:NREND):: RN_NEU

      NS=NS_NSA(NSA)

      log10_neu0=LOG10(RN_NEU0)
      log10_neus=LOG10(RN_NEUS)

!     neutral gas profile
      SELECT CASE(MODEL_CX_LOSS)
      CASE(1)
         alpha=2.5D0
         beta=0.8D0
         DO NR=NRSTART, NREND
            RN_NEU(NR)=10**( (log10_neu0-log10_neus)*(1.D0-RM(NR)**alpha)**beta+log10_neus )
         END DO
      END SELECT

      IF(MODEL_DELTA_F(NS).eq.0)THEN
         DO NR=NRSTART, NREND
            DO NP=NPSTART, NPEND
               k_energy = 0.5D0*(PTFP0(NSA)*PM(NP,NS))**2/(AMFP(NSA)*AEE)
               log_energy = dlog10(k_energy)
!     sigma_cx in literature is written in [cm^2]
               sigma_cx = 0.6937D-14*(1.D0-0.155D0*log_energy)**2/ &
                    (1.D0+0.1112D-14*k_energy**3.3D0)*1.D-4
!            WRITE(*,'(I5,1P5E14.6)') NP, PM(NP,NSA), k_energy, k_energy*AEE, sigma_cx, sigma_cx*VTFP0(NSA)*PM(NP,NSA)
               DO NTH=1, NTHMAX
                  SPPL_CX(NTH,NP,NR,NSA) = -RN_NEU(NR)*1.D20 &
                       *sigma_cx*VTFP0(NSA)*PM(NP,NS)*FNSP(NTH,NP,NR,NSA)
               END DO
            END DO
         END DO
      ELSEIF(MODEL_DELTA_F(NS).eq.1)THEN
         DO NR=NRSTART, NREND
            DO NP=NPSTART, NPEND
               k_energy = 0.5D0*(PTFP0(NSA)*PM(NP,NS))**2/(AMFP(NSA)*AEE)
               log_energy = dlog10(k_energy)
!     sigma_cx in literature is written in [cm^2]
               sigma_cx = 0.6937D-14*(1.D0-0.155D0*log_energy)**2/ &
                    (1.D0+0.1112D-14*k_energy**3.3D0)*1.D-4
!            WRITE(*,'(I5,1P5E14.6)') NP, PM(NP,NSA), k_energy, k_energy*AEE, sigma_cx, sigma_cx*VTFP0(NSA)*PM(NP,NSA)
               DO NTH=1, NTHMAX
                  SPPL_CX(NTH,NP,NR,NSA) = -RN_NEU(NR)*1.D20 &
                       *sigma_cx*VTFP0(NSA)*PM(NP,NS)*FNSP_DEL(NTH,NP,NR,NSA)
               END DO
            END DO
         END DO


      END IF

      
!      DO NR=NRSTART, NREND
!         DO NP=NPSTART, NPEND
!            DO NTH=1, NTHMAX
!               SPPL(NTH,NP,NR,NSA) = SPPL(NTH,NP,NR,NSA) + SPPL_CX(NTH,NP,NR,NSA)
!            END DO
!         END DO
!      END DO

      END SUBROUTINE CX_LOSS_TERM
!-------------------------------------------------------------
      SUBROUTINE update_fnsb_maxwell

      IMPLICIT NONE
      INTEGER:: NS, NR, NP, NTH
      double precision:: FL

      DO NS=1, NSMAX
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPEND
               FL=FPMXWL(PM(NP,NS),NR,NS)
               DO NTH=1,NTHMAX
                  FNSB(NTH,NP,NR,NS)=FL
               END DO
            ENDDO
         END DO
      END DO

      END SUBROUTINE update_fnsb_maxwell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END MODULE fpcoef
