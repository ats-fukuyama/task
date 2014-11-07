!     $Id$

! *****************
!     MAIN LOOP
! *****************

      MODULE fploop

      USE fpcomm
      use fpexec
      use fpdrexec
      use fpcoef
      use fpsave
      use libmpi
      use fpmpi
      use fpcaleind
      use fpdisrupt

      contains

!-----------------------------

      SUBROUTINE FP_LOOP

      USE libmtx
      USE FPMPI
      USE fpprep, only: Coulomb_log 
      IMPLICIT NONE
      real(kind8),dimension(NRSTART:NREND,NSAMAX):: RJNS
      real(kind8),dimension(NRSTART:NREND):: RJN,RJ3,E3,DELE
      real(kind8),dimension(NSAMAX)::RSUMF,RSUMF0,RSUM_SS
      real(kind8):: RSUMF_, RSUMF0_

      integer:: NT, NR, NP, NTH, NSA, NTI, NSBA, NTE
      integer:: L, IERR, I
      real(kind8):: RSUM, DELEM, RJNL, dw, RSUM1, RSUM2
      real(4):: gut, gut1, gut2, gut3, gut4, gut5, gut6, gut7
      real(4):: gut_ex, gut_coef, gut_1step, gut_conv, gut_cale
      real(kind8),DIMENSION(nsize):: RSUMA
      integer,dimension(NSAMAX)::NTCLSTEP2
      real(kind8):: DEPS_MAX, DEPS, DEPS1, DEPSE
      real(kind8),dimension(NSASTART:NSAEND):: DEPS_MAXVL, DEPSV, DEPS_SS_LOCAL
      real(kind8),dimension(NSAMAX):: DEPS_MAXV
      integer,dimension(NSASTART:NSAEND):: ILOCL
      integer,dimension(NSAMAX):: ILOC
      integer:: NSTEST
      real(kind8):: temp_send, temp_recv
      character:: fmt*40
      integer:: modela_temp, NSW, NSWI,its
      integer:: ILOC1, nsend,j, ISW_D, NDIMPL, IMPL_criterion
      real(8):: sigma, ip_all, ip_ohm, ip_run, DEPS_E, jbs, IP_bs, l_ind, IP_prim, DEPS_E2
      real(8),dimension(NRSTART:NREND):: E_SIGMA
      real(8):: IP_all_FP, FPL
!      real(8),dimension(NPSTART:NPEND):: DCPP_L1, DCPP_L2, FCPP_L1, FCPP_L2, DCTT_L1, DCTT_L2, DPP_L, FPP_L, DTT_L
!      real(8),dimension(NPMAX):: DCPP_1, DCPP_2, FCPP_1, FCPP_2, DCTT_1, DCTT_2, DPP_A, FPP_A, DTT_A

!     +++++ Time loop +++++

      DO NT=1,NTMAX
         IF(TIMEFP.le.time_quench_start)THEN
            IMPL_criterion=2
         ELSE
            IMPL_criterion=20
         END IF

         N_IMPL=0
         DEPS=1.D0
         DO NSA=NSASTART,NSAEND
            NSBA=NSB_NSA(NSA)
            DO NR=NRSTART-1,NREND+1 ! local
               IF(NR.ge.1.and.NR.le.NRMAX)THEN
                  DO NP=NPSTARTW,NPENDWM
                     DO NTH=1,NTHMAX
!                        FNSP(NTH,NP,NR,NSBA)=FNS(NTH,NP,NR,NSBA)  ! new step: variant each N_IMPL ! for fp_load
                        FNSM(NTH,NP,NR,NSBA)=FNSP(NTH,NP,NR,NSBA) ! old step: invariant during N_IMPL 
                     END DO
                  END DO
               END IF
            END DO
         END DO
 
         gut_EX = 0.D0
         gut_conv=0.D0
         gut_COEF= 0.D0
         gut_1step= 0.D0
         gut_cale=0.D0
         CALL GUTIME(gut5)

         IF(MODEL_DISRUPT.ne.0)THEN
            DO NR=NRSTART,NREND
               EM(NR)=E1(NR)
               EP(NR)=E1(NR)
               SIGMA_SPM(NR)=conduct_sp(NR)
            END DO
            DO NR=NRSTART-1, NREND+1
               IF(NR.eq.0)THEN
                  EM_W(NR)=0.D0
               ELSEIF(NR.ne.NRMAX+1)THEN
                  EM_W(NR)=E1(NR)
               END IF
            END DO
            IF(TIMEFP.ge.time_quench_start) CALL Tquench_trans
            CALL set_post_disrupt_Clog
            DO NR=NRSTART,NREND
               CALL SPITZER_SIGMA(NR,sigma)
               SIGMA_SPP(NR)=sigma
               RJ_bsm(NR)=RJ_bs(NR)
            END DO
            CALL mtx_set_communicator(comm_nr)
            call mtx_allgather_real8(SIGMA_SPP,NREND-NRSTART+1,conduct_sp)
            CALL mtx_reset_communicator
!            IF(MODEL_BS.ge.1)THEN
!               DO NR=1,NRMAX
!                  call bootstrap_current(NR,jbs)
!                  RJ_bs(NR)=jbs
!               END DO
!            END IF
            DO NSA=NSASTART, NSAEND
               DO NR=1,NRMAX
                  RJS_M(NR,NSA)=RJS(NR,NSA)
               END DO
            END DO
!            IF(MODEL_jfp.eq.1)THEN
               DO NSA=NSASTART,NSAEND
                  IF (MOD(NT,NTCLSTEP).EQ.0) CALL FP_COEF(NSA)
               END DO
!            END IF
            DO NR=1,NRMAX
               RN_runaway_M(NR)=RN_runaway(NR)
            END DO
         END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         nsw = NSAEND-NSASTART+1
         DO WHILE(N_IMPL.le.LMAXFP) ! start do while
            N_IMPL=N_IMPL+1
            DO NSA=NSASTART,NSAEND 
               NSBA=NSB_NSA(NSA)
               DO NR=NRSTART,NREND
               DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX
                  F(NTH,NP,NR)=FNSP(NTH,NP,NR,NSBA) ! used at fpweight only!
               END DO
               END DO
               END DO

               ISW_D=1
               CALL GUTIME(gut1)
               IF(ISW_D.eq.0)THEN ! separate p, r
                  modeld_temp=modeld
                  modeld=0
                  CALL fp_exec(NSA,IERR,its) ! F1 and FNS0 changed
                  modeld=modeld_temp
                  IF(MODELD.ge.1)THEN
                     CALL fp_drexec(NSA,IERR,its)
                  END IF
               ELSEIF(ISW_D.eq.1)THEN ! 
                  CALL fp_exec(NSA,IERR,its) ! F1 and FNS0 changed
                  IERR=0
               END IF
               nt_init=1

               IF(IERR.NE.0) GOTO 250
               CALL GUTIME(gut2)
               GUT_EX = GUT_EX + (gut2-gut1)

               RSUMF(NSA)=0.D0
               RSUMF0(NSA)=0.D0
               RSUM_SS(NSA)=0.D0
               DO NR=NRSTART,NREND
               DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX
                  RSUMF(NSA)=RSUMF(NSA) &
                         +ABS(FNSP(NTH,NP,NR,NSBA)-FNS0(NTH,NP,NR,NSBA) )**2
                  RSUMF0(NSA)=RSUMF0(NSA) &
                         +ABS(FNSM(NTH,NP,NR,NSBA))**2
               ENDDO
               ENDDO
               ENDDO
               DO NR=NRSTARTW,NRENDWM
                  DO NP=NPSTARTW,NPENDWM
                     DO NTH=1,NTHMAX
                        FNSP(NTH,NP,NR,NSBA)=FNS0(NTH,NP,NR,NSBA)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO ! END OF NSA
!---------- convergence criterion
            CALL mtx_set_communicator(comm_np) 
            DO NSA=NSASTART,NSAEND
               RSUMF_=RSUMF(NSA)
               RSUMF0_=RSUMF0(NSA)
               CALL p_theta_integration(RSUMF_) 
               CALL p_theta_integration(RSUMF0_) 
               RSUMF(NSA)=RSUMF_
               RSUMF0(NSA)=RSUMF0_
            END DO

            DEPS=0.D0
            DO NSA=NSASTART,NSAEND
               DEPSV(NSA) = RSUMF(NSA)/RSUMF0(NSA)
               DEPS1 = DEPSV(NSA)
               DEPS=MAX(DEPS,DEPS1)
            END DO
            DEPS_MAX=0.D0
            CALL mtx_set_communicator(comm_nsa) 
            CALL mtx_reduce1_real8(DEPS,1,DEPS_MAX,ILOC1) ! MAX DEPS among NSA
            DEPS = DEPS_MAX
            CALL mtx_set_communicator(comm_nr) 
            CALL mtx_reduce1_real8(DEPS,1,DEPS_MAX,ILOC1) ! MAX DEPS among NR
            CALL mtx_reset_communicator
            DEPS = DEPS_MAX

!            IF(NRANK.eq.0.and.DEPS.le.EPSFP.and.N_IMPL.ge.2)THEN
!            IF(NRANK.eq.0.and.DEPS.le.EPSFP.and.N_IMPL.ge.IMPL_criterion)THEN
!               N_IMPL=1+LMAXFP ! exit dowhile
!            ENDIF
!            CALL mtx_broadcast1_integer(N_IMPL)

            CALL mtx_set_communicator(comm_nr) !3D
            CALL mtx_allreduce_real8(DEPSV,NSW,4,DEPS_MAXVL,ILOCL) ! the peak DEPSV for each NSA

            CALL mtx_set_communicator(comm_nsa) !3D
            CALL mtx_gather_real8(DEPS_MAXVL,nsw,DEPS_MAXV) 
            CALL mtx_gather_integer(ILOCL,nsw,ILOC) 
            CALL mtx_reset_communicator

            IF (MOD(NT,NTG1STEP).EQ.0) THEN
            IF(nrank.eq.0) THEN
               WRITE(fmt,'(a16,I1,a6,I1,a3)') &
                    '(A,1PE12.4,I4,1P',NSAMAX,'E12.4,',NSAMAX,'I4)'
               WRITE(6,fmt) 'DEPS',&
                    DEPS,N_IMPL,(DEPS_MAXV(NSA),NSA=1,NSAMAX) &
                    ,(ILOC(NSA),NSA=1,NSAMAX)
            ENDIF
            END IF
!---------- end of convergence 
            CALL GUTIME(gut3)
            gut_conv = gut_conv + (gut3-gut2)

            CALL fusion_source_init

!           update FNSB 
            IF(MODELC.ge.2)THEN
               CALL mtx_set_communicator(comm_nsa)
               CALL update_fnsb
               CALL mtx_reset_communicator
            END IF
!           end of update FNSB
!
            NDIMPL=0
            DEPS_E=0.D0
            IF(MODEL_DISRUPT.ne.0.and.MODEL_jfp.eq.0)THEN
               call update_rns_rjs(IP_all_FP)
               DO while(NDIMPL.le.2)
                  NDIMPL=NDIMPL+1
                  call calculation_runaway_rate
                  CALL AVALANCHE
                  IF(TIMEFP.ge.time_quench_start) CALL E_IND_IMPLICIT
                  DEPS_E=(E1(NRSTART)-EP(NRSTART))**2/EM(NRSTART)**2

                  DO NR=NRSTART,NREND
                     EP(NR)=E1(NR)
                  END DO
                  DO NSA=NSASTART, NSAEND
                     CALL FP_CALE(NSA)
                  END DO
                  CALL update_fpp
               END DO
               IF (MOD(NT,NTG1STEP).EQ.0.and.NRANK.eq.0) &
                    WRITE(6,'(A,E14.6)') "CALE_CONVERSION = ", DEPS_E
            END IF
            DEPS_E=0.D0
            DEPS_E2=0.D0
!
            IF(MODEL_DISRUPT.ne.0.and.MODEL_jfp.eq.1)THEN
               call update_rns_rjs(IP_all_FP)
               DO NR=NRSTART, NREND
                  E_SIGMA(NR)=EP(NR)
               END DO
               DO while(NDIMPL.le.LMAXFP)
                  NDIMPL=NDIMPL+1
                  call calculation_runaway_rate
                  CALL AVALANCHE
                  IF(TIMEFP.ge.time_quench_start) CALL E_IND_IMPLICIT_FP(E_SIGMA)
!                  IF(TIMEFP.ge.time_quench_start) CALL E_IND_CRANK_NICOLSON_FP
!                  IF(TIMEFP.ge.time_quench_start) CALL E_IND_PSI_IMPLICIT_FP
                  DEPS_E=(E1(NRSTART)-EP(NRSTART))**2/EM(NRSTART)**2

                  CALL mtx_set_communicator(comm_nr) 
                  CALL mtx_reduce1_real8(DEPS_E,1,DEPS_MAX,ILOC1) ! MAX DEPS among NR
                  CALL mtx_reset_communicator
                  DEPS_E = DEPS_MAX
                  IF(NRANK.eq.0.and.DEPS_E.le.EPSFP.and.NDIMPL.ge.2)THEN
                     NDIMPL=NDIMPL+LMAXFP ! exit dowhile
                  ENDIF
                  CALL mtx_broadcast1_integer(NDIMPL)
                  DO NR=NRSTART,NREND
                     EP(NR)=E1(NR)
                  END DO
                  DO NSA=NSASTART, NSAEND
                     CALL FP_CALE(NSA)
                  END DO
                  CALL update_fpp
               END DO

               DEPS_E2=(EP(NRSTART)-E_SIGMA(NRSTART))**2/E_SIGMA(NRSTART)**2
               CALL mtx_set_communicator(comm_nr) 
               CALL mtx_reduce1_real8(DEPS_E2,1,DEPS_MAX,ILOC1) ! MAX DEPS among NR
               CALL mtx_reset_communicator
               DEPS_E2 = DEPS_MAX

               IF (MOD(NT,NTG1STEP).EQ.0.and.NRANK.eq.0) &
                    WRITE(6,'(A,2E14.6,I3,A,2E14.6,A,2E14.6)') "CALE_CONVERSION = ", &
                    DEPS_E, DEPS_E2, NDIMPL," RJS ", RJS(1,1),RJS(NRMAX,1), &
                    " E1 ",E1(1), E1(NRMAX)
            END IF
!                  CALL djdt

!            IF(NRANK.eq.0.and.DEPS.le.EPSFP.and.DEPS_E2.le.EPSFP.and.N_IMPL.ge.2)THEN
            IF(NRANK.eq.0.and.DEPS.le.EPSFP.and.DEPS_E2.le.EPSFP)THEN
               N_IMPL=1+LMAXFP ! exit dowhile
            ENDIF
            CALL mtx_broadcast1_integer(N_IMPL)
            
            CALL GUTIME(gut7)
            gut_cale = gut_cale + gut7-gut3
!            
            IF(MODEL_LNL.eq.0) CALL Coulomb_log
!            IF(MODEL_jfp.ne.1)THEN
!            DO NSA=NSASTART,NSAEND
!               IF (MOD(NT,NTCLSTEP).EQ.0) CALL FP_COEF(NSA)
!            END DO
!            END IF

!           sum up SPPF
!            CALL mtx_set_communicator(comm_nr) !2D
!            CALL mtx_set_communicator(comm_nrnp) !3D
!            CALL source_allreduce(SPPF)
!            CALL mtx_reset_communicator
!           end of sum up SPPF

            CALL GUTIME(gut4)
            GUT_COEF = GUT_COEF + (gut4-gut7)
         END DO ! END OF DOWHILE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO NSA=NSASTART,NSAEND
            DO NR=NRSTART,NREND
               DO NP=NPSTARTW,NPENDWM
                  DO NTH=1,NTHMAX
                     F(NTH,NP,NR)=FNSP(NTH,NP,NR,NSBA) ! used at fpweight only!
                  END DO
               END DO
            END DO
            CALL FPWEIGHT(NSA,IERR) ! SET FPWEIGHT FOR FPSAVE
         END DO

         CALL GUTIME(gut6)
         GUT_1step = gut6-gut5

         IF (MOD(NT,NTG1STEP).EQ.0) THEN
         IF(NRANK.eq.0) &
              WRITE(6,'(A,E12.4, A,E12.4, A,E12.4, A,E12.4, A,E12.4, A,E12.4)') &
              " GUT_EXEC= ", GUT_EX,   " GUT_CONV= ",GUT_conv, &
              " GUT_CALE= ", gut_cale, &
              " GUT_COEF= ", GUT_COEF, " GUT_1step= ", GUT_1step, &
              " EXEC_RATIO = ", GUT_EX/GUT_1step
         END IF
!
  250    CONTINUE
!
!         CALL FPGRAC('F -2',F,4)
!         CALL FPGRAC('F1-2',F1,4)

!     +++++ calculate and save global data +++++
!         CALL mtx_set_communicator(comm_nr)
!         call mtx_allgather_real8(EP,NREND-NRSTART+1,E1)
!         CALL mtx_reset_communicator

         CALL GUTIME(gut1)
         TIMEFP=TIMEFP+DELT

         ISAVE=0
         CALL FPSSUB
         IF (MOD(NT,NTG1STEP).EQ.0) THEN
            IF(NRANK.EQ.0) THEN
               CALL FPSGLB
               CALL FPWRTGLB
            ENDIF
         ENDIF
         IF (MOD(NT,NTG2STEP).EQ.0) THEN
            IF(NRANK.EQ.0) THEN
               CALL FPSPRF
               CALL FPWRTPRF
            ENDIF
         ENDIF
         CALL mtx_broadcast_real8(RT_T,NRMAX*NSAMAX)
         CALL mtx_broadcast_real8(RNS,NRMAX*NSAMAX)
         CALL mtx_broadcast1_integer(NTG1)
         CALL mtx_broadcast1_integer(NTG2)
!         IF(NRANK.EQ.0.AND.NTG1.GT.0) call FPWRTSNAP
         CALL GUTIME(gut2)
         IF (MOD(NT,NTG1STEP).EQ.0) THEN
            IF(NRANK.eq.0) WRITE(*,'(A,E14.6)') "--------SAVE_TIME=",gut2-gut1
         END IF

         IF(MODEL_DISRUPT.ne.0)THEN
            CALL update_disruption_quantities(IP_all, IP_ohm, IP_run, IP_prim, IP_bs, l_ind)
            IF (MOD(NT,NTG1STEP).EQ.0) THEN
            IF(NRANK.eq.0)THEN
               WRITE(6,'(7A14)') "IP_all","IP_ohm", "IP_run", "IP_primary", "IP_bs", "DEL_RJS", "sigma*E"
               WRITE(6,'(1P7E14.6)') IP_all, IP_ohm, IP_run, IP_prim, IP_bs, (RJS(1,1)-RJS_M(1,1))*1.D6, SIGMA_SPP(1)*EP(1)-SIGMA_SPM(1)*EM(1)
               WRITE(6,'(A,1PE14.6)') "Zeff= ", ZEFF
            END IF
            END IF
         END IF

         IF(IERR.NE.0) RETURN
         
         IF(NRANK.eq.0.and.model_disrupt.ne.0)THEN ! time_evol.dat
               WRITE(10,'(I6,1P30E16.7e3)') NT, TIMEFP, RT_T(1,1), &
                    RT_quench(1), RJ_ohm(1), &
                    RJ_runaway(1), RN_disrupt(1), RN_runaway(1), &
                    RN_drei(1), &
                    IP_all, ip_ohm, ip_run, ip_bs, &
                    l_ind, ZEFF, IP_all_FP
         END IF

         IF(NRANK.eq.0.and.model_disrupt.ne.0)THEN
            WRITE(11,'(1P128E14.6)') TIMEFP, (E1(NR), NR=1,NRMAX), (ER_drei(NR), NR=1,NRMAX)
!            WRITE(12,'(1P128E15.6e3)') TIMEFP, (Rconner(NR)*RN_disrupt(NR)*1.d20, NR=1,NRMAX), &
            IF(MODEL_conner_fp.eq.1)THEN
               WRITE(12,'(1P128E15.6e3)') TIMEFP, (RFP(NR)*RN_disrupt(NR)*1.d20, NR=1,NRMAX), &
                    (RFP_ava(NR)*RN_disrupt(NR)*1.d20, NR=1,NRMAX)
            ELSE
               WRITE(12,'(1P128E15.6e3)') TIMEFP, (Rconner(NR)*RN_disrupt(NR)*1.d20, NR=1,NRMAX), &
                    (RFP_ava(NR)*RN_disrupt(NR)*1.d20, NR=1,NRMAX)
            END IF

            WRITE(15,'(1P128E15.6e3)') TIMEFP, (RE_PITCH(NTH), NTH=1,NTHMAX)

            IF (MOD(NT,NTG1STEP).EQ.0) THEN ! radial.dat
               WRITE(13,'(A, 1PE15.6e3, i7)') "# TIME ", TIMEFP, N_f1
               WRITE(13,'(A)') "# RM, RN, RT, RN_P, RN_Sm J_ohm, J_run, J_bs, E1, E_drei, dndt, dndt_ava,RJS " 
               IF(MODEL_conner_fp.eq.1)THEN
                  DO NR=1,NRMAX
                     WRITE(13,'(1P13E17.8e3)'), RM(NR), RN_disrupt(NR), RT_quench(NR), & 
                       RN_runaway(NR), RN_runaway(NR)-RN_drei(NR), &
                       RJ_ohm(NR), RJ_runaway(NR), RJ_bs(NR), &
                       E1(NR), ER_drei(NR), &
                       RFP(NR)*RN_disrupt(NR)*1.D20, RFP_AVA(NR)*RN_disrupt(NR)*1.D20, &
                       RJS(NR,1)
                  END DO
               ELSE
                  DO NR=1,NRMAX
                     WRITE(13,'(1P13E17.8e3)'), RM(NR), RN_disrupt(NR), RT_quench(NR), & 
                       RN_runaway(NR), RN_runaway(NR)-RN_drei(NR), &
                       RJ_ohm(NR), RJ_runaway(NR), RJ_bs(NR), &
                       E1(NR), ER_drei(NR), &
                       Rconner(NR)*RN_disrupt(NR)*1.D20, RFP_AVA(NR)*RN_disrupt(NR)*1.D20, &
                       RJS(NR,1)
                  END DO
               END IF
               WRITE(13,'(A)') " "
               WRITE(13,'(A)') " "

               WRITE(14,'(A, 1PE15.6e3, i7)') "# TIME ", TIMEFP, N_f1 ! RE_pitch.dat
               DO NTH=1,NTHMAX
                  WRITE(14,'(I4, 1P13E17.8e3)'), NTH, RE_PITCH(NTH)
               END DO
               WRITE(14,'(A)') " "
               WRITE(14,'(A)') " "
               N_f1=N_f1+1
            END IF
         END IF
      ENDDO ! END OF NT LOOP
!      CLOSE(10)
!      CLOSE(9)

!     +++++ end of time loop +++++

!         DO NP=NPSTART,NPEND
!            DCPP_L1(NP)=DCPP2(1,NP,1,1,1)
!            DCPP_L2(NP)=DCPP2(1,NP,1,2,1)
!            FCPP_L1(NP)=FCPP2(1,NP,1,1,1)
!            FCPP_L2(NP)=FCPP2(1,NP,1,2,1)
!            DCTT_L1(NP)=DCTT2(1,NP,1,1,1)
!            DCTT_L2(NP)=DCTT2(1,NP,1,2,1)            
!            DPP_L(NP)=DPP(1,NP,1,1)
!            DTT_L(NP)=DTT(1,NP,1,1)
!            FPP_L(NP)=FPP(1,NP,1,1)
!         END DO
!         CALL mtx_set_communicator(comm_np)
!         call mtx_allgather_real8(DCPP_L1,NPEND-NPSTART+1,DCPP_1)
!         call mtx_allgather_real8(DCPP_L2,NPEND-NPSTART+1,DCPP_2)
!         call mtx_allgather_real8(FCPP_L1,NPEND-NPSTART+1,FCPP_1)
!         call mtx_allgather_real8(FCPP_L2,NPEND-NPSTART+1,FCPP_2)
!         call mtx_allgather_real8(DCTT_L1,NPEND-NPSTART+1,DCTT_1)
!         call mtx_allgather_real8(DCTT_L2,NPEND-NPSTART+1,DCTT_2)
!         call mtx_allgather_real8(DTT_L,NPEND-NPSTART+1,DTT_A)
!         call mtx_allgather_real8(DPP_L,NPEND-NPSTART+1,DPP_A)
!         call mtx_allgather_real8(FPP_L,NPEND-NPSTART+1,FPP_A)
!         CALL mtx_reset_communicator

!
      IF(NRANK.eq.0)THEN
!         open(10,file='time_evol.dat')
!         DO NTI=1,NTG1
!            WRITE(10,'(I4,F12.3,1P8E16.8)') NTI, PTG(NTI)*1000, &

!               WRITE(9,'(I4,F12.3,14E16.8)') NTI, PTG(NTI)*1000 &
!                    ,PTT2(1,NTI),PTT2(2,NTI),PIT(1,NTI),PIT(2,NTI) &
!                    ,EPTR(1,NTI),EPTR(NRMAX,NTI) &
!                    ,RTT(1,1,NTI),RTT(1,2,NTI),RTT(NRMAX,1,NTI),RTT(NRMAX,2,NTI) &
!                    ,RJT(1,1,NTI),RJT(1,2,NTI),RJT(NRMAX,1,NTI),RJT(NRMAX,2,NTI)
!               WRITE(9,645) NTI, PTG(NTI)*1000 &
!                    ,PPCT2(1,1,NTI),PPCT2(2,1,NTI),PPCT2(3,1,NTI),PPCT2(4,1,NTI),PPCT(1,NTI) &
!                    ,PPCT2(1,2,NTI),PPCT2(2,2,NTI),PPCT2(3,2,NTI),PPCT2(4,2,NTI),PPCT(2,NTI) &
!                    ,PPCT2(1,3,NTI),PPCT2(2,3,NTI),PPCT2(3,3,NTI),PPCT2(4,3,NTI),PPCT(3,NTI) &
!                    ,PPCT2(1,4,NTI),PPCT2(2,4,NTI),PPCT2(3,4,NTI),PPCT2(4,4,NTI),PPCT(4,NTI) &
!                    ,PPWT(1,NTI),PPWT(2,NTI),PPWT(3,NTI),PPWT(4,NTI) &
!                    ,PDR(1,NTI),PDR(2,NTI),PDR(3,NTI),PDR(4,NTI) &
!                    ,PWT(1,NTI),PWT(2,NTI),PWT(3,NTI),PWT(4,NTI) &
!                    ,PNT(1,NTI),PNT(2,NTI),PNT(3,NTI),PNT(4,NTI) &
!                    ,PTT2(1,NTI),PTT2(2,NTI),PTT2(3,NTI),PTT2(4,NTI) &
!                    ,PTT_BULK(1,NTI),PTT_BULK(2,NTI),PTT_BULK(3,NTI),PTT_BULK(4,NTI) &
!                    ,PSPBT(2,NTI),PSPFT(2,NTI),PSPFT(3,NTI),PSPFT(4,NTI) &
!                   ,PECT(1,NTI)
!         END DO
!         close(10)
      END IF



      CALL GUTIME(gut1)
      CALL update_fns
      CALL GUTIME(gut2)
      IF(NRANK.eq.0) WRITE(6,'(A,E14.6)') "---------TIME UPDATE FNS =",gut2-gut1

      IF(NRANK.eq.0)THEN
         DO NP=1,NPMAX
            WRITE(9,'(1PE12.4,I6,1P20E17.8e3)') PTG(NTG1)*1000, NP, PM(NP,1), &
                 PM(NP,1)*PTFP0(1)/AMFP(1)/VC/SQRT(1.D0+PM(NP,1)**2*THETA0(1)), &
                 PM(NP,1)**2, &
                 PTFP0(1)**2*PM(NP,1)**2/(AEE*AMFP(1)*1.D3), FNS(1,NP,1,1), FNS(NTHMAX,NP,1,1), &
                 FNS(NTHMAX/2,NP,1,1)
!                 , DCPP_1(NP), DCPP_2(NP), FCPP_1(NP), FCPP_2(NP), DCTT_1(NP), DCTT_2(NP), &
!                 DPP_A(NP), DTT_A(NP), FPP_A(NP), (1.D0-SQRT(1.D0+PM(NP,1)**2*THETA0(1)))/THETA0(1)
         END DO
         WRITE(9,*) " "
         WRITE(9,*) " "
!         WRITE(16,'(A, 1PE15.6e3)') "# TIME ", TIMEFP
!         DO NP=1, NPMAX
!            DO NTH=1,NTHMAX
!               FPL=LOG10(ABS(FNS(NTH,NP,1,1)))
!               IF(FPL.le.-70.D0) FPL=-70.D0
!               WRITE(16,'(4E15.6e3,2E16.7e3)') PM(NP,1), THM(NTH), PM(NP,1)*COSM(NTH), PM(NP,1)*SINM(NTH), FNS(NTH, NP, 1, 1), FPL
!            END DO
!         END DO
!         WRITE(16,*) " "
!         WRITE(16,*) " "
      END IF

!      IF(NRANK.eq.0)THEN
!         DO NSA=1,NSAMAX
!            SELECT CASE(NSA)
!            CASE(1)
!               open(9,file='fns_e.dat')
!            CASE(2)
!               open(9,file='fns_D.dat')
!            CASE(3)
!               open(9,file='fns_T.dat')
!            END SELECT
!            IF(NSA.LE.3) THEN
!               DO NR=1,NRMAX
!                  DO NP=1,NPMAX
!                     DO NTH=1,NTHMAX
!                        WRITE(9,'(3I4,3E17.8)') NR, NP, NTH, &
!                             FNS(NTH,NP,NR,NSA), &
!                             PM(NP,NSA)*COSM(NTH), PM(NP,NSA)*SINM(NTH)
!                     END DO
!                  END DO
!                  WRITE(9,*) " "
!                  WRITE(9,*) " "
!               END DO
!               close(9)
!            ENDIF
!         END DO
!      END IF

      RETURN
      END SUBROUTINE FP_LOOP

! ************************************
!     PREDICTION OF ELECTRIC FIELD
! ************************************

!      SUBROUTINE FPNEWE
!
!      IMPLICIT NONE
!      integer:: NR
!
!      DO NR=2,NRMAX
!         BP(NR)=BP(NR)+(E1(NR)-E1(NR-1))*DELT/(RA*DELR)
!      ENDDO
!
!      DO NR=NRSTART,NREND
!         RJ2(NR)=(RG(NR+1)*BP(NR+1)-RG(NR)*BP(NR)) &
!                 /(RMU0*RM(NR)*DELR*RA)
!      ENDDO
!
!      DO NR=NRSTART,NREND
!         E2(NR)=RJ2(NR)*E1(NR)/RJ1(NR)
!      ENDDO
!
!      RETURN
!      END SUBROUTINE FPNEWE

!------------------------------------
      END MODULE fploop
