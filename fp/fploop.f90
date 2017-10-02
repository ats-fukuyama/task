!     $Id: fploop.f90,v 1.40 2013/02/08 07:36:24 nuga Exp $

! *****************
!     MAIN LOOP
! *****************

      MODULE fploop

      USE fpcomm
      use fpexec
!      use fpdrexec
      use fpcoef
      use fpsave
      use libmpi
      use fpmpi
      use fpdisrupt
      use eg_read
      use fpoutdata

      contains

!-----------------------------

      SUBROUTINE FP_LOOP

      USE libmtx
      USE plprof
      USE FPMPI
      USE fpprep, only: Coulomb_log 
      USE fpnfrr
      IMPLICIT NONE
      real(kind8),dimension(NSAMAX)::RSUMF,RSUMF0,RSUM_SS
      real(kind8):: RSUMF_, RSUMF0_!, RGAMA, FACTP, FACTR, diff_c, tau_dB

      integer:: NT, NR, NP, NTH, NSA, NS
      integer:: IERR, I
      real(4):: gut_exe1, gut_exe2, gut_conv3, gut_coef1
      real(4):: gut_loop1, gut_loop2, gut_cale7, gut_coef2, gut1, gut2
      real(4):: gut_out1, gut_out2, gut_out
      real(4):: gut_ex, gut_coef, gut_1step, gut_conv, gut_cale
      real(kind8):: DEPS_MAX, DEPS, DEPS1
      real(kind8),dimension(NSASTART:NSAEND):: DEPS_MAXVL, DEPSV
      real(kind8),dimension(NSAMAX):: DEPS_MAXV
      integer,dimension(NSASTART:NSAEND):: ILOCL
      integer,dimension(NSAMAX):: ILOC
      character:: fmt*40
      integer:: NSW, its
      integer:: ILOC1, ISW_D
      real(8):: DEPS_E2
      real(8):: IP_all_FP, pitch_angle_av, FL
!      real(8),dimension(NPSTART:NPEND):: DCPP_L1, DCPP_L2, FCPP_L1, FCPP_L2, DCTT_L1, DCTT_L2, DPP_L, FPP_L, DTT_L
!      real(8),dimension(NPMAX):: DCPP_1, DCPP_2, FCPP_1, FCPP_2, DCTT_1, DCTT_2, DPP_A, FPP_A, DTT_A
      real(8),dimension(NRMAX,NSMAX):: tempt, tempn
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      real(8):: RHON
!     +++++ Time loop +++++

      DO NT=1,NTMAX
         N_IMPL=0
         DEPS=1.D0
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)
!            DO NR=NRSTART-1,NREND+1 ! local
            DO NR=NRSTARTW,NRENDWM ! local
!               IF(NR.ge.1.and.NR.le.NRMAX)THEN
                  DO NP=NPSTARTW,NPENDWM
                     DO NTH=1,NTHMAX
!                        FNSP(NTH,NP,NR,NSA)=FNS(NTH,NP,NR,NSA)  ! new step: variant in each N_IMPL ! for fp_load
                        FNSM(NTH,NP,NR,NSA)=FNSP(NTH,NP,NR,NSA) ! minus step: invariant during N_IMPL 
                        IF(MODEL_DELTA_F(NS).eq.1)THEN
                           FNSM(NTH,NP,NR,NSA)=FNSP_DEL(NTH,NP,NR,NSA) ! minus step: invariant during N_IMPL 
                        END IF
                     END DO
                  END DO
!               END IF
            END DO
         END DO
 
         gut_EX = 0.0
         gut_conv=0.0
         gut_COEF= 0.0
         gut_1step= 0.0
         gut_cale=0.0
         CALL GUTIME(gut_loop1)

         IF(MODEL_DISRUPT.ne.0)THEN
            CALL TOP_OF_TIME_LOOP_DISRUPT(NT) ! include fpcoef
         END IF
         CALL GUTIME(gut_coef2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         nsw = NSAEND-NSASTART+1
         DO WHILE(N_IMPL.le.LMAXFP) ! start do while
            N_IMPL=N_IMPL+1
            DO NSA=NSASTART,NSAEND 
               NS=NS_NSA(NSA)
               DO NR=NRSTART,NREND
               DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX
                  F(NTH,NP,NR)=FNSP(NTH,NP,NR,NSA) ! used at fpweight only!
               END DO
               END DO
               END DO

               IF(MODEL_DELTA_F(NS).eq.1)THEN
                  DO NR=NRSTART,NREND
                     DO NP=NPSTARTW,NPENDWM
                        DO NTH=1,NTHMAX
                           FNSP(NTH,NP,NR,NSA)=FNSP_DEL(NTH,NP,NR,NSA)
                        END DO
                     END DO
                  END DO
               END IF

               ISW_D=1
               CALL GUTIME(gut_exe1)
               IF(ISW_D.eq.0)THEN 
                  CALL fp_exec(NSA,IERR,its) ! F1 and FNS0 changed
               ELSEIF(ISW_D.eq.1)THEN !
                  IF(MODEL_connor_fp.eq.1.or.MODEL_DISRUPT.eq.0)THEN ! Connor model doesn't require f evolution
                     CALL fp_exec(NSA,IERR,its) ! F1 and FNS0 are changed
                  END IF
                  IERR=0
               END IF
               nt_init=1

               IF(MODEL_DELTA_F(NS).eq.1)THEN
                  DO NR=NRSTART,NREND
                     DO NP=NPSTARTW,NPENDWM
                        DO NTH=1,NTHMAX
                           FNSP(NTH,NP,NR,NSA)=FNSP_DEL(NTH,NP,NR,NSA)+FNSP_MXWL(NTH,NP,NR,NSA)
                           FNSP_DEL(NTH,NP,NR,NSA)=FNS0(NTH,NP,NR,NSA)
                           FNS0(NTH,NP,NR,NSA)=FNSP_DEL(NTH,NP,NR,NSA)+FNSP_MXWL(NTH,NP,NR,NSA)
                        END DO
                     END DO
                  END DO
               END IF

!--------------convergence check of f
               IF(IERR.NE.0) GOTO 250
               CALL GUTIME(gut_exe2)
               GUT_EX = GUT_EX + (gut_exe2-gut_exe1)

               RSUMF(NSA)=0.D0
               RSUMF0(NSA)=0.D0
               RSUM_SS(NSA)=0.D0
               DO NR=NRSTART,NREND
               DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX
                  RSUMF(NSA)=RSUMF(NSA) &
                         +ABS(FNSP(NTH,NP,NR,NSA)-FNS0(NTH,NP,NR,NSA) )**2
                  RSUMF0(NSA)=RSUMF0(NSA) &
                         +ABS(FNSP(NTH,NP,NR,NSA))**2
               ENDDO
               ENDDO
               ENDDO
               DO NR=NRSTARTW,NRENDWM
                  DO NP=NPSTARTW,NPENDWM
                     DO NTH=1,NTHMAX
                        FNSP(NTH,NP,NR,NSA)=FNS0(NTH,NP,NR,NSA)
                     ENDDO
                  ENDDO
               ENDDO
               IF(MODELD_boundary.eq.1.and.NREND.eq.NRMAX)THEN
                  CALL update_radial_f_boundary(NSA)
               END IF
            ENDDO ! END OF NSA
!!---------- convergence criterion
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

            CALL mtx_set_communicator(comm_nr) !3D
            CALL mtx_allreduce_real8(DEPSV,NSW,4,DEPS_MAXVL,ILOCL) 
                                                ! the peak DEPSV for each NSA

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
!---------- end of convergence check
            CALL GUTIME(gut_conv3)
            gut_conv = gut_conv + (gut_conv3-gut_exe2)


            DEPS_E2=0.D0 ! changed
            IF(MODEL_DISRUPT.eq.1)THEN ! E field evolution for DISRUPT
               CALL E_FIELD_EVOLUTION_DISRUPT(NT,IP_all_FP,DEPS_E2)
            END IF
!                  CALL djdt

            IF(NRANK.eq.0.and.DEPS.le.EPSFP.and.DEPS_E2.le.EPSFP)THEN
               N_IMPL=1+LMAXFP ! exit dowhile
            ENDIF
            CALL mtx_broadcast1_integer(N_IMPL)
            
            CALL GUTIME(gut_cale7)
            gut_cale = gut_cale + gut_cale7-gut_conv3
!            
            CALL update_RN_RT ! update RN_TEMP, RT_TEMP for Coulomb log and fpcalcnr
            IF(MODEL_LNL.eq.0) CALL Coulomb_log ! update coulomb log

            CALL fusion_source_init
!           update FNSB (fnsb is required by NL collsion and NF reaction)
            IF(MODELC.ge.4.or.MODELS.ge.2)THEN
               CALL mtx_set_communicator(comm_nsa)
               CALL update_fnsb_maxwell
               CALL update_fnsb 
               CALL mtx_reset_communicator
            END IF
!           end of update FNSB

! IF MODEL_DISRUPT=1, FP_COEF is already called in TOP_OF_TIME_LOOP_DISRUPT
            IF(MODEL_DISRUPT.eq.0)THEN 
               IF (MOD(NT,NTSTEP_COEF).EQ.0) CALL FP_COEF(NT)
            END IF

!           sum up SPPF
            IF(MODELS.ne.0)THEN
               CALL mtx_set_communicator(comm_nsa) 
               CALL source_allreduce(SPPF)
               CALL mtx_reset_communicator
            END IF
!           end of sum up SPPF

            CALL GUTIME(gut_coef1)
            GUT_COEF = GUT_COEF + (gut_coef1-gut_cale7) + (gut_coef2-gut_loop1)

         END DO ! END OF DOWHILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF(MODEL_BULK_CONST.ge.1)THEN
!        Bulk f is replaced by Maxwellian
            DO NSA=1, NSAMAX
               DO NR=1, NRMAX
                  NS=NS_NSA(NSA)
                  tempn(NR,NS)=RN_TEMP(NR,NS)
                  tempt(NR,NS)=RT_TEMP(NR,NS)
               END DO
            END DO

            DO NSA=NSASTART, NSAEND
               NS=NS_NSA(NSA)
               DO NR=NRSTARTW, NRENDWM
                  RHON=RM(NR)
                  CALL PL_PROF(RHON,PLF) ! bulk values are fixed to initial values
                  RN_TEMP(NR,NS)=PLF(NS)%RN
                  RT_TEMP(NR,NS)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
                  IF(MODEL_DELTA_F(NS).eq.0)THEN
                     DO NP=NPSTARTW, NPENDWM
                        IF(NP.le.NP_BULK(NR,NSA))THEN
                           DO NTH=1, NTHMAX
                              FL=FPMXWL_EXP(PM(NP,NS),NR,NS)
                              FNS0(NTH,NP,NR,NSA)=FL
                              FNSP(NTH,NP,NR,NSA)=FL
                           END DO
                        END IF
                     END DO
                  ELSE ! delta_f mode: update f_M and delta_f
                     DO NP=NPSTARTW, NPENDWM
                        IF(NP.le.NP_BULK(NR,NSA).and.MODEL_BULK_CONST.eq.1)THEN ! eliminate delta f in bulk region
                           DO NTH=1, NTHMAX
                              FNSP_DEL(NTH,NP,NR,NSA)=0.D0
                           END DO
                        END IF
                        DO NTH=1, NTHMAX
                           FL=FPMXWL_EXP(PM(NP,NS),NR,NS)
                           FNSP_MXWL(NTH,NP,NR,NSA)=FL
                        END DO
                        DO NTH=1, NTHMAX
                           FNS0(NTH,NP,NR,NSA)=FNSP_DEL(NTH,NP,NR,NSA)+FNSP_MXWL(NTH,NP,NR,NSA)
                           FNSP(NTH,NP,NR,NSA)=FNSP_DEL(NTH,NP,NR,NSA)+FNSP_MXWL(NTH,NP,NR,NSA)
                        END DO
                     END DO
                  END IF

               END DO
            END DO

            DO NSA=1, NSAMAX
               DO NR=1, NRMAX
                  NS=NS_NSA(NSA)
                  RN_TEMP(NR,NS)=tempn(NR,NS)
                  RT_TEMP(NR,NS)=tempt(NR,NS)
               END DO
            END DO
         END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF(MODELS.ne.0)THEN
            CALL ALLREDUCE_NF_RATE
            IF(NRANK.eq.0)THEN
!               WRITE(25,'(A,2I5)') "# ", NRANK, NSASTART
               WRITE(25,'(20E14.6)') TIMEFP+DELT, (RATE_NF(1,i),i=1,6), (RATE_NF_BB(1,i),i=1,6)
            END IF
         END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO NSA=NSASTART,NSAEND
            DO NR=NRSTART,NREND
               DO NP=NPSTARTW,NPENDWM
                  DO NTH=1,NTHMAX
                     F(NTH,NP,NR)=FNSP(NTH,NP,NR,NSA) ! used at fpweight only!
                  END DO
               END DO
            END DO
            CALL FPWEIGHT(NSA,IERR) ! SET FPWEIGHT FOR FPSAVE
         END DO

         CALL GUTIME(gut_loop2)
         GUT_1step = gut_loop2-gut_loop1

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

         CALL GUTIME(gut1)
         TIMEFP=TIMEFP+DELT

         CALL GUTIME(gut_out1)
         IF(MODEL_DISRUPT.eq.1)THEN
            CALL FILE_OUTPUT_DISRUPT(NT,IP_all_FP)
         END IF
         CALL GUTIME(gut_out2)
         gut_out=gut_out2-gut_out1
         IF (MOD(NT,NTG1STEP).EQ.0) THEN
            IF(NRANK.eq.0) WRITE(*,'(A,E14.6)') "--------FILE_OUTPUT_TIME in NT LOOP=",gut_out
         END IF

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
!         CALL mtx_broadcast_real8(RT_T,NRMAX*NSAMAX)
         CALL mtx_broadcast_real8(RNS,NRMAX*NSAMAX)
         CALL mtx_broadcast1_integer(NTG1)
         CALL mtx_broadcast1_integer(NTG2)


         DO NSA=1,NSAMAX
            DO NR=1,NRMAX
               IF(RNS(NR,NSA).lt.0)THEN
                  ierr_g = ierr_g + 1
                  IF(NRANK.eq.0)THEN
                     WRITE(*,'(A,2I5)') "NEGATIVE DENS. at NR= ", NR, NSA 
                  END IF
               END IF
            END DO
         END DO

!         IF(NRANK.EQ.0.AND.NTG1.GT.0) call FPWRTSNAP
         CALL GUTIME(gut2)
         IF (MOD(NT,NTG1STEP).EQ.0) THEN
            IF(NRANK.eq.0) WRITE(*,'(A,E14.6)') "--------SAVE_TIME=",gut2-gut1
         END IF

         IF(IERR.NE.0) RETURN

         IF(ierr_g.ne.0)THEN
            call mtx_abort(ierr_g)
         END IF

      ENDDO ! END OF NT LOOP

!     +++++ end of time loop +++++


!
      CALL GUTIME(gut1)
      CALL update_fns
      CALL GUTIME(gut2)
!      IF(MODEL_DISRUPT.eq.1) CALL FLUXS_PTH
      IF(NRANK.eq.0) WRITE(6,'(A,E14.6)') "---------TIME UPDATE FNS =",gut2-gut1

      IF(NRANK.eq.0)THEN
         IF(OUTPUT_TXT_F1.eq.1) CALL OUT_TXT_F1
      END IF

      RETURN
      END SUBROUTINE FP_LOOP

!-------------------------------------------------------------

      SUBROUTINE update_radial_f_boundary(NSA)

      IMPLICIT NONE
      integer,intent(in):: NSA
      integer:: NTH, NP

      DO NP=NPSTARTW, NPENDWM
         DO NTH=1, NTHMAX
            FS2(NTH,NP,NSA) = 2.D0*FS1(NTH,NP,NSA) - FNSP(NTH,NP,NRMAX,NSA)
         END DO
      END DO

      END SUBROUTINE update_radial_f_boundary

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
