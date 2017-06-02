!     $Id: fploop.f90,v 1.40 2013/02/08 07:36:24 nuga Exp $

! *****************
!     MAIN LOOP
! *****************

      MODULE fploop_exp

      USE fpcomm
      use fpexec
!      use fpdrexec
      use fpcoef
      use fpsave
      use libmpi
      use fpmpi
      use fpdisrupt
      USE EG_READ

      contains

!-----------------------------

      SUBROUTINE FP_LOOP_EXP

      USE libmtx
      USE fploop
      USE FPMPI
      USE fpprep, only: Coulomb_log 
      IMPLICIT NONE
      real(kind8),dimension(NSAMAX)::RSUMF,RSUMF0,RSUM_SS
      real(kind8):: RSUMF_, RSUMF0_, FL

      integer:: NT, NR, NP, NTH, NSA, NSBA, NS
      integer:: IERR
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
!      real(8):: sigma, ip_all, ip_ohm, ip_run, jbs, IP_bs, l_ind, IP_prim, 
      real(8):: pitch_angle_av, beam_peak_value
      integer:: NP_2e_h, NP_2e_l, NP_1e_h, NP_1e_l, NP_half_l, NP_half_h, NP_B_peak

!     +++++ Time loop +++++

      DO NT=1,NTMAX
         CALL MAKE_EXP_PROF(timefp)

         N_IMPL=0
         DEPS=1.D0
         DO NSA=NSASTART,NSAEND
            NSBA=NSB_NSA(NSA)
            DO NR=NRSTART-1,NREND+1 ! local
               IF(NR.ge.1.and.NR.le.NRMAX)THEN
                  DO NP=NPSTARTW,NPENDWM
                     IF(MODEL_DELTA_F(NSA).eq.0)THEN
                        DO NTH=1,NTHMAX
                           FNSM(NTH,NP,NR,NSBA)=FNSP(NTH,NP,NR,NSBA) ! minus step: invariant during N_IMPL 
                        END DO
                     ELSEIF(MODEL_DELTA_F(NSA).eq.1)THEN
                        DO NTH=1,NTHMAX
                           FNSM(NTH,NP,NR,NSBA)=FNSP_DEL(NTH,NP,NR,NSBA) ! minus step: invariant during N_IMPL 
                        END DO
                     END IF
                  END DO
               END IF
            END DO
         END DO
 
         gut_EX = 0.0
         gut_conv=0.0
         gut_COEF= 0.0
         gut_1step= 0.0
         gut_cale=0.0
         CALL GUTIME(gut_loop1)

         CALL GUTIME(gut_coef2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

               IF(MODEL_DELTA_F(NSA).eq.1)THEN 
                  DO NR=NRSTART,NREND
                     DO NP=NPSTARTW,NPENDWM
                        DO NTH=1,NTHMAX
                           FNSP(NTH,NP,NR,NSBA)=FNSP_DEL(NTH,NP,NR,NSBA)
                        END DO
                     END DO
                  END DO
               END IF

               ISW_D=1
               CALL GUTIME(gut_exe1)
               IF(ISW_D.eq.0)THEN ! separate p, r ! usually not use this
                  CALL fp_exec(NSA,IERR,its) ! F1 and FNS0 changed
               ELSEIF(ISW_D.eq.1)THEN !
                  IF(MODEL_connor_fp.eq.1.or.MODEL_DISRUPT.eq.0)THEN ! Connor model doesn't require f evolution
                     CALL fp_exec(NSA,IERR,its) ! F1 and FNS0 changed
                  END IF
                  IERR=0
               END IF
               nt_init=1

               IF(MODEL_DELTA_F(NSA).eq.1)THEN
                  DO NR=NRSTART,NREND
                     DO NP=NPSTARTW,NPENDWM
                        DO NTH=1,NTHMAX
                           FNSP(NTH,NP,NR,NSBA)=FNSP_DEL(NTH,NP,NR,NSBA)+FNSP_MXWL(NTH,NP,NR,NSBA)
                           FNSP_DEL(NTH,NP,NR,NSBA)=FNS0(NTH,NP,NR,NSBA)
                           FNS0(NTH,NP,NR,NSBA)=FNSP_DEL(NTH,NP,NR,NSBA)+FNSP_MXWL(NTH,NP,NR,NSBA)
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
                         +ABS(FNSP(NTH,NP,NR,NSBA)-FNS0(NTH,NP,NR,NSBA) )**2
                  RSUMF0(NSA)=RSUMF0(NSA) &
                         +ABS(FNSP(NTH,NP,NR,NSBA))**2
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
!---------- end of convergence check
            CALL GUTIME(gut_conv3)
            gut_conv = gut_conv + (gut_conv3-gut_exe2)

            IF(NRANK.eq.0.and.DEPS.le.EPSFP)THEN
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
!        Bulk f is replaced by Maxwellian
         CALL Define_Bulk_NP
!        temporary
!         DO NSA=NSASTART, NSAEND
!            NS=NS_NSA(NSA)
!            DO NR=NRSTART, NREND
!               RN_READ(NR,NS)=RN_TEMP(NR,NS)
!               RT_READ(NR,NS)=RT_TEMP(NR,NS)
!!               WRITE(*,'(A,2I5,2E14.6)') "TEST READ ", NR, NS, RN_READ(NR,NS), RT_READ(NR,NS)
!            END DO
!         END DO
!
         DO NSA=NSASTART, NSAEND
            NS=NS_NSB(NSA)
            IF(MODEL_DELTA_F(NSA).eq.0)THEN
               DO NR=NRSTART, NREND
                  DO NP=NPSTARTW, NPENDWM
                     IF(NP.le.NP_BULK(NR,NSA))THEN
                        DO NTH=1, NTHMAX
                           FL=FPMXWL_EXP(PM(NP,NSA),NR,NS)
                           FNSP(NTH,NP,NR,NSA)=FL
                        END DO
                     END IF
                  END DO
               END DO
            ELSEIF(MODEL_DELTA_F(NSA).eq.1)THEN
               DO NR=NRSTART, NREND
                  DO NP=NPSTARTW, NPENDWM
                     IF(NP.le.NP_BULK(NR,NSA))THEN
                        DO NTH=1, NTHMAX
                           FNSP_DEL(NTH,NP,NR,NSA)=0.D0
                        END DO
                     ELSE
                        DO NTH=1, NTHMAX
                           FL=FPMXWL_EXP(PM(NP,NSA),NR,NS)
                           FNSP_MXWL(NTH,NP,NR,NSA)=FL
                        END DO
                     END IF
                     DO NTH=1, NTHMAX
                        FNS0(NTH,NP,NR,NS)=FNSP_DEL(NTH,NP,NR,NS)+FNSP_MXWL(NTH,NP,NR,NS)
                        FNSP(NTH,NP,NR,NS)=FNSP_DEL(NTH,NP,NR,NS)+FNSP_MXWL(NTH,NP,NR,NS)  
                     END DO
                  END DO
               END DO
            END IF
         END DO
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
         CALL mtx_broadcast_real8(RT_T,NRMAX*NSAMAX)
         CALL mtx_broadcast_real8(RNS,NRMAX*NSAMAX)
         CALL mtx_broadcast1_integer(NTG1)
         CALL mtx_broadcast1_integer(NTG2)


         DO NSA=1,NSAMAX
            DO NR=1,NRMAX
               IF(RNS(NR,NSA).lt.0)THEN
                  ierr_g = ierr_g + 1
                  IF(NRANK.eq.0)THEN
                     WRITE(*,'(A,3I5,E14.6)') "NEGATIVE DENS. at NR= ", NR, NSA, ierr_g, RNS(NR,NSA)
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
!      IF(MODEL_DISRUPT.eq.1) CALL FLUXS_PTH
      IF(NRANK.eq.0) WRITE(6,'(A,E14.6)') "---------TIME UPDATE FNS =",gut2-gut1

!      IF(NRANK.eq.0.and.MODEL_DISRUPT.ne.0)THEN
      IF(NRANK.eq.0)THEN
         NR=1
         DO NP=1,NPMAX
!pitch angle average
            pitch_angle_av = 0.D0
            DO NTH=1,NTHMAX
               pitch_angle_av = pitch_angle_av + FNS(NTH,NP,NR,NS_F1)*SINM(NTH)*0.5D0
            END DO
            pitch_angle_av = pitch_angle_av/NTHMAX
            WRITE(9,'(1PE12.4,I6,1P30E17.8e3)') PTG(NTG1)*1000, NP, PM(NP,NS_F1), &
                 PM(NP,NS_F1)*PTFP0(NS_F1)/AMFP(NS_F1)/VC/SQRT(1.D0+PM(NP,NS_F1)**2*THETA0(NS_F1)), &
                 PM(NP,NS_F1)**2, &
                 0.5D0*PTFP0(NS_F1)**2*PM(NP,NS_F1)**2/(AEE*AMFP(NS_F1)*1.D3), FNS(1,NP,NR,NS_F1), FNS(NTH_F1,NP,NR,NS_F1), &
                 FNS(NTHMAX/2,NP,NR,NS_F1), pitch_angle_av!, &
         END DO
         WRITE(9,*) " "
         WRITE(9,*) " "

         NP_2e_l=1
         NP_2e_h=1
         NP_1e_l=1
         NP_1e_h=1
         NP_half_l=1
         NP_half_h=1
         beam_peak_value=0.D0
         DO NP=1, NPMAX-1
            IF(FNS(NTH_F1,NP,NR,NS_F1).le.2.D-5.and.FNS(NTH_F1,NP+1,NR,NS_F1).ge.2.D-5)THEN
               NP_2e_l=NP
            ELSEIF(FNS(NTH_F1,NP,NR,NS_F1).ge.2.D-5.and.FNS(NTH_F1,NP+1,NR,NS_F1).le.2.D-5)THEN
               NP_2e_h=NP
            END IF
            IF(FNS(NTH_F1,NP,NR,NS_F1).le.5.D-6.and.FNS(NTH_F1,NP+1,NR,NS_F1).ge.5.D-6)THEN
               NP_1e_l=NP
            ELSEIF(FNS(NTH_F1,NP,NR,NS_F1).ge.5.D-6.and.FNS(NTH_F1,NP+1,NR,NS_F1).le.5.D-6)THEN
               NP_1e_h=NP
            END IF
         END DO

         DO NP=NP_BULK(NR,NS_F1), NPMAX-1 ! half value width
            IF(FNS(NTH_F1,NP-1,NR,NS_F1).le.FNS(NTH_F1,NP,NR,NS_F1).and.FNS(NTH_F1,NP,NR,NS_F1).ge.FNS(NTH_F1,NP+1,NR,NS_F1))THEN
               beam_peak_value=FNS(NTH_F1,NP,NR,NS_F1)*0.5D0
               NP_B_peak=NP
            END IF
         END DO
         DO NP=1, NPMAX-1
            IF(FNS(NTH_F1,NP,NR,NS_F1).le.beam_peak_value.and.FNS(NTH_F1,NP+1,NR,NS_F1).ge.beam_peak_value)THEN
               NP_half_l=NP
            ELSEIF(FNS(NTH_F1,NP,NR,NS_F1).ge.beam_peak_value.and.FNS(NTH_F1,NP+1,NR,NS_F1).le.beam_peak_value)THEN
               NP_half_h=NP
            END IF
         END DO

         WRITE(24,'(E14.6,3I5,20E14.6)') TIMEFP, &
!              NP_2e_l, NP_2e_h, &
              NP_half_l, NP_half_h, &
              NP_B_peak, &
!              PTFP0(2)**2*PM(NP_2e_l,2)**2/(AEE*AMFP(2)*1.D3)*0.5D0, &
!              PTFP0(2)**2*PM(NP_2e_h,2)**2/(AEE*AMFP(2)*1.D3)*0.5D0, &
!              FNS(NTH_F1,NP_2e_l,NR,2), FNS(NTH_F1,NP_2e_h,NR,2), &
!              PTFP0(2)**2*PM(NP_1e_l,2)**2/(AEE*AMFP(2)*1.D3)*0.5D0, &
!              PTFP0(2)**2*PM(NP_1e_h,2)**2/(AEE*AMFP(2)*1.D3)*0.5D0, &
!              FNS(NTH_F1,NP_1e_l,NR,2), FNS(NTH_F1,NP_1e_h,NR,2), &
              PTFP0(NS_F1)**2*PM(NP_half_l,NS_F1)**2/(AEE*AMFP(NS_F1)*1.D3)*0.5D0, &
              PTFP0(NS_F1)**2*PM(NP_half_h,NS_F1)**2/(AEE*AMFP(NS_F1)*1.D3)*0.5D0, &
              PTFP0(NS_F1)**2*PM(NP_B_peak,NS_F1)**2/(AEE*AMFP(NS_F1)*1.D3)*0.5D0, &
              FNS(NTH_F1,NP_half_l,NR,NS_F1), FNS(NTH_F1,NP_half_h,NR,NS_F1), beam_peak_value*2.D0
         WRITE(*,'(A,E14.6,6I5,E14.6,I5)') "BEAM_HALF ",TIMEFP, NP_2e_l, NP_2e_h, NP_1e_l, NP_1e_h, NP_half_l, NP_half_h, beam_peak_value, NP_BULK(NR,NS_F1)
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
      END SUBROUTINE FP_LOOP_EXP

! ! ****************************************
! !     MAXWELLIAN VELOCITY DISTRIBUTION
! ! ****************************************

!       FUNCTION FPMXWL_EXP(PML,NR,NS)

!       USE plprof
!       implicit none
!       integer :: NR, NS
!       real(kind8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
!       real(kind8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
!       TYPE(pl_plf_type),DIMENSION(NSMAX):: plf
!       real(kind8):: FPMXWL_EXP

!       AMFDL=PA(NS)*AMP
!       AEFDL=PZ(NS)*AEE
!       RNFD0L=PN(NS)
!       RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
!       PTFD0L=SQRT(RTFD0L*1.D3*AEE*AMFDL)

!       IF(NR.eq.0)THEN
!          RL=0.D0
!          RHON=ABS(RL)
!       ELSEIF(NR.EQ.NRMAX+1) THEN
!          RL=RM(NRMAX)+DELR
!          RHON=MIN(RL,1.D0)
!       ELSE
!          RL=RM(NR)
!          RHON=RL
!       ENDIF
!       CALL PL_PROF(RHON,PLF)

!       RNFDL=RN_TEMP(NR,NS)/RNFD0L
!       RTFDL=RT_TEMP(NR,NS)

!       IF(MODELR.EQ.0) THEN
!          FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0L)**3
!          EX=PML**2/(2.D0*RTFDL/RTFD0L)
!          FPMXWL_EXP=FACT*EXP(-EX)
!       ELSE
!          THETA0L=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
!          THETAL=THETA0L*RTFDL/RTFD0L
!          Z=1.D0/THETAL
!          DKBSL=BESEKN(2,Z)
!          FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL) &
!               *RTFD0L
!          EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
!          FPMXWL_EXP=FACT*EXP(EX)
!       END IF

!       RETURN
!       END FUNCTION FPMXWL_EXP
!------------------------------------
      END MODULE fploop_exp
