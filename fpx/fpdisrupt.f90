! fpdisrupt.f90

MODULE fpdisrupt

     USE fpcomm

      REAL(rkind),parameter:: Z_ef=3.D0
!      REAL(rkind),parameter:: Z_ef=5.D0
!      REAL(rkind),parameter:: Z_ef=1.D0
!      REAL(rkind),parameter:: Z_ef=3.D0

!      REAL(rkind),parameter:: IP_init=1.8D6 ! [A] ! JET
!      REAL(rkind),parameter:: IP_init=1.9D6 ! [A] ! JET 2008
!     REAL(rkind),parameter:: IP_init=1.052D6 ! [A] !60 ITB
!     REAL(rkind),parameter:: IP_init=1.5D6 ! [A] !60
!     REAL(rkind),parameter:: IP_init=0.8D6 ! [A] !60 ITB
!     REAL(rkind),parameter:: IP_init=0.5D6 ! [A] !60 ITB
!     REAL(rkind),parameter:: IP_init=1.04D6 ! [A] !60 49255
!      REAL(rkind),parameter:: IP_init=2.D6 ! [A] ! 60 typical
!      REAL(rkind),parameter:: IP_init=15.D6 ! [A] ! ITER
!     REAL(rkind),parameter:: IP_init=1.0D3 ! [A] !test
!     REAL(rkind):: IP_init=RIP*1.D6 ! [A] !60 ITB

      integer,parameter:: ISW_Z=1 !0=Z_ef, 1=ZEFF

      REAL(rkind),dimension(:),pointer:: rt_init
!      REAL(rkind),parameter:: lnL_ED=10.D0
      REAL(rkind),parameter:: lnL_ED=18.D0
      integer,parameter:: ISW_Q=2 ! 2=exponential, 3=2step quench, 0=linear
      contains

! ------------------------------------------

      SUBROUTINE set_initial_disrupt_param 

      USE plprof
      USE fpmpi
      USE libmpi
      IMPLICIT NONE

      INTEGER:: NR, NSA, NSB, NS, NP, NSW, N, NTH
      TYPE(pl_prf_type),DIMENSION(NSMAX):: PLF
      REAL(rkind),dimension(NRSTART:NREND):: &
           RN1_temp, RN2_temp, RJ1_temp, RJ2_temp
      REAL(rkind),dimension(NRSTART:NREND):: RT1_temp, RT2_temp
      REAL(rkind):: RHON, SUMIP, IP_norm, jbs, SUMIP_BS, sumip_ohm, T0_init
      INTEGER:: VLOC, n_ite
      REAL(rkind):: Erunm, Erunp, minimum_RE_E

      allocate(rt_init(NRSTART:NREND))
!------critical momentum
!      IF(MODELR.eq.0)THEN
!         pc_runaway=1.D0/SQRT(E0)
!      ELSEIF(MODELR.eq.1)THEN
!         pc_runaway=1.D0/SQRT(E0-THETA0(1))
!      END IF
!      IF(pc_runaway.ge.PG(NPMAX+1,1))THEN
!         NPC_runaway=NPMAX+1
!      ELSE
!         DO NP=NPMAX+1,1,-1
!            IF(PG(NP,1).ge.pc_runaway) NPC_runaway=NP
!         END DO
!      END IF
!----- SET NPC_runaway = 1MeV
      IF(MODEL_RE_pmax.eq.1)THEN
         CALL mtx_set_communicator(comm_np)
         minimum_RE_E=1.D3 ! keV
         NPC_runaway=0
         DO NP=NPSTART, NPEND
            Erunm=PTFP0(1)**2*PM(NP,1)**2/(AEE*AMFP(1)*1.D3)
            IF(Erunm.le.minimum_RE_E)THEN
               Erunp=PTFP0(1)**2*PM(NP+1,1)**2/(AEE*AMFP(1)*1.D3)
               IF(Erunp.ge.minimum_RE_E)THEN
                  NPC_runaway = NP + 1
               END IF
            END IF
         END DO
         CALL mtx_allreduce1_integer(NPC_runaway,3,NPC_runaway,vloc)
         CALL mtx_reset_communicator
         IF(NRANK.eq.0) &
              WRITE(6,'(A,I5,A,E14.6)') "RE boundary NP=",NPC_runaway, " PM(NP,1)=",PM(NPC_runaway,1)
      END IF
!----- SET END

      DO NR=1,NRMAX
         RN_disrupt(NR)=0.D0
         RN_runaway(NR)=0.D0
         RN_drei(NR)=0.D0
         Rj_ohm(NR)=0.D0
         RJ_runaway(NR)=0.D0
         RJ_bs(NR)=0.D0
      END DO
      RJ_bsm(:)=0.D0
!-----initial profile of n, j, T
      SUMIP=0.D0
      SUMIP_BS=0.D0
      DO NR=NRSTART,NREND
         RN1_TEMP(NR)=RNFP(NR,1)
         RN2_TEMP(NR)=0.D0

         IF(MODEL_jfp.eq.0)THEN
!            RJ1_temp(NR)=(1.D0-RM(NR)**0.95)**3 ! SMITH2006 ! SWITCH
!            RJ1_temp(NR)=(1.D0-RM(NR)**2) ! parabora
!            RJ1_temp(NR)=(1.D0-RM(NR)**1.74D0)**3.23D0 ! matsuyama 60U
            RJ1_temp(NR)=(1.D0-RM(NR)**RJPROF1)**RJPROF2 ! RJPROF2
         ELSE
            RJ1_temp(NR)=(1.D0-RM(NR)**1.5)**3 ! matsuyama 60U
         END IF
         RJ2_temp(NR)=0.D0
         SUMIP=RJ1_temp(NR)*VOLR(NR)/(2.D0*PI*RR) + SUMIP

         RHON=RM(NR)
         CALL PL_PROF(RHON,PLF)
         T0_init=2.D0
         RT1_temp(NR)=RTFD(NR,1)
         RT_INIT(NR)=RTFD(NR,1)
!         RT2_temp(NR)=T0_quench ! flat profile
         RT2_temp(NR)=T0_quench*(1.D0-0.9D0*RM(NR)**2) ! smith
!         RT2_temp(NR)=T0_quench*(1.D0-0.1D0*RM(NR)**2) ! nearly flat
      END DO

      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(RN1_temp,NREND-NRSTART+1,RN_disrupt)
      call mtx_allgather_real8(RN2_temp,NREND-NRSTART+1,RN_runaway)

      CALL mtx_allreduce1_real8(SUMIP,3,IP_norm,vloc)
      DO NR=NRSTART,NREND
         RJ1_temp(NR)=RJ1_temp(NR)*RIP*1.D6/IP_norm*1.D-6
      END DO
      call mtx_allgather_real8(RJ1_temp,NREND-NRSTART+1,Rj_ohm)
      call mtx_allgather_real8(RJ2_temp,NREND-NRSTART+1,RJ_runaway)
!
      call mtx_allgather_real8(RT1_temp,NREND-NRSTART+1,RT_quench)
      call mtx_allgather_real8(RT2_temp,NREND-NRSTART+1,RT_quench_f)

      CALL mtx_reset_communicator

      DO NSA=1,NSAMAX
         E_drei0(NSA)=PTFP0(NSA)/(ABS(AEFP(NSA))*tau_ta0(NSA) ) 
         E_crit0(NSA)=E_drei0(NSA)*THETA0(NSA) 
      END DO

      DO NR=1,NRMAX
         Rconnor(NR)=0.D0
         RFP_AVA(NR)=0.D0
         RFP(NR)=0.D0
      END DO
      DO NR=NRSTART,NREND
         RFPL(NR)=0.D0
      END DO

      IF(MODEL_BS.ge.1)THEN
         n_ite=0
         DO while(N_ite.le.100)
            SUMIP_BS=0.D0
            SUMIP=0.D0
            SUMIP_ohm=0.D0
            DO NR=1,NRMAX
               call bootstrap_current(NR,jbs)
!               RJ_bs(NR)=jbs
               RJ_bs(NR)=jbs*0.D0
               SUMIP_ohm= RJ_ohm(NR)*VOLR(NR)/(2.D0*PI*RR) + SUMIP_ohm
               SUMIP_BS=RJ_bs(NR)*VOLR(NR)/(2.D0*PI*RR) + SUMIP_BS
            END DO
            DO NR=1,NRMAX
!               RJ_ohm(NR) = (RIP*1.D6/SUMIP_BS*1.D-6 -1.D0)* RJ_bs(NR) ! change j_ohm prof. to j_bs
               RJ_ohm(NR) = (RIP -SUMIP_BS)/(SUMIP_ohm)*RJ_ohm(NR) ! keep j_ohm prof. ! SWITCH
               SUMIP= (RJ_ohm(NR)+RJ_bs(NR))*VOLR(NR)/(2.D0*PI*RR) + SUMIP
               QLM(NR)=2*PI*RM(NR)**2*RA**2*BB/(RMU0*SUMIP*1.D6*RR)
            END DO
            n_ite=n_ite+1
         END DO
      END IF
      DO NTH=1,NTHMAX
         RE_PITCH(NTH)=0.D0
      END DO
      IF(NRANK.eq.0) &
           WRITE(6,'(A,1P3E14.6)') "initial IP_all, IP_ohm, IP_bs=", SUMIP, SUMIP-SUMIP_BS, SUMIP_BS


      END SUBROUTINE set_initial_disrupt_param
!**************************************************
      SUBROUTINE set_post_disrupt_Clog_f

      IMPLICIT NONE
      INTEGER:: NR, NSA, NSB, NS, NSFP, NSFD, ISW_CLOG
      REAL(rkind):: RNA,RNB,RTA,RTB,RLNRL
      REAL(rkind),dimension(NSAMAX):: PTFP0_f

      DO NR=NRSTART, NRENDWM
!      DO NR=1, NRMAX

!-----   post disrupt Coulomb log
         ISW_CLOG=0 ! =0 Wesson, =1 NRL
         DO NSA=1,NSAMAX
            NSFP=NS_NSB(NSA)
            IF(MODEL_IMPURITY.eq.0)THEN
               RNA=RNFP(NR,NSA)
            ELSE
!               RNA=RNFP(NR,NSA)+PZ(N_impu)*SPITOT*RNFP(NR,NSA)/RNFP0(NSA)
               RNA=SPITOT*RNFP(NR,NSA)/RNFP0(NSA)
            END IF
            RTA=RT_quench_f(NR)
            PTFP0_f(NSA)=SQRT(T0_quench*1.D3*AEE*AMFP(NSA))
            DO NSB=1,NSBMAX
               NSFD=NS_NSB(NSB)
               IF(MODEL_IMPURITY.eq.0)THEN
                  RNB=RNFD(NR,NSB)
               ELSE
                  IF(NSB.ne.N_impu)THEN
                     RNB=RNFD(NR,NSB)
                  ELSE
!                     RNB=RNFD(NR,NSB)+SPITOT*RNFD(NR,NSB)/RNFD0(NSB)
                     RNB=((SPITOT-RNFP0(1))*RNFP(NR,1)/RNFP0(1))/zeff_imp
                  END IF
               END IF
               RTB=RT_quench_f(NR)
               IF(ISW_CLOG.eq.0)THEN
                  IF(PZ(NSFP).eq.-1.and.PZ(NSFD).eq.-1) THEN !e-e
                     RLNRL=14.9D0-0.5D0*LOG(RNA)+LOG(RTA) 
                  ELSEIF(PZ(NSFP).eq.-1.OR.PZ(NSFD).eq.-1) THEN
                     RLNRL=15.2D0-0.5D0*LOG(RNA)+LOG(RTA) ! e-i T>10eV
                  ELSE
                     RLNRL=17.3D0-0.5D0*LOG(RNA)+1.5D0*LOG(RTB) ! i-i T < m_i/m_p*10 keV, single charge
                  ENDIF
               ELSEIF(ISW_CLOG.eq.1)THEN
                  IF(PZ(NSFP).eq.-1.and.PZ(NSFD).eq.-1) THEN !e-e
                     RLNRL=23.5D0-LOG(SQRT(RNA*1.D14)*(RTA*1.D3)**(-1.25D0) ) - &
                          SQRT(1.D-5+(LOG(RTA*1.D3)-2.D0)**2/16.D0 )
                  ELSEIF(PZ(NSFP).eq.-1.OR.PZ(NSFD).eq.-1) THEN
                     RLNRL=24.D0-LOG(SQRT(RNA*1.D14)*(RTA*1.D3)**(-1) ) ! Ti*(me/mi) < 10eV < Te
                  ELSE
                     RLNRL=23.D0-LOG(PZ(NSA)*PZ(NSB)*(PA(NSA)+PA(NSB)) &
                          /(PA(NSA)*(RTB*1.D3)+PA(NSB)*(RTA*1.D3))* &
                          SQRT((RNA*1.D14)*PZ(NSA)**2/(RTA*1.D3) + (RNB*1.D14)*PZ(NSB)**2/(RTB*1.D3) ) )
                  ENDIF
               END IF
               POST_LNLAM_f(NR,NSB,NSA)=RLNRL
            END DO ! NSB
            IF(NRANK.eq.0)THEN
               IF(MODEL_IMPURITY.eq.0)THEN
                  POST_tau_ta0_f(NSA)=(4.D0*PI*EPS0**2)*PTFP0_f(NSA)**3 &
                       /( AMFP(NSA)*AEFP(NSA)**4*POST_LNLAM_f(1,NSA,NSA)*RNFP0(NSA)*1.D20 )
               ELSE
                  POST_tau_ta0_f(NSA)=(4.D0*PI*EPS0**2)*PTFP0_f(NSA)**3 &
                       /( AMFP(NSA)*AEFP(NSA)**4*POST_LNLAM_f(1,NSA,NSA)*(SPITOT)*1.D20 )
               END IF
            END IF
         END DO ! NSA

      END DO ! NR

      END SUBROUTINE set_post_disrupt_Clog_f
!**************************************************
      SUBROUTINE set_post_disrupt_Clog

      USE fpmpi
      USE libmpi
      IMPLICIT NONE
      INTEGER:: NR, NSA, NSB, NS, NSFP, NSFD, ISW_CLOG, N, NSW
      REAL(rkind):: RNA,RNB,RTA,RTB,RLNRL,P_thermal, rtheta
      REAL(rkind),dimension(NSAMAX):: PTFP0_f
      REAL(rkind),dimension(NRSTART:NREND,NSAMAX):: tau_ta
      REAL(rkind),dimension(NRSTART:NREND):: LNL_l

      IF(MODEL_LNL.eq.0)THEN ! update Coulomb Log with changing T
      DO NR=NRSTART, NREND
!-----   post disrupt Coulomb log
         ISW_CLOG=0 ! =0 Wesson, =1 NRL
         DO NSA=1,NSAMAX
            NSFP=NS_NSB(NSA)
            IF(MODEL_IMPURITY.eq.0)THEN
               RNA=RNFP(NR,NSA)
            ELSE
               RNA=RN_MGI(NR,NSA)
            END IF
            RTA=RT_quench(NR)
            DO NSB=1,NSBMAX
               NSFD=NS_NSB(NSB)
               IF(MODEL_IMPURITY.eq.0)THEN
                  RNB=RNFD(NR,NSB)
               ELSE
                  RNB=RN_MGI(NR,NSB)
               END IF
               RTB=RT_quench(NR)
               IF(ISW_CLOG.eq.0)THEN
                  IF(PZ(NSFP).eq.-1.and.PZ(NSFD).eq.-1) THEN !e-e
                     RLNRL=14.9D0-0.5D0*LOG(RNA)+LOG(RTA) 
                  ELSEIF(PZ(NSFP).eq.-1.OR.PZ(NSFD).eq.-1) THEN
                     RLNRL=15.2D0-0.5D0*LOG(RNA)+LOG(RTA) ! e-i T>10eV
                  ELSE
                     RLNRL=17.3D0-0.5D0*LOG(RNA)+1.5D0*LOG(RTB) ! i-i T < m_i/m_p*10 keV, single charge
                  ENDIF
               ELSEIF(ISW_CLOG.eq.1)THEN
                  IF(PZ(NSFP).eq.-1.and.PZ(NSFD).eq.-1) THEN !e-e
                     RLNRL=23.5D0-LOG(SQRT(RNA*1.D14)*(RTA*1.D3)**(-1.25D0) ) - &
                          SQRT(1.D-5+(LOG(RTA*1.D3)-2.D0)**2/16.D0 )
                  ELSEIF(PZ(NSFP).eq.-1.OR.PZ(NSFD).eq.-1) THEN
                     RLNRL=24.D0-LOG(SQRT(RNA*1.D14)*(RTA*1.D3)**(-1) ) ! Ti*(me/mi) < 10eV < Te
                  ELSE
                     RLNRL=23.D0-LOG(PZ(NSA)*PZ(NSB)*(PA(NSA)+PA(NSB)) &
                          /(PA(NSA)*(RTB*1.D3)+PA(NSB)*(RTA*1.D3))* &
                          SQRT((RNA*1.D14)*PZ(NSA)**2/(RTA*1.D3) + (RNB*1.D14)*PZ(NSB)**2/(RTB*1.D3) ) )
                  ENDIF
               END IF
               POST_LNLAM(NR,NSB,NSA)=RLNRL
            END DO ! NSB
         END DO ! NSA
      END DO ! NR
      ELSEIF(MODEL_LNL.eq.1)THEN ! fix initial Clog
         DO NSB=1, NSBMAX
            DO NSA=1, NSAMAX
               DO NR=NRSTART,NREND
                  POST_LNLAM(NR,NSB,NSA)=LNLAM(NR,NSB,NSA)
               END DO
            END DO
         END DO
      ELSEIF(MODEL_LNL.eq.2)THEN ! fix post-quech Clog
         DO NSB=1, NSBMAX
            DO NSA=1, NSAMAX
               DO NR=NRSTART,NREND
                  POST_LNLAM(NR,NSB,NSA)=POST_LNLAM_f(NR,NSB,NSA)
               END DO
            END DO
         END DO
      END IF

! e-e collision time
      DO NSA=1,NSAMAX
         DO NR=NRSTART, NREND
            IF(MODEL_IMPURITY.eq.0)THEN
               RNA=RNFP(NR,NSA)
            ELSE
               RNA=RN_MGI(NR,NSA)
            END IF
            RTA=RT_quench(NR)
            P_thermal=SQRT(RTA*1.D3*AEE*AMFP(NSA))
            tau_ta(NR,NSA)=(4.D0*PI*EPS0**2)*P_thermal**3 &
                 /( AMFP(NSA)*AEFP(NSA)**4*POST_LNLAM(NR,NSA,NSA)*RNA*1.D20 )
         END DO
      END DO
      DO NR=NRSTART, NREND
         lnl_l(NR)=POST_LNLAM(NR,1,1)
      END DO

      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(LNL_L,NREND-NRSTART+1,LNL_GL)

      CALL mtx_set_communicator(comm_nsanr)
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL fp_gatherv_real8_sav(tau_ta,SAVLEN(NRANK+1),POST_tau_ta,N,NSA)
      END DO
      CALL mtx_reset_communicator
      CALL mtx_broadcast_real8(POST_tau_ta,NRMAX*NSAMAX)

      DO NR=1,NRMAX
         ER_drei(NR)=SQRT(RT_quench(NR)*1.D3*AEE*AMFP(1))/(ABS(AEFP(1))*POST_tau_ta(NR,1) )
         ER_crit(NR)=ER_drei(NR)*RT_quench(NR)*1.D3*AEE/(AMFP(1)*VC*VC) &
              /lnl_gl(NR)*lnL_ED  ! lnL for knock on collision

!         ER_drei(NR)=SQRT(RT_quench(NR)*1.D3*AEE*AMFP(1))/(ABS(AEFP(1))*POST_tau_ta(NR,1) ) &
!              /lnl_gl(NR)*lnL_ED ! lnL for fixed value
!         ER_crit(NR)=ER_drei(NR)*RT_quench(NR)*1.D3*AEE/(AMFP(1)*VC*VC) ! SWITCH
      END DO

      END SUBROUTINE set_post_disrupt_Clog
!*****************************************
      SUBROUTINE display_disrupt_initials

      IMPLICIT NONE
      REAL(rkind):: alp, z_i, h_alpha_z, lambda_alpha, gamma_alpha_z, G_connor
      REAL(rkind):: G_connor_nr, G_connor_lm, tau_th, v_thermal, tau_rela
      integer:: j

      WRITE(6,*)" ---- DISRUPTION PARAM ------"
      WRITE(6,'(A,1PE14.6)') "Pre  disruption Te0 [keV]=", RTFP0(1)
      WRITE(6,'(A,1P4E14.6)') &
             "Pre  disrupt tau_ta0[sec]=", (tau_ta0(j), j=1,NSAMAX)

      IF(ISW_Z.eq.0)THEN
         z_i = Z_ef 
      ELSEIF(ISW_Z.eq.1)THEN
         z_i = ZEFF
      END IF
      IF(MODELR.eq.0)THEN
         WRITE(6,'(A,1P3E14.6,A,1PE14.6)') &
              "E0*E_drei0=E1[V/m] ->",E0, E_drei0(1), E1(1)&
              ," p_{c_runaway}=",1.D0/SQRT(E0)

!         v_thermal=SQRT( RT_quench(1)*1.D3*AEE/AMFD(1)) 
!         tau_th=(4.D0*PI*EPS0**2)*AMFP(1)**2*v_thermal**3/ &
!              ( AEFP(1)**4*POST_LNLAM(1,1,1)*RNFP0(1)*1.D20 )
         G_connor_nr=0.35D0*E0**(-3.D0*(Z_i+1.D0)/16.D0)&
              *EXP(-0.25D0/E0 -SQRT( (1.D0+z_i)/E0 ) )
         WRITE(6,'(A,E14.6)') "RE gen rate_KB=", G_connor_nr!/tau_th
      ELSE
         WRITE(6,'(A,1P4E14.6)') &
              "E0*E_drei0=E1[V/m], E_crit0 ->",E0, E_drei0(1), E1(1), E_crit0(1)
         WRITE(6,'(A,1PE14.6,A,1PE14.6)') &
              " p_{c_runaway}=",1.D0/SQRT(E0)           &
              ," p_{cr_runaway}=",1.D0/SQRT(E0-THETA0(1))
         
         alp = E0/THETA0(1)
         h_alpha_z=( alp*(z_i+1.D0) - z_i + 7.D0 + &
              2.D0*SQRT(alp/(alp-1.D0))*(1.D0+z_i)*(alp-2.D0) ) &
              /( 16.D0*(alp-1.D0) )
         lambda_alpha=8.D0*alp*(alp-0.5D0-SQRT(alp*(alp-1.D0)) )
         gamma_alpha_z=SQRT( (1.0+z_i)*alp**2/(8.D0*(alp-1.D0)) )&
              *(0.5D0*PI-ASIN(1.D0-2.D0/alp))
         
         G_connor=0.35D0* E0**(-h_alpha_z)*EXP(-0.25D0*lambda_alpha/E0 &
              -SQRT(2.D0/E0)*gamma_alpha_z )
         
         G_connor_nr=0.35D0*E0**(-3.D0*(Z_i+1.D0)/16.D0)&
              *EXP(-0.25D0/E0 -SQRT( (1.D0+z_i)/E0 ) )
         G_connor_lm=G_connor_nr &
              *EXP(-THETA0(1)*(0.125D0/E0**2 + 2.D0/3.D0/SQRT(E0**3)*SQRT(1.D0+z_i) ) )
         
         WRITE(6,'(A,1PE14.6,A,1PE14.6,A,1PE14.6)') " alpha = ", alp, " G_connor= ", G_connor, &
              " THETA0", THETA0(1)
         WRITE(6,'(A,1PE14.6,A,1PE14.6)') " G_Connor_nr = ", &
              G_connor_nr, " G_connor_lm = ", G_connor_lm
      END IF

      IF(MODEL_IMPURITY.eq.0)THEN
         tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
              ( AEFP(1)**4*POST_LNLAM(1,1,1)*RNFP0(1)*1.D20 )
      ELSE
         tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
              ( AEFP(1)**4*POST_LNLAM(1,1,1)*RN0_MGI(1)*1.D20 )
      END IF

      WRITE(6,*) "-----------------------"
      WRITE(6,'(A,1PE14.6)') "decay time [msec] = ", tau_quench*1000
      WRITE(6,*) "-----------------------"

      WRITE(6,'(A,1PE14.6)') "Post disruption Te0 [keV]=", T0_quench
      WRITE(6,'(A,1P4E14.6)') &
             "Post disrupt tau_rela[sec]=", tau_rela
      WRITE(6,'(A,1P4E14.6)') &
             "Post disrupt tau_ta0[sec]=", (POST_tau_ta0_f(j), j=1,NSAMAX)
      WRITE(6,'(A,1P4E14.6)') &
             "Post disrupt tau_ta[sec]=", POST_tau_ta(1,1), ER_drei(1), ER_crit(1)

      IF(MODEL_IMPURITY.eq.1)THEN
         WRITE(6,'(A)') "  MGI PARAM "
         WRITE(6,'(A,I3)') "     MGI ion species = ", N_impu
         WRITE(6,'(A,1PE14.6)') "     target Zeff = ", target_zeff
         WRITE(6,'(A,1PE14.6)') "     target electron dens. = ", SPITOT
      END IF

      WRITE(6,*)" ---- END OF DISRUPTION PARAM ------"
      

      END SUBROUTINE display_disrupt_initials
! ****************************************************
      SUBROUTINE SPITZER_SIGMA(NR,Sigma)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      REAL(rkind),INTENT(OUT):: Sigma
      INTEGER:: NSA, NSB
      REAL(rkind):: fact, taue_col, RTE, RNE, RTI, RNI, Z_i
      REAL(rkind):: neoc,  phi, f_T, C_, tau_rela, theta_l, C_log

      NSA=1
      NSB=2

      IF(ISW_Z.eq.0)THEN
         z_i = Z_ef 
      ELSEIF(ISW_Z.eq.1)THEN
         z_i = ZEFF
      END IF

      IF(MODEL_DISRUPT.eq.0)THEN
         RTE=RTFP(NR,NSA)
         RNE=RNFP(NR,NSA)
         C_log = LNLAM(NR,NSA,NSA)
      ELSE
         IF(NT_init.eq.0)THEN
            RTE=RT_INIT(NR)
            RNE=RNFP(NR,NSA)
         ELSE
            RTE=RT_quench(NR)
            IF(MODEL_IMPURITY.eq.0)THEN
               RNE=RNS(NR,1)
            ELSE
               RNE=RN_MGI(NR,1)
            END IF
         END IF
         C_log = POST_LNLAM(NR,NSA,NSA)
      END IF

!      FACT=AEFP(NSA)**2*AEFD(NSB)**2*POST_LNLAM(NR,NSB,NSA)*RNE*1.D20
      FACT=Z_i*AEFP(NSA)**4*C_log*RNE*1.D20
      taue_col=3.D0*SQRT((2.D0*PI)**3)/FACT &
           *SQRT(AMFP(1)*(AEE*RTE*1.D3)**3)*(EPS0**2)
      sigma=1.96D0*RNE*1.D20*AEFP(NSA)**2*taue_col/AMFP(NSA) ! Wesson P. 174, 737
!      sigma= ! P. 71

!!     NEO-CLASSICAL CORRECTION
!      neoc=(1.D0-SQRT(invasp))**2 ! P. 174     
      IF(MODEL_DISRUPT.eq.0)THEN
         theta_l=THETA0(1)*RTFP(NR,1)/RTFP0(1)
         tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
              ( AEFP(1)**4*LNLAM(NR,1,1)*RNFP(NR,NSA)*1.D20 )
      ELSE
         theta_l=THETA0(1)*RT_quench(NR)/RTFP0(1)
         IF(MODEL_IMPURITY.eq.0)THEN
            tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
                 ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RNFP(NR,NSA)*1.D20 )
         ELSE
            tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
                 ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RN_MGI(NR,NSA)*1.D20 )
         END IF
      END IF

      C_ = 0.56D0/Z_i*(3.0D0-Z_i)/(3.D0+Z_i)
      f_t=1.D0 -(1.D0-EPSRM2(NR))**2 &
               /( SQRT(1.D0-EPSRM2(NR)**2)*(1.D0+1.46D0*SQRT(EPSRM2(NR))) )
!      phi = f_t/(1.D0 + (0.58D0+0.2D0*Z_i)*theta_l**2 &
!                       *(2.D0*RR*QLM(NR)*EPSRM2(NR)**(-1.5D0) ) &
      !                       /(3.D0*SQRT(2.D0*PI)*VC*tau_rela) )
!      IF(NRANK.EQ.0) write(6,'(A,I5,1P5E12.4)') &
!           'NR,RR,QLM,EPSRM2,tau_rela,theta_l', &
!           NR,RR,QLM(NR),EPSRM2(NR),tau_rela,theta_l

      phi = f_t/(1.D0 + (0.58D0+0.2D0*Z_i) &
                       *(2.D0*RR*QLM(NR)*EPSRM2(NR)**(-1.5D0) )*theta_l**2 &
                       /(3.D0*SQRT(2.D0*PI)*VC*tau_rela))
! AF modified 190313
!                       *(2.D0*RR*QLM(NR)*EPSRM2(NR)**(-1.5D0) ) &
!                       /(3.D0*SQRT(2.D0*PI)*VC*tau_rela)/theta_l**2 )
! /AF

      neoc=(1.D0-phi)*(1.D0-C_*phi)*(1.D0+0.47D0*(Z_i-1.D0))/ &
           (Z_i*(1.D0+0.27D0*(Z_i-1.D0)) )

!      neoc=( 1.D0-SQRT(EPSRM2(NR)) )**2 ! simple correction
      sigma=sigma*neoc
!      sigma=sigma*1.D0
!      WRITE(*,'(2I,6E14.6)') NRANK, NR, SIGMA, neoc, Z_i, phi, tau_rela, POST_LNLAM(NR,1,1)

      END SUBROUTINE SPITZER_SIGMA
! *******************************************************
      SUBROUTINE Tquench_trans

      USE fpmpi
      USE libmpi
      IMPLICIT NONE
      INTEGER:: NR, NSA, N, NSW
      REAL(rkind):: k, T0, Ts, Tf, r_tauq
      REAL(rkind),dimension(NRSTART:NREND):: RT1_temp
      REAL(rkind),save:: T_switch, time_switch

      IF(ISW_Q.eq.0)THEN ! linear 
         IF(TIMEFP+DELT.le.5.D0*tau_quench)THEN
            DO NR=NRSTART,NREND
!               k=( RT_quench_f(NR)-RTFP(NR,1) )/(5.D0*tau_quench)
               k=( RT_quench_f(NR)-RT_init(NR) )/(5.D0*tau_quench)
               RT1_temp(NR)=RT_init(NR) + k*(TIMEFP+DELT)
            END DO
         ELSE
            DO NR=NRSTART,NREND
               RT1_temp(NR)=RT_quench_f(NR)
            END DO
         END IF
      ELSEIF(ISW_Q.eq.1)THEN !
         IF(TIMEFP+DELT.le.tau_quench)THEN
            DO NR=NRSTART,NREND
               k=-LOG(RT_quench_f(NR)/RTFP(NR,1))
               RT1_temp(NR)=RTFP(NR,1)*EXP(-k*(TIMEFP+DELT)/tau_quench)
            END DO
         ELSE
            DO NR=NRSTART,NREND
               RT1_temp(NR)=RT_quench_f(NR)
            END DO
         END IF
      ELSEIF(ISW_Q.eq.2)THEN ! SMITH
         DO NR=NRSTART,NREND
            RT1_temp(NR)=RT_quench_f(NR)+ &
                 (RT_init(NR)-RT_quench_f(NR))*EXP(-(TIMEFP+DELT-time_quench_start)/tau_quench)
         END DO
      ELSEIF(ISW_Q.eq.3)THEN ! ITB break
         IF(TIMEFP.le.tau_quench*6.D0)THEN
            DO NR=NRSTART,NREND
               T0=1.D-1
               Ts=1.D-2
               Tf=(T0-Ts)*(1.D0-RM(NR)**2)+Ts
               RT1_temp(NR)=Tf + &
                    (RTFP(NR,1)-Tf)*EXP(-(TIMEFP+DELT-time_quench_start)/(tau_quench))
            END DO
            T_switch=RT1_temp(NRSTART)
            time_switch=TIMEFP+DELT
         ELSE
            DO NR=NRSTART,NREND
               RT1_temp(NR)=RT_quench_f(NR)+ &
                    (T_switch-RT_quench_f(NR))*EXP(-(TIMEFP+DELT-time_quench_start-time_switch)/tau_quench )
            END DO
         END IF
      ELSEIF(ISW_Q.eq.4)THEN ! ITB break and loss from edge
         IF(TIMEFP.le.tau_quench*1.D0)THEN
            DO NR=NRSTART,NREND
               T0=0.5D0*RTFP0(1)
               Ts=RTFPS(1)
               Tf=(T0-Ts)*(1.D0-RM(NR)**2)+Ts
               RT1_temp(NR)=Tf + &
                    (RTFP(NR,1)-Tf)*EXP(-(TIMEFP+DELT)/(tau_quench*0.5D0))
            END DO
            T_switch=RT1_temp(NRSTART)
            time_switch=TIMEFP+DELT
         ELSE 
            DO NR=NRSTART,NREND
               RT1_temp(NR)=RT_quench_f(NR)+ &
                    (T_switch-RT_quench_f(NR))*EXP(-(TIMEFP+DELT-time_switch)/(SQRT(1.D0-RM(NR))*tau_quench) )
            END DO
         END IF
      ELSEIF(ISW_Q.eq.5)THEN ! loss from edge
         DO NR=NRSTART,NREND
            r_tauq=tau_quench*(1.D0-0.9*RM(NR))
            RT1_temp(NR)=RT_quench_f(NR)+ &
                 (RTFP(NR,1)-RT_quench_f(NR))*EXP(-(TIMEFP+DELT)/r_tauq )
         END DO
      END IF

      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(RT1_temp,NREND-NRSTART+1,RT_quench)
      CALL mtx_reset_communicator

      END SUBROUTINE Tquench_trans
! *******************************************************
      SUBROUTINE update_rns_rjs(IP_all_FP)

      USE fpmpi
      USE libmpi
      IMPLICIT NONE
      integer:: NTH, NP, NR, NSA, N, NSW, NS, NSB
      integer:: nps, npe, nite
      REAL(rkind):: FLUXS_PMAX, FFP, RSUM1, FACT, SUMZ, RSUM2, PV
      REAL(rkind),intent(out):: IP_all_FP
      REAL(rkind):: RSUMP1, RSUMP2, zeff_imp2
      REAL(rkind),dimension(NRSTART:NREND):: temp_r
      
      CALL mtx_set_communicator(comm_np)
      DO NR=NRSTART, NREND
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)
            RSUM1=0.D0
            RSUM2=0.D0
            FLUXS_PMAX=0.D0

! for check
!            IF(NRSTART.eq.1.and.NPSTART.eq.1)THEN
!               open(17,file='p-n.dat')
!            END IF
            RSUMP2=0.D0
            DO NP=NPSTART, NPEND
               DO NTH=1, NTHMAX
                  RSUMP2 = RSUMP2 + VOLP(NTH,NP,NS)*FNSP(NTH,NP,NR,NSA)
               END DO
            END DO
            CALL p_theta_integration(RSUMP2)
            DO NITE=1,NPMAX
               RSUMP1=0.D0
               IF(NPSTART.ge.NITE)THEN
                  NPE=NPSTART
               ELSEIF(NPSTART.le.NITE.and.NPEND.ge.NITE)THEN
                  NPE=NITE
               ELSE
                  NPE=NPEND
               END IF

               DO NP=NPSTART, NPE
                  DO NTH=1, NTHMAX
                     RSUMP1 = RSUMP1 + VOLP(NTH,NP,NS)*FNSP(NTH,NP,NR,NSA)
                  END DO
               END DO
               CALL p_theta_integration(RSUMP1)
!               IF(NRSTART.eq.1.and.NPSTART.eq.1)THEN
!                  WRITE(17,'(I4,5E14.6)') NITE, PM(NITE,1), RSUMP2, RSUMP1, RSUMP1/RSUMP2, 1-RSUMP1/RSUMP2
!               END IF
            END DO
!            IF(NRSTART.eq.1.and.NPSTART.eq.1)THEN
!               close(17)
!            END IF
! end 

            IF(MODELA.eq.0)THEN
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NS)*FNSP(NTH,NP,NR,NSA)
                  END DO
               ENDDO
            ELSE
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NS)*FNSP(NTH,NP,NR,NSA) &
                          *RLAMDA(NTH,NR)
                  END DO
               ENDDO
            END IF

            IF(MODELA.eq.0) THEN
               IF(MODELR.EQ.0) THEN
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NS)*FNSP(NTH,NP,NR,NSA) &
                             *PM(NP,NS)*COSM(NTH)
                     END DO
                  ENDDO
               ELSE
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NS)**2)
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NS)*FNSP(NTH,NP,NR,NSA) &
                             *PM(NP,NS)*COSM(NTH)/PV
                     END DO
                  END DO
               ENDIF
            END IF
            
            CALL p_theta_integration(RSUM1)
            CALL p_theta_integration(RSUM2)
            FACT=RNFP0(NSA)*1.D20/RFSADG(NR)
            RNSL(NR,NSA) = RSUM1*FACT*1.D-20
            RJSL(NR,NSA) = RSUM2*FACT*AEFP(NSA)*PTFP0(NSA) &
                           /AMFP(NSA)*1.D-6
         END DO
      END DO

      CALL mtx_set_communicator(comm_nsanr) 
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL fp_gatherv_real8_sav(RNSL,SAVLEN(NRANK+1),RNS,N,NSA)
         CALL fp_gatherv_real8_sav(RJSL,SAVLEN(NRANK+1),RJS,N,NSA)
      END DO
      CALL mtx_reset_communicator
      CALL mtx_broadcast_real8(RNS,NRMAX*NSAMAX)
      CALL mtx_broadcast_real8(RJS,NRMAX*NSAMAX)

      IP_all_FP=0.D0
      DO NR=1,NRMAX
         IP_all_FP = IP_all_FP + RJS(NR,1)*VOLR(NR)
      END DO
      IP_ALL_FP = IP_ALL_FP/(2.D0*PI*RR)

      IF(MODEL_IMPURITY.eq.1.and.MODEL_connor_fp.eq.0)THEN
         DO NR=NRSTART, NREND
            temp_r(NR)=RN_MGI(NR,1)
         END DO

         CALL mtx_set_communicator(comm_nr)
         CALL mtx_allgather_real8(temp_r,NREND-NRSTART+1,RN_disrupt)
         CALL mtx_reset_communicator
      ELSE
         DO NR=1,NRMAX
            RN_disrupt(NR)=RNS(NR,1)
         END DO
      END IF

      SUMZ=0.D0
      IF(MODEL_IMPURITY.eq.0)THEN
         IF(NRANK.eq.0)THEN
            DO NSB=2,NSBMAX
               IF(NSB.le.NSAMAX)THEN
                  SUMZ=SUMZ+RNS(1,NSB)*PZ(NSB)**2
               ELSE
                  SUMZ=SUMZ+RNFD(1,NSB)*PZ(NSB)**2
               END IF
            END DO
            ZEFF = SUMZ/RNS(1,1)
         END IF
         CALL mtx_broadcast1_real8(ZEFF)
      ELSE
         IF(NRANK.eq.0)THEN
!         WRITE(*,*) "TEST ZEFF", zeff_imp2, target_zeff, spitot
            DO NSB=2,NSBMAX
               IF(NSB.ne.N_impu)THEN
                  SUMZ=SUMZ+RN_MGI(1,NSB)*PZ(NSB)**2
               ELSE
                  SUMZ=SUMZ+RN_MGI(1,NSB)*zeff_imp**2
               END IF
            END DO
            ZEFF = SUMZ / RN_MGI(1,1)
         END IF
         CALL mtx_broadcast1_real8(ZEFF)
      END IF

      END SUBROUTINE update_rns_rjs
! *******************************************************
      SUBROUTINE calculation_runaway_rate

      USE fpmpi
      USE libmpi
      IMPLICIT NONE
      integer:: NTH, NP, NR, NSA, N, NSW, NS, NSB
      REAL(rkind):: FLUXS_PMAX, FFP, RSUM1, FACT, SUMZ, RSUM2, PV, PITCH, v_thermal, tau_th
      REAL(rkind),dimension(NRSTART:NREND):: Rconnor_l
      REAL(rkind):: alp, z_i, h_alpha_z, lambda_alpha, gamma_alpha_z, theta_l, E00, tau_rela
      REAL(rkind),dimension(NTHMAX):: RE_PITCH_L
      
      CALL mtx_set_communicator(comm_np)
      DO NR=NRSTART, NREND
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)
            FLUXS_PMAX=0.D0
            RE_PITCH_L(:)=0.D0
            IF(MODEL_RE_pmax.eq.0)THEN
! FLUX S_p crossing pmax for runaway rate
               IF(NPEND.eq.NPMAX)THEN
                  DO NTH=1,NTHMAX
                     FFP=    PG(NPMAX+1,NS)*FNSP(NTH,NPMAX,NR,NSA)
                     
                     FLUXS_PMAX = FLUXS_PMAX +  &
                          FPP(NTH,NPMAX+1,NR,NSA)*FFP  &
                          * PG(NPMAX+1,NS)*SINM(NTH)!*tau_ta0(NSA) ! gamma
                     RE_PITCH_L(NTH)= &
                          FPP(NTH,NPMAX+1,NR,NSA)*FFP  & 
                          * PG(NPMAX+1,NS)*SINM(NTH)!*tau_ta0(NSA) ! gamma
                  END DO
               END IF
            ELSEIF(MODEL_RE_pmax.eq.1)THEN
! FLUX S_p across NPC_runaway for runaway rate
               IF(NPSTART.le.NPC_runaway.and.NPEND.ge.NPC_runaway)THEN
                  DO NTH=1,NTHMAX
                     FFP=    PG(NPC_runaway+1,NS)*FNSP(NTH,NPC_runaway,NR,NSA)
                     
                     FLUXS_PMAX = FLUXS_PMAX +  &
                          FPP(NTH,NPC_runaway+1,NR,NSA)*FFP  &
                          * PG(NPC_runaway+1,NS)*SINM(NTH)!*tau_ta0(NSA) ! gamma
                     RE_PITCH_L(NTH)= &
                          FPP(NTH,NPC_runaway+1,NR,NSA)*FFP  & 
                          * PG(NPC_runaway+1,NS)*SINM(NTH)!*tau_ta0(NSA) ! gamma
                  END DO
               END IF
            END IF
            
            CALL p_theta_integration(FLUXS_PMAX) ! same as allgather
            DO NTH=1,NTHMAX
               PITCH=RE_PITCH_L(NTH)
               CALL p_theta_integration(PITCH)
               RE_PITCH_L(NTH)=PITCH
            END DO
            
            IF(NSA.eq.1) THEN ! time integration
               FACT=RNFP0(NSA)*1.D20/RFSADG(NR)
               RSUM1=RNS(NR,1)/FACT*1.D20
               RFPL(NR)=2.D0*PI*DELTH* FLUXS_PMAX / RSUM1
               
               DO NTH=1,NTHMAX
                  RE_PITCH(NTH)=RE_PITCH(NTH)+ &
                       2.D0*PI*RE_PITCH_L(NTH)*FACT/(RNS(NR,1)*1.D20)&
                       *DELT*RN_disrupt(NR)
               END DO
            END IF
         END DO
      END DO
      CALL mtx_set_communicator(comm_nsa) 
      CALL mtx_broadcast_real8(RFPL,NREND-NRSTART+1)
      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(RFPL,NREND-NRSTART+1,RFP)
         
      IF(ISW_Z.eq.0)THEN
         z_i = Z_ef 
      ELSEIF(ISW_Z.eq.1)THEN
         z_i = ZEFF
      END IF
      DO NR=NRSTART,NREND
         IF(EP(NR).ge.ER_crit(NR))THEN
            E00 = EP(NR)/ER_drei(NR)
            IF(MODELR.eq.0)THEN
               v_thermal=SQRT( RT_quench(NR)*1.D3*AEE/AMFD(1)) 
               tau_th=(4.D0*PI*EPS0**2)*AMFP(1)**2*v_thermal**3/ &
                    ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RNFP0(1)*1.D20 )
               Rconnor_l(NR)=0.35D0*E00**(-3.D0*(Z_i+1.D0)/16.D0)&
                    *EXP(-0.25D0/E00 -SQRT( (1.D0+z_i)/E00 ) )!/tau_th
            ELSE
               theta_l=THETA0(1)*RT_quench(NR)/RTFP0(1)
               alp = E00/theta_l
               h_alpha_z=( alp*(z_i+1.D0) - z_i + 7.D0 + &
                    2.D0*SQRT(alp/(alp-1.D0))*(1.D0+z_i)*(alp-2.D0) ) &
                    /( 16.D0*(alp-1.D0) )
               lambda_alpha=8.D0*alp*(alp-0.5D0-SQRT(alp*(alp-1.D0)) )
               gamma_alpha_z=SQRT( (1.0+z_i)*alp**2/(8.D0*(alp-1.D0)) )&
                    *(0.5D0*PI-ASIN(1.D0-2.D0/alp))
               Rconnor_l(NR)= E00**(-h_alpha_z) &
                    *EXP(-0.25D0*lambda_alpha/E00-SQRT(2.D0/E00)*gamma_alpha_z )
               
               IF(MODEL_IMPURITY.eq.0)THEN
!                  tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
!                       ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RNFP0(1)*1.D20 )
                  tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
                       ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RNFP(NR,1)*1.D20 )
               ELSE
!                  tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
!                       ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RN0_MGI(1)*1.D20 )
                  tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
                       ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RN_MGI(NR,1)*1.D20 )
               END IF
               Rconnor_l(NR)=Rconnor_l(NR)/(tau_rela*SQRT(2*theta_l)**3)
!               Rconnor_l(NR)=Rconnor_l(NR)/(tau_rela*SQRT(theta_l)**3)
            END IF
         ELSE
            Rconnor_l(NR)=0.D0
         END IF
      END DO
      CALL mtx_set_communicator(comm_nr) 
      call mtx_allgather_real8(Rconnor_l,NREND-NRSTART+1,Rconnor)
      CALL mtx_reset_communicator


      END SUBROUTINE calculation_runaway_rate
! *******************************************************
      SUBROUTINE update_disruption_quantities(IP_all,IP_ohm,IP_run,IP_prim,IP_BS,l_ind)

      USE fpmpi
      USE libmpi
      USE plprof
      IMPLICIT NONE
      INTEGER:: NR, NSA
      REAL(rkind),intent(out):: IP_all, IP_ohm, IP_run, IP_BS, l_ind, IP_prim
      REAL(rkind),dimension(NRSTART:NREND):: &
           RN1_temp, RN2_temp, RJ1_temp, QLM_L, RN2P_temp
      REAL(rkind),dimension(NRMAX):: SUMIP4
      REAL(rkind):: SUMIP, SUMIP_ohm, SUMIP_run, SUMIP_BS, rhon, BP_a, SUMIP_prim, sumip_integ
      REAL(rkind),save:: integ

      IF(TIMEFP.eq.DELT) integ=0.D0
      IF(TIMEFP.eq.DELT) previous_rate(:)=0.D0
      IF(TIMEFP.eq.DELT) previous_rate_p(:)=0.D0

      DO NR=NRSTART,NREND
         IF(MODEL_Connor_FP.eq.0)THEN
            RN1_temp(NR)=-DELT*(Rconnor(NR)+RFP_ava(NR))*RN_disrupt(NR)+RN_disrupt(NR)
         ELSEIF(MODEL_Connor_FP.eq.1)THEN
            RN1_temp(NR)=-0.5D0*DELT* &
                 ( (RFP(NR)+RFP_ava(NR))*RN_disrupt(NR) + previous_rate(NR) ) &
                 +RN_disrupt(NR)
            previous_rate(NR)=(RFP(NR)+RFP_ava(NR))*RN_disrupt(NR)
            previous_rate_p(NR)=(RFP(NR))*RN_disrupt(NR)
         END IF
         RJ1_temp(NR)=SIGMA_SPP(NR)*EP(NR)*1.D-6
      END DO
      
      IF(MODELD_n_RE.eq.1.and.RN_runaway(1).ge.1.D-40)THEN
         CALL N_RE_transport ! RN_runaway, RN_drei updated
      ELSE
         DO NR=NRSTART,NREND
            IF(MODEL_Connor_FP.eq.0)THEN
               RN2_temp(NR)=DELT*(Rconnor(NR)+RFP_ava(NR))*RN_disrupt(NR)+RN_runaway(NR)
               RN2P_temp(NR)=DELT*(Rconnor(NR))*RN_disrupt(NR)+RN_drei(NR)
            ELSEIF(MODEL_Connor_FP.eq.1)THEN
               RN2_temp(NR)=0.5D0*DELT* &
                    ( (RFP(NR)+RFP_ava(NR))*RN_disrupt(NR) + previous_rate(NR) ) &
                    +RN_runaway(NR)
               RN2P_temp(NR)=0.5D0*DELT* &
                    ( (RFP(NR))*RN_disrupt(NR) + previous_rate_p(NR) ) &
                    +RN_drei(NR)
            END IF
         END DO
         CALL mtx_set_communicator(comm_nr)
         call mtx_allgather_real8(RN2_temp,NREND-NRSTART+1,RN_runaway)
         call mtx_allgather_real8(RN2P_temp,NREND-NRSTART+1,RN_drei)
      END IF

      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(RN1_temp,NREND-NRSTART+1,RN_disrupt)
      call mtx_allgather_real8(RJ1_temp,NREND-NRSTART+1,Rj_ohm)
      SUMIP=0.D0
      SUMIP_ohm=0.D0
      SUMIP_run=0.D0
      SUMIP_prim=0.D0
      SUMIP_BS=0.D0
      SUMIP_integ=0.D0
      SUMIP4(:)=0.D0
      DO NR=1,NRMAX
         RJ_runaway(NR)=AEE*v_RE*VC*RN_runaway(NR)*1.D20*1.D-6
         SUMIP = (RJ_ohm(NR)+RJ_runaway(NR)+RJ_bs(NR)) &
              *VOLR(NR)/(2.D0*PI*RR) + SUMIP
         SUMIP_ohm = (RJ_ohm(NR)) &
              *VOLR(NR)/(2.D0*PI*RR) + SUMIP_ohm
         SUMIP_run =  (RJ_runaway(NR)) &
              *VOLR(NR)/(2.D0*PI*RR) + SUMIP_run
         SUMIP_BS =  (RJ_bs(NR)) &
              *VOLR(NR)/(2.D0*PI*RR) + SUMIP_BS
         SUMIP_prim = AEE*v_RE*VC*RN_drei(NR)*1.D20*1.D-6 &
              *VOLR(NR)/(2.D0*PI*RR) + SUMIP_prim
         SUMIP_integ=(RJ_ohm(NR) +RJ_runaway(NR) +RJ_bs(NR))*VOLR(NR)/(2.D0*PI*RR) + SUMIP_integ
         SUMIP4(NR)=SUMIP_integ
      END DO
      IP_all=SUMIP
      IP_ohm=SUMIP_ohm
      IP_run=SUMIP_run
      IP_prim=SUMIP_prim
      IP_BS=SUMIP_BS
!
      SUMIP=0.D0
      l_ind=0.D0
      DO NR=1,NRMAX
         SUMIP=SUMIP+SUMIP4(NR)
         QLM(NR)=2*PI*RM(NR)**2*RA**2*BB/(RMU0*SUMIP*1.D6*RR)
         rhon=RM(NR)
         BP_a = (RSRHON(RHON)*BB/(RR*QLM(NR)) )**2
         l_ind = l_ind + RM(NR)*DELR*BP_a
      END DO
      l_ind = 2*l_ind/BP_a
      CALL mtx_reset_communicator

      END SUBROUTINE update_disruption_quantities
! *******************************************************
      SUBROUTINE AVALANCHE

      USE fpmpi
      USE libmpi
      IMPLICIT NONE
      integer:: NR
      REAL(rkind),dimension(NRSTART:NREND):: rate_ava_l
      REAL(rkind):: E_hat, tau_rela, phi, theta_l, Z_i

      IF(ISW_Z.eq.0)THEN
         Z_i=Z_ef
      ELSEIF(ISW_Z.eq.1)THEN
         Z_i=ZEFF
      END IF
      DO NR=NRSTART, NREND
         IF(MODEL_IMPURITY.eq.0)THEN
!            tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
!                 ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RNFP0(1)*1.D20 )
            tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
                 ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RNFP(NR,1)*1.D20 )
         ELSE
!            tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
!                 ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RN0_MGI(1)*1.D20 )
            tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
                 ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RN_MGI(NR,1)*1.D20 )
         END IF
!         theta_l=THETA0(1)*RT_quench(NR)/RTFP0(1)
         E_hat=EP(NR)/ER_crit(NR)
         phi=1.D0-1.46D0*SQRT(EPSRM2(NR))+1.72D0*EPSRM2(NR)
         IF(E_hat.ge.1)THEN
            IF(MODEL_Connor_FP.eq.0)THEN
               rate_ava_l(NR)= &
                    RN_runaway(NR)/RN_disrupt(NR)* &
                    (E_hat-1.D0)/ &
                    (tau_rela*POST_LNLAM(NR,1,1))* &
                    SQRT(PI*phi/(3.D0*(Z_i+5.D0)))/ &
                    SQRT(1.D0-1.D0/E_hat+ (4.D0*PI*(Z_i+1.D0)**2)/ &
                    (3.D0*phi*(Z_i+5.D0)*(E_hat**2+4.D0/phi**2-1.D0)) )
            ELSEIF(MODEL_Connor_FP.eq.1)THEN
               rate_ava_l(NR)= &
                    RN_runaway(NR)/RN_disrupt(NR)* &
                    (E_hat-1.D0)/ &
                    (tau_rela*POST_LNLAM(NR,1,1))* &
                    SQRT(PI*phi/(3.D0*(Z_i+5.D0)))/ &
                    SQRT(1.D0-1.D0/E_hat+ (4.D0*PI*(Z_i+1.D0)**2)/ &
                    (3.D0*phi*(Z_i+5.D0)*(E_hat**2+4.D0/phi**2-1.D0)) )
            END IF
         ELSE
            rate_ava_l(NR)=0.D0
         END IF
      END DO

!      rate_ava_l(NRSTART)=0.D0 ! for given E
      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(rate_ava_l,NREND-NRSTART+1,RFP_ava)
!
      CALL mtx_reset_communicator

      END SUBROUTINE AVALANCHE
! *******************************************************
      SUBROUTINE update_fpp

      IMPLICIT NONE
      integer:: nr, nsa, nth, np

      DO NSA=NSASTART, NSAEND
         DO NR= NRSTART, NREND
            DO NP=NPSTART, NPENDWG
               DO NTH=1, NTHMAX
                  FPP(NTH,NP,NR,NSA)= &
                       FEPP(NTH,NP,NR,NSA)+FCPP(NTH,NP,NR,NSA)
               END DO
            END DO

            DO NP=NPSTARTW, NPENDWM
               DO NTH=1, NTHMAX+1
                  FTH(NTH,NP,NR,NSA)= &
                       FETH(NTH,NP,NR,NSA)+FCTH(NTH,NP,NR,NSA)
               END DO
            END DO

            IF(NPENDWG.eq.NPMAX+1)THEN 
               DO NTH=1,NTHMAX
                  FPP(NTH,NPMAX+1,NR,NSA)=max(0.D0,FPP(NTH,NPMAX+1,NR,NSA))
               END DO
            END IF
         END DO
      END DO

      END SUBROUTINE update_fpp
!**********************************************
      SUBROUTINE bootstrap_current(NR,jbs)

      USE plprof
      IMPLICIT NONE
      INTEGER,intent(in):: NR
      REAL(rkind),intent(out):: jbs
      REAL(rkind):: BP, RNE, RNI, RTE, RTI, dndr, dtedr, dtidr, rhon

      RTE=RT_quench(NR)*1.D3*AEE
      RTI=RT_quench(NR)*1.D3*AEE
      RNE=RN_disrupt(NR)*1.D20
      IF(NRMAX.ne.1)THEN
         IF(NR.eq.1)THEN
            dndr=(RN_disrupt(NR+1)-RN_disrupt(NR))/(DELR*RA)*1.D20
            dtedr=(RT_quench(NR+1)-RT_quench(NR))/(DELR*RA)*1.D3*AEE
            dtidr=dtedr
         ELSEIF(NR.eq.NRMAX)THEN
            dndr=(RN_disrupt(NRMAX)-RN_disrupt(NR-1))/(DELR*RA)*1.D20
            dtedr=(RT_quench(NRMAX)-RT_quench(NR-1))/(DELR*RA)*1.D3*AEE
            dtidr=dtedr
         ELSE
            dndr=(RN_disrupt(NR+1)-RN_disrupt(NR-1))/(2.D0*DELR*RA)*1.D20
            dtedr=(RT_quench(NR+1)-RT_quench(NR-1))/(2.D0*DELR*RA)*1.D3*AEE
            dtidr=dtedr            
         END IF
      ELSE
         dndr=-PROFN1(1)*R1**(PROFN1(1)-1.D0) &
              *(1-R1**PROFN1(1))**(PROFN2(1)-1.D0)*PROFN2(1) &
              *(PN(1)-PNS(1))*1.D20/RA
         dtedr=-PROFT1(1)*R1**(PROFT1(1)-1.D0) &
              *(1-R1**PROFT1(1))**(PROFT2(1)-1.D0)*PROFT2(1) &
              *(PTPR(1)-PTS(1))*1.D3*AEE/RA
         dtidr=dtedr
      END IF

      RHON=RM(NR)
      BP= RSRHON(RHON)*BB/(RR*QLM(NR))

!     Wesson P. 173
      jbs=-SQRT(EPSRM2(NR))*RNE/BP* &
           (2.44D0*(RTE+RTI)*dndr/RNE + &
           0.69D0*dtedr - 0.42D0*dtidr) *1.D-6! *1.D-10

!      IF(NRANK.eq.0.and.nr.eq.1)THEN
!         WRITE(6,'(A,7E14.6)') "TEST_BS=",jbs, bp, rte, rti, rne, dndr, dtedr
!      END IF

      END SUBROUTINE bootstrap_current
!**********************************************
      SUBROUTINE djdt

      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      INTEGER:: NSA, NR, NP, NTH, NPS, NS
      REAL(rkind):: PV, WPL, WPM, WPP, DFP, DFT, FFP
      REAL(rkind):: P2_S_P, P2_S_T, RSUMP, RSUMT, FACT
      REAL(rkind),DIMENSION(NRSTART:NREND):: R_djdtl

      CALL mtx_set_communicator(comm_np)
      DO NSA=NSASTART, NSAEND
         NS=NS_NSA(NSA)
         DO NR=NRSTART, NREND
            RSUMP=0.D0
            RSUMT=0.D0
            IF(NPSTART.eq.1)THEN
               NPS=2
            ELSE
               NPS=NPSTART
            END IF
            DO NP=NPS,NPEND
               PV=SQRT(1.D0+THETA0(NSA)*PG(NP,NS)**2)
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
                  DFP=    PG(NP,NS) &
                       /DELP(NS)*(FNSP(NTH,NP,NR,NSA)-FNSP(NTH,NP-1,NR,NSA))
                  IF(NTH.EQ.1) THEN
                     DFT=1.D0/DELTH                             &
                         *(                                     &
                            ((1.D0-WPP)*FNSP(NTH+1,NP  ,NR,NSA)   &
                                  +WPP *FNSP(NTH+1,NP-1,NR,NSA))&
                           -                                    &
                            ((1.D0-WPM)*FNSP(NTH,NP  ,NR,NSA)     &
                                  +WPM *FNSP(NTH,NP-1,NR,NSA))&
                          )

                  ELSE IF(NTH.EQ.NTHMAX) THEN
                     DFT=    1.D0/DELTH                         & 
                         *(-                                    &
                            ((1.D0-WPM)*FNSP(NTH-1,NP  ,NR,NSA)   &
                                  +WPM *FNSP(NTH-1,NP-1,NR,NSA))&
                          +                                     &
                            ((1.D0-WPP)*FNSP(NTH,NP  ,NR,NSA)     &
                                  +WPP *FNSP(NTH,NP-1,NR,NSA))&
                          )
                  ELSE
                     DFT=    1.D0/(2.D0*DELTH)                  &
                         *(                                     &
                            ((1.D0-WPP)*FNSP(NTH+1,NP  ,NR,NSA)   &
                                  +WPP *FNSP(NTH+1,NP-1,NR,NSA))&
                           -                                    &
                            ((1.D0-WPM)*FNSP(NTH-1,NP  ,NR,NSA)   &
                                  +WPM *FNSP(NTH-1,NP-1,NR,NSA))&
                                  )
                  ENDIF
                  FFP=    PG(NP,NS)                           &
                         *((1.D0-WPL)*FNSP(NTH  ,NP  ,NR,NSA)  &
                                +WPL *FNSP(NTH  ,NP-1,NR,NSA))

                  P2_S_P = PG(NP,NS)*( &
                       -DPP(NTH,NP,NR,NSA)*DFP &
                       -DPT(NTH,NP,NR,NSA)*DFT &
                       +FPP(NTH,NP,NR,NSA)*FFP )
                  
                  RSUMP = RSUMP + &
                       SINM(NTH)*( &
                       P2_S_P*COSM(NTH)/PV**3 )
               END DO
            END DO

            DO NP=NPSTART,NPEND
               PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NS)**2)
               DO NTH=2,NTHMAX
                  WPL=WEIGHT(NTH  ,NP,NR,NSA)
                  IF(NP.EQ.1) THEN
                     WPM=0.D0
                  ELSE
                     WPM=WEIGHT(NTH,NP-1,NR,NSA)
                  ENDIF
                  IF(NP.EQ.NPMAX) THEN
                     WPP=0.D0
                  ELSE
                     WPP=WEIGHT(NTH,NP+1,NR,NSA)
                  ENDIF
!                     DFP = PM(NP,NS)*0.5D0/DELP*( &
!                          (1.D0-WPP)*FNSP(NTH+1,NP+1,NR,NSA) &
!                          +     WPP *FNSP(NTH,  NP+1,NR,NSA) - &
!                          (1.D0-WPM)*FNSP(NTH+1,NP-1,NR,NSA) &
!                          -     WPM *FNSP(NTH,  NP-1,NR,NSA) &
!                          )
                  DFT = 1.D0/DELTH &
                       *(FNSP(NTH,NP,NR,NSA)-FNSP(NTH-1,NP,NR,NSA))
                  FFP = PM(NP,NS)*( &
                       (1.D0-WPL)*FNSP(NTH,  NP,NR,NSA) &
                       +     WPL *FNSP(NTH-1,NP,NR,NSA) )

                  P2_S_T = PM(NP,NS)*( &
!                       -DTP(NTH,NP,NR,NSA)*DFP &
                       -DTT(NTH,NP,NR,NSA)*DFT &
                       +FTH(NTH,NP,NR,NSA)*FFP )

                  RSUMT = RSUMT + &
                       SING(NTH)**2*P2_S_T/PV
               END DO
            END DO
            CALL p_theta_integration(RSUMP)
            CALL p_theta_integration(RSUMT)

            FACT=RNFP0(NSA)*1.D20/RFSADG(NR)
            R_djdtl(NR) = (RSUMP-RSUMT)*FACT*AEFP(NSA)*PTFP0(NSA) &
                 /AMFP(NSA)*1.D-6
         END DO
      END DO

      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(R_djdtl,NREND-NRSTART+1,R_djdt)
      CALL mtx_reset_communicator

      IF(NRANK.eq.0)THEN
         WRITE(6,'(A, 2E14.6)') "DELJ = ", R_djdt(1)*DELT, (RJS(1,1)-RJS_M(1,1))
      END IF

      END SUBROUTINE djdt
!**********************************************
      SUBROUTINE N_RE_transport

      USE libmpi
      USE libmtx
      IMPLICIT NONE
      INTEGER:: NR
      integer:: imtxstart1,imtxend1,its
      REAL(rkind):: Dr_coef, factp, factr
      REAL(rkind):: Aij, Aij_m1, Aij_p1, RHS, b_wall, E_e, a_e, coef_ln

      CALL mtx_set_communicator(comm_nr) 
      FACTP=PI*RR*VC/SQRT(2.D0)

!=RN_runaway
      CALL mtx_setup(NRMAX,imtxstart1,imtxend1)
!---- DIAGONAL TERM
      DO NR=NRSTART,NREND
         FACTR=QLM(NR)*deltaB_B**2
         DR_coef=FACTP*FACTR/(RA**2)
         Aij = 1.D0 + 2.D0*DR_coef*DELT/DELR**2
         CALL mtx_set_matrix(NR,NR,Aij)
         CALL mtx_set_vector(NR,RN_runaway(NR))
      END DO
!---- OFF DIAGONAL
      DO NR=NRSTART,NREND
         FACTR=QLM(NR)*deltaB_B**2
         DR_coef=FACTP*FACTR/(RA**2)
         Aij_m1=-DR_coef*DELT/(DELR**2)*(1.D0-0.5D0*DELR/RM(NR))
         Aij_p1=-DR_coef*DELT/(DELR**2)*(1.D0+0.5D0*DELR/RM(NR))
         IF(NR.ne.1)THEN
            CALL mtx_set_matrix(NR,NR-1,Aij_m1)
         END IF
         IF(NR.ne.NRMAX)THEN
            CALL mtx_set_matrix(NR,NR+1,Aij_p1)
         END IF
      END DO
!---- RIGHT HAND SIDE
      DO NR=NRSTART,NREND
         IF(MODEL_Connor_FP.eq.0)THEN
            RHS = RN_runaway(NR) + (Rconnor(NR)+RFP_ava(NR))*RN_disrupt(NR)*DELT
         ELSE
            RHS = RN_runaway(NR) + (RFP(NR)+RFP_ava(NR))*RN_disrupt(NR)*DELT
         END IF
         CALL mtx_set_source(NR,RHS)
      END DO
!---- SOLVE
      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) 
      CALL mtx_gather_vector(RN_runaway)
      CALL mtx_cleanup
      CALL mtx_reset_communicator

      CALL mtx_set_communicator(comm_nr) 
!=RN_drei
      CALL mtx_setup(NRMAX,imtxstart1,imtxend1)
!---- DIAGONAL TERM
      DO NR=NRSTART,NREND
         FACTR=QLM(NR)*deltaB_B**2
         DR_coef=FACTP*FACTR/(RA**2)
         Aij = 1.D0 + 2.D0*DR_coef*DELT/(DELR)**2
         CALL mtx_set_matrix(NR,NR,Aij)
         CALL mtx_set_vector(NR,RN_drei(NR))
      END DO
!---- OFF DIAGONAL
      DO NR=NRSTART,NREND
         FACTR=QLM(NR)*deltaB_B**2
         DR_coef=FACTP*FACTR/(RA**2)
         Aij_m1=-DR_coef*DELT/(DELR**2)*(1.D0-0.5D0*DELR/RM(NR))
         Aij_p1=-DR_coef*DELT/(DELR**2)*(1.D0+0.5D0*DELR/RM(NR))
         IF(NR.ne.1)THEN
            CALL mtx_set_matrix(NR,NR-1,Aij_m1)
         END IF
         IF(NR.ne.NRMAX)THEN
            CALL mtx_set_matrix(NR,NR+1,Aij_p1)
         END IF
      END DO
!---- RIGHT HAND SIDE
      DO NR=NRSTART,NREND
         IF(MODEL_Connor_FP.eq.0)THEN
            RHS = RN_drei(NR) + (Rconnor(NR))*RN_disrupt(NR)*DELT
         ELSE
            RHS = RN_drei(NR) + (RFP(NR))*RN_disrupt(NR)*DELT
         END IF
         CALL mtx_set_source(NR,RHS)
      END DO
!---- SOLVE
      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) 
      CALL mtx_gather_vector(RN_drei)
      CALL mtx_cleanup
      CALL mtx_reset_communicator

      END SUBROUTINE N_RE_transport
!**********************************************
!      SUBROUTINE N_RE_transport_EXPLICIT
!       
!      IMPLICIT NONE
!      integer:: NR
!      REAL(rkind):: factr, factp, dr_coef
!
!      DO NR=1,NRMAX
!         FACTP=PI*RR*VC/SQRT(2.D0)
!         FACTR=QLM(NR)*deltaB_B**2
!         DR_coef=FACTP*FACTR/(RA**2)
!
!         RN_runaway(NR) + DR_coef*DELT/(RM(NR)*DELR**2)*( &
!              RG(NR+1)*RN_runaway(NR+1) &
!              -2.D0*RM(NR)*RN_runaway(NR) &
!              +RG(NR)*RN_runaway(NR-1) ) &
!              +(Rconnor(NR)+RFP_ava(NR))*RN_disrupt(NR)
!
!      END DO
!
!      END SUBROUTINE N_RE_transport_EXPLICIT
!**********************************************
      SUBROUTINE FLUXS_PTH

      USE libmpi
      IMPLICIT NONE
      INTEGER:: NTH, NP, NR, NSA, NPS, NPE, NS, NTHSTEP
      REAL(rkind):: WPL, WPM, WPP, WTP, DFT, DFP, FFP, FFT
      REAL(rkind), dimension(NTHMAX,NPSTART:NPEND+1):: FLUXS_PL_L
      REAL(rkind), dimension(NTHMAX+1,NPSTART:NPEND):: FLUXS_TL_L
      REAL(rkind), dimension(NTHMAX,NPSTART:NPEND):: FLUXS_PL, FLUXS_TL!, dfdt_l
      REAL(rkind), dimension(NTHMAX,NPMAX):: FLUXS_P, FLUXS_T!, dfdt

      IF(NPSTART.eq.1)THEN
         NPS=2
      ELSE
         NPS=NPSTART
      END IF
      IF(NPEND.eq.NPMAX)THEN
         NPE = NPEND
      ELSE
         NPE = NPEND+1
      END IF

      NSA=1
      FLUXS_PL_L(:,:)=0.D0
      NR=NRSTART
      NS=NS_NSA(NSA)
!FLUXS_P
      DO NP=NPS,NPE
         DO NTH=1,NTHMAX
            WPL=WEIGHP(NTH  ,NP,NR,NSA)
            DFP=    1.D0 &
                 /DELP(NS)*(FNSP(NTH,NP,NR,NSA)-FNSP(NTH,NP-1,NR,NSA))
            FFP=    1.D0                           &
                 *((1.D0-WPL)*FNSP(NTH  ,NP  ,NR,NSA)  &
                 +WPL *FNSP(NTH  ,NP-1,NR,NSA))
            FLUXS_PL_L(NTH,NP)= &
                 -DPP(NTH,NP,NR,NSA)*DFP+FPP(NTH,NP,NR,NSA)*FFP
         END DO
      END DO
      
      DO NP=NPSTART,NPEND
         DO NTH=1,NTHMAX
            IF(FNSP(NTH,NP,NR,NSA).le.1.D-43)THEN
               FLUXS_PL(NTH,NP)= 0.D0
            ELSE
               FLUXS_PL(NTH,NP)= &
                    0.5D0*( FLUXS_PL_L(NTH,NP)+FLUXS_PL_L(NTH,NP+1) )
            END IF
         END DO
      END DO
!FLUXS_T
      DO NP=NPSTART,NPEND
         DO NTH=1,NTHMAX+1
            WTP = WEIGHT(NTH,NP,NR,NSA)
            IF(NTH.eq.1)THEN
               FFT= FNSP(1,NP,NR,NSA) 
               DFT = 0.D0
            ELSEIF(NTH.eq.NTHMAX+1)THEN
               FFT= FNSP(NTHMAX,NP,NR,NSA) 
               DFT = 0.D0
            ELSE
               FFT= (1.D0-WTP)*FNSP(NTH,NP,NR,NSA)+ &
                    WTP * FNSP(NTH,NP,NR,NSA) 
               DFT = 1.D0/DELTH*&
                    (FNSP(NTH,NP,NR,NSA)-FNSP(NTH-1,NP,NR,NSA))
            END IF
            FLUXS_TL_L(NTH,NP) = -DFT*DCTT(NTH,NP,NR,NSA)/PM(NP,NSA) &
                 + FFT*FTH(NTH,NP,NR,NSA)
         END DO
      END DO
      DO NP=NPSTART,NPEND
         DO NTH=1,NTHMAX
            IF(FNSP(NTH,NP,NR,NSA).le.1.D-43)THEN ! p~15
               FLUXS_TL(NTH,NP)= 0.D0
            ELSE
               FLUXS_TL(NTH,NP)= &
                 0.5D0*( FLUXS_TL_L(NTH+1,NP)+FLUXS_TL_L(NTH,NP) )
            END IF
         END DO
      END DO
! dfdt_l
!      DO NP=NPSTART,NPEND
!         DO NTH=1,NTHMAX
!            dfdt_l(NTH,NP)=(FNSP(NTH,NP,1,1)-FNSM(NTH,NP,1,1) )/DELT
!         END DO
!      END DO

      CALL mtx_set_communicator(comm_np)
      call mtx_allgather_real8(FLUXS_PL,(NPEND-NPSTART+1)*NTHMAX,FLUXS_P)
      call mtx_allgather_real8(FLUXS_TL,(NPEND-NPSTART+1)*NTHMAX,FLUXS_T)
!      call mtx_allgather_real8(DFDT_L,(NPEND-NPSTART+1)*NTHMAX,DFDT)
      CALL mtx_reset_communicator 


      IF(NSASTART.eq.1.and.NRSTART.eq.1.and.NRANK.eq.0)THEN      

      WRITE(17,'(A, E14.6)') "# X, Y, Z, vx, vy, vz, time=", TIMEFP
      DO NP=2, NPMAX, N_partition_p
         IF(PM(NP,NSA).le.2.D0)THEN
            NTHSTEP=8*NTHMAX/64
         ELSEIF(PM(NP,NSA).le.5.D0)THEN
            NTHSTEP=4*NTHMAX/64
         ELSE
            NTHSTEP=2*NTHMAX/64
         END IF
         DO NTH=1,NTHMAX,NTHSTEP
            WRITE(17,'(20E16.8)') &
                 PM(NP,NSA)*COSM(NTH), PM(NP,NSA)*SINM(NTH), PM(NP,NSA)*0.D0, &
                 FLUXS_P(NTH,NP)*COSM(NTH)-FLUXS_T(NTH,NP)*SINM(NTH), &
                 FLUXS_P(NTH,NP)*SINM(NTH)+FLUXS_T(NTH,NP)*COSM(NTH), &
                 FLUXS_P(NTH,NP)*0.D0, FLUXS_P(NTH,NP), FLUXS_T(NTH,NP), FNS(NTH,NP,1,1)!, &
!                 DFDT(NTH,NP)*PM(NP,NSA)*COSM(NTH), DFDT(NTH,NP)*PM(NP,NSA)*SINM(NTH)
         END DO
      END DO
      WRITE(17,*) " "
      WRITE(17,*) " "
      END IF

      END SUBROUTINE FLUXS_PTH
!**********************************************
      SUBROUTINE TOP_OF_TIME_LOOP_DISRUPT(NT)

      USE libmpi
      USE libmtx
      USE fpcoef
      use fpcaleind
      IMPLICIT NONE
      INTEGER:: NTH, NP, NR, NSA
      INTEGER,intent(in):: NT
      REAL(rkind):: sigma

      DO NR=NRSTART,NREND
         EM(NR)=E1(NR)
         EP(NR)=E1(NR)
         SIGMA_SPM(NR)=conduct_sp(NR)
      END DO
      DO NR=NRSTART-1, NREND+1 ! for crank-nicolson
         IF(NR.eq.0)THEN
            EM_W(NR)=0.D0
         ELSEIF(NR.ne.NRMAX+1)THEN
            EM_W(NR)=E1(NR)
         END IF
      END DO

      IF(TIMEFP.ge.time_quench_start) CALL Tquench_trans ! evolution of temperature
      IF(MODEL_IMPURITY.eq.1.and.TIMEFP.le.tau_mgi)THEN ! evolution of electron density
         CALL MGI_DENSITY 
      END IF

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
            IF (MOD(NT,NTSTEP_COEF).EQ.0) CALL FP_COEF(NT)
!            END IF
         DO NR=1,NRMAX
            RN_runaway_M(NR)=RN_runaway(NR)
         END DO
!         IF(TIMEFP.ge.time_quench_start) CALL READ_E_FIELD  !given

      END SUBROUTINE TOP_OF_TIME_LOOP_DISRUPT
!**********************************************
      SUBROUTINE E_FIELD_EVOLUTION_DISRUPT(NT,IP_all_FP,DEPS_E2)

      USE libmpi
      USE libmtx
      USE fpcoef
      use fpcaleind
      IMPLICIT NONE
      INTEGER,intent(in):: NT
      INTEGER:: NR, NSA
      INTEGER:: NDIMPL, ILOC1
      REAL(rkind):: DEPS_E, DEPS_MAX
      REAL(rkind),intent(out):: IP_all_FP, DEPS_E2
      REAL(rkind),dimension(NRSTART:NREND):: E_SIGMA 

      NDIMPL=0
      DEPS_E=0.D0
      IF(MODEL_jfp.eq.0)THEN
         call update_rns_rjs(IP_all_FP)
         DO while(NDIMPL.le.2)
            NDIMPL=NDIMPL+1
            call calculation_runaway_rate
            CALL AVALANCHE
            IF(TIMEFP.ge.time_quench_start) CALL E_IND_IMPLICIT !non-given
            DEPS_E=(E1(NRSTART)-EP(NRSTART))**2/EM(NRSTART)**2
            
            DO NR=NRSTART,NREND
               EP(NR)=E1(NR)
            END DO
            CALL FP_CALE
            CALL update_fpp
         END DO
         IF (MOD(NT,NTG1STEP).EQ.0.and.NRANK.eq.0) &
              WRITE(6,'(A,E14.6)') "CALE_CONVERSION ON AXIS= ", DEPS_E
      END IF
      DEPS_E=0.D0
      DEPS_E2=0.D0

!
      IF(MODEL_jfp.eq.1)THEN
         call update_rns_rjs(IP_all_FP)
         DO NR=NRSTART, NREND
            E_SIGMA(NR)=EP(NR)
         END DO
         DO while(NDIMPL.le.LMAXFP)
            NDIMPL=NDIMPL+1
            call calculation_runaway_rate
            CALL AVALANCHE
            IF(TIMEFP.ge.time_quench_start) CALL E_IND_IMPLICIT_FP(E_SIGMA)
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
            CALL FP_CALE
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


      END SUBROUTINE E_FIELD_EVOLUTION_DISRUPT
!**********************************************
      SUBROUTINE FILE_OUTPUT_DISRUPT(NT,IP_all_FP)

      IMPLICIT NONE
      INTEGER,intent(in):: NT
      INTEGER:: NR, NSA, NP, NTH

      REAL(rkind):: ip_all, ip_ohm, ip_run, jbs, IP_bs
      REAL(rkind):: IP_all_FP, l_ind, IP_prim, Z_i

      CALL update_disruption_quantities(IP_all, IP_ohm, IP_run, IP_prim, IP_bs, l_ind)

      IF(NRANK.eq.0)THEN

      IF(ISW_Z.eq.0)THEN
         Z_i=Z_ef
      ELSEIF(ISW_Z.eq.1)THEN
         Z_i=ZEFF 
      END IF

      IF(MOD(NT,NTG1STEP).EQ.0) THEN
         WRITE(6,'(7A14)') "IP_all","IP_ohm", "IP_run", "IP_primary", "IP_bs", "RE gen.rate"
         WRITE(6,'(1P7E14.6)') IP_all, IP_ohm, IP_run, IP_prim, IP_bs, RFP(1)*tau_ta0(1)
         WRITE(6,'(A,1PE14.6,A,1PE14.6)') "Zeff= ", ZEFF, " l_ind= ", l_ind
         IF(MODEL_IMPURITY.eq.1)THEN
            WRITE(6,'(A,1P3E14.6)') "RN_MGI = ", RN_MGI(1,1), RN_MGI(1,2), RN_MGI(1,3)
         END IF
! radial.dat
         WRITE(13,'(A, 1PE15.6e3, i7)') "# TIME ", TIMEFP, N_f1
         WRITE(13,'(A)') "# RM, RN, RT, RN_r, RN_p, RN_s J_ohm, J_run, J_bs, E1, E_drei, dndt_p, dndt_s, RJS " 
         IF(MODEL_connor_fp.eq.1)THEN
            DO NR=1,NRMAX
               WRITE(13,'(1P100E17.8e3)') &
                    RM(NR), RN_disrupt(NR), RT_quench(NR), & 
                    RN_runaway(NR), RN_drei(NR), RN_runaway(NR)-RN_drei(NR), &
                    RJ_ohm(NR), RJ_runaway(NR), RJ_bs(NR), &
                    E1(NR), ER_drei(NR), &
                    RFP(NR)*RN_disrupt(NR)*1.D20, RFP_AVA(NR)*RN_disrupt(NR)*1.D20, &
                    RJS(NR,1),ER_crit(NR)
            END DO
         ELSE
            DO NR=1,NRMAX
               WRITE(13,'(1P14E17.8e3)') &
                    RM(NR), RN_disrupt(NR), RT_quench(NR), & 
                    RN_runaway(NR), RN_drei(NR), RN_runaway(NR)-RN_drei(NR), &
                    RJ_ohm(NR), RJ_runaway(NR), RJ_bs(NR), &
                    E1(NR), ER_drei(NR), &
                    Rconnor(NR)*RN_disrupt(NR)*1.D20, RFP_AVA(NR)*RN_disrupt(NR)*1.D20, &
                    RJS(NR,1),ER_crit(NR)
            END DO
         END IF
         WRITE(13,'(A)') " "
         WRITE(13,'(A)') " "
! nth-re
         WRITE(14,'(A, 1PE15.6e3, i7)') "# TIME ", TIMEFP, N_f1 ! RE_pitch.dat
         DO NTH=1,NTHMAX
            WRITE(14,'(I4, 1P13E17.8e3)') &
                 NTH, RE_PITCH(NTH), SINM(NTH), THM(NTH)
         END DO
         WRITE(14,'(A)') " "
         WRITE(14,'(A)') " "

!         WRITE(19,'(A, 1PE15.6e3, i7)') "# TIME ", TIMEFP, N_f1
!         DO NP=NPSTART, NPEND
!            WRITE(19,'(I5,1P7E14.6)') &
!         NP, PM(NP,1), DCPP2(1,NP,1,1,1), DCPP2(1,NP,1,2,1) &
!                 , DCTT2(1,NP,1,1,1), DCTT2(1,NP,1,2,1) &
!                 , FCPP2(1,NP,1,1,1), FCPP2(1,NP,1,2,1)
!         END DO
!         WRITE(19,*) " "
!         WRITE(19,*) " "

         N_f1=N_f1+1
      END IF

! time_evol.dat
      WRITE(10,'(I6,1P30E16.7e3)') NT, TIMEFP, RT_TEMP(1,1), &
           RT_quench(1), RJ_ohm(1), &
           RJ_runaway(1), RN_disrupt(1), RN_runaway(1), &
           RN_drei(1), &
           IP_all, ip_ohm, ip_run, ip_prim, &
           l_ind, ZEFF, IP_all_FP, RFP(1)*tau_ta0(1),RN_MGI(1,1)

! efield_e1
      WRITE(11,'(1P256E14.6)') TIMEFP, (E1(NR), NR=1,NRMAX), (ER_drei(NR), NR=1,NRMAX), (ER_crit(NR), NR=1,NRMAX)

! dndt
      IF(MODEL_connor_fp.eq.1)THEN
         WRITE(12,'(1P128E15.6e3)') TIMEFP, (RFP(NR)*RN_disrupt(NR)*1.d20, NR=1,NRMAX), &
              (RFP_ava(NR)*RN_disrupt(NR)*1.d20, NR=1,NRMAX)
      ELSE
         WRITE(12,'(1P128E15.6e3)') TIMEFP, (Rconnor(NR)*RN_disrupt(NR)*1.d20, NR=1,NRMAX), &
              (RFP_ava(NR)*RN_disrupt(NR)*1.d20, NR=1,NRMAX)
      END IF
      
! re_pitch
      WRITE(15,'(1P128E15.6e3)') TIMEFP, (RE_PITCH(NTH), NTH=1,NTHMAX)

! e-e collision time: taus = 0.5*tau_ee
      WRITE(23,'(1P128E15.6e3)') TIMEFP, (POST_tau_ta(NR,1), NR=1,NRMAX)


      END IF
      
      END SUBROUTINE FILE_OUTPUT_DISRUPT
!**********************************************
      SUBROUTINE MGI_DENSITY

      IMPLICIT NONE
      INTEGER:: NR, NSB
      REAL(rkind):: FACTZ, tau_imp, imp_charge
      
!      tau_imp=3.D0*tau_quench
      tau_imp=tau_mgi
      imp_charge=PZ(n_impu)

      DO NR=NRSTART,NREND
!      DO NR=1,NRMAX
         DO NSB=1,NSBMAX
            IF(NSB.eq.1)THEN
               FACTZ=zeff_imp
               RN_MGI(NR,NSB) = RN_MGI(NR,NSB) + (SPITOT-RNFP0(1))/tau_imp*DELT*RNFD(NR,NSB)/RNFD0(NSB)
            ELSEIF(NSB.eq.N_impu)THEN
               FACTZ=1.D0/zeff_imp
               RN_MGI(NR,NSB) = RN_MGI(NR,NSB) + FACTZ*(SPITOT-RNFP0(1))/tau_imp*DELT*RNFD(NR,NSB)/RNFD0(NSB)
            ELSE
               FACTZ=0.D0
            END IF
         END DO
      END DO
      
      DO NSB=1,NSBMAX
         IF(NSB.eq.1)THEN
            FACTZ=1
         ELSEIF(NSB.eq.N_impu)THEN
            FACTZ=1.D0/zeff_imp
         ELSE
            FACTZ=0.D0
         END IF
         RN0_MGI(NSB) = RN0_MGI(NSB) + FACTZ*(SPITOT-RNFP0(1))/tau_imp*DELT
      END DO
      

!      WRITE(*,'(I3,6E14.6)') NRANK, RN0_MGI(1), RN0_MGI(2), RN0_MGI(3), RN_MGI(NRSTART,1), RN_MGI(NRSTART,2), RN_MGI(NRSTART,3)
    
      END SUBROUTINE MGI_DENSITY
!**********************************************

      END MODULE fpdisrupt
