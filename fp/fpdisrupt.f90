      MODULE fpdisrupt

      USE fpcomm

      real(8),parameter:: Z_ef=3.D0
!      real(8),parameter:: Z_ef=1.D0
      real(8),parameter:: IP_init=1.8D6 ! [A]
!      real(8),parameter:: IP_init=15.D6 ! [A]

!      integer,parameter:: ISW_R=0
      integer,parameter:: ISW_R=1
      contains

! ------------------------------------------

      SUBROUTINE set_initial_disrupt_param 

      USE plprof
      USE fpmpi
      USE libmpi
      IMPLICIT NONE

      INTEGER:: NR, NSA, NSB, NS, NP, NSW, N
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      REAL(8),dimension(NRSTART:NREND):: &
           RN1_temp, RN2_temp, RJ1_temp, RJ2_temp
      REAL(8),dimension(NRSTART:NREND):: RT1_temp, RT2_temp
      real(8):: RHON, SUMIP, IP_norm
      INTEGER:: VLOC

!------critical momentum
      IF(MODELR.eq.0)THEN
         pc_runaway=1.D0/SQRT(E0)
      ELSEIF(MODELR.eq.1)THEN
         pc_runaway=1.D0/SQRT(E0-THETA0(1))
      END IF
      IF(pc_runaway.ge.PG(NPMAX+1,1))THEN
         NPC_runaway=NPMAX+1
      ELSE
         DO NP=NPMAX+1,1,-1
            IF(PG(NP,1).ge.pc_runaway) NPC_runaway=NP
         END DO
      END IF

!-----initial profile of n, j, T
      SUMIP=0.D0
      DO NR=NRSTART,NREND
!      DO NR=1,NRMAX
         RN1_TEMP(NR)=RNFP(NR,1)
         RN2_TEMP(NR)=0.D0

!         RJ1_temp(NR)=(1.D0-RM(NR)**2)**2
         RJ1_temp(NR)=(1.D0-RM(NR)**0.95)**3
         RJ2_temp(NR)=0.D0
         SUMIP=RJ1_temp(NR)*VOLR(NR)/(2.D0*PI*RR) + SUMIP

         RHON=RM(NR)
         CALL PL_PROF(RHON,PLF)
!            RT_quench_f(NR)=PLF(NS)%RTPR*T0_quench
         RT1_temp(NR)=RTFD(NR,1)
!            RT2_temp(NR,NSA)=PLF(NS)%RTPR/RTFP0(NSA)*T0_quench ! sustain profile
!            RT2_temp(NR,NSA)=T0_quench ! flat profile
         RT2_temp(NR)=T0_quench*(1.D0-0.9D0*RM(NR)**2) ! smith
      END DO

      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(RN1_temp,NREND-NRSTART+1,RN_disrupt)
      call mtx_allgather_real8(RN2_temp,NREND-NRSTART+1,RN_runaway)

      CALL mtx_allreduce1_real8(SUMIP,3,IP_norm,vloc)
      DO NR=NRSTART,NREND
         RJ1_temp(NR)=RJ1_temp(NR)*IP_init/IP_norm*1.D-6
      END DO
      call mtx_allgather_real8(RJ1_temp,NREND-NRSTART+1,RJ_disrupt)
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
         Rconner(NR)=0.D0
         RFP_AVA(NR)=0.D0
         RFP(NR)=0.D0
      END DO
      DO NR=NRSTART,NREND
         RFPL(NR)=0.D0
      END DO

      END SUBROUTINE set_initial_disrupt_param
!**************************************************
      SUBROUTINE set_post_disrupt_Clog_f

      IMPLICIT NONE
      INTEGER:: NR, NSA, NSB, NS, NSFP, NSFD, ISW_CLOG
      real(8):: RNA,RNB,RTA,RTB,RLNRL
      real(8),dimension(NSAMAX):: PTFP0_f

      DO NR=NRSTART, NRENDWM
!      DO NR=1, NRMAX

!-----   post disrupt Coulomb log
         ISW_CLOG=0 ! =0 Wesson, =1 NRL
         DO NSA=1,NSAMAX
            NSFP=NS_NSB(NSA)
            RNA=RNFP(NR,NSA)
            RTA=RT_quench_f(NR)
            PTFP0_f(NSA)=SQRT(T0_quench*1.D3*AEE*AMFP(NSA))
            DO NSB=1,NSBMAX
               NSFD=NS_NSB(NSB)
               RNB=RNFD(NR,NSB)
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
!               POST_LNLAM_f(NR,NSB,NSA)=18.D0
            END DO ! NSB
            IF(NRANK.eq.0)THEN
               POST_tau_ta0_f(NSA)=(4.D0*PI*EPS0**2)*PTFP0_f(NSA)**3 &
                    /( AMFP(NSA)*AEFP(NSA)**4*POST_LNLAM_f(1,NSA,NSA)*RNFP0(NSA)*1.D20 )
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
      real(8):: RNA,RNB,RTA,RTB,RLNRL,P_thermal, rtheta
      real(8),dimension(NSAMAX):: PTFP0_f
      real(8),dimension(NRSTART:NREND,NSAMAX):: tau_ta
      real(8),dimension(NRSTART:NREND):: LNL_l

      DO NR=NRSTART, NREND
!      DO NR=1, NRMAX

!-----   post disrupt Coulomb log
         ISW_CLOG=0 ! =0 Wesson, =1 NRL
         DO NSA=1,NSAMAX
            NSFP=NS_NSB(NSA)
            RNA=RNFP(NR,NSA)
            RTA=RT_quench(NR)
            DO NSB=1,NSBMAX
               NSFD=NS_NSB(NSB)
               RNB=RNFD(NR,NSB)
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
!               POST_LNLAM(NR,NSB,NSA)=18.D0
            END DO ! NSB
            P_thermal=SQRT(RTA*1.D3*AEE*AMFP(NSA))
            tau_ta(NR,NSA)=(4.D0*PI*EPS0**2)*P_thermal**3 &
                 /( AMFP(NSA)*AEFP(NSA)**4*POST_LNLAM(NR,NSA,NSA)*RNA*1.D20 )
         END DO ! NSA
         lnl_l(NR)=POST_LNLAM(NR,1,1)
      END DO ! NR

      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(LNL_L,NREND-NRSTART+1,LNL_G)

      CALL mtx_set_communicator(comm_nsanr)
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL fp_gatherv_real8_sav(tau_ta,SAVLEN(NRANK+1),POST_tau_ta,N,NSA)
      END DO
      CALL mtx_reset_communicator
      CALL mtx_broadcast_real8(POST_tau_ta,NRMAX*NSAMAX)

      DO NR=1,NRMAX
!         ER_drei(NR)=SQRT(RT_quench(NR)*1.D3*AEE*AMFP(1))/(ABS(AEFP(1))*POST_tau_ta(NR,1) )
!         ER_crit(NR)=ER_drei(NR)*RT_quench(NR)*1.D3*AEE/(AMFP(1)*VC*VC)

         ER_drei(NR)=SQRT(RT_quench(NR)*1.D3*AEE*AMFP(1))/(ABS(AEFP(1))*POST_tau_ta(NR,1) ) &
              /lnL_G(NR)*18.D0
         ER_crit(NR)=ER_drei(NR)*RT_quench(NR)*1.D3*AEE/(AMFP(1)*VC*VC) &
              /lnL_G(NR)*18.D0
      END DO

      END SUBROUTINE set_post_disrupt_Clog
!*****************************************
      SUBROUTINE display_disrupt_initials

      IMPLICIT NONE
      real(8):: alp, z_i, h_alpha_z, lambda_alpha, gamma_alpha_z, G_conner
      real(8):: G_conner_nr, G_conner_lm
      integer:: j

      WRITE(6,*)" ---- DISRUPTION PARAM ------"
      WRITE(6,'(A,1PE14.6)') "Pre  disruption Te0 [keV]=", RTFP0(1)
      WRITE(6,'(A,1P4E14.6)') &
             "Pre  disrupt tau_ta0[sec]=", (tau_ta0(j), j=1,NSAMAX)

      IF(MODELR.eq.0)THEN
         WRITE(6,'(A,1P3E14.6,A,1PE14.6)') &
              "E0*E_drei0=E1[V/m] ->",E0, E_drei0(1), E1(1)&
              ," p_{c_runaway}=",1.D0/SQRT(E0)
      ELSE
         WRITE(6,'(A,1P4E14.6)') &
              "E0*E_drei0=E1[V/m], E_crit0 ->",E0, E_drei0(1), E1(1), E_crit0(1)
         WRITE(6,'(A,1PE14.6,A,1PE14.6)') &
              " p_{c_runaway}=",1.D0/SQRT(E0)           &
              ," p_{cr_runaway}=",1.D0/SQRT(E0-THETA0(1))
         
         alp = E0/THETA0(1)
!         z_i = PZ(2)
         z_i = Z_ef
         h_alpha_z=( alp*(z_i+1.D0) - z_i + 7.D0 + &
              2.D0*SQRT(alp/(alp-1.D0))*(1.D0+z_i)*(alp-2.D0) ) &
              /( 16.D0*(alp-1.D0) )
         lambda_alpha=8.D0*alp*(alp-0.5D0-SQRT(alp*(alp-1.D0)) )
         gamma_alpha_z=SQRT( (1.0+z_i)*alp**2/(8.D0*(alp-1.D0)) )&
              *(0.5D0*PI-ASIN(1.D0-2.D0/alp))
         
         G_conner=0.35D0* E0**(-h_alpha_z)*EXP(-0.25D0*lambda_alpha/E0 &
              -SQRT(2.D0/E0)*gamma_alpha_z )
         
         G_conner_nr=0.35D0*E0**(-3.D0*(Z_i+1.D0)/16.D0)&
              *EXP(-0.25D0/E0 -SQRT( (1.D0+z_i)/E0 ) )
         G_conner_lm=G_conner_nr &
              *EXP(-THETA0(1)*(0.125D0/E0**2 + 2.D0/3.D0/SQRT(E0**3)*SQRT(1.D0+z_i) ) )
         
         WRITE(6,'(A,1PE14.6,A,1PE14.6)') " alpha = ", alp, " G_conner= ", G_conner
         WRITE(6,'(A,1PE14.6,A,1PE14.6)') " G_Conner_nr = ", &
              G_conner_nr, " G_conner_lm = ", G_conner_lm
      END IF

      WRITE(6,*) "-----------------------"
      WRITE(6,'(A,1PE14.6)') "decay time [msec] = ", tau_quench*1000
      WRITE(6,*) "-----------------------"

      WRITE(6,'(A,1PE14.6)') "Post disruption Te0 [keV]=", T0_quench
      WRITE(6,'(A,1P4E14.6)') &
             "Post disrupt tau_ta0[sec]=", (POST_tau_ta0_f(j), j=1,NSAMAX)
      WRITE(6,*)" ---- END OF DISRUPTION PARAM ------"
      

      END SUBROUTINE display_disrupt_initials
! ****************************************************
!     only for 2 species plasma
      SUBROUTINE SPITZER_SIGMA(NR,Sigma)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      REAL(8),INTENT(OUT):: Sigma
      INTEGER:: NSA, NSB
      real(8):: fact, taue_col, RTE, RNE, RTI, RNI
      real(8):: neoc,  phi, f_T, C_, tau_rela, theta_l

      NSA=1
      NSB=2
      FACT=AEFP(NSA)**2*AEFD(NSB)**2*POST_LNLAM(NR,NSB,NSA)/(EPS0**2)
      IF(NTG2.eq.0)THEN
         RTE=RTFP(NR,NSA)
         RNE=RNFP(NR,NSA)
!         RNI=RNFD(NR,NSB)
      ELSE
         RTE=RT_quench(NR)
         RNE=RNFP(NR,NSA)
!         RNI=RNFD(NR,NSB)
      END IF
      taue_col=3.D0*SQRT((2.D0*PI)**3)/FACT &
           *SQRT(AMFP(1)*(AEE*RTE*1.D3)**3)/(RNE*1.D20)
      sigma=1.96D0*RNE*1.D20*AEFP(NSA)**2*taue_col/AMFP(NSA) ! Wesson P. 174
!      sigma= ! P. 71

!      neoc=(1.D0-SQRT(invasp))**2 ! P. 174
!      Z_EF=3.D0
      theta_l=THETA0(1)*RT_quench(NR)/RTFP0(1)
      tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
           ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RNFP0(NSA)*1.D20 )
      C_ = 0.56D0/Z_ef*(3.0D0-Z_ef)/(3.D0+Z_ef)
      f_t=1.D0 -(1.D0-EPSRM(NR))**2/ ( SQRT(1.D0-EPSRM(NR)**2)*(1.D0+1.46D0*SQRT(EPSRM(NR))) )
      phi = f_t/(1.D0 + (0.58D0+0.2D0*Z_ef)*(2.D0*RR*QLM(NR)*EPSRM(NR)**(-1.5D0) )/ &
           (3.D0*SQRT(2.D0*PI)*VC*tau_rela)/theta_l**2 )
      neoc=(1.D0-phi)*(1.D0-C_*phi)*(1.D0+0.47D0*(Z_ef-1.D0))/ &
           (Z_ef*(1.D0+0.27D0*(Z_ef-1.D0)) )

      sigma=sigma*neoc

      END SUBROUTINE SPITZER_SIGMA
! *******************************************************
      SUBROUTINE Tquench_trans

      USE fpmpi
      USE libmpi
      IMPLICIT NONE
      INTEGER:: NR, NSA, N, NSW, ISW_Q
      real(8):: k
      real(8),dimension(NRSTART:NREND):: RT1_temp

      ISW_Q=2

      IF(ISW_Q.eq.0)THEN ! linear 
         IF(TIMEFP+DELT.le.tau_quench)THEN
            DO NR=NRSTART,NREND
               k=( RT_quench_f(NR)-RTFP(NR,1) )/tau_quench
               RT1_temp(NR)=RTFP(NR,1) + k*(TIMEFP+DELT)
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
                 (RTFP(NR,1)-RT_quench_f(NR))*EXP(-(TIMEFP+DELT)/tau_quench)
         END DO
      END IF

      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(RT1_temp,NREND-NRSTART+1,RT_quench)
      CALL mtx_reset_communicator

      END SUBROUTINE Tquench_trans
! *******************************************************
      SUBROUTINE calculation_runaway_rate

      USE fpmpi
      USE libmpi
      IMPLICIT NONE
      integer:: NTH, NP, NR, NSA, N, NSW, NS, NSBA
      real(8):: FLUXS_PMAX, FFP, RSUM1, FACT
      real(8),dimension(NRSTART:NREND):: Rconner_l
      real(8):: alp, z_i, h_alpha_z, lambda_alpha, gamma_alpha_z, theta_l, E00, tau_rela

      CALL mtx_set_communicator(comm_np)
      DO NR=NRSTART, NRENDX
         DO NSA=NSASTART,NSAEND
            NSBA=NSB_NSA(NSA)
            NS=NS_NSA(NSA)
            RSUM1=0.D0
            FLUXS_PMAX=0.D0
               
            IF(MODELA.eq.0)THEN
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)
                  END DO
               ENDDO
            ELSE
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                          *RLAMDA(NTH,NR)
                  END DO
               ENDDO
            END IF
               
! FLUX S_p across pmax for runaway rate
            IF(NPEND.eq.NPMAX)THEN
               DO NTH=1,NTHMAX
                  FFP=    PG(NPMAX+1,NSBA)*FNSP(NTH,NPMAX,NR,NSBA)
                  
                  FLUXS_PMAX = FLUXS_PMAX +  &
                       FPP(NTH,NPMAX+1,NR,NSA)*FFP  &
                       * PG(NPMAX+1,NSBA)*SINM(NTH)!*tau_ta0(NSA)
               END DO
            END IF
            
            CALL p_theta_integration(RSUM1)
            CALL p_theta_integration(FLUXS_PMAX)
            FACT=RNFP0(NSA)*1.D20/RFSADG(NR)
            RNSL(NR,NSA) = RSUM1*FACT*1.D-20!*RCOEFNG(NR)
            IF(NSA.eq.1) THEN
               RFPL(NR)=2.D0*PI*DELTH* FLUXS_PMAX / RSUM1
            END IF
         END DO
      END DO
      CALL mtx_set_communicator(comm_nsa) 
      CALL mtx_broadcast_real8(RFPL,NREND-NRSTART+1)
      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(RFPL,NREND-NRSTART+1,RFP)

      CALL mtx_set_communicator(comm_nsanr) 
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL fp_gatherv_real8_sav(RNSL,SAVLEN(NRANK+1),RNS,N,NSA)
      END DO
      CALL mtx_reset_communicator
      CALL mtx_broadcast_real8(RNS,NRMAX*NSAMAX)


      DO NR=NRSTART,NREND
         IF(EP(NR).ge.ER_crit(NR))THEN
            E00 = EP(NR)/ER_drei(NR)
            theta_l=THETA0(1)*RT_quench(NR)/RTFP0(1)
            alp = E00/theta_l
!            z_i = PZ(2)
            z_i = Z_ef
            h_alpha_z=( alp*(z_i+1.D0) - z_i + 7.D0 + &
                 2.D0*SQRT(alp/(alp-1.D0))*(1.D0+z_i)*(alp-2.D0) ) &
                 /( 16.D0*(alp-1.D0) )
            lambda_alpha=8.D0*alp*(alp-0.5D0-SQRT(alp*(alp-1.D0)) )
            gamma_alpha_z=SQRT( (1.0+z_i)*alp**2/(8.D0*(alp-1.D0)) )&
                 *(0.5D0*PI-ASIN(1.D0-2.D0/alp))
            
!            Rconner_l(NR)=0.35D0* E00**(-h_alpha_z) &
            Rconner_l(NR)= E00**(-h_alpha_z) &
                 *EXP(-0.25D0*lambda_alpha/E00-SQRT(2.D0/E00)*gamma_alpha_z )
! non-rela 
            Rconner(NR)=E00**(-3.D0*(Z_i+1.D0)/16.D0)&
                 *EXP(-0.25D0/E00 -SQRT( (1.D0+z_i)/E00 ) )

            tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
                 ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RNFP0(NSA)*1.D20 )
            Rconner_l(NR)=Rconner_l(NR)/(tau_rela*SQRT(2*theta_l)**3)
            
         ELSE
            Rconner_l(NR)=0.D0
         END IF
      END DO
      CALL mtx_set_communicator(comm_nr) 
      call mtx_allgather_real8(Rconner_l,NREND-NRSTART+1,Rconner)
      CALL mtx_reset_communicator


      END SUBROUTINE calculation_runaway_rate
! *******************************************************
      SUBROUTINE update_disruption_quantities(IP_all,IP_ohm,IP_run)

      USE fpmpi
      USE libmpi
      IMPLICIT NONE
      INTEGER:: NR, NSA
      real(8),intent(out):: IP_all, IP_ohm, IP_run
      REAL(8),dimension(NRSTART:NREND):: &
           RN1_temp, RN2_temp, RJ1_temp, RJ2_temp, QLM_L, SUMIP4_L
      REAL(8),dimension(NRSTART:NREND,NSAMAX):: RT1_temp, RT2_temp
      real(8),dimension(NRMAX):: SUMIP4
      REAL(8):: SUMIP, SUMIP2, SUMIP3
      real(8),save:: integ
!      real(8),dimension(nrstart:nrend),save:: previous_rate
      INTEGER:: VLOC

      IF(TIMEFP.eq.DELT) integ=0.D0
      IF(TIMEFP.eq.DELT) previous_rate(:)=0.D0

      SUMIP=0.D0
      SUMIP3=0.D0
      SUMIP2=0.D0
      SUMIP4(:)=0.D0
      SUMIP4_L(:)=0.D0
      DO NR=NRSTART,NREND
         IF(ISW_R.eq.0)THEN
            RN2_temp(NR)=DELT*(Rconner(NR)+RFP_ava(NR))*RN_disrupt(NR)+RN_runaway(NR)
            RN1_temp(NR)=-DELT*(Rconner(NR)+RFP_ava(NR))*RN_disrupt(NR)+RN_disrupt(NR)
!            integ = (Rconner(NR)+RFP_ava(NR))*DELT + integ
!            integ = 0.5D0*( (Rconner(NR)+RFP_ava(NR))+previous_rate )*DELT + integ
!            previous_rate = (Rconner(NR)+RFP_ava(NR))
!            RN1_temp(NR) = RNFP(NR,1)*EXP(-integ)
!            RN2_temp(NR) = RNFP(NR,1)*( 1.D0- EXP(-integ) )
         ELSEIF(ISW_R.eq.1)THEN
!            RN2_temp(NR)=DELT*(RFP(NR)+RFP_ava(NR))*RN_disrupt(NR)+RN_runaway(NR)
!            RN1_temp(NR)=-DELT*(RFP(NR)+RFP_ava(NR))*RN_disrupt(NR)+RN_disrupt(NR)

            RN2_temp(NR)=0.5D0*DELT* &
                 ( (RFP(NR)+RFP_ava(NR))*RN_disrupt(NR) + previous_rate(NR) ) &
                 +RN_runaway(NR)
            RN1_temp(NR)=-0.5D0*DELT* &
                 ( (RFP(NR)+RFP_ava(NR))*RN_disrupt(NR) + previous_rate(NR) ) &
                 +RN_disrupt(NR)

            previous_rate(NR)=(RFP(NR)+RFP_ava(NR))*RN_disrupt(NR)
         END IF

         RJ2_temp(NR)=AEE*VC*RN2_temp(NR)*1.D20*1.D-6
         RJ1_temp(NR)=SIGMA_SPP(NR)*EP(NR)*1.D-6 !- RJ2_temp(NR) 
         SUMIP=(RJ1_temp(NR) +RJ2_temp(NR) )*VOLR(NR)/(2.D0*PI*RR) + SUMIP
         SUMIP4_L(NR)=SUMIP
         SUMIP2=RJ1_temp(NR)*VOLR(NR)/(2.D0*PI*RR) + SUMIP2
         SUMIP3=RJ2_temp(NR)*VOLR(NR)/(2.D0*PI*RR) + SUMIP3
      END DO
      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(RN1_temp,NREND-NRSTART+1,RN_disrupt)
      call mtx_allgather_real8(RN2_temp,NREND-NRSTART+1,RN_runaway)
      call mtx_allgather_real8(RJ1_temp,NREND-NRSTART+1,RJ_disrupt)
      call mtx_allgather_real8(RJ2_temp,NREND-NRSTART+1,RJ_runaway)
      call mtx_allgather_real8(SUMIP4_L,NREND-NRSTART+1,SUMIP4)
!
      CALL mtx_allreduce1_real8(SUMIP,3,IP_all,vloc)
      CALL mtx_allreduce1_real8(SUMIP2,3,IP_ohm,vloc)
      CALL mtx_allreduce1_real8(SUMIP3,3,IP_run,vloc)
      SUMIP=0.D0
      DO NR=1,NRMAX
         SUMIP=SUMIP+SUMIP4(NR)
         QLM(NR)=2*PI*RM(NR)**2*RA**2*BB/(RMU0*SUMIP*1.D6*RR)
      END DO
      CALL mtx_reset_communicator

!      IF(NRANK.eq.0)THEN
!         WRITE(*,*) "IP_post_disruption [MA] = ", IP_norm*1.D-6
!      END IF

      END SUBROUTINE update_disruption_quantities
! *******************************************************
      SUBROUTINE AVALANCHE

      USE fpmpi
      USE libmpi
      IMPLICIT NONE
      integer:: NR
      real(8),dimension(NRSTART:NREND):: rate_ava_l
      real(8):: E_hat, tau_rela, phi, theta_l

      DO NR=NRSTART, NREND
         tau_rela=(4.D0*PI*EPS0**2)*AMFP(1)**2*VC**3/ &
              ( AEFP(1)**4*POST_LNLAM(NR,1,1)*RNFP0(1)*1.D20 )
!         theta_l=THETA0(1)*RT_quench(NR)/RTFP0(1)
         E_hat=EP(NR)/ER_crit(NR)
         phi=1.D0-1.46D0*SQRT(EPSRM(NR))+1.72D0*EPSRM(NR)
         IF(E_hat.ge.1)THEN
            IF(ISW_R.eq.0)THEN
               rate_ava_l(NR)= &
                    RN_runaway(NR)/RN_disrupt(NR)* &
                    (E_hat-1.D0)/ &
                    (tau_rela*POST_LNLAM(NR,1,1))* &
                    SQRT(PI*phi/(3.D0*(Z_ef+5.D0)))/ &
                    SQRT(1.D0-1.D0/E_hat+ (4.D0*PI*(Z_ef+1.D0)**2)/ &
                    (3.D0*phi*(Z_ef+5.D0)*(E_hat**2+4.D0/phi**2-1.D0)) )
            ELSEIF(ISW_R.eq.1)THEN
               rate_ava_l(NR)= &
                    RN_runaway(NR)/RN_disrupt(NR)* &
                    (E_hat-1.D0)/ &
                    (tau_rela*POST_LNLAM(NR,1,1))* &
                    SQRT(PI*phi/(3.D0*(Z_ef+5.D0)))/ &
                    SQRT(1.D0-1.D0/E_hat+ (4.D0*PI*(Z_ef+1.D0)**2)/ &
                    (3.D0*phi*(Z_ef+5.D0)*(E_hat**2+4.D0/phi**2-1.D0)) )
            END IF
         ELSE
            rate_ava_l(NR)=0.D0
         END IF
      END DO

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
            DO NP=NPSTART, NPEND
               DO NTH=1, NTHMAX
                  FPP(NTH,NP,NR,NSA)= &
                       FEPP(NTH,NP,NR,NSA)+FCPP(NTH,NP,NR,NSA)
               END DO
            END DO
         END DO
      END DO

      END SUBROUTINE update_fpp
!**********************************************

      END MODULE fpdisrupt
