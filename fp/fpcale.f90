!
! *************************
!     SAVE DATA ROUTINE
! *************************
!
      MODULE fpcaleind

      USE fpcomm
      USE libmpi
      USE libmtx
      USE plprof
      USE equnit_mod

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE E_IND_IMPLICIT

      IMPLICIT NONE
      integer:: NR
      integer:: imtxstart1,imtxend1,its
      real(4):: gut1, gut2
      real(8):: Aij, Aij_m1, Aij_p1, RHS, b_wall, E_e, a_e, coef_ln
      real(8):: factr, factp, DR_coef, RN_runm
      real(8),dimension(NRMAX):: dndt

      CALL GUTIME(gut1)

      IF(MODELD_n_RE.eq.1.and.RN_runaway(1).ge.1.D-40)THEN
         CALL N_RE_TRANSPORT_dndt(dndt)
      END IF

      CALL mtx_set_communicator(comm_nr) !3D

      CALL mtx_setup(NRMAX,imtxstart1,imtxend1)
!---- DIAGONAL TERM
      DO NR=NRSTART,NREND
         Aij=SIGMA_SPP(NR)+2.D0*DELT/(RMU0*DELR**2*RA**2)
         CALL mtx_set_matrix(NR,NR,Aij)
         CALL mtx_set_vector(NR,EM(NR))
      END DO
!---- OFF DIAGONAL
      DO NR=NRSTART,NREND
         Aij_m1=-DELT/(RMU0*DELR**2*RA**2)*(1.D0-0.5D0*DELR/RM(NR))
         Aij_p1=-DELT/(RMU0*DELR**2*RA**2)*(1.D0+0.5D0*DELR/RM(NR))
         IF(NR.ne.1)THEN
            CALL mtx_set_matrix(NR,NR-1,Aij_m1)
         END IF
         IF(NR.ne.NRMAX)THEN
            CALL mtx_set_matrix(NR,NR+1,Aij_p1)
         END IF
      END DO
!---- RIGHT HAND SIDE
      DO NR=NRSTART,NREND
         IF(NR.ne.NRMAX)THEN
            IF(MODELD_n_RE.eq.1.and.RN_runaway(1).ge.1.D-40)THEN
               FACTR=QLM(NR)*deltaB_B**2
               FACTP=PI*RR*VC/SQRT(2.D0)
               DR_coef=FACTP*FACTR/(RA**2)
               IF(NR.eq.1)THEN
                  rn_runm=0.D0
               ELSE
                  rn_runm=rn_runaway(NR)
               END IF
               IF(MODEL_Conner_FP.eq.0)THEN
                  RHS=SIGMA_SPM(NR)*EM(NR) &
                       -&
!                       ( &
!                       DR_coef/(RM(NR)*DELR**2)*( &
!                       RG(NR+1)*RN_runaway(NR+1)-2.D0*RM(NR)*RN_runaway(NR)+RG(NR)*RN_runm&
!                       )*1.D20 &
!                       + ((Rconner(NR)+RFP_ava(NR))*RN_disrupt(NR)*1.D20)&
!                       )&
                       dndt(NR)*1.D20 &
                       *AEE*VC*DELT &
                       - (RJ_bs(NR)-RJ_bsm(NR))*1.d6
               ELSE
                  RHS=SIGMA_SPM(NR)*EM(NR) &
                       - &
!                       ( &
!                       DR_coef/(RM(NR)*DELR**2)*( &
!                       RG(NR+1)*RN_runaway(NR+1)-2.D0*RM(NR)*RN_runaway(NR)+RG(NR)*RN_runm &
!                       )*1.D20 &
!                       + ((RFP(NR)+RFP_ava(NR))*RN_disrupt(NR)*1.D20)&
!                       ) &
                       dndt(NR)*1.D20 &
                       *AEE*VC*DELT &
                       - (RJ_bs(NR)-RJ_bsm(NR))*1.d6
               END IF
            ELSE ! no RE transport
               IF(MODEL_Conner_FP.eq.0)THEN ! default
!               IF(MODEL_Conner_FP.eq.1)THEN ! for given E prof. 
                  RHS=SIGMA_SPM(NR)*EM(NR) &
                       - AEE*v_RE*VC*((Rconner(NR)+RFP_ava(NR))*RN_disrupt(NR)*1.D20)*DELT &
                       - (RJ_bs(NR)-RJ_bsm(NR))*1.d6
               ELSE
                  RHS=SIGMA_SPM(NR)*EM(NR) &
                       - AEE*v_RE*VC*((RFP(NR)+RFP_ava(NR))*RN_disrupt(NR)*1.D20)*DELT &
                       - (RJ_bs(NR)-RJ_bsm(NR))*1.d6
               END IF
            END IF
         ELSE
! ------ boundary condition
            Aij_p1=-DELT/(RMU0*DELR**2*RA**2)*(1.D0+0.5D0*DELR/RM(NR))
            b_wall=RB
            a_e=RA
            coef_ln=-LOG(b_wall/(a_e*(1.D0+DELR*0.5D0) ) )*(b_wall/(a_e*DELR))
            E_e=-coef_ln*EP(NRMAX)/(1.D0-coef_ln)
            IF(MODELD_n_RE.eq.1.and.RN_runaway(1).ge.1.D-40)THEN
               FACTR=QLM(NR)*deltaB_B**2
               FACTP=PI*RR*VC/SQRT(2.D0)
               DR_coef=FACTP*FACTR/(RA**2)
               IF(MODEL_Conner_FP.eq.0)THEN
                  RHS=SIGMA_SPM(NR)*EM(NR) -Aij_p1*E_e &
                       - &
!                       ( &
!                       DR_coef/(RM(NR)*DELR**2)*( &
!                       -2.D0*RM(NR)*RN_runaway(NR)+RG(NR)*RN_runaway(NR-1)&
!                       )*1.D20 &
!                       + ((Rconner(NR)+RFP_ava(NR))*RN_disrupt(NR)*1.D20)&
!                       ) &
                       dndt(NR)*1.D20 &
                       *AEE*VC*DELT &
                       - (RJ_bs(NR)-RJ_bsm(NR))*1.d6
               ELSE
                  RHS=SIGMA_SPM(NR)*EM(NR) -Aij_p1*E_e &
                       - &
!                       ( &
!                       DR_coef/(RM(NR)*DELR**2)*( &
!                       -2.D0*RM(NR)*RN_runaway(NR)+RG(NR)*RN_runaway(NR-1) &
!                       )*1.D20 &
!                       + ((RFP(NR)+RFP_ava(NR))*RN_disrupt(NR)*1.D20)&
!                       ) &
                       dndt(NR)*1.D20 &
                       *AEE*VC*DELT &
                       - (RJ_bs(NR)-RJ_bsm(NR))*1.d6
               END IF
            ELSE ! no RE transport
               IF(MODEL_Conner_FP.eq.0)THEN ! default
!               IF(MODEL_Conner_FP.eq.1)THEN ! for given E prof.
                  RHS=SIGMA_SPM(NR)*EM(NR) -Aij_p1*E_e &
                       - AEE*v_RE*VC*((Rconner(NR)+RFP_ava(NR))*RN_disrupt(NR)*1.D20)*DELT &
                       - (RJ_bs(NR)-RJ_bsm(NR))*1.d6
               ELSE
                  RHS=SIGMA_SPM(NR)*EM(NR) -Aij_p1*E_e &
                       - AEE*v_RE*VC*((RFP(NR)+RFP_ava(NR))*RN_disrupt(NR)*1.D20)*DELT &
                       - (RJ_bs(NR)-RJ_bsm(NR))*1.d6
               END IF
            END IF
         END IF
         CALL mtx_set_source(NR,RHS)
      END DO
!---- SOLVE
      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) 
      CALL mtx_gather_vector(E1)
      CALL mtx_cleanup
      CALL mtx_reset_communicator

      CALL GUTIME(gut2)

      END SUBROUTINE E_IND_IMPLICIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE E_IND_EXPLICIT

      IMPLICIT NONE
      integer:: NR
      real(4):: gut1, gut2
      real(8):: coef, term1, term2, term3, term4
      real(8),dimension(NRSTART:NREND):: EP_local

      CALL GUTIME(gut1)

      CALL mtx_set_communicator(comm_nr) !3D

      DO NR=NRSTART,NREND
         coef=DELT/(RMU0*RM(NR)*DELR**2)
         term2=(SIGMA_SPM(NR)-2*DELT/(RMU0*DELR**2) )*EM(NR)
         term4=-AEE*VC*(RFP(NR)*RN_disrupt(NR)*1.D20)*DELT 
         IF(NR.eq.1)THEN
            term1=coef*(RM(NR)+DELR)*EM(NR+1)
            term3=0.D0
         ELSEIF(NR.eq.NRMAX)THEN
            term1=0.D0
            term3=coef*(RM(NR)-DELR)*EM(NR-1)
         ELSE
            term1=coef*(RM(NR)+DELR)*EM(NR+1)
            term3=coef*(RM(NR)-DELR)*EM(NR-1)
         END IF
         EP_local(NR)=(term1+term3+term4+term2)/SIGMA_SPP(NR)
      END DO

!      WRITE(*,'(A,I3,1P4E14.6)') "EXPLICIT", NRANK, term1, term2, term3, term4
      call mtx_allgather_real8(EP_local,NREND-NRSTART+1,EP)
      CALL mtx_reset_communicator

      CALL GUTIME(gut2)
      IF(NRANK.eq.0) WRITE(*,*) "GUT_EVOL_ITERATION = ",gut2-gut1

      END SUBROUTINE E_IND_EXPLICIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      SUBROUTINE E_IND_IMPLICIT_FP(E_SIGMA) ! MODEL_jfp=1

      IMPLICIT NONE
      integer:: NR
      integer:: imtxstart1,imtxend1,its
      real(8),dimension(NRSTART:NREND),INTENT(IN):: E_SIGMA
      real(4):: gut1, gut2
      real(8):: Aij, Aij_m1, Aij_p1, RHS, b_wall, E_e, a_e, coef_ln
      real(8):: RJ_P

      CALL GUTIME(gut1)
      CALL mtx_set_communicator(comm_nr) !3D

      CALL mtx_setup(NRMAX,imtxstart1,imtxend1)
!---- DIAGONAL TERM
      DO NR=NRSTART,NREND
         RJ_P=RJS(NR,1)*1.D6
!         Aij=2.D0*DELT/(RMU0*DELR**2*RA**2)
!         Aij=( 2.D0*DELT/(RMU0*DELR**2*RA**2) +RJ_P/EP(NR) )
         Aij=( 2.D0*DELT/(RMU0*DELR**2*RA**2) +RJ_P/E_SIGMA(NR) )
         CALL mtx_set_matrix(NR,NR,Aij)
!         CALL mtx_set_vector(NR,EM(NR))
         CALL mtx_set_vector(NR,EP(NR))
      END DO
!---- OFF DIAGONAL
      DO NR=NRSTART,NREND
         Aij_m1=-DELT/(RMU0*DELR**2*RA**2)*(1.D0-0.5D0*DELR/RM(NR))
         Aij_p1=-DELT/(RMU0*DELR**2*RA**2)*(1.D0+0.5D0*DELR/RM(NR))
         IF(NR.ne.1)THEN
            CALL mtx_set_matrix(NR,NR-1,Aij_m1)
         END IF
         IF(NR.ne.NRMAX)THEN
            CALL mtx_set_matrix(NR,NR+1,Aij_p1)
         END IF
      END DO
!---- RIGHT HAND SIDE
      DO NR=NRSTART,NREND
         IF(NR.ne.NRMAX)THEN
            IF(MODEL_Conner_FP.eq.0)THEN
!               RHS= -( RJ_P-RJS_M(NR,1) )*1.D6 &
               RHS= -( -RJS_M(NR,1) )*1.D6 &
                    - AEE*VC*((Rconner(NR)+RFP_ava(NR))*RNS(NR,1)*1.D20)*DELT &
                    - (RJ_bs(NR)-RJ_bsm(NR))*1.d6
            ELSEIF(MODEL_Conner_FP.eq.1)THEN
!               RHS= -( RJ_P-RJS_M(NR,1) )*1.D6 &
               RHS= -( -RJS_M(NR,1) )*1.D6 &
                    - AEE*VC*((RFP(NR)+RFP_ava(NR))*RNS(NR,1)*1.D20)*DELT &
                    - (RJ_bs(NR)-RJ_bsm(NR))*1.d6
            END IF
         ELSE
! ------ boundary condition
            Aij_p1=-DELT/(RMU0*DELR**2*RA**2)*(1.D0+0.5D0*DELR/RM(NR))
            b_wall=RB
            a_e=RA
            coef_ln=-LOG(b_wall/(a_e*(1.D0+DELR*0.5D0) ) )*(b_wall/(a_e*DELR))
!            E_e=-coef_ln*EP(NRMAX)/(1.D0-coef_ln)
            E_e=-coef_ln*E_SIGMA(NRMAX)/(1.D0-coef_ln)
            IF(MODEL_Conner_FP.eq.0)THEN
!               RHS= -( RJ_P-RJS_M(NR,1) )*1.D6 &
               RHS= -( -RJS_M(NR,1) )*1.D6 &
                    - Aij_p1*E_e &
                    - AEE*VC*((Rconner(NR)+RFP_ava(NR))*RNS(NR,1)*1.D20)*DELT &
                    - (RJ_bs(NR)-RJ_bsm(NR))*1.d6
            ELSEIF(MODEL_Conner_FP.eq.1)THEN
!               RHS= -( RJ_P-RJS_M(NR,1) )*1.D6 &
               RHS= -( -RJS_M(NR,1) )*1.D6 &
                    - Aij_p1*E_e &
                    - AEE*VC*((RFP(NR)+RFP_ava(NR))*RNS(NR,1)*1.D20)*DELT &
                    - (RJ_bs(NR)-RJ_bsm(NR))*1.d6
            END IF
         END IF
         CALL mtx_set_source(NR,RHS)
      END DO
!---- SOLVE
      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) 
!      IF(NRANK.eq.0) write(6,*) 'E_IND_EVOL, Number of iterations    =',its
      CALL mtx_gather_vector(E1)
      CALL mtx_cleanup
      CALL mtx_reset_communicator

      CALL GUTIME(gut2)
!      IF (MOD(NT,NTG1STEP).EQ.0) THEN
!         IF(NRANK.eq.0) WRITE(*,*) "GUT_EVOL_ITERATION = ",gut2-gut1
!      END IF

      END SUBROUTINE E_IND_IMPLICIT_FP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      SUBROUTINE E_IND_CRANK_NICOLSON_FP

      IMPLICIT NONE
      integer:: NR
      integer:: imtxstart1,imtxend1,its
      real(4):: gut1, gut2
      real(8):: Aij, Aij_m1, Aij_p1, RHS, b_wall, E_e, a_e, coef_ln, RJ_P

      CALL GUTIME(gut1)
      CALL mtx_set_communicator(comm_nr) !3D

      CALL mtx_setup(NRMAX,imtxstart1,imtxend1)
!---- Boundary condition
      IF(NREND.eq.NRMAX)THEN
         Aij_p1=-DELT/(RMU0*DELR**2*RA**2)*(1.D0+0.5D0*DELR/RM(NRMAX))
         b_wall=RB
         a_e=RA
         coef_ln=-LOG(b_wall/(a_e*(1.D0+DELR*0.5D0) ) )*(b_wall/(a_e*DELR))
         E_e=-coef_ln*EP(NRMAX)/(1.D0-coef_ln)
!         coef_ln=-LOG(b_wall/(1.D0+DELR*0.5D0))
!         E_e=-coef_ln/DELR*EP(NRMAX)/(1.D0-coef_ln/DELR)
         EM_W(NRMAX+1)=-coef_ln/DELR*EM(NRMAX)/(1.D0-coef_ln/DELR)
      END IF

!---- DIAGONAL TERM
      DO NR=NRSTART,NREND
         RJ_P=RJS(NR,1)*1.D6
!         Aij=1.D0*DELT/(RMU0*DELR**2*RA**2)
         Aij=( 1.D0*DELT/(RMU0*DELR**2*RA**2) +RJ_P/EP(NR) )
         CALL mtx_set_matrix(NR,NR,Aij)
         CALL mtx_set_vector(NR,EM(NR))
      END DO
!---- OFF DIAGONAL
      DO NR=NRSTART,NREND
         Aij_m1=-0.5D0*DELT/(RMU0*DELR**2*RA**2)*(1.D0-0.5D0*DELR/RM(NR))
         Aij_p1=-0.5D0*DELT/(RMU0*DELR**2*RA**2)*(1.D0+0.5D0*DELR/RM(NR))
         IF(NR.ne.1)THEN
            CALL mtx_set_matrix(NR,NR-1,Aij_m1)
         END IF
         IF(NR.ne.NRMAX)THEN
            CALL mtx_set_matrix(NR,NR+1,Aij_p1)
         END IF
      END DO
!---- RIGHT HAND SIDE
      DO NR=NRSTART,NREND
         Aij=1.D0*DELT/(RMU0*DELR**2*RA**2)
         Aij_m1=0.5D0*DELT/(RMU0*DELR**2*RA**2)*(1.D0-0.5D0*DELR/RM(NR))
         Aij_p1=0.5D0*DELT/(RMU0*DELR**2*RA**2)*(1.D0+0.5D0*DELR/RM(NR))
         IF(NR.ne.NRMAX)THEN
            IF(MODEL_Conner_FP.eq.0)THEN
               RHS= -( -RJS_M(NR,1) )*1.D6 &
                    - AEE*VC*((Rconner(NR)+RFP_ava(NR))*RNS(NR,1)*1.D20)*DELT &
                    - (RJ_bs(NR)-RJ_bsm(NR))*1.d6 &
                    + Aij_m1*EM_W(NR-1) - Aij*EM_W(NR) + Aij_p1*EM_W(NR+1)
            ELSEIF(MODEL_Conner_FP.eq.1)THEN
               RHS= -( -RJS_M(NR,1) )*1.D6 &
                    - AEE*VC*((RFP(NR)+RFP_ava(NR))*RNS(NR,1)*1.D20)*DELT &
                    - (RJ_bs(NR)-RJ_bsm(NR))*1.d6 &
                    + Aij_m1*EM_W(NR-1) - Aij*EM_W(NR) + Aij_p1*EM_W(NR+1)
            END IF
         ELSE
! ------ boundary condition
            IF(MODEL_Conner_FP.eq.0)THEN
               RHS= -( -RJS_M(NR,1) )*1.D6 &
                    - Aij_p1*E_e &
                    - AEE*VC*((Rconner(NR)+RFP_ava(NR))*RNS(NR,1)*1.D20)*DELT &
                    - (RJ_bs(NR)-RJ_bsm(NR))*1.d6 &
                    + Aij_m1*EM_W(NR-1) - Aij*EM_W(NR) + Aij_p1*EM_W(NR+1)
            ELSEIF(MODEL_Conner_FP.eq.1)THEN
               RHS= -( -RJS_M(NR,1) )*1.D6 &
                    - Aij_p1*E_e &
                    - AEE*VC*((RFP(NR)+RFP_ava(NR))*RNS(NR,1)*1.D20)*DELT &
                    - (RJ_bs(NR)-RJ_bsm(NR))*1.d6 &
                    + Aij_m1*EM_W(NR-1) - Aij*EM_W(NR) + Aij_p1*EM_W(NR+1)
            END IF
         END IF
         CALL mtx_set_source(NR,RHS)
      END DO
!---- SOLVE
      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) 
!      IF(NRANK.eq.0) write(6,*) 'E_IND_EVOL, Number of iterations    =',its
      CALL mtx_gather_vector(E1)
      CALL mtx_cleanup
      CALL mtx_reset_communicator

      CALL GUTIME(gut2)
!      IF (MOD(NT,NTG1STEP).EQ.0) THEN
!         IF(NRANK.eq.0) WRITE(*,*) "GUT_EVOL_ITERATION = ",gut2-gut1
!      END IF

      END SUBROUTINE E_IND_CRANK_NICOLSON_FP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      SUBROUTINE E_IND_PSI_IMPLICIT_FP

      IMPLICIT NONE
      integer:: NR
      integer:: imtxstart1,imtxend1,its
      real(8):: Aij, Aij_m1, Aij_p1, RHS, b_wall, E_e, a_e, coef_ln
      real(8):: RJ_P, SUMB_pol_p, SUMB_POL_M, SUM_PSI_P, SUM_PSI_M 
      real(8),dimension(NRMAX):: RB_pol_P, RB_pol_M
      real(8):: RPSI_POL_P, RPSI_POL_M

      SUMB_pol_P=0.D0
      SUMB_pol_M=0.D0
      DO NR=1,NRMAX
         SUMB_pol_P = SUMB_pol_P + (RJS(NR,1)*1.D6+AEE*VC*RN_runaway(NR)*1.D20)*VOLR(NR)/(2.D0*PI*RR)
         RB_pol_P(NR) = SUMB_pol_P*RMU0/(RM(NR)*RA)
         SUMB_pol_M = SUMB_pol_M + (RJS_M(NR,1)*1.D6+AEE*VC*RN_runaway_M(NR)*1.D20)*VOLR(NR)/(2.D0*PI*RR)
         RB_pol_M(NR) = SUMB_pol_M*RMU0/(RM(NR)*RA)
      END DO

      SUM_PSI_P=0.D0
      SUM_PSI_M=0.D0
      DO NR=1, NRSTART
         SUM_PSI_P = SUM_PSI_P + RB_POL_P(NR)*DELR*RA
         SUM_PSI_M = SUM_PSI_M + RB_POL_M(NR)*DELR*RA
      END DO
      RPSI_POL_P = SUM_PSI_P
      RPSI_POL_M = SUM_PSI_M

      CALL mtx_set_communicator(comm_nr) !3D
      CALL mtx_setup(NRMAX,imtxstart1,imtxend1)
!---- DIAGONAL TERM
      DO NR=NRSTART,NREND
         Aij=1.D0
         CALL mtx_set_matrix(NR,NR,Aij)
         CALL mtx_set_vector(NR,EM(NR))
      END DO
!---- OFF DIAGONAL
      DO NR=NRSTART,NREND
         Aij_p1=-1.D0
         IF(NR.ne.NRMAX)THEN
            CALL mtx_set_matrix(NR,NR+1,Aij_p1)
         END IF
      END DO
!---- RIGHT HAND SIDE
      DO NR=NRSTART,NREND
         IF(NR.ne.NRMAX)THEN
            IF(MODEL_Conner_FP.eq.0)THEN
               RHS= -( RB_POL_P(NR+1)*DELR*RA - RB_POL_M(NR+1)*DELR*RA)/DELT
            ELSEIF(MODEL_Conner_FP.eq.1)THEN
               RHS= -( RB_POL_P(NR+1)*DELR*RA - RB_POL_M(NR+1)*DELR*RA)/DELT
            END IF
         ELSE
! ------ boundary condition
            Aij_p1=-DELT/(RMU0*DELR**2*RA**2)*(1.D0+0.5D0*DELR/RM(NR))
            b_wall=RB
            a_e=RA
            coef_ln=-LOG(b_wall/(a_e*(1.D0+DELR*0.5D0) ) )*(b_wall/(a_e*DELR))
            E_e=-coef_ln*EP(NRMAX)/(1.D0-coef_ln)
            IF(MODEL_Conner_FP.eq.0)THEN
               RHS= -( RB_POL_P(NRMAX)*RM(NRMAX)/(RM(NRMAX)+DELR)*DELR &
                    -RB_POL_M(NRMAX)*RM(NRMAX)/(RM(NRMAX)+DELR)*DELR )/DELT +E_e
            ELSEIF(MODEL_Conner_FP.eq.1)THEN
               RHS= -( RB_POL_P(NRMAX)*RM(NRMAX)/(RM(NRMAX)+DELR)*DELR &
                    -RB_POL_M(NRMAX)*RM(NRMAX)/(RM(NRMAX)+DELR)*DELR )/DELT +E_e
            END IF
         END IF
         CALL mtx_set_source(NR,RHS)
      END DO
!---- SOLVE
      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) 
!      IF(NRANK.eq.0) write(6,*) 'E_IND_EVOL, Number of iterations    =',its
      CALL mtx_gather_vector(E1)
      CALL mtx_cleanup
      CALL mtx_reset_communicator


      END SUBROUTINE E_IND_PSI_IMPLICIT_FP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      SUBROUTINE N_RE_TRANSPORT_dndt(dndt)

      USE libmpi
      USE libmtx
      IMPLICIT NONE
      INTEGER:: NR
      real(8),dimension(NRMAX),intent(out):: dndt
      integer:: imtxstart1,imtxend1,its
      real(8):: Dr_coef, factp, factr
      real(8):: Aij, Aij_m1, Aij_p1, RHS, b_wall, E_e, a_e, coef_ln
      real(8),dimension(NRMAX):: RN_new
      real(8):: sum1, sum2

      CALL mtx_set_communicator(comm_nr) 
      FACTP=PI*RR*VC/SQRT(2.D0)
      RN_new(:)=0.D0

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
         IF(MODEL_Conner_FP.eq.0)THEN
            RHS = RN_runaway(NR) + (Rconner(NR)+RFP_ava(NR))*RN_disrupt(NR)*DELT
         ELSE
            RHS = RN_runaway(NR) + (RFP(NR)+RFP_ava(NR))*RN_disrupt(NR)*DELT
         END IF
         CALL mtx_set_source(NR,RHS)
      END DO
!---- SOLVE
      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) 
      CALL mtx_gather_vector(RN_new)
      CALL mtx_cleanup
      CALL mtx_reset_communicator

      sum1=0.D0
      sum2=0.D0
      DO NR=1,NRMAX
         dndt(NR)=(RN_new(NR)-RN_runaway(NR))/DELT
         IF(NRANK.eq.0) WRITE(*,'(A,I3,4E14.6)') "test NR",NR ,dndt(NR), &
              (Rconner(NR)+RFP_ava(NR))*RN_disrupt(NR), &
              RN_new(NR), RN_runaway(NR)+(Rconner(NR)+RFP_ava(NR))*RN_disrupt(NR)*DELT
         sum1 = sum1 + dndt(NR)*VOLR(NR)
         sum2 = sum2 +  (Rconner(NR)+RFP_ava(NR))*RN_disrupt(NR)*VOLR(NR)
      END DO
      IF(NRANK.eq.0) WRITE(*,*) "sum1 sum2 ", sum1, sum2
      IF(NRANK.eq.0) WRITE(*,*) " "

      END SUBROUTINE N_RE_TRANSPORT_dndt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      SUBROUTINE READ_E_FIELD

      USE libmpi
      IMPLICIT NONE
      INTEGER:: i,j
      double precision:: time_read
      double precision,dimension(NRMAX):: read_E

      CALL mtx_reset_communicator

      time_read=0.D0
      IF(NRANK.eq.0)THEN
         open(18,file='efield_ref.dat',status='old')

         DO while(timefp.ge.time_read)
            read(18,'(33E14.6)') time_read, (read_E(j),j=1,NRMAX)
!            WRITE(*,'(A,3E14.6)') "TEST",timefp, time_read, read_E(1)
         END DO
         DO i=1,NRMAX
            E1(i)=read_E(i)
         END DO
         close(18)
      END IF

      CALL mtx_broadcast_real8(E1,NRMAX)

      END SUBROUTINE READ_E_FIELD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      end MODULE fpcaleind
