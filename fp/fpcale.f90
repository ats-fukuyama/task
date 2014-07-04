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

!      integer,parameter:: isw_e=0
      integer,parameter:: isw_e=1

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE E_IND_IMPLICIT

      IMPLICIT NONE
      integer:: NR
      integer:: imtxstart1,imtxend1,its
      real(4):: gut1, gut2
      real(8):: Aij, Aij_m1, Aij_p1, RHS, b_wall, E_e, a_e, coef_ln

      CALL GUTIME(gut1)
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
            IF(ISW_E.eq.0)THEN
               RHS=SIGMA_SPM(NR)*EM(NR) &
                    - AEE*VC*((Rconner(NR)+RFP_ava(NR))*RN_disrupt(NR)*1.D20)*DELT
            ELSEIF(ISW_E.eq.1)THEN
               RHS=SIGMA_SPM(NR)*EM(NR) &
                    - AEE*VC*((RFP(NR)+RFP_ava(NR))*RN_disrupt(NR)*1.D20)*DELT
            END IF
         ELSE
            Aij_p1=-DELT/(RMU0*DELR**2*RA**2)*(1.D0+0.5D0*DELR/RM(NR))
            b_wall=1.2D0
            a_e=1.D0
            coef_ln=-LOG(b_wall/(1.D0+DELR*0.5D0))
            E_e=-coef_ln/DELR*EP(NRMAX)/(1.D0-coef_ln/DELR)
            IF(ISW_E.eq.0)THEN
               RHS=SIGMA_SPM(NR)*EM(NR) -Aij_p1*E_e &
                    - AEE*VC*((Rconner(NR)+RFP_ava(NR))*RN_disrupt(NR)*1.D20)*DELT
            ELSEIF(ISW_E.eq.1)THEN
               RHS=SIGMA_SPM(NR)*EM(NR) -Aij_p1*E_e &
                    - AEE*VC*((RFP(NR)+RFP_ava(NR))*RN_disrupt(NR)*1.D20)*DELT
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

      end MODULE fpcaleind
