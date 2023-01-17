      MODULE fpcalj

      USE fpcomm
      USE fpcalc
      USE libmpi
      USE libmtx
      IMPLICIT NONE

      contains
!**********************************************
      SUBROUTINE bootstrap_current(NR,RTE_keV,RTI_keV,RNE_norm,jbs)

      USE plprof
      IMPLICIT NONE
      INTEGER,intent(in):: NR
      double precision,dimension(NRMAX),intent(in):: RTE_keV, RTI_keV, RNE_norm
      real(8),intent(out):: jbs ! [MA]
      real(8):: BP, RNE, RNI, RTE, RTI, dndr, dtedr, dtidr, rhon

      RTE=RTE_keV(NR)*1.D3*AEE
      RTI=RTI_keV(NR)*1.D3*AEE
      RNE=RNE_norm(NR)*1.D20
      IF(NRMAX.ne.1)THEN
         IF(NR.eq.1)THEN
            dndr=(RNE_norm(NR+1)-RNE_norm(NR))/(DELR*RA)*1.D20
            dtedr=(RTE_keV(NR+1)-RTE_keV(NR))/(DELR*RA)*1.D3*AEE
            dtidr=(RTI_keV(NR+1)-RTI_keV(NR))/(DELR*RA)*1.D3*AEE
         ELSEIF(NR.eq.NRMAX)THEN
            dndr=(RNE_norm(NRMAX)-RNE_norm(NR-1))/(DELR*RA)*1.D20
            dtedr=(RTE_keV(NRMAX)-RTE_keV(NR-1))/(DELR*RA)*1.D3*AEE
            dtidr=(RTI_keV(NRMAX)-RTI_keV(NR-1))/(DELR*RA)*1.D3*AEE
         ELSE
            dndr=(RNE_norm(NR+1)-RNE_norm(NR-1))/(2.D0*DELR*RA)*1.D20
            dtedr=(RTE_keV(NR+1)-RTE_keV(NR-1))/(2.D0*DELR*RA)*1.D3*AEE
            dtidr=(RTI_keV(NR+1)-RTI_keV(NR-1))/(2.D0*DELR*RA)*1.D3*AEE
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
!      BP= RSRHON(RHON)*BB/(RR*QLM(NR))
      BP= RA*RM(NR)*BB/(RR*QLM(NR))

!     Wesson P. 173
      jbs=-SQRT(EPSRM2(NR))*RNE/BP* &
           (2.44D0*(RTE+RTI)*dndr/RNE + &
           0.69D0*dtedr - 0.42D0*dtidr) *1.D-6! *1.D-10

      jbs = ABS(jbs)

!      IF(NRANK.eq.0.and.nr.eq.1)THEN
!         WRITE(6,'(A,7E14.6)') "TEST_BS=",jbs, bp, rte, rti, rne, dndr, dtedr
!      END IF

      END SUBROUTINE bootstrap_current
!**********************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      SUBROUTINE UPDATE_QLM(NR_IN)

      IMPLICIT NONE
      INTEGER,intent(IN):: NR_IN
      INTEGER:: NR, NSA
      double precision:: RSUM_Ip, RSUM_j, Bp, Bp_init
      
      RSUM_Ip=0.D0
      IF(NR_IN.ne.1)THEN
         DO NR=1, NR_IN-1
            RSUM_j=0.D0
            DO NSA=1, NSAMAX
               RSUM_j = RSUM_j + RJS(NR,NSA) + RJES(NR,NSA)
            END DO
            RSUM_j = RSUM_j + RJ_BS(NR) + RJ_IND(NR)
            RSUM_Ip = RSUM_Ip + RSUM_j*VOLR(NR)/(2.D0*PI*RR)
         END DO
      END IF
      RSUM_j=0.D0
      DO NSA=1, NSAMAX
         RSUM_j = RSUM_j + RJS(NR_IN,NSA) + RJES(NR_IN,NSA)
      END DO
      RSUM_j = RSUM_j + RJ_BS(NR_IN) + RJ_IND(NR_IN)
      RSUM_Ip = RSUM_Ip + RSUM_j*VOLR(NR_IN)/(2.D0*PI*RR)*0.5D0

      Bp_init = RA*RM(NR_IN)*BB/(RR*QLM_INIT(NR_IN))
      Bp = RMU0*RSUM_Ip*1.D6/(2.D0*PI*RA*RM(NR_IN)) + Bp_init

      QLM(NR_IN) = RA*RM(NR_IN)*BB/(RR*Bp)

!      IF(NPSTART.eq.1) &
!           WRITE(*,'(A,I4,4E14.6)') "QLMupdate, ", NR_IN, RSUM_Ip, Bp_init, Bp, QLM(NR_IN)

      END SUBROUTINE UPDATE_QLM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---- for (NB, RF) current drive including inductive current
!-----Need MODEL_CD=1
      SUBROUTINE TOP_OF_TIME_LOOP_CD(RJ_NIND_M)

      IMPLICIT NONE
      integer:: NR, NSA
      double precision:: SUMJ, sigma, para2toro, jbs
      double precision,dimension(NRMAX),intent(out):: RJ_NIND_M
      double precision,dimension(NRMAX):: RTE_TEMP, RTI_TEMP, RNE_TEMP

      DO NR=1, NRMAX
         RTE_TEMP(NR)=RT_BULK(NR,1)
         RTI_TEMP(NR)=RT_BULK(NR,2)
         RNE_TEMP(NR)=RN_TEMP(NR,1)
      END DO

!     SUM non-inductive current
      RJ_BS(:) = 0.D0
      RJ_IND(:) = 0.D0
      RJ_NIND_M(:) = 0.D0
      DO NR=1, NRMAX
         para2toro=1.D0/SQRT(1.D0+ (RM(NR)*RA/(QLM(NR)*RR) )**2 )
         CALL bootstrap_current(NR,RTE_TEMP,RTI_TEMP,RNE_TEMP,jbs)
         RJ_BS(NR)=jbs*para2toro
         SUMJ=0.D0
         DO NSA=1,NSAMAX 
            SUMJ = SUMJ + RJS(NR,NSA) + RJES(NR,NSA)
         END DO
         RJ_NIND_M(NR) = SUMJ + RJ_BS(NR)
      END DO

      DO NR=NRSTART, NREND
         CALL SPITZER_SIGMA(NR,RT_BULK(NR,1),RN_TEMP(NR,1),given_zeff,sigma)
         para2toro=1.D0/SQRT(1.D0+ (RM(NR)*RA/(QLM(NR)*RR) )**2 )
!         SIGMA_SPM(NR)=sigma*para2toro
!         SIGMA_SPM(NR)=sigma*SQRT(para2toro**2+0.25D0*(1.D0-para2toro**2))
         SIGMA_SPM(NR)=sigma*SQRT(para2toro**2+(1.D0-para2toro**2)/1.96D0)
      END DO
      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(SIGMA_SPM,NREND-NRSTART+1,conduct_sp) 
      CALL mtx_reset_communicator

      DO NR=NRSTART, NREND
         EM(NR)=E1(NR)
         EP(NR)=E1(NR)
!         IF(E1(NR).eq.0.D0)THEN
!            EM(NR)=1.D-12
!            EP(NR)=1.D-12
!         ELSE
            EM(NR)=E1(NR)
            EP(NR)=E1(NR)
!         END IF
      END DO

      IF(NRANK.eq.0)THEN
         DO NR=NRSTART, NREND
            para2toro=1.D0/SQRT(1.D0+ (RM(NR)*RA/(QLM(NR)*RR) )**2 )
!            WRITE(*,'(A, 5E14.6)') "SIGMA_SPM, ", SIGMA_SPM(NR), conduct_sp(NR), &
!                 EM(NR), EP(NR), para2toro
         END DO
      END IF

      END SUBROUTINE TOP_OF_TIME_LOOP_CD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE AFTER_DO_WHILE_CD(RJ_NIND_M, NT)

      USE fpcaleind
      IMPLICIT NONE
      integer,intent(in):: NT
      double precision,dimension(NRMAX),intent(in):: RJ_NIND_M
      double precision,dimension(NRMAX):: EP_M, EM_M
      integer:: NR, NSA, N_IMPL_E
      double precision,dimension(NRSTART:NREND):: djdt_nind
      double precision:: RSUM, DEPS, sigma, para2toro, jbs
      double precision,dimension(NRMAX):: RTE_TEMP, RTI_TEMP, RNE_TEMP

      DO NR=1, NRMAX
         RTE_TEMP(NR)=RT_BULK(NR,1)
         RTI_TEMP(NR)=RT_BULK(NR,2)
         RNE_TEMP(NR)=RN_TEMP(NR,1)
      END DO
      DO NR=1, NRMAX
         para2toro=1.D0/SQRT(1.D0+ (RM(NR)*RA/(QLM(NR)*RR) )**2 )
         CALL bootstrap_current(NR,RTE_TEMP,RTI_TEMP,RNE_TEMP,jbs)
         RJ_BS(NR)=jbs*para2toro
      END DO

      DO NR=NRSTART, NREND
         CALL SPITZER_SIGMA(NR,RT_BULK(NR,1),RN_TEMP(NR,1),given_zeff,sigma)
         para2toro=1.D0/SQRT(1.D0+ (RM(NR)*RA/(QLM(NR)*RR) )**2 )
!         SIGMA_SPP(NR)=sigma*para2toro
!         SIGMA_SPP(NR)=sigma*SQRT(para2toro**2+0.25D0*(1.D0-para2toro**2))
         SIGMA_SPP(NR)=sigma*SQRT(para2toro**2+(1.D0-para2toro**2)/1.96D0)
!         SIGMA_SPP(NR)=sigma
      END DO
      CALL mtx_set_communicator(comm_nr)
      call mtx_allgather_real8(SIGMA_SPP,NREND-NRSTART+1,conduct_sp) 
      CALL mtx_reset_communicator

      RSUM=0.D0
      DO NR=NRSTART, NREND
         DO NSA=1, NSAMAX
            RSUM = RSUM + RJS(NR,NSA) + RJES(NR,NSA)
         END DO
         djdt_nind(NR)= (RSUM + RJ_BS(NR) - RJ_NIND_M(NR) )/DELT ! MA/(m2*s)
      END DO
!      WRITE(*,'(A,I4,5E14.6)') "RJ_NIND, ", NRSTART, RSUM, RJ_NIND_M(NRSTART), djdt_nind(NRSTART)

      DEPS=1.D0
      N_IMPL_E=0
      DO NR=1, NRMAX
         IF(E1(NR).eq.0.D0)THEN
            EM_M(NR)=1.D-12
         ELSE
            EM_M(NR)=E1(NR)
         END IF
      END DO
      DO WHILE(N_IMPL_E.le.10.and.DEPS.ge.1.D-12)
!      DO WHILE(N_IMPL_E.le.5)
         DEPS=0.D0
         N_IMPL_E = N_IMPL_E + 1
         CALL E_IND_IMPLICIT_CD(djdt_nind) ! E1(1:NRMAX) is updated
         DO NR=1, NRMAX
            DEPS = DEPS + (E1(NR)-EP_M(NR))**2/EM_M(NR)**2
            EP_M(NR)=E1(NR)
         END DO
         DO NR=NRSTART, NREND
            EP(NR)=E1(NR)
         END DO
         DEPS = DEPS / NRMAX
      END DO

!      DO NR=NRSTART, NREND
!         WRITE(*,'(A, 5E14.6)') "SIGMA_SPP, ", SIGMA_SPP(NR), conduct_sp(NR), EM(NR), EP(NR)
!      END DO

!      IF (NRANK.eq.0) &
!           WRITE(6,'(A,I3,E14.6)') "CALE_CONVERSION ON AXIS= ", N_IMPL_E, DEPS

!      IF(NRANK.eq.0)THEN
!         WRITE(*,'(A,99E14.6)') "TEST E_IND", djdt_nind(1)*1.D6, E1(1), sigma_spp(1)
!      END IF

      DO NR=NRSTART, NREND
         EP(NR)=E1(NR)
      END DO

      DO NR=1, NRMAX
         RJ_IND(NR)=conduct_sp(NR)*E1(NR)*1.D-6
      END DO

      DO NR=1, NRMAX
         CALL UPDATE_QLM(NR)
      END DO

      END SUBROUTINE AFTER_DO_WHILE_CD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      end MODULE fpcalj
