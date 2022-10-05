!

! ***********************************************************
!       Nuclear Fusion Reaction Rate with Legendre Expansion
! ***********************************************************

      MODULE fpnflg

      USE fpcomm
      USE fpinit
      USE libspf, ONLY: dpleg

      PUBLIC NF_LG_FUNCTION
      PUBLIC NF_REACTION_RATE_LG

      PRIVATE

      contains

      SUBROUTINE NF_LG_FUNCTION

      USE fpcomm

      INTEGER::NTH, IERR, L
      REAL(rkind),DIMENSION(0:LLMAX_NF):: PL

      !     ---- legendre polynomials
      DO NTH=1,NTHMAX
         CALL DPLEG(COSM(NTH),LLMAX_NF,PL,IERR)
         DO L=0,LLMAX_NF
            PL_NF(L,NTH)=PL(L)
         END DO
      END DO

      END SUBROUTINE NF_LG_FUNCTION

      SUBROUTINE NF_REACTION_RATE_LG(NR,ID)

      USE fpcomm
      USE libmpi
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR,ID
      INTEGER:: L, K, L2, NSB1, NSB2, NSA1, NSA2, NP1, NP2, NTH, VLOC
      REAL(rkind):: RSUM, FACT, RSUM2, FACT1, FACT2
      REAL(rkind):: RSUM_SUM, double_count, RSUM_B
      REAL(rkind),DIMENSION(0:LLMAX_NF,NPSTART:NPEND):: F_LG1, F_LG2, R_1
      REAL(rkind),DIMENSION(0:LLMAX_NF,NPMAX):: F_LG2_ALL, R_2
      REAL(rkind),DIMENSION(NTHMAX,NPMAX):: FNSB_B2_VOLP
      REAL(rkind),DIMENSION(NTHMAX,NPSTART:NPEND):: FNSB_temp

      NSB1=NSB1_NF(ID)
      NSB2=NSB2_NF(ID)
      NSA1=NSA1_NF(ID)
      NSA2=NSA2_NF(ID)

      double_count=1.D0
      IF(ID.eq.1.or.ID.eq.2.or.ID.eq.5) double_count=0.5D0

      FACT = RNFD0(NSB1)*RNFD0(NSB2)*1.D20*double_count

!      Decomposing the distribution function into Legendre harmonics

      DO NP1=NPSTART,NPEND
         FACT1=2.D0*PI*DELP(NSB1)*PM(NP1,NSB1)**2
         DO L=0,LLMAX_NF
            F_LG1(L,NP1)=0.D0
            DO NTH=1,NTHMAX
               F_LG1(L,NP1)=F_LG1(L,NP1) &
                    +FNSB(NTH,NP1,NR,NSB1) &
                    *PL_NF(L,NTH)*SINM(NTH)*DELTH
            END DO
            F_LG1(L,NP1)=FACT1*F_LG1(L,NP1)
         END DO
      END DO

      DO NP1=NPSTART,NPEND
         FACT2=2.D0*PI*DELP(NSB2)*PM(NP1,NSB2)**2
         DO L=0,LLMAX_NF
            F_LG2(L,NP1)=0.D0
            DO NTH=1,NTHMAX
               F_LG2(L,NP1)=F_LG2(L,NP1) &
                    +FNSB(NTH,NP1,NR,NSB2) &
                    *PL_NF(L,NTH)*SINM(NTH)*DELTH
            END DO
            F_LG2(L,NP1)=FACT2*F_LG2(L,NP1)
         END DO
      END DO
      CALL mtx_set_communicator(comm_np)
      CALL mtx_allgather_real8(F_LG2,(LLMAX_NF+1)*(NPEND-NPSTART+1),F_LG2_ALL)

      DO NP1=NPSTART,NPEND
         DO NTH=1,NTHMAX
            RATE_NF_D1(NTH,NP1,NR,ID)=0.D0
            RATE_NF_D2(NTH,NP1,NR,ID)=0.D0
         END DO
         DO L=0,LLMAX_NF
            R_1(L,NP1)=0.D0
         END DO
      END DO

      DO NP2=1,NPMAX
         DO L=0,LLMAX_NF
            R_2(L,NP2)=0.D0
         END DO
      END DO
!      Calculating the fusion reaction rate

      RSUM=0.D0
      RSUM2=0.D0

      DO NP2=1,NPMAX
         DO K=0,LLMAX_NF
            DO NP1=NPSTART,NPEND
               DO L=0,LLMAX_NF
                  R_1(L,NP1)=R_1(L,NP1)+FACT*F_LG2_ALL(K,NP2) &
                       *SIGMAV_LG(L,NP1,K,NP2,ID)

                  IF(PM(NP1,NSB1).ge.pmax_bb(NSB1).and.PM(NP2,NSB2).ge.pmax_bb(NSB2))THEN
                     RSUM2=RSUM2+FACT*F_LG1(L,NP1)*F_LG2_ALL(K,NP2) &
                          *SIGMAV_LG(L,NP1,K,NP2,ID)
                  END IF
               END DO
            END DO
         END DO
      END DO

      DO NP2=1,NPMAX
         DO K=0,LLMAX_NF
            RSUM_B=0.D0
            DO NP1=NPSTART,NPEND
               DO L=0,LLMAX_NF
                  RSUM_B=RSUM_B+FACT*F_LG1(L,NP1) &
                       *SIGMAV_LG(L,NP1,K,NP2,ID)
               END DO
            END DO
            CALL mtx_allreduce1_real8(RSUM_B,3,RSUM_SUM,VLOC)
            R_2(K,NP2)=RSUM_SUM
         END DO
      END DO

      DO NP1=NPSTART,NPEND
         DO L=0,LLMAX_NF
            RSUM=RSUM+F_LG1(L,NP1)*R_1(L,NP1)
         END DO
      END DO

      CALL mtx_allreduce1_real8(RSUM,3,RSUM_SUM,VLOC) ! RATE_NF
      RSUM=RSUM_SUM
      CALL mtx_allreduce1_real8(RSUM2,3,RSUM_SUM,VLOC) ! integrate np1
      RSUM2=RSUM_SUM

      DO NP1=NPSTART,NPEND
         DO NTH=1,NTHMAX
            RSUM_B=0.D0
            DO L=0,LLMAX_NF
               RATE_NF_D1(NTH,NP1,NR,ID)= &
                    RATE_NF_D1(NTH,NP1,NR,ID)+R_1(L,NP1)*PL_NF(L,NTH)*FNSB(NTH,NP1,NR,NSB1)
            END DO
         END DO
      END DO

!     comm fnsb for nsb2
      DO NP1=NPSTART, NPEND
         DO NTH=1,NTHMAX
            FNSB_temp(NTH,NP1)=FNSB(NTH,NP1,NR,NSB2)
         END DO
      END DO
      CALL mtx_set_communicator(comm_np)
      CALL mtx_allgather_real8(FNSB_temp,NTHMAX*(NPEND-NPSTART+1),FNSB_B2_VOLP)


      DO NP2=1,NPMAX ! not devide
         DO NTH=1,NTHMAX
            RSUM_B=0.D0
            DO L=0,LLMAX_NF
               RSUM_B=RSUM_B+R_2(L,NP2)*PL_NF(L,NTH)*FNSB_B2_VOLP(NTH,NP2)
            END DO
            RATE_NF_D2(NTH,NP2,NR,ID)=RSUM_B
         END DO
      END DO

      CALL mtx_reset_communicator

      IF((N_IMPL.eq.0.or.N_IMPL.gt.LMAXFP).and.NR.eq.1.and.NPSTART.eq.1)THEN
         WRITE(6,'(A,3I4)') '|-NF_REACTION_RATE:', nrank, comm_nr%rank, comm_nsa%rank
         WRITE(6,'(A,I5,A,2I5,A,2I5)') '   |-ID,NSB1,NSB2 -> NSA1,NSA2=' &
              ,ID,':  ',NSB1,NSB2,' -> ',NSA1,NSA2
         WRITE(6, *) "  |-ID,  NR,NSB1,NSB2,  <sigma*v>,      ENG1_NF,      RATE_NF,   RATE_beam-beam"
         WRITE(6,'("  ",4I5,1PE14.6,1PE12.4, 1P2E14.6)') ID,NR,NSB1,NSB2,RSUM/FACT,ENG1_NF(ID), RSUM, RSUM2
      END IF

      RATE_NF(NR,ID)=RSUM

      END SUBROUTINE NF_REACTION_RATE_LG

      END MODULE fpnflg
