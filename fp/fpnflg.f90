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
      REAL(8),DIMENSION(0:LLMAX_NF):: PL

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
      REAL(8):: RSUM, FACT, RSUM2, FACT1, FACT2, FACT3, FACT4, FACT5
      REAL(8):: RSUM_SUM, double_count
      REAL(8),DIMENSION(0:LLMAX_NF,NPSTART:NPEND):: F_LG1, F_LG2, R1, R2
      REAL(8),DIMENSION(0:LLMAX_NF,NPMAX):: F_LG1_ALL, F_LG2_ALL, R1_SUM, R2_SUM

      NSB1=NSB1_NF(ID)
      NSB2=NSB2_NF(ID)
      NSA1=NSA1_NF(ID)
      NSA2=NSA2_NF(ID)

      double_count=1.D0
      IF(ID.eq.1.or.ID.eq.2.or.ID.eq.5) double_count=0.5D0

      FACT = RNFD0(NSB1)*RNFD0(NSB2)*1.D20*double_count

!      Decomposing the distribution function into Legendre harmonics

      DO NP1=NPSTART,NPEND
         DO L=0,LLMAX_NF
            F_LG1(L,NP1)=0.D0
            F_LG2(L,NP1)=0.D0
            DO NTH=1,NTHMAX
               F_LG1(L,NP1)=F_LG1(L,NP1) &
                    +FNSB(NTH,NP1,NR,NSB1) &
                    *PL_NF(L,NTH)*SINM(NTH)*DELTH
               F_LG2(L,NP1)=F_LG2(L,NP1) &
                    +FNSB(NTH,NP1,NR,NSB2) &
                    *PL_NF(L,NTH)*SINM(NTH)*DELTH
            END DO
            F_LG1(L,NP1)=0.5D0*(2.D0*L+1.D0)*F_LG1(L,NP1)
            F_LG2(L,NP1)=0.5D0*(2.D0*L+1.D0)*F_LG2(L,NP1)
         END DO
      END DO
      CALL mtx_set_communicator(comm_np)
      CALL mtx_allgather_real8(F_LG1,(LLMAX_NF+1)*(NPEND-NPSTART+1),F_LG1_ALL)
      CALL mtx_allgather_real8(F_LG2,(LLMAX_NF+1)*(NPEND-NPSTART+1),F_LG2_ALL)

      DO NP1=NPSTART,NPEND
         DO NTH=1,NTHMAX
            RATE_NF_D1(NTH,NP1,NR,ID)=0.D0
            RATE_NF_D2(NTH,NP1,NR,ID)=0.D0
         END DO
         DO L=0,LLMAX_NF
            R1(L,NP1)=0.D0
            R2(L,NP1)=0.D0
         END DO
      END DO

!      Calculating the fusion reaction rate

      RSUM=0.D0

      DO L=0,LLMAX_NF
         FACT1=2.D0/(2.D0*L+1)
         DO K=0,LLMAX_NF
            FACT2=2.D0/(2.D0*K+1)
            DO NP1=NPSTART,NPEND
               FACT3=2.D0*PI*DELP(NSB1)*PM(NP1,NSB1)**2
               DO NP2=1,NPMAX
                  FACT4=2.D0*PI*DELP(NSB2)*PM(NP2,NSB2)**2
                  FACT5=2.D0*PI*DELP(NSB1)*PM(NP2,NSB1)**2
                  RSUM=RSUM+FACT*FACT1*FACT2*FACT3*FACT4 &
                       *F_LG1_ALL(L,NP1)*F_LG2_ALL(K,NP2) &
                       *SIGMAV_LG(L,K,NP1,NP2,ID)

                  R1(L,NP1)=R1(L,NP1)+FACT*FACT2*FACT4 &
                       *SIGMAV_LG(L,K,NP1,NP2,ID) &
                       *F_LG2_ALL(K,NP2)

                  R2(K,NP1)=R2(K,NP1)+FACT*FACT1*FACT5 &
                       *SIGMAV_LG(L,K,NP2,NP1,ID) &
                       *F_LG1_ALL(L,NP2)
               END DO
            END DO
         END DO
      END DO

      CALL mtx_allreduce1_real8(RSUM,3,RSUM_SUM,VLOC)
      RSUM=RSUM_SUM

      DO NTH=1,NTHMAX
         DO NP1=NPSTART,NPEND
            DO L=0,LLMAX_NF
               RATE_NF_D1(NTH,NP1,NR,ID)= &
                    RATE_NF_D1(NTH,NP1,NR,ID) &
                    +R1(NP1,L)*PL_NF(L,NTH)*FNSB(NTH,NP1,NR,NSB1)
               RATE_NF_D2(NTH,NP1,NR,ID)= &
                    RATE_NF_D2(NTH,NP1,NR,ID) &
                    +R2(NP1,L)*PL_NF(L,NTH)*FNSB(NTH,NP1,NR,NSB2)
            END DO
         END DO
      END DO

      CALL mtx_reset_communicator

      IF((N_IMPL.eq.0.or.N_IMPL.gt.LMAXFP).and.NR.eq.1.and.NPSTART.eq.1)THEN
         WRITE(6,'(A,3I4)') '|-NF_REACTION_RATE:', nrank, comm_nr%rank, comm_nsa%rank
         WRITE(6,'(A,I5,A,2I5,A,2I5)') '   |-ID,NSB1,NSB2 -> NSA1,NSA2=' &
              ,ID,':  ',NSB1,NSB2,' -> ',NSA1,NSA2
         WRITE(6, *) "  |-ID,  NR,NSB1,NSB2,  <sigma*v>,      ENG1_NF,      RATE_NF"
         WRITE(6,'("  ",4I5,1PE14.6,1PE12.4, 1P2E14.6)') ID,NR,NSB1,NSB2,RSUM/FACT,ENG1_NF(ID), RSUM
      END IF

      RATE_NF(NR,ID)=RSUM

      WRITE(6,*) 'RATE_NF', RATE_NF(NR,ID), 'ID', ID

      END SUBROUTINE NF_REACTION_RATE_LG

      END MODULE fpnflg
