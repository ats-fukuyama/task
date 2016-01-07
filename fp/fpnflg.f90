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
      INTEGER:: L, K, L2, NSB1, NSB2, NSA1, NSA2, NP1, NP2, NTH
      REAL(8):: RSUM, FACT, RSUM2, FACT1, FACT2, FACT3, FACT4
      REAL(8),DIMENSION(NPMAX,0:LLMAX_NF):: R1, R2       
      REAL(8),dimension(0:LLMAX_NF,NPMAX,NSAMAX):: F_LG1, F_LG2

      NSB1=NSB1_NF(ID)
      NSB2=NSB2_NF(ID)
      NSA1=NSA1_NF(ID)
      NSA2=NSA2_NF(ID)

!      Decomposing the distribution function into Legendre harmonics

      DO NP1=1,NPMAX
         DO L=0,LLMAX_NF
            F_LG1(L,NP1,NSB1)=0.D0
            F_LG2(L,NP1,NSB2)=0.D0
            DO NTH=1,NTHMAX
               F_LG1(L,NP1,NSB1)=F_LG1(L,NP1,NSB1) &
                    +FNSB(NTH,NP1,NR,NSB1) &
                    *PL_NF(L,NTH)*SINM(NTH)*DELTH
               F_LG2(L,NP1,NSB2)=F_LG2(L,NP1,NSB2) &
                    +FNSB(NTH,NP1,NR,NSB2) &
                    *PL_NF(L,NTH)*SINM(NTH)*DELTH
            END DO
            F_LG1(L,NP1,NSB1)=0.5D0*(2.D0*L+1.D0)*F_LG1(L,NP1,NSB1)
            F_LG2(L,NP1,NSB2)=0.5D0*(2.D0*L+1.D0)*F_LG2(L,NP1,NSB2)
         END DO
      END DO

      FACT = RNFD0(NSB1)*RNFD0(NSB2)*1.D20

      DO NP1=1,NPMAX
         DO NTH=1,NTHMAX
            RATE_NF_D1(NTH,NP1,NR,ID)=0.D0
            RATE_NF_D2(NTH,NP1,NR,ID)=0.D0
         END DO
         DO L=0,LLMAX_NF
            R1(NP1,L)=0.D0
            R2(NP1,L)=0.D0
         END DO
      END DO

!      Calculating the fusion reaction rate

      RSUM=0.D0

      DO L=0,LLMAX_NF
         FACT1=2.D0/(2.D0*L+1)
         DO K=0,LLMAX_NF
            FACT2=2.D0/(2.D0*K+1)
            DO NP1=1,NPMAX
               FACT3=2.D0*PI*DELP(NSB1)*PM(NP1,NSB1)**2
               DO NP2=1,NPMAX
                  FACT4=2.D0*PI*DELP(NSB2)*PM(NP2,NSB2)**2
                  RSUM=RSUM+FACT*FACT1*FACT2*FACT3*FACT4 &
                       *F_LG1(L,NP1,NSB1)*F_LG2(K,NP2,NSB2) &
                       *SIGMAV_LG(L,K,NP1,NP2,ID)

                  R1(NP1,L)=R1(NP1,L)+FACT*FACT2*FACT4 &
                       *SIGMAV_LG(L,K,NP1,NP2,ID) &
                       *F_LG2(K,NP2,NSB2)

                  R2(NP2,K)=R2(NP2,K)+FACT*FACT1*FACT3 &
                       *SIGMAV_LG(L,K,NP1,NP2,ID) &
                       *F_LG1(L,NP1,NSB1)
               END DO
            END DO
         END DO
      END DO

      DO NTH=1,NTHMAX
         DO NP1=1,NPMAX
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

      RATE_NF(NR,ID)=RSUM
      IF((N_IMPL.eq.0.or.N_IMPL.gt.LMAXFP).and.NR.eq.1)THEN
         WRITE(6,*) '|-NF_REACTION_RATE:', nrank, comm_nr%rank
         WRITE(6,'(A,3I5,A,2I5)') '   |-ID,NSB1,NSB2 -> NSA1,NSA2=' &
              ,ID,NSB1,NSB2,' -> ',NSA1,NSA2
         WRITE(6, *) "  |-ID,  NR,NSB1,NSB2,  <sigma*v>,      ENG1_NF"
         WRITE(6,'("  ",4I5,1PE14.6,1PE12.4)') ID,NR,NSB1,NSB2,RATE_NF(NR,ID)/FACT,ENG1_NF(ID)
      END IF

      WRITE(6,*) 'RATE_NF', RATE_NF(NR,ID), 'ID', ID

      END SUBROUTINE NF_REACTION_RATE_LG

      END MODULE fpnflg
