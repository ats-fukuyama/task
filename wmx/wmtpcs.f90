! wmtpcs.f90

!***********************************************************************

!     Multi-mode calculation for TOPICS

!***********************************************************************

MODULE wmtpcs
  
  PRIVATE
  PUBLIC wm_topics

CONTAINS

  SUBROUTINE wm_topics_TOPICS(IERR)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    REAL(rkind),ALLOCATABLE:: PABSR_SV(:,:),PABST_SV(:),PABS_SV(:,:,:,:)
    INTEGER:: NPH0_SV,NS,NR,NHH,NTH

    ALLOCATE(PABSR_SV(nrmax+1,nsmax),PABST_SV(nsmax))
    ALLOCATE(PABS_SV(nthmax,nhhmax,nrmax+1,nsmax)

    NPH0_SV  = NPH0

    DO NS = 1, NSMAX
       PABST_SV(NS)=0.D0
       DO NR = 1, NRMAX
          PABSR_SV(NR,NS) = 0.D0
          DO NHH=1,NDSIZ
             DO NTH=1,MDSIZ
                PABS_SV(NTH,NHH,NR,NS) = 0.0
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    WRITE(6,*) "========= MULTI-MODE CALCULATION START ========="
    DO NPH = 1, NPHMAX
       NPH0 = NPH0_LOOP(NPH)
       PRFIN = PRFIN_SV * PFRACL(NPH)
       WRITE(6,*)
       WRITE(6,'(A,I4)') "== Toroidal mode number : ",NPH0_LOOP(NPH)
       
         CALL WMEXEC(IERR)
         CALL mtx_barrier
         IF(IERR.NE.0) EXIT

         PABSTT_SV=PABSTT_SV+PABSTT
         DO NS=1,NSMAX
            PABST_SV(NS)=PABST_SV(NS)+PABST(NS)
         DO NR = 1, NRMAX
            PABSR_SV(NR,NS) = PABSR_SV(NR,NS) + PABSR(NR,NS)
         DO NHH=1,NDSIZ
         DO NTH=1,MDSIZ
            PABS_SV(NTH,NHH,NR,NS) = PABS_SV(NTH,NHH,NR,NS) 
     &          + PABS(NTH,NHH,NR,NS)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

      ENDDO

      IF(IERR.EQ.0.AND.NRANK.EQ.0) THEN
         PABSTT=PABSTT_SV
         DO NS=1,NSMAX
            PABST(NS)=PABST_SV(NS)
         DO NR = 1, NRMAX
            PABSR(NR,NS) = PABSR_SV(NR,NS)
         DO NHH=1,NDSIZ
         DO NTH=1,MDSIZ
            PABS(NTH,NHH,NR,NS) = PABS_SV(NTH,NHH,NR,NS) 
         ENDDO
         ENDDO
         ENDDO
         ENDDO

         CALL WM_TOPICS_OUT

      ENDIF

      NPH0  = NPH0_SV

      IF(IERR.NE.0) THEN
         WRITE(6,*)
     &        "======= MULTI-MODE CALCULATION ABNORMAL END ======="
         RETURN
      ENDIF

      WRITE(6,*)
      WRITE(6,'(A,1PE12.4)')  "TOT. PABS=",PABSTT
      WRITE(6,'(A,1P6E12.4)') "PABS(NS) =",(PABST(NS),NS=1,NSMAX)
      WRITE(6,*)

      WRITE(6,*) "======= MULTI-MODE CALCULATION NORMAL END ======="
      WRITE(6,*)

      return
      end

!***********************************************************************

!     Make Pabs(r,s) output file for TOPICS

!***********************************************************************

      SUBROUTINE WM_TOPICS_OUT

      INCLUDE './wmcomm.inc'

      CHARACTER KNAMWT*80
      DATA KNAMWT /'wm_topics.out'/

!     Output: XRHO(NR), PABSR(NR,NS)

      ntopics=21
      CALL FWOPEN(ntopics,KNAMWT,1,MODEFW,'WM',IERR)
      IF(IERR.NE.0) RETURN

      REWIND(ntopics)
      WRITE(ntopics,'(1P7E15.7)') (XRHO(NR),(PABSR(NR,NS),NS=1,6),
     &     NR=1,NRMAX)
      CLOSE(ntopics)

      RETURN
      END
