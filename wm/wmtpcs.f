C***********************************************************************
C
C     Multi-mode calculation for TOPICS
C
C***********************************************************************
c
      SUBROUTINE WM_TOPICS(IERR)
c
      include '../wm/wmcomm.inc'
      DIMENSION PABSR_SV(NRM,NSM,MNPHMX)
c
      NPH0_SV  = NPH0
      PRFIN_SV = PRFIN
C
      WRITE(6,*) "========== MULTI-MODE CALCULATION =========="
      DO NM = 1, MNPHMAX
         NPH0 = MNPH0(NM)
         PRFIN = PRFIN_SV * PFRAC(NM)
         WRITE(6,*)
         WRITE(6,*) "== Toroidal mode number : ",MNPH0(NM)
         WRITE(6,*) "== Power fraction       : ",PFRAC(NM)
C
         CALL WMEXEC(IERR)
         CALL MPSYNC
         IF(IERR.NE.0) EXIT
C
         DO NS = 1, NSMAX
            DO NR = 1, NRMAX
               PABSR_SV(NR,NS,NM) = PABSR(NR,NS)
            ENDDO
         ENDDO
      ENDDO
C
      IF(IERR.EQ.0.AND.MYRANK.EQ.0) THEN
         DO NS = 1, NSMAX
            DO NR = 1, NRMAX
               PABSRM = 0.D0
               DO NM = 1, MNPHMAX
                  PABSRM = PABSRM + PABSR_SV(NR,NS,NM)
               ENDDO
               PABSR(NR,NS) = PABSRM
            ENDDO
         ENDDO
         CALL WM_TOPICS_OUT
      ENDIF
C
      NPH0  = NPH0_SV
      PRFIN = PRFIN_SV
C
      WRITE(6,*) "======= MULTI-MODE CALCULATION NORMAL END ======="
      WRITE(6,*)
C
      return
      end
C
C***********************************************************************
C
C     Make Pabs(r,s) output file for TOPICS
C
C***********************************************************************
C
      SUBROUTINE WM_TOPICS_OUT
C
      INCLUDE '../wm/wmcomm.inc'
C
      CHARACTER KNAMWT*80
      DATA KNAMWT /'wm_topics.out'/
C
C     Output: XRHO(NR), PABSR(NR,NS)
C
      ntopics=21
      CALL FWOPEN(ntopics,KNAMWT,1,MODEFW,'WM',IERR)
      IF(IERR.NE.0) RETURN
C
      REWIND(ntopics)
      WRITE(ntopics,'(1P7E15.7)') (XRHO(NR),(PABSR(NR,NS),NS=1,6),
     &     NR=1,NRMAX)
      CLOSE(ntopics)
C
      RETURN
      END
