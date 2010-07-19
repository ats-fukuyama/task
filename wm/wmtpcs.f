C***********************************************************************
C
C     Multi-mode calculation for TOPICS
C
C***********************************************************************
c
      SUBROUTINE WM_TOPICS(IERR)
c
      include '../wm/wmcomm.inc'
      DIMENSION PABSR_SV(NRM,NSM),PABST_SV(NSM)
      DIMENSION PABS_SV(MDM,NDM,NRM,NSM)
c
      NPH0_SV  = NPH0
      PRFIN_SV = PRFIN
C
      DO NS = 1, NSMAX
         PABST_SV(NS)=0.D0
      DO NR = 1, NRMAX
         PABSR_SV(NR,NS) = 0.D0
      DO NPH=1,NDSIZ
      DO NTH=1,MDSIZ
         PABS_SV(NTH,NPH,NR,NS) = 0.0
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      WRITE(6,*) "========= MULTI-MODE CALCULATION START ========="
      DO NM = 1, MNPHMAX
         NPH0 = MNPH0(NM)
         PRFIN = PRFIN_SV * PFRAC(NM)
         WRITE(6,*)
         WRITE(6,'(A,I4)') "== Toroidal mode number : ",MNPH0(NM)
         WRITE(6,'(A,F8.6)') "== Power fraction       : ",PFRAC(NM)
C
         CALL WMEXEC(IERR)
         CALL MPSYNC
         IF(IERR.NE.0) EXIT
C
         PABSTT_SV=PABSTT_SV+PABSTT
         DO NS=1,NSMAX
            PABST_SV(NS)=PABST_SV(NS)+PABST(NS)
         DO NR = 1, NRMAX
            PABSR_SV(NR,NS) = PABSR_SV(NR,NS) + PABSR(NR,NS)
         DO NPH=1,NDSIZ
         DO NTH=1,MDSIZ
            PABS_SV(NTH,NPH,NR,NS) = PABS_SV(NTH,NPH,NR,NS) 
     &          + PABS(NTH,NPH,NR,NS)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
C
      ENDDO
C
      IF(IERR.EQ.0.AND.MYRANK.EQ.0) THEN
         PABSTT=PABSTT_SV
         DO NS=1,NSMAX
            PABST(NS)=PABST_SV(NS)
         DO NR = 1, NRMAX
            PABSR(NR,NS) = PABSR_SV(NR,NS)
         DO NPH=1,NDSIZ
         DO NTH=1,MDSIZ
            PABS(NTH,NPH,NR,NS) = PABS_SV(NTH,NPH,NR,NS) 
         ENDDO
         ENDDO
         ENDDO
         ENDDO

         CALL WM_TOPICS_OUT

      ENDIF
C
      NPH0  = NPH0_SV
      PRFIN = PRFIN_SV
C
      IF(IERR.NE.0) THEN
         WRITE(6,*)
     &        "======= MULTI-MODE CALCULATION ABNORMAL END ======="
         RETURN
      ENDIF
C
      WRITE(6,*)
      WRITE(6,'(A,1PE12.4)')  "TOT. PABS=",PABSTT
      WRITE(6,'(A,1P6E12.4)') "PABS(NS) =",(PABST(NS),NS=1,NSMAX)
      WRITE(6,*)
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
