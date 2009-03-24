C***********************************************************************
C
C     Multi-mode calculation for TOPICS
C
C***********************************************************************
c
      SUBROUTINE WM_TOPICS(IERR)
c
      include '../wm/wmcomm.inc'
      DIMENSION PABSR_SV(NRM,NSM,MNPHMX), PABSTL(NSM), FACT(MNPHMX)
c
      NPH0_SV  = NPH0
      PRFIN_SV = PRFIN
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
         DO NS = 1, NSMAX
            PABSTL(NS) = 0.D0
            DO NR = 1, NRMAX
               PABSR_SV(NR,NS,NM) = PABSR(NR,NS)
               PABSTL(NS) = PABSTL(NS) + PABSR(NR,NS)
            ENDDO
         ENDDO
C
         PABSTTL = 0.D0
         DO NS = 1, NSMAX
            PABSTTL = PABSTTL + PABSTL(NS)
         ENDDO
C
         IF(PRFIN.GT.0.D0.AND.PABSTT.GT.0.D0) THEN
            FACT(NM)=PRFIN/PABSTTL
         ELSE
            FACT(NM)=1.D0
         ENDIF
      ENDDO
C
      IF(IERR.EQ.0.AND.MYRANK.EQ.0) THEN
         PABSTTM = 0.D0
         DO NS = 1, NSMAX
            PABSTL(NS) = 0.D0
            DO NR = 1, NRMAX
               PABSRM = 0.D0
               DO NM = 1, MNPHMAX
                  PABSRM = PABSRM + FACT(NM) * PABSR_SV(NR,NS,NM)
               ENDDO
               PABSR(NR,NS) = PABSRM
               PABSTL(NS) = PABSTL(NS) + PABSRM
            ENDDO
            PABSTTM = PABSTTM + PABSTL(NS)
         ENDDO
         CALL WM_TOPICS_OUT
      ENDIF
C
      NPH0  = NPH0_SV
      PRFIN = PRFIN_SV
C
      WRITE(6,*)
      WRITE(6,'(A,1PE12.4)')  "TOT. PABS=",PABSTTM
      WRITE(6,'(A,1P6E12.4)') "PABS(NS) =",(PABSTL(NS),NS=1,NSMAX)
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
