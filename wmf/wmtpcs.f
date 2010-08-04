c$$$C***********************************************************************
c$$$C
c$$$C     Multi-mode calculation for TOPICS
c$$$C
c$$$C***********************************************************************
c$$$c
c$$$      SUBROUTINE WM_TOPICS(IERR)
c$$$c
c$$$      include 'wmcomm.inc'
c$$$      DIMENSION PABSR_SV(NRM,NSM,NPHSM), PABSTL(NSM), FACT(NPHSM)
c$$$      DIMENSION GPABS_SV(MDM,NDM,NRM,NSM)
c$$$c
c$$$      NPH0_SV  = NPH0
c$$$      PRFIN_SV = PRFIN
c$$$C
c$$$      DO NS=1,NSMAX
c$$$      DO NR=1,NRMAX+1
c$$$      DO NPH=1,NDSIZ
c$$$      DO NTH=1,MDSIZ
c$$$         GPABS_SV(NTH,NPH,NR,NS) = 0.0
c$$$      ENDDO
c$$$      ENDDO
c$$$      ENDDO
c$$$      ENDDO
c$$$C
c$$$      WRITE(6,*) "========= MULTI-MODE CALCULATION START ========="
c$$$      DO NPHS = 1, NPHSMAX
c$$$         NPH0 = NPH0S(NPHS)
c$$$         PRFIN = PRFIN_SV * PFRACS(NPHS)
c$$$         WRITE(6,*)
c$$$         WRITE(6,'(A,I4)') "== Toroidal mode number : ",NPH0S(NPHS)
c$$$         WRITE(6,'(A,F8.6)') "== Power fraction       : ",PFRACS(NPHS)
c$$$C
c$$$         CALL WMEXEC(IERR)
c$$$         CALL MPSYNC
c$$$         IF(IERR.NE.0) EXIT
c$$$C
c$$$         IF(MYRANK.EQ.0) THEN
c$$$            CALL WMPOUT(NPH0)
c$$$            IF(MODELW.EQ.1) CALL WMDOUT(IERR)
c$$$         ENDIF
c$$$C
c$$$         DO NS = 1, NSMAX
c$$$            PABSTL(NS) = 0.D0
c$$$            DO NR = 1, NRMAX
c$$$               PABSR_SV(NR,NS,NPHS) = PABSR(NR,NS)
c$$$               PABSTL(NS) = PABSTL(NS) + PABSR(NR,NS)
c$$$            ENDDO
c$$$         ENDDO
c$$$C
c$$$         PABSTTL = 0.D0
c$$$         DO NS = 1, NSMAX
c$$$            PABSTTL = PABSTTL + PABSTL(NS)
c$$$         ENDDO
c$$$C
c$$$         IF(PRFIN.GT.0.D0.AND.PABSTT.GT.0.D0) THEN
c$$$            FACT(NM)=PRFIN/PABSTTL
c$$$         ELSE
c$$$            FACT(NM)=1.D0
c$$$         ENDIF
c$$$C
c$$$C     === For Graphics: PP* ===
c$$$C
c$$$         DO NS=1,NSMAX
c$$$         DO NR=1,NRMAX+1
c$$$         DO NPH=1,NDSIZ
c$$$         DO NTH=1,MDSIZ
c$$$            GPABS_SV(NTH,NPH,NR,NS) = GPABS_SV(NTH,NPH,NR,NS) 
c$$$     &          + REAL(FACT(NM)) * GUCLIP(PABS(NTH,NPH,NR,NS))
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$C
c$$$      ENDDO
c$$$C
c$$$      IF(IERR.EQ.0.AND.MYRANK.EQ.0) THEN
c$$$         PABSTTM = 0.D0
c$$$         DO NS = 1, NSMAX
c$$$            PABST(NS) = 0.D0
c$$$            DO NR = 1, NRMAX
c$$$               PABSRM = 0.D0
c$$$               DO NM = 1, MNPHMAX
c$$$                  PABSRM = PABSRM + FACT(NM) * PABSR_SV(NR,NS,NM)
c$$$               ENDDO
c$$$               PABSR(NR,NS) = PABSRM
c$$$               PABST(NS) = PABST(NS) + PABSRM
c$$$            ENDDO
c$$$            PABSTTM = PABSTTM + PABST(NS)
c$$$         ENDDO
c$$$         CALL WM_TOPICS_OUT
c$$$C
c$$$C     === For Graphics: PP* ===
c$$$C
c$$$         DO NS=1,NSMAX
c$$$         DO NR=1,NRMAX+1
c$$$         DO NPH=1,NDSIZ
c$$$         DO NTH=1,MDSIZ
c$$$            PABS(NTH,NPH,NR,NS) = DBLE(GPABS_SV(NTH,NPH,NR,NS))
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO
c$$$         ENDDO         
c$$$      ENDIF
c$$$C
c$$$      NPH0  = NPH0_SV
c$$$      PRFIN = PRFIN_SV
c$$$C
c$$$      IF(IERR.NE.0) THEN
c$$$         WRITE(6,*)
c$$$     &        "======= MULTI-MODE CALCULATION ABNORMAL END ======="
c$$$         RETURN
c$$$      ENDIF
c$$$C
c$$$      WRITE(6,*)
c$$$      WRITE(6,'(A,1PE12.4)')  "TOT. PABS=",PABSTTM
c$$$      WRITE(6,'(A,1P6E12.4)') "PABS(NS) =",(PABST(NS),NS=1,NSMAX)
c$$$      WRITE(6,*)
c$$$C
c$$$      WRITE(6,*) "======= MULTI-MODE CALCULATION NORMAL END ======="
c$$$      WRITE(6,*)
c$$$C
c$$$      return
c$$$      end
C
C***********************************************************************
C
C     Make P_abs(r,s)/J_CD(r) output file for TOPICS/ACCOME
C
C***********************************************************************
C
      SUBROUTINE WM_TOPICS_OUT(PABSTS,IERR)
C
      INCLUDE 'wmcomm.inc'
C
      real(8), dimension(NSMAX,NPHSMAX), intent(in) :: PABSTS
      integer(4), intent(out) :: IERR
C
      CHARACTER KNAMWT(1:2)*80
      DATA KNAMWT /'wm_pwrcur.out','wm_antenna.out'/
C
C     Output: XRHO(NR), PABSR(NR,NS), PCURR(NR)
C             8 columns
C
      nout=21
      CALL FWOPEN(nout,KNAMWT(1),1,MODEFW,'WM',IERR)
      IF(IERR.NE.0) RETURN
C
      REWIND(nout)
      WRITE(nout,'(A1,2(2X,A))') '# units:','PABSR(NS) [W/m^3],',
     &     'PCURR [A/m^2]'
      WRITE(nout,'(A1,6X,A3,11X,A4,10X,6(A8,7X),A5)') '#','RHO','PSIN',
     &     'PABSR(1)','PABSR(2)','PABSR(3)','PABSR(4)',
     &     'PABSR(5)','PABSR(6)','PCURR'
      WRITE(nout,'(1P9E15.7)') (XRHO(NR),FNPSIN(XRHO(NR)),
     &     (PABSR(NR,NS),NS=1,6),
     &     PCURR(NR),NR=1,NRMAX)
      CLOSE(nout)
C
C     Power absorption ratio per toroidal mode
C     Output: NPH0S(NPHS), Power fraction(NPHS), PABSTS(NS,NPHS)
C             max 8 columns
C
      nout=21
      CALL FWOPEN(nout,KNAMWT(2),1,MODEFW,'WM',IERR)
      IF(IERR.NE.0) RETURN
C
      REWIND(nout)
      WRITE(nout,'(3A)') '# NPH','    PFRACS ', '   PABSTS(NS) [W]'
      DO NPHS=1,NPHSMAX
         WRITE(nout,'(I5,1P6E12.4)') 
     &        NPH0S(NPHS),
     &        SUM(PABSTS(1:MIN(NSMAX,6),NPHS))/PABSTT,
     &        (PABSTS(NS,NPHS),NS=1,MIN(NSMAX,6))
      ENDDO
      WRITE(nout,'(A)') 
     &     '--------------------------------------------'//
     &     '---------------------------------'
      WRITE(nout,'(17X,1P5E12.4)') (PABST(NS),NS=1,MIN(NSMAX,6))
      WRITE(nout,'(A8,1P1E11.4)') 'TOTAL = ',PABSTT
      CLOSE(nout)
C
      RETURN
      END
