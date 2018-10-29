C***********************************************************************
C
C     Multi-mode calculation
C
C***********************************************************************
c
      SUBROUTINE WM_LOOP(IERR)
c
      include 'wmcomm.inc'
      DIMENSION PABSR_SV(NRM,NSM),PABST_SV(NSM)
      CHARACTER KNHC*2, KNPH0*4
c
      IERR=0
      MODELEG=0
      MODELK=1
C     MODELK=0

      CALL WMSETG(IERR)
      IF(IERR.NE.0) RETURN

      DO NNR=1,NRMAX
         IF(XRHO(NNR)>1.d0) EXIT
      ENDDO
      NR_S=NNR-1
      IF(nrank.EQ.0) WRITE(6,'(A,I5)') 'NR_S=',NR_S
      
      CALL DPCHEK(NTHMAX,NRMAX+1,XRHO(1),XRHO(NRMAX+1),RR,IERR)
      IF(IERR.NE.0) RETURN
      CALL WMSOLV_PREP

      NPH0_SV  = NPH0
      NPHMAX_SV = NPHMAX
      NHHMAX_SV = NHHMAX

      IF(NPHMAX.EQ.1) THEN
         CALL WMEXEC(IERR)
         RETURN
      ENDIF

      IF(NPHMAX.eq.NHHMAX .and. NHC .gt. 1) THEN
         NHHMAX=NPHMAX/NHC
         NPHMAX=NHC
      END IF
C
C     Axisymmetric case : NHHMAX=1
C
      DO NS = 1, NSMAX
         PABST_SV(NS)=0.D0
      DO NR = 1, NRMAX
         PABSR_SV(NR,NS) = 0.D0
      ENDDO
      ENDDO

      DO NR=1,NRMAX+1
         DO NPH=1,NPHMAX
            DO NTH=1,NTHMAX
               CEFLD3D(1,NTH,NPH,NR)=(0.D0,0.D0)
               CEFLD3D(2,NTH,NPH,NR)=(0.D0,0.D0)
               CEFLD3D(3,NTH,NPH,NR)=(0.D0,0.D0)
            END DO
         END DO
      END DO
      DO NS=1,NSMAX
         DO NR=1,NRMAX+1
            DO NHH=1,NHHMAX !NPHMAX_SV
               DO NTH=1,NTHMAX
                  CPABS3D(NTH,NHH,NR,NS)=0.D0
                   PABS3D(NTH,NHH,NR,NS)=0.D0
               END DO
            END DO
         END DO
      END DO
      
      CALL WMPOUT_INIT
      DO NPH = 1,NPHMAX
         NPH0 = NPH-NPHMAX/2-1 + NPH0_SV
         WRITE(6,'(A,I4)') "== Toroidal mode number : ",NPH0
C
         CALL WMEXEC(IERR)
         CALL mtx_barrier
         IF(IERR.NE.0) EXIT

         DO NS=1,NSMAX
            DO NR=1,NRMAX+1
               DO NHH=1,NHHMAX
                  NPH1=NPH + NHC*(NHH-1)
                  DO NTH=1,NTHMAX
                  END DO
               END DO
               PABSR_SV(NR,NS)=PABSR_SV(NR,NS)
     &                        + PABSR(NR,NS)
            END DO
         END DO

         DO NS=1,NSMAX
            PABST3D(NPH,NS)=PABST(NS)
            PABST_SV(NS)=PABST_SV(NS)
     &                  + PABST(NS)
         END DO
         PABSTT3D(NPH)=PABSTT

C         write(6,'(A,1P2E12.4)') 
C     &        'CEFLD3D(1,5,NPH,5)=',CEFLD3D(1,5,NPH,5)
C
         IF(KNAMWM /=' ')THEN
           IF(NRANK.EQ.0) THEN
             KNAMWM_SAVE=KNAMWM
             IKNWM =index(KNAMWM_SAVE,' ') -1
             WRITE(KNHC,"(i2.2)")NHC
             IF (NPH0 .ge. 0)THEN
               WRITE(KNPH0,"(A1,i3.3)")"p",NPH0
             ELSE
               WRITE(KNPH0,"(A1,i3.3)")"m",-NPH0
             ENDIF
C             KNAMWM=''//KNAMWM_SAVE(1:IKNWM)//'_'//KNPH0//'_'//KNHC//''
             KNAMWM=''//KNAMWM_SAVE(1:IKNWM)//'_'//KNPH0//''
             CALL WMSAVE
             KNAMWM=KNAMWM_SAVE
           ENDIF
         ENDIF
      ENDDO
C
      NPH0  = NPH0_SV
      CALL WMEFLD_POST
      CALL WMPABS_POST
      CALL WMPOUT

!      DO NS=1,NSMAX
!         DO NR=1,NRMAX+1
!         PABSR(NR,NS) = PABSR_SV(NR,NS)
!         ENDDO
!         PABST(NS)=PABST_SV(NS)
!      ENDDO
!      DO NR=1,NRMAX+1
!         write(3133,*)"#",NR
!         DO NHH=1,NHHMAX
!            DO NTH=1,NTHMAX
!              write(3133,"(2(i6,1x),3(es15.8,1x))")NTH,NHH,
!     &                                PABS3D(NTH,NHH,NR,1:3)
!            ENDDO
!            write(3133,*)
!         ENDDO
!         write(3133,*)
!         write(3133,*)
!      ENDDO
C
      return
      end
