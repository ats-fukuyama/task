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

      IF(NPHMAX.EQ.1) THEN
         CALL WMEXEC(IERR)
         RETURN
      ENDIF

      NPH0_SV  = NPH0
!      NPHMAX_SV = NPHMAX
!      NHHMAX_SV = NHHMAX
!
!      IF(NHHMAX.GT.1) THEN
!         NHHMAX=NPHMAX/NHC
!         NPHMAX=NHC
!      ELSE
!         NHC=1
!      END IF
   
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
      
      CALL WMSETG(IERR)
      CALL WMPOUT_INIT
      DO NPH = 1,NPHMAX
         NPH0 = NPH-NPHMAX/2-1 + NPH0_SV
         WRITE(6,'(A,I4)') "== Toroidal mode number : ",NPH0
C
         CALL WMEXEC(IERR)
         CALL mtx_barrier
         IF(IERR.NE.0) EXIT

C         DO NR=1,NRMAX+1
C            DO NHH=1,NHHMAX
C               DO NTH=1,NTHMAX
C               CEFLD3D(1,NTH,NHH,NR)=CEFLD(1,NTH,NHH,NR)
C               CEFLD3D(2,NTH,NHH,NR)=CEFLD(2,NTH,NHH,NR)
C               CEFLD3D(3,NTH,NHH,NR)=CEFLD(3,NTH,NHH,NR)
C               END DO
C            END DO
C         END DO
C
         DO NS=1,NSMAX
            DO NR=1,NRMAX+1
               DO NHH=1,NHHMAX
                  NPH1=NPH + NHC*(NHH-1)
                  DO NTH=1,NTHMAX
CC                  CPABS3D(NTH,NPH1,NR,NS)=CPABS(NTH,NHH,NR,NS)
CC                   PABS3D(NTH,NPH1,NR,NS)= PABS(NTH,NHH,NR,NS)
CC                  CPABS3D(NTH,NHH,NR,NS)=CPABS3D(NTH,NHH,NR,NS)
CC     &                          + CPABS(NTH,NHH,NR,NS)
CC                   PABS3D(NTH,NHH,NR,NS)= PABS3D(NTH,NHH,NR,NS)
CC     &                          + PABS(NTH,NHH,NR,NS)
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
