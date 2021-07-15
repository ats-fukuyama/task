! wmloop.f90

MODULE wmloop

  PRIVATE
  PUBLIC wm_loop

CONTAINS

!***********************************************************************

!     Multi-mode calculation

!***********************************************************************

  SUBROUTINE WM_LOOP(IERR)

    USE wmcomm
    USE wmsetg
    USE wmexec
    USE wmpout
    USE wmfile
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NPH0_SV,NPH,NR,NTH,NS,NHH,NPH1
    CHARACTER(LEN=80):: KNAMWM_SAVE
    CHARACTER(LEN=2):: KNHC
    CHARACTER(LEN=4):: KNPH0

    NPH0_SV  = NPH0

    IF(NPHMAX.EQ.1) THEN
       NPH_LOOP(1)=NPH0
    ELSE
       DO NPH=1,NPHMAX
          NPH_LOOP(NPH)=NPH-NPHMAX/2-1
       END DO
    END IF

    DO NR=1,NRMAX+1
       DO NPH=1,NPHMAX
          DO NTH=1,NTHMAX
             CEFLD3D(1,NTH,NPH,NR)=(0.D0,0.D0)
             CEFLD3D(2,NTH,NPH,NR)=(0.D0,0.D0)
             CEFLD3D(3,NTH,NPH,NR)=(0.D0,0.D0)
             CBFLD3D(1,NTH,NPH,NR)=(0.D0,0.D0)
             CBFLD3D(2,NTH,NPH,NR)=(0.D0,0.D0)
             CBFLD3D(3,NTH,NPH,NR)=(0.D0,0.D0)
             CEFLDK3D(1,NTH,NPH,NR)=(0.D0,0.D0)
             CEFLDK3D(2,NTH,NPH,NR)=(0.D0,0.D0)
             CEFLDK3D(3,NTH,NPH,NR)=(0.D0,0.D0)
             CBFLDK3D(1,NTH,NPH,NR)=(0.D0,0.D0)
             CBFLDK3D(2,NTH,NPH,NR)=(0.D0,0.D0)
             CBFLDK3D(3,NTH,NPH,NR)=(0.D0,0.D0)
          END DO
       END DO
    END DO
    DO NS=1,NSMAX
       DO NR=1,NRMAX
          DO NPH=1,NPHMAX
             DO NTH=1,NTHMAX
                CPABS3D(NTH,NPH,NR,NS)=0.D0
                CPABSK3D(NTH,NPH,NR,NS)=0.D0
             END DO
          END DO
       END DO
    END DO
      
    CALL wm_setg(IERR)
    CALL wm_pout_init
      
    DO NPH = 1,NPHMAX
       NPH0 = NPH_LOOP(NPH)
       IF(nrank.EQ.0) WRITE(6,'(A,I4)') "== Toroidal mode number : ",NPH0

       CALL wm_exec(IERR)
       IF(IERR.NE.0) EXIT

       DO NR=1,NRMAX+1
          DO NHH=1,NHHMAX
             NPH1=NPH+NHC*(NHH-1)
             DO NTH=1,NTHMAX
                CEFLD3D(1,NTH,NPH1,NR)=CEFLD3D(1,NTH,NPH1,NR) &
                     +CEFLD(1,NTH,NHH,NR)
                CEFLD3D(2,NTH,NPH1,NR)=CEFLD3D(2,NTH,NPH1,NR) &
                     +CEFLD(2,NTH,NHH,NR) 
                CEFLD3D(3,NTH,NPH1,NR)=CEFLD3D(3,NTH,NPH1,NR) &
                     +CEFLD(3,NTH,NHH,NR)
                CBFLD3D(1,NTH,NPH1,NR)=CBFLD3D(1,NTH,NPH1,NR) &
                     +CBFLD(1,NTH,NHH,NR)
                CBFLD3D(2,NTH,NPH1,NR)=CBFLD3D(2,NTH,NPH1,NR) &
                     +CBFLD(2,NTH,NHH,NR)
                CBFLD3D(3,NTH,NPH1,NR)=CBFLD3D(3,NTH,NPH1,NR) &
                     +CBFLD(3,NTH,NHH,NR)
                CEFLDK3D(1,NTH,NPH1,NR)=CEFLDK3D(1,NTH,NPH1,NR) &
                     +CEFLDK(1,NTH,NHH,NR)
                CEFLDK3D(2,NTH,NPH1,NR)=CEFLDK3D(2,NTH,NPH1,NR) &
                     +CEFLDK(2,NTH,NHH,NR)
                CEFLDK3D(3,NTH,NPH1,NR)=CEFLDK3D(3,NTH,NPH1,NR) &
                     +CEFLDK(3,NTH,NHH,NR)
                CBFLDK3D(1,NTH,NPH1,NR)=CBFLDK3D(1,NTH,NPH1,NR) &
                     +CBFLDK(1,NTH,NHH,NR)
                CBFLDK3D(2,NTH,NPH1,NR)=CBFLDK3D(2,NTH,NPH1,NR) &
                     +CBFLDK(2,NTH,NHH,NR)
                CBFLDK3D(3,NTH,NPH1,NR)=CBFLDK3D(3,NTH,NPH1,NR) &
                     +CBFLDK(3,NTH,NHH,NR)
             END DO
          END DO
       END DO

       DO NS=1,NSMAX
          DO NR=1,NRMAX
             DO NHH=1,NHHMAX
                NPH1=NPH + NHC*(NHH-1)
                DO NTH=1,NTHMAX
                   CPABS3D(NTH,NPH1,NR,NS)=CPABS3D(NTH,NPH1,NR,NS) &
                        +CPABS(NTH,NHH,NR,NS)
                   CPABSK3D(NTH,NPH1,NR,NS)=CPABSK3D(NTH,NPH1,NR,NS) &
                        +CPABSK(NTH,NHH,NR,NS)
                END DO
             END DO
          END DO
       END DO

       DO NS=1,NSMAX
          DO NR=1,NRMAX
             PABSR3D(NPH,NR,NS)=PABSR3D(NPH,NR,NS)+PABSR(NR,NS)
          END DO
       END DO

       DO NS=1,NSMAX
          PABST3D(NPH,NS)=PABST3D(NPH,NS)+PABST(NS)
       END DO

       PABSTT3D(NPH)=PABSTT

!         write(6,'(A,1P2E12.4)') 
!     &        'CEFLD3D(1,5,NPH,5)=',CEFLD3D(1,5,NPH,5)

       IF(KNAMWM /=' ')THEN
          IF(NRANK.EQ.0) THEN
             KNAMWM_SAVE=KNAMWM
             WRITE(KNHC,"(i2.2)") NHC
             IF (NPH0 .ge. 0)THEN
                WRITE(KNPH0,"(A1,i3.3)")"p",NPH0
             ELSE
                WRITE(KNPH0,"(A1,i3.3)")"m",-NPH0
             ENDIF
             KNAMWM=''//TRIM(KNAMWM_SAVE)//'_'//KNPH0//''
             CALL wm_save(ierr)
             KNAMWM=KNAMWM_SAVE
          ENDIF
       ENDIF
    ENDDO

    WRITE(6,*) '@@@ point 1'
    CALL wm_pout_sum
    WRITE(6,*) '@@@ point 2'
    CALL wm_pout
    WRITE(6,*) '@@@ point 3'

    NPH0  = NPH0_SV

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

    RETURN
  END SUBROUTINE WM_LOOP
END MODULE wmloop
