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
    INTEGER:: NPH0_save,NPH,NR,NTH,NS,NHH,NPOW,NPP
    CHARACTER(LEN=80):: KNAMWM_SAVE
    CHARACTER(LEN=2):: KNHC
    CHARACTER(LEN=4):: KNPH0

    NPH0_save  = NPH0

    IF(nhhmax.EQ.1) THEN ! tokamak
       IF(nphmax.EQ.1) THEN ! single mode
          nph0_npp(1)=nph0
          nph_nhh_npp(1,1)=nph0
       ELSE                 ! multi mode
          DO npp=1,nppmax
             nph0_npp(npp)=npp-nppmax/2-1        ! -NPHMAX/2:NPHMAX/2-1
             nph_nhh_npp(1,npp)=npp-nppmax/2-1   ! -NPHMAX/2:NPHMAX/2-1
          END DO
       END IF
    ELSE                 ! helical
       IF(nppmax.EQ.1) THEN ! single mode-group
          nph0_npp(1)=nph0
          DO nhh=1,nhhmax
             nph_nhh_npp(nhh,1)=nph0+nhc*(nhh-nhhmax/2-1)
          END DO
       ELSE                 ! multi mode-group
          DO NPP=1,NPPMAX
             nph0_npp(npp)=npp-nppmax/2-1
             DO nhh=1,nhhmax
                nph_nhh_npp(nhh,npp)=npp-nppmax/2-1+nhc*(nhh-nhhmax/2-1)
             END DO
          END DO
       END IF
    END IF

    IF(nrank.EQ.0) WRITE(6,'(A,4I8)') &
         '== nhhmax,nppmax,nphmax,nphtot=', &
         nhhmax,nppmax,nphmax,nphtot

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
      
    DO NPP = 1,NPPMAX
       NPH0 = nph0_npp(NPP)
       IF(nrank.EQ.0) WRITE(6,'(A,I4)') &
            "== Toroidal mode number : ",NPH0

       CALL wm_exec(IERR)
       IF(IERR.NE.0) EXIT

       DO NR=1,NRMAX+1
          DO NHH=1,NHHMAX
             NPH=nph_nhh_npp(nhh,npp)
             DO NTH=1,NTHMAX
                CEFLD3D(1,NTH,NPH,NR)=CEFLD3D(1,NTH,NPH,NR) &
                     +CEFLD(1,NTH,NHH,NR)
                CEFLD3D(2,NTH,NPH,NR)=CEFLD3D(2,NTH,NPH,NR) &
                     +CEFLD(2,NTH,NHH,NR) 
                CEFLD3D(3,NTH,NPH,NR)=CEFLD3D(3,NTH,NPH,NR) &
                     +CEFLD(3,NTH,NHH,NR)
                CBFLD3D(1,NTH,NPH,NR)=CBFLD3D(1,NTH,NPH,NR) &
                     +CBFLD(1,NTH,NHH,NR)
                CBFLD3D(2,NTH,NPH,NR)=CBFLD3D(2,NTH,NPH,NR) &
                     +CBFLD(2,NTH,NHH,NR)
                CBFLD3D(3,NTH,NPH,NR)=CBFLD3D(3,NTH,NPH,NR) &
                     +CBFLD(3,NTH,NHH,NR)
                CEFLDK3D(1,NTH,NPH,NR)=CEFLDK3D(1,NTH,NPH,NR) &
                     +CEFLDK(1,NTH,NHH,NR)
                CEFLDK3D(2,NTH,NPH,NR)=CEFLDK3D(2,NTH,NPH,NR) &
                     +CEFLDK(2,NTH,NHH,NR)
                CEFLDK3D(3,NTH,NPH,NR)=CEFLDK3D(3,NTH,NPH,NR) &
                     +CEFLDK(3,NTH,NHH,NR)
                CBFLDK3D(1,NTH,NPH,NR)=CBFLDK3D(1,NTH,NPH,NR) &
                     +CBFLDK(1,NTH,NHH,NR)
                CBFLDK3D(2,NTH,NPH,NR)=CBFLDK3D(2,NTH,NPH,NR) &
                     +CBFLDK(2,NTH,NHH,NR)
                CBFLDK3D(3,NTH,NPH,NR)=CBFLDK3D(3,NTH,NPH,NR) &
                     +CBFLDK(3,NTH,NHH,NR)
             END DO
          END DO
       END DO

       DO NS=1,NSMAX
          DO NR=1,NRMAX
             DO NHH=1,NHHMAX
                NPH=nph_nhh_npp(nhh,npp)
                DO NTH=1,NTHMAX
                   CPABS3D(NTH,NPH,NR,NS)=CPABS3D(NTH,NPH,NR,NS) &
                        +CPABS(NTH,NHH,NR,NS)
                   CPABSK3D(NTH,NPH,NR,NS)=CPABSK3D(NTH,NPH,NR,NS) &
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

       IF(idebuga(71).NE.0) THEN
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
       END IF
    END DO ! npp

    CALL wm_pout_sum
    CALL wm_pout

    NPH0  = NPH0_save

    RETURN
  END SUBROUTINE WM_LOOP
END MODULE wmloop
