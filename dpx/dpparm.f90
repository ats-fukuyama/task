MODULE DPPARM

  PRIVATE
  PUBLIC dp_parm,dp_chek,dpprep,dpprep_local,dp_view

CONTAINS

!     ****** INPUT PARAMETERS ******

  SUBROUTINE DP_PARM(MODE,KIN,IERR)

!     MODE=0 : standard namelinst input
!     MODE=1 : namelist file input
!     MODE=2 : namelist line input

!     IERR=0 : normal end
!     IERR=1 : namelist standard input error
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE
!     IERR=10X : input parameter out of range

    IMPLICIT NONE
    INTEGER,INTENT(IN):: MODE
    CHARACTER(LEN=*),INTENT(IN):: KIN
    INTEGER,INTENT(OUT):: IERR

    IERR=0

1   CALL TASK_PARM(MODE,'DP',KIN,DPNLIN,DPPLST,IERR)
    IF(IERR.NE.0) RETURN

    CALL EQCHEK(IERR)
    CALL DP_CHEK(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
    CALL DPPREP_LOCAL(IERR)
    RETURN
  END SUBROUTINE DP_PARM

!     ****** INPUT NAMELIST ******

  SUBROUTINE DPNLIN(NID,IST,IERR)

    USE dpcomm_parm,DPX=>DP
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NID
    INTEGER,INTENT(OUT):: IST,IERR
    INTEGER:: NS

    NAMELIST /DP/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
                  NSMAX,NPA,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS, &
                  ID_NS,KID_NS, &
                  PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
                  RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB, &
                  MODELG,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF, &
                  RHOGMN,RHOGMX, &
                  KNAMEQ,KNAMWR,KNAMFP,MODEFR,MODEFW,IDEBUG, &
                  MODELP,MODELV,NDISP1,NDISP2,NS_NSA,PMAX,EMAX, &
                  RF0,RFI0,RKX0,RKY0,RKZ0,RX0,RY0,RZ0, &
                  RF1,RFI1,RKX1,RKY1,RKZ1,RX1,RY1,RZ1, &
                  RF2,RFI2,RKX2,RKY2,RKZ2,RX2,RY2,RZ2, &
                  NGXMAX,NGYMAX,NGPMAX,EPSRT,LMAXRT, &
                  NPMAX,NTHMAX,NRMAX,NSAMAX,RMIN,RMAX

    READ(NID,DP,IOSTAT=IST,ERR=9800,END=9900)
    IF(MODEL_PROF.EQ.0) THEN
       DO NS=1,NSMAX
          PROFN1(NS)=PROFN1(1)
          PROFN2(NS)=PROFN2(1)
          PROFT1(NS)=PROFT1(1)
          PROFT2(NS)=PROFT2(1)
          PROFU1(NS)=PROFU1(1)
          PROFU2(NS)=PROFU2(1)
       END DO
    END IF
    IERR=0
    RETURN

9800 IERR=8
    RETURN
9900 IERR=9
    RETURN
  END SUBROUTINE DPNLIN

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE DPPLST

    WRITE(6,601)
    RETURN

601 FORMAT(' ','# &DP : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'/ &
            9X,'NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,'/ &
            9X,'ID_NS,KID_NS,'/ &
            9X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'/ &
            9X,'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,'/ &
            9X,'MODELG,MODELN,MODELQ,MODEFA,'/ &
            9X,'KNAMEQ,KNAMWR,KNAMFP,MODEFR,MODEFW,IDEBUG,'/ &
            9X,'MODELP,MODELV,NDISP1,NDISP2,PMAX,'/ &
            9X,'RF0,RFI0,RKX0,RKY0,RKZ0,RX0,RY0,RZ0,'/ &
            9X,'RF1,RFI1,RKX1,RKY1,RKZ1,RX1,'/ &
            9X,'RF2,RFI2,RKX2,RKY2,RKZ2,RX2,'/ &
            9X,'NGXMAX,NGYMAX,NGPMAX,EPSRT,LMAXRT,'/ &
            9X,'NPMAX,NTHMAX,NRMAX,NSAMAX,RMIN,RMAX')
  END SUBROUTINE DPPLST

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE DP_CHEK(IERR)

    USE dpcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR

    IERR=0

    RETURN
  END SUBROUTINE DP_CHEK

!     ***** Setup velocity distribution function *****

  SUBROUTINE DPPREP(NTHMAX_1,NRMAX_1,RMIN_1,RMAX_1,IERR)

    USE dpcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NTHMAX_1,NRMAX_1,RMIN_1,RMAX_1
    INTEGER,INTENT(OUT):: IERR

    NTHMAX=NTHMAX_1
    NRMAX=NRMAX_1
    RMIN=RMIN_1
    RMAX=RMAX_1
    CALL DPPREP_LOCAL(IERR)
    RETURN
  END SUBROUTINE DPPREP

  SUBROUTINE DPPREP_LOCAL(IERR)

    USE dpcomm
    USE dpfpin
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NS,IND

    IERR=0
    IND=0
    DO NS=1,NSMAX
       IF(MODELV(NS).EQ.2.OR. &
          MODELV(NS).EQ.4.OR. &
          MODELV(NS).EQ.9) IND=1
    END DO

    IF(IND.EQ.1) THEN
       CALL DPLDFP(IERR)
       IF(IERR.NE.0) RETURN
    END IF

    DO NS=1,NSMAX
       IF(MODELV(NS).EQ.1) THEN
          CALL DPLDFM(NS,0,IERR)
          IF(IERR.NE.0) RETURN
       END IF
       IF(MODELV(NS).EQ.3) THEN
          CALL DPLDFM(NS,1,IERR)
          IF(IERR.NE.0) RETURN
       END IF
    END DO
    RETURN
  END SUBROUTINE DPPREP_LOCAL

!     ****** SHOW PARAMETERS ******

  SUBROUTINE DP_VIEW

    USE dpcomm_parm
    IMPLICIT NONE
    INTEGER:: NS

    WRITE(6,100)
    DO NS=1,NSMAX
       WRITE(6,110) &
            NS,MODELP(NS),MODELV(NS),NDISP1(NS),NDISP2(NS),PMAX(NS)
    ENDDO
    WRITE(6,601) 'RMIN  ',RMIN  ,'RMAX  ',RMAX
    WRITE(6,602) 'NPMAX ',NPMAX ,'NTHMAX',NTHMAX, &
                 'NRMAX ',NRMAX ,'NSAMAX',NSAMAX

    RETURN

100 FORMAT(1H ,'NS    MODELP  MODELV  NDISP1  NDISP2  PMAX')
110 FORMAT(1H ,I2,' ',4I8,1PE12.4)                               
601 FORMAT(' ',A6,'=',1PE11.3 :2X,A6,'=',1PE11.3: &
            2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
602 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
            2X,A6,'=',I7,4X  :2X,A6,'=',I7)
  END SUBROUTINE DP_VIEW
END MODULE DPPARM
