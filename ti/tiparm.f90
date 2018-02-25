! MODULE tiparm

MODULE tiparm

  PRIVATE
  PUBLIC ti_parm,ti_view

CONTAINS

!     ***********************************************************

!           PARAMETER INPUT

!     ***********************************************************

  SUBROUTINE ti_parm(MODE,KIN,IERR)

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

      USE ticomm,ONLY: NSMAX,PT,PTPR,PTPP
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: MODE
      CHARACTER(LEN=*),INTENT(IN)::  KIN
      INTEGER(4),INTENT(OUT):: IERR
      INTEGER:: NS

      DO NS=1,NSMAX
         PT(NS)=0.D0
      END DO

    1 CALL TASK_PARM(MODE,'TI',KIN,ti_nlin,ti_plst,IERR)

      DO NS=1,NSMAX
         IF(PT(NS).EQ.0.D0) THEN
            PT(NS)=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
         ELSE
            PTPR(NS)=PT(NS)
            PTPP(NS)=PT(NS)
         END IF
      END DO

      IF(IERR.NE.0) RETURN

      CALL ti_check(IERR)
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100

      RETURN
  END SUBROUTINE ti_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE ti_nlin(NID,IST,IERR)

      USE ticomm_parm
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NID
      INTEGER(4),INTENT(OUT):: IST, IERR

      NAMELIST /TI/ &
           RR,RA,RKAP,RDLT,BB,RIP,NSMAX, &
           PM,PZ,PN,PNS,PT,PTS,PU,PUS, &
           ID_NS,NZMIN_NS,NZMAX_NS,NZINI_NS,KID_NS, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           MODELG,MODELN,MODELQ,MODEL_NPROF, &
           KNAMEQ,KNAMEQ2,KNAMTR, &
           PT,PROFJ1,PROFJ2,DT, &
           NRMAX,NTMAX,NTSTEP,NGTSTEP,NGRSTEP, &
           NGTMAX_INIT,NGRMAX_INIT,NTSTEP_COEF,EPSLTI,LMAXTI, &
           MODEL_EQB,MODEL_EQN,MODEL_EQT,MODEL_EQU, &
           MODEL_KAI,MODEL_DRR,MODEL_VR,MODEL_NC, &
           MODEL_NF,MODEL_NB,MODEL_EC,MODEL_LH,MODEL_IC, &
           MODEL_CD,MODEL_SYNC,MODEL_PEL,MODEL_PSC, &
           AK0,AD0,AV0,DK0,DKS, &
           PNBIN,PNBR0,PNBRW,PNBENG,PNBRTG, &
           PECIN,PECR0,PECRW,PECTOE,PECNPR, &
           PLHIN,PLHR0,PLHRW,PLHTOE,PLHNPR, &
           PICIN,PICR0,PICRW,PICTOE,PICNPR, &
           PNBCD,PECCD,PLHCD,PICCD,PBSCD, &
           PELIN,PELR0,PELRW,PELPAT, &
           PELTS,PELDT,PELTE,PELRAD,PELVEL, &
           PSCIN,PSCR0,PSCRW,PSCPAT, &
           SYNC_WALL,SYNC_CONV, &
           KNAMLOG
      
      IERR=0

      READ(NID,TI,IOSTAT=IST,ERR=9800,END=9900)
      IF(IST.NE.0) THEN
         WRITE(6,'(A,I5)') 'XX ti_nlin: READ ERROR: IOSTAT=',IST
         IERR=IST
         RETURN
      END IF

      RETURN

 9800 IERR=8
      WRITE(6,'(A)') 'XX ti_nlin: READ ERROR'
      RETURN
 9900 IERR=9
      WRITE(6,'(A)') 'XX ti_nlin: READ END OF FILE'
      RETURN
    END SUBROUTINE ti_nlin

!     ***** INPUT PARAMETER LIST *****

    SUBROUTINE ti_plst

      WRITE(6,'(A)') '# &TI : RR,RA,RKAP,RDLT,BB,RIP,NSMAX,'
      WRITE(6,'(A)') '        PM,PZ,PN,PNS,PT,PTS,PU,PUS,'
      WRITE(6,'(A)') '        ID_NS,NZMIN_NS,NZMAX_NS,NZINI_NS,KID_NS,'
      WRITE(6,'(A)') '        PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'
      WRITE(6,'(A)') '        MODELG,MODELN,MODELQ,MODEL_NPROF,'
      WRITE(6,'(A)') '        KNAMEQ,KNAMEQ2,KNAMTR,'
      WRITE(6,'(A)') '        PT,PROFJ1,PROFJ2,DT,'
      WRITE(6,'(A)') '        NRMAX,NTMAX,NTSTEP,NGTSTEP,NGRSTEP,NTSTEP_CONF,'
      WRITE(6,'(A)') '        NGTMAX_INIT,NGRMAX_INIT,EPSLTI,LMAXTI,'
      WRITE(6,'(A)') '        MODEL_EQB,MODEL_EQN,MODEL_EQT,MODEL_EQU,'
      WRITE(6,'(A)') '        MODEL_KAI,MODEL_DRR,MODEL_VR,MODEL_NC,'
      WRITE(6,'(A)') '        MODEL_NF,MODEL_NB,MODEL_EC,MODEL_LH,MODEL_IC,'
      WRITE(6,'(A)') '        MODEL_CD,MODEL_SYNC,MODEL_PEL,MODEL_PSC,'
      WRITE(6,'(A)') '        AK0,AD0,AV0,DK0,DKS,'
      WRITE(6,'(A)') '        PNBIN,PNBR0,PNBRW,PNBENG,PNBRTG,'
      WRITE(6,'(A)') '        PECIN,PECR0,PECRW,PECTOE,PECNPR,'
      WRITE(6,'(A)') '        PLHIN,PLHR0,PLHRW,PLHTOE,PLHNPR,'
      WRITE(6,'(A)') '        PICIN,PICR0,PICRW,PICTOE,PICNPR,'
      WRITE(6,'(A)') '        PNBCD,PECCD,PLHCD,PICCD,PBSCD,'
      WRITE(6,'(A)') '        PELIN,PELR0,PELRW,PELPAT(NSM),'
      WRITE(6,'(A)') '        PELTS,PELDT,PELTE,PELRAD,PELVEL,'
      WRITE(6,'(A)') '        PSCIN,PSCR0,PSCRW,PSCPAT(NSM),'
      WRITE(6,'(A)') '        SYNC_WALL,SYNC_CONV,KNAMLOG'
      RETURN
    END SUBROUTINE ti_plst

!     ***** CHECK INPUT PARAMETERS *****

    SUBROUTINE ti_check(IERR)

      USE ticomm_parm
      IMPLICIT NONE
      INTEGER, INTENT(OUT):: IERR

      IERR=0


      IF(NRMAX.LT.1) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NRMAX'
         WRITE(6,*) '                  NRMAX =',NRMAX
         IERR=1
      ENDIF

      RETURN
    END SUBROUTINE ti_check

!     ***********************************************************

!           VIEW INPUT PARAMETER

!     ***********************************************************

    SUBROUTINE ti_view

      USE ticomm_parm
      IMPLICIT NONE
      INTEGER :: NS,NPSC

      WRITE(6,601) 'RR          ',RR          , &
                   'RA          ',RA          , &
                   'RKAP        ',RKAP
      WRITE(6,601) 'RDLT        ',RDLT        , &
                   'BB          ',BB          , &
                   'RIP         ',RIP
      WRITE(6,*) 'NS  ','PM          ','PZ          ', &
                 'ID  ','MIN ','MAX ','INI ','KID '
      WRITE(6,*) '    ','PN          ','PNS         ', &
                        'PT          ','PTS         ', &
                        'PU          ','PUS         '

      DO NS=1,NSMAX
         WRITE(6,611) NS,PM(NS),PZ(NS),ID_NS(NS),NZMIN_NS(NS),NZMAX_NS(NS), &
                      NZINI_NS(NS),KID_NS(NS)
         WRITE(6,612) PN(NS),PNS(NS),PT(NS),PTS(NS),PU(NS),PUS(NS)
      END DO
      
      WRITE(6,601) 'PROFN1      ',PROFN1, &
                   'PROFT1      ',PROFT1, &
                   'PROFU1      ',PROFU1
      WRITE(6,601) 'PROFN2      ',PROFN2, &
                   'PROFT2      ',PROFT2, &
                   'PROFU2      ',PROFU2
      WRITE(6,601) 'DT          ',DT, &
                   'PROFJ1      ',PROFJ1, &
                   'PROFJ2      ',PROFJ2
      WRITE(6,601) 'EPSLTI      ',EPSLTI
      WRITE(6,602) 'MODELG      ',MODELG, &
                   'MODELN      ',MODELN, &
                   'MODELQ      ',MODELQ, &
                   'MODEL_NPROF ',MODEL_NPROF
      WRITE(6,602) 'NTMAX       ',NTMAX, &
                   'NTSTEP      ',NTSTEP, &
                   'NGTSTEP     ',NGTSTEP, &
                   'NGRSTEP     ',NGRSTEP
      WRITE(6,602) 'NRMAX       ',NRMAX, &
                   'NTSTEP_COEF ',NTSTEP_COEF, &
                   'NGTMAX_INIT ',NGTMAX_INIT, &
                   'NGRMAX_INIT ',NGRMAX_INIT
      WRITE(6,602) 'MODEL_EQB   ',MODEL_EQB, &
                   'MODEL_EQN   ',MODEL_EQN, &
                   'MODEL_EQT   ',MODEL_EQT, &
                   'MODEL_EQU   ',MODEL_EQU
      WRITE(6,602) 'MODEL_KAI   ',MODEL_KAI, &
                   'MODEL_DRR   ',MODEL_DRR, &
                   'MODEL_VR    ',MODEL_VR, &
                   'MODEL_NC    ',MODEL_NC
      WRITE(6,602) 'MODEL_NF    ',MODEL_NF, &
                   'MODEL_NB    ',MODEL_NB, &
                   'MODEL_EC    ',MODEL_EC, &
                   'MODEL_LH    ',MODEL_LH
      WRITE(6,602) 'MODEL_IC    ',MODEL_IC, &
                   'MODEL_CD    ',MODEL_CD, &
                   'MODEL_PEL   ',MODEL_PEL, &
                   'MODEL_PSC   ',MODEL_PSC
      WRITE(6,602) 'MODEL_SYNC  ',MODEL_SYNC, &
                   'LMAXTI      ',LMAXTI
      WRITE(6,*)
      WRITE(6,601) 'AK0         ',AK0, &
                   'AD0         ',AD0, &
                   'AV0         ',AV0
      WRITE(6,601) 'DK0         ',DK0, &
                   'DKS         ',DKS
      WRITE(6,601) 'PNBIN       ',PNBIN, &
                   'PNBR0       ',PNBR0, &
                   'PNBRW       ',PNBRW
      WRITE(6,601) 'PNBENG      ',PNBENG, &
                   'PNBRTG      ',PNBRTG
      WRITE(6,601) 'PECIN       ',PECIN, &
                   'PECR0       ',PECR0, &
                   'PECRW       ',PECRW
      WRITE(6,601) 'PECTOE      ',PECTOE, &
                   'PECNPR      ',PECNPR
      WRITE(6,601) 'PLHIN       ',PLHIN, &
                   'PLHR0       ',PLHR0, &
                   'PLHRW       ',PLHRW
      WRITE(6,601) 'PLHTOE      ',PLHTOE, &
                   'PLHNPR      ',PLHNPR
      WRITE(6,601) 'PICIN       ',PICIN, &
                   'PICR0       ',PICR0, &
                   'PICRW       ',PICRW
      WRITE(6,601) 'PICTOE      ',PICTOE, &
                   'PICNPR      ',PICNPR
      WRITE(6,601) 'PNBCD       ',PNBCD, &
                   'PECCD       ',PECCD, &
                   'PLHCD       ',PLHCD
      WRITE(6,601) 'PICCD       ',PICCD, &
                   'PBSCD       ',PBSCD
      WRITE(6,601) 'PELIN       ',PELIN, &
                   'PELR0       ',PELR0, &
                   'PELRW       ',PELRW
      WRITE(6,601) 'PELTS       ',PELTS, &
                   'PELDT       ',PELDT, &
                   'PELTE       ',PELTE
      WRITE(6,601) 'PELRAD      ',PELRAD, &
                   'PELVEL      ',PELVEL
      WRITE(6,601) 'PSCIN       ',PSCIN, &
                   'PSCR0       ',PSCR0, &
                   'PSCRW       ',PSCRW
      WRITE(6,601) 'SYNC_WALL   ',SYNC_WALL, &
                   'SYNC_CONV   ',SYNC_CONV

      IF(PELIN.GT.0.D0.OR.PSCIN.GT.0.D0) THEN
         WRITE(6,'(A4,2A12)') 'NS  ','PELPAT      ','PSCPAT      '
         DO NS=1,NSMAX
            WRITE(6,613) NS,PELPAT(NS),PSCPAT(NS)
         END DO
      END IF

      WRITE(6,'(A,A60)') 'KNAMEQ =',KNAMEQ
      WRITE(6,'(A,A60)') 'KNAMEQ2=',KNAMEQ2
      WRITE(6,'(A,A60)') 'KNAMTR =',KNAMTR
      WRITE(6,'(A,A60)') 'KNAMLOG=',KNAMLOG

      WRITE(6,'(A)')  ' *** END OF INPUT PARAMETER LIST *** '
      RETURN

  601 FORMAT(A12,'=',1PE12.4:2X,A12,'=',1PE12.4:2X,A12,'=',1PE12.4)
  602 FORMAT(A12,'=',I5:2X,A12,'=',I5:2X,A12,'=',I5,:2X,A12,'=',I5)
  611 FORMAT(I4,2(1PE12.4:),4I4,A)
  612 FORMAT(4X,6(1PE12.4:))
  613 FORMAT(I4,2(1PE12.4:))

    END SUBROUTINE ti_view
  END MODULE tiparm
