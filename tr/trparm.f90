! trparm.f90

MODULE trparm

  PRIVATE
  PUBLIC tr_parm,tr_nlin

CONTAINS

!     ***********************************************************

!           PARAMETER INPUT

!     ***********************************************************

  SUBROUTINE tr_parm(MODE,KIN,IERR)

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

      USE TRCOMM
      USE libkio
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: MODE
      CHARACTER(LEN=*),INTENT(IN)::  KIN
      INTEGER,INTENT(OUT):: IERR

    1 CALL TASK_PARM(MODE,'TR',KIN,tr_nlin,trplst,IERR)
      IF(IERR.NE.0) RETURN

      CALL TRCHEK(IERR)
      NTMAX_SAVE=NTMAX
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100

      RETURN
    END SUBROUTINE tr_parm

!     ****** INPUT NAMELIST ******

    SUBROUTINE tr_nlin(NID,IST,IERR)

      USE trcomm
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NID
      INTEGER,INTENT(OUT):: IST, IERR

      NAMELIST /TR/ RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RHOA, &
                    PA,PZ,PN,PNS,PT,PTS,PNC,PNFE,PNNU,PNNUS, &
                    PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
                    PROFJ1,PROFJ2,ALP,AD0,AV0,CNP,CNH,CDP,CDH,CNN,CDW, &
                    CWEB,CALF,CNB,CSPRS, &
                    MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,MDLTPF, &
                    DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST, &
                    EPSLTR,LMAXTR,CHP,CK0,CK1,CKALFA,CKBETA,CKGUMA, &
                    TPRST,CDW, &
                    MDLST,MDLNF,IZERO,MODELG,NTEQIT,MDEDGE,MDLIMP, &
                    MDLXP,MDLUF,MDNCLS,MDLWLD,MDLFLX,MDLER,MDCD05, &
                    PNBTOT,PNBR0,PNBRW,PNBVY,PNBVW,PNBENG,PNBRTG,MDLNB, &
                    PECTOT,PECR0,PECRW,PECTOE,PECNPR,MDLEC, &
                    PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,MDLLH, &
                    PICTOT,PICR0,PICRW,PICTOE,PICNPR,MDLIC, &
                    PNBCD,PECCD,PLHCD,PICCD,PBSCD,MDLCD, &
                    PELTOT,PELR0,PELRW,PELRAD,PELVEL,MDLPEL, &
                    MDLPR,SYNCABS,SYNCSELF, &
                    PELTIM,PELPAT,KNAMEQ,KNAMEQ2,KNAMTR,KFNLOG,KFNTXT,KFNCVS, &
                    MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE, &
                    MDLEOI,NSMAX,NSZMAX,NSNMAX, &
                    KUFDIR,KUFDEV,KUFDCG,TIME_INT,MODEP,MDNI,MDLJQ,MDTC, &
                    MDLPCK,MDLPSC,NPSCMAX,PSCTOT,PSCR0,PSCRW,NSPSC


      IF(NID.GE.0) THEN
         READ(NID,TR,IOSTAT=IST,ERR=9800,END=9900)
         NTMAX_SAVE=NTMAX
      ENDIF
      IST=0
      IERR=0
      RETURN

 9800 IERR=8
      RETURN
 9900 IERR=9
      RETURN
    END SUBROUTINE tr_nlin

!     ***** INPUT PARAMETER LIST *****

    SUBROUTINE trplst

      WRITE(6,601)
      RETURN

  601 FORMAT(' ','# &TR : RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RHOA'/ &
             ' ',8X,'(PA,PZ,PN,PNS,PT,PTS:NSM)'/ &
             ' ',8X,'PNC,PNFE,PNNU,PNNUS'/ &
             ' ',8X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2'/ &
             ' ',8X,'PROFJ1,PROFJ2,ALP'/ &
             ' ',8X,'CK0,CK1,CNP,CNH,CDP,CDH,CNN,CDW,CNB,CSPRS'/ &
             ' ',8X,'CWEB,CALF,CKALFA,CKBETA,MDLKNC,MDLTPF'/ &
             ' ',8X,'AD0,CHP,MDLAD,MDLAVK,CKGUMA,MDLKAI,MDLETA,MDLJBS'/ &
             ' ',8X,'DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST'/ &
             ' ',8X,'EPSLTR,LMAXTR,PRST,MDLST,MDLNF,IZERO,PBSCD,MDLCD'/ &
             ' ',8X,'PNBTOT,PNBR0,PNBRW,PNBVY,PNBVW,PNBENG,PNBRTG'/ &
             ' ',8X,'PNBCD,MDLNB'/ &
             ' ',8X,'PECTOT,PECR0,PECRW,PECTOE,PECNPR,PECCD,MDLEC'/ &
             ' ',8X,'PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,PLHCD,MDLLH'/ &
             ' ',8X,'PICTOT,PICR0,PICRW,PICTOE,PICNPR,PICCD,MDLIC'/ &
             ' ',8X,'PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM,MDLPEL'/ &
             ' ',8X,'PELTIM,PELPAT,MDLPR,SYNCABS,SYNCSELF,MODELG,NTEQIT'/&
             ' ',8X,'MDEDGE,MDLIMP,'/ &
             ' ',8X,'MDLXP,MDLUF,MDNCLS,MDLWLD,MDLFLX,MDLER,MDCD05'/ &
             ' ',8X,'MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0'/ &
             ' ',8X,'MDLEQE,MDLEOI,NSMAX,NSZMAX,NSNMAX,KUFDIR,KUFDEV,KUFDCG'/ &
             ' ',8X,'TIME_INT,MODEP,MDNI,MDLJQ,MDTC,MDLPCK'/ &
             ' ',8X,'KNAMEQ,KNAMEQ2,KNAMTR,KFNLOG,KFNTXT,KFNCVS,'/ &
             ' ',8X,'MDLPSC,NPSCMAX,PSCTOT,PSCR0,PSCRW,NSPSC')
    END SUBROUTINE trplst

!     ***** CHECK INPUT PARAMETERS *****

    SUBROUTINE trchek(IERR)

      USE TRCOMM, ONLY : NGTSTP, NRMAX, NTM, NTMAX
      IMPLICIT NONE
      INTEGER, INTENT(OUT):: IERR


      IERR=0

      IF(NRMAX.LT.1) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NRMAX'
         WRITE(6,*) '                  NRMAX =',NRMAX
         IERR=1
      ENDIF

      IF(NTMAX.LT.0.OR.NTMAX/NGTSTP.GT.NTM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NTMAX'
         WRITE(6,*) '                  NTMAX,NTM =',NTMAX,NTM
         NTMAX=NTM*NGTSTP
         IERR=1
      ENDIF

      RETURN
    END SUBROUTINE trchek

END MODULE trparm
