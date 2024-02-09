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

      USE trcomm
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

      NAMELIST /TR/ &
           RR,RA,RB,RKAP,RDLT,BB,RIPS,RIPE,RHOA, &
           NSMAX,NSZMAX,NSNMAX, &
           PA,PZ,PN,PNS,PT,PTS,PU,PUS, &
           MDLIMP,PNC,PNFE,PNNU,PNNUS, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           PROFNU1,PROFNU2,PROFJ1,PROFJ2,ALP, &
           model_prof,knam_prof,knam_nfixed,knam_tfixed, &
           AD0,AV0,CNP,CNH,CDP,CDH,CNN,CDW, &
           MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,MDLTPF, &
           MDLWLD,MDLDW,MDLCD05, &
           CHP,CK0,CK1,CWEB,CALF,CKALFA,CKBETA,CKGUMA, &
           DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST, &
           EPSLTR,LMAXTR,TPRST,MDLST,IZERO,MDLNF, &
           NNBMAX,NECMAX,NLHMAX,NICMAX,NPELMAX,NPSCMAX, &
           MDLNB,MDLEC,MDLLH,MDLIC,MDLPEL,MDLPSC, &
           PNBIN,PNBR0,PNBRW,PNBCD,PNBVY,PNBVW,PNBENG,PNBRTG,NSPNB,NRNBMAX, &
           PECIN,PECR0,PECRW,PECCD,PECTOE,PECNPR, &
           PLHIN,PLHR0,PLHRW,PLHCD,PLHTOE,PLHNPR, &
           PICIN,PICR0,PICRW,PICCD,PICTOE,PICNPR, &
           PELIN,PELR0,PELRW,PELRAD,PELVEL,PELTIM, &
           pellet_time_start,pellet_time_interval,PELPAT, &
           number_of_pellet_repeat, &
           PSCIN,PSCR0,PSCRW,NSPSC, &
           MDLCD,PBSCD, &
           MDLPR,SYNCABS,SYNCSELF, &
           MDLEDGE,CSPRS, &
           MODELG,NTEQIT,MODEP,MDLJQ,MDLFLX, &
           MDLXP,MDLUF, &
           MDLNI, &
           MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE,MDLEOI, &
           MDLER,MDLNCL,NSLMAX, &
           MDLELM,ELMWID,ELMDUR,ELMNRD,ELMTRD,ELMENH, &
           MDLTC,MDLPCK,model_nfixed,model_tfixed, &
           KUFDIR,KUFDev,KUFDCG, &
           KNAMEQ,KNAMEQ2,KNAMTR,KFNLOG,KFNTXT,KFNCVS

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

  601 FORMAT(' ','# &TR : RR,RA,RB,RKAP,RDLT,BB,RIPS,RIPE,RHOA'/ &
             ' ',8X,'(PA,PZ,PN,PNS,PT,PTS:NSM)'/ &
             ' ',8X,'PNC,PNFE,PNNU,PNNUS'/ &
             ' ',8X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2'/ &
             ' ',8X,'PROFJ1,PROFJ2,ALP'/ &
             ' ',8X,'CK0,CK1,CNP,CNH,CDP,CDH,CNN,CDW,CSPRS'/ &
             ' ',8X,'CWEB,CALF,CKALFA,CKBETA,MDLKNC,MDLTPF'/ &
             ' ',8X,'AD0,CHP,MDLAD,MDLAVK,CKGUMA,MDLKAI,MDLETA,MDLJBS'/ &
             ' ',8X,'DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST'/ &
             ' ',8X,'EPSLTR,LMAXTR,PRST,MDLST,MDLNF,IZERO,PBSCD,MDLCD'/ &
             ' ',8X,'PNBTOT,PNBR0,PNBRW,PNBVY,PNBVW,PNBENG,PNBRTG'/ &
             ' ',8X,'PNBCD,MDLNB,NSPNB,NRNBMAX'/ &
             ' ',8X,'PECIN,PECR0,PECRW,PECTOE,PECNPR,PECCD,MDLEC'/ &
             ' ',8X,'PLHIN,PLHR0,PLHRW,PLHTOE,PLHNPR,PLHCD,MDLLH'/ &
             ' ',8X,'PICIN,PICR0,PICRW,PICTOE,PICNPR,PICCD,MDLIC'/ &
             ' ',8X,'PELIN,PELR0,PELRW,PELRAD,PELVEL,PELTIM,MDLPEL'/ &
             ' ',8X,'PELTIM,PELPAT'/ &
             ' ',8X,'pellet_time_start,pellet_time_interval'/ &
             ' ',8X,'number_of_pellet_repeat'/ &
             ' ',8X,'MDLPR,SYNCABS,SYNCSELF,MODELG,NTEQIT'/&
             ' ',8X,'MDEDGE,MDLIMP,model_prof,knam_prof,'/ &
             ' ',8X,'MDLXP,MDLUF,MDLNCL,MDLWLD,MDLFLX,MDLER,MDCD05'/ &
             ' ',8X,'MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0'/ &
             ' ',8X,'MDLEQE,MDLEOI,NSMAX,NSZMAX,NSNMAX,KUFDIR,KUFDEV,KUFDCG'/ &
             ' ',8X,'TIME_INT,MODEP,MDNI,MDLJQ,MDLTC,MDLPCK'/ &
             ' ',8X,'KNAMEQ,KNAMEQ2,KNAMTR,KFNLOG,KFNTXT,KFNCVS,'/ &
             ' ',8X,'MDLPSC,NPSCMAX,PSCIN,PSCR0,PSCRW,NSPSC'/ &
             ' ',8X,'knam_nfixed,knam_tfixed'/ &
             ' ',8X,'model_nfixed,model_tfixed')
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
