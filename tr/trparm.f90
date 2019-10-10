! trparm.f90

MODULE trparm

  PRIVATE
  PUBLIC tr_parm,tr_nlin,tr_view

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
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: MODE
      CHARACTER(LEN=*),INTENT(IN)::  KIN
      INTEGER(4),INTENT(OUT):: IERR

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
      INTEGER(4),INTENT(IN) :: NID
      INTEGER(4),INTENT(OUT):: IST, IERR

      NAMELIST /TR/ RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RHOA, &
                    PA,PZ,PN,PNS,PT,PTS,PNC,PNFE,PNNU,PNNUS, &
                    PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
                    PROFJ1,PROFJ2,ALP,AD0,AV0,CNP,CNH,CDP,CDH,CNN,CDW, &
                    CWEB,CALF,CNB,CSPRS, &
                    MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,MDLTPF, &
                    DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST, &
                    EPSLTR,LMAXTR,CHP,CK0,CK1,CKALFA,CKBETA,CKGUMA, &
                    TPRST,CDW, &
                    MDLST,MDLNF,IZERO,MODELG,NTEQIT,MDEDGE, &
                    MDLXP,MDLUF,MDNCLS,MDLWLD,MDLFLX,MDLER,MDCD05, &
                    PNBTOT,PNBR0,PNBRW,PNBVY,PNBVW,PNBENG,PNBRTG,MDLNB, &
                    PECTOT,PECR0,PECRW,PECTOE,PECNPR,MDLEC, &
                    PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,MDLLH, &
                    PICTOT,PICR0,PICRW,PICTOE,PICNPR,MDLIC, &
                    PNBCD,PECCD,PLHCD,PICCD,PBSCD,MDLCD, &
                    PELTOT,PELR0,PELRW,PELRAD,PELVEL,MDLPEL, &
                    MDLPR,SYNCABS,SYNCSELF, &
                    PELTIM,PELPAT,KNAMEQ,KNAMEQ2,KNAMTR,KFNLOG, &
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
             ' ',8X,'MDEDGE'/ &
             ' ',8X,'MDLXP,MDLUF,MDNCLS,MDLWLD,MDLFLX,MDLER,MDCD05'/ &
             ' ',8X,'MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0'/ &
             ' ',8X,'MDLEQE,MDLEOI,NSMAX,NSZMAX,NSNMAX,KUFDIR,KUFDEV,KUFDCG'/ &
             ' ',8X,'TIME_INT,MODEP,MDNI,MDLJQ,MDTC,MDLPCK'/ &
             ' ',8X,'KNAMEQ,KNAMEQ2,KNAMTR,KFNLOG'/ &
             ' ',8X,'MDLPSC,NPSCMAX,PSCTOT,PSCR0,PSCRW,NSPSC')
    END SUBROUTINE trplst

!     ***** CHECK INPUT PARAMETERS *****

    SUBROUTINE trchek(IERR)

      USE TRCOMM, ONLY : NGTSTP, NRMAX, NTM, NTMAX
      IMPLICIT NONE
      INTEGER(4), INTENT(OUT):: IERR


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

!     ***********************************************************

!           VIEW INPUT PARAMETER

!     ***********************************************************

    SUBROUTINE tr_view(ID)

      USE trcomm

      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ID
      INTEGER(4) :: NS,NPSC


      WRITE(6,*) '** TRANSPORT **'
      WRITE(6,602) 'MDLEQB',MDLEQB,'MDLEQN',MDLEQN,'MDLEQT',MDLEQT,'MDLEQU',MDLEQU
      WRITE(6,602) 'MDLEQZ',MDLEQZ,'MDLEQ0',MDLEQ0,'MDLEQE',MDLEQE,'MDLEOI',MDLEOI
      WRITE(6,602) 'NSMAX ',NSMAX, 'NSZMAX',NSZMAX,'NSNMAX',NSNMAX
      WRITE(6,601) 'RR    ',RR,    'RA    ',RA,    'RKAP  ',RKAP,  'RDLT  ',RDLT
      WRITE(6,601) 'RIPS  ',RIPS,  'RIPE  ',RIPE,  'BB    ',BB

      WRITE(6,611)
  611 FORMAT(' ','NS',2X,'PA           PZ      PN(E20)  PNS(E20) ','PT(KEV)  PTS(KEV) PELPAT')
      DO NS=1,NSMAX
         WRITE(6,612) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PT(NS),PTS(NS),PELPAT(NS)
  612    FORMAT(' ',I2,1PD12.4,0P,F8.3,5F9.4)
      ENDDO

      WRITE(6,601) 'PNC   ',PNC,   'PNFE  ',PNFE,  'PNNU  ',PNNU,  'PNNUS ',PNNUS
      WRITE(6,601) 'PROFN1',PROFN1,'PROFT1',PROFT1,'PROFU1',PROFU1,'PROFJ1',PROFJ1
      WRITE(6,601) 'PROFN2',PROFN2,'PROFT2',PROFT2,'PROFU2',PROFU2,'PROFJ2',PROFJ2
      WRITE(6,601) 'ALP(1)',ALP(1),'ALP(2)',ALP(2),'ALP(3)',ALP(3),'PBSCD ',PBSCD
      WRITE(6,602) 'MDLKAI',MDLKAI,'MDLETA',MDLETA,'MDLAD ',MDLAD, 'MDLAVK',MDLAVK
      WRITE(6,602) 'MDLJBS',MDLJBS,'MDLKNC',MDLKNC,'MDLTPF',MDLTPF,'MDNCLS',MDNCLS
      WRITE(6,604) 'MDLUF ',MDLUF, 'KUFDEV',KUFDEV,'KUFDCG',KUFDCG,'MDNI  ',MDNI
      WRITE(6,605) 'MDLJQ ',MDLJQ, 'MDLFLX',MDLFLX,'MDTC  ',MDTC,  'RHOA  ',RHOA
      WRITE(6,602) 'MDLWLD',MDLWLD,'MDLER ',MDLER, 'MODELG',MODELG,'NTEQIT',NTEQIT
      WRITE(6,603) 'MDCD05',MDCD05,'CK0   ',CK0,   'CK1   ',CK1
      WRITE(6,603) 'MDEDGE',MDEDGE,'CSPRS ',CSPRS, 'CNN   ',CNN
      WRITE(6,601) 'CNP   ',CNP,   'CNH   ',CNH,   'CDP   ',CDP,   'CDH   ',CDH
      WRITE(6,601) 'AD0   ',AD0,   'CHP   ',CHP,   'CWEB  ',CWEB,  'CALF  ',CALF
      IF((MDLKAI.GE.1.AND.MDLKAI.LT.10).OR.ID.EQ.1) &
         WRITE(6,601) 'CKALFA',CKALFA,'CKBETA',CKBETA,'CKGUMA',CKGUMA

      IF((MDLKAI.GE.10.AND.MDLKAI.LT.20).OR.ID.EQ.1)  &
         WRITE(6,613) CDW(1),CDW(2),CDW(3),CDW(4),CDW(5),CDW(6),CDW(7),CDW(8)
  613 FORMAT(' ','    AKDW(E) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/ &
             ' ','    AKDW(D) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/ &
             ' ','    AKDW(T) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/ &
             ' ','    AKDW(A) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW')

      WRITE(6,601) 'DT    ',DT,    'EPSLTR',EPSLTR,'TSST  ',TSST,  'TPRST ',TPRST
      WRITE(6,602) 'LMAXTR',LMAXTR,'NRMAX ',NRMAX, 'NTMAX ',NTMAX, 'NTSTEP',NTSTEP
      WRITE(6,602) 'NGRSTP',NGRSTP,'NGTSTP',NGTSTP,'NGPST ',NGPST, 'IZERO ',IZERO
      WRITE(6,602) 'MDLST ',MDLST, 'MDLCD ',MDLCD, 'MDLNF ',MDLNF

      IF((PNBTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PNBTOT',PNBTOT,'PNBR0 ',PNBR0,'PNBRW ',PNBRW,'PNBENG',PNBENG
         WRITE(6,603) 'MDLNB ',MDLNB, 'PNBRTG',PNBRTG,'PNBCD ',PNBCD,'PNBVY ',PNBVY
         WRITE(6,601) 'PNBVW ',PNBVW
      ENDIF

      IF((PECTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PECTOT',PECTOT, 'PECR0 ',PECR0,'PECRW ',PECRW,'PECTOE',PECTOE
         WRITE(6,603) 'MDLEC ',MDLEC,  'PECNPR',PECNPR,'PECCD ',PECCD
      ENDIF

      IF((PLHTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PLHTOT',PLHTOT, 'PLHR0 ',PLHR0,'PLHRW ',PLHRW, 'PLHTOE',PLHTOE
         WRITE(6,603) 'MDLLH ',MDLLH,  'PLHNPR',PLHNPR,'PLHCD ',PLHCD
      ENDIF

      IF((PICTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PICTOT',PICTOT, 'PICR0 ',PICR0, 'PICRW ',PICRW,'PICTOE',PICTOE
         WRITE(6,603) 'MDLIC ',MDLIC,  'PICNPR',PICNPR,'PICCD ',PICCD
      ENDIF

      IF((PELTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PELTOT',PELTOT,'PELR0 ',PELR0,'PELRW ',PELRW
         WRITE(6,603) 'MDLPEL',MDLPEL,'PELRAD',PELRAD,'PELVEL',PELVEL,'PELTIM',PELTIM
      ENDIF
      IF(MDLPSC.GT.0) THEN
         WRITE(6,622) 'MDLPSC  ',MDLPSC,'NPSCMAX ',NPSCMAX
         DO NPSC=1,NPSCMAX
            WRITE(6,603) 'NSPSC ',NSPSC(NPSC) ,'PSCTOT',PSCTOT(NPSC), &
                         'PSCR0 ',PSCR0(NPSC) ,'PSCRW ',PSCRW(NPSC)
         END DO
      END IF

      IF((MDLPR.GE.1).OR.(ID.EQ.1)) THEN
         WRITE(6,623) 'MDLPR   ',MDLPR,   'SYNCABS ',SYNCABS, &
                      'SYNCSELF',SYNCSELF
      ENDIF

      WRITE(6,601) 'CNB   ',CNB
      RETURN

  601 FORMAT(' ',A6,'=',1PE11.3 :2X,A6,'=',1PE11.3: &
              2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
  602 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',I7,4X  : &
              2X,A6,'=',I7,4X   :2X,A6,'=',I7)
  603 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',1PE11.3: &
              2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
  604 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',1X,A6,4X: &
              2X,A6,'=',1X,A6,4X:2X,A6,'=',I7)
  605 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',I7,4X  : &
              2X,A6,'=',I7,4X   :2X,A6,'=',1PE11.3)
  622 FORMAT(' ',A8,'=',I5,4X   :2X,A8,'=',I5,4X  : &
              2X,A8,'=',I5,4X   :2X,A8,'=',I5)
  623 FORMAT(' ',A8,'=',I7,4X   :2X,A8,'=',1PE11.3: &
              2X,A8,'=',1PE11.3)
    END SUBROUTINE tr_view
END MODULE trparm
