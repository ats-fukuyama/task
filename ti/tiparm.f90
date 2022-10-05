! MODULE tiparm

MODULE tiparm

  PRIVATE
  PUBLIC ti_parm,ti_view,ti_broadcast

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
      USE libkio
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: MODE
      CHARACTER(LEN=*),INTENT(IN)::  KIN
      INTEGER,INTENT(OUT):: IERR
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
      INTEGER,INTENT(IN) :: NID
      INTEGER,INTENT(OUT):: IST, IERR

      NAMELIST /TI/ &
           RR,RA,RKAP,RDLT,BB,RIP,NSMAX, &
           NPA,PM,PZ,PN,PNS,PT,PTS,PU,PUS, &
           ID_NS,NZMIN_NS,NZMAX_NS,NZINI_NS,KID_NS, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           MODELG,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF, &
           KNAMEQ,KNAMEQ2,KNAMTR, &
           PT,PROFJ1,PROFJ2,DT, &
           NRMAX,NTMAX,NTSTEP,NGTSTEP,NGRSTEP, &
           ngt_allocate_step,ngr_allocate_step, &
           EPSLOOP,MAXLOOP,EPSMAT,MATTYPE, &
           MODEL_EQB,MODEL_EQN,MODEL_EQT,MODEL_EQU, &
           MODEL_KAI,MODEL_DRR,MODEL_VR,MODEL_NC, &
           MODEL_NF,MODEL_NB,MODEL_EC,MODEL_LH,MODEL_IC, &
           MODEL_CD,MODEL_SYNC,MODEL_PEL,MODEL_PSC, &
           MODEL_BND,BND_VALUE, &
           DN0,DT0,DU0,VDN0,VDT0,VDU0,DR0,DRS, &
           DN0_NS,DT0_NS,DU0_NS,VDN0_NS,VDT0_NS,VDU0_NS, &
           PNBIN,PNBR0,PNBRW,PNBENG,PNBRTG, &
           PECIN,PECR0,PECRW,PECTOE,PECNPR, &
           PLHIN,PLHR0,PLHRW,PLHTOE,PLHNPR, &
           PICIN,PICR0,PICRW,PICTOE,PICNPR, &
           PNBCD,PECCD,PLHCD,PICCD,PBSCD, &
           PELIN,PELR0,PELRW,PELPAT, &
           PELTS,PELDT,PELTE,PELRAD,PELVEL, &
           PSCIN,PSCR0,PSCRW,PSCPAT, &
           SYNC_WALL,SYNC_CONV, &
           KNAMLOG,IDEBUG,glog_min, &
           adpost_dir,adpost_filename,adas_adf11_dir,adas_adf11_filename
      
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
      WRITE(6,'(A)') '        NPA,PM,PZ,PN,PNS,PT,PTS,PU,PUS,'
      WRITE(6,'(A)') '        ID_NS,NZMIN_NS,NZMAX_NS,NZINI_NS,KID_NS,'
      WRITE(6,'(A)') '        PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'
      WRITE(6,'(A)') '        MODELG,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF,'
      WRITE(6,'(A)') '        KNAMEQ,KNAMEQ2,KNAMTR,'
      WRITE(6,'(A)') '        PT,PROFJ1,PROFJ2,DT,'
      WRITE(6,'(A)') '        NRMAX,NTMAX,NTSTEP,NGTSTEP,NGRSTEP,'
      WRITE(6,'(A)') '        ngt_allocate_step,ngr_allocate_step,'
      WRITE(6,'(A)') '        EPSLOOP,MAXLOOP,EPSMAT,MATTYPE,'
      WRITE(6,'(A)') '        MODEL_EQB,MODEL_EQN,MODEL_EQT,MODEL_EQU,'
      WRITE(6,'(A)') '        MODEL_KAI,MODEL_DRR,MODEL_VR,MODEL_NC,'
      WRITE(6,'(A)') '        MODEL_NF,MODEL_NB,MODEL_EC,MODEL_LH,MODEL_IC,'
      WRITE(6,'(A)') '        MODEL_CD,MODEL_SYNC,MODEL_PEL,MODEL_PSC,'
      WRITE(6,'(A)') '        MODEL_BND,BND_VALUE,'
      WRITE(6,'(A)') '        DN0,DT0,DU0,VDN0,VDT0,VDU0,DR0,DRS,'
      WRITE(6,'(A)') '        DN0_NS,DT0_NS,DU0_NS,VDN0_NS,VDT0_NS,VDU0_NS,'
      WRITE(6,'(A)') '        PNBIN,PNBR0,PNBRW,PNBENG,PNBRTG,'
      WRITE(6,'(A)') '        PECIN,PECR0,PECRW,PECTOE,PECNPR,'
      WRITE(6,'(A)') '        PLHIN,PLHR0,PLHRW,PLHTOE,PLHNPR,'
      WRITE(6,'(A)') '        PICIN,PICR0,PICRW,PICTOE,PICNPR,'
      WRITE(6,'(A)') '        PNBCD,PECCD,PLHCD,PICCD,PBSCD,'
      WRITE(6,'(A)') '        PELIN,PELR0,PELRW,PELPAT(NSM),'
      WRITE(6,'(A)') '        PELTS,PELDT,PELTE,PELRAD,PELVEL,'
      WRITE(6,'(A)') '        PSCIN,PSCR0,PSCRW,PSCPAT(NSM),'
      WRITE(6,'(A)') '        SYNC_WALL,SYNC_CONV,KNAMLOG,IEBUG,'
      WRITE(6,'(A)') '        glog_min,'
      WRITE(6,'(A)') '        adpost_dir,adpost_filename,'
      WRITE(6,'(A)') '        adas_adf11_dir,adas_adf11_filename'
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
      INTEGER :: NS,NV

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
      WRITE(6,*) '    ','PROFN1      ','PROFN2      ', &
                        'PROTT1      ','PROFT2      ', &
                        'PROFU1      ','PROFU2      '

      DO NS=1,NSMAX
         WRITE(6,611) NS,PM(NS),PZ(NS),ID_NS(NS),NZMIN_NS(NS),NZMAX_NS(NS), &
                      NZINI_NS(NS),KID_NS(NS)
         WRITE(6,612) PN(NS),PNS(NS),PT(NS),PTS(NS),PU(NS),PUS(NS)
         WRITE(6,612) PROFN1(NS),PROFN2(NS),PROFT1(NS),PROFT2(NS),&
                      PROFU1(NS),PROFU2(NS)
      END DO
      
      WRITE(6,601) 'DT          ',DT, &
                   'PROFJ1      ',PROFJ1, &
                   'PROFJ2      ',PROFJ2
      WRITE(6,601) 'EPSLOOP     ',EPSLOOP, &
                   'EPSMAT      ',EPSMAT
      WRITE(6,602) 'MODELG      ',MODELG, &
                   'MODELN      ',MODELN, &
                   'MODELQ      ',MODELQ, &
                   'MODEL_NPROF ',MODEL_NPROF
      WRITE(6,602) 'NTMAX       ',NTMAX, &
                   'NTSTEP      ',NTSTEP, &
                   'NGTSTEP     ',NGTSTEP, &
                   'NGRSTEP     ',NGRSTEP
      WRITE(6,602) 'NRMAX       ',NRMAX, &
                   'ngt_alc_step',ngt_allocate_step, &
                   'ngr_alc_step',ngr_allocate_step
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
                   'MAXLOOP     ',MAXLOOP, &
                   'MATTYPE     ',MATTYPE, &
                   'IDEBUG      ',IDEBUG

      WRITE(6,'(A)') 'NS,NV,MODEL_BNDA,,BND_VALUE='
      DO NS=1,NSMAX
         WRITE(6,'(I5,3(2I5,1PE12.4))') &
              NS,(NV,MODEL_BND(NV,NS),BND_VALUE(NV,NS),NV=1,3)
      END DO

      WRITE(6,*)
      WRITE(6,601) 'DN0         ',DN0, &
                   'DT0         ',DT0, &
                   'DU0         ',DU0
      WRITE(6,601) 'VDN0        ',VDN0, &
                   'VDT0        ',VDT0, &
                   'VDU0        ',VDU0
      WRITE(6,'(A,A)') '   NS      DN0_NS      DT0_NS      DU0_NS', &
                          '     VDN0_NS     VDT0_NS     VDU0_NS'
      DO NS=1,NSMAX
         WRITE(6,'(I5,1P6E12.4)') &
              NS,DN0_NS(NS),DT0_NS(NS),DU0_NS(NS), &
                 VDN0_NS(NS),VDT0_NS(NS),VDU0_NS(NS)
      END DO

      WRITE(6,601) 'DR0         ',DR0, &
                   'DRS         ',DRS
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
                   'SYNC_CONV   ',SYNC_CONV, &
                   'glog_min    ',glog_min

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
  602 FORMAT(A12,'=',I5:2X,A12,'=',I5:2X,A12,'=',I5:2X,A12,'=',I5)
  611 FORMAT(I4,2(1PE12.4:),4I4,A)
  612 FORMAT(4X,6(1PE12.4:))
  613 FORMAT(I4,2(1PE12.4:))

    END SUBROUTINE ti_view

!     ***** BROADCAST INPUT PARAMETERS *****

    SUBROUTINE ti_broadcast

      USE ticomm_parm
      USE libmpi
      IMPLICIT NONE
      INTEGER,DIMENSION(99):: idata
      REAL(rkind),DIMENSION(99):: rdata
      INTEGER:: NS

!----- PL input parameters -----     

      idata( 1)=NSMAX
      idata( 2)=MODELG
      idata( 3)=MODELN
      idata( 4)=MODELQ
      idata( 5)=IDEBUG
      idata( 6)=MODEFR
      idata( 7)=MODEFW
      idata( 8)=MODEL_PROF
      idata( 9)=MODEL_NPROF

      CALL mtx_broadcast_integer(idata,9)
      NSMAX =idata( 1)
      MODELG=idata( 2)
      MODELN=idata( 3)
      MODELQ=idata( 4)
      IDEBUG=idata( 5)
      MODEFR=idata( 6)
      MODEFW=idata( 7)
      MODEL_PROF=idata( 8)
      MODEL_NPROF=idata( 9)

      rdata( 1)=RR
      rdata( 2)=RA
      rdata( 3)=RB
      rdata( 4)=RKAP
      rdata( 5)=RDLT
      rdata( 6)=BB
      rdata( 7)=Q0
      rdata( 8)=QA
      rdata( 9)=RIP
      rdata(10)=PROFJ
      rdata(11)=RHOMIN
      rdata(12)=QMIN
      rdata(13)=RHOEDG
      rdata(14)=RHOGMN
      rdata(15)=RHOGMX
      rdata(16)=PPN0
      rdata(17)=PTN0
      rdata(18)=RF_PL

      CALL mtx_broadcast_real8(rdata,18)
      RR    =rdata( 1)
      RA    =rdata( 2)
      RB    =rdata( 3)
      RKAP  =rdata( 4)
      RDLT  =rdata( 5)
      BB    =rdata( 6)
      Q0    =rdata( 7)
      QA    =rdata( 8)
      RIP   =rdata( 9)
      PROFJ =rdata(10)
      RHOMIN=rdata(11)
      QMIN  =rdata(12)
      RHOEDG=rdata(13)
      RHOGMN=rdata(14)
      RHOGMX=rdata(15)
      PPN0=rdata(16)
      PTN0=rdata(17)
      RF_PL=rdata(18)

      CALL mtx_broadcast_real8(PM,NSMAX)
      CALL mtx_broadcast_real8(PZ,NSMAX)
      CALL mtx_broadcast_integer(NPA,NSMAX)
      CALL mtx_broadcast_real8(PN,NSMAX)
      CALL mtx_broadcast_real8(PNS,NSMAX)
      CALL mtx_broadcast_real8(PTPR,NSMAX)
      CALL mtx_broadcast_real8(PTPP,NSMAX)
      CALL mtx_broadcast_real8(PTS,NSMAX)
      CALL mtx_broadcast_real8(PU,NSMAX)
      CALL mtx_broadcast_real8(PUS,NSMAX)
      CALL mtx_broadcast_real8(RHOITB,NSMAX)
      CALL mtx_broadcast_real8(PNITB,NSMAX)
      CALL mtx_broadcast_real8(PTITB,NSMAX)
      CALL mtx_broadcast_real8(PUITB,NSMAX)
      CALL mtx_broadcast_real8(PROFN1,NSMAX)
      CALL mtx_broadcast_real8(PROFN2,NSMAX)
      CALL mtx_broadcast_real8(PROFT1,NSMAX)
      CALL mtx_broadcast_real8(PROFT2,NSMAX)
      CALL mtx_broadcast_real8(PROFU1,NSMAX)
      CALL mtx_broadcast_real8(PROFU2,NSMAX)
      CALL mtx_broadcast_real8(PZCL,NSMAX)

      CALL mtx_broadcast_integer(ID_NS,NSMAX)
      DO NS=1,NSMAX
         CALL mtx_broadcast_character(KID_NS(NS),2)
      END DO

      CALL mtx_broadcast_real8(r_corner,3)
      CALL mtx_broadcast_real8(z_corner,3)
      CALL mtx_broadcast_real8(br_corner,3)
      CALL mtx_broadcast_real8(bz_corner,3)
      CALL mtx_broadcast_real8(bt_corner,3)
      DO NS=1,NSMAX
         CALL mtx_broadcast_real8(pn_corner(1:3,NS),3)
         CALL mtx_broadcast_real8(ptpr_corner(1:3,NS),3)
         CALL mtx_broadcast_real8(ptpp_corner(1:3,NS),3)
      END DO
      CALL mtx_broadcast_character(KNAMEQ,80)
      CALL mtx_broadcast_character(KNAMWR,80)
      CALL mtx_broadcast_character(KNAMFP,80)
      CALL mtx_broadcast_character(KNAMWM,80)
      CALL mtx_broadcast_character(KNAMPF,80)
      CALL mtx_broadcast_character(KNAMFO,80)
      CALL mtx_broadcast_character(KNAMTR,80)
      CALL mtx_broadcast_character(KNAMEQ2,80)

!----- TI input parameters -----

      idata( 1)=NRMAX
      idata( 2)=NTMAX
      idata( 3)=NTSTEP
      idata( 4)=NGTSTEP
      idata( 5)=NGRSTEP
      idata( 6)=ngt_allocate_step
      idata( 7)=ngr_allocate_step
      idata( 8)=0.D0
      idata( 9)=MAXLOOP
      idata(10)=MATTYPE
      idata(11)=MODEL_EQB
      idata(12)=MODEL_EQN
      idata(13)=MODEL_EQT
      idata(14)=MODEL_EQU
      idata(15)=MODEL_KAI
      idata(16)=MODEL_DRR
      idata(17)=MODEL_VR
      idata(18)=MODEL_NC
      idata(19)=MODEL_NF
      idata(20)=MODEL_NB
      idata(21)=MODEL_EC
      idata(22)=MODEL_LH
      idata(23)=MODEL_IC
      idata(24)=MODEL_CD
      idata(25)=MODEL_SYNC
      idata(26)=MODEL_PEL
      idata(27)=MODEL_PSC
      idata(28)=IDEBUG

      CALL mtx_broadcast_integer(idata,28)
      NRMAX=idata( 1)
      NTMAX=idata( 2)
      NTSTEP=idata( 3)
      NGTSTEP=idata( 4)
      NGRSTEP=idata( 5)
      ngt_allocate_step=idata( 6)
      ngr_allocate_step=idata( 7)
!      =idata( 8)
      MAXLOOP=idata( 9)
      MATTYPE=idata(10)
      MODEL_EQB=idata(11)
      MODEL_EQN=idata(12)
      MODEL_EQT=idata(13)
      MODEL_EQU=idata(14)
      MODEL_KAI=idata(15)
      MODEL_DRR=idata(16)
      MODEL_VR=idata(17)
      MODEL_NC=idata(18)
      MODEL_NF=idata(19)
      MODEL_NB=idata(20)
      MODEL_EC=idata(21)
      MODEL_LH=idata(22)
      MODEL_IC=idata(23)
      MODEL_CD=idata(24)
      MODEL_SYNC=idata(25)
      MODEL_PEL=idata(26)
      MODEL_PSC=idata(27)
      IDEBUG=idata(28)

      rdata( 1)=PROFJ1
      rdata( 2)=PROFJ2
      rdata( 3)=DT
      rdata( 4)=EPSLOOP
      rdata( 5)=EPSMAT
      rdata( 6)=DN0
      rdata( 7)=DT0
      rdata( 8)=DU0
      rdata( 9)=VDN0
      rdata(10)=VDT0
      rdata(11)=VDU0
      rdata(12)=DR0
      rdata(13)=DRS
      rdata(14)=PNBIN
      rdata(15)=PNBR0
      rdata(16)=PNBRW
      rdata(17)=PNBENG
      rdata(18)=PNBRTG
      rdata(19)=PECIN
      rdata(20)=PECR0
      rdata(21)=PECRW
      rdata(22)=PECTOE
      rdata(23)=PECNPR
      rdata(24)=PLHR0
      rdata(25)=PLHRW
      rdata(26)=PLHTOE
      rdata(27)=PLHNPR
      rdata(28)=PICIN
      rdata(29)=PICR0
      rdata(30)=PICRW
      rdata(31)=PICTOE
      rdata(32)=PICNPR
      rdata(33)=PNBCD
      rdata(34)=PECCD
      rdata(35)=PLHCD
      rdata(36)=PICCD
      rdata(37)=PBSCD
      rdata(38)=PELIN
      rdata(39)=PELR0
      rdata(40)=PELRW
      rdata(41)=PELDT
      rdata(42)=PELTE
      rdata(43)=PELRAD
      rdata(44)=PELVEL
      rdata(45)=PSCIN
      rdata(46)=PSCR0
      rdata(47)=PSCRW
      rdata(48)=SYNC_WALL
      rdata(49)=SYNC_CONV
      rdata(50)=glog_min

      CALL mtx_broadcast_real8(rdata,50)

      PROFJ1=rdata( 1)
      PROFJ2=rdata( 2)
      DT=rdata( 3)
      EPSLOOP=rdata( 4)
      EPSMAT=rdata( 5)
      DN0=rdata( 6)
      DT0=rdata( 7)
      DU0=rdata( 8)
      VDN0=rdata( 9)
      VDT0=rdata(10)
      VDU0=rdata(11)
      DR0=rdata(12)
      DRS=rdata(13)
      PNBIN=rdata(14)
      PNBR0=rdata(15)
      PNBRW=rdata(16)
      PNBENG=rdata(17)
      PNBRTG=rdata(18)
      PECIN=rdata(19)
      PECR0=rdata(20)
      PECRW=rdata(21)
      PECTOE=rdata(22)
      PECNPR=rdata(23)
      PLHR0=rdata(24)
      PLHRW=rdata(25)
      PLHTOE=rdata(26)
      PLHNPR=rdata(27)
      PICIN=rdata(28)
      PICR0=rdata(29)
      PICRW=rdata(30)
      PICTOE=rdata(31)
      PICNPR=rdata(32)
      PNBCD=rdata(33)
      PECCD=rdata(34)
      PLHCD=rdata(35)
      PICCD=rdata(36)
      PBSCD=rdata(37)
      PELIN=rdata(38)
      PELR0=rdata(39)
      PELRW=rdata(40)
      PELDT=rdata(41)
      PELTE=rdata(42)
      PELRAD=rdata(43)
      PELVEL=rdata(44)
      PSCIN=rdata(45)
      PSCR0=rdata(46)
      PSCRW=rdata(47)
      SYNC_WALL=rdata(48)
      SYNC_CONV=rdata(49)
      glog_min=rdata(50)

      CALL mtx_broadcast_integer(NZMIN_NS,NSMAX)
      CALL mtx_broadcast_integer(NZMAX_NS,NSMAX)
      CALL mtx_broadcast_integer(NZINI_NS,NSMAX)
      CALL mtx_broadcast_real8(PT,NSMAX)
      CALL mtx_broadcast_real8(PELPAT,NSMAX)
      CALL mtx_broadcast_real8(PSCPAT,NSMAX)
      CALL mtx_broadcast_real8(DN0_NS,NSMAX)
      CALL mtx_broadcast_real8(DT0_NS,NSMAX)
      CALL mtx_broadcast_real8(DU0_NS,NSMAX)
      CALL mtx_broadcast_real8(VDN0_NS,NSMAX)
      CALL mtx_broadcast_real8(VDT0_NS,NSMAX)
      CALL mtx_broadcast_real8(VDU0_NS,NSMAX)

      DO NS=1,NSMAX
         CALL mtx_broadcast_integer(MODEL_BND(1:3,NS),3)
         CALL mtx_broadcast_real8(BND_VALUE(1:3,NS),3)
      END DO

    END SUBROUTINE ti_broadcast

  END MODULE tiparm
