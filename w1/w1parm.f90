Module w1parm

  USE w1comm_parm
  NAMELIST /W1/ &
       BB,RR,RZ,RA,RD,RB,RF,WALLR,APRFPN,APRFTR,APRFTP, &
       RKZ,DRF,DRKZ,DXFACT,DXWDTH,NAMAX, &
       AJYH,AJZH,APYH,APZH,ALZH,APHH,AJYL,AJZL,APYL,APZL,ALZL,APHL,&
       PA,PZ,PN,PTPP,PTPR,PU,PNS,PTS,PZCL,NSMAX, &
       NXMAX,NZMAX,NPRINT,NFILE,NGRAPH,NLOOP,NSYM, &
       NMODEL,NALPHA,NDMAX,XDMAX,IHARM,NSYS,NGDSP,MODELN, &
       EPSH,ZEFF,WVYSIZ,NCDTYP,NXABS,IELEC, &
       MDLWG,MDLWGS,WGZ1,WGZ2,WGAMP,WGNZ,job_id,nfile_data
  
  PRIVATE
  PUBLIC w1_parm
  PUBLIC w1_view
  PUBLIC w1_save_parm
  PUBLIC w1_load_parm

CONTAINS

!     ****** PARAMETER INPUT ******

  SUBROUTINE w1_parm(MODE,KIN,IERR)

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

    USE libkio
    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode
    CHARACTER(LEN=*),INTENT(IN)::  kin
    INTEGER,INTENT(OUT):: ierr

1   CALL TASK_PARM(MODE,'w1',KIN,w1_nlin,w1_plst,IERR)
    IF(IERR.NE.0 .AND. IERR.NE.2) RETURN

    CALl w1_check(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
    IF(IERR.NE.0) IERR=IERR+100

    RETURN
  END SUBROUTINE w1_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE w1_nlin(NID,IST,IERR)

    USE w1comm_parm

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NID
    INTEGER,INTENT(OUT) :: IST,IERR


    RZ=0.D0

    READ(NID,W1,IOSTAT=IST,ERR=9800,END=9900)

    IF(RZ.EQ.0.D0) RZ=2.D0*PI*(RR+RB)

    IERR=0
    RETURN

9800 IERR=8
    RETURN
9900 IERR=9
    RETURN
  END SUBROUTINE w1_nlin

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE w1_plst

    IMPLICIT NONE
    WRITE(6,'(A)') &
         '# &W1 : BB,RR,RZ,RA,RD,RB,RF,WALLR,APRFPN,APRFTR,APRFTP,',&
         '        RKZ,DRF,DRKZ,DXFACT,DXWDTH,NXMAX, ',&
         '        AJYH,AJZH,APYH,APZH,ALZH,APHH, ',&
         '        AJYL,AJZL,APYL,APZL,ALZL,APHL,', &
         '        PA,PZ,PN,PTPP,PTPR,PU,PNS,PTS,PZCL,NSMAX,', &
         '        NXMAX,NZMAX,NPRINT,NFILE,NGRAPH,NLOOP,NSYM,', &
         '        NMODEL,NALPHA,NDMAX,XDMAX,IHARM,NSYS,NGDSP,MODELN,', &
         '        EPSH,ZEFF,WVYSIZ,NCDTYP,NXABS,IELEC', &
         '        MDLWG,MDLWGS,WGZ1,WGZ2,WGAMP,WGNZ,job_id,nfile_data'
    RETURN

  END SUBROUTINE w1_plst

!     ****** CHECK INPUT PARAMETER ******

  SUBROUTINE w1_check(IERR)

    USE w1comm_parm
    IMPLICIT NONE
    INTEGER:: IERR

    IERR=0

    IF(NXMAX < 0) THEN
       WRITE(6,'(A,I8)') 'W1 w1_check: INVALID nxmax: nxmax=',nxmax
       IERR=1
    ENDIF

    IF(NAMAX.GT.1.AND.NZMAX.EQ.1) THEN
       WRITE(6,'(A,I8)') 'W1 w1_check: NAMAX shoule be 1 for NZMAX=1',namax
       IERR=2
    ENDIF

    RETURN
  END SUBROUTINE w1_check

!     ****** SHOW PARAMETERS ******

  SUBROUTINE w1_view

    use w1comm_parm
    implicit none
    INTEGER:: NS,NA

    WRITE(6,601) 'BB    ',BB    ,'RR    ',RR, &
                 'RZ    ',RZ    ,'RA    ',RA, &
                 'RD    ',RD    ,'RB    ',RB, &
                 'RF    ',RF    ,'WALLR ',WALLR, &
                 'APRFPN',APRFPN,'APRFTR',APRFTR, &
                 'APRFTP',APRFTP,'RKZ   ',RKZ, &
                 'DRF   ',DRF   ,'DRKZ  ',DRKZ, &
                 'DXFACT',DXFACT,'DXWDTH',DXWDTH, &
                 'EPSH  ',EPSH  ,'ZEFF  ',ZEFF, &
                 'WVYSIZ',WVYSIZ,'XDMAX ',XDMAX, &
                 'WGZ1  ',WGZ1  ,'WGZ1  ',WGZ2,  &
                 'WGAMP ',WGAMP ,'WGNZ  ',WGNZ

    WRITE(6,'(A)') '     ', &
          '  NS/IELEC   PA/PZ      PN/PNS    PTPR/PTS     PTPP/PU    PZCL'
    DO NS=1,NSMAX
       WRITE(6,'(A,I6,1P5E12.4)') &
            '     ',NS,PA(NS),PN(NS),PTPR(NS),PTPP(NS),PZCL(NS)
       WRITE(6,'(A,I6,1P4E12.4)') &
            '     ',IELEC(NS),PZ(NS),PNS(NS),PTS(NS),PU(NS)
    END DO

    WRITE(6,'(A,A,A)') '   ', &
          '    NA AJYH/AJYL   AJZH/AJZL   APYH/APYL   APZH/APZL', &
              '   ALZH/ALZL   APHH/APHL'
    DO NA=1,NAMAX
       IF(AJYH(NA).NE.0.D0) &
       WRITE(6,'(A,I2,1P6E12.4)') &
            ' ANT-H',NA,AJYH(NA),AJZH(NA),APYH(NA),APZH(NA),ALZH(NA),APHH(NA)
       IF(AJYL(NA).NE.0.D0) &
       WRITE(6,'(A,I2,1P6E12.4)') &
            ' ANT-L',NA,AJYL(NA),AJZL(NA),APYL(NA),APZL(NA),ALZL(NA),APHL(NA)
    END DO

    WRITE(6,602) 'NXMAX ',NXMAX , &
                 'NZMAX ',NZMAX ,'NSMAX ',NSMAX, &
                 'NAMAX ',NAMAX ,'NDMAX ',NDMAX, &
                 'NPRINT',NPRINT,'NFILE ',NFILE, &
                 'NGRAPH',NGRAPH,'NLOOP ',NLOOP, &
                 'NSYM  ',NSYM  ,'NMODEL',NMODEL, &
                 'NALPHA',NALPHA,'MODELN',MODELN, &
                 'NSYS  ',NSYS  ,'NGDSP ',NGDSP, &
                 'NCDTYP',NCDTYP,'NXABS ',NXABS, &
                 'MDLWG ',MDLWG ,'MDLWGS',MDLWGS
    WRITE(6,'(A)') TRIM(job_id)
    WRITE(6,'(A,I6)') 'nfile_data = ',nfile_data
    RETURN

601 FORMAT(    A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
            2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
602 FORMAT(    A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
            2X,A6,'=',I7,4X  :2X,A6,'=',I7)
  END SUBROUTINE w1_view

!     ****** save NAMELIST w1 ******

  SUBROUTINE w1_save_parm(nfc)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: nfc

    WRITE(nfc,w1)
    RETURN
  END SUBROUTINE w1_save_parm

!     ****** load NAMELIST w1 ******

  SUBROUTINE w1_load_parm(nfc,ierr)

    USE w1comm_parm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nfc
    INTEGER,INTENT(OUT) :: ierr
    INTEGER:: ist

    CALL w1_nlin(nfc,ist,ierr)
    IF(ierr.EQ.8) GOTO 9800
    IF(ierr.EQ.9) GOTO 9900
    ierr=0
    WRITE(6,'(A,I5)') '## w1_load_parm: w1.namelist loaded.'
    RETURN

9800 CONTINUE
    ierr=1
    WRITE(6,'(A,I5)') 'XX w1_load_parm: read error: ist=',ist
    RETURN

9900 CONTINUE
    ierr=2
    WRITE(6,'(A)') 'XX w1_load_parm: file open error:'
    RETURN
  END SUBROUTINE w1_load_parm
    
END MODULE w1parm
  
