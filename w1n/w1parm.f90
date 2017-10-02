Module w1parm

  PRIVATE

  PUBLIC w1_parm, w1_view

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

    NAMELIST /W1/ BB,RR,RZ,RA,RD,RB,RF,WALLR,APRFPN,APRFTR,APRFTP, &
                  RKZ,DRF,DRKZ,DXFACT,DXWDTH, &
                  AJYH,ALYH,APYH,AJYL,ALYL,APYL,AJZH,AJZL,NAMAX, &
                  PA,PZ,PN,PTPP,PTPR,PU,PNS,PTS,PZCL,NSMAX, &
                  NXPMAX,NXVMAX,NZPMAX,NPRINT,NFILE,NGRAPH,NLOOP,NSYM, &
                  NMODEL,NALPHA,NDMAX,XDMAX,IHARM,NSYS,NDISP, &
                  EPSH,ZEFF,WVYSIZ,NCDTYP,NXABS,IELEC

    READ(NID,W1,IOSTAT=IST,ERR=9800,END=9900)

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
         '        RKZ,DRF,DRKZ,DXFACT,DXWDTH, ',&
         '        AJYH,ALYH,APYH,AJYL,ALYL,APYL,AJZH,AJZL,NAMAX,', &
         '        PA,PZ,PN,PTPP,PTPR,PU,PNS,PTS,PZCL,NSMAX,', &
         '        NXPMAX,NXVMAX,NZPMAX,NPRINT,NFILE,NGRAPH,NLOOP,NSYM,', &
         '        NMODEL,NALPHA,NDMAX,XDMAX,IHARM,NSYS,NDISP,', &
         '        EPSH,ZEFF,WVYSIZ,NCDTYP,NXABS,IELEC'
    RETURN

  END SUBROUTINE w1_plst

!     ****** CHECK INPUT PARAMETER ******

  SUBROUTINE w1_check(IERR)

    USE w1comm_parm
    IMPLICIT NONE
    INTEGER:: IERR

    IERR=0

    IF(NXPMAX < 0) THEN
       WRITE(6,'(A,I8)') 'W1 w1_check: INVALID nxmax: nxmax=',nxpmax
       IERR=1
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
                 'WVYSIZ',WVYSIZ,'XDMAX ',XDMAX

    WRITE(6,'(A)') '     ', &
          '  NS/IELEC   PA/PZ      PN/PNS    PTPR/PTS     PTPP/PU    PZCL'
    DO NS=1,NSMAX
       WRITE(6,'(A,I6,1P5E12.4)') &
            '     ',NS,PA(NS),PN(NS),PTPR(NS),PTPP(NS),PZCL(NS)
       WRITE(6,'(A,I6,1P4E12.4)') &
            '     ',IELEC(NS),PZ(NS),PNS(NS),PTS(NS),PU(NS)
    END DO

    WRITE(6,'(A)') '     ', &
          '  NA     AJYH/AJYL   AJZH/AJZL   ALYH/ALYL   APYH/APYL'
    DO NA=1,NAMAX
       WRITE(6,'(A,I6,1P4E12.4)') &
            'ANT-H',NA,AJYH(NA),AJZH(NA),ALYH(NA),APYH(NA)
       WRITE(6,'(A,I6,1P4E12.4)') &
            'ANT-L',NA,AJYL(NA),AJZL(NA),ALYL(NA),APYL(NA)
    END DO

    WRITE(6,602) 'NXPMAX',NXPMAX,'NXVMAX',NXVMAX, &
                 'NZPMAX',NZPMAX,'NSMAX ',NSMAX, &
                 'NAMAX ',NAMAX ,'NDMAX ',NDMAX, &
                 'NPRINT',NPRINT,'NFILE ',NFILE, &
                 'NGRAPH',NGRAPH,'NLOOP ',NLOOP, &
                 'NSYM  ',NSYM  ,'NMODEL',NMODEL, &
                 'NALPHA',NALPHA,'NHARM ',NHARM, &
                 'NSYS  ',NSYS  ,'NDISP ',NDISP, &
                 'NCDTYP',NCDTYP,'NXABS ',NXABS
    RETURN

601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
            2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
602 FORMAT(' ',A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
            2X,A6,'=',I7,4X  :2X,A6,'=',I7)
  END SUBROUTINE w1_view

END MODULE w1parm
  
