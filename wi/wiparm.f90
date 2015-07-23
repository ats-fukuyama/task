!     $Id$

Module wiparm

  PRIVATE

  PUBLIC wi_parm, wi_view

CONTAINS

!     ****** PARAMETER INPUT ******

  SUBROUTINE wi_parm(MODE,KIN,IERR)

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

1   CALL TASK_PARM(MODE,'WI',KIN,wi_nlin,wi_plst,IERR)
    IF(IERR.NE.0 .AND. IERR.NE.2) RETURN

    CALl wi_check(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
    IF(IERR.NE.0) IERR=IERR+100

    RETURN
  END SUBROUTINE wi_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE wi_nlin(NID,IST,IERR)

    USE wicomm, ONLY: modelg,xmin,xmax,dx0,xwint,pn0,alfa,any,beta,cfyn, &
                      ntaumax,taumin,taumax,nalfamax,alfamin,alfamax, &
                      modelp,pnu,dxmin,xwmin, &
                      modewi,kfscan,idebug,rkind

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NID
    INTEGER,INTENT(OUT) :: IST,IERR
    REAL(rkind),SAVE:: pn0_save

    NAMELIST /WI/ modelg,xmin,xmax,dx0,xwint,pn0,alfa,any,beta,cfyn, &
                  ntaumax,taumin,taumax,nalfamax,alfamin,alfamax, &
                  modelp,pnu,dxmin,xwmin, &
                  modewi,kfscan,idebug

    IF(modewi.EQ.1) THEN
       pn0=pn0_save
       alfa=alfa/beta
       xmin=xmin*beta
       xmax=xmax*beta
       dx0=dx0*beta
       dxmin=dxmin*beta
       xwmin=xwmin*beta
    ENDIF

    READ(NID,WI,IOSTAT=IST,ERR=9800,END=9900)

    IF(modewi.EQ.1) THEN
       pn0=1.D0
       alfa=alfa*beta
       xmin=xmin/beta
       xmax=xmax/beta
       dx0=dx0/beta
       dxmin=dxmin/beta
       xwmin=xwmin/beta
    ENDIF

    IERR=0
    RETURN

9800 IERR=8
    RETURN
9900 IERR=9
    RETURN
  END SUBROUTINE wi_nlin

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE wi_plst

    IMPLICIT NONE
    WRITE(6,'(A)') '# &WI : modelg,xmin,xmax,dx0,xwint,pn0,alfa,any,beta,cfyn,'
    WRITE(6,'(A)') '        ntaumax,taumin,taumax,nalfamax,alfamin,alfamax,'
    WRITE(6,'(A)') '        modelp,pnu,dxmin,xwmin,'
    WRITE(6,'(A)') '        modewi,kfscan,idebug'
    RETURN

  END SUBROUTINE wi_plst

!     ****** CHECK INPUT PARAMETER ******

  SUBROUTINE wi_check(IERR)

    USE wicomm,ONLY: modelg,modelp,xmin,xmax,taumin,taumax
    IMPLICIT NONE
    INTEGER:: IERR

    IERR=0

    IF(modelg /= 0) THEN
       WRITE(6,'(A,I6)') 'XX wi_check: INVALID modelg: modelg=',modelg
       IERR=1
    ENDIF
    IF(modelp < 0 .OR. modelp > 2) THEN
       WRITE(6,'(A,I6)') 'XX wi_check: INVALID modelp: modelp=',modelp
       IERR=1
    ENDIF
    IF(xmax <= xmin) THEN
       WRITE(6,'(A,1PE12.4)') 'XX wi_check: INVALID xmax: xmax=',xmax
       IERR=1
    ENDIF
    IF(taumax <= taumin) THEN
       WRITE(6,'(A,A,1P2E12.4)') 'XX wi_check: INVALID taumin,taumax: ',&
                                 'taumin,taumax =',taumin,taumax
       IERR=1
    ENDIF
    RETURN
  END SUBROUTINE wi_check

!     ****** SHOW PARAMETERS ******

  SUBROUTINE wi_view

    use wicomm
    implicit none

    WRITE(6,602) 'modelg  ',modelg, 'modelp  ',modelp, &
                 'ntaumax ',ntaumax,'nalfamax',nalfamax
    WRITE(6,602) 'modewi  ',modewi, 'idebug  ',idebug
    IF(modewi.eq.0) THEN
       WRITE(6,601) 'xmin  ',xmin  ,'xmax  ',xmax  , &
                    'dx0   ',dx0   ,'alfa  ',alfa
       WRITE(6,601) 'dxmin ',dxmin ,'xwmin ',xwmin
    ELSE
       WRITE(6,601) 'xmin  ',xmin*beta  ,'xmax  ',xmax*beta  , &
                    'dx0   ',dx0*beta   ,'alfa  ',alfa/beta
       WRITE(6,601) 'dxmin ',dxmin*beta ,'xwmin ',xwmin*beta
    END IF
       
    WRITE(6,601) 'pn0   ',pn0   ,'any   ',any, &
                 'xwint ',xwint 
    WRITE(6,601) 'beta  ',beta  ,'pnu   ',pnu, &
                 'taumin',taumin,'taumax',taumax
    WRITE(6,603) 'cfyn  ',cfyn
    RETURN

601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
            2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
602 FORMAT(' ',A8,'=',I7,2X  :2X,A8,'=',I7,2X  : &
            2X,A8,'=',I7,2X  :2X,A8,'=',I7)
603 FORMAT(' ',A6,'=',1P2E11.3:2X,A6,'=',1P2E11.3)

  END SUBROUTINE wi_view

END MODULE wiparm
  
