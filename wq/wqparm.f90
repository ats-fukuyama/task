! wqparm.f90

MODULE wqparm

  PRIVATE
  PUBLIC wq_parm,wq_broadcast

CONTAINS

!     ****** INPUT PARAMETERS ******

  SUBROUTINE wq_parm(mode,kin,ierr)
    
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
    CHARACTER(LEN=*),INTENT(IN):: kin
    INTEGER,INTENT(OUT):: ierr
    
1   CALL TASK_PARM(mode,'WQ',kin,wq_nlin,wq_plst,ierr)
    IF(ierr.NE.0) RETURN

    CALl wq_check(ierr)
    IF(mode.EQ.0.AND.ierr.NE.0) GOTO 1
    IF(ierr.NE.0) ierr=ierr+100

    RETURN
  END SUBROUTINE wq_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE wq_nlin(nid,ist,ierr)

    USE wqcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nid
    INTEGER,INTENT(OUT):: ist,ierr

    NAMELIST /WQ/ &
         FREQ,dtfactor,dxfactor,dyfactor,nufactor, &
         B0,RR,RA,q0,qa,n0,ntmax,nxmax,nymax,INMODE,TMN

    ierr=0

    READ(nid,WQ,IOSTAT=ist,ERR=9800,END=9900)
    IF(ist.NE.0) THEN
       WRITE(6,'(A,I5)') 'XX wq_nlin: READ ERROR: IOSTAT=',ist
       ierr=ist
       RETURN
    END IF
    RETURN

9800 IERR=8
    WRITE(6,'(A)') 'XX wq_nlin: READ ERROR'
    RETURN
9900 IERR=9
    WRITE(6,'(A)') 'XX wq_nlin: READ END OF FILE'
    RETURN
  END SUBROUTINE wq_nlin

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE wq_plst

    WRITE(6,'(A)') &
         '# &WQ: FREQ,dtfactor,dxfactor,dyfactor,nufactor,', &
         '       B0,RR,RA,q0,qa,n0,ntmax,nxmax,nymax,INMODE,TMN'
    RETURN
  END SUBROUTINE wq_plst

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE wq_check(ierr)

    USE wqcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    ierr=0

    RETURN
  END SUBROUTINE wq_check

!     ***** BROADCAST INPUT PARAMETERS *****

  SUBROUTINE wq_broadcast

    USE wqcomm_parm
    USE libmpi
    IMPLICIT NONE
    INTEGER,DIMENSION(99):: idata
    REAL(rkind),DIMENSION(99):: rdata

! --- WQ specific input parameters ---

    idata( 1)=ntmax
    idata( 2)=nxmax
    idata( 3)=nymax
    idata( 4)=INMODE

    CALL mtx_broadcast_integer(idata,4)

    ntmax=idata( 1)
    nxmax=idata( 2)
    nymax=idata( 3)
    INMODE=idata( 4)

    rdata( 1)=FREQ
    rdata( 2)=dtfactor
    rdata( 3)=dxfactor
    rdata( 4)=dyfactor
    rdata( 5)=nufactor
    rdata( 6)=B0
    rdata( 7)=RR
    rdata( 8)=RA
    rdata( 9)=q0
    rdata(10)=qa
    rdata(11)=n0
    rdata(12)=TMN
    
    CALL mtx_broadcast_real8(rdata,12)

    FREQ=rdata( 1)
    dtfactor=rdata( 2)
    dxfactor=rdata( 3)
    dyfactor=rdata( 4)
    nufactor=rdata( 5)
    B0=rdata( 6)
    RR=rdata( 7)
    RA=rdata( 8)
    q0=rdata( 9)
    qa=rdata(10)
    n0=rdata(11)
    TMN=rdata(12)

  END SUBROUTINE wq_broadcast
END MODULE wqparm
