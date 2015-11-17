!  ***** TASK/PIC PARAMETER *****

Module picparm

  PRIVATE
  PUBLIC pic_parm, pic_view,pic_broadcast

CONTAINS

!     ****** PARAMETER INPUT ******

  SUBROUTINE pic_parm(MODE,KIN,IERR)

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

1   CALL TASK_PARM(MODE,'PIC',kin,pic_nlin,pic_plst,IERR)
    IF(IERR.NE.0 .AND. IERR.NE.2) RETURN

    CALl pic_check(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
    IF(IERR.NE.0) IERR=IERR+100

    RETURN
  END SUBROUTINE pic_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE pic_nlin(NID,IST,IERR)

    USE piccomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NID
    INTEGER,INTENT(OUT) :: IST,IERR

    NAMELIST /PIC/ npx,npy,nx,ny,iend,nhmod, &
         me,mi,chrge,chrgi,te,ti,dt,eps,bxmin,bxmax,bymin,bymax,bzmin,bzmax,&
         c,omega,jxbg,jybg,jzbg,f,thetax,thetay,thetaz

    READ(NID,PIC,IOSTAT=IST,ERR=9800,END=9900)

    IERR=0
    RETURN

9800 IERR=8
    RETURN
9900 IERR=9
    RETURN
  END SUBROUTINE pic_nlin

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE pic_plst

    IMPLICIT NONE
    WRITE(6,'(A/)') '# &PIC : npx,npy,nx,ny,iend,nhmod,', '        &
         & me,mi,chrge,chrgi,te,ti,dt,eps,', '         bxbg,bybg,bzbg&
         &,c,omega,jxbg,jybg,jzbg,f', '         thetax,thetay,thetaz'
    RETURN

  END SUBROUTINE pic_plst

!     ****** CHECK INPUT PARAMETER ******

  SUBROUTINE pic_check(IERR)

    USE piccomm_parm
    IMPLICIT NONE
    INTEGER:: IERR

    IERR=0

    IF(dt <= 0.D0) THEN
       WRITE(6,'(A,1PE12.4)') 'XX pic_check: INVALID dt: dt=',dt
       IERR=1
    ENDIF
    RETURN
  END SUBROUTINE pic_check

!     ***** BROADCAST INPUT PARAMETERS *****

  SUBROUTINE pic_broadcast

    USE piccomm,ONLY: myid
    USE piccomm_parm
    USE libmpi
    IMPLICIT NONE
    integer:: idata(6)
    REAL(8):: ddata(23)

    IF(myid == 0) THEN
       idata(1)=npx
       idata(2)=npy
       idata(3)=nx
       idata(4)=ny
       idata(5)=iend
       idata(6)=nhmod
    END IF
    CALL mtx_broadcast_integer(idata,11)
       npx=idata(1)
       npy=idata(2)
       nx=idata(3)
       ny=idata(4)
       iend=idata(5)
       nhmod=idata(6)

    IF(myid == 0) THEN
       ddata(1)=me
       ddata(2)=mi
       ddata(3)=chrge
       ddata(4)=chrgi
       ddata(5)=te
       ddata(6)=ti
       ddata(7)=dt
       ddata(8)=eps
       ddata(9)=bxmin
       ddata(10)=bxmax
       ddata(11)=bymin
       ddata(12)=bymax
       ddata(13)=bzmin
       ddata(14)=bzmax
       ddata(15)=c
       ddata(16)=omega
       ddata(17)=jxbg
       ddata(18)=jybg
       ddata(19)=jzbg
       ddata(20)=f
       ddata(21)=thetax
       ddata(22)=thetay
       ddata(23)=thetaz
    END IF
    CALL mtx_broadcast_real8(ddata,8)
       me=ddata(1)
       mi=ddata(2)
       chrge=ddata(3)
       chrgi=ddata(4)
       te=ddata(5)
       ti=ddata(6)
       dt=ddata(7)
       eps=ddata(8)
       bxmin=ddata(9)
       bxmax=ddata(10)
       bymin=ddata(11)
       bymax=ddata(12)
       bzmin=ddata(13)
       bzmax=ddata(14)
       c=ddata(15)
       omega=ddata(16)
       jxbg=ddata(17)
       jybg=ddata(18)
       jzbg=ddata(19)
       f=ddata(20)
       thetax=ddata(21)
       thetay=ddata(22)
       thetaz=ddata(23)
    RETURN

  END SUBROUTINE pic_broadcast

!     ****** SHOW PARAMETERS ******

  SUBROUTINE pic_view

    use piccomm_parm
    implicit none

    WRITE(6,601) 'npx   ',npx   ,'npy   ',npy
    WRITE(6,601) 'nx    ',nx    ,'ny    ',ny
    WRITE(6,601) 'iend  ',iend  ,'nhmod ',nhmod
    WRITE(6,602) 'me    ',me    ,'mi    ',mi   , 'chrge ',chrge &
         ,'chrgi ',chrgi
    WRITE(6,602) 'te    ',te    ,'ti    ',ti    ,'dt    ',dt    ,'eps&
            ',eps
    WRITE(6,602) 'bxmin ',bxmin ,'bxmax ',bxmax ,'bymin ',bymin &
        ,'bymax ',bymax ,'bzmin ',bzmax &
        ,'c    ',c     ,'omega ',omega , 'jxbg  ',jxbg  ,'jybg  ',jybg  &
        ,'jzbg  ',jzbg  , 'f    ',f    , 'thetax',thetax,'thetay'&
        ,thetay,'thetaz',thetaz
    RETURN

601 FORMAT(' ',A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
            2X,A6,'=',I7,4X  :2X,A6,'=',I7)
602 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
            2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  END SUBROUTINE pic_view

END MODULE picparm

