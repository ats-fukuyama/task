!  ***** TASK/PIC PARAMETER *****

Module picparm

  PRIVATE
  PUBLIC pic_parm,pic_view,pic_broadcast

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

    NAMELIST /PIC/ npxmax,npymax,nxmax,nymax,ntmax,ntstep,ntgstep,ntpstep, &
         me,mi,chrge,chrgi,te,ti,dt,eps,bxmin,bxmax,bymin,bymax,bzmin,bzmax,&
         vcfact,omega,jxant,jyant,jzant,phxant,phyant,phzant, &
         model_boundary,model_antenna,model_wg, &
         xmin_wg,xmax_wg,ymin_wg,ymax_wg,amp_wg, ph_wg,rot_wg,eli_wg

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
    WRITE(6,'(A/)') &
         '# &PIC : npxmax,npymax,nxmax,nymax,ntmax,ntstep,ntgstep,ntpstep,', &
         '         me,mi,chrge,chrgi,te,ti,dt,eps,vcfact,bxbg,bybg,bzbg,', &
         '         omega,jxant,jyant,jzant,phxant,phyant,phzant', &
         '         xmin_wg,xmax_wg,ymin_wg,ymax_wg,amp_wg, ph_wg,rot_wg,eli_wg'
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
    integer,parameter:: nint=11
    integer,parameter:: ndbl=30
    integer:: idata(nint)
    REAL(8):: ddata(ndbl)

    IF(myid == 0) THEN
       idata(1)=npxmax
       idata(2)=npymax
       idata(3)=nxmax
       idata(4)=nymax
       idata(5)=ntmax
       idata(6)=ntstep
       idata(7)=ntgstep
       idata(8)=ntpstep
       idata(9)=model_boundary
       idata(10)=model_antenna
       idata(11)=model_wg
    END IF
    CALL mtx_broadcast_integer(idata,nint)
       npxmax=idata(1)
       npymax=idata(2)
       nxmax=idata(3)
       nymax=idata(4)
       ntmax=idata(5)
       ntstep=idata(6)
       ntgstep=idata(7)
       ntpstep=idata(8)
       model_boundary=idata(9)
       model_antenna=idata(10)
       model_wg=idata(11)

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
       ddata(15)=vcfact
       ddata(16)=omega
       ddata(17)=jxant
       ddata(18)=jyant
       ddata(19)=jzant
       ddata(20)=phxant
       ddata(21)=phyant
       ddata(22)=phzant
       ddata(23)=xmin_wg
       ddata(24)=xmax_wg
       ddata(25)=ymin_wg
       ddata(26)=ymax_wg
       ddata(27)=amp_wg
       ddata(28)=ph_wg
       ddata(29)=rot_wg
       ddata(30)=eli_wg
    END IF
    CALL mtx_broadcast_real8(ddata,ndbl)
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
       vcfact=ddata(15)
       omega=ddata(16)
       jxant=ddata(17)
       jyant=ddata(18)
       jzant=ddata(19)
       phxant=ddata(20)
       phyant=ddata(21)
       phzant=ddata(22)
       xmin_wg=ddata(23)
       xmax_wg=ddata(24)
       ymin_wg=ddata(25)
       ymax_wg=ddata(26)
       amp_wg=ddata(27)
       ph_wg=ddata(28)
       rot_wg=ddata(29)
       eli_wg=ddata(30)
    RETURN

  END SUBROUTINE pic_broadcast

!     ****** SHOW PARAMETERS ******

  SUBROUTINE pic_view

    use piccomm_parm
    implicit none

    WRITE(6,601) 'npxmax',npxmax,'npymax',npymax
    WRITE(6,601) 'nxmax ',nxmax ,'nymax ',nymax
    WRITE(6,601) 'ntmax ',ntmax ,'ntstep',ntstep
    WRITE(6,611) 'ntgstep   ',ntgstep,'ntpstep   ',ntpstep
    WRITE(6,602) 'me    ',me    ,'mi    ',mi    , &
                 'chrge ',chrge ,'chrgi ',chrgi
    WRITE(6,602) 'te    ',te    ,'ti    ',ti    , &
                 'dt    ',dt    ,'eps   ',eps
    WRITE(6,602) 'bxmin ',bxmin ,'bxmax ',bxmax , &
                 'bymin ',bymin ,'bymax ',bymax , &
                 'bzmin ',bzmax
    WRITE(6,602) 'vcfact',vcfact,'omega ',omega 
    WRITE(6,602) 'jxant ',jxant ,'jyant ',jyant , &
                 'jzant ',jzant , &
                 'phxant',phxant,'phyant',phyant, &
                 'phzant',phzant
    WRITE(6,611) 'model_boundary ',model_boundary
    WRITE(6,611) 'model_antenna ',model_antenna
    WRITE(6,611) 'model_wg ',model_wg
    WRITE(6,602) 'xmin_wg ',xmin_wg,'xmax_wg ',xmax_wg,&
                 'ymin_wg ',ymin_wg,'ymax_wg ',ymax_wg,'amp_wg ',amp_wg,&
                 'ph_wg ',ph_wg,'rot_wg ',rot_wg,'eli_wg ',eli_wg
    RETURN

601 FORMAT(' ',A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
            2X,A6,'=',I7,4X  :2X,A6,'=',I7)
602 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
            2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
611 FORMAT(' ',A12,'=',I7,4X  :2X,A12,'=',I7,4X  : &
            2X,A12,'=',I7)
  END SUBROUTINE pic_view

END MODULE picparm

