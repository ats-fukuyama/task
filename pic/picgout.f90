!  ***** TASK/XX GRAPHIC OUTPU *****

Module picgout
  PRIVATE
  PUBLIC pic_gout
 
CONTAINS

  SUBROUTINE pic_gout

    USE piccomm
    USE picparm,ONLY: pic_parm
    USE grd2dp_mod,ONLY: grd2dp
    USE libgrf
    IMPLICIT NONE
    INTEGER:: kid,mode,err,ich1,ich2,nx,ny,np
    CHARACTER(LEN=2):: kch
    CHARACTER(LEN=80):: line
    REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: work
    REAL(rkind),ALLOCATABLE,DIMENSION(:):: x,y,vtote,vtoti
    REAL(rkind):: WORK_RGB(3,1)
    REAL(rkind):: aspect,vtotemax,vtotimax

    aspect=DBLE(nymax)/DBLE(nxmax)

    ALLOCATE(x(nxmax1),y(nymax1))
    ALLOCATE(vtote(npmax),vtoti(npmax))
    DO nx=1,nxmax1
       x(nx)=DBLE(nx-1)/DBLE(nxmax)
    END DO
    DO ny=1,nymax1
       y(ny)=DBLE(ny-1)/DBLE(nymax)
    END DO
    vtotemax=0.D0
    vtotimax=0.D0
    DO np=1,npmax
       vtote(np)=SQRT(vparae(np)**2+vperpe(np)**2)
       vtoti(np)=SQRT(vparai(np)**2+vperpi(np)**2)
       vtotemax=MAX(vtotemax,vtote(np))
       vtotimax=MAX(vtotimax,vtoti(np))
    END DO

1   CONTINUE
    err=0
    WRITE(6,'(A)') &
         '#### PIC GOUT: T1 E1-5 F1 X/exit'
    CALL TASK_KLIN(line,kid,mode,pic_parm)
    IF(mode == 2 .OR. mode == 3) GOTO 1

    ICH1=ICHAR(LINE(1:1))
    IF(ICH1.GE.97.AND.ICH1.LE.122) ICH1=ICH1-32
    ICH2=ICHAR(LINE(2:2))
    IF(ICH2.GE.97.AND.ICH2.LE.122) ICH2=ICH2-32
    KCH=CHAR(ICH1)//CHAR(ICH2)

    SELECT CASE(kch)
    CASE('T1')
       ALLOCATE(work(ntgmax,5))
       work(1:ntgmax,1)=akinet(1:ntgmax)
       work(1:ntgmax,2)=akinit(1:ntgmax)
       work(1:ntgmax,3)=aktott(1:ntgmax)
       work(1:ntgmax,4)=apott(1:ntgmax)
       work(1:ntgmax,5)=atott(1:ntgmax)
       CALL PAGES
       CALL GRD1D(0,timet,work,ntgmax,ntgmax,5,'@Ke,Ki,Ktot,Ptot,Wtot vs t@')
       CALL PAGEE
       DEALLOCATE(work)
    CASE('E1')
       CALL PAGES
       CALL GRD2D(1,x,y,ex,nxmax1,nxmax1,nymax1,'@Ex@',ASPECT=aspect)
       CALL GRD2D(2,x,y,ey,nxmax1,nxmax1,nymax1,'@Ey@',ASPECT=aspect)
       CALL GRD2D(3,x,y,rho,nxmax1,nxmax1,nymax1,'@rho@',ASPECT=aspect)
       CALL GRD2D(4,x,y,phi,nxmax1,nxmax1,nymax1,'@phi@',ASPECT=aspect)
       CALL PAGEE
    CASE('E2')
       CALL PAGES
       CALL GRD2D(1,x,y,bx,nxmax1,nxmax1,nymax1,'@Bx@',ASPECT=aspect)
       CALL GRD2D(2,x,y,by,nxmax1,nxmax1,nymax1,'@By@',ASPECT=aspect)
       CALL GRD2D(3,x,y,bz,nxmax1,nxmax1,nymax1,'@Bz@',ASPECT=aspect)
       CALL GRD2D(4,x,y,bb,nxmax1,nxmax1,nymax1,'@BB@',ASPECT=aspect)
       CALL PAGEE
    CASE('E3')
       CALL PAGES
       CALL GRD2D(1,x,y,Ax,nxmax1,nxmax1,nymax1,'@Ax@',ASPECT=aspect)
       CALL GRD2D(2,x,y,Ay,nxmax1,nxmax1,nymax1,'@Ay@',ASPECT=aspect)
       CALL GRD2D(3,x,y,Az,nxmax1,nxmax1,nymax1,'@Az@',ASPECT=aspect)
       CALL GRD2D(4,x,y,AA,nxmax1,nxmax1,nymax1,'@AA@',ASPECT=aspect)
       CALL PAGEE
    CASE('E4')
       CALL PAGES
       CALL GRD2D(1,x,y,jx,nxmax1,nxmax1,nymax1,'@jx@',ASPECT=aspect)
       CALL GRD2D(2,x,y,jy,nxmax1,nxmax1,nymax1,'@jy@',ASPECT=aspect)
       CALL GRD2D(3,x,y,jz,nxmax1,nxmax1,nymax1,'@jz@',ASPECT=aspect)
       CALL GRD2D(4,x,y,rho,nxmax1,nxmax1,nymax1,'@rho@',ASPECT=aspect)
       CALL PAGEE
    CASE('E5')
       CALL PAGES
       CALL GRD2D( 5,x,y,ex,nxmax1,nxmax1,nymax1,'@Ex@',ASPECT=aspect)
       CALL GRD2D( 6,x,y,ey,nxmax1,nxmax1,nymax1,'@Ey@',ASPECT=aspect)
       CALL GRD2D( 7,x,y,ez,nxmax1,nxmax1,nymax1,'@Ez@',ASPECT=aspect)
       CALL GRD2D( 8,x,y,esx,nxmax1,nxmax1,nymax1,'@ESx@',ASPECT=aspect)
       CALL GRD2D( 9,x,y,esy,nxmax1,nxmax1,nymax1,'@ESy@',ASPECT=aspect)
       CALL GRD2D(10,x,y,esz,nxmax1,nxmax1,nymax1,'@ESz@',ASPECT=aspect)
       CALL GRD2D(11,x,y,emx,nxmax1,nxmax1,nymax1,'@EMx@',ASPECT=aspect)
       CALL GRD2D(12,x,y,emy,nxmax1,nxmax1,nymax1,'@EMy@',ASPECT=aspect)
       CALL GRD2D(13,x,y,emz,nxmax1,nxmax1,nymax1,'@EMz@',ASPECT=aspect)
       CALL PAGEE
    CASE('F1')
       CALL PAGES
       WORK_RGB(1,1)=0.D0
       WORK_RGB(2,1)=0.D0
       WORK_RGB(3,1)=1.D0
       CALL GRD2DP(1,vparae,vperpe,vtote,npmax,'@vpara-vperp:e@',NMMAX=1, &
                   MARK_RGB=WORK_RGB,XMIN=-vtotemax,XMAX=vtotemax, &
                   YMIN=0.D0,YMAX=vtotemax,ASPECT=0.5D0)
       CALL GRD2DP(3,xe,ye,ze,npmax,'@x-y-z:e@',NMMAX=1, &
                   MARK_RGB=WORK_RGB,ASPECT=aspect)
       WORK_RGB(1,1)=1.D0
       WORK_RGB(2,1)=0.D0
       WORK_RGB(3,1)=0.D0
       CALL GRD2DP(2,vparai,vperpi,vtoti,npmax,'@vpara-vperp:i@',NMMAX=1, &
                   MARK_RGB=WORK_RGB,XMIN=-vtotimax,XMAX=vtotimax, &
                   YMIN=0.D0,YMAX=vtotimax,ASPECT=0.5D0)
       CALL GRD2DP(4,xi,yi,zi,npmax,'@x-y-z:i@',NMMAX=1, &
                   MARK_RGB=WORK_RGB,ASPECT=aspect)
       CALL PAGEE
    CASE('X') 
       GO TO 9000
    END SELECT
    GO TO 1

9000 CONTINUE
    DEALLOCATE(x,y)
    RETURN
  END SUBROUTINE pic_gout
END Module picgout
