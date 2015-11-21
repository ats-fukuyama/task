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
    INTEGER:: kid,mode,err,ich1,ich2,i
    CHARACTER(LEN=2):: kch
    CHARACTER(LEN=80):: line
    REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: work
    REAL(rkind),ALLOCATABLE,DIMENSION(:):: x,y,vtote,vtoti
    REAL(rkind):: WORK_RGB(3,1)

    ALLOCATE(x(nxmax+1),y(nymax+1))
    ALLOCATE(vtote(npmax),vtoti(npmax))
    DO i=1,nxmax+1
       x(i)=DBLE(i-1)/DBLE(nxmax)
    END DO
    DO i=1,nymax+1
       y(i)=DBLE(i-1)/DBLE(nymax)
    END DO
    DO i=1,npmax
       vtote(npmax)=SQRT(vparae(npmax)**2+vperpe(npmax)**2)
       vtoti(npmax)=SQRT(vparai(npmax)**2+vperpi(npmax)**2)
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
       ALLOCATE(work(ienemax,5))
       work(1:ienemax,1)=akinet(1:ienemax)
       work(1:ienemax,2)=akinit(1:ienemax)
       work(1:ienemax,3)=aktott(1:ienemax)
       work(1:ienemax,4)=apott(1:ienemax)
       work(1:ienemax,5)=atott(1:ienemax)
       CALL PAGES
       CALL GRD1D(0,timet,work,ienemax,ienemax,5,'@Ke,Ki,Ktot,Ptot,Wtot vs t@')
       CALL PAGEE
       DEALLOCATE(work)
    CASE('E1')
       CALL PAGES
       CALL GRD2D(1,x,y,ex,nxmax+1,nxmax+1,nymax+1,'@Ex@')
       CALL GRD2D(2,x,y,ey,nxmax+1,nxmax+1,nymax+1,'@Ey@')
       CALL GRD2D(3,x,y,rho,nxmax+1,nxmax+1,nymax+1,'@rho@')
       CALL GRD2D(4,x,y,phi,nxmax+1,nxmax+1,nymax+1,'@phi@')
       CALL PAGEE
    CASE('E2')
       CALL PAGES
       CALL GRD2D(1,x,y,bx,nxmax+1,nxmax+1,nymax+1,'@Bx@')
       CALL GRD2D(2,x,y,by,nxmax+1,nxmax+1,nymax+1,'@By@')
       CALL GRD2D(3,x,y,bz,nxmax+1,nxmax+1,nymax+1,'@Bz@')
       CALL GRD2D(4,x,y,bb,nxmax+1,nxmax+1,nymax+1,'@BB@')
       CALL PAGEE
    CASE('E3')
       CALL PAGES
       CALL GRD2D(1,x,y,Ax,nxmax+1,nxmax+1,nymax+1,'@Ax@')
       CALL GRD2D(2,x,y,Ay,nxmax+1,nxmax+1,nymax+1,'@Ay@')
       CALL GRD2D(3,x,y,Az,nxmax+1,nxmax+1,nymax+1,'@Az@')
       CALL GRD2D(4,x,y,AA,nxmax+1,nxmax+1,nymax+1,'@AA@')
       CALL PAGEE
    CASE('E4')
       CALL PAGES
       CALL GRD2D(1,x,y,jx,nxmax+1,nxmax+1,nymax+1,'@jx@')
       CALL GRD2D(2,x,y,jy,nxmax+1,nxmax+1,nymax+1,'@jy@')
       CALL GRD2D(3,x,y,jz,nxmax+1,nxmax+1,nymax+1,'@jz@')
       CALL GRD2D(4,x,y,rho,nxmax+1,nxmax+1,nymax+1,'@rho@')
       CALL PAGEE
    CASE('E5')
       CALL PAGES
       CALL GRD2D( 5,x,y,ex,nxmax+1,nxmax+1,nymax+1,'@Ex@')
       CALL GRD2D( 6,x,y,ey,nxmax+1,nxmax+1,nymax+1,'@Ey@')
       CALL GRD2D( 7,x,y,ez,nxmax+1,nxmax+1,nymax+1,'@Ez@')
       CALL GRD2D( 8,x,y,esx,nxmax+1,nxmax+1,nymax+1,'@ESx@')
       CALL GRD2D( 9,x,y,esy,nxmax+1,nxmax+1,nymax+1,'@ESy@')
       CALL GRD2D(10,x,y,esz,nxmax+1,nxmax+1,nymax+1,'@ESz@')
       CALL GRD2D(11,x,y,emx,nxmax+1,nxmax+1,nymax+1,'@EMx@')
       CALL GRD2D(12,x,y,emy,nxmax+1,nxmax+1,nymax+1,'@EMy@')
       CALL GRD2D(13,x,y,emz,nxmax+1,nxmax+1,nymax+1,'@EMz@')
       CALL PAGEE
    CASE('F1')
       CALL PAGES
       WORK_RGB(1,1)=0.D0
       WORK_RGB(2,1)=0.D0
       WORK_RGB(3,1)=1.D0
       CALL GRD2DP(1,vparae,vperpe,vtote,npmax,'@vpara-vperp:e@',NMMAX=1, &
                   MARK_RGB=WORK_RGB)
       CALL GRD2DP(3,xe,ye,ze,npmax,'@x-y-z:e@',NMMAX=1, &
                   MARK_RGB=WORK_RGB)
       WORK_RGB(1,1)=1.D0
       WORK_RGB(2,1)=0.D0
       WORK_RGB(3,1)=0.D0
       CALL GRD2DP(2,vparai,vperpi,vtoti,npmax,'@vpara-vperp:i@',NMMAX=1, &
                   MARK_RGB=WORK_RGB)
       CALL GRD2DP(4,xi,yi,zi,npmax,'@x-y-z:i@',NMMAX=1, &
                   MARK_RGB=WORK_RGB)
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
