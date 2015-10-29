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

    ALLOCATE(x(nx+1),y(ny+1))
    ALLOCATE(vtote(np),vtoti(np))
    DO i=1,nx+1
       x(i)=DBLE(i-1)/DBLE(nx)
    END DO
    DO i=1,ny+1
       y(i)=DBLE(i-1)/DBLE(ny)
    END DO
    DO i=1,np
       vtote(np)=SQRT(vparae(np)**2+vperpe(np)**2)
       vtoti(np)=SQRT(vparai(np)**2+vperpi(np)**2)
    END DO

1   CONTINUE
    err=0
    WRITE(6,'(A)') &
         '#### PIC GOUT: T1 E1 E2 F1 X/exit'
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
       CALL GRD2D(1,x,y,ex,nx+1,nx+1,ny+1,'@Ex@')
       CALL GRD2D(2,x,y,ey,nx+1,nx+1,ny+1,'@Ey@')
       CALL GRD2D(3,x,y,rho,nx+1,nx+1,ny+1,'@rho@')
       CALL GRD2D(4,x,y,phi,nx+1,nx+1,ny+1,'@phi@')
       CALL PAGEE
    CASE('E2')
       CALL PAGES
       CALL GRD2D(1,x,y,Ax,nx+1,nx+1,ny+1,'@Ax@')
       CALL GRD2D(2,x,y,Ay,nx+1,nx+1,ny+1,'@Ay@')
       CALL GRD2D(3,x,y,Az,nx+1,nx+1,ny+1,'@Az@')
       CALL GRD2D(4,x,y,phi,nx+1,nx+1,ny+1,'@phi@')
       CALL PAGEE
    CASE('F1')
       CALL PAGES
       WORK_RGB(1,1)=0.D0
       WORK_RGB(2,1)=0.D0
       WORK_RGB(3,1)=1.D0
       CALL GRD2DP(1,vparae,vperpe,vtote,np,'@vpara-vperp:e@',NMMAX=1, &
                   MARK_RGB=WORK_RGB)
       CALL GRD2DP(3,xe,ye,ze,np,'@x-y-z:e@',NMMAX=1, &
                   MARK_RGB=WORK_RGB)
       WORK_RGB(1,1)=1.D0
       WORK_RGB(2,1)=0.D0
       WORK_RGB(3,1)=0.D0
       CALL GRD2DP(2,vparai,vperpi,vtoti,np,'@vpara-vperp:i@',NMMAX=1, &
                   MARK_RGB=WORK_RGB)
       CALL GRD2DP(4,xi,yi,zi,np,'@x-y-z:i@',NMMAX=1, &
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
