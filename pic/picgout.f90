!  ***** TASK/XX GRAPHIC OUTPU *****

Module picgout
  PRIVATE
  PUBLIC pic_gout
 
CONTAINS

  SUBROUTINE pic_gout

    USE piccomm
    USE picparm,ONLY: pic_parm
    USE libgrf
    IMPLICIT NONE
    INTEGER:: kid,mode,err,ich1,ich2,i
    CHARACTER(LEN=2):: kch
    CHARACTER(LEN=80):: line
    REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: work
    REAL(rkind),ALLOCATABLE,DIMENSION(:):: x,y

    ALLOCATE(x(nx+1),y(ny+1))
    DO i=1,nx+1
       x(i)=DBLE(i-1)/DBLE(nx)
    END DO
    DO i=1,ny+1
       y(i)=DBLE(i-1)/DBLE(ny)
    END DO

1   CONTINUE
    err=0
    WRITE(6,'(A)') &
         '#### PIC GOUT: T1 E1 F1 X/exit'
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
    CASE('X') 
       GO TO 9000
    END SELECT
    GO TO 1

9000 CONTINUE
    DEALLOCATE(x,y)
    RETURN
  END SUBROUTINE pic_gout
END Module picgout
