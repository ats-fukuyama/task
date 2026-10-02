! wftest.f90

MODULE wftest
  PRIVATE
  PUBLIC wf_test

CONTAINS

  SUBROUTINE wf_test
    USE wfcomm
    USE feminterpolate
    USE libgrf
    IMPLICIT NONE
    INTEGER:: init=0
    INTEGER:: mode_test,nxmax,nymax
    REAL(rkind):: xmin,xmax,ymin,ymax
    REAL(rkind),ALLOCATABLE:: xg(:),yg(:),fg(:,:)
    INTEGER:: nx,ny,nelm
    REAL(rkind):: delx,dely,x,y
    COMPLEX(rkind):: ce
    CHARACTER(LEN=80):: title

    IF(init.EQ.0) THEN
       mode_test=1
       xmin=bdrmin
       xmax=bdrmax
       ymin=bdzmin
       ymax=bdzmax
       nxmax=101
       nymax=101
       init=1
    END IF
       
1   CONTINUE
    WRITE(6,'(A)') &
         '## INPUT: mode_test(1-6,9),xmin,xmax,ymin,ymax,nxmax,nymax:'
    WRITE(6,'(11X,I4,4ES12.4,2I6)') mode_test,xmin,xmax,ymin,ymax,nxmax,nymax
    READ(5,*,ERR=1,END=9000) mode_test,xmin,xmax,ymin,ymax,nxmax,nymax
    IF(mode_test.EQ.9) GOTO 9000
    IF(mode_test.LT.1.OR.mode_test.GT.6) GOTO 1

    ALLOCATE(xg(nxmax),yg(nymax),fg(nxmax,nymax))

    delx=(xmax-xmin)/(nxmax-1)
    dely=(ymax-ymin)/(nymax-1)
    DO ny=1,nymax
       yg(ny)=ymin+dely*(ny-1)
    END DO
    DO nx=1,nxmax
       xg(nx)=xmin+delx*(nx-1)
    END DO
    
    DO ny=1,nymax
       y=yg(ny)
       DO nx=1,nxmax
          x=xg(nx)
          CALL fem_find_nelm_for_xy(x,y,nelm)
          SELECT CASE(mode_test)
          CASE(1,2)
             CALL FIELDCR(nelm,x,y,cesd,ce)
          CASE(3,4)
             CALL FIELDCZ(nelm,x,y,cesd,ce)
          CASE(5,6)
             CALL FIELDCP(nelm,x,y,cend,ce)
          END SELECT
          SELECT CASE(mode_test)
          CASE(1,3,5)
             fg(nx,ny)=REAL(ce)
          CASE(2,4,6)
             fg(nx,ny)=AIMAG(ce)
          END SELECT
          WRITE(24,'(A,2I6,3ES12.4)') 'test:',nx,ny,xg(nx),yg(ny),fg(nx,ny)
       END DO
    END DO

    CALL PAGES
    SELECT CASE(mode_test)
    CASE(1)
       title='@Ex real(x,y)@'
    CASE(2)
       title='@Ex imag(x,y)@'
    CASE(3)
       title='@Ey real(x,y)@'
    CASE(4)
       title='@Ey imag(x,y)@'
    CASE(5)
       title='@Ez real(x,y)@'
    CASE(6)
       title='@Ez imag(x,y)@'
    END SELECT
    CALL GRD2D(0,xg,yg,fg,nxmax,nxmax,nymax,title)
    CALL PAGEE

    DEALLOCATE(xg,yg,fg)
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wf_test
END MODULE wftest
    
