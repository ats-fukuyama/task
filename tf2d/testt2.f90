!   $Id$

PROGRAM testt2

  USE bpsd,ONLY: rkind,pi
  USE LIBT2,ONLY: integ_f,func_1,x_ab
  IMPLICIT NONE
  REAL(rkind):: xmax,dx
  INTEGER:: nxmax,nx
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: x,f,g,err

  xmax=10.D0
  nxmax=101

1 CONTINUE
  
  WRITE(6,'(A,1PE12.4,I10)') '## INPUT xmax,nxmax = ',xmax,nxmax
  READ(5,*,ERR=1,END=9000) xmax,nxmax

  ALLOCATE(x(nxmax),f(nxmax),g(nxmax),err(nxmax))

  dx=xmax/(nxmax-1)
  DO nx=1,nxmax
     x(nx)=dx*(nx-1)
     x_ab=x(nx)
     IF(x_ab.LE.0.D0) THEN
        f(nx)=4.D0/(3.D0*SQRT(pi))
        err(nx)=0.D0
        g(nx)=4.D0/(3.D0*SQRT(pi))
     ELSE
        CALL integ_f(FUNC_1,f(nx),err(nx),EPS=1.D-8,ILST=1,KID='FUNC_1')
        g(nx)=(4.D0/(3.D0*SQRT(pi))) &
             *(SQRT(1.D0+x_ab**2)+x_ab**2*LOG(x_ab/(1.D0+SQRT(1+x_ab**2))))
     END IF
  END DO

  WRITE(6,'(1P5E12.4)') (x(nx),f(nx),g(nx),f(nx)-g(nx),err(nx),nx=1,nxmax)

  DEALLOCATE(x,f,g,err)

  GO TO 1

9000 CONTINUE
  STOP
END PROGRAM testt2
