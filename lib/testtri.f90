! testtri.f90

PROGRAM testtri

  USE task_kinds,ONLY: dp
  USE libderiv
  IMPLICIT NONE
  REAL(dp):: x(3),f(3)
  REAL(dp):: x0,f0,df0,ddf0
  INTEGER:: ierr

1 CONTINUE
  WRITE(6,'(A)') '## Input: x1,x2,x3'
  READ(5,*,ERR=1,END=9000) x(1),x(2),x(3)
2 CONTINUE
  WRITE(6,'(A)') '## Input: f1,f2,f3'
  READ(5,*,ERR=2,END=1) f(1),f(2),f(3)
3 CONTINUE
  WRITE(6,'(A)') '## Input: x0'
  READ(5,*,ERR=3,END=2) x0

  CALL tri_deriv(x,f,x0,f0,df0,ddf0,ierr)

  IF(ierr.NE.0) GO TO 3
  WRITE(6,'(A,4ES12.4)') 'x0,f0,df0,ddf0=',x0,f0,df0,ddf0
  WRITE(6,'(A,4ES12.4)') '   f1,f0+df0* =', &
       f(1),f0+df0*(x(1)-x0)+0.5D0*ddf0*(x(1)-x0)**2
  WRITE(6,'(A,4ES12.4)') '   f2,f0+df0* =', &
       f(2),f0+df0*(x(2)-x0)+0.5D0*ddf0*(x(2)-x0)**2
  WRITE(6,'(A,4ES12.4)') '   f3,f0+df0* =', &
       f(3),f0+df0*(x(3)-x0)+0.5D0*ddf0*(x(3)-x0)**2
  GO TO 3

9000 CONTINUE
  STOP
END PROGRAM testtri
  
