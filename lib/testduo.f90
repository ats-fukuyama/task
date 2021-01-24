! testduo.f90

PROGRAM testduo

  USE task_kinds,ONLY: dp
  USE libderiv
  IMPLICIT NONE
  REAL(dp):: x(2),f(2)
  REAL(dp):: x0,f0,df0,ddf0
  INTEGER:: id,ierr

1 CONTINUE
  WRITE(6,'(A)') '## Input: x1,x2'
  READ(5,*,ERR=1,END=9000) x(1),x(2)
2 CONTINUE
  WRITE(6,'(A)') '## Input: f1,f2'
  READ(5,*,ERR=2,END=1) f(1),f(2)
3 CONTINUE
  WRITE(6,'(A)') '## Input: id,x0'
  READ(5,*,ERR=3,END=2) id,x0

  CALL duo_deriv(x,f,x0,f0,df0,ddf0,id,ierr)

  IF(ierr.NE.0) GO TO 3
  WRITE(6,'(A,4ES12.4)') 'x0,f0,df0,ddf0=',x0,f0,df0,ddf0
  WRITE(6,'(A,4ES12.4)') '   f1,f0+df0* =', &
       f(1),f0+df0*(x(1)-x0)+0.5D0*ddf0*(x(1)-x0)**2
  WRITE(6,'(A,4ES12.4)') '   f2,f0+df0* =', &
       f(2),f0+df0*(x(2)-x0)+0.5D0*ddf0*(x(2)-x0)**2
  GO TO 3

9000 CONTINUE
  STOP
END PROGRAM testduo
  
