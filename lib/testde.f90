! testde.f90

PROGRAM testde

  USE task_kinds,ONLY: dp
  USE libde
  USE libgrf
  
  IMPLICIT NONE

  REAL(dp):: a,b,h0,eps,cs,es
  INTEGER:: ilst
  
  a=0.D0
  b=1.D0
  h0=0.01D0
  eps=1.D-6
  ilst=0

1 CONTINUE

  WRITE(6,'(A,2ES12.4)') '## Input A and B: ',A,B
  READ(5,*,ERR=1,END=9000) A,B
  
  CALL DEFTAB(A,B,CS,ES,H0,EPS,ILST,f1,'f1')
  WRITE(6,'(A,ES16.8,ES12.4)') '## f1: integral and error:',CS,ES
  CALL DEFTAB(A,B,CS,ES,H0,EPS,ILST,f2,'f2')
  WRITE(6,'(A,ES16.8,ES12.4)') '## f2: integral and error:',CS,ES
  CALL DEFTAB(A,B,CS,ES,H0,EPS,ILST,f3,'f3')
  WRITE(6,'(A,ES16.8,ES12.4)') '## f3: integral and error:',CS,ES
  GO TO 1

9000 CONTINUE
  STOP

CONTAINS

  FUNCTION f1(y,by,ya)
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: y,by,ya
    REAL(dp):: f1
    REAL(dp):: dummy

    dummy=y
    dummy=by
    dummy=ya
    f1=1.D0
    RETURN
  END FUNCTION f1

  FUNCTION f2(y,by,ya)
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: y,by,ya
    REAL(dp):: f2
    REAL(dp):: dummy

    dummy=by
    dummy=ya
    f2=y
    RETURN
  END FUNCTION f2

  FUNCTION f3(y,by,ya)
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: y,by,ya
    REAL(dp):: f3
    REAL(dp):: dummy

    dummy=y
    f3=1.D0/SQRT(by)+1.D0/SQRT(ya)
    RETURN
  END FUNCTION f3

END PROGRAM testde
