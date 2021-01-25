! testrk.f90

MODULE testrk_local

  USE task_kinds,ONLY: dp
  INTEGER,PARAMETER:: neqmax=1
  INTEGER:: model
  REAL(dp):: parm1,parm2

CONTAINS

  SUBROUTINE SUB(X,Y,F)
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: X,Y(neqmax)
    REAL(dp),INTENT(OUT):: F(neqmax)

    F(1)=Y(1)
    RETURN
  END SUBROUTINE SUB

END MODULE testrk_local

PROGRAM testrk

  USE task_kinds,ONLY: dp
  USE testrk_local
  USE librkfn
  IMPLICIT NONE
  EXTERNAL ODERK,ODERKN,RKF

  INTEGER:: method,ntmax
  REAL(dp):: abserr,relerr
  REAL(dp):: t0,f0(neqmax),t1,f1(neqmax),esterr(neqmax)
  REAL(dp):: t0x,f0x(neqmax)
  INTEGER:: init,nit,ierr
  REAL(dp):: htemp,y0temp(neqmax),y1temp(neqmax),errtemp(neqmax)
  REAL(dp):: work2(neqmax,2),work11(neqmax,11)

  model=1
  method=1
  ntmax=100
  abserr=1.D-6
  relerr=1.D-6
  t0=0.D0
  t1=1.D0
  f0(1)=1.D0
  parm1=1.D0
  parm2=0.D0
  
1 REWIND(5)
2 CONTINUE
  WRITE(6,'(A,I3,2ES12.4)') &
       '## INPUT model,parm1,parm2=', &
       model,parm1,parm2
  READ(5,*,ERR=2,END=9000) model,parm1,parm2

3 REWIND(5)
4 CONTINUE
  WRITE(6,'(A,I3,I6,2ES12.4)') &
       '## INPUT method,ntmax,abserr,relerr=', &
       method,ntmax,abserr,relerr
  READ(5,*,ERR=4,END=1) method,ntmax,abserr,relerr

5 CONTINUE
  WRITE(6,'(A,3ES12.4)') &
       '## INPUT t0,f0,t1=',t0,f0(1),t1
  READ(5,*,ERR=5,END=3) t0,f0(1),t1

  SELECT CASE(method)
  CASE(1)
     t0x=t0
     f0x(1:neqmax)=f0(1:neqmax)
     CALL ODERK(neqmax,SUB,t0x,t1,ntmax,f0x,f1,work2)
     esterr=0.D0
     nit=ntmax
     ierr=0
  CASE(2)
     init=1
     t0x=t0
     f0x(1:neqmax)=f0(1:neqmax)
     CALL RKFN(neqmax,SUB,t0x,t1,f0x,init,relerr,abserr, &
              f1,esterr,nit,ierr,htemp,y0temp,y1temp,errtemp,work11)
  CASE(3)
     init=1
     t0x=t0
     f0x(1:neqmax)=f0(1:neqmax)
     CALL RKF(neqmax,SUB,t0x,t1,f0x,init,relerr,abserr, &
              f1,esterr,nit,ierr,htemp,y0temp,y1temp,errtemp,work11)
  END SELECT
  WRITE(6,'(A,ES12.4,ES18.10,I6,ES12.4)') &
       't1,f1,nit,esterr=',t1,f1(1),nit,esterr(1)
  GO TO 5

9000 CONTINUE
  STOP
END PROGRAM testrk
