! testrk.f90

MODULE testrk_local

  USE task_kinds,ONLY: dp
  INTEGER,PARAMETER:: neqmax=2
  INTEGER:: model
  REAL(dp):: parm1,parm2

CONTAINS

  SUBROUTINE SUB(X,Y,F)
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: X,Y(neqmax)
    REAL(dp),INTENT(OUT):: F(neqmax)

    SELECT CASE(MODEL)
    CASE(1)
       F(1)= Y(1)
       F(2)=-Y(2)
    CASE(2)
       F(1)=-Y(2)
       F(2)= Y(1)
    END SELECT
    RETURN
  END SUBROUTINE SUB

END MODULE testrk_local

PROGRAM testrk

  USE task_kinds,ONLY: dp
  USE testrk_local
  USE librkf
  IMPLICIT NONE
  EXTERNAL ODERK,ODERKN

  INTEGER:: method,ntmax
  REAL(dp):: abserr,relerr
  REAL(dp):: t0,f0(neqmax),t1,f1(neqmax),esterr(neqmax)
  REAL(dp):: t0x,f0x(neqmax),f1x,f2x
  INTEGER:: init,nit,ierr
  REAL(dp):: htemp,y0temp(neqmax),y1temp(neqmax),errtemp(neqmax)
  REAL(dp):: work2(neqmax,2),work11(neqmax,11)

  model=1
  method=1
  ntmax=1
  abserr=1.D-6
  relerr=1.D-6
  t0=0.D0
  t1=1.D0
  f0(1)=1.D0
  f0(2)=1.D0
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
  WRITE(6,'(A,4ES12.4)') &
       '## INPUT t0,f1,f2,t1=',t0,f0(1),f0(2),t1
  READ(5,*,ERR=5,END=3) t0,f0(1),f0(2),t1

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
     CALL RKF(neqmax,SUB,t0x,t1,f0x,init,relerr,abserr, &
              f1,esterr,nit,ierr,htemp,y0temp,y1temp,errtemp,work11)
  END SELECT
  SELECT CASE(model)
  CASE(1)
     f1x=EXP( t1)
     f2x=EXP(-t1)
  CASE(2)
     f1x=COS(t1)
     f2x=SIN(t1)
  END SELECT
  WRITE(6,'(A,ES12.4,ES18.10,I6,2ES12.4)') &
       't1,f1,nit,esterr=',t1,f1(1),nit,esterr(1),f1(1)-f1x
  WRITE(6,'(17X,ES12.4,ES18.10,I6,2ES12.4)') &
                           t1,f1(2),nit,esterr(2),f1(2)-f2x
  GO TO 5

9000 CONTINUE
  STOP
END PROGRAM testrk
