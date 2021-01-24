! libderiv.pf90  

MODULE libderiv

  PRIVATE
  PUBLIC tri_deriv,duo_deriv

CONTAINS

! *** calculate function and derivatives, f0, df0, ddf0, at x0
! *** from three points data f(1),f(2),f(3) at x(1),x(2),x(3)
  
  SUBROUTINE tri_deriv(x,f,x0,f0,df0,ddf0,ierr)
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: x(3),f(3),x0
    REAL(dp),INTENT(OUT):: f0,df0,ddf0
    INTEGER,INTENT(OUT):: ierr
    REAL(dp):: dx1,dx2,dx3,addf0,adf0

    ierr=0
    dx1=x(1)-x0
    dx2=x(2)-x0
    dx3=x(3)-x0

    addf0=0.5D0*((dx2**2-dx1**2)*(dx3-dx2)-(dx3**2-dx2**2)*(dx2-dx1))
    IF(ABS(addf0).LE.0.D0) THEN
       WRITE(6,'(A)') 'XX tri_deriv: addf0=0.D0 : INCONSISTENT POSITION x'
       WRITE(6,'(A,3ES12.4)') '   x(1),x(2),x(3)=',x(1),x(2),x(3)
       ierr=1
       RETURN
    END IF
    ddf0=((f(2)-f(1))*(dx3-dx2)-(f(3)-f(2))*(dx2-dx1))/addf0
    adf0=dx2-dx1
    df0=(f(2)-f(1)-0.5D0*ddf0*(dx2**2-dx1**2))/adf0
    f0=f(1)-df0*dx1-0.5D0*ddf0*dx1**2
    RETURN
  END SUBROUTINE tri_deriv

! *** calculate function and derivatives, f0, df0, ddf0, at x0
! *** either of f0, df0, ddf0 is zero
! *** from three points data f(1),f(2) at x(1),x(2)
  
  SUBROUTINE duo_deriv(x,f,x0,f0,df0,ddf0,id,ierr)
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id     ! id=0:f0=0, id=1:df0=0, id=2:ddf0=0
    REAL(dp),INTENT(IN):: x(2),f(2),x0
    REAL(dp),INTENT(OUT):: f0,df0,ddf0
    INTEGER,INTENT(OUT):: ierr
    REAL(dp):: dx1,dx2,det

    ierr=0
    dx1=x(1)-x0
    dx2=x(2)-x0

    det=0.5D0*dx1*dx2*(dx2-dx1)
    IF(ABS(det).LE.0.D0) THEN
       WRITE(6,'(A)') 'XX duo_deriv: det=0.D0 : INCONSISTENT POSITION x'
       WRITE(6,'(A,2ES12.4)') '   x(1),x(2)=',x(1),x(2)
       ierr=1
       RETURN
    END IF
    
    SELECT CASE(id)
    CASE(0)
       f0=0.D0
       df0=0.5D0*(f(1)*dx2**2-f(2)*dx1**2)/det
       ddf0=     (-f(1)*dx2+f(2)*dx1)/det
    CASE(1)
       ddf0=2.D0*(f(2)-f(1))/(dx2**2-dx1**2)
       df0=0.D0
       f0=f(1)-0.5D0*ddf0*dx1**2
    CASE(2)
       df0=(f(2)-f(1))/(dx2-dx1)
       ddf0=0.D0
       f0=f(1)-df0*dx1
    CASE DEFAULT
       WRITE(6,'(A,I3)') 'XX duo_deriv: undeifned id: id=',id
       ierr=2
       RETURN
    END SELECT
    RETURN
  END SUBROUTINE duo_deriv
END MODULE libderiv
  
