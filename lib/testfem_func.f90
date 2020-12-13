PROGRAM testfem_func

  USE task_kinds,ONLY: dp
  USE libgrf
  USE libfem

  IMPLICIT NONE
  INTEGER:: type

  CALL gsopen
  type=0

1 CONTINUE
  WRITE(6,*) '## input type 0:l 1:p 2:h 3:g '
  READ(5,*,ERR=1,END=9000) type

  SELECT CASE(type)
  CASE(0)
     CALL plot_fem_func(fem_func_l,2)
  CASE(1)
     CALL plot_fem_func(fem_func_p,2)
  CASE(2)
     CALL plot_fem_func(fem_func_h,4)
  CASE(3)
     CALL plot_fem_func(fem_func_g,4)
  END SELECT
  GO TO 1

9000 CONTINUE
  CALL gsclos
  STOP

CONTAINS

  SUBROUTINE plot_fem_func(func,idmax)
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    INTEGER,INTENT(IN):: idmax
    INTEGER:: id,nx
    INTEGER,PARAMETER:: nxmax=101
    REAL(dp):: x,dx
    REAL(dp),DIMENSION(nxmax):: gx
    REAL(dp),DIMENSION(nxmax,2*idmax):: gy1,gy2,gy3
    INTERFACE 
       FUNCTION func(x,id,mode)
         USE task_kinds,ONLY: dp
         REAL(dp),INTENT(IN):: x
         INTEGER,INTENT(IN):: id,mode
         REAL(dp):: func
       END FUNCTION func
    END INTERFACE

    dx=1.D0/(nxmax-1)
    DO nx=1,nxmax
       x=dx*(nx-1)
       gx(nx)=x
       DO id=1,idmax
          gy1(nx,id)=func(x,id,0)
          gy1(nx,id+idmax)=(func(x+dx,id,2)-func(x-dx,id,2))/(2.D0*dx)
          gy2(nx,id)=func(x,id,1)
          gy2(nx,id+idmax)=(func(x+dx,id,0)-func(x-dx,id,0))/(2.D0*dx)
          gy3(nx,id)=func(x,id,2)
       END DO
    END DO

    CALL PAGES
    CALL GRD1D(1,gx,gy1,nxmax,nxmax,2*idmax,'@Function@')
    CALL GRD1D(2,gx,gy2,nxmax,nxmax,2*idmax,'@Derivative@')
    CALL GRD1D(3,gx,gy3,nxmax,nxmax,  idmax,'@Integral@')
    CALL PAGEE
    RETURN
  END SUBROUTINE plot_fem_func
END PROGRAM testfem_func
