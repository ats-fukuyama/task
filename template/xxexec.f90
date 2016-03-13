!    $Id$

Module xxexec
  PRIVATE
  PUBLIC xx_exec
 
CONTAINS

  SUBROUTINE xx_exec(ierr)

    USE xxcomm,ONLY: nxmax,x0,dx,am,x,y,xx_allocate
    
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nx

    ierr=0
    DO nx=1,nxmax
       x(nx)=x0+dx*(nx-1)
       y(nx)=am*sin(x(nx))
    END DO

  END SUBROUTINE xx_exec
END Module xxexec
