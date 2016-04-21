!  ***** TASK/XXC PREPARATION *****

Module xxprep
  PRIVATE
  PUBLIC xx_prep
 
CONTAINS

  SUBROUTINE xx_prep(ierr)

    USE xxcomm,ONLY: nxmax,x0,dx,am,x,y,xx_allocate
    
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nx

    ierr=0
    DO nx=1,nxmax
       x(nx)=x0+dx*(nx-1)
       y(nx)=am*sin(x(nx))
    END DO

  END SUBROUTINE xx_prep
END Module xxprep
