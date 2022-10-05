  ! wimeasure.f90

Module wimeasure

  PRIVATE
  PUBLIC wi_measure_pabs

CONTAINS

  ! *** measure max position and half width ***

  SUBROUTINE wi_measure_pabs

    USE wicomm
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nxmax+1):: x_nx,s_nx
    REAL(rkind):: x_max,s_max,d,x0,x1,x2,s0,s1,s2,xp,xm
    INTEGER:: nx,nxpos
    
    DO nx=1,nxmax+1 
       x_nx(nx)=xgrid(nx-1)
       s_nx(nx)=REAL(CPOWER(NX-1))
    END DO

    nxpos=1
    s_max=s_nx(nxpos)
    DO nx=2,nxmax+1
       IF(s_nx(nx).GT.s_max) THEN
          nxpos=nx
          s_max=s_nx(nxpos)
       END IF
    END DO
    IF(nxpos.NE.1.AND.nxpos.NE.nxmax+1) THEN
       ! s=s_max-d*(x-x_max)^2
       ! s0=s_max-d*(x0-x_max)^2
       ! s1=s_max-d*(x1-x_max)^2
       ! s2=s_max-d*(x2-x_max)^2
       ! s0-s1=-d*(x0-x_max)^2+d*(x1-x_max)^2=d*(x1^2-x0^2)+2*d*x_max*(x0-x1)
       ! s2-s1=-d*(x2-x_max)^2+d*(x1-x_max)^2=d*(x1^2-x2^2)+2*d*x_max*(x2-x1)
       ! (s0-s1)*(x2-x1)=d*(x1^2-x0^2)*(x2-x1)+2*d*x_max*(x0-x1)*(x2-x1)
       ! (s2-s1)*(x0-x1)=d*(x1^2-x2^2)*(x0-x1)+2*d*x_max*(x2-x1)*(x0-x1)
       ! d=((s0-s1)*(x2-x1)-(s2-s1)*(x0-x1)) &
       !   /((x1^2-x0^2)*(x2-x1)-(x1^2-x2^2)*(x0-x1)
       ! x_max=(s0-s1-d*(x1^2-x0^2))/(2*d*(x0-x1))
       ! s_max=s1+d*(x1-x_max)^2
       x0=x_nx(nxpos-1)
       x1=x_nx(nxpos)
       x2=x_nx(nxpos+1)
       s0=s_nx(nxpos-1)
       s1=s_nx(nxpos)
       s2=s_nx(nxpos+1)
       d=((s0-s1)*(x2-x1)-(s2-s1)*(x0-x1)) &
            /((x1**2-x0**2)*(x2-x1)-(x1**2-x2**2)*(x0-x1))
       x_max=(s0-s1-d*(x1**2-x0**2))/(2*d*(x0-x1))
       s_max=s1+d*(x1-x_max)**2
    ELSE
       x_max=x_nx(nxpos)
       s_max=s_nx(nxpos)
    END IF

    WRITE(6,'(A,3ES12.4)') &
         '## xpos:', x0,x1,x2
    WRITE(6,'(A,3ES12.4)') &
         '## spos:', s0,s1,s2
    WRITE(6,'(A,3ES12.4)') &
         '## x_max, s_max, 2nd deriv:', x_max,s_max,d

    ! position xp and xm where s=0.5D0*s_max
    !   x1   xp    x2
    !   s1 s_max/2 s2
    !   (xp-x1)/(x2-x1)=(s_max/2-s1)/(s2-s1)
    !   xp=x1+(x2-x1)*(s_max/2-s1)/(s2-s1)

    s0=0.5D0*s_max
    xp=x_nx(nxpos)
    DO nx=nxpos,nxmax+1
       IF(s_nx(nx).LT.s0) THEN
          x1=x_nx(nx-1)
          x2=x_nx(nx)
          s1=s_nx(nx-1)
          s2=s_nx(nx)
          xp=x1+(x2-x1)*(s0-s1)/(s2-s1)
          EXIT
       END IF
    END DO
    DO nx=nxpos,1,-1
       IF(s_nx(nx).LT.s0) THEN
          x1=x_nx(nx-1)
          x2=x_nx(nx)
          s1=s_nx(nx-1)
          s2=s_nx(nx)
          xm=x1+(x2-x1)*(s0-s1)/(s2-s1)
          EXIT
       END IF
    END DO
    WRITE(6,'(A,3ES12.4)') &
         '## width,xp,xm :', xp-xm,xp,xm
    
          RETURN
  END SUBROUTINE wi_measure_pabs
END Module wimeasure
