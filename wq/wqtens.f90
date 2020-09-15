! wqtens

MODULE wqtens

  PRIVATE
  PUBLIC wq_tens

CONTAINS
  
subroutine wq_tens(omega_,CD_)

  use wqcomm
  implicit none

  REAL(rkind),INTENT(IN):: omega_
  COMPLEX(rkind),INTENT(OUT):: CD_(3,3,nxmax,nymax)
  COMPLEX(8) :: Comega,COEF
  integer(4) :: nx,ny
  REAL(8)    :: OCEX(nxmax,nymax),OCEY(nxmax,nymax),OCEZ(nxmax,nymax)
  REAL(8)    :: BX(nxmax,nymax),BY(nxmax,nymax),BZ(nxmax,nymax),r,rx,ry,q
  
  DO ny=1,nymax
     DO nx=1,nxmax
       
        rx  = (DBLE(nx)-DBLE(nxmax/2))*dx
        ry  = (DBLE(ny)-DBLE(nymax/2))*dy
        r   =  SQRT(rx**2+ry**2)

        IF(r.le.RA)THEN

           q=q0+(qa-q0)*(r/RA)**2

           BX(nx,ny)=-B0*ry/(RR*q)
           BY(nx,ny)= B0*rx/(RR*q)
           BZ(nx,ny)= B0*RR/(RR+rx)
           ne(nx,ny)= n0*(1.d0-(r/RA)**2)

        ELSE

           q=qa*RA/r

           BX(nx,ny)=-B0*ry/(RR*q)
           BY(nx,ny)= B0*rx/(RR*q)
           BZ(nx,ny)= B0*RR/(RR+rx)
           ne(nx,ny)= 0.0d0

        END IF

     END DO
  END DO

!  DO ny=1,nymax
!    DO nx=1,nxmax
!        BX(nx,ny) = 0.d0
!        BY(nx,ny) = 0.d0
!        BZ(nx,ny) = B0
!        ne(nx,ny) = n0*(DBLE(nxmax-nx)/DBLE(nxmax) &
!                       +DBLE(nymax-ny)/DBLE(nymax)
!     END DO
!  END DO

  DO ny=1,nymax
     DO nx=1,nxmax

        OPE  (nx,ny) = sqrt(ne(nx,ny)*AEE**2/(AME*EPS0))
        OCEX (nx,ny) = AEE*BX(nx,ny)/AME
        OCEY (nx,ny) = AEE*BY(nx,ny)/AME
        OCEZ (nx,ny) = AEE*BZ(nx,ny)/AME
        OCE  (nx,ny) = sqrt(OCEX(nx,ny)**2+OCEY(nx,ny)**2+OCEZ(nx,ny)**2)
        OUH  (nx,ny) = sqrt(OCE(nx,ny)**2+OPE(nx,ny)**2)
        OR   (nx,ny) = 0.5d0*(sqrt(OCE(nx,ny)**2+4.d0*OPE(nx,ny)**2) &
                     +ABS(OCE(nx,ny)))
        OL   (nx,ny) = 0.5d0*(sqrt(OCE(nx,ny)**2+4.d0*OPE(nx,ny)**2) &
                     -ABS(OCE(nx,ny)))

     END DO
  END DO

  Comega=-CI*omega_+nu !collision version of omega

  DO ny=1,nymax
     DO nx=1,nxmax
        
        COEF        = ne(nx,ny)*AME**2/(AME*Comega*(OCE(nx,ny)**2+Comega**2))

        CD_(1,1,nx,ny) = COEF*(    Comega**2    +OCEX(nx,ny)**2       ) 
        CD_(1,2,nx,ny) = COEF*( OCEZ(nx,ny)*Comega+OCEX(nx,ny)*OCEY(nx,ny))
        CD_(1,3,nx,ny) = COEF*(-OCEY(nx,ny)*Comega+OCEZ(nx,ny)*OCEX(nx,ny))
        CD_(2,1,nx,ny) = COEF*(-OCEZ(nx,ny)*Comega+OCEX(nx,ny)*OCEY(nx,ny))
        CD_(2,2,nx,ny) = COEF*(    Comega**2    +OCEY(nx,ny)**2       )
        CD_(2,3,nx,ny) = COEF*( OCEX(nx,ny)*Comega+OCEY(nx,ny)*OCEZ(nx,ny))
        CD_(3,1,nx,ny) = COEF*( OCEY(nx,ny)*Comega+OCEZ(nx,ny)*OCEX(nx,ny))
        CD_(3,2,nx,ny) = COEF*(-OCEX(nx,ny)*Comega+OCEY(nx,ny)*OCEZ(nx,ny))
        CD_(3,3,nx,ny) = COEF*(    Comega**2    +OCEZ(nx,ny)**2       )

    END DO
  END DO

  RETURN
END SUBROUTINE wq_tens

END MODULE wqtens
