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
  COMPLEX(rkind):: Comega,COEF,sigma,omega_coll
  INTEGER:: nx,ny,medium
  REAL(rkind):: OCEX(nxmax,nymax),OCEY(nxmax,nymax),OCEZ(nxmax,nymax)
  REAL(rkind):: BX(nxmax,nymax),BY(nxmax,nymax),BZ(nxmax,nymax)
  REAL(rkind):: x,y,r,q
  
  CD_(1:3,1:3,1:nxmax,1:nymax) = (0.D0,0.D0) ! vacuum

  DO ny=1,nymax
     DO nx=1,nxmax
        medium=medium_nx_ny(nx,ny)
        IF(medium.EQ.0) CYCLE
        SELECT CASE(id_medium(medium))
        CASE(1) ! plasma
           x  = xg(nx)*dxfactor
           y  = yg(ny)*dyfactor
           r   =  SQRT(x**2+y**2)

           IF(r.le.RA)THEN
              q=q0+(qa-q0)*(r/RA)**2
              ne(nx,ny)= density_medium(medium)*(1.d0-(r/RA)**2)
           ELSE
              q=qa*RA/r
              ne(nx,ny)= 0.0d0
           END IF
           BX(nx,ny)=-B0*y/(RR*q)
           BY(nx,ny)= B0*x/(RR*q)
           BZ(nx,ny)= B0*RR/(RR+x)

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
           
           Comega=-CI*omega_*(1.D0+CI*collision_medium(medium))

           COEF = ne(nx,ny)*AME**2/(AME*Comega*(OCE(nx,ny)**2+Comega**2))
           CD_(1,1,nx,ny) = COEF*(    Comega**2    +OCEX(nx,ny)**2       ) 
           CD_(1,2,nx,ny) = COEF*( OCEZ(nx,ny)*Comega+OCEX(nx,ny)*OCEY(nx,ny))
           CD_(1,3,nx,ny) = COEF*(-OCEY(nx,ny)*Comega+OCEZ(nx,ny)*OCEX(nx,ny))
           CD_(2,1,nx,ny) = COEF*(-OCEZ(nx,ny)*Comega+OCEX(nx,ny)*OCEY(nx,ny))
           CD_(2,2,nx,ny) = COEF*(    Comega**2    +OCEY(nx,ny)**2       )
           CD_(2,3,nx,ny) = COEF*( OCEX(nx,ny)*Comega+OCEY(nx,ny)*OCEZ(nx,ny))
           CD_(3,1,nx,ny) = COEF*( OCEY(nx,ny)*Comega+OCEZ(nx,ny)*OCEX(nx,ny))
           CD_(3,2,nx,ny) = COEF*(-OCEX(nx,ny)*Comega+OCEY(nx,ny)*OCEZ(nx,ny))
           CD_(3,3,nx,ny) = COEF*(    Comega**2    +OCEZ(nx,ny)**2       )
        CASE(2) ! dielectric
           sigma=omega_/(CI*VC**2*RMU0) &
                *(dielectric_medium(medium)-1.D0) &
                *(1.D0+CI*collision_medium(medium))
           CD(1,1,nx,ny)=sigma
           CD(2,2,nx,ny)=sigma
           CD(3,3,nx,ny)=sigma
        CASE(3) ! dielectric with frequency resonance
           sigma=omega_/(CI*VC**2*RMU0) &
                *(dielectric_medium(medium)-1.D0)
           omega_coll=omega*res_freq_medium(medium) &
                *(1.D0-CI*res_coll_medium(medium))
           sigma=sigma*omega_/(omega_-omega_coll)
           CD(1,1,nx,ny)=sigma
           CD(2,2,nx,ny)=sigma
           CD(3,3,nx,ny)=sigma
        END SELECT

    END DO
  END DO

  RETURN
END SUBROUTINE wq_tens

END MODULE wqtens
