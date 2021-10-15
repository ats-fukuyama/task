! wqsolv.f90

MODULE wqsolv

  PRIVATE
  PUBLIC wq_solv

CONTAINS
  
subroutine wq_solv

  use wqcomm
  implicit none

  integer(4) :: nx,ny
  real(8)    :: kz
  complex(8) :: EXNEXT(nxmax,nymax),EYNEXT(nxmax,nymax),EZNEXT(nxmax,nymax)
  complex(8) :: BV(3),P,Q,R1xy,R2xx,R2yy,R3x,R3y
  COMPLEX(rkind):: CM(3,3,-1:1,-1,1)

  kz = TMN/RR

  P  = omega**2/VC**2
  Q  = CI*omega*RMU0
  R1xy = 0.25d0/(dx*dy)
  R2xx = 1.00d0/dx**2
  R2yy = 1.00d0/dy**2
  R3x = 0.50d0*CI*kz/dx
  R3y = 0.50d0*CI*kz/dy

  do ny=1,nymax
     do nx=1,nxmax
        EXNEXT(nx,ny)=0.0d0
        EYNEXT(nx,ny)=0.0d0
        EZNEXT(nx,ny)=0.0d0
     end do
  end do

!$omp do        
  do ny=2,nymax-1
     do nx=2,nxmax-1
        
        BV(1)= R1xy *(EY(nx+1,ny+1)-     EY(nx+1,ny-1)  &
                     -EY(nx-1,ny+1)+     EY(nx-1,ny-1)) &
              +R3x  *(EZ(nx+1,ny  )-     EZ(nx-1,ny  )              ) &
              -R2yy *(EX(nx  ,ny+1)-2.d0*EX(nx  ,ny  )+EX(nx  ,ny-1)) &
              +kz**2* EX(nx  ,ny  ) &
              -P    * EX(nx  ,ny  ) &
              -Q    *(CD(1,1,nx,ny)*EX(nx,ny) &
                     +CD(1,2,nx,ny)*EY(nx,ny) &
                     +CD(1,3,nx,ny)*EZ(nx,ny))
        CM(1,2,+1,+1)= R1xy
        CM(1,2,+1,-1)=-R1xy
        CM(1,2,-1,+1)=-R1xy
        CM(1,2,-1,-1)= R1xy
        CM(1,3,+1, 0)= R3x
        CM(1,3,-1, 0)=-R3x
        CM(1,1, 0,+1)=-     R2yy
        CM(1,1, 0, 0)=+2.D0*R2yy+kz**2-P-Q*CD(1,1,nx,ny)
        CM(1,1, 0,+1)=-     R2yy
        CM(1,2, 0, 0)=                  -Q*CD(1,2,nx,ny)
        CM(1,3, 0, 0)=                  -Q*CD(1,3,nx,ny)
        
        BV(2)= R3y  *(EZ(nx  ,ny+1)-     EZ(nx  ,ny-1)) &
              +R1xy *(EX(nx+1,ny+1)-     EX(nx+1,ny-1)  &
                     -EX(nx-1,ny+1)+     EX(nx-1,ny-1)) &
              +kz**2* EY(nx  ,ny  ) &
              -R2xx *(EY(nx+1,ny  )-2.d0*EY(nx  ,ny  )+EY(nx-1,ny  )) &
              -P    * EY(nx  ,ny  ) &
              -Q    *(CD(2,1,nx,ny)*EX(nx,ny) &
                     +CD(2,2,nx,ny)*EY(nx,ny) &
                     +CD(2,3,nx,ny)*EZ(nx,ny))

        CM2(3,nx  ,ny+1)= R3y
        CM2(3,nx  ,ny-1)=-R3y
        CM2(1,nx+1,ny+1)= R1xy
        CM2(1,nx+1,ny-1)=-R1xy
        CM2(1,nx-1,ny+1)=-R1xy
        CM2(1,nx-1,ny-1)= R1xy
        CM2(2,nx+1,ny  )=-     R2xx
        CM2(2,nx  ,ny  )=+2.D0*R2xx+kz**2-P-Q*CD(2,2,nx,ny)
        CM2(2,nx-1  ,ny)=-     R2xx
        CM2(1,nx  ,ny  )=                  -Q*CD(2,1,nx,ny)
        CM2(3,nx  ,ny  )=                  -Q*CD(2,3,nx,ny)
        
        BV(3)= R3x  *(EX(nx+1,ny  )-     EX(nx-1,ny  )) &
              +R3y  *(EY(nx  ,ny+1)-     EY(nx  ,ny-1)) &
              -R2xx *(EZ(nx+1,ny  )-2.d0*EZ(nx  ,ny  )+EZ(nx-1,ny  )) &
              -R2yy *(EZ(nx  ,ny+1)-2.d0*EZ(nx  ,ny  )+EZ(nx  ,ny-1)) &
              -P    * EZ(nx  ,ny  ) &
              -Q    *(CD(3,1,nx,ny)*EX(nx,ny) &
                     +CD(3,2,nx,ny)*EY(nx,ny) &
                     +CD(3,3,nx,ny)*EZ(nx,ny))

        CM3(1,nx+1,ny  )= R3x
        CM3(1,nx-1,ny  )=-R3x
        CM3(2,nx  ,ny+1)= R3y
        CM3(2,nx  ,ny-1)=-R3y
        CM3(3,nx+1,ny  )=-     R2xx
        CM3(3,nx  ,ny+1)=-               R2yy
        CM3(3,nx  ,ny  )=+2.D0*R2xx+2.D0*R2yy-P-Q*CD(3,3,nx,ny)
        CM3(3,nx-1,ny  )=-     R2xx
        CM3(3,nx  ,ny-1)=-               R2yy
        CM3(1,nx  ,ny  )=                      -Q*CD(3,1,nx,ny)
        CM3(2,nx  ,ny  )=                      -Q*CD(3,2,nx,ny)

        EXNEXT(nx,ny)=dt*(Ainv(1,1,nx,ny)*BV(1) &
                         +Ainv(1,2,nx,ny)*BV(2) &
                         +Ainv(1,3,nx,ny)*BV(3))+EX(nx,ny)
        EYNEXT(nx,ny)=dt*(Ainv(2,1,nx,ny)*BV(1) &
                         +Ainv(2,2,nx,ny)*BV(2) &
                         +Ainv(2,3,nx,ny)*BV(3))+EY(nx,ny)
        EZNEXT(nx,ny)=dt*(Ainv(3,1,nx,ny)*BV(1) &
                         +Ainv(3,2,nx,ny)*BV(2) &
                         +Ainv(3,3,nx,ny)*BV(3))+EZ(nx,ny)

     end do
  end do

  do ny=1,nymax
     do nx=1,nxmax

        if(nx.eq.1)then
           EX(nx,ny)=EXNEXT(nx+1,ny  )
        else if(nx.eq.nxmax)then
           EX(nx,ny)=EXNEXT(nx-1,ny  )
        else if(ny.eq.1)then
           EY(nx,ny)=EYNEXT(nx  ,ny+1)
        else if(ny.eq.nymax)then
           EY(nx,ny)=EYNEXT(nx  ,ny-1)
        else
           EX(nx,ny)=EXNEXT(nx,ny)
           EY(nx,ny)=EYNEXT(nx,ny)
           EZ(nx,ny)=EZNEXT(nx,ny)
        end if
        
     end do
  end do

  return
end subroutine wq_solv
END MODULE wqsolv
