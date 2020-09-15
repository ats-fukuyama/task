subroutine wmqsolv(EX,EY,EZ,W,A,CD,RR,dt,dx,omega,TMN)

  !$ use omp_lib
  use bpsd
  implicit none
  integer(4),intent(in)   :: W
  complex(8),intent(inout):: EX(W,W),EY(W,W),EZ(W,W)
  complex(8),intent(in)   :: A(W,W,3,3),CD(W,W,3,3)
  real(8)   ,intent(in)   :: RR,dt,dx,omega,TMN

  integer(4) :: i,j
  real(8)    :: kz
  complex(8) :: EXNEXT(W,W),EYNEXT(W,W),EZNEXT(W,W),B(3),P,Q,R1,R2,R3

  kz = TMN/RR

  P  = omega**2/VC**2
  Q  = CI*omega*RMU0
  R1 = 0.25d0/dx**2
  R2 = 1.00d0/dx**2
  R3 = 0.50d0*CI*kz/dx

  do j=1,W
     do i=1,W
        EXNEXT(i,j)=0.0d0
        EYNEXT(i,j)=0.0d0
        EZNEXT(i,j)=0.0d0
     end do
  end do

!$omp do        
  do j=2,W-1
     do i=2,W-1
        
        B(1)= R1   *(EY(i+1,j+1)-     EY(i+1,j-1)-EY(i-1,j+1)+EY(i-1,j-1))&
             +R3   *(EZ(i+1,j  )-     EZ(i-1,j  )                        )&
             -R2   *(EX(i  ,j+1)-2.d0*EX(i  ,j  )+EX(i  ,j-1)            )&
             +kz**2* EX(i  ,j  )                                          &
             -P    * EX(i  ,j  )                                          &
             -Q    *(CD(i,j,1,1)*EX(i,j)+CD(i,j,1,2)*EY(i,j)+CD(i,j,1,3)*EZ(i,j))
        
        B(2)= R3   *(EZ(i  ,j+1)-     EZ(i  ,j-1)                        )&
             +R1   *(EX(i+1,j+1)-     EX(i+1,j-1)-EX(i-1,j+1)+EX(i-1,j-1))&
             +kz**2* EY(i  ,j  )                                          &
             -R2   *(EY(i+1,j  )-2.d0*EY(i  ,j  )+EY(i-1,j  )            )&
             -P    * EY(i  ,j  )                                          &
             -Q    *(CD(i,j,2,1)*EX(i,j)+CD(i,j,2,2)*EY(i,j)+CD(i,j,2,3)*EZ(i,j))

        B(3)= R3   *(EX(i+1,j  )-     EX(i-1,j  )                        )&
             +R3   *(EY(i  ,j+1)-     EY(i  ,j-1)                        )&
             -R2   *(EZ(i+1,j  )-2.d0*EZ(i  ,j  )+EZ(i-1,j  )            )&
             -R2   *(EZ(i  ,j+1)-2.d0*EZ(i  ,j  )+EZ(i  ,j-1)            )&
             -P    * EZ(i  ,j  )                                          &
             -Q    *(CD(i,j,3,1)*EX(i,j)+CD(i,j,3,2)*EY(i,j)+CD(i,j,3,3)*EZ(i,j))

        EXNEXT(i,j)=dt*(A(i,j,1,1)*B(1)+A(i,j,1,2)*B(2)+A(i,j,1,3)*B(3))+EX(i,j)
        EYNEXT(i,j)=dt*(A(i,j,2,1)*B(1)+A(i,j,2,2)*B(2)+A(i,j,2,3)*B(3))+EY(i,j)
        EZNEXT(i,j)=dt*(A(i,j,3,1)*B(1)+A(i,j,3,2)*B(2)+A(i,j,3,3)*B(3))+EZ(i,j)

     end do
  end do

  do j=1,W
     do i=1,W

        if(i.eq.1)then
           EX(i,j)=EXNEXT(i+1,j  )
        else if(i.eq.W)then
           EX(i,j)=EXNEXT(i-1,j  )
        else if(j.eq.1)then
           EY(i,j)=EYNEXT(i  ,j+1)
        else if(j.eq.W)then
           EY(i,j)=EYNEXT(i  ,j-1)
        else
           EX(i,j)=EXNEXT(i,j)
           EY(i,j)=EYNEXT(i,j)
           EZ(i,j)=EZNEXT(i,j)
        end if
        
     end do
  end do

  return
end subroutine wmqsolv
