!!! Miscellaneous libraries related to the physics or physical model

!***************************************************************
!
!   Collision frequency between i and j species
!
!     Inputs (integer*4): NR : radial grid number
!            (integer*4):  i : species number
!            (integer*4):  j : species number
!            (real*8) optional : eps : inverse aspect ratio, output effictive coll. freq.
!     Output (real*8)   : coll_freq : collision frequency [1/s]
!
!***************************************************************

real(8) function coll_freq(NR,i,j,eps)
  use tx_commons, only : pi, aee, eps0, amp, rkev, achg, amas, rr, q, Var
  use tx_interface, only : coulog
  implicit none
  integer(4), intent(in) :: NR, i, j
  real(8), intent(in), optional :: eps
  real(8) :: const, sqeps3, vtm

  const = AEE**4 * 1.d20 / (eps0**2 * 6.d0 * pi * sqrt(2.d0 * pi * amp * rKeV) * rKeV)

  coll_freq = Var(NR,j)%n * (achg(i)*achg(j))**2 / ( sqrt(amas(i) * Var(NR,i)%T) * Var(NR,i)%T )  &
       &    * coulog(1.d0,Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amas(i),achg(i),amas(j),achg(j)) &
       &    * const

  ! Effective collision frequency
  if(present(eps)) then
     sqeps3 = sqrt(eps) * eps
     if( sqeps3 == 0.d0 ) then
        coll_freq = 0.d0 ! Note that it must be huge in nature.
     else
        vtm = sqrt(2.d0 * Var(NR,i)%T * rKeV / (amas(i) * amp))
        coll_freq = coll_freq * (rr * q(NR) / (sqeps3 * vtm))
     end if
  end if
    
end function coll_freq

!***************************************************************
!
!   Correction factor for resistivity
!     (Hirshman and Sigmar, (1981), Eq. (7.36))
!
!***************************************************************

pure REAL(8) FUNCTION CORR(X)
  ! X is the effective charge number
  real(8), intent(in) :: X

  CORR = (1.D0 + (1.198D0 + 0.222D0 * X) * X) * X &
  &    / (1.D0 + (2.966D0 + 0.753D0 * X) * X)

END FUNCTION CORR

!***************************************************************
!
!   Trapped particle fraction
!     (Y. B. Kim, et al., Phys. Fluids B 3 (1990) 2050)
!
!***************************************************************

elemental real(8) function ftfunc(x)
  ! x is the inverse aspect ratio ; x is either scaler or array.
  real(8), intent(in) :: x

  ftfunc = (1.46d0 - 0.46d0 * x) * sqrt(x)

end function ftfunc
