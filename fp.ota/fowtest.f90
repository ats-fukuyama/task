subroutine fow_forbitten_boundary(upper_boundary)
  use fowcomm
  use fpcomm

  real(rkind),intent(out) :: upper_boundary(:,:,:)
  integer :: np,nth,nr,nsa
  integer :: ir
  real(rkind) :: v_stagnation_orbit
  real(rkind),allocatable :: B_m(:,:), Gm(:,:,:), dBmdpsi(:,:), dFdpsi(:,:)


  lorentzFactor = sqrt(1.d0+theta0(nsa)*pm(np,nsa)**2)
  beta = sqrt(1.d0-1.d0/lorentzFactor**2)

  
  allocate(B_m(nthmax,nrmax),Gm(nthmax,nrmax,nsamax),dBmdpsi(nthmax,nrmax), dFdpsi(nthmax,nrmax))

  do nr = 1, nrmax
    do nth = 1, nthmax
      if ( xi(nth) <= 0.d0 ) then
        B_m(nth,nr) = Bin(nr)
      else
        B_m(nth,nr) = Bout(nr)
      end if
    end do
  end do

  do nsa = 1, nsamax
    do nr = 1, nrmax
      do nth = 1, nthmax
        Gm(nth,nr,nsa) = aefp(nsa)*B_m(nth,nr)*psim(nr)/(amfp(nsa)*vc*Fpsi(nr))
      end do
    end do
  end do

  do nr = 1, nrmax
    do nth = 1, nthmax
      if ( xi(nth)>=0.d0 ) then
        if ( nr==1 ) then
          dBmdpsi(nth,nr) = (-3*Bout(1)+4*Bout(2)-Bout(3))/(-3*psim(1)+4*psim(2)-psim(3))
          dFdpsi(nth,nr) = (-3*Fpsi(1)+4*Fpsi(2)-Fpsi(3))/(-3*psim(1)+4*psim(2)-psim(3))
        else if ( nr==nrmax ) then
          dBmdpsi(nth,nr) = (3*Bout(nrmax)-4*Bout(nrmax-1)+Bout(nrmax-2))&
                  /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
          dFdpsi(nth,nr) = (3*Fpsi(nrmax)-4*Fpsi(nrmax-1)+Fpsi(nrmax-2))&
                  /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
        else
          dBmdpsi(nth,nr) = (Bout(nr+1)-Bout(nr-1))/(psim(nr+1)-psim(nr-1))
          dFdpsi(nth,nr) = (Fpsi(nr+1)-Fpsi(nr-1))/(psim(nr+1)-psim(nr-1))
        end if
    
      else if ( xi(nth)<0.d0) then
        if ( nr==1 ) then
          dBmdpsi(nth,nr) = (-3*Bin(1)+4*Bin(2)-Bin(3))/(-3*psim(1)+4*psim(2)-psim(3))
          dFdpsi(nth,nr) = (-3*Fpsi(1)+4*Fpsi(2)-Fpsi(3))/(-3*psim(1)+4*psim(2)-psim(3))
        else if ( nr==nrmax ) then
          dBmdpsi(nth,nr) = (3*Bin(nrmax)-4*Bin(nrmax-1)+Bin(nrmax-2))&
                  /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
          dFdpsi(nth,nr) = (3*Fpsi(nrmax)-4*Fpsi(nrmax-1)+Fpsi(nrmax-2))&
                  /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
        else
          dBmdpsi(nth,nr) = (Bin(nr+1)-Bin(nr-1))/(psim(nr+1)-psim(nr-1))
          dFdpsi(nth,nr) = (Fpsi(nr+1)-Fpsi(nr-1))/(psim(nr+1)-psim(nr-1))
        end if
      end if    
    end do
  end do

  do nsa = 1, nsamax
    do nr = 1, nrmax
      do nth = 1, nthmax
        v_stagnation_orbit = Gm(nth,nr,nsa)*xi(nth)/(xi(nth)**2*dFdpsi(nth,nr,nsa)/Fpsi(nr)&
                                      -0.5d0*(1.d0+xi(nth)**2)*dBmdpsi(nth,nr)/Bm(nth,nr)) ! LHS = gamma*beta 
        v_stagnation_orbit = vc*sqrt(v_stagnation_orbit**2/(1.d0+v_stagnation_orbit**2)) ! LHS = velocity of stagnation orbit

        upper_boundary(nth,nr,nsa) = amfp(nsa)*v_stagnation_orbit/(1.d0-v_stagnation_orbit**2/vc**2) ! mmomentum of stagnation orbit
        upper_boundary(nth,nr,nsa) = upper_boundary(nth,nr,nsa)/ptfp0(nsa) ! normalize
      end do
    end do
  end do

end subroutine fow_forbitten_boundary