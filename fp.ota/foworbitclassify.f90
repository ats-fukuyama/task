module foworbitclassify
  use fowcomm
  use fpcomm
  implicit none

  private

  public :: trapped, passing, forbitten

contains

  function trapped(np,nth,nr,nsa,mode)
    ! mode(1) = 0 :  np is mesh point
    !           1 :  np is grid point
    ! mode(2) = 0 : nth is mesh point
    !           1 : nth is grid point
    ! mode(3) = 0 : nr is mesh point
    !           1 : nr is grid point

    logical :: trapped
    integer,intent(in) :: np, nth, nr,nsa,mode(3)
    integer :: ir
    real(rkind) :: xil, Bml, psiml, D_orbit, Gm, bn, psin, dnr_bn, lorentzFactor, beta

    select case(mode(1))
    case(0)
      lorentzFactor = sqrt(1.d0+theta0(nsa)*pm(np,nsa)**2)
      beta = sqrt(1.d0-1.d0/lorentzFactor**2)
    case(1)
      lorentzFactor = sqrt(1.d0+theta0(nsa)*pg(np,nsa)**2)
      beta = sqrt(1.d0-1.d0/lorentzFactor**2)
    end select

    select case(mode(2))
    case(0)
      xil = xi(nth)
    case(1)
      xil = xig(nth)
    end select

    select case(mode(3))
    case(0)
      if ( xil>0.d0 ) then
        Bml = Bout(nr)
      else if (xil<0.d0) then
        Bml = Bin(nr)
      else
        trapped = .false.
        return
      end if
      psiml = psim(nr)
      Gm = aefp(nsa)*Bml*psiml/(amfp(nsa)*vc*Fpsi(nr))
      Bn = Bml/(1.d0-xil**2)

      if( Bn<=Boutg(1) ) then
        do ir = 1, nrmax
          if( Boutg(ir+1)<=Bn ) then
            dnr_bn = (Bn-Boutg(ir))/(Boutg(ir+1)-Boutg(ir))
            psin = psimg(ir)+(psimg(ir+1)-psimg(ir))*dnr_bn
            exit
          else if ( Bn<Boutg(nrmax+1) ) then
            dnr_bn = (Bn-Boutg(nrmax+1))/(Boutg(nrmax+1)-Boutg(nrmax))
            psin = psimg(nrmax+1)+(psimg(nrmax+1)-psimg(nrmax))*dnr_bn
            exit
          end if
        end do
      else
        do ir = 1, nrmax
          if( Bn<=Bing(ir+1) ) then
            dnr_bn = (Bn-Bing(ir))/(Bing(ir+1)-Bing(ir))
            psin = psimg(ir)+(psimg(ir+1)-psimg(ir))*dnr_bn
            exit
          else if ( Bn>Bing(nrmax+1) ) then
            dnr_bn = (Bn-Bing(nrmax+1))/(Bing(nrmax+1)-Bing(nrmax))
            psin = psimg(nrmax+1)+(psimg(nrmax+1)-psimg(nrmax))*dnr_bn
            exit
          end if
        end do
      end if

    case(1)
      if ( xil>0.d0 ) then
        Bml = Boutg(nr)
      else if ( xil<0.d0 ) then
        Bml = Bing(nr)
      else
        trapped = .false.
        return
      end if
      psiml = psimg(nr)
      Gm = aefp(nsa)*Bml*psiml/(amfp(nsa)*vc*Fpsig(nr))
      Bn = Bml/(1.d0-xil**2)

      if ( Bn<=Boutg(1) ) then
        do ir = 1,nrmax
          if ( Bout(ir+1)<=Bn )then
            dnr_bn = (Bn-Boutg(ir))/(Boutg(ir+1)-Boutg(ir))
            psin = psimg(ir)+(psimg(ir+1)-psimg(ir))*dnr_bn
            exit
          else if ( Bn<Boutg(nrmax+1) ) then
            dnr_bn = (Bn-Boutg(nrmax+1))/(Boutg(nrmax+1)-Boutg(nrmax))
            psin = psimg(nrmax+1)+(psimg(nrmax+1)-psimg(nrmax))*dnr_bn
            exit
          end if
        end do
      else
        do ir = 1, nrmax
          if ( Bn <= Bing(ir+1) ) then
            dnr_bn = (Bn-Bing(ir))/(Bing(ir+1)-Bing(ir))
            psin = psimg(ir)+(psimg(ir+1)-psimg(ir))*dnr_bn
            exit
          else if ( Bn>Bing(nrmax+1) ) then
            dnr_bn = (Bn-Bing(nrmax+1))/(Bing(nrmax+1)-Bing(nrmax))
            psin = psim(nrmax+1)+(psimg(nrmax)-psimg(nrmax))*dnr_bn
            exit
          end if
        end do
      end if

    end select

    if( Boutg(nrmax+1)<=Bn .and. Bn<=Bing(nrmax+1) .and. psin<=psiml ) then
      D_orbit = 1.d0/sqrt(1.d0+xil**2/(Gm*(1.d0-psin/psiml))**2)
    else
      D_orbit = 0.d0
    end if

    if( beta<=D_orbit ) then
      trapped = .true.
    else
      trapped = .false.
    end if

  end function  trapped

  function passing(np,nth,nr,nsa,mode)
    ! mode(1) = 0 :  np is mesh point
    !           1 :  np is grid point
    ! mode(2) = 0 : nth is mesh point
    !           1 : nth is grid point
    ! mode(3) = 0 : nr is mesh point
    !           1 : nr is grid point

    logical :: passing
    integer,intent(in) :: np,nth,nr,nsa,mode(3)

    passing = (.not. trapped(np,nth,nr,nsa,mode)).and.(.not.forbitten(np,nth,nr,nsa,mode))

  end function passing

  function forbitten(np,nth,nr,nsa,mode)
    ! mode(1) = 0 :  np is mesh point
    !           1 :  np is grid point
    ! mode(2) = 0 : nth is mesh point
    !           1 : nth is grid point
    ! mode(3) = 0 : nr is mesh point
    !           1 : nr is grid point
    
    logical :: forbitten
    integer,intent(in) :: np,nth,nr,nsa,mode(3)
    integer :: ir
    real(rkind) :: xil, Bml, psiml, Fpsil, Gm, lorentzFactor, beta&
                  , dFdpsi, dBmdpsi, stagnation_orbit

    select case(mode(1))
    case(0)
      lorentzFactor = sqrt(1.d0+theta0(nsa)*pm(np,nsa)**2)
      beta = sqrt(1.d0-1.d0/lorentzFactor**2)
    case(1)
      lorentzFactor = sqrt(1.d0+theta0(nsa)*pg(np,nsa)**2)
      beta = sqrt(1.d0-1.d0/lorentzFactor**2)
    end select


    select case(mode(2))
    case(0)
      xil = xi(nth)
    case(1)
      xil = xig(nth)
    end select

    select case(mode(3))
    case(0)
      if ( xil>=0.d0 ) then
        Bml = Bout(nr)
        if ( nr==1 ) then
          dBmdpsi = (-3*Bout(1)+4*Bout(2)-Bout(3))/(-3*psim(1)+4*psim(2)-psim(3))
          dFdpsi = (-3*Fpsi(1)+4*Fpsi(2)-Fpsi(3))/(-3*psim(1)+4*psim(2)-psim(3))
        else if ( nr==nrmax ) then
          dBmdpsi = (3*Bout(nrmax)-4*Bout(nrmax-1)+Bout(nrmax-2))&
                  /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
          dFdpsi = (3*Fpsi(nrmax)-4*Fpsi(nrmax-1)+Fpsi(nrmax-2))&
                  /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
        else
          dBmdpsi = (Bout(nr+1)-Bout(nr-1))/(psim(nr+1)-psim(nr-1))
          dFdpsi = (Fpsi(nr+1)-Fpsi(nr-1))/(psim(nr+1)-psim(nr-1))
        end if

      else if ( xil<0.d0) then
        Bml = Bin(nr)
        if ( nr==1 ) then
          dBmdpsi = (-3*Bin(1)+4*Bin(2)-Bin(3))/(-3*psim(1)+4*psim(2)-psim(3))
          dFdpsi = (-3*Fpsi(1)+4*Fpsi(2)-Fpsi(3))/(-3*psim(1)+4*psim(2)-psim(3))
        else if ( nr==nrmax ) then
          dBmdpsi = (3*Bin(nrmax)-4*Bin(nrmax-1)+Bin(nrmax-2))&
                  /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
          dFdpsi = (3*Fpsi(nrmax)-4*Fpsi(nrmax-1)+Fpsi(nrmax-2))&
                  /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
        else
          dBmdpsi = (Bin(nr+1)-Bin(nr-1))/(psim(nr+1)-psim(nr-1))
          dFdpsi = (Fpsi(nr+1)-Fpsi(nr-1))/(psim(nr+1)-psim(nr-1))
        end if

      end if
      psiml = psim(nr)
      Fpsil = Fpsi(nr)
      Gm = aefp(nsa)*Bml*psiml/(amfp(nsa)*vc*Fpsil)

    case(1)
      if ( xil>=0.d0 ) then
        Bml = Boutg(nr)
        if ( nr==1 ) then
          dBmdpsi = (-3*Boutg(1)+4*Boutg(2)-Boutg(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
          dFdpsi = (-3*Fpsig(1)+4*Fpsig(2)-Fpsig(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
        else if ( nr==nrmax+1 ) then
          dBmdpsi = (3*Boutg(nrmax+1)-4*Boutg(nrmax)+Boutg(nrmax-1))&
                 /(3*psimg(nrmax+1)-4*psimg(nrmax)+psimg(nrmax-1))
          dFdpsi = (3*Fpsig(nrmax+1)-4*Fpsig(nrmax)+Fpsig(nrmax-1))&
                /(3*psimg(nrmax+1)-4*psimg(nrmax)+psimg(nrmax-1))
        else
          dBmdpsi = (Boutg(nr+1)-Boutg(nr-1))/(psimg(nr+1)-psimg(nr-1))
          dFdpsi = (Fpsig(nr+1)-Fpsig(nr-1))/(psimg(nr+1)-psimg(nr-1))
        end if

      else
        Bml = Bing(nr)
        if ( nr==1 ) then
          dBmdpsi = (-3*Bing(1)+4*Bing(2)-Bing(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
          dFdpsi = (-3*Fpsig(1)+4*Fpsig(2)-Fpsig(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
        else if ( nr==nrmax+1 ) then
          dBmdpsi = (3*Bing(nrmax+1)-4*Bing(nrmax)+Bing(nrmax-1))&
                 /(3*psimg(nrmax+1)-4*psimg(nrmax)+psimg(nrmax-1))
          dFdpsi = (3*Fpsig(nrmax+1)-4*Fpsig(nrmax)+Fpsig(nrmax-1))&
                /(3*psimg(nrmax+1)-4*psimg(nrmax)+psimg(nrmax-1))
        else
          dBmdpsi = (Bing(nr+1)-Bing(nr-1))/(psimg(nr+1)-psimg(nr-1))
          dFdpsi = (Fpsig(nr+1)-Fpsig(nr-1))/(psimg(nr+1)-psimg(nr-1))
        end if
      end if
      psiml = psimg(nr)
      Fpsil = Fpsig(nr)
      Gm = aefp(nsa)*Bml*psiml/(amfp(nsa)*vc*Fpsil)

    end select

    stagnation_orbit = Gm*xil/(xil**2*dFdpsi/Fpsil-0.5d0*(1.d0+xil**2)*dBmdpsi/Bml)
    stagnation_orbit = (1.d0+stagnation_orbit**(-2))**(-0.5d0)

    if ( beta<stagnation_orbit ) then
      forbitten = .false.
    else
      forbitten = .true.
    end if

  end function forbitten

end module foworbitclassify
