module foworbitclassify
  use fowcomm
  use fpcomm
  implicit none

  private

  public :: trapped, passing, forbitten

contains

  function trapped(np,nze,nps,nsa,mode)
    ! mode(1) = 0 :  np is mesh point
    !           1 :  np is grid point
    ! mode(2) = 0 : nze is mesh point
    !           1 : nze is grid point
    ! mode(3) = 0 : nps is mesh point
    !           1 : nps is grid point

    logical :: trapped
    integer,intent(in) :: np, nze, nps,nsa,mode(3)
    integer :: ips
    real(rkind) :: zetal, Bml, psiml, D_orbit, Gm, bn, psin, dnps_bn, lorentzFactor, beta

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
      zetal = zeta(nze)
    case(1)
      zetal = zetag(nze)
    end select

    select case(mode(3))
    case(0)
      if ( zetal>0.d0 ) then
        Bml = Bout(nps)
      else if (zetal<0.d0) then
        Bml = Bin(nps)
      else
        trapped = .false.
        return
      end if
      psiml = psim(nps)
      Gm = aefp(nsa)*Bml*psiml/(amfp(nsa)*vc*Fpsi(nps))
      Bn = Bml/(1.d0-zetal**2)

      if( Bn<=Boutg(1) ) then
        do ips = 1, npsmax
          if( Boutg(ips+1)<=Bn ) then
            dnps_bn = (Bn-Boutg(ips))/(Boutg(ips+1)-Boutg(ips))
            psin = psimg(ips)+(psimg(ips+1)-psimg(ips))*dnps_bn
            exit
          else if ( Bn<Boutg(npsmax+1) ) then
            dnps_bn = (Bn-Boutg(npsmax+1))/(Boutg(npsmax+1)-Boutg(npsmax))
            psin = psimg(npsmax+1)+(psimg(npsmax+1)-psimg(npsmax))*dnps_bn
            exit
          end if
        end do
      else
        do ips = 1, npsmax
          if( Bn<=Bing(ips+1) ) then
            dnps_bn = (Bn-Bing(ips))/(Bing(ips+1)-Bing(ips))
            psin = psimg(ips)+(psimg(ips+1)-psimg(ips))*dnps_bn
            exit
          else if ( Bn>Bing(npsmax+1) ) then
            dnps_bn = (Bn-Bing(npsmax+1))/(Bing(npsmax+1)-Bing(npsmax))
            psin = psimg(npsmax+1)+(psimg(npsmax+1)-psimg(npsmax))*dnps_bn
            exit
          end if
        end do
      end if

    case(1)
      if ( zetal>0.d0 ) then
        Bml = Boutg(nps)
      else if ( zetal<0.d0 ) then
        Bml = Bing(nps)
      else
        trapped = .false.
        return
      end if
      psiml = psimg(nps)
      Gm = aefp(nsa)*Bml*psiml/(amfp(nsa)*vc*Fpsig(nps))
      Bn = Bml/(1.d0-zetal**2)

      if ( Bn<=Boutg(1) ) then
        do ips = 1,npsmax
          if ( Bout(ips+1)<=Bn )then
            dnps_bn = (Bn-Boutg(ips))/(Boutg(ips+1)-Boutg(ips))
            psin = psimg(ips)+(psimg(ips+1)-psimg(ips))*dnps_bn
            exit
          else if ( Bn<Boutg(npsmax+1) ) then
            dnps_bn = (Bn-Boutg(npsmax+1))/(Boutg(npsmax+1)-Boutg(npsmax))
            psin = psimg(npsmax+1)+(psimg(npsmax+1)-psimg(npsmax))*dnps_bn
            exit
          end if
        end do
      else
        do ips = 1, npsmax
          if ( Bn <= Bing(ips+1) ) then
            dnps_bn = (Bn-Bing(ips))/(Bing(ips+1)-Bing(ips))
            psin = psimg(ips)+(psimg(ips+1)-psimg(ips))*dnps_bn
            exit
          else if ( Bn>Bing(npsmax+1) ) then
            dnps_bn = (Bn-Bing(npsmax+1))/(Bing(npsmax+1)-Bing(npsmax))
            psin = psim(npsmax+1)+(psimg(npsmax)-psimg(npsmax))*dnps_bn
            exit
          end if
        end do
      end if

    end select

    if( Boutg(npsmax+1)<=Bn .and. Bn<=Bing(npsmax+1) .and. psin<=psiml ) then
      D_orbit = 1.d0/sqrt(1.d0+zetal**2/(Gm*(1.d0-psin/psiml))**2)
    else
      D_orbit = 0.d0
    end if

    if( beta<=D_orbit ) then
      trapped = .true.
    else
      trapped = .false.
    end if

  end function  trapped

  function passing(np,nze,nps,nsa,mode)
    ! mode(1) = 0 :  np is mesh point
    !           1 :  np is grid point
    ! mode(2) = 0 : nze is mesh point
    !           1 : nze is grid point
    ! mode(3) = 0 : nps is mesh point
    !           1 : nps is grid point

    logical :: passing
    integer,intent(in) :: np,nze,nps,nsa,mode(3)

    passing = (.not. trapped(np,nze,nps,nsa,mode)).and.(.not.forbitten(np,nze,nps,nsa,mode))

  end function passing

  function forbitten(np,nze,nps,nsa,mode)
    ! mode(1) = 0 :  np is mesh point
    !           1 :  np is grid point
    ! mode(2) = 0 : nze is mesh point
    !           1 : nze is grid point
    ! mode(3) = 0 : nps is mesh point
    !           1 : nps is grid point
    
    logical :: forbitten
    integer,intent(in) :: np,nze,nps,nsa,mode(3)
    integer :: ips
    real(rkind) :: zetal, Bml, psiml, Fpsil, Gm, lorentzFactor, beta&
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
      zetal = zeta(nze)
    case(1)
      zetal = zetag(nze)
    end select

    select case(mode(3))
    case(0)
      if(zetal>=0.d0)then
        Bml = Bout(nps)
        dBmdpsi = (Boutg(nps+1)-Boutg(nps))/(psimg(nps+1)-psimg(nps))
      else
        Bml = Bin(nps)
        dBmdpsi = (Bing(nps+1)-Bing(nps))/(psimg(nps+1)-psimg(nps))
      end if
      psiml = psim(nps)
      Fpsil = Fpsi(nps)
      dFdpsi = (Fpsig(nps+1)-Fpsig(nps))/(psimg(nps+1)-psimg(nps))
      Gm = aefp(nsa)*Bml*psiml/(amfp(nsa)*vc*Fpsi(nps))

    case(1)
      if(zetal>=0.d0)then
        Bml = Boutg(nps)
        if(nps==1)then
          dBmdpsi = (-3*Boutg(1)+4*Boutg(2)-Boutg(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
          dFdpsi = (-3*Fpsig(1)+4*Fpsig(2)-Fpsig(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
        else if(nps==npsmax+1)then
          dBmdpsi = (3*Boutg(npsmax+1)-4*Boutg(npsmax)+Boutg(npsmax-1))&
                 /(3*psimg(npsmax+1)-4*psimg(npsmax)+psimg(npsmax-1))
          dFdpsi = (3*Fpsig(npsmax+1)-4*Fpsig(npsmax)+Fpsig(npsmax-1))&
                /(3*psimg(npsmax+1)-4*psimg(npsmax)+psimg(npsmax-1))
        else
          dBmdpsi = (Bout(nps)-Bout(nps-1))/(psim(nps)-psim(nps-1))
          dFdpsi = (Fpsi(nps)-Fpsi(nps-1))/(psim(nps)-psim(nps-1))
        end if

      else
        Bml = Bing(nps)
        if(nps==1)then
          dBmdpsi = (-3*Bing(1)+4*Bing(2)-Bing(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
          dFdpsi = (-3*Fpsig(1)+4*Fpsig(2)-Fpsig(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
        else if(nps==npsmax+1)then
          dBmdpsi = (3*Bing(npsmax+1)-4*Bing(npsmax)+Bing(npsmax-1))&
                 /(3*psimg(npsmax+1)-4*psimg(npsmax)+psimg(npsmax-1))
          dFdpsi = (3*Fpsig(npsmax+1)-4*Fpsig(npsmax)+Fpsig(npsmax-1))&
                /(3*psimg(npsmax+1)-4*psimg(npsmax)+psimg(npsmax-1))
        else
          dBmdpsi = (Bin(nps)-Bin(nps-1))/(psim(nps)-psim(nps-1))
          dFdpsi = (Fpsi(nps)-Fpsi(nps-1))/(psim(nps)-psim(nps-1))
        end if
      end if

      psiml = psimg(nps)
      Fpsil = Fpsig(nps)
      Gm = aefp(nsa)*Bml*psiml/(amfp(nsa)*vc*Fpsig(nps))
    end select

    stagnation_orbit = Gm*zetal/(zetal**2*dFdpsi/Fpsil-0.5*(1.d0+zetal**2)*dBmdpsi/Bml)

    if(lorentzFactor*beta<stagnation_orbit)then
      forbitten = .false.
    else
      forbitten = .true.
    end if

  end function forbitten

end module foworbitclassify
