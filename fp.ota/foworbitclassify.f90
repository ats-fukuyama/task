module foworbitclassify

  implicit none

  private

  public :: trapped, passing, forbitten, fow_trapped_boundary, fow_forbitten_boundary, fow_stagnation_type, fow_test

contains

  function trapped(np,nth,nr,nsa,mode)
    ! mode(1) = 0 :  np is mesh point
    !           1 :  np is grid point
    ! mode(2) = 0 : nth is mesh point
    !           1 : nth is grid point
    ! mode(3) = 0 : nr is mesh point
    !           1 : nr is grid point
    use fowcomm
    use fpcomm  

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
    use fowcomm
    use fpcomm
  
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
    use fowcomm
    use fpcomm
  
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
    ! WRITE(6,'(A,2ES12.4)') 'xil,stg_ob:', xil,stagnation_orbit
    stagnation_orbit = (1.d0+stagnation_orbit**(-2))**(-0.5d0)

    if ( beta<stagnation_orbit ) then
      forbitten = .false.
    else
      forbitten = .true.
    end if

  end function forbitten
  
  subroutine fow_trapped_boundary(upper_boundary)
    ! upper_boundary is maximum momentum of trapped particles for given psi_m, xi and particle species
    use fowcomm
    use fpcomm
  
    real(rkind) :: upper_boundary(:,:,:)
    integer :: np, nth, nr, nsa, ir
    real(rkind) :: dnr_bn, v_D_orbit
    real(rkind),allocatable :: psin(:,:), B_m(:,:), Gm(:,:,:), Bn(:,:)
  
    allocate(psin(nthmax,nrmax),B_m(nthmax,nrmax),Bn(nthmax,nrmax),Gm(nthmax,nrmax,nsamax))
  
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
        Bn(nth,nr) = B_m(nth,nr)/(1.d0-xi(nth)**2)
  
        if( Bn(nth,nr)<=Boutg(1) ) then
          do ir = 1, nrmax
            if( Boutg(ir+1)<=Bn(nth,nr) ) then
              dnr_bn = (Bn(nth,nr)-Boutg(ir))/(Boutg(ir+1)-Boutg(ir))
              psin(nth,nr) = psimg(ir)+(psimg(ir+1)-psimg(ir))*dnr_bn
              exit
            else if ( Bn(nth,nr)<Boutg(nrmax+1) ) then
              dnr_bn = (Bn(nth,nr)-Boutg(nrmax+1))/(Boutg(nrmax+1)-Boutg(nrmax))
              psin(nth,nr) = psimg(nrmax+1)+(psimg(nrmax+1)-psimg(nrmax))*dnr_bn
              exit
            end if
          end do
        else
          do ir = 1, nrmax
            if( Bn(nth,nr)<=Bing(ir+1) ) then
              dnr_bn = (Bn(nth,nr)-Bing(ir))/(Bing(ir+1)-Bing(ir))
              psin(nth,nr) = psimg(ir)+(psimg(ir+1)-psimg(ir))*dnr_bn
              exit
            else if ( Bn(nth,nr)>Bing(nrmax+1) ) then
              dnr_bn = (Bn(nth,nr)-Bing(nrmax+1))/(Bing(nrmax+1)-Bing(nrmax))
              psin(nth,nr) = psimg(nrmax+1)+(psimg(nrmax+1)-psimg(nrmax))*dnr_bn
              exit
            end if
          end do
        end if
    
      end do
    end do
  
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do nth = 1, nthmax
          if( Boutg(nrmax+1)<=Bn(nth,nr) .and. Bn(nth,nr)<=Bing(nrmax+1) .and. psin(nth,nr)<=psim(nr) ) then
            v_D_orbit = vc/sqrt(1.d0+xi(nth)**2/(Gm(nth,nr,nsa)*(1.d0-psin(nth,nr)/psim(nr)))**2)
            upper_boundary(nth,nr,nsa) = amfp(nsa)*v_D_orbit/sqrt(1.d0-v_D_orbit**2/vc**2)
            upper_boundary(nth,nr,nsa) = upper_boundary(nth,nr,nsa)/ptfp0(nsa) ! normalize
          else
            upper_boundary(nth,nr,nsa) = 0.d0
          end if
        end do
      end do
    end do
  
  
  end subroutine fow_trapped_boundary

  subroutine fow_forbitten_boundary(upper_boundary)
    ! upper_boundary is maximum momentum of not-forbitten particles for given psi_m, xi and particle species
    use fowcomm
    use fpcomm
  
    real(rkind),intent(out) :: upper_boundary(:,:,:)
    integer :: np,nth,nr,nsa
    integer :: ir
    real(rkind) :: v_stagnation_orbit
    real(rkind),allocatable :: B_m(:,:), Gm(:,:,:), dBmdpsi(:,:), dFdpsi(:)
      
    allocate(B_m(nthmax,nrmax),Gm(nthmax,nrmax,nsamax),dBmdpsi(nthmax,nrmax), dFdpsi(nrmax))
  
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

    call first_order_derivative(dFdpsi, Fpsi, psim)

    do nth = 1, nthmax
      call first_order_derivative(dBmdpsi(nth,:), B_m(nth,:), psim)
    end do
    
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do nth = 1, nthmax
          v_stagnation_orbit = Gm(nth,nr,nsa)*xi(nth)/(xi(nth)**2*dFdpsi(nr)/Fpsi(nr)&
                                        -0.5d0*(1.d0+xi(nth)**2)*dBmdpsi(nth,nr)/B_m(nth,nr)) ! LHS = gamma*beta 
          v_stagnation_orbit = vc*sqrt(v_stagnation_orbit**2/(1.d0+v_stagnation_orbit**2)) ! LHS = velocity of stagnation orbit
  
          upper_boundary(nth,nr,nsa) = amfp(nsa)*v_stagnation_orbit/(1.d0-v_stagnation_orbit**2/vc**2) ! momentum of stagnation orbit
          upper_boundary(nth,nr,nsa) = upper_boundary(nth,nr,nsa)/ptfp0(nsa) ! normalize
        end do
      end do
    end do
  
  end subroutine fow_forbitten_boundary

  subroutine fow_stagnation_type(xi_Xtype_boundary)
    use fpcomm
    use fowcomm
    use fpwrite
    implicit none
    ! F_p : ( dF / d(psi_p/psi_m) )/F 
    ! F_pp : ( d^2F / d(psi_p/psi_m)^2 )/F
    real(rkind), intent(out) :: xi_Xtype_boundary(:,:,:)
    integer :: nth, np, nr, nsa, ir
    real(rkind) :: C(3), F_p, F_pp, B_p, B_pp
    complex(rkind) :: z(2)
    real(rkind),allocatable :: dFdpsi(:), d2Fdpsi(:), dBmdpsi(:,:), d2Bmdpsi(:,:)&
                              ,F_(:), F__(:), B_(:,:), B__(:,:)

    allocate(dFdpsi(nrmax), d2Fdpsi(nrmax),dBmdpsi(nrmax,2), d2Bmdpsi(nrmax,2))

    call first_order_derivative(dFdpsi, Fpsi, psim)
    call second_order_derivative(d2Fdpsi, Fpsi, psim)
    
    call first_order_derivative(dBmdpsi(:,1), Bout, psim)
    call first_order_derivative(dBmdpsi(:,2), Bin, psim)

    call second_order_derivative(d2Bmdpsi(:,1), Bout, psim)
    call second_order_derivative(d2Bmdpsi(:,2), Bin, psim)

    call fpcsv1D(Fpsi,"./csv/Fpsi.csv")
    call fpcsv1D(dFdpsi,"./csv/dFdpsi.csv")
    call fpcsv1D(d2Fdpsi,"./csv/d2Fdpsi.csv")
    call fpcsv1D(Bout,"./csv/Bout.csv")
    call fpcsv1D(Bin,"./csv/Bin.csv")
    call fpcsv1D(dBmdpsi(:,1),"./csv/dBmdpsi_out.csv")
    call fpcsv1D(dBmdpsi(:,2),"./csv/dBmdpsi_in.csv")
    call fpcsv1D(d2Bmdpsi(:,1),"./csv/d2Bmdpsi_out.csv")
    call fpcsv1D(d2Bmdpsi(:,2),"./csv/d2Bmdpsi_in.csv")

    allocate(F_(nrmax),F__(nrmax),B_(nrmax,2), B__(nrmax,2))
    do nr = 1, nrmax
      F_(nr)=dFdpsi(nr)*psim(nr)/Fpsi(nr)
      F__(nr)=d2Fdpsi(nr)*psim(nr)**2/Fpsi(nr)
      B_(nr,1) = dBmdpsi(nr,1)*psim(nr)/Bout(nr)
      B__(nr,1) = d2Bmdpsi(nr,1)*psim(nr)**2/Bout(nr)
      B_(nr,2) = dBmdpsi(nr,2)*psim(nr)/Bin(nr)
      B__(nr,2) = d2Bmdpsi(nr,2)*psim(nr)**2/Bin(nr)
    end do

    call fpcsv1D(F_,"./csv/F_.csv")
    call fpcsv1D(F__,"./csv/F__.csv")
    call fpcsv1D(B_(:,1),"./csv/B_o.csv")
    call fpcsv1D(B__(:,1),"./csv/B__o.csv")
    call fpcsv1D(B_(:,2),"./csv/B_i.csv")
    call fpcsv1D(B__(:,2),"./csv/B__i.csv")

    do nr = 1, nrmax
      F_p = dFdpsi(nr)*psim(nr)/Fpsi(nr)
      F_pp = d2Fdpsi(nr)*psim(nr)**2/Fpsi(nr)

      ! calculate for xi > 0
      B_p = dBmdpsi(nr,1)*psim(nr)/Bout(nr)
      B_pp = d2Bmdpsi(nr,1)*psim(nr)**2/Bout(nr)

      C(3) = 4.d0*F_pp-4.d0*B_p*F_p-2.d0*B_pp+3.d0*B_p**2
      C(2) = 6.d0*B_p**2-2.d0*B_pp-4.d0*B_p*F_p
      C(1) = -1.d0*B_p**2
      ! write(*,*)"+nr",nr,F_p,F_pp,B_p,B_pp

      call solve_quadratic_equation(z, C)

      ! if ( aimag(z(2)) == 0.d0 ) then
      !   if ( real(z(2)) >= 1.d0 ) then
      !     xi_Xtype_boundary(nr,1) = 1.d0
      !   else if ( real(z(2)) <= 0.d0 ) then
      !     xi_Xtype_boundary(nr,1) = 0.d0
      !   else
      !     xi_Xtype_boundary(nr,1) = sqrt(real(z(2)))
      !   end if
      ! else
      !   xi_Xtype_boundary(nr,1) = 0.d0
      ! end if
      if ( aimag(z(2)) == 0.d0 ) then
        ! if ( real(z(2)) >= 1.d0 ) then
        !   xi_Xtype_boundary(nr,1,2) = 1.d0
        if ( real(z(2)) <= 0.d0 ) then
          xi_Xtype_boundary(nr,1,2) = 0.d0
        else
          xi_Xtype_boundary(nr,1,2) = sqrt(real(z(2)))
        end if

        ! if ( real(z(1)) >= 1.d0 ) then
        !   xi_Xtype_boundary(nr,1,1) = 1.d0
        if ( real(z(1)) <= 0.d0 ) then
          xi_Xtype_boundary(nr,1,1) = 0.d0
        else
          xi_Xtype_boundary(nr,1,1) = sqrt(real(z(1)))
        end if

      else
        xi_Xtype_boundary(nr,1,1) = 0.d0
        xi_Xtype_boundary(nr,1,2) = 0.d0
      end if

      ! calculate for xi < 0
      B_p = dBmdpsi(nr,2)*psim(nr)/Bin(nr)
      B_pp = d2Bmdpsi(nr,2)*psim(nr)**2/Bin(nr)

      C(3) = 4.d0*F_pp-4.d0*B_p*F_p-2.d0*B_pp+3.d0*B_p**2
      C(2) = 6.d0*B_p**2-2.d0*B_pp-4.d0*B_p*F_p
      C(1) = -1.d0*B_p**2

      call solve_quadratic_equation(z, C)

      if ( aimag(z(2)) == 0.d0 ) then
        ! if ( real(z(2)) >= 1.d0 ) then
        !   xi_Xtype_boundary(nr,2,2) = -1.d0
        if ( real(z(2)) <= 0.d0 ) then
          xi_Xtype_boundary(nr,2,2) = 0.d0
        else
          xi_Xtype_boundary(nr,2,2) = -sqrt(real(z(2)))
        end if

        ! if ( real(z(1)) >= 1.d0 ) then
        !   xi_Xtype_boundary(nr,2,1) = -1.d0
        if ( real(z(1)) <= 0.d0 ) then
          xi_Xtype_boundary(nr,2,1) = 0.d0
        else
          xi_Xtype_boundary(nr,2,1) = -sqrt(real(z(1)))
        end if

      else
        xi_Xtype_boundary(nr,2,1) = 0.d0
        xi_Xtype_boundary(nr,2,2) = 0.d0
      end if

    end do

  end subroutine fow_stagnation_type

  subroutine solve_quadratic_equation(z,C)
    ! solve C(3)*z**2+C(2)*z+C(1) = 0
    ! z(1) : -sqrt(D)
    ! z(2) : +sqrt(D)
    use fpcomm,only:rkind

    implicit none
    real(rkind),intent(in) :: C(3)
    complex(rkind),intent(out) :: z(2)
    real(rkind) :: D
    complex(rkind) :: ei=(0.d0,1.d0)

    z(1) = (0.d0, 0.d0)
    z(2) = (0.d0, 0.d0)

    D = C(2)**2-4.d0*C(1)*C(3)

    if ( D >= 0.d0 ) then
      z(1) = (-C(2)-sqrt(D))/(2.d0*C(3))+0.d0*ei
      z(2) = (-C(2)+sqrt(D))/(2.d0*C(3))+0.d0*ei
    else
      z(1) = (-C(2)-ei*sqrt(-D))/(2.d0*C(3))
      z(2) = (-C(2)+ei*sqrt(-D))/(2.d0*C(3))
    end if

  end subroutine solve_quadratic_equation

  subroutine fow_test
    use fpwrite
    implicit none
    INTEGER,PARAMETER :: imax=300
    real(8) :: y(imax), x(imax), dy(imax), d2y(imax)
    integer :: i

    do i = 1, imax
      x(i) = (dble(i)-0.5d0)/imax*9
      y(i) = log(x(i))
      dy(i) = 1.d0/x(i)
      d2y(i) = -1.d0/x(i)**2
    end do
    call fpcsv1D(y,"./csv/y_r.csv")
    call fpcsv1D(dy,"./csv/dy_r.csv")
    call fpcsv1D(d2y,"./csv/d2y_r.csv")
    call first_order_derivative(dy,y,x)
    call second_order_derivative(d2y,y,x)
    call fpcsv1D(dy,"./csv/dy_c.csv")
    call fpcsv1D(d2y,"./csv/d2y_c.csv")
    write(*,*)"kaijou",func_kaijou(2),func_kaijou(3),func_kaijou(4),func_kaijou(5)

  end subroutine fow_test

  subroutine nth_order_derivative(dfdx, f, x, n)
    ! calcurate dfdx, the n-th order derivative of f , n=1, 2, 3, or 4
    ! error order is O(dx**2)
    use fpcomm,only:rkind

    implicit none
    real(rkind),intent(out) :: dfdx(:)
    real(rkind),intent(in)  :: f(:), x(:)
    integer, intent(in) :: n
    real(rkind) :: s(5), A(4,4), B(4), C(5,5), D(5)
    integer :: imax, i, k, l

    imax = size(f)

    do i = 1, imax
      if ( i+2 <= imax .and. i-2 >= 1) then
        s(1) = x(i+2)-x(i)
        s(2) = x(i+1)-x(i)
        s(3) = x(i-1)-x(i)
        s(4) = x(i-2)-x(i)

        do k = 1, 4
          do l = 1, 4
            A(k,l) = s(k)**l/func_kaijou(l)
          end do
        end do
        B(1) = f(i+2)-f(i)
        B(2) = f(i+1)-f(i)
        B(3) = f(i-1)-f(i)
        B(4) = f(i-2)-f(i)

        call gauss_jordan(A, B, 4)

        dfdx(i) = B(n)

      else if ( i == 1 .or. i==2 ) then
        s(1) = x(i+1)-x(i)
        s(2) = x(i+2)-x(i)
        s(3) = x(i+3)-x(i)
        s(4) = x(i+4)-x(i)
        s(5) = x(i+5)-x(i)

        do k = 1, 5
          D(k) = f(i+k)-f(i)
          do l = 1, 5
            C(k,l) = s(k)**l/func_kaijou(l)
          end do
        end do

        call gauss_jordan(C, D, 5)

        dfdx(i) = B(n)

      else if ( i == imax .or. i==imax-1 ) then
        s(1) = x(i)-x(i-1)
        s(2) = x(i)-x(i-2)
        s(3) = x(i)-x(i-3)
        s(4) = x(i)-x(i-4)
        s(5) = x(i)-x(i-5)

        do k = 1, 5
          D(k) = f(i-k)-f(i)
          do l = 1, 5
            C(k,l) = (-1.d0*s(k))**l/func_kaijou(l)
          end do
        end do

        call gauss_jordan(C, D, 5)

        dfdx(i) = B(n)
      end if
    end do

  end subroutine nth_order_derivative

  subroutine first_order_derivative(dfdx, f, x)
    use fpcomm,only:rkind
    implicit none
    real(rkind),intent(out) :: dfdx(:)
    real(rkind),intent(in)  :: f(:), x(:)

    call nth_order_derivative(dfdx, f, x, 1)

  end subroutine

  subroutine second_order_derivative(dfdx, f, x)
    use fpcomm,only:rkind
    implicit none
    real(rkind),intent(out) :: dfdx(:)
    real(rkind),intent(in)  :: f(:), x(:)

    call nth_order_derivative(dfdx, f, x, 2)

  end subroutine


  ! subroutine first_order_derivative(dfdx, f, x)
  !   ! calcurate dfdx, the first order derivative of f
  !   ! error order is O(dx**2)
  !   use fpcomm,only:rkind

  !   implicit none
  !   real(rkind),intent(out) :: dfdx(:)
  !   real(rkind),intent(in)  :: f(:), x(:)
  !   real(rkind) :: s, t
  !   integer :: imax, i

  !   imax = size(f)

  !   do i = 1, imax
  !     if ( i /= imax .and. i /= 1) then
  !       t = x(i+1)-x(i)
  !       s = x(i)-x(i-1)
  !       dfdx(i) = (s**2*f(i+1)+(t**2-s**2)*f(i)-t**2*f(i-1))/(s*t*(s+t))
  !     else if ( i == 1 ) then
  !       s = x(2)-x(1)
  !       t = x(3)-x(2)
  !       dfdx(1) = ((s**2-(s+t)**2)*f(1)+(s+t)**2*f(2)-s**2*f(3))/(s*t*(s+t))
  !     else if ( i == imax ) then
  !       t = x(imax)-x(imax-1)
  !       s = x(imax-1)-x(imax-2)
  !       dfdx(imax) = (((s+t)**2-t**2)*f(imax)-(s+t)**2*f(imax-1)+t**2*f(imax-2))/(s*t*(s+t))
  !     end if
  !   end do

  ! end subroutine first_order_derivative

  ! subroutine second_order_derivative(d2fdx2, f, x)
  !   ! calcurate d2fdx2, the second order derivative of f
  !   ! error order is O(dx**2)

  !   use fpcomm,only:rkind

  !   implicit none
  !   real(rkind),intent(out) :: d2fdx2(:)
  !   real(rkind),intent(in)  :: f(:), x(:)
  !   real(rkind) :: s, t, r, v, w, A(3,3), B(3)
  !   integer :: imax, i

  !   imax = size(f)

  !   do i = 1, imax
  !     if ( i /= imax .and. i /= 1) then
  !       t = x(i+1)-x(i)
  !       s = x(i)-x(i-1)
  !       d2fdx2(i) = (s*f(i+1)-(s+t)*f(i)+t*f(i-1))/(0.5d0*s*t*(s+t))
  !     else if ( i == 1 ) then
  !       r = x(2)-x(1)
  !       s = x(3)-x(2)
  !       t = x(4)-x(3)
  !       v = r+s
  !       w = r+s+t

  !       A(1,1) = r
  !       A(1,2) = r**2/2
  !       A(1,3) = r**3/6

  !       A(2,1) = v
  !       A(2,2) = v**2/2
  !       A(2,3) = v**3/6

  !       A(3,1) = w
  !       A(3,2) = w**2/2
  !       A(3,3) = w**3/6

  !       B(1) = f(2)-f(1)
  !       B(2) = f(3)-f(1)
  !       B(3) = f(4)-f(1)

  !       call gauss_jordan(A, B, 3)

  !       d2fdx2(1) = B(2)

  !     else if ( i == imax ) then
  !       r = x(imax-2)-x(imax-3)
  !       s = x(imax-1)-x(imax-2)
  !       t = x(imax)-x(imax-1)
  !       v = s+t
  !       w = r+s+t

  !       A(1,1) = -1.d0*t
  !       A(1,2) = t**2/2
  !       A(1,3) = -1.d0*t**3/6

  !       A(2,1) = -1.d0*v
  !       A(2,2) = v**2/2
  !       A(2,3) = -1.d0*v**3/6

  !       A(3,1) = -1.d0*w
  !       A(3,2) = w**2/2
  !       A(3,3) = -1.d0*w**3/6

  !       B(1) = f(imax-1)-f(imax)
  !       B(2) = f(imax-2)-f(imax)
  !       B(3) = f(imax-3)-f(imax)

  !       call gauss_jordan(A, B, 3)

  !       d2fdx2(imax) = B(2)
        
  !     end if
  !   end do

  ! end subroutine second_order_derivative

  subroutine gauss_jordan(A, B, n)
    use fpcomm,only:rkind

    implicit none
    real(rkind) :: A(n,n), B(n)
    integer n,i,j,k
  
    do k = 1, n
      do j = k + 1, n
        a(k,j) = a(k,j) / a(k,k)
      end do
      b(k) = b(k) / a(k,k)
  
      do i = 1, n
        if ( i .ne. k ) then
          do j = k + 1, n
            a(i,j) = a(i,j) - a(i,k) * a(k,j)
          end do
          b(i) = b(i) - a(i,k) * b(k)
        end if
      end do
  
    end do
  
  end subroutine gauss_jordan
  
  recursive function func_kaijou(n) result(m)
    implicit none
    integer,intent(in) :: n
    integer :: m

    if(n == 1) then
      m = 1
    else
      m = n*func_kaijou(n-1)
    end if

  end function func_kaijou
  
end module foworbitclassify
