module fowbounceaverage
  implicit none
  private
contains

  ! subroutine fow_calculate_jacobian
  !   use fpcomm
  !   use fowcomm
  ! end subroutine fow_calculate_jacobian

  subroutine fow_transformation_matrix(dIdu, th_in, p_in, B_in, F_in, np, nth, nr, nsa, mode)

    use fpcomm
    use fowcomm

    real(rkind),intent(out) :: dIdu(3,3)
    integer,intent(in) :: np, nsa, nth, nr, mode(3)
    real(rkind),intent(in) :: th_in, p_in, B_in, F_in

    ! elements of transformation matrix, dIdu
    real(rkind) ::  dpdp,  dzedp,  dpsdp,&
                dpdth, dzedth, dpsdth,&
                dpdr,  dzedr,  dpsdr
    real(rkind) :: A,B,C,D,E
    real(rkind) :: xil, Fpsil, Bml, dBmdpsi, dFoBdpsi ! dFoBdpsi=d(Fpsil/Bml)/dpsi_m

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
          dFoBdpsi = (-3*Fpsi(1)/Bout(1)+4*Fpsi(2)/Bout(2)-Fpsi(3)/Bout(3))/(-3*psim(1)+4*psim(2)-psim(3))
        else if ( nr==nrmax ) then
          dBmdpsi = (3*Bout(nrmax)-4*Bout(nrmax-1)+Bout(nrmax-2))&
                  /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
          dFoBdpsi = (3*Fpsi(nrmax)/Bout(nrmax)-4*Fpsi(nrmax-1)/Bout(nrmax-1)+Fpsi(nrmax-2)/Bout(nrmax-2))&
                  /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
        else
          dBmdpsi = (Bout(nr+1)-Bout(nr-1))/(psim(nr+1)-psim(nr-1))
          dFoBdpsi = (Fpsi(nr+1)/Bout(nr+1)-Fpsi(nr-1)/Bout(nr-1))/(psim(nr+1)-psim(nr-1))
        end if

      else if ( xil<0.d0) then
        Bml = Bin(nr)
        if ( nr==1 ) then
          dBmdpsi = (-3*Bin(1)+4*Bin(2)-Bin(3))/(-3*psim(1)+4*psim(2)-psim(3))
          dFoBdpsi = (-3*Fpsi(1)/Bin(1)+4*Fpsi(2)/Bin(2)-Fpsi(3)/Bin(3))/(-3*psim(1)+4*psim(2)-psim(3))
        else if ( nr==nrmax ) then
          dBmdpsi = (3*Bin(nrmax)-4*Bin(nrmax-1)+Bin(nrmax-2))&
                  /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
          dFoBdpsi = (3*Fpsi(nrmax)/Bin(nrmax)-4*Fpsi(nrmax-1)/Bin(nrmax-1)+Fpsi(nrmax-2)/Bin(nrmax-2))&
                  /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
        else
          dBmdpsi = (Bin(nr+1)-Bin(nr-1))/(psim(nr+1)-psim(nr-1))
          dFoBdpsi = (Fpsi(nr+1)/Bin(nr+1)-Fpsi(nr-1)/Bin(nr-1))/(psim(nr+1)-psim(nr-1))
        end if

      end if
      Fpsil = Fpsi(nr)

    case(1)
      if ( xil>=0.d0 ) then
        Bml = Boutg(nr)
        if ( nr==1 ) then
          dBmdpsi = (-3*Boutg(1)+4*Boutg(2)-Boutg(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
          dFoBdpsi = (-3*Fpsig(1)/Boutg(1)+4*Fpsig(2)/Boutg(2)-Fpsig(3)/Boutg(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
        else if ( nr==nrmax+1 ) then
          dBmdpsi = (3*Boutg(nrmax+1)-4*Boutg(nrmax)+Boutg(nrmax-1))&
                  /(3*psimg(nrmax+1)-4*psimg(nrmax)+psimg(nrmax-1))
          dFoBdpsi = (3*Fpsig(nrmax+1)/Boutg(nrmax+1)-4*Fpsig(nrmax)/Boutg(nrmax)+Fpsig(nrmax-1)/Boutg(nrmax-1))&
                /(3*psimg(nrmax+1)-4*psimg(nrmax)+psimg(nrmax-1))
        else
          dBmdpsi = (Boutg(nr+1)-Boutg(nr-1))/(psimg(nr+1)-psimg(nr-1))
          dFoBdpsi = (Fpsig(nr+1)/Boutg(nr+1)-Fpsig(nr-1)/Boutg(nr-1))/(psimg(nr+1)-psimg(nr-1))
        end if

      else
        Bml = Bing(nr)
        if ( nr==1 ) then
          dBmdpsi = (-3*Bing(1)+4*Bing(2)-Bing(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
          dFoBdpsi = (-3*Fpsig(1)/Bing(1)+4*Fpsig(2)/Bing(2)-Fpsig(3)/Bing(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
        else if ( nr==nrmax+1 ) then
          dBmdpsi = (3*Bing(nrmax+1)-4*Bing(nrmax)+Bing(nrmax-1))&
                  /(3*psimg(nrmax+1)-4*psimg(nrmax)+psimg(nrmax-1))
          dFoBdpsi = (3*Fpsig(nrmax+1)/Bing(nrmax+1)-4*Fpsig(nrmax)/Bing(nrmax)+Fpsig(nrmax-1)/Bing(nrmax-1))&
                /(3*psimg(nrmax+1)-4*psimg(nrmax)+psimg(nrmax-1))
        else
          dBmdpsi = (Bing(nr+1)-Bing(nr-1))/(psimg(nr+1)-psimg(nr-1))
          dFoBdpsi = (Fpsig(nr+1)/Bing(nr+1)-Fpsig(nr-1)/Bing(nr-1))/(psimg(nr+1)-psimg(nr-1))
        end if
      end if
      Fpsil = Fpsig(nr)

    end select

    A = F_in/B_in*cos(th_in)-Fpsil/Bml*xil
    B = -(1.d0-xil**2)/(2.d0*xil)*dBmdpsi/Bml
    C = dFoBdpsi*p_in*xil-AEFP(nsa)
    D = p_in*sin(th_in)/B_in*(Fpsil*cos(th_in)/xil-F_in)
    E = Bml/B_in*sin(th_in)*cos(th_in)/xil

    dpdp =1.d0
    dpdth=0.d0
    dpdr =0.d0
    dzedr=0.d0
    dpsdr=0.d0
    dpsdp = A/(B*dFoBdpsi*p_in+C)
    dzedp = A/(B*dFoBdpsi*p_in+C)*B
    dpsdth= D/(B*dFoBdpsi*p_in+C)
    dzedth= D/(B*dFoBdpsi*p_in+C)*B-E

    dIdu(1,1) = dpdp
    dIdu(1,2) = dzedp
    dIdu(1,3) = dpsdp
    dIdu(2,1) = dpdth
    dIdu(2,2) = dzedth
    dIdu(2,3) = dpsdth
    dIdu(3,1) = dpdr
    dIdu(3,2) = dzedr
    dIdu(3,3) = dpsdr

  end subroutine fow_transformation_matrix

end module fowbounceaverage