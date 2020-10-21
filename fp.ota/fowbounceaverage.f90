module fowbounceaverage
  private
contains

  subroutine fow_dcoef_baverage(D_out, D_in, orbit_in, variable)

    use fpcomm
    use fowcomm
    use foworbitclassify

    implicit none
    real(rkind),intent(out) :: D_out(:,:,:,:,:)
    real(rkind),intent(in) ::  D_in(:,:,:,:,:)
    type(orbit),intent(in) :: orbit_in(:,:,:,:)
    character(*) :: variable
    integer :: nth, np, nr, nsa, mode(3), nstp, nstpmax
    real(rkind) :: dIdu(3,3), sign_xi
    real(rkind),allocatable :: dBmdpsi(:,:), dFdpsi(:)

    ! determine the row of diffusion tensor
    if ( variable == 'xi' ) then
      mode=[1,0,0]
    else if ( variable == 'p' ) then
      mode = [0,1,0]
    else if ( variable == 'psim' ) then
      mode = [0,0,1]
    end if

    ! if xi > 0, then second dimension of dBmdpsi = 1 and Bm = Bout
    ! if xi > 0, then second dimension of dBmdpsi = 2 and Bm = Bin
    allocate(dBmdpsi(nrmax+mode(3),2), dFdpsi(nrmax+mode(3)))
    select case(mode(3))
    case(0)
      call first_order_derivative(dFdpsi, Fpsi, psim)
      call first_order_derivative(dBmdpsi(:,1), Bout, psim)
      call first_order_derivative(dBmdpsi(:,2), Bin, psim)
    case(1)
      call first_order_derivative(dFdpsi, Fpsi, psim)
      call first_order_derivative(dBmdpsi(:,1), Bout, psim)
      call first_order_derivative(dBmdpsi(:,2), Bin, psim)
    end select

    ! execute integral D_in over orbit periods
    ! D_out = \int_0^taup D_in dt
    do nsa = 1, nsamax
      do nr = 1, nrmax+mode(3)
        do np = 1, npmax+mode(2)
          do nth = 1, nthmax+mode(1)

            select case(mode(1))
            case(0)
              if ( xi(nth) /= 0.d0 ) then
                sign_xi = xi(nth)/abs(xi(nth))
              else
                sign_xi = 0.d0
              end if
            case(1)
              if ( xig(nth) /= 0.d0 ) then
                sign_xi = xig(nth)/abs(xig(nth))
              else
                sign_xi = 0.d0
              end if
            end select

            if ( forbitten(nth,np,nr,nsa,mode) ) then
              continue
            else
              D_out(nth,np,nr,nsa,1) = 0.d0 ! column for xi
              D_out(nth,np,nr,nsa,2) = 0.d0 ! column for p
              D_out(nth,np,nr,nsa,3) = 0.d0 ! column for psi_m

              nstpmax = orbit_in(nth,np,nr,nsa)%nstp_max
              do nstp = 2, nstpmax

                if ( sign_xi > 0.d0) then
                  call fow_transformation_matrix&
                  (dIdu, orbit_in(nth,np,nr,nsa), nth, np, nr, nsa, nstp, mode, dBmdpsi(:,1), dFdpsi)
                else 
                  call fow_transformation_matrix&
                  (dIdu, orbit_in(nth,np,nr,nsa), nth, np, nr, nsa, nstp, mode, dBmdpsi(:,2), dFdpsi)
                end if

                D_out(nth,np,nr,nsa,1) = D_out(nth,np,nr,nsa,1) +
                D_out(nth,np,nr,nsa,1) = D_out(nth,np,nr,nsa,2) +
                D_out(nth,np,nr,nsa,1) = D_out(nth,np,nr,nsa,3) +
                
              end do
            end if
          end do
        end do
      end do
    end do

  end subroutine fow_dcoef_baverage

  subroutine fow_transformation_matrix(dIdu, orbit_in, nth_in, np_in, nr_in, nsa_in, nstp_in, mode, dBmdpsi, dFdpsi)

    use fpcomm
    use fowcomm

    implicit none

    real(rkind),intent(out) :: dIdu(3,3)
    type(orbit),intent(in) :: orbit_in
    integer,intent(in) :: nth_in, np_in, nr_in, nsa_in, nstp_in, mode(3)
    real(rkind),intent(in) :: dBmdpsi(:), dFdpsi(:)

    ! elements of transformation matrix, dIdu
    real(rkind) ::  dpdp,  dzedp,  dpsdp,&
                dpdth, dzedth, dpsdth,&
                dpdr,  dzedr,  dpsdr
    real(rkind) :: A,B,C,D,E
    real(rkind) :: pitch_angle, Fl, Bl, xil, pl, dFoBdpsi

    select case(mode(1))
    case(0)
      xil = xi(nth_in)
    case(1)
      xil = xig(nth_in)
    end select

    select case(mode(2))
    case(0)
      pl = pm(np_in)
    case(1)
      pl = pg(np_in)
    end select

    select case(mode(3))
    case(0)
      Fl = Fpsi(nr_in)
      if ( xil >= 0.d0 ) then
        Bl = Bout(nr_in)
      else
        Bl = Bin(nr_in)
      end if
    case(1)
      Fl = Fpsig(nr_in)
      if ( xil >= 0.d0 ) then
        Bl = Boutg(nr_in)
      else
        Bl = Bing(nr_in)
      end if
    end select

    dFoBdpsi = (dFdpsi(nr_in)*Bl-Fl*dBmdpsi(nr_in))/Bl**2

    pitch_angle = acos( orbit_in%vpara(nstp_in)/sqrt(orbit_in%vpara(nstp_in)**2+orbit_in%vperp(nstp_in)**2))

    A = orbit_in%RB_phi(nstp_in)/orbit_in%Babs(nstp_in)*cos(pitch_angle)-Fl/Bl*xil
    B = -(1.d0-xil**2)/(2.d0*xil)*dBmdpsi(nr_in)/Bl
    C = dFoBdpsi*pl*xil-AEFP(nsa_in)
    D = pl*sin(pitch_angle)/orbit_in%babs(nstp_in)*(Fl*cos(pitch_angle)/xil-orbit_in%RB_phi(nstp_in))
    E = Bl/orbit_in%babs(nstp_in)*sin(pitch_angle)*cos(pitch_angle)/xil

    dpdp =1.d0
    dpdth=0.d0
    dpdr =0.d0
    dzedr=0.d0
    dpsdr=0.d0
    dpsdp = A/(B*dFoBdpsi*pl+C)
    dzedp = A/(B*dFoBdpsi*pl+C)*B
    dpsdth= D/(B*dFoBdpsi*pl+C)
    dzedth= D/(B*dFoBdpsi*pl+C)*B-E

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