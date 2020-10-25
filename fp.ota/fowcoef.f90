module fowcoef
  private
contains

  subroutine fow_calulate_coef
  end subroutine fow_calulate_coef

  subroutine fow_bounce_average

    use fpcomm
    use fowcomm
    use fowprep,only:first_order_derivative

    implicit none
    integer :: nth, np, nr, nsa, mode(3), nstp, nstpmax
    real(rkind),allocatable :: dBmdpsi(:,:), dFdpsi(:), dBmgdpsi(:,:), dFgdpsi(:)
    real(rkind),allocatable :: dIdu(:,:,:), D_in(:,:), F_in(:)
    real(rkind),allocatable,dimension(:,:,:,:) :: Dppfow, Dptfow, Dprfow,&
                                                  Dtpfow, Dttfow, Dtrfow,&
                                                  Drpfow, Drtfow, Drrfow,&
                                                  Fpfow,  Ftfow,  Frfow

    real(rkind),allocatable,dimension(:) :: Dpporg, Dptorg, Fporg,&
                                            Dtporg, Dttorg, Ftorg

    ! if xi > 0, then second dimension of dBmdpsi = 1 and Bm = Bout
    ! if xi > 0, then second dimension of dBmdpsi = 2 and Bm = Bin
    allocate(dBmdpsi(nrmax,2), dFdpsi(nrmax))
    allocate(dBmgdpsi(nrmax+1,2), dFgdpsi(nrmax+1))

    call first_order_derivative(dFdpsi, Fpsi, psim)
    call first_order_derivative(dBmdpsi(:,1), Bout, psim)
    call first_order_derivative(dBmdpsi(:,2), Bin, psim)
    call first_order_derivative(dFgdpsi, Fpsig, psimg)
    call first_order_derivative(dBmdpsi(:,1), Boutg, psimg)
    call first_order_derivative(dBmdpsi(:,2), Bing, psimg)

    allocate(Dppfow(nthmax,npmax+1,nrmax,nsamax),Dptfow(nthmax,npmax+1,nrmax,nsamax),Dprfow(nthmax,npmax+1,nrmax,nsamax))
    allocate(Dtpfow(nthmax+1,npmax,nrmax,nsamax),Dttfow(nthmax+1,npmax,nrmax,nsamax),Dtrfow(nthmax+1,npmax,nrmax,nsamax))
    allocate(Drpfow(nthmax,npmax,nrmax+1,nsamax),Drtfow(nthmax,npmax,nrmax+1,nsamax),Drrfow(nthmax,npmax,nrmax+1,nsamax))
    allocate(Fpfow(nthmax,npmax+1,nrmax,nsamax),Ftfow(nthmax,npmax+1,nrmax,nsamax),Frfow(nthmax,npmax+1,nrmax,nsamax))


    ! calculate Dpp, Dpt, Dpr, Fp
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax+1
          do nth = 1, nthmax

            nstpmax = orbit_p(nth,np,nr,nsa)%nstp_max

            allocate(Dpporg(nstpmax),Dptorg(nstpmax),Fporg(nstpmax))
            allocate(Dtporg(nstpmax),Dttorg(nstpmax),Ftorg(nstpmax))

            allocate(dIdu(3,3,nstpmax), D_in(3,nstpmax), F_in(nstpmax))

            if ( thetam(nth,np,nr,nsa) <= pi/2.d0 ) then
              call transformation_matrix(dIdu, orbit_p(nth,np,nr,nsa), nth, np, nr, nsa, [0,1,0], dBmdpsi(:,1), dFdpsi)
            else
              call transformation_matrix(dIdu, orbit_p(nth,np,nr,nsa), nth, np, nr, nsa, [0,1,0], dBmdpsi(:,2), dFdpsi)
            end if
            
            do nstp = 1, nstpmax
              D_in(1, nstp) = Dpporg(nstp)
              D_in(2, nstp) = Dptorg(nstp)*dIdu(1,2,nstp)+Dptorg(nstp)*dIdu(2,2,nstp)
              D_in(3, nstp) = Dpporg(nstp)*dIdu(1,3,nstp)+Dptorg(nstp)*dIdu(2,3,nstp)
              F_in(nstp) = Fporg(nstp)
            end do

            ! Dpp
            call orbit_integral(Dppfow(nth,np,nr,nsa),D_in(1,:),orbit_p(nth,np,nr,nsa))
            ! Dpt
            call orbit_integral(Dptfow(nth,np,nr,nsa),D_in(2,:),orbit_p(nth,np,nr,nsa))
            ! Dpr
            call orbit_integral(Dprfow(nth,np,nr,nsa),D_in(3,:),orbit_p(nth,np,nr,nsa))
            ! Fp
            call orbit_integral(Fpfow(nth,np,nr,nsa),F_in,orbit_p(nth,np,nr,nsa))

            deallocate(Dpporg, Dptorg, Dtporg, Dttorg, Fporg, Ftorg)
            deallocate(dIdu, D_in, F_in)

          end do
        end do
      end do
    end do


  end subroutine fow_bounce_average

  subroutine calculate_ZOW_coef(Dppzow, Dptzow, Dtpzow, Dttzow, Fpzow, Ftzow, nth_in, np_in, nr_in, nsa_in, mode)

    use fpcomm
    use fowcomm

    implicit none
        
    real(rkind),intent(out) :: Dppzow(:), Dptzow(:), Dtpzow(:), Dttzow(:), Fpzow(:), Ftzow(:)
    integer,intent(in) :: nth_in, np_in, nr_in, nsa_in, mode(3)

  end subroutine

  subroutine transformation_matrix(dIdu, orbit_in, nth_in, np_in, nr_in, nsa_in, mode, dBmdpsi, dFdpsi)

    use fpcomm
    use fowcomm
    use foworbit,only:func_orbit_F

    implicit none

    real(rkind),intent(out) :: dIdu(:,:,:)
    type(orbit),intent(in) :: orbit_in
    integer,intent(in) :: nth_in, np_in, nr_in, nsa_in, mode(3)
    real(rkind),intent(in) :: dBmdpsi(:), dFdpsi(:)

    ! elements of transformation matrix, dIdu
    real(rkind) ::  dpdp,  dxidp,  dpsdp,&
                dpdth, dxidth, dpsdth,&
                dpdr,  dxidr,  dpsdr
    real(rkind) :: A,B,C,D,E
    real(rkind) :: Fl, Bl, xil, pl, dFoBdpsi, Fob, Bob, thetaob
    integer :: nstp, nstpmax

    select case(mode(1))
    case(0)
      xil = cos(thetam(nth_in,np_in,nr_in,nsa_in))
    case(1)
      xil = cos(thetamg(nth_in,np_in,nr_in,nsa_in))
    end select

    select case(mode(2))
    case(0)
      pl = pm(np_in,nsa_in)
    case(1)
      pl = pg(np_in,nsa_in)
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

    nstpmax = orbit_in%nstp_max

    do nstp = 1, nstpmax
      Fob = func_orbit_F(orbit_in, nstp, nr_in)
      Bob = orbit_in%Babs(nstp)
      thetaob = orbit_in%theta(nstp)
  
      dFoBdpsi = (dFdpsi(nr_in)*Bl-Fl*dBmdpsi(nr_in))/Bl**2
  
      A =Fob/Bob*cos(thetaob)-Fl/Bl*xil
      B = -(1.d0-xil**2)/(2.d0*xil)*dBmdpsi(nr_in)/Bl
      C = dFoBdpsi*pl*xil-AEFP(nsa_in)
      D = pl*sin(thetaob)/Bob*(Fl*cos(thetaob)/xil-Fob)
      E = Bl/Bob*sin(thetaob)*cos(thetaob)/xil
  
      dpdp =1.d0
      dpdth=0.d0
      dpdr =0.d0
      dxidr=0.d0
      dpsdr=0.d0
      dpsdp = A/(B*dFoBdpsi*pl+C)
      dxidp = A/(B*dFoBdpsi*pl+C)*B
      dpsdth= D/(B*dFoBdpsi*pl+C)
      dxidth= D/(B*dFoBdpsi*pl+C)*B-E
  
      dIdu(1,1,nstp) = dpdp
      dIdu(1,2,nstp) = dxidp  / sqrt(1.d0-xil**2)*(-1.d0) ! convert dxi/dp  to dthetam/dp
      dIdu(1,3,nstp) = dpsdp
      dIdu(2,1,nstp) = dpdth
      dIdu(2,2,nstp) = dxidth / sqrt(1.d0-xil**2)*(-1.d0) ! convert dxi/dth to dthetam/dth
      dIdu(2,3,nstp) = dpsdth
      dIdu(3,1,nstp) = dpdr
      dIdu(3,2,nstp) = dxidr
      dIdu(3,3,nstp) = dpsdr  
      
    end do

  end subroutine transformation_matrix

  subroutine orbit_integral(D_out, D_in, orbit_in)
    ! D_out = int_0^tau_p D_in / tau_p dt
    ! tau_p is poloidal period of orbit_in
    use fpcomm
    use fowcomm

    real(rkind),intent(out) :: D_out
    real(rkind),intent(in) :: D_in(:) ! nstpmax
    type(orbit),intent(in) :: orbit_in

    integer :: nstp, nstpmax
    real(rkind) :: dt

    nstpmax = orbit_in%nstp_max

    D_out = 0.d0
    do nstp = 2, nstpmax
      dt = orbit_in%time(nstp)-orbit_in%time(nstp-1)
      D_out = D_out + D_in(nstp)*dt
    end do

    D_out = D_out/orbit_in%time(nstpmax)

  end subroutine orbit_integral

end module fowcoef