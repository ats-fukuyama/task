module check_neoclass
  private

  public :: output_neoclass

contains

  subroutine output_neoclass
    use fpcomm
    use fowcomm
    use fpwrite

    implicit none
    double precision,dimension(nrmax,nsamax) :: nu_D_theory, nu_star_theory, &
                                                Drhmrhm, nu_star
                                                 
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: wb, wb_model, omega_b, omega_b_model
    double precision,dimension(nrmax) :: nu_ei, nu_p, nu_b, q_ion, q_neo
    double precision,dimension(nrmax,nsamax) :: D_p, D_b, Drw
    double precision,dimension(nrmax,nsbmax,nsamax) :: nu_D
    double precision,dimension(3,3,nthmax,npmax,nrmax,nsamax) :: dIdu_m
    double precision,dimension(3,3,nthmax+1,npmax,nrmax,nsamax) :: dIdu_th
    double precision,dimension(3,3,nthmax,npmax+1,nrmax,nsamax) :: dIdu_p
    double precision,dimension(3,3,nthmax,npmax,nrmax+1,nsamax) :: dIdu_r
    integer nth, np, nr, nsa, nsb, mode(3)

    call system('mkdir -p dat')

    call integral_Drr(Drhmrhm)
    ! call D_neoclass(D_b, D_p, nu_ei, nu_p, nu_b)
    call D_random_walk_baverage(Drw)

    call banana_width_and_omega_bounce(wb, omega_b)
    call banana_width_and_omega_bounce_model(wb_model, omega_b_model)

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do nsb = 1, nsbmax
          nu_D(nr,nsb,nsa) = rnud(nr,nsb,nsa)*rnsl(nr,nsb)/rnfp0(nsb)
        end do
      end do
    end do

    call fptxt3D(nu_D,"dat/nud.txt")

    call fptxt2D(Drhmrhm,"dat/Drhmrhm.txt")

    call fptxt4D(wb,"dat/banana_width.txt")
    call fptxt4D(wb_model,"dat/banana_width_model.txt")
    call fptxt4D(omega_b,"dat/omega_bounce.txt")
    call fptxt4D(omega_b_model,"dat/omega_bounce_model.txt")

    ! call fptxt2D(D_p,"dat/D_plateau.txt")
    ! call fptxt2D(D_b,"dat/D_banana.txt")
    call fptxt2D(Drw,"dat/D_random_walk.txt")
    ! call fptxt1D(nu_ei,"dat/nu_Spitz.txt")
    ! call fptxt1D(nu_p,"dat/nu_p.txt")
    ! call fptxt1D(nu_b,"dat/nu_b.txt")

    mode = [0,0,0]
    call mean_transmatrix(dIdu_m,mode)
    call fptxt4D(dIdu_m(1,1,:,:,:,:),"dat/dpdp_m.txt")
    call fptxt4D(dIdu_m(2,1,:,:,:,:),"dat/dpdth_m.txt")
    call fptxt4D(dIdu_m(3,1,:,:,:,:),"dat/dpdr_m.txt")
    call fptxt4D(dIdu_m(1,2,:,:,:,:),"dat/dthdp_m.txt")
    call fptxt4D(dIdu_m(2,2,:,:,:,:),"dat/dthdth_m.txt")
    call fptxt4D(dIdu_m(3,2,:,:,:,:),"dat/dthdr_m.txt")
    call fptxt4D(dIdu_m(1,3,:,:,:,:),"dat/drdp_m.txt")
    call fptxt4D(dIdu_m(2,3,:,:,:,:),"dat/drdth_m.txt")
    call fptxt4D(dIdu_m(3,3,:,:,:,:),"dat/drdr_m.txt")
    mode = [1,0,0]
    call mean_transmatrix(dIdu_th,mode)
    call fptxt4D(dIdu_th(1,1,:,:,:,:),"dat/dpdp_th.txt")
    call fptxt4D(dIdu_th(2,1,:,:,:,:),"dat/dpdth_th.txt")
    call fptxt4D(dIdu_th(3,1,:,:,:,:),"dat/dpdr_th.txt")
    call fptxt4D(dIdu_th(1,2,:,:,:,:),"dat/dthdp_th.txt")
    call fptxt4D(dIdu_th(2,2,:,:,:,:),"dat/dthdth_th.txt")
    call fptxt4D(dIdu_th(3,2,:,:,:,:),"dat/dthdr_th.txt")
    call fptxt4D(dIdu_th(1,3,:,:,:,:),"dat/drdp_th.txt")
    call fptxt4D(dIdu_th(2,3,:,:,:,:),"dat/drdth_th.txt")
    call fptxt4D(dIdu_th(3,3,:,:,:,:),"dat/drdr_th.txt")
    mode = [0,0,1]
    call mean_transmatrix(dIdu_r,mode)
    call fptxt4D(dIdu_r(1,1,:,:,:,:),"dat/dpdp_r.txt")
    call fptxt4D(dIdu_r(2,1,:,:,:,:),"dat/dpdth_r.txt")
    call fptxt4D(dIdu_r(3,1,:,:,:,:),"dat/dpdr_r.txt")
    call fptxt4D(dIdu_r(1,2,:,:,:,:),"dat/dthdp_r.txt")
    call fptxt4D(dIdu_r(2,2,:,:,:,:),"dat/dthdth_r.txt")
    call fptxt4D(dIdu_r(3,2,:,:,:,:),"dat/dthdr_r.txt")
    call fptxt4D(dIdu_r(1,3,:,:,:,:),"dat/drdp_r.txt")
    call fptxt4D(dIdu_r(2,3,:,:,:,:),"dat/drdth_r.txt")
    call fptxt4D(dIdu_r(3,3,:,:,:,:),"dat/drdr_r.txt")
    mode = [0,1,0]
    call mean_transmatrix(dIdu_p,mode)
    call fptxt4D(dIdu_p(1,1,:,:,:,:),"dat/dpdp_p.txt")
    call fptxt4D(dIdu_p(2,1,:,:,:,:),"dat/dpdth_p.txt")
    call fptxt4D(dIdu_p(3,1,:,:,:,:),"dat/dpdr_p.txt")
    call fptxt4D(dIdu_p(1,2,:,:,:,:),"dat/dthdp_p.txt")
    call fptxt4D(dIdu_p(2,2,:,:,:,:),"dat/dthdth_p.txt")
    call fptxt4D(dIdu_p(3,2,:,:,:,:),"dat/dthdr_p.txt")
    call fptxt4D(dIdu_p(1,3,:,:,:,:),"dat/drdp_p.txt")
    call fptxt4D(dIdu_p(2,3,:,:,:,:),"dat/drdth_p.txt")
    call fptxt4D(dIdu_p(3,3,:,:,:,:),"dat/drdr_p.txt")

  end subroutine output_neoclass

  subroutine integral_Drr(Drr_out)
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out) :: Drr_out
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom
    double precision Drr_, dVI, sumVI, JIl_p, JIl_m
    integer nth,np,nr,nsa,ns

    do nsa = 1, nsamax
      do np = 1, npmax
        do nth = 1, nthmax
          call first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
        end do
      end do
    end do

    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        Drr_out(nr,nsa) = 0.d0
        sumVI = 0.d0
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit
          do nth = 1, nthmax

            if ( nr == 1 ) then
              JIl_m = JI(nth,np,1,nsa)
              JIl_p = ( JI(nth,np,2,nsa)+JI(nth,np,1,nsa) )*0.5d0
            else if ( nr == nrmax ) then
              JIl_m = ( JI(nth,np,nrmax,nsa)+JI(nth,np,nrmax-1,nsa) )*0.5d0
              JIl_p = JI(nth,np,nrmax,nsa)
            else
              JIl_m = ( JI(nth,np,nr,nsa)+JI(nth,np,nr-1,nsa) )*0.5d0
              JIl_p = ( JI(nth,np,nr,nsa)+JI(nth,np,nr+1,nsa) )*0.5d0
            end if

            dVI = delp(ns)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)
            Drr_ = ( Drrfow(nth,np,nr+1,nsa)/JIl_p+Drrfow(nth,np,nr,nsa)/JIl_m )*0.5d0

            Drr_out(nr,nsa) = Drr_out(nr,nsa) + Drr_*dfdrhom(nth,np,nr,nsa)*dVI
            sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
          end do
        end do
        Drr_out(nr,nsa) = Drr_out(nr,nsa)/sumVI
      end do
    end do

  end subroutine integral_Drr

  subroutine ionHeatFlux(q_ion)
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nrmax),intent(out) :: q_ion
    double precision,dimension(nthmax,npmax,nrmax) :: dfdrhom
    double precision :: Drr_, dVI, K, pv
    integer :: nth,np,nr,nsa,ns

    do np = 1, npmax
      do nth = 1, nthmax
        call first_order_derivative(dfdrhom(nth,np,:), fnsp(nth,np,:,2), rm)
        do nr = 1, nrmax
          dfdrhom(nth,np,nr) = dfdrhom(nth,np,nr)/RA
        end do
      end do
    end do

    do nr = 1, nrmax
      q_ion(nr) = 0.d0
      do np = 1, npmax
        if ( pm(np,2) > fact_bulk ) exit
        do nth = 1, nthmax
          dVI = delp(2)*delthm(nth,np,nr,2)*JIR(nth,np,nr,2)
          Drr_ = ( Drrfow(nth,np,nr+1,2)+Drrfow(nth,np,nr,2) )*0.5d0/JI(nth,np,nr,2)
          pv = SQRT(1.d0+theta0(2)*pm(np,2)**2)
          K = amfp(2)*vc**2*(pv-1.d0)

          q_ion(nr) = q_ion(nr) - K*Drr_*dfdrhom(nth,np,nr)*dVI*1.d20
        end do
      end do
      write(*,*)"ion",q_ion(nr)
    end do

  end subroutine ionHeatFlux

  subroutine D_neoclass(D_banana, D_plateau, nu_ei, nu_p, nu_b)
    use fowcomm
    use fpcomm

    implicit none
    double precision,dimension(nrmax),intent(out) :: D_banana, D_plateau, nu_ei, nu_p, nu_b
    integer :: nr, nsb, nssb
    double precision :: dens ,Te, Vt, omega_e, Baxis, rho_e, epst
    double precision :: factSpitz

    Baxis = Bing(1)
    omega_e = aee*Baxis/ame
    factSpitz = aee**4/( 93.d0*eps0**2*SQRT(ame) )
    do nr = 1, nrmax
      Te = rwsl(nr,1)*1.d6/( 1.5d0*rnsl(nr,1)*1.d20 ) ! [J]
      nu_ei(nr) = 0.d0
      do nsb = 1, nsbmax
        nssb = ns_nsb(nsb)
        if ( nssb == 1 ) cycle
        dens = rnsl(nr,nsb)*1.d20
        nu_ei(nr) = nu_ei(nr) + pz(nssb)**2*dens*lnlam(nr,nsb,1)
      end do
      nu_ei(nr) = nu_ei(nr)*factSpitz/Te**1.5d0
      ! nu_ei = rnud(nr,2,1)

      Vt = SQRT( Te/ame )
      rho_e = Vt/omega_e
      epst = rm(nr)*RA/RR
      nu_p(nr) = Vt/( safety_factor(nr)*RR )
      nu_b(nr) = epst**1.5d0*nu_p(nr)

      D_banana (nr) = ( safety_factor(nr)*rho_e )**2 * nu_ei(nr) / epst**1.5d0
      D_plateau(nr) = ( safety_factor(nr)*rho_e )**2 * nu_p(nr)
    end do

  end subroutine D_neoclass

  subroutine D_random_walk(Drw)
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: Drw
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: w_b, omega_b, dfdrhom
    double precision,dimension(npmax,nrmax,nsamax)        :: nud
    double precision step_len, step_time, Drwlocal
    double precision dVI, sumVI
    integer nth, np, nr, nsa

    call banana_width_and_omega_bounce(w_b, omega_b)
    call nu_deflection(nud)

    do nsa = 1, nsamax
      do np = 1, npmax
        do nth = 1, nthmax
          call first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
        end do
      end do
    end do

    do nsa = 1, nsamax
      do nr = 1, nrmax
        Drw(nr,nsa) = 0.d0
        sumVI = 0.d0

        do np = 1, npmax
          step_time = 1.d0/nud(np,nr,nsa)
          do nth = 1, nthmax
            step_len = w_b(nth,np,nr,nsa)/RA
            Drwlocal = step_len**2/step_time
            
            dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)
            sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
            Drw(nr,nsa) = Drw(nr,nsa) + Drwlocal*dfdrhom(nth,np,nr,nsa)*dVI

          end do
        end do

        Drw(nr,nsa) = Drw(nr,nsa)/sumVI

      end do
    end do

  end subroutine D_random_walk

  subroutine D_random_walk_baverage(Drw)
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: Drw
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: w_b, omega_b, dfdrhom
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
    double precision,dimension(npmax,nrmax,nsamax) :: nud
    double precision step_len, step_time
    double precision dVI, sumVI
    integer nth, np, nr, nsa

    call banana_width_and_omega_bounce(w_b, omega_b)
    call nu_deflection(nud)

    do nsa = 1, nsamax
      do np = 1, npmax
        do nth = 1, nthmax
          call first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
        end do
      end do
    end do

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          step_time = 1.d0/nud(np,nr,nsa)
          do nth = 1, nthmax
            step_len = w_b(nth,np,nr,nsa)
            Drwlocal(nth,np,nr,nsa) = step_len**2/step_time
          end do
        end do
      end do
    end do

    call bounce_average_for_Drw(Drwba, Drwlocal)

    do nsa = 1, nsamax
      do nr = 1, nrmax
        Drw(nr,nsa) = 0.d0
        sumVI = 0.d0
        do np = 1, npmax
          do nth = 1, nthmax
            dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)

            Drw(nr,nsa) = Drw(nr,nsa)&
                          + Drwba(nth,np,nr,nsa)*dfdrhom(nth,np,nr,nsa)*dVI
            sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
          end do
        end do
        Drw(nr,nsa) = Drw(nr,nsa)/sumVI
      end do
    end do

  end subroutine D_random_walk_baverage

  subroutine nu_deflection(nud)
    use fpcomm
    use fowcomm
    
    implicit none
    double precision,dimension(npmax,nrmax,nsamax),intent(out) :: nud
    double precision x, fx, gx
    double precision nudb, fact, vthb, Tb, va, Ta, pv
    double precision nuBra, factBra, numerator, denominator
    integer nth,np,nr,nsa,nsb,nssa,nssb

    fact = 3.d0/2.d0*SQRT( pi/2.d0 )
    factBra = 1.d0/( 12.d0*pi**1.5d0*eps0**2 )

    do nsa = 1, nsamax
      nssa = ns_nsa(nsa)
      do nr = 1, nrmax
        do np = 1, npmax
          nud(np,nr,nsa) = 0.d0
          pv = SQRT(1.d0+theta0(nsa)*pm(np,nsa)**2)
          va = vc * SQRT( 1.d0-pv**(-2) )
          do nsb = 1, nsbmax
            nssb = ns_nsb(nsa)
            Tb = rwsl(nr,nsb)*1.d6/( 1.5d0*rnsl(nr,nsb)*1.d20 )
            vthb = SQRT( 2.d0*Tb/amfp(nsb) )

            numerator  = rnsl(nr,nsb)*1.d20*( aefp(nsa)*aefp(nsb) )**2 * lnlam(nr,nsb,nsa)
            denominator= SQRT( amfp(nsb) )*Tb**1.5d0*( amfp(nsa)/amfp(nsb) )**2
            nuBra = factBra*numerator/denominator

            x = va/vthb
            call gosakn(x,fx,gx)

            nudb = fact*nuBra*(fx-gx)/x**3

            nud(np,nr,nsa) = nud(np,nr,nsa)+nudb

          end do
        end do
      end do
    end do

  end subroutine

  subroutine gosakn(x,fx,gx)
    use fpcomm,only:pi
    implicit none
    real(8),intent(in)  :: x
    real(8),intent(out) :: fx,gx
    real(8) :: f,fx1,fx2
    real(8) :: rh
    DATA RH / 0.70710678118654752440D+00/

    f= DERF(X)*0.5D0
    fx=2.d0*f
    fx1=2.d0/sqrt(pi)*exp(-x**2)
    fx2=-2.d0*x*fx1
    gx=0.d0
    if(abs(x).gt.1.e-10) then
        gx=(fx-x*fx1)/(2.d0*x**2)
    end if

    return

  end subroutine gosakn

  subroutine bounce_average_for_Drw(Dout, Din)
    use fpcomm
    use fowcomm
    use fowcoef
    implicit none
    double precision,dimension(nthmax,npmax,nrmax,nsamax),intent(out) :: Dout
    double precision,dimension(nthmax,npmax,nrmax,nsamax),intent(in)  :: Din
    double precision,dimension(3,3,max_stp) :: dIdul
    double precision,allocatable :: U(:,:,:,:,:,:,:,:), Drwl(:,:,:,:,:)
    double precision dt, taup, cpitch_ob, psip_ob, thetap_ob, Drw_ob
    type(orbit) ob
    integer nth, np, nr, nsa, nstp, nstpmax, nthp, mode(3)
    write(*,*)"b"

    allocate(U(4,4,4,nthmax,nrmax,nthpmax,npmax,nsamax))
    allocate(Drwl(nthmax,npmax,nrmax,nthpmax,nsamax))

    mode = [0,0,0]

    do nsa = 1, nsamax
      do nthp = 1, nthpmax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              Drwl(nth,np,nr,nthp,nsa) = Din(nth,np,nr,nsa)
            end do
          end do
        end do
      end do
    end do

    do nsa = 1, nsamax
      call make_U_Dxy(U, Drwl, 'm', nsa)
    end do

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            ob = orbit_m(nth,np,nr,nsa)
            nstpmax = ob%nstp_max
            taup = ob%time(nstpmax)
        
            call transformation_matrix(dIdul, ob, nth, np, nr, nsa, mode)

            Dout(nth,np,nr,nsa) = 0.d0

            do nstp = 2, nstpmax
              dt = ob%time(nstp)-ob%time(nstp-1)
              cpitch_ob = ob%costh(nstp)
              psip_ob   = ob%psip(nstp)
              thetap_ob = ob%thetap(nstp)
              call interpolate_D_unlessZero(Drw_ob, U(:,:,:,:,:,:,np,nsa), &
                    1.d0, cpitch_ob, psip_ob, thetap_ob)
              
              Dout(nth,np,nr,nsa) = Dout(nth,np,nr,nsa)&
                                  + Drw_ob*dIdul(3,3,nstp)**2*dt
            end do

            Dout(nth,np,nr,nsa) = Dout(nth,np,nr,nsa)/taup
        
          end do
        end do
      end do
    end do


  end subroutine

  subroutine banana_width_and_omega_bounce(w_b, omega_b)
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nthmax,npmax,nrmax,nsamax),intent(out) :: w_b, omega_b
    double precision,allocatable :: drdt(:), U(:,:)
    type(orbit) ob
    integer nth, np, nr, nsa, nstp, nstpmax, ierr
    double precision r_min

    ierr = 0

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            ob = orbit_m(nth,np,nr,nsa)
            nstpmax = ob%nstp_max
            omega_b(nth,np,nr,nsa) = 2.d0*pi/ob%time(nstpmax)
            allocate(drdt(nstpmax))
            allocate(U(4, nstpmax))
            call spl1D(ob%time, ob%r, drdt, U, nstpmax, 0, ierr)
            call spl1DF(ob%time(nstpmax)*0.5d0, r_min, ob%time, U ,nstpmax, ierr)
            w_b(nth,np,nr,nsa) = ( rm(nr)-r_min )*RA
            deallocate(drdt)
            deallocate(U)
          end do
        end do
      end do
    end do

  end subroutine banana_width_and_omega_bounce

  subroutine banana_width_and_omega_bounce_model(w_b, omega_b)
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nthmax,npmax,nrmax,nsamax),intent(out) :: w_b, omega_b
    double precision :: rho, eps,sign_prev,sign,Baxis,vl,pv,B_theta,thetapl
    type(orbit) :: ob
    logical :: isTrapped
    integer :: nth,np,nr,nsa,ns,nstp,nstpmax

    Baxis = Boutg(1)

    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        eps = rm(nr)*RA/RR
        do np = 1, npmax
          pv = SQRT(1.d0+theta0(nsa)*pm(np,nsa)**2)
          vl = vc*SQRT(1.d0-1.d0/pv**2)
          do nth = 1, nthmax
            ob = orbit_m(nth,np,nr,nsa)
            nstpmax = ob%nstp_max

            isTrapped = .false.
            if ( ob%costh(1) < 0.d0 ) then
              sign_prev = -1.d0
            else if ( ob%costh(1) > 0.d0 ) then
              sign_prev = 1.d0
            else
              sign_prev = 0.d0
            end if

            do nstp = 2, nstpmax
              if ( ob%costh(nstp) /= 0.d0 ) then
                sign = ob%costh(nstp) / ABS( ob%costh(nstp) )
              else
                sign = 0.d0
              end if

              if ( sign_prev*sign > 0.d0 ) then
                cycle
              else if ( sign_prev*sign < 0.d0 ) then
                isTrapped = .True.
              else
                sign_prev = sign
              end if
            end do

            if ( isTrapped ) then
              rho = pm(np,ns)*ptfp0(ns)/ABS( pz(ns)*aee*Baxis )
              w_b(nth,np,nr,nsa) = safety_factor(nr)*rho/SQRT(eps)
              omega_b(nth,np,nr,nsa) = SQRT(eps)*vl/( safety_factor(nr)*RR )
            else
              if ( pz(ns)*COS( thetam(nth,np,nr,nsa) ) > 0.d0 ) then
                thetapl = 0.d0
              else
                thetapl = pi
              end if
              B_theta = safety_factor(nr)/( RR+RA*rm(nr)*COS( thetapl ) )*dpsimdr(nr)
              w_b(nth,np,nr,nsa) = 0.d0
              omega_b(nth,np,nr,nsa) = B_theta*vl*COS( thetam(nth,np,nr,nsa) )&
                                      /( Baxis*rm(nr)*RA )
              omega_b(nth,np,nr,nsa) = ABS( omega_b(nth,np,nr,nsa) )
              ! omega_b(nth,np,nr,nsa) = 0.d0
            end if

          end do
        end do
      end do
    end do

  end subroutine banana_width_and_omega_bounce_model

  subroutine mean_transmatrix(dIdu,mode)
    use fpcomm
    use fowcomm
    use fowcoef
    implicit none
    double precision,dimension(:,:,:,:,:,:),intent(out) :: dIdu
    integer,dimension(3),intent(in) :: mode
    double precision,dimension(3,3,max_stp) :: dIdul
    double precision :: dt, taup
    integer :: nth,np,nr,nsa,nstp,nstpmax,i,j
    type(orbit) :: ob

    do nsa = 1, nsamax
      do nr = 1, nrmax+mode(3)
        do np = 1, npmax+mode(2)
          do nth = 1, nthmax+mode(1)
            do j = 1, 3
              do i = 1, 3
                dIdu(i,j,nth,np,nr,nsa) = 0.d0
              end do
            end do
          end do
        end do
      end do
    end do

    do nsa = 1, nsamax
      do nr = 1+mode(3), nrmax+mode(3)
        do np = 1+mode(2), npmax+mode(2)
          do nth = 1, nthmax+mode(1)

            if ( mode(1) == 0 .and. mode(2) == 0 .and. mode(3) == 0 ) then
              ob = orbit_m(nth,np,nr,nsa)
            else if ( mode(1) == 1 .and. mode(2) == 0 .and. mode(3) == 0 ) then
              ob = orbit_th(nth,np,nr,nsa)
              if ( nth == nth_stg(nsa) ) cycle
            else if ( mode(1) == 0 .and. mode(2) == 1 .and. mode(3) == 0 ) then
              ob = orbit_p(nth,np,nr,nsa)
            else if ( mode(1) == 0 .and. mode(2) == 0 .and. mode(3) == 1 ) then
              ob = orbit_r(nth,np,nr,nsa)
            end if

            nstpmax = ob%nstp_max
            taup = ob%time(nstpmax)
            call transformation_matrix(dIdul, ob, nth, np, nr, nsa, mode)

            do nstp = 2, nstpmax
              dt = ( ob%time(nstp)-ob%time(nstp-1) )/taup
              do j = 1, 3
                do i = 1, 3
                  dIdu(i,j,nth,np,nr,nsa) = dIdu(i,j,nth,np,nr,nsa) + dIdul(i,j,nstp)*dt
                end do
              end do

            end do

          end do
        end do
      end do
    end do

  end subroutine  mean_transmatrix

end module check_neoclass
