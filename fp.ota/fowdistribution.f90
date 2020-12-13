module fowdistribution
  implicit none
  private
  public :: fI_Maxwellian,fI_Maxwellian_sub
  public :: convert_fI_to_fu
  public :: moment_0th_order_COM, moment_2nd_order_COM
  public :: total_N, radial_particle_flux, effective_diffusion_cosfficient


contains

  function fI_Maxwellian(nth,np,nr,nsa) result(fI)
    use plprof
    use fpcomm
    use fowcomm
    use foworbit

    implicit none
    real(rkind) :: fI
    integer,intent(in) :: nth, np, nr, nsa
    real(rkind) :: ra_Bmin, theta_Bmin, rtfd0l, ptfd0l, fact, ex
    real(rkind) :: rnfd0l, rnfdl, rtfdl
    integer :: ns
    type(pl_plf_type),dimension(nsmax) :: plf

    ns = ns_nsa(nsa)

    call quantities_at_Bminimum(ra_Bmin, theta_Bmin, orbit_m(nth,np,nr,nsa))

    rtfd0l = ( ptpr(ns)+2.d0*ptpp(ns) )/3.d0
    ptfd0l = SQRT( rtfd0l*1.d3*aee*amfp(ns) )
    rnfd0l = pn(ns)


    if( model_ex_read_tn == 0 ) then
      call pl_prof(ra_Bmin, plf)
      rnfdl = plf(ns)%rn/rnfd0l
      rtfdl = ( plf(ns)%rtpr+2.d0*plf(ns)%rtpp)/3.d0

    else
      rnfdl=rn_temp(nr,ns)/rnfd0l
      rtfdl=rt_temp(nr,ns)

    end if

    fact = rnfdl / SQRT(2.d0*pi*rtfdl/rtfd0l)**3
    ex = pm(np,nsa)**2/(2.d0*rtfdl/rtfd0l)
    fI = fact*EXP( -1.d0*ex )! / JI(nth,np,nr,nsa)

  end function

  subroutine fI_Maxwellian_sub(fI, nsa)
    use plprof
    use fpcomm
    use fpsub
    use fowcomm
    use foworbit

    implicit none
    real(rkind),intent(out) :: fI(:,:,:)
    integer,intent(in) :: nsa
    real(rkind) :: ra_Bmin, theta_Bmin, rtfd0l, ptfd0l, fact, ex
    real(rkind) :: rnfd0l, rnfdl, rtfdl
    integer :: nth, np, nr, ns
    real(rkind) :: sumI, sumu
    type(pl_plf_type),dimension(nsmax) :: plf

    ns = ns_nsa(nsa)

    sumI = 0.d0
    sumu = 0.d0
    do nr = 1, nrmax
      do np = 1, npmax
        do nth = 1, nthmax

          call quantities_at_Bminimum(ra_Bmin, theta_Bmin, orbit_m(nth,np,nr,nsa))

          rtfd0l = ( ptpr(ns)+2.d0*ptpp(ns) )/3.d0
          ptfd0l = SQRT( rtfd0l*1.d3*aee*amfp(ns) )
          rnfd0l = pn(ns)
      
      
          if( model_ex_read_tn == 0 ) then
            call pl_prof(ra_Bmin, plf)
            rnfdl = plf(ns)%rn/rnfd0l
            rtfdl = ( plf(ns)%rtpr+2.d0*plf(ns)%rtpp)/3.d0
      
          else
            rnfdl=rn_temp(nr,ns)/rnfd0l
            rtfdl=rt_temp(nr,ns)
      
          end if
      
          fact = rnfdl / SQRT(2.d0*pi*rtfdl/rtfd0l)**3
          ex = pm(np,nsa)**2/(2.d0*rtfdl/rtfd0l)
          fI(nth, np, nr) = fact*EXP( -1.d0*ex )! / JI(nth,np,nr,nsa)
      
        end do
      end do
    end do

    do nr = 1, nrmax
      sumI = 0.d0
      do np = 1, npmax
        do nth = 1, nthmax
          sumI = sumI + fI(nth,np,nr)*JIR(nth,np,nr,nsa)*delp(nsa)*delthm(nth,np,nr,nsa)
        end do
      end do

      do np = 1, npmax
        do nth = 1, nthmax
          fI(nth,np,nr) = fI(nth,np,nr)*rn_temp(nr,ns)/rnfd0l/sumI
        end do
      end do

    end do

  end subroutine

  subroutine convert_fI_to_fu(fu_l, fI)

    use fpcomm
    use fowcomm

    implicit none
    real(rkind),dimension(nthmax,npmax,nrmax,nsamax),intent(in) :: fI
    real(rkind),dimension(nthmax,npmax,nrmax,nthpmax,nsamax),intent(out) :: fu_l

    integer :: nth,np,nr,nthp,nsa
    integer :: ir, ith, irr, irl, ithr, ithl
    real(rkind) :: psiml, thetaml, suml, sumI, volI
    real(rkind) :: dxm, dxp, dym, dyp, f11, f12, f21, f22, dd, vr, vp

    do nsa = 1, nsamax

      suml = 0.d0
      do nthp = 1,nthpmax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              if ( time_loss(nth,np,nr,nthp,nsa) /= 0 ) then ! loss region
                fu_l(nth,np,nr,nthp,nsa) = 0.d0

              else
                thetaml = thetam_local(nth,np,nr,nthp,nsa)
                psiml   = psim_local(nth,np,nr,nthp,nsa)
  
                ! search the cell include thetaml, psiml
                do ir = 1, nrmax-1
                  if ( psim(ir) <= psiml .and. psiml <= psim(ir+1) ) then
                    exit
                  end if
                end do
  
                do ith = 1, nthmax-1
                  if ( thetam(ith,np,ir,nsa) <= thetaml &
                      .and. thetaml <= thetam(ith+1,np,ir,nsa) ) then
                    exit
                  end if
                end do

                ! use bilinear intepolation
                ithr = min( ith+1, nthmax )
                ithl = ithr-1
                irr  = min( ir+1, nrmax )
                irl  = irr-1

                f11 = fI(ithl,np,irl,nsa)*JI(ithl,np,irl,nsa) &
                      *delthm(ithl,np,irl,nsa)*delp(nsa)*delps(irl)

                f12 = fI(ithl,np,irr,nsa)*JI(ithl,np,irr,nsa) &
                      *delthm(ithl,np,irr,nsa)*delp(nsa)*delps(irr)

                f21 = fI(ithr,np,irl,nsa)*JI(ithr,np,irl,nsa) &
                      *delthm(ithr,np,irl,nsa)*delp(nsa)*delps(irl)

                f22 = fI(ithr,np,irr,nsa)*JI(ithr,np,irr,nsa) &
                      *delthm(ithr,np,irr,nsa)*delp(nsa)*delps(irr)

                dxm = thetaml - thetam(ithl,np,irl,nsa)
                dxp = thetam(ithr,np,irl,nsa) - thetaml 
                dym = psiml - psim(irl)
                dyp = psim(irr) - psiml
                dd = (delps(irl)+delps(irr))/2.d0 * (delthm(ithl,np,irl,nsa)+delthm(ithr,np,irl,nsa))/2.d0

                fu_l(nth,np,nr,nthp,nsa) = ( dyp*dxp*f11 + dyp*dxm*f21 + dym*dxp*f12 + dym*dxm*f22 )/dd * fnorm(nsa)

                vr = twopi*rr*(rg(nr+1)-rg(nr))*rm(nr)*twopi/dble(nthpmax)

                fu_l(nth,np,nr,nthp,nsa) = fu_l(nth,np,nr,nthp,nsa) / (vr*volp(nth,np,nsa)) / RNFP0(NSA)*1.D20

                suml = suml + vr*volp(nth,np,nsa)*fu_l(nth,np,nr,nthp,nsa)
  
              end if
              
            end do
          end do
        end do
      end do

      sumI = 0.d0
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            sumI = sumI + delps(nr)*delp(nsa)*delthm(nth,np,nr,nsa) &
                            *JI(nth,np,nr,nsa)*fI(nth,np,nr,nsa)
          end do
        end do
      end do

      do nthp = 1, nthpmax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              vr = twopi*rr*(rg(nr+1)-rg(nr))*rm(nr)*twopi/dble(nthpmax)
              fu_l(nth,np,nr,nthp,nsa) = fu_l(nth,np,nr,nthp,nsa) * sumI / suml 
            end do
          end do
        end do
      end do

    end do

  end subroutine convert_fI_to_fu

  subroutine moment_0th_order_COM(M0, fI_in)
    use fpcomm
    use fowcomm
    use foworbit
    implicit none
    real(rkind),dimension(nthmax,npmax,nrmax,nsamax),intent(in) :: fI_in
    real(rkind),dimension(nrmax,nsamax),intent(out) :: M0
    integer :: nth, np, nr, nsa, nthp, ir
    real(rkind) :: sum_f, dVr, mean_ra, dVu, dVI


    do nsa = 1, nsamax
      do nr = 1, nrmax

        M0(nr,nsa) = 0.d0
        do np = 1, npmax
          do nth = 1, nthmax
            dVI = delthm(nth,np,nr,nsa)*delp(nsa)
            M0(nr,nsa) = M0(nr,nsa) + fI_in(nth,np,nr,nsa) * dVI * JIR(nth,np,nr,nsa)
          end do
        end do

      end do
    end do

  end subroutine

  subroutine moment_2nd_order_COM(M2, fI_in)
    use fpcomm
    use fowcomm
    use foworbit
    implicit none
    real(rkind),dimension(nthmax,npmax,nrmax,nsamax),intent(in) :: fI_in
    real(rkind),dimension(nrmax,nsamax),intent(out) :: M2
    integer :: nth, np, nr, nsa
    real(rkind) :: dVI, K, PV


    do nsa = 1, nsamax
      do nr = 1, nrmax

        M2(nr,nsa) = 0.d0
        do np = 1, npmax
          PV = sqrt(1.d0+theta0(nsa)*pm(np,nsa)**2)
          do nth = 1, nthmax
            K = (PV-1.d0)*amfp(nsa)*vc**2
            dVI = delthm(nth,np,nr,nsa)*delp(nsa)
            M2(nr,nsa) = M2(nr,nsa) + K*fI_in(nth,np,nr,nsa) * dVI * JIR(nth,np,nr,nsa)
          end do
        end do

      end do
    end do

  end subroutine

  subroutine radial_particle_flux(gamma_r, nsa)
    use fpcomm
    use fowcomm
    implicit none
    real(rkind),intent(out) :: gamma_r(:)
    integer,intent(in) :: nsa
    integer :: nr
    real(rkind) :: sumg, dndt, vprimem, vprimeg
    real(rkind),allocatable :: rnsl_prev(:,:)

    allocate(rnsl_prev(nrmax,nsamax))
    call moment_0th_order_COM(rnsl_prev, fnsm)

    sumg = 0.d0
    do nr = 2, nrmax
      dndt = ( rnsl(nr,nsa)-rnsl_prev(nr,nsa) )/delt
      vprimem=4.d0*pi**2*rr*(dble(nr)-0.5d0)/dble(nrmax)*ra
      vprimeg=4.d0*pi**2*rr*dble(nr)/dble(nrmax)*ra
      sumg = sumg + (-1.d0*dndt)*vprimem*delr*ra
      gamma_r(nr) = sumg/vprimeg
    end do

  end subroutine radial_particle_flux

  subroutine effective_diffusion_cosfficient(Deff, nsa)
    use fpcomm
    use fowcomm
    implicit none
    real(rkind),intent(out) :: Deff(:)
    integer,intent(in) :: nsa
    real(rkind),dimension(nrmax) :: g_r, dndr
    integer :: nr


    call radial_particle_flux(g_r, nsa)

    do nr = 1, nrmax
      if ( nr /= nrmax ) then
        dndr(nr) = (rnsl(nr+1,nsa) - rnsl(nr,nsa))/delr
      else
        dndr(nr) = (3.d0*rnsl(nr,nsa)-4.d0*rnsl(nr-1,nsa)+rnsl(nr-2,nsa))/delr/2.d0
      end if
      Deff(nr) = -1.d0*g_r(nr)/dndr(nr)
    end do


  end subroutine effective_diffusion_cosfficient

  subroutine total_N(N,fI,nsa)
    use fpcomm
    use fowcomm

    implicit none
    real(rkind),intent(in) :: fi(:,:,:,:)
    integer,intent(in) :: nsa
    real(rkind),intent(out) :: N
    integer :: nth,np,nr
    real(rkind) :: sumI, dVI

    sumI = 0.d0
    do nr = 1, nrmax
      do np = 1, npmax
        do nth = 1, nthmax
          dVI = delp(nsa)*delthm(nth,np,nr,nsa)*delps(nr)! * 2.d0*pi*pm(np,nsa)**2
          sumI = sumI+FI(nth,np,nr,nsa)*dVI*JI(nth,np,nr,nsa)
        end do
      end do
    end do

    N = sumI * fnorm(nsa)
    
  end subroutine total_N

  subroutine moment_0th_order_COMl(M0, fI_in)
    use fpcomm
    use fowcomm
    use foworbit
    implicit none
    real(rkind),dimension(nthmax,npmax,nrmax,nsamax),intent(in) :: fI_in
    real(rkind),dimension(nrmax,nsamax),intent(out) :: M0
    integer :: nth, np, nr, nsa, nthp, ir
    real(rkind) :: sum_f, dVr, mean_ra, dVu, dVI
    real(rkind),allocatable :: fu_l(:,:,:,:,:) 

    allocate(fu_l(nthmax,npmax,nrmax,nthpmax,nsamax))

    call convert_fI_to_fu(fu_l, fI_in)

    do nsa = 1, nsamax

      do ir = 1, nrmax
        M0(ir,nsa) = 0.d0
      end do

      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            mean_ra = get_mean_ra(orbit_m(nth,np,nr,nsa))
            
            do ir = 1, nrmax
              if ( rg(nr) <= mean_ra .and. mean_ra <= rg(nr+1) ) then
                dVI = JI(nth,np,nr,nsa)*delp(nsa)*delthm(nth,np,nr,nsa)*delps(nr) * fnorm(nsa)
                dVu = volp(nth,np,nsa)*twopi*rr*(rg(ir+1)-rg(ir))*rm(ir)*twopi/dble(nthpmax)! * RNFP0(NSA)*1.0d20

                dVr = twopi*RR*twopi*rm(nr)
                M0(ir,nsa) = M0(ir,nsa) + fI_in(nth,np,nr,nsa)*dVI/dVu * dVr * volp(nth,np,nsa)
              end if
            end do

          end do
        end do
      end do

    end do

    deallocate(fu_l)

  end subroutine

  subroutine moment_2nd_order_COMl(M2, fI_in)
    use fpcomm
    use fowcomm
    implicit none
    real(rkind),dimension(nthmax,npmax,nrmax,nsamax),intent(in) :: fI_in
    real(rkind),dimension(nrmax,nsamax),intent(out) :: M2 
    integer :: nth, np, nr, nsa, nthp
    real(rkind) :: sum_f, dV_r, PV
    real(rkind),allocatable :: fu_l(:,:,:,:,:) 

    allocate(fu_l(nthmax,npmax,nrmax,nthpmax,nsamax))

    call convert_fI_to_fu(fu_l, fI_in)

    do nsa = 1, nsamax
      do nr = 1, nrmax

        sum_f = 0.d0
        dV_r = 2.d0*pi*RR*2.d0*pi*rm(nr)
        do nthp = 1, nthpmax
          do np = 1, npmax

            PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSA)**2)
            do nth = 1, nthmax
              sum_f = sum_f + (amfp(nsa)*vc**2*(PV-1.d0))&
                              *fu_l(nth,np,nr,nthp,nsa)*dV_r*VOLP(nth,np,nsa)
            end do

          end do
        end do
        M2(nr,nsa) = sum_f*fnorm(nsa) ! [J]
        M2(nr,nsa) =  M2(nr,nsa)*1.D-6      ! [MJ]
        ! M2(nr,nsa) = sum_f*RNFP0(NSA)*1.D20/(aee*1.d-3) ! [keV]

      end do
    end do

    deallocate(fu_l)

  end subroutine

end module fowdistribution