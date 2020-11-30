module fowdistribution
  implicit none
  private
  public :: fow_Maxwellian_COM
  public :: convert_fI_to_fu
  public :: moment_0th_order_COM, moment_2nd_order_COM


contains

  subroutine fow_Maxwellian_COM(fm_I)
    ! calculate Maxwellian in COM space using the mean psi_p of orbits. fm_I is value on mesh points.
      use fowcomm
      use fpcomm
      use fpsub
      use fpwrite
      use foworbit

      real(rkind),intent(out) :: fm_I(:,:,:,:)
      real(rkind),allocatable :: fm_u(:,:)
      real(rkind) :: mean_ra, dr, s, t
      integer :: np, nth, nr, nsa, ierr, ir
      real(rkind) :: dVp, dVI, sumI, sumu

      ierr = 0

      allocate(fm_u(npmax,nrmax))

      ! calculate Maxwellian in (p,r/a) space
      sumu = 0.d0
      do nsa = 1, nsamax

        do nr = 1, nrmax
          do np = 1, npmax

            dVp = 0.d0
            do nth = 1, nthmax
              dVp = dVp+volp(nth,np,nsa)*ptfp0(nsa)**3
            end do

            fm_u(np,nr) = FPMXWL(PM(NP,NSA),NR,NSA) * dVp * volr(nr)*ra**2
            sumu = sumu + fm_u(np,nr)

          end do
        end do

        sumI = 0.d0
        do np = 1, npmax
          do nth = 1, nthmax
            do nr = 1, nrmax

              ! fm_I(nth,np,nr,nsa) = 0.d0 in forbitten region
              if ( nth == nth_forbitten(nsa) ) then
                fm_I(nth,np,nr,nsa) = 0.d0
                cycle
              end if

              ! calculate mean minor radius of an orbit
              mean_ra = get_mean_ra( orbit_m(nth,np,nr,nsa) )

              ! linear-interpolartion to mean_ra
              do ir = 2, nrmax
                  if( mean_ra <= rm(ir) .and. ir == 1 )then
                    dr = rm(1)-mean_ra
                    fm_I(nth,np,nr,nsa) = fm_u(np,1)-dr*(fm_u(np,1)-fm_u(np,2))/(rm(1)-rm(2))
                    exit
                  else if ( mean_ra <= rm(ir) .and. ir/=1 ) then
                    dr = rm(ir)-mean_ra
                    fm_I(nth,np,nr,nsa) = fm_u(np,ir)-dr*(fm_u(np,ir)-fm_u(np,ir-1))/(rm(ir)-rm(ir-1))
                    exit
                  else if( mean_ra >= rm(nrmax))then
                    fm_I(nth,np,nr,nsa) = fm_u(np,nrmax)
                  exit
                end if
              end do

              dVI = Jacobian_I(nth,np,nr,nsa) * delthm(nth,np,nr,nsa) * delp(nsa) * delps(nr) * ptfp0(nsa)
              fm_I(nth,np,nr,nsa) = fm_I(nth,np,nr,nsa) / dVI

              sumI = sumI + fm_I(nth,np,nr,nsa)*dVI

            end do
          end do
        end do

        write(*,*)"sumu/sumI",sumu/sumI
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              fm_I(nth,np,nr,nsa) = fm_I(nth,np,nr,nsa) *sumu/sumI
            end do
          end do
        end do

      end do ! nsa loop

  end subroutine fow_Maxwellian_COM

  subroutine convert_fI_to_fu(fu_l, fI)

    use fpcomm
    use fowcomm

    implicit none
    real(rkind),dimension(nthmax,npmax,nrmax,nsamax),intent(in) :: fI
    real(rkind),dimension(nthmax,npmax,nrmax,nthpmax,nsamax),intent(out) :: fu_l

    integer :: nth,np,nr,nthp,nsa
    integer :: ir, ith, irr, irl, ithr, ithl
    real(rkind) :: psiml, thetaml, sum_ful, sum_fI, volI
    real(rkind) :: dxm, dxp, dym, dyp, f11, f12, f21, f22, dd, vr, vp

    do nsa = 1, nsamax

      sum_ful = 0.d0
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

                f11 = fI(ithl,np,irl,nsa)*Jacobian_I(ithl,np,irl,nsa) &
                      *delthm(ithl,np,irl,nsa)*delp(nsa)*delps(irl)

                f12 = fI(ithl,np,irr,nsa)*Jacobian_I(ithl,np,irr,nsa) &
                      *delthm(ithl,np,irr,nsa)*delp(nsa)*delps(irr)

                f21 = fI(ithr,np,irl,nsa)*Jacobian_I(ithr,np,irl,nsa) &
                      *delthm(ithr,np,irl,nsa)*delp(nsa)*delps(irl)

                f22 = fI(ithr,np,irr,nsa)*Jacobian_I(ithr,np,irr,nsa) &
                      *delthm(ithr,np,irr,nsa)*delp(nsa)*delps(irr)

                dxm = thetaml - thetam(ithl,np,irl,nsa)
                dxp = thetam(ithr,np,irl,nsa) - thetaml 
                dym = psiml - psim(irl)
                dyp = psim(irr) - psiml
                dd = (delps(irl)+delps(irr))/2.d0 * (delthm(ithl,np,irl,nsa)+delthm(ithr,np,irl,nsa))/2.d0

                fu_l(nth,np,nr,nthp,nsa) = ( dyp*dxp*f11 + dyp*dxm*f21 + dym*dxp*f12 + dym*dxm*f22 )/dd * ptfp0(nsa)

                if ( nthp == 1 ) then
                  vr = (rg(nr+1)-rg(nr))*ra * rm(nr)*theta_p(2)/2.d0 * twopi*rr
                else if ( nthp == nthpmax ) then
                  vr = (rg(nr+1)-rg(nr))*ra * rm(nr)*(twopi-0.5d0*theta_p(nthp)-0.5d0*theta_p(nthp-1))/2.d0 * twopi*rr
                else
                  vr = (rg(nr+1)-rg(nr))*ra * rm(nr)*(theta_p(nthp+1)-theta_p(nthp-1))/2.d0 * twopi*rr
                end if

                fu_l(nth,np,nr,nthp,nsa) = fu_l(nth,np,nr,nthp,nsa) / (vr*volp(nth,np,nsa)*ptfp0(nsa)**3)

                sum_ful = sum_ful + vr*volp(nth,np,nsa)*ptfp0(nsa)**3*fu_l(nth,np,nr,nthp,nsa)
  
              end if
              
            end do
          end do
        end do
      end do

      sum_fI = 0.d0
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            sum_fI = sum_fI + delps(nr)*delp(nsa)*ptfp0(nsa)*delthm(nth,np,nr,nsa) &
                            *Jacobian_I(nth,np,nr,nsa)*fI(nth,np,nr,nsa)
          end do
        end do
      end do

      do nthp = 1, nthpmax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              fu_l(nth,np,nr,nthp,nsa) = fu_l(nth,np,nr,nthp,nsa) * sum_fI / sum_ful
            end do
          end do
        end do
      end do

    end do

  end subroutine convert_fI_to_fu

  subroutine moment_0th_order_COM(M0, fI_in)
    use fpcomm
    use fowcomm
    implicit none
    real(rkind),dimension(nthmax,npmax,nrmax,nsamax),intent(in) :: fI_in
    real(rkind),dimension(nrmax,nsamax),intent(out) :: M0
    integer :: nth, np, nr, nsa, nthp
    real(rkind) :: sum_f, dV_r
    real(rkind),allocatable :: fu_l(:,:,:,:,:) 

    allocate(fu_l(nthmax,npmax,nrmax,nthpmax,nsamax))

    call convert_fI_to_fu(fu_l, fI_in)

    do nsa = 1, nsamax
      do nr = 1, nrmax

        sum_f = 0.d0
        dV_r = 2.d0*pi/dble(nthpmax)
        do nthp = 1, nthpmax
          do np = 1, npmax
            do nth = 1, nthmax
              sum_f = sum_f + fu_l(nth,np,nr,nthp,nsa)*dV_r*VOLP(nth,np,nsa)
            end do
          end do
        end do
        M0(nr,nsa) = sum_f*RNFP0(NSA)*1.D20*ptfp0(nsa)**3
        M0(nr,nsa) = M0(nr,nsa)*1.0d-20

      end do
    end do

    deallocate(fu_l)

  end subroutine

  subroutine moment_2nd_order_COM(M2, fI_in)
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
        dV_r = 2.d0*pi*RR*(rg(nr+1)-rg(nr))*rm(nr)*2.d0*pi/dble(nthpmax)
        do nthp = 1, nthpmax
          do np = 1, npmax

            PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSA)**2)
            do nth = 1, nthmax
              sum_f = sum_f + (amfp(nsa)*vc**2*(PV-1.d0))&
                              *fu_l(nth,np,nr,nthp,nsa)*dV_r*VOLP(nth,np,nsa)
            end do

          end do
        end do
        M2(nr,nsa) = sum_f*RNFP0(NSA)*1.D20*ptfp0(nsa)**3 ! [J]
        M2(nr,nsa) =  M2(nr,nsa)*1.D-6      ! [MJ]
        ! M2(nr,nsa) = sum_f*RNFP0(NSA)*1.D20/(aee*1.d-3) ! [keV]

      end do
    end do

    deallocate(fu_l)

  end subroutine

end module fowdistribution