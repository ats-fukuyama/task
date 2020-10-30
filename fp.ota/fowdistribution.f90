module fowdistribution
  implicit none
  private
  public :: fow_distribution_maxwellian_inCOM
  public :: convert_fI_to_fu
  public :: moment_0th_order_inCOM, moment_2nd_order_inCOM


contains

  subroutine fow_distribution_maxwellian_inCOM(fm_I)
    ! calculate Maxwellian in COM space using the mean psi_p of orbits. fm_I is value on mesh points.
      use fowcomm
      use fpcomm
      use fpsub
      use fpwrite
      use foworbit,only:func_orbit_mean_ra

      real(rkind),intent(out) :: fm_I(:,:,:,:)
      real(rkind),allocatable :: U_fm_u(:,:), fm_u(:,:,:), fm_ux(:), U_RM(:,:), fx(:)
      real(rkind) :: mean_ra, dr, s, t
      integer :: np, nth, nr, nsa, ierr, ir
      real(rkind) :: deltaps, deltap, deltath, dV_p, sumI, sumu

      ierr = 0

      allocate(fm_u(npmax,nrmax,nsamax),U_fm_u(4,nrmax),fm_ux(nrmax))

      ! prepare spl-coef for func_orbit_mean_ra
      allocate(fx(nrmax),U_RM(4,nrmax))
      s = psim(2)-psim(1)
      t = psim(3)-psim(2)
      fx(1) = ((s**2-(s+t)**2)*rm(1)+(s+t)**2*rm(2)-s**2*rm(3))/(s*t*(s+t)) ! calculate edge dericvative
      t = psim(nrmax)-psim(nrmax-1)
      s = psim(nrmax-1)-psim(nrmax-2)
      fx(nrmax) = (((s+t)**2-t**2)*rm(nrmax)-(s+t)**2*rm(nrmax-1)+t**2*rm(nrmax-2))/(s*t*(s+t)) ! calculate edge dericvative
      call SPL1D(psim, rm, fx, U_RM,nrmax, 3, ierr)

      ! calculate Maxwellian in (p,r/a) space
      do nsa = 1, nsamax
        do nr = 1, nrmax
          do np = 1, npmax
            fm_u(np,nr,nsa) = FPMXWL(PM(NP,NSA),NR,NSA)
          end do
        end do
      end do

      do nsa = 1, nsamax
        do np = 1, npmax
          do nth = 1, nthmax
            do nr = 1, nrmax

              ! calculate mean minor radius of an orbit
              mean_ra = func_orbit_mean_ra(orbit_m(nth,np,nr,nsa),U_RM)

              ! calculate fm_I using linear interpolartion
              do ir = 2, nrmax
                  if( mean_ra <= rm(ir) .and. ir == 1 )then
                    dr = rm(1)-mean_ra
                    fm_I(nth,np,nr,nsa) = fm_u(np,1,nsa)-dr*(fm_u(np,1,nsa)-fm_u(np,2,nsa))/(rm(1)-rm(2))
                    exit
                  else if ( mean_ra <= rm(ir) .and. ir/=1 ) then
                  dr = rm(ir)-mean_ra
                  fm_I(nth,np,nr,nsa) = fm_u(np,ir,nsa)-dr*(fm_u(np,ir,nsa)-fm_u(np,ir-1,nsa))/(rm(ir)-rm(ir-1))
                  exit
                else if( mean_ra >= rm(nrmax))then
                  fm_I(nth,np,nr,nsa) = fm_u(np,nrmax,nsa)
                  exit
                end if
              end do

            end do
          end do
        end do
      end do

      do nsa = 1, nsamax

        sumI = 0.d0
        do nr = 1, nrmax
          deltaps = psimg(nr+1)-psimg(nr)
          do np = 1, npmax
            deltap = (pg(np+1,nsa)-pg(np,nsa))*ptfp0(nsa)
            do nth = 1, nthmax
              deltath = thetamg(nth+1,np,nr,nsa)-thetamg(nth,np,nr,nsa)
              sumI = sumI + fm_I(nth,np,nr,nsa)*deltap*deltath*deltaps*Jacobian_I(nth,np,nr,nsa)
            end do
          end do
        end do

        sumu = 0.d0
        do nr = 1, nrmax
          do np = 1, npmax

            dV_p = 0.d0
            do nth = 1, nthmax
              dV_p = dV_p + VOLP(nth,np,nsa)
            end do
            sumu = sumu + fm_u(np,nr,nsa)*Volr(nr)*dV_p

          end do
        end do
  
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              fm_I(nth,np,nr,nsa) = fm_I(nth,np,nr,nsa)*sumu/sumI
            end do
          end do
        end do
  
      end do  

  end subroutine fow_distribution_maxwellian_inCOM

  subroutine convert_fI_to_fu(fu_l, fI)

    use fpcomm
    use fowcomm

    implicit none
    real(rkind),dimension(nthmax,npmax,nrmax,nsamax),intent(in) :: fI
    real(rkind),dimension(nthmax,npmax,nrmax,nsamax,nthpmax),intent(out) :: fu_l

    integer :: nth,np,nr,nsa,nthp, ir, ir_min, ith,ip,isa, ith_min
    real(rkind) :: thetam_l, psim_l
    real(rkind) :: g, gmin, mu_l, Pzeta_l, pl, sign_thetam, difth_min, difth, Babsl
    real(rkind) :: sumI, sumu, deltap, deltaps, deltath, dV_r


    do nthp = 1, nthpmax
      do nsa = 1, nsamax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              if ( nthp == nthmax ) then
                Babsl = 0.5d0*Babs(nr,nthpmax)+0.5d0*Babs(nr,1)
              else
                Babsl = 0.5d0*Babs(nr,nthp)+0.5d0*Babs(nr,nthp+1)
              end if
              pl = pm(np,nsa)*ptfp0(nsa) 
              mu_l = pl**2*sinm(nth)**2/(2.d0*amfp(nsa)*Babsl)
              Pzeta_l = Fpsi(nr)/Babsl*pl*cosm(nth)-aefp(nsa)*psim(nr)
              gmin = 1.0d8

              do ir = 1, nrmax
                g = mu_l - pl**2/(2.d0*amfp(nsa)*Bout(ir))*(1.d0-(Bout(ir)/Fpsi(ir)/pl)**2*(Pzeta_l+aefp(nsa)*psim(ir))**2)
                if ( abs(g) <= abs(gmin) ) then
                  gmin = g
                  ir_min = ir
                  sign_thetam = 1.0
                end if

                g = mu_l - pl**2/(2.d0*amfp(nsa)*Bin(ir))*(1.d0-(Bin(ir)/Fpsi(ir)/pl)**2*(Pzeta_l+aefp(nsa)*psim(ir))**2)
                if ( abs(g) <= abs(gmin) ) then
                  gmin = g
                  ir_min = ir
                  sign_thetam = -1.0
                end if

              end do

              psim_l = psim(ir_min)
              if ( sign_thetam >= 0.d0 ) then
                thetam_l = acos((Pzeta_l+aefp(nsa)*psim(ir_min))/(Fpsi(ir_min)/Bout(ir_min)*pl))
              else
                thetam_l = acos((Pzeta_l+aefp(nsa)*psim(ir_min))/(Fpsi(ir_min)/Bin(ir_min)*pl))
              end if

              difth_min = 1.0d8

              do ith = 1, nthmax
                difth = abs(thetam_l-thetam(ith,np,ir_min,nsa))
                if ( difth <= difth_min ) then
                  difth_min = difth
                  ith_min = ith
                end if
              end do

              fu_l(nth,np,nr,nsa,nthp) = fI(ith_min, np, ir_min, nsa)

            end do
          end do
        end do
      end do
    end do


    do nsa = 1, nsamax

      sumI = 0.d0
      do nr = 1, nrmax
        deltaps = psimg(nr+1)-psimg(nr)
        do np = 1, npmax
          deltap = (pg(np+1,nsa)-pg(np,nsa))*ptfp0(nsa)
          do nth = 1, nthmax
            deltath = thetamg(nth+1,np,nr,nsa)-thetamg(nth,np,nr,nsa)
            sumI = sumI + fI(nth,np,nr,nsa)*deltap*deltath*deltaps*Jacobian_I(nth,np,nr,nsa)
          end do
        end do
      end do

      sumu = 0.d0
      do nthp = 1, nthpmax
        do nr = 1, nrmax

          dV_r = 2.d0*pi*RR*(rg(nr+1)-rg(nr))*rm(nr)*2.d0*pi/dble(nthpmax)
          do np = 1, npmax
            do nth = 1, nthmax
              sumu = sumu + fu_l(nth,np,nr,nsa,nthp)*dV_r*VOLP(nth,np,nsa)
            end do
          end do

        end do
      end do

      do nthp = 1, nthpmax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              fu_l(nth,np,nr,nsa,nthp) = fu_l(nth,np,nr,nsa,nthp)*sumI/sumu
            end do
          end do
        end do
      end do

    end do


  end subroutine convert_fI_to_fu

  subroutine moment_0th_order_inCOM(fI_in, M0)
    use fpcomm
    use fowcomm
    implicit none
    real(rkind),dimension(nthmax,npmax,nrmax,nsamax),intent(in) :: fI_in
    real(rkind),dimension(nrmax,nsamax),intent(out) :: M0 ! [/m^3]
    integer :: nth, np, nr, nsa, nthp
    real(rkind) :: sum_f, dV_r
    real(rkind),allocatable :: fu_l(:,:,:,:,:) 

    allocate(fu_l(nthmax,npmax,nrmax,nsamax,nthpmax))

    call convert_fI_to_fu(fu_l, fI_in)

    do nsa = 1, nsamax
      do nr = 1, nrmax

        sum_f = 0.d0
        dV_r = 2.d0*pi*RR*(rg(nr+1)-rg(nr))*rm(nr)*2.d0*pi/dble(nthpmax)
        do nthp = 1, nthpmax
          do np = 1, npmax
            do nth = 1, nthmax
              sum_f = sum_f + fu_l(nth,np,nr,nsa,nthp)*dV_r*VOLP(nth,np,nsa)
            end do
          end do
        end do
        M0(nr,nsa) = sum_f*RNFP0(NSA)*1.D20

      end do
    end do

    deallocate(fu_l)

  end subroutine

  subroutine moment_2nd_order_inCOM(fI_in, M2)
    use fpcomm
    use fowcomm
    implicit none
    real(rkind),dimension(nthmax,npmax,nrmax,nsamax),intent(in) :: fI_in
    real(rkind),dimension(nrmax,nsamax),intent(out) :: M2 ! [J]
    integer :: nth, np, nr, nsa, nthp
    real(rkind) :: sum_f, dV_r, PV
    real(rkind),allocatable :: fu_l(:,:,:,:,:) 

    allocate(fu_l(nthmax,npmax,nrmax,nsamax,nthpmax))

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
                              *fu_l(nth,np,nr,nsa,nthp)*dV_r*VOLP(nth,np,nsa)
            end do

          end do
        end do
        M2(nr,nsa) = sum_f*RNFP0(NSA)*1.D20

      end do
    end do

    deallocate(fu_l)

  end subroutine


end module fowdistribution