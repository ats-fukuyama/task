! check_neoclass.f90
! [2022/12/9]
! ***************************************
!  Calculation of diffusion coefficients
! ***************************************
! made by ota / modified by anzai
! ver.4

module check_neoclass
  private

  public :: output_neoclass

contains

!===================================================
! COMMON modules
!===================================================
  subroutine output_neoclass
  ! !----------------------------------------
  ! module for file out put as *.txt
  !----------------------------------------
    use fpcomm
    use fowcomm
    use fpwrite

    implicit none
    double precision,dimension(nrmax,nsamax) :: Dr, Sr_part
    !double precision,dimension(nrmax,nsamax) :: Drw,Drwav,Drweff,Drweffav
    double precision,dimension(nrmax,nsamax) :: Dnba,Dnpla
    double precision,dimension(nrmax,nsamax) :: Sr_ba,Sr_pla
    double precision,dimension(nrmax,nsamax) :: Sr_ch
    double precision,dimension(nrmax,nsamax) :: heatfow!, heatrw, heatrwav
    double precision,dimension(nrmax,nsamax) :: Sr_hh, Dr_hh, heat_hh, chi_hh
    double precision,dimension(nrmax,nsamax) :: chi_a
    double precision,dimension(nrmax,nsamax) :: hfowout_r, hfowout_p
    double precision,dimension(nrmax,nsamax) :: hfowout_t, hfowout_f
    double precision,dimension(nrmax,nsamax) :: chi_Dp,chi_Dt
    double precision,dimension(nrmax,nsamax) :: chi_Dr, chi_Fr
    double precision,dimension(nrmax,nsamax) :: chi_ch
    double precision,dimension(nrmax,nsamax) :: heatba,heatpla
    double precision,dimension(nrmax,nsamax) :: chi_neo_ba,chi_neo_pla
    double precision,dimension(nrmax,nsamax) :: cyclo_rho, dTadr
    double precision,dimension(nrmax,nsamax) :: Ta, spVI, spVI2,sumdf, sTI
    double precision,dimension(nrmax,nsamax) :: nu_eff
    double precision,dimension(nrmax,nsamax) :: pla_parm, inv_asp
    double precision,dimension(nrmax) :: eps_t,tau_ele,tau_i
    integer nth, np, nr, nsa, nsb, mode(3), nstpmax, nstp

    ! do nr = 1, nrmax
    ! write(*,*)"Bin:",Bin(nr)
    ! end do

    ! do nsa = 1, nsamax
      do nr = 1, nrmax
        ! do np = 1, npmax
          ! do nth = 1, nthmax
            ! nstpmax = orbit_m(nth,np,nr,nsa)%nstp_max
            ! do nstp = 1, nstpmax
          ! write(*,*)"thetap_pnc:", NO_PINCH_ORBIT!nsa,nth_pnc_tg(np,nr,nsa)
          ! if (IBCflux_ratio(np,nr,nsa)==0.d0)write(*,*)"nsa,IBC:",nsa, IBCflux_ratio(np,nr,nsa)
          ! write(*,*)"nsa,IBC:",nsa, IBCflux_ratio(np,nr,nsa)
          ! write(*,*)"nsa,nr_pinch:",nsa, theta_pnc(np,nr,nsa),nr_rhom_pinch(np,nr,nsa)
          ! write(*,*)"nsa,IBC:",nsa, thetam(nth,np,nr,nsa)
          ! write(*,*)"nsa,stg:",nsa, theta_co_stg(np,nr,nsa), theta_cnt_stg(np,nr,nsa)
          ! write(*,*)"nsa,F:",nsa, orbit_m(nth,np,nr,nsa)%F(nstp)
          write(*,*)"nr,psi:",nr, psim(nr), Bout(nr)
          ! write(*,*)"pinc:",nth_pnc(nsa),nth_pnc_tg(np,nr,nsa), nth_pnc_pg(np,nr,nsa), nth_pnc_rg(np,nr,nsa)
          ! if(theta_pnc_pg(np,nr,nsa)/= NO_PINCH_ORBIT)write(*,*)"pinc:",nth_pnc_pg(np,nr,nsa), theta_cnt_stg_pg(np,nr,nsa)
          ! write(*,*)"pinc:",nth_pnc(nsa),theta_pnc_rg(np,nr,nsa)
          ! write(*,*)"pism_rg:",nr,psim_rg(nr),psim(nr)
          ! write(*,*)"pism_rg:",thetam_rg(nth,np,nr,nsa)
          ! write(*,*)"pinc:",nsa,nr_cgt_stg_pg(np,nr,nsa),nr_cgt_stg_tg(np,nr,nsa),nr_cgt_stg_rg(np,nr,nsa)
          ! write(*,*)"pinc:",nsa,nr_stg_inv_pg(np,nr,nsa),nr_stg_inv_tg(np,nr,nsa),nr_stg_inv_rg(np,nr,nsa)
          ! if(nr == 2)write(*,*)"pm_stg:",nsa,pm_stg_pg(nth,np,nr,nsa),pm_stg_tg(nth,np,nr,nsa),pm_stg_rg(nth,np,nr,nsa)
          ! write(*,*)"pm_stg:",shape(thetam)
          ! write(*,*)"pm_stg:",PTFP()：w
          ! write(*,*)"nsa,ptfp,amfp:",nsa,ptfp0(nsa),AMFP(nsa)
          ! write(*,*)"pm_stg:",thetam_pg(nth,np,nr,nsa)
          ! write(*,*)"pm_stg:",Fpsi(nr), (RR-RA*rm(nr))*Bin(nr), &
          ! (RR+RA*rm(nr))*Bout(nr)
          ! B_p = rm(nr)*RA*BB/(safety_factor(nr)*RR)
          ! if(pm_stg_pg(nth,np,nr,nsa)>1.d0)write(*,*)"pm_stg:",nsa,pm_stg_pg(nth,np,nr,nsa)
            ! end do
          ! end do
        ! end do
      end do
    ! end do
          ! do nth = 1, nthmax
          ! write(*,*)"dble:",dble(nth)
          ! end do

    !**** make data folder
    write(*,*)"------------------- check neoclass  start---------------------"
    call system('mkdir -p dat')

    !**** calculation of physical values
    call integral_ParticleDif(Sr_part,Dr)
    write(*,*)"------------------- check neoclass 1 particle ---------------------"
    ! call Particleneoclass(Sr_ba,Sr_pla,Dnba,Dnpla) ![2022/3/7] editted by anzai
    !**
    call integral_Heatdiff(heatfow,chi_a) ![2022/6/6] editted by anzai
    write(*,*)"------------------- check neoclass 2 heat ---------------------"
    ! call Rosenbluth_Hazeltine_Hinton_Neoclass(heatba,heatpla,chi_neo_ba,chi_neo_pla) ![2022/3/14] editted by anzai
    call Rosenbluth_Hazeltine_Hinton_Neoclass(heatba,heatpla,chi_neo_ba,chi_neo_pla,Sr_ba,Sr_pla, Dnba, Dnpla) ![2022/9/5] editted by anzai
    write(*,*)"------------------- check neoclass 3 RHH ---------------------"
    call Chang_Hinton_neoclass(Sr_ch, chi_ch) ![2022/10/26] editted by anzai
    write(*,*)"------------------- check neoclass 4 CH---------------------"
    call Hinton_Hazeltine_neoclass(Sr_hh, Dr_hh, heat_hh, chi_hh) ![2022/10/27] added by anzai
    !** For Factor check
    ! call check_factorneo(eps_t, cyclo_rho, dTadr,Ta,tau_ele,tau_i, spVI,spVI2,sTI) ![2022/2/27] editted by anzai
    call ch_intHD(hfowout_r,hfowout_p,hfowout_t,hfowout_f,chi_Dr,chi_Dp,chi_Dt,chi_Fr) ![2022/6/6]
    write(*,*)"------------------- check neoclass 5 HH---------------------"

    !**** write txt file
    !** For Particle Neoclass
    call fptxt2D(Sr_part,"dat/Sr_part.txt")
    call fptxt2D(Dr,"dat/Dr.txt")
    call fptxt2D(Sr_ba, "dat/Sr_ba.txt") ![2022/7/25] edited by anzai
    call fptxt2D(Sr_pla, "dat/Sr_pla.txt") ![2022/7/25] editted by anzai
    call fptxt2D(Dnba,"dat/Dnba.txt")
    call fptxt2D(Dnpla,"dat/Dnpla.txt")
    ! call fptxt2D(Drw,"dat/Drw.txt")
    ! call fptxt2D(Drwav,"dat/Drwav.txt")
    ! call fptxt2D(Drweff,"dat/Drweff.txt") ![2022/3/2] editted by anzai
    ! call fptxt2D(Drweffav,"dat/Drweffav.txt") ![2022/3/2] editted by anzai
    !** For Heat Neoclass
    call fptxt2D(heatfow, "dat/heatfow.txt")![2022/2/19] editted by anzai
    call fptxt2D(chi_a, "dat/chi_a.txt")![2022/2/19] editted by anzai
    call fptxt2D(heatba, "dat/heatba.txt")![2022/3/4] editted by anzai
    call fptxt2D(heatpla, "dat/heatpla.txt")![2022/3/4] editted by anzai
    call fptxt2D(chi_neo_ba, "dat/chi_neo_ba.txt")![2022/2/19] editted by anzai
    call fptxt2D(chi_neo_pla, "dat/chi_neo_pla.txt")![2022/2/19] editted by anzai
    call fptxt2D(Sr_ch, "dat/Sr_ch.txt")![2022/10/27] editted by anzai
    call fptxt2D(chi_ch, "dat/chi_ch.txt")![2022/10/27] editted by anzai
    call fptxt2D(Sr_hh, "dat/Sr_hh.txt")![2022/10/27] editted by anzai
    call fptxt2D(heat_hh, "dat/heat_hh.txt")![2022/10/27] editted by anzai
    call fptxt2D(Dr_hh, "dat/Dr_hh.txt")![2022/10/27] editted by anzai
    call fptxt2D(chi_hh, "dat/chi_hh.txt")![2022/10/27] editted by anzai
    ! call fptxt2D(heatrw, "dat/heatrw.txt")![2022/2/23] editted by anzai
    ! call fptxt2D(heatrwav, "dat/heatrwav.txt")![2022/2/23] editted by anzai

    !***** For factor check
    ! call fptxt1D(eps_t, "dat/eps_t.txt") ![2022/2/27] editted by anzai
    ! call fptxt2D(cyclo_rho, "dat/cyclo_rho.txt") ![2022/2/27] editted by anzai
    ! call fptxt2D(spVI, "dat/spVI.txt") ![2022/3/4] editted by anzai
    ! call fptxt2D(spVI2, "dat/spVI2.txt") ![2022/3/4] editted by anzai
    ! call fptxt2D(sumdf, "dat/sumdf.txt") ![2022/3/4] editted by anzai
    ! call fptxt2D(dTadr, "dat/dTadr.txt") ![2022/2/27] editted by anzai
    ! call fptxt2D(Ta, "dat/Ta.txt") ![2022/2/27] editted by anzai
    ! call fptxt2D(sTI, "dat/sTI.txt") ![2022/6/10] editted by anzai
    ! call fptxt1D(tau_ele, "dat/tau_ele.txt") ![2022/2/27] editted by anzai
    ! call fptxt1D(tau_i, "dat/tau_i.txt") ![2022/2/27] editted by anzai
    ! call fptxt2D(nu_eff, "dat/nu_eff.txt") ![2022/3/24] editted by anzai
    ! call fptxt2D(hfowout_r, "dat/hfow_r.txt")![2022/6/6] editted by anzai
    ! call fptxt2D(hfowout_p, "dat/hfow_p.txt")![2022/6/6] editted by anzai
    ! call fptxt2D(hfowout_t, "dat/hfow_t.txt")![2022/6/6] editted by anzai
    ! call fptxt2D(hfowout_f, "dat/hfow_f.txt")![2022/6/6] editted by anzai
    call fptxt2D(chi_Dr, "dat/chi_Dr.txt")![2022/2/19] editted by anzai
    call fptxt2D(chi_Dp, "dat/chi_Dp.txt")![2022/2/19] editted by anzai
    call fptxt2D(chi_Dt, "dat/chi_Dt.txt")![2022/2/19] editted by anzai
    call fptxt2D(chi_Fr, "dat/chi_Fr.txt")![2022/2/19] editted by anzai
    call fptxt2D(pla_parm, "dat/pla_parm.txt")![2022/2/19] editted by anzai
    call fptxt2D(inv_asp, "dat/inv_asp.txt")![2022/2/19] editted by anzai

  end subroutine output_neoclass

!!!****==================================
  ! Modules of collisional flecuency
!!!****==================================
  subroutine cal_nu_ei(nu)
  !--------------------------------------
  ![2022/2/8] Modified by anzai
  !--------------------------------------
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out) :: nu
    double precision fact, Te, tau
    integer nsa, nr

    fact = 12.d0*pi**1.5d0/sqrt(2.d0)*eps0**2/aee**2*SQRT(AME)
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Te = rwsl(nr,1)*1.d6/( 1.5d0*rnsl(nr,1)*1.d20 )
        tau = fact*Te**1.5d0/( lnlam(nr,2,1)*rnsl(nr,2)*1.d20*aefp(2)**2 )
        nu(nr,nsa) = 1.d0/tau !** [J]
      end do
    end do

  end subroutine

  subroutine nu_effective(nu_eff)
  !------------------------------------------------------------------
  ! Calculation of effective collisional flecency[keV] in "Tokamaks"
  ! Formulae from "Tokamaks fourth edition"
  !------------------------------------------------------------------
    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out) :: nu_eff
    double precision,dimension(nrmax,nsamax) :: Ta
    double precision :: fact, tau_e, tau_i
    integer nr, nsa

    !**** initialization
    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)

    !**** Temperature and flecency make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/( 1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
        !**** Ta(temperature)[keV] /RKEV =*1/( AEE*1.D3)
        if (nsa == 1) then
          nu_eff(nr,nsa) = (1.D0/6.4D14)*rnsl(nr,nsa)*1.D20/(Ta(nr,nsa)**1.5D0)
          nu_eff(nr,nsa) = nu_eff(nr, nsa)/(rm(nr)*RA/RR)
        else
          nu_eff(nr,nsa) = (1.D0/6.6D17)*rnsl(nr,nsa) &
                         * 1.D20/(Ta(nr,nsa)**1.5D0) &
                         * sqrt(AMP/AMFP(nsa)) &
                         * lnlam(nr,nsa,nsa)
                         !** AMP : proton mass
          nu_eff(nr,nsa) = nu_eff(nr, nsa)/(rm(nr)*RA/RR)
        end if
        !**** Effective collisional flecuency nu
      end do
    end do
  end subroutine nu_effective

  subroutine nu_deflection(nud)
  !----------------------------------
  ! modified by anzai [2022/3/2]
  !----------------------------------
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
          !** particle verocity [m/s]

          va = vc * SQRT( 1.d0-pv**(-2) )
          !** vc is speed of light[m/s]

          do nsb = 1, nsbmax
            nssb = ns_nsb(nsa)
            Tb = rwsl(nr,nsb)*1.d6/( 1.5d0*rnsl(nr,nsb)*1.d20 )
            !** Tb[J]
            vthb = SQRT( 2.d0*Tb/amfp(nsb) )

            numerator  = rnsl(nr,nsb)*1.d20 &
                 *AEFP(nsa)**2*AEFP(nsb)**2 * lnlam(nr,nsb,nsa)
            denominator= SQRT( amfp(nsb) )*Tb**1.5d0*SQRT(amfp(nsa)/amfp(nsb))
            nuBra = factBra*numerator/denominator

            x = va/vthb
            call gosakn(x,fx,gx)

            nudb = fact*nuBra*(fx-gx)/x**3
            nud(np,nr,nsa) = nud(np,nr,nsa)+nudb

          end do
        end do
      end do
    end do

  end subroutine nu_deflection

!!!****======================
  ! Other needed subroutine
!!!****======================
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

    if(abs(x) .gt. 1.e-10) then
    !** .gt. is >,1.e0 is single precision

         gx=(fx-x*fx1)/(2.d0*x**2)
    end if

    return

  end subroutine gosakn

  ! subroutine bounce_average_for_Drw(Dout, Din)
 
  !   use fpcomm
  !   use fowcomm
  !   use fowcoef
 
  !   implicit none
  !   double precision,dimension(nthmax,npmax,nrmax,nsamax),intent(out) :: Dout
  !   double precision,dimension(nthmax,npmax,nrmax,nsamax),intent(in)  :: Din
  !   double precision,dimension(3,3,max_stp) :: dIdul
  !   double precision,allocatable :: U(:,:,:,:,:,:,:,:), Drwl(:,:,:,:,:)
  !   double precision dt, taup, cpitch_ob, psip_ob, thetap_ob, Drw_ob
  !   type(orbit) ob
  !   integer nth, np, nr, nsa, nstp, nstpmax, nthp, mode(3)
 
  !   allocate(U(4,4,4,nthmax,nrmax,nthpmax,npmax,nsamax))
  !   allocate(Drwl(nthmax,npmax,nrmax,nthpmax,nsamax))

  !   mode = [0,0,0]
 
  !   do nsa = 1, nsamax
  !     do nthp = 1, nthpmax
  !       do nr = 1, nrmax
  !         do np = 1, npmax
  !           do nth = 1, nthmax
  !             Drwl(nth,np,nr,nthp,nsa) = Din(nth,np,nr,nsa)
  !           end do
  !         end do
  !       end do
  !     end do
  !   end do

  !   do nsa = 1, nsamax
  !     call make_U_Dxy(U, Drwl, 'm', nsa)
  !   end do
 
  !   do nsa = 1, nsamax
  !     do nr = 1, nrmax
  !       do np = 1, npmax
  !         do nth = 1, nthmax
  !           ob = orbit_m(nth,np,nr,nsa)
  !           nstpmax = ob%nstp_max
  !           taup = ob%time(nstpmax)
 
  !           call transformation_matrix(dIdul, ob, nth, np, nr, nsa, mode)
 
  !           Dout(nth,np,nr,nsa) = 0.d0

  !           do nstp = 2, nstpmax
  !             dt = ob%time(nstp)-ob%time(nstp-1)
  !             cpitch_ob = ob%costh(nstp)
  !             psip_ob   = ob%psip(nstp)
  !             thetap_ob = ob%thetap(nstp)
  !             call interpolate_D_unlessZero(Drw_ob, U(:,:,:,:,:,:,np,nsa), &
  !                   1.d0, cpitch_ob, psip_ob, thetap_ob)
 
  !             Dout(nth,np,nr,nsa) = Dout(nth,np,nr,nsa)&
  !                                 + Drw_ob*dIdul(3,3,nstp)**2*dt
  !           end do
  !           Dout(nth,np,nr,nsa) = Dout(nth,np,nr,nsa)/taup
  !         end do
  !       end do
  !     end do
  !   end do
  ! end subroutine

!===================================================
! Subroutines for particle diffusion coef
!===================================================
  subroutine integral_ParticleDif(Sr_part,Dr)
  !---------------------------------------------
  ! calculation of FOW diff coef
  ! 1/J * dfdx_i * dV(volume)
  !---------------------------------------------
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out) :: Sr_part, Dr
    double precision,dimension(nrmax,nsamax) :: Na, dNadr
    double precision,dimension(3,nthmax,npmax,nrmax,nsamax) ::moment_vector ! by anzai
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom,dfdp, divDdfdp
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: temp_dfdthm,temp_dfdp
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdthm, divDdfdthm
    double precision,dimension(nrmax,nsamax) :: add_sr
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: fnsp_l
    real(rkind) :: sum_, sum2_
    double precision Drr_, dVI, Frr_, JIp, JIm
    integer nth,np,nr,nsa,ns, nr_temp

    !**** initialization ****
    Sr_part(:,:) = 0.d0
    Dr(:,:) = 0.d0
    divDdfdp(:,:,:,:) = 0.d0
    divDdfdthm(:,:,:,:) = 0.d0
    add_sr(:,:) = 0.d0
    !****temperature make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Na(nr,nsa) = rnsl(nr,nsa)*1.d20
      end do
    end do
    !****first order derivation
    do nsa = 1, nsamax
      call first_order_derivative(dNadr(:,nsa), Na(:,nsa), rm)
    end do
    ! write(*,*)"------------------- check 1 ---------------------"

    !********** for new dfdrhom module [2022/3/4]
    fnsp_l(:,:,:,:)=0.d0
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr= 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit
          do nth = 1, nthmax
          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa)
          end do
        end do
      end do
    end do

    !**** first order derivative
    do nsa = 1, nsamax
      do np = 1, npmax
        do nth = 1, nthmax
          call first_order_derivative(dfdrhom(nth,np,:,nsa), &
                 fnsp_l(nth,np,:,nsa)*1.d20, rm)
        end do
      end do
    end do
    ! write(*,*)"------------------- check 2 ---------------------"
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          call first_order_derivative(dfdthm(:,np,nr,nsa), &
                 fnsp_l(:,np,nr,nsa)*1.d20, thetam(:,np,nr,nsa))
        end do
      end do
    end do
    ! write(*,*)"------------------- check 3 ---------------------"
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do nth = 1, nthmax
          call first_order_derivative(dfdp(nth,:,nr,nsa), &
                 fnsp_l(nth,:,nr,nsa)*1.d20, pm(:,nsa))
        end do
      end do
    end do
    ! write(*,*)"------------------- check 4 ---------------------"
    ! !**** Making additional source term
    ! do nsa = 1, nsamax
    !   do nr = 1, nrmax
    !     do np = 1, npmax
    !       do nth = 1, nthmax
    !         temp_dfdp(nth,np,nr,nsa) = 0.d0 &
    !                                 !  + Dppfow(nth,np,nr,nsa) * dfdp(nth,np,nr,nsa) &
    !                                 !  + Dptfow(nth,np,nr,nsa) * dfdthm(nth,np,nr,nsa) &
    !                                  + Dprfow(nth,np,nr,nsa) * dfdrhom(nth,np,nr,nsa)! &
    !                                 !  - Fppfow(nth,np,nr,nsa) * fnsp_l(nth,np,nr,nsa)*1.d20
    !         temp_dfdthm(nth,np,nr,nsa) = 0.d0 &
    !                                   !  + Dtpfow(nth,np,nr,nsa) * dfdp(nth,np,nr,nsa) &
    !                                   !  + Dttfow(nth,np,nr,nsa) * dfdthm(nth,np,nr,nsa) &
    !                                    + Dtrfow(nth,np,nr,nsa) * dfdrhom(nth,np,nr,nsa)! &
    !                                   !  - Fthfow(nth,np,nr,nsa) * fnsp_l(nth,np,nr,nsa)*1.d20
    !       end do
    !     end do
    !   end do
    ! end do
    ! do nsa = 1, NSAMAX
    !   do nr = 1, NRMAX
    !     do nth = 1, NTHMAX
    !       call first_order_derivative(divDdfdp(nth,:,nr,nsa), &
    !               temp_dfdp(nth,:,nr,nsa), pm(:,nsa))
    !     end do
    !   end do
    ! end do
    ! do nsa = 1, nsamax
    !   do nr = 1, nrmax
    !     do np = 1, npmax
    !       call first_order_derivative(divDdfdthm(:,np,nr,nsa), &
    !               temp_dfdthm(:,np,nr,nsa), thetam(:,np,nr,nsa))
    !     end do
    !   end do
    ! end do 
    ! do nr_temp = 1, nrmax !** for place integral 
    !   do nsa = 1, nsamax
    !     do nr = 1, nr_temp
    !     ! do nr = 1, nrmax
    !       do np = 1, npmax
    !         do nth = 1, nthmax
    !           add_sr(nr_temp,nsa) = -( 0.d0&
    !           ! add_sr(nr,nsa) = -( 0.d0&
    !                                      + divDdfdp(nth,np,nr,nsa)/JIR(nth,np,nr,nsa) &
    !                                     !  + divDdfdthm(nth,np,nr,nsa)/JIR(nth,np,nr,nsa)&
    !                                     ) &
    !                                      * delp(nsa) * delthm(nth,np,nr,nsa) * delr &
    !                                      * JIR(nth,np,nr,nsa)!*1.d20
    !         end do !** np
    !       end do !** nth
    !     end do !** nr
    !   end do !** nsa
    ! end do !** nr_temp

    !**** end of new dfdrhom module
    ! !**** making radial mean vector
    ! call make_moment_vector(moment_vector)
    ! do nsa = 1, nsamax
    !   ns = ns_nsa(nsa)
    !   do nr = 1, nrmax
    !     do np = 1, npmax
    !       if ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
    !       do nth = 1, nthmax
    !        !**** making mean vector
    !         Sr_part(nr,nsa) = moment_vector(3,nth,np,nr,nsa)*ptfp0(nsa)/AMFP(nsa) &
    !                             * fnsp_l(nth,np,nr,nsa)*1.d20 &
    !                             * JIR(nth,np,nr,nsa) &
    !                             / JIR(nth,np,nr,nsa) &
    !                             * delp(ns) * delthm(nth,np,nr,nsa)
    !         sum_ = sum_ + fnsp_l(nth,np,nr,nsa)*1.d20 &
    !             * ( &
    !             ! moment_vector(1,nth,np,nr,nsa)**2.0 &
    !             ! + moment_vector(2,nth,np,nr,nsa)**2.0 &
    !             + moment_vector(3,nth,np,nr,nsa)**2.0 &
    !             ! + moment_vector(2,nth,np,nr,nsa)*moment_vector(3,nth,np,nr,nsa)*2.0 &
    !             ! + moment_vector(1,nth,np,nr,nsa)*moment_vector(2,nth,np,nr,nsa)*2.0 &
    !             ! + moment_vector(1,nth,np,nr,nsa)*moment_vector(3,nth,np,nr,nsa)*2.0 &
    !             ) /(pm(np,nsa)**2.0) &
    !                             * JIR(nth,np,nr,nsa) &
    !                             *delp(ns) * delthm(nth,np,nr,nsa)
    !         sum2_ = sum2_ + fnsp_l(nth,np,nr,nsa)*1.d20 &
    !                             * JIR(nth,np,nr,nsa) &
    !                             *delp(ns) * delthm(nth,np,nr,nsa)
    !       end do !** nth
    !     end do !** np
    !     ! radial_mean_moment_vector(nr,nsa) = radial_mean_moment_vector(nr,nsa)/sum_
    !     ! mean_moment_vector(nr,nsa) = mean_moment_vector(nr,nsa)/sum_
    !     ! write(*,*)radial_mean_moment_vector(nr,nsa),mean_moment_vector(nr,nsa)
    !     Dr(nr,nsa) = -Sr_part(nr,nsa)/dNadr(nr,nsa)
    !     write(*,*)'sum_:', sum_/sum2_
    !     sum_ = 0.0
    !     sum2_ = 0.0
    !   end do !** nr
    ! end do !** nsa

    !**** Source calculation
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit
          do nth = 1, nthmax

            if (nr /= 1 .and. nr /= nrmax) then
              JIp = (JIR(nth,np,nr+1,nsa)+JIR(nth,np,nr,nsa))*0.5d0
              JIm = (JIR(nth,np,nr,nsa)+JIR(nth,np,nr-1,nsa))*0.5d0
              Drr_ = ((Drrfow(nth,np,nr+1,nsa)+Drrfow(nth,np,nr,nsa))/JIp*0.5d0 &
                   + (Drrfow(nth,np,nr,nsa)+Drrfow(nth,np,nr-1,nsa))/JIm*0.5d0)*0.5d0
              Frr_ = ((Frrfow(nth,np,nr+1,nsa)+Frrfow(nth,np,nr,nsa))/JIp*0.5d0 &
                   + (Frrfow(nth,np,nr,nsa)+Frrfow(nth,np,nr-1,nsa))/JIm*0.5d0)*0.5d0
            else if (nr == 1) then
              JIp = (JIR(nth,np,nr+1,nsa)+JIR(nth,np,nr,nsa))*0.5d0
              JIm = (JIR(nth,np,nr,nsa)+JIR(nth,np,nr,nsa))*0.5d0
              Drr_ = ((Drrfow(nth,np,nr+1,nsa)+Drrfow(nth,np,nr,nsa))/JIp*0.5d0 &
                   + (Drrfow(nth,np,nr,nsa)+Drrfow(nth,np,nr,nsa))/JIm*0.5d0)*0.5d0
              Frr_ = ((Frrfow(nth,np,nr+1,nsa)+Frrfow(nth,np,nr,nsa))/JIp*0.5d0 &
                   + (Frrfow(nth,np,nr,nsa)+Frrfow(nth,np,nr,nsa))/JIm*0.5d0)*0.5d0
            else
              JIp = (JIR(nth,np,nr,nsa)+JIR(nth,np,nr,nsa))*0.5d0
              JIm = (JIR(nth,np,nr-1,nsa)+JIR(nth,np,nr,nsa))*0.5d0
              Drr_ = ((Drrfow(nth,np,nr+1,nsa)+Drrfow(nth,np,nr,nsa))/JIp*0.5d0 &
                   + (Drrfow(nth,np,nr-1,nsa)+Drrfow(nth,np,nr,nsa))/JIm*0.5d0)*0.5d0
              Frr_ = ((Frrfow(nth,np,nr+1,nsa)+Frrfow(nth,np,nr,nsa))/JIp*0.5d0 &
                   + (Frrfow(nth,np,nr-1,nsa)+Frrfow(nth,np,nr,nsa))/JIm*0.5d0)*0.5d0
            end if

            Sr_part(nr,nsa) = Sr_part(nr,nsa) &
                            - 1.d0&
                            !* (Drr_j(nth,np,nr,nsa)&
                            * ( 0.d0 &
                            ! + Drrfow(nth,np,nr,nsa)&
                            + Drr_ &
                            * dfdrhom(nth,np,nr,nsa)&!*1.d20 &
                            ! + Drpfow(nth,np,nr,nsa) &
                            ! * dfdp(nth,np,nr,nsa)&!*1.d20 &
                            ! + Drtfow(nth,np,nr,nsa) &
                            ! * dfdthm(nth,np,nr,nsa)&!*1.d20&
                            ) &
                            ! * JIR(nth,np,nr,nsa) &
                            ! / JIR(nth,np,nr,nsa) &
                            * delp(ns) * delthm(nth,np,nr,nsa) &

                            - 1.d0&
                            !* Frr_j(nth,np,nr,nsa) &
                            ! * Frrfow(nth,np,nr,nsa) &
                            * fnsp_l(nth,np,nr,nsa)*1.d20 &
                          !  * JIR(nth,np,nr,nsa) &
                            * Frr_ &
                            ! / JIR(nth,np,nr,nsa) &
                            * delp(ns) * delthm(nth,np,nr,nsa) &
                            + 0.d0
    ! write(*,*)"------------------- check tt ---------------------"

                            ! + add_sr(nr,nsa)
                            ! write(*,*)"add:",add_sr(nr,nsa)
          end do
        end do
    ! write(*,*)"------------------- check error ---------------------"
        Sr_part(nr,nsa) = Sr_part(nr,nsa) !- add_sr(nr,nsa)
                            ! write(*,*)"add:",add_sr(nr,nsa)
    ! write(*,*)"------------------- check error 2 ---------------------"
        Dr(nr,nsa) = -Sr_part(nr,nsa)/dNadr(nr,nsa)
    ! write(*,*)"------------------- check error 3 ---------------------"
      end do
    end do
    ! write(*,*)"------------------- check 5 ---------------------"

  end subroutine integral_ParticleDif

!   subroutine D_random_walk_baverage(Drwav, Drw, Drweff, Drweffav)
!   !-----------------------------------------------------
!   ! For particle diffusion in banana region
!   ! modified and added by anzai [2022/3/2]
!   !-----------------------------------------------------
!     use fpcomm
!     use fowcomm
!     implicit none
!     double precision,dimension(nrmax,nsamax),intent(out)  :: Drw, Drwav
!     double precision,dimension(nrmax,nsamax),intent(out)  :: Drweff, Drweffav
!     double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom
!     double precision,dimension(nthmax,npmax,nrmax,nsamax) :: fnsp_l
!     double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
!     double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drweff_l, Drweff_ba
!     double precision,dimension(npmax,nrmax,nsamax) :: nu
!     double precision,dimension(nrmax,nsamax) :: Ta
!     double precision,dimension(nrmax,nsamax) :: nu_eff
!     double precision step_len, step_time, step_time_eff, eps_t
!     double precision rho_theta, Baxis, pv, B_theta, rho
!     double precision dVI, sumVI
!     integer nth, np, nr, nsa, ns

! !     !**** initialization ****
!     Baxis = Bin_rg(1)
!     Drw(:,:) = 0.d0
!     Drwav(:,:) = 0.d0
!     sumVI = 0.d0
!     call nu_deflection(nu)
!
! !********** for new dfdrhom module [2022/3/4]
!     fnsp_l(:,:,:,:)=0.d0
!     do nsa = 1, nsamax
!       ns = ns_nsa(nsa)
!       do nr= 1, nrmax
!         do np = 1, npmax
!           if ( pm(np,ns) > fact_bulk ) exit
!           do nth = 1, nthmax
!           fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa)
!           end do
!         end do
!       end do
!     end do
!
!     !**** first order derivative
!     do nsa = 1, nsamax
!       do np = 1, npmax
!         do nth = 1, nthmax
!           call first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp_l(nth,np,:,nsa), rm)
!         end do
!       end do
!     end do
! !!**** end of new dfdrhom module
!
!     !************************by anzai[2022/2/15]
!     do nsa = 1, nsamax
!       do nr = 1, nrmax
!         Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
!         !** Ta(temperature)[keV]
!       end do
!     end do
!
!     call nu_effective(nu_eff)
!
!     do nsa = 1, nsamax
!       do nr = 1, nrmax
!         eps_t = rm(nr)*RA/RR
!         do np = 1, npmax
!
!           step_time = 1.d0/nu(np, nr, nsa)
!           step_time_eff = 1.d0/nu_eff(nr, nsa)
!           rho = pm(np,nsa)*ptfp0(nsa)/(aefp(nsa)*Baxis)
!           !** [keV]verocity's length
!
!           do nth = 1, nthmax
!             step_len = rho * safety_factor(nr)/sqrt(eps_t)
!             Drwlocal(nth,np,nr,nsa) = step_len**2/step_time*sqrt(eps_t)
!             Drweff_l(nth,np,nr,nsa) = step_len**2/step_time_eff*sqrt(eps_t)
!           end do
!         end do
!       end do
!     end do
!
!     call bounce_average_for_Drw(Drwba, Drwlocal)
!     call bounce_average_for_Drw(Drweff_ba, Drweff_l)
!
!     do nsa = 1, nsamax
!       do nr = 1, nrmax
!         do np = 1, npmax
!           do nth = 1, nthmax
!             dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)
!
!             Drw(nr,nsa) = Drw(nr,nsa) - Drwlocal(nth,np,nr,nsa) &
!                         * dfdrhom(nth,np,nr,nsa)*dVI
!             Drwav(nr,nsa) = Drwav(nr,nsa) - Drwba(nth,np,nr,nsa) &
!                           * dfdrhom(nth,np,nr,nsa)*dVI
!             Drweff(nr,nsa) = Drweff(nr,nsa) - Drweff_l(nth,np,nr,nsa) &
!                            * dfdrhom(nth,np,nr,nsa)*dVI
!             Drweffav(nr,nsa) = Drweffav(nr,nsa) - Drweff_ba(nth,np,nr,nsa) &
!                              * dfdrhom(nth,np,nr,nsa)*dVI
!             sumVI = sumVI - dfdrhom(nth,np,nr,nsa)*dVI
!           end do
!         end do
!         Drw(nr,nsa) = Drw(nr,nsa)/sumVI
!         Drwav(nr,nsa) = Drwav(nr,nsa)/sumVI
!         Drweff(nr,nsa) = Drweff(nr,nsa)/sumVI
!         Drweffav(nr,nsa) = Drweffav(nr,nsa)/sumVI
!       end do
!     end do

  ! end subroutine D_random_walk_baverage

  ! subroutine Particleneoclass(Sr_ba,Sr_pla,Dnba,Dnpla)
    !------------------------------------------------------------
    ! Subroutine for neoclassical particle flux in banana region
    ! Calculate particle flux and diffusion coefficient
    ! Formulae from "Tokamaks fourth edition"
    !------------------------------------------------------------

  !   use fpcomm
  !   use fowcomm

  !   implicit none
  !   double precision,dimension(nrmax,nsamax),intent(out)  :: Dnba, Dnpla
  !   double precision,dimension(nrmax,nsamax),intent(out)  :: Sr_ba, Sr_pla
  !       double precision,dimension(nrmax,nsamax)  :: Ta,dTadr,dndr
  !   double precision fact,fact_ba,fact_pla,tau_ele,rho_e,eps_t,Baxis,B_p
  !   integer nr, nsa

  !   fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
  !   Baxis = Bin_rg(1) ! approximation on B by Baxis

  !   !****temperature make
  !   do nsa = 1, nsamax
  !     do nr = 1, nrmax
  !       Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)
  !       !Ta(temperature)[J]
  !     end do
  !   end do

  !   !****first order derivation
  !   do nsa = 1, nsamax
  !     call first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
  !     call first_order_derivative(dndr(:,nsa), rnsl(:, nsa), rm)
  !   end do

  !   !****calculate flux and diffusion coefficient
  !   do nsa = 1, nsamax
  !     do nr = 1, nrmax
  !       tau_ele = fact*sqrt(AMFP(1))*((Ta(nr,1))**1.5d0)/ &
  !            (rnsl(nr,2)*1.d20*AEFP(2)**2*lnlam(nr,2,1)*AEE**2)
  !       rho_e = sqrt(2*Ta(nr,1)/AMFP(1))*AMFP(1)/(AEE*Baxis)
  !       eps_t = rm(nr)*RA/RR
  !       B_p = rm(nr)*RA*BB/(safety_factor(nr)*RR)
  !       fact_ba = (safety_factor(nr)**2)*1.d20*rnsl(nr,nsa) &
  !            *(rho_e**2)/((eps_t**1.5)*tau_ele)*RKEV ! modified by anzai[*RKEV]
  !       fact_pla = -sqrt(pi)*(eps_t**2)*Ta(nr,1)*RKEV &! modified by anzai [*RKEV]
  !            *rho_e*1.d20*rnsl(nr,nsa)/(2*AEE*B_p*rm(nr)*RA)

  !       Sr_ba(nr,nsa) = fact_ba*(-1.22d0*(1+Ta(nr,2)/Ta(nr,1)) &
  !                     * 1.d20 * dndr(nr, nsa) &
  !                     + 4.3d-1*dTadr(nr,1)/Ta(nr,1) &
  !                     + 1.9d-1*dTadr(nr,2)/Ta(nr,1) )
  !       !**** neglect E_para term
  !       Sr_pla(nr,nsa) = fact_pla*((1+Ta(nr,2)/Ta(nr,1)) &
  !                      * 1.d20 * dndr(nr,nsa) &
  !                      / (1.d20*rnsl(nr,nsa)) &
  !                      + 1.5d0*dTadr(nr,1)/Ta(nr,1) + 1.5d0*dTadr(nr,2)/Ta(nr,1) )
  !       !**** neglect B term

  !       Dnba(nr,nsa) = - Sr_ba(nr, nsa)/(1.d20*dndr(nr, nsa))
  !       Dnpla(nr,nsa) = - Sr_pla(nr, nsa)/(1.d20*dndr(nr, nsa))
  !    end do
  !  end do

  ! end subroutine Particleneoclass

!======================================================================
!subroutines for heat fluxes
!======================================================================

  ! subroutine integral_Heatdiff(heatfow_out,chi_a)
  ! !--------------------------------------------
  ! ! Subroutine for FOW heat diffuion coefficient
  ! ! calculate only heat diffusion coefficient
  ! ! 1/J * dfdx_i dV(volume)
  ! !--------------------------------------------

  !   use fpcomm
  !   use fowcomm

  !   implicit none
  !   double precision,dimension(nrmax,nsamax),intent(out) :: heatfow_out
  !   double precision,dimension(nrmax,nsamax),intent(out) :: chi_a
  !   double precision,dimension(nrmax,nsamax) :: heatfow_l,sumVI_l,sumVI
  !   double precision,dimension(nrmax,nsamax) :: Ta, dTadr
  !   double precision,dimension(nthmax,npmax,nrmax,nsamax) ::fnsp_l, dfdrhom
  !   double precision,dimension(nthmax,npmax,nrmax,nsamax) ::dfdthm, dfdp
  !   double precision,dimension(3,nthmax,npmax,nrmax,nsamax) ::moment_vector
  !   double precision,dimension(nrmax,nsamax) :: radial_mean_moment_vector, mean_moment_vector
  !   double precision :: K, PV, sum_
  !   integer nth,np,nr,nsa,ns

  !   !**** initialization ****
  !   heatfow_out(:,:) = 0.d0
  !   sumVI_l(:,:) = 0.d0
  !   fnsp_l(:,:,:,:)=0.d0
  !   sum_ = 0.d0

  !   !****temperature make
  !   do nsa = 1, nsamax
  !     do nr = 1, nrmax
  !       Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
  !       !Ta(temperature)[keV] rwsl[1/m^3*MJ]
  !     end do !** nr
  !   end do !** nsa
  !   !****first order derivation
  !   do nsa = 1, nsamax
  !     call first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
  !   end do !** nsa
  !   !********** for new dfdrhom module [2022/3/4]
  !   do nsa = 1, nsamax
  !     ns = ns_nsa(nsa)
  !     do nr= 1, nrmax
  !       do np = 1, npmax
  !         if ( pm(np,ns) > fact_bulk ) exit
  !         do nth = 1, nthmax
  !         fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa)
  !         end do !** nth
  !       end do !** np
  !     end do !** nr
  !   end do ! ** nsa
  !   !**** first order derivative
  !   do nsa = 1, nsamax
  !     do np = 1, npmax
  !       do nth = 1, nthmax
  !         call first_order_derivative(dfdrhom(nth,np,:,nsa), &
  !                fnsp_l(nth,np,:,nsa), rm)
  !       end do !** nth
  !     end do !** np
  !   end do !**nsa
  !   do nsa = 1, nsamax
  !     do nr = 1, nrmax
  !       do np = 1, npmax
  !         call first_order_derivative(dfdthm(:,np,nr,nsa), &
  !                fnsp_l(:,np,nr,nsa), thetam(:,np,nr,nsa))
  !       end do !** np
  !     end do !** nr
  !   end do !** nsa
  !   do nsa = 1, nsamax
  !     do nr = 1, nrmax
  !       do nth = 1, nthmax
  !         call first_order_derivative(dfdp(nth,:,nr,nsa), &
  !                                     fnsp_l(nth,:,nr,nsa), pm(:,nsa))
  !       end do !** nth
  !     end do !** nr
  !   end do !**nsa
  !   !!**** end of new dfdrhom module
  !   !**** making radial mean vector
  !   call make_moment_vector(moment_vector)
  !   do nsa = 1, nsamax
  !     ns = ns_nsa(nsa)
  !     do nr = 1, nrmax
  !       do np = 1, npmax
  !         if ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
  !         do nth = 1, nthmax
  !          !**** making mean vector
  !           radial_mean_moment_vector(nr,nsa) = moment_vector(3,nth,np,nr,nsa)*ptfp0(nsa) &
  !                               * fnsp_l(nth,np,nr,nsa)*1.d20 &
  !                               * JIR(nth,np,nr,nsa) &
  !                               / JIR(nth,np,nr,nsa) &
  !                               * delp(ns) * delthm(nth,np,nr,nsa)
  !           mean_moment_vector(nr,nsa) = pm(np,nsa)*ptfp0(nsa) &
  !                               * fnsp_l(nth,np,nr,nsa)*1.d20 &
  !                               * JIR(nth,np,nr,nsa) &
  !                               / JIR(nth,np,nr,nsa) &
  !                               * delp(ns) * delthm(nth,np,nr,nsa)
  !           sum_ = sum_ + fnsp_l(nth,np,nr,nsa)*1.d20 &
  !                               * JIR(nth,np,nr,nsa) &
  !                               *delp(ns) * delthm(nth,np,nr,nsa)
  !         end do !** nth
  !       end do !** np
  !       radial_mean_moment_vector(nr,nsa) = radial_mean_moment_vector(nr,nsa)/sum_
  !       mean_moment_vector(nr,nsa) = mean_moment_vector(nr,nsa)/sum_
  !       ! write(*,*)'r_meanv,mean_v',radial_mean_moment_vector(nr,nsa),mean_moment_vector(nr,nsa)
  !     end do !** nr
  !   end do !** nsa
  !   !****Integration over moment(np) and pitch angle(nth)
  !   do nsa = 1, nsamax
  !     ns = ns_nsa(nsa)
  !     do nr = 1, nrmax
  !       do np = 1, npmax
  !         if ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
  !         do nth = 1, nthmax
  !           ! PV = sqrt(1.d0+theta0(nsa)*pm(np,nsa)**2)
  !           ! K = (PV-1.d0)*AMFP(nsa)*vc**2
  !           K = ((pm(np,nsa)-mean_moment_vector(nr,nsa))*ptfp0(nsa))**2/AMFP(nsa)/2.d0 &
  !             / ((moment_vector(3,nth,np,nr,nsa)-radial_mean_moment_vector(nr,nsa))*ptfp0(nsa)/AMFP(nsa))/AEE/1.d3
  !           ! write(*,*)'Kinetic E:',K
  !           !**** for dfdrhom is made of fnsp*dVI
  !           heatfow_out(nr,nsa) = heatfow_out(nr,nsa) &
  !                               !- (pm(np,nsa)*ptfp0(nsa))**2 &
  !                               !/ (2*AMFP(nsa)) &
  !                               + K &
  !                               * fnsp_l(nth,np,nr,nsa)*1.d20 &
  !                               * JIR(nth,np,nr,nsa) &
  !                               / JIR(nth,np,nr,nsa) &
  !                               * delp(ns) * delthm(nth,np,nr,nsa) 
  !         end do !** nth
  !       end do !** np
  !       chi_a(nr,nsa) = -heatfow_out(nr,nsa)/(dTadr(nr,nsa)*rnsl(nr,nsa)*1.d20)
  !     end do !** nr
  !   end do !** nse
  ! end subroutine integral_Heatdiff

  subroutine integral_Heatdiff(heatfow_out,chi_a)
  !--------------------------------------------
  ! Subroutine for FOW heat diffuion coefficient
  ! calculate only heat diffusion coefficient
  ! 1/J * K * dfdx_i dV(volume)
  !--------------------------------------------

    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out) :: heatfow_out
    double precision,dimension(nrmax,nsamax),intent(out) :: chi_a
    double precision,dimension(nrmax,nsamax) :: heatfow_l,sumVI_l,sumVI
    double precision,dimension(nrmax,nsamax) :: Ta, dTadr, add_heat
    double precision,dimension(nthmax,npmax,nrmax,nsamax) ::fnsp_l, dfdrhom
    double precision,dimension(nthmax,npmax,nrmax,nsamax) ::dfdthm, dfdp
    double precision,dimension(3,nthmax,npmax,nrmax,nsamax) ::moment_vector
    double precision,dimension(nrmax,nsamax) :: radial_mean_moment, mean_moment
    double precision :: K, PV, sum_, Drr_, Frr_, Drp_, Drt_, JIp, JIm
    integer nth,np,nr,nsa,ns



    !**** initialization ****
    heatfow_out(:,:) = 0.d0
    sumVI_l(:,:) = 0.d0
    fnsp_l(:,:,:,:)=0.d0
    radial_mean_moment(:,:) = 0.d0
    add_heat(:,:) = 0.d0
    sum_ = 0.d0

    !****temperature make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
        !Ta(temperature)[keV] rwsl[1/m^3*MJ]
      end do
    end do
    !****first order derivation
    do nsa = 1, nsamax
      call first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
    end do
!********** for new dfdrhom module [2022/3/4]
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr= 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit
          do nth = 1, nthmax
          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa)
          end do
        end do
      end do
    end do

    !**** first order derivative
    do nsa = 1, nsamax
      do np = 1, npmax
        do nth = 1, nthmax
          call first_order_derivative(dfdrhom(nth,np,:,nsa), &
                 fnsp_l(nth,np,:,nsa), rm)
        end do
      end do
    end do
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          call first_order_derivative(dfdthm(:,np,nr,nsa), &
                 fnsp_l(:,np,nr,nsa), thetam(:,np,nr,nsa))
        end do
      end do
    end do
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do nth = 1, nthmax
          call first_order_derivative(dfdp(nth,:,nr,nsa), &
                 fnsp_l(nth,:,nr,nsa), pm(:,nsa))
        end do
      end do
    end do
!!**** end of new dfdrhom module
    !**** making radial mean vector
    ! call make_moment_vector(moment_vector)
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
          do nth = 1, nthmax
           !**** making mean vector
            mean_moment(nr,nsa) = pm(np,nsa)*PTFP0(nsa) &
                                * fnsp_l(nth,np,nr,nsa)*1.d20 &
                                * JIR(nth,np,nr,nsa) &
                                * delp(ns) * delthm(nth,np,nr,nsa)
            sum_ = sum_ + fnsp_l(nth,np,nr,nsa)*1.d20 &
                                * JIR(nth,np,nr,nsa) &
                                *delp(ns) * delthm(nth,np,nr,nsa)
          end do
        end do
        mean_moment(nr,nsa) = mean_moment(nr,nsa)/sum_
        sum_ = 0.d0
        ! write(*,*)"mean_moment:",mean_moment(nr,nsa)
      end do
    end do

    !**** Making mean radial moment
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
          do nth = 1, nthmax
           !**** radial mean velocity
            radial_mean_moment(nr,nsa) = radial_mean_moment(nr,nsa) &
                                -( 0.d0 &
                                ! + Drrfow(nth,np,nr,nsa)&
                                ! * dfdrhom(nth,np,nr,nsa)*1.d20 &
                                ! + Drpfow(nth,np,nr,nsa) &
                                ! * dfdp(nth,np,nr,nsa)*1.d20 &
                                ! + Drtfow(nth,np,nr,nsa) &
                                ! * dfdthm(nth,np,nr,nsa)*1.d20&
                                ) &
                                ! * JIR(nth,np,nr,nsa) &
                                / JIR(nth,np,nr,nsa) &
                                * delp(ns) * delthm(nth,np,nr,nsa) &

                                ! unit converter [J] to [keV] &
                                !* Frr_j(nth,np,nr,nsa) &
                                + Frrfow(nth,np,nr,nsa) &
                                * fnsp_l(nth,np,nr,nsa)*1.d20 &
                                !* JIR(nth,np,nr,nsa) &
                                / JIR(nth,np,nr,nsa) &
                                * delp(ns) * delthm(nth,np,nr,nsa) &
                                !** unit [keV] [22/6/6]
                                + 0.d0
          end do
        end do
        ! write(*,*)radial_mean_moment(nr,nsa)
        ! radial_mean_moment(nr,nsa) = radial_mean_moment(nr,nsa)/(rnsl(nr,nsa)*1.d20)
        ! write(*,*)radial_mean_moment(nr,nsa)
      end do
    end do
    !****Integration over moment(np) and pitch angle(nth)
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
          do nth = 1, nthmax
            ! PV = sqrt(1.d0+theta0(nsa)*pm(np,nsa)**2)
            PV = sqrt(1.d0+theta0(nsa)*(pm(np,nsa)-mean_moment(nr,nsa)/PTFP0(nsa))**2)
            K = (PV-1.d0)*AMFP(nsa)*vc**2 &
            ! K = (pm(np,nsa)*PTFP0(nsa)-mean_moment(nr,nsa))**2.d0/AMFP(nsa)/2.d0 &
            ! K = (pm(np,nsa)*PTFP0(nsa))**2.d0/AMFP(nsa)/2.d0 &
                                / (AEE*1.D3)*2.d0/3.d0 * 1.d0

            if (nr /= 1 .and. nr /= nrmax) then
              JIp = (JIR(nth,np,nr+1,nsa)+JIR(nth,np,nr,nsa))*0.5d0
              JIm = (JIR(nth,np,nr,nsa)+JIR(nth,np,nr-1,nsa))*0.5d0
              Drr_ = ((Drrfow(nth,np,nr+1,nsa)+Drrfow(nth,np,nr,nsa))/JIp*0.5d0 &
                   + (Drrfow(nth,np,nr,nsa)+Drrfow(nth,np,nr-1,nsa))/JIm*0.5d0)*0.5d0
              Frr_ = ((Frrfow(nth,np,nr+1,nsa)+Frrfow(nth,np,nr,nsa))/JIp*0.5d0 &
                   + (Frrfow(nth,np,nr,nsa)+Frrfow(nth,np,nr-1,nsa))/JIm*0.5d0)*0.5d0
            else if (nr == 1) then
              JIp = (JIR(nth,np,nr+1,nsa)+JIR(nth,np,nr,nsa))*0.5d0
              JIm = (JIR(nth,np,nr,nsa)+JIR(nth,np,nr,nsa))*0.5d0
              Drr_ = ((Drrfow(nth,np,nr+1,nsa)+Drrfow(nth,np,nr,nsa))/JIp*0.5d0 &
                   + (Drrfow(nth,np,nr,nsa)+Drrfow(nth,np,nr,nsa))/JIm*0.5d0)*0.5d0
              Frr_ = ((Frrfow(nth,np,nr+1,nsa)+Frrfow(nth,np,nr,nsa))/JIp*0.5d0 &
                   + (Frrfow(nth,np,nr,nsa)+Frrfow(nth,np,nr,nsa))/JIm*0.5d0)*0.5d0
            else
              JIp = (JIR(nth,np,nr,nsa)+JIR(nth,np,nr,nsa))*0.5d0
              JIm = (JIR(nth,np,nr-1,nsa)+JIR(nth,np,nr,nsa))*0.5d0
              Drr_ = ((Drrfow(nth,np,nr+1,nsa)+Drrfow(nth,np,nr,nsa))/JIp*0.5d0 &
                   + (Drrfow(nth,np,nr-1,nsa)+Drrfow(nth,np,nr,nsa))/JIm*0.5d0)*0.5d0
              Frr_ = ((Frrfow(nth,np,nr+1,nsa)+Frrfow(nth,np,nr,nsa))/JIp*0.5d0 &
                   + (Frrfow(nth,np,nr-1,nsa)+Frrfow(nth,np,nr,nsa))/JIm*0.5d0)*0.5d0
            end if
              ! Drr_ = Drrfow(nth,np,nr,nsa)/JIR(nth,np,nr,nsa)
              ! Frr_ = Frrfow(nth,np,nr,nsa)/JIR(nth,np,nr,nsa)

            !**** for dfdrhom is made of fnsp*dVI
            heatfow_out(nr,nsa) = heatfow_out(nr,nsa) &
                                !- (pm(np,nsa)*ptfp0(nsa))**2 &
                                !/ (2*AMFP(nsa)) &
                                - k &
                                *( 0.d0&
                                + Drr_&
                                * dfdrhom(nth,np,nr,nsa)*1.d20 &
                                ! + Drpfow(nth,np,nr,nsa) &
                                ! * dfdp(nth,np,nr,nsa)*1.d20 &
                                ! + Drtfow(nth,np,nr,nsa) &
                                ! * dfdthm(nth,np,nr,nsa)*1.d20&
                                ) &
                                ! / JIR(nth,np,nr,nsa)**2.d0 &
                                ! / JIR(nth,np,nr,nsa) &
                                * delp(ns) * delthm(nth,np,nr,nsa) &

                                !+ (pm(np,nsa)*ptfp0(nsa))**2 &
                                !/ (2*AMFP(nsa) )&
                                
                                - K &
                                ! unit converter [J] to [keV] &
                                !* Frr_j(nth,np,nr,nsa) &
                                * Frr_ &
                                * fnsp_l(nth,np,nr,nsa)*1.d20 &
                                ! / JIR(nth,np,nr,nsa)**2.d0 &
                                ! / JIR(nth,np,nr,nsa) &
                                * delp(ns) * delthm(nth,np,nr,nsa) &
                                !** unit [keV] [22/6/6]
                                + 0.d0
              ! add_heat(nr,nsa)  = K &
              !                   * fnsp_l(nth,np,nr,nsa) * 1.d20 &
              !                   ! * radial_mean_moment(nr,nsa)&!/AMFP(nsa) &
              !                   * JIR(nth,np,nr,nsa) * delp(ns) * delthm(nth,np,nr,nsa) 
! write(*,*)"add term"
          end do
        end do
        heatfow_out(nr,nsa) = ( 0.d0&
                            + heatfow_out(nr,nsa) &
                            ! -5.d0/2.d0*Ta(nr,nsa)*radial_mean_moment(nr,nsa)&
                            +0.d0 )
        ! if(nsa ==1)heatfow_out(nr,nsa) = heatfow_out(nr,nsa)*sqrt(AMFP(1)/AMFP(2))
        ! heatfow_out(nr,nsa) = heatfow_out(nr,nsa)-Ta(nr,nsa)*rnsl(nr,nsa)*radial_mean_moment(nr,nsa)*1.d20!&
        ! write(*,*)"ta_radial" ,TA(nr,nsa)*radial_mean_moment(nr,nsa)
                            ! + 2.5d0*Ta(nr,nsa)*(rnsl(nr,nsa)*1.d20)*radial_mean_moment(nr,nsa)
        ! write(*,*)"add heat:",add_heat(nr,nsa)*radial_mean_moment(nr,nsa)
        if(nsa ==1)heatfow_out(nr,nsa) = heatfow_out(nr,nsa)*sqrt(AMFP(1)/AMFP(2))
        chi_a(nr,nsa) = -heatfow_out(nr,nsa)/(dTadr(nr,nsa)*rnsl(nr,nsa)*1.d20)
      end do
    end do

    ! !**** making radial mean vector
    ! call make_moment_vector(moment_vector)
    ! do nsa = 1, nsamax
    !   ns = ns_nsa(nsa)
    !   do nr = 1, nrmax
    !     do np = 1, npmax
    !       if ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
    !       do nth = 1, nthmax
    !        !**** making mean vector
    !         radial_mean_moment_vector(nr,nsa) = moment_vector(3,nth,np,nr,nsa) &
    !                             * fnsp_l(nth,np,nr,nsa)*1.d20 &
    !                             * JIR(nth,np,nr,nsa) &
    !                             * delp(ns) * delthm(nth,np,nr,nsa)
    !         sum_ = sum_ + fnsp_l(nth,np,nr,nsa)*1.d20 &
    !                             * JIR(nth,np,nr,nsa) &
    !                             *delp(ns) * delthm(nth,np,nr,nsa)
    !       end do
    !     end do
    !     radial_mean_moment_vector(nr,nsa) = radial_mean_moment_vector(nr,nsa)/sum_
    !     write(*,*)radial_mean_moment_vector(nr,nsa)
    !   end do
    ! end do

  end subroutine integral_Heatdiff

  ! subroutine Heat_rw_baverage(heatrwav, heatrw)
  !---------------------------------------------
  ! heat diffusion coef calcul module
  ! This module corresponds to integral_Heatdiff
  !---------------------------------------------

!     use fpcomm
!     use fowcomm

!     implicit none
!     double precision,dimension(nrmax,nsamax),intent(out) :: heatrw, heatrwav
!     double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom, fnsp_l
!     double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
!     double precision,dimension(nrmax,nsamax) :: nu, Ta, Drw, vth
!     double precision step_len, step_time, eps_t, rho_theta, Baxis, pv, rho
!     double precision dVI, sumVI
!     integer nth, np, nr, nsa, ns

!     !**** initialization ****
!     Baxis = Bin_rg(1)
!     heatrw(:,:) = 0.d0
!     heatrwav(:,:) = 0.d0
!     sumVI = 0.d0
!  !**** for new dfdrhom module [2022/3/4]
!      fnsp_l(:,:,:,:)=0.d0
!      do nsa = 1, nsamax
!        ns = ns_nsa(nsa)
!        do nr= 1, nrmax
!          do np = 1, npmax
!            if ( pm(np,ns) > fact_bulk ) exit
!            do nth = 1, nthmax
!            fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa)
!            end do
!          end do
!        end do
!      end do

!      !**** first order derivative
!      do nsa = 1, nsamax
!        do np = 1, npmax
!          do nth = 1, nthmax
!            call first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp_l(nth,np,:    ,nsa), rm)
!          end do
!        end do
!      end do
!  !**** end of new dfdrhom module

!     !**** Temperature and flecency make
!     do nsa = 1, nsamax
!       do nr = 1, nrmax
!         Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)
!         !**** Ta(temperature)[J]
!         ! /RKEV = 1/(AEE*1.D3) unit converter [keV]to [J]
!       end do
!     end do

!    call nu_effective(nu)

!     !**** Calculation of particle diff coef
!     do nsa = 1, nsamax
!       do nr = 1, nrmax
!         eps_t = rm(nr)*RA/RR
!         do np = 1, npmax
!           step_time = 1.d0/nu(nr,nsa)
!           !** thermal rho
! !          vth(nr,nsa) = SQRT( 2.d0*Ta(nr,nsa)/AMFP(nsa)*RKEV)
! !          rho = vth(nr,nsa)*AMFP(nsa)/(aefp(nsa)*Baxis)
!           !** original rho
!           rho = (pm(np,nsa)*ptfp0(nsa) ) &!/sqrt(AEE*1.d3)) &
!               / (aefp(nsa)*Baxis)
!           !****Unit [J]

!           do nth = 1, nthmax
!             step_len = rho * safety_factor(nr)/sqrt(eps_t)
!             Drwlocal(nth,np,nr,nsa) = step_len**2/step_time*sqrt(eps_t)
!           end do
!         end do
!       end do
!     end do

!     call bounce_average_for_Drw(Drwba, Drwlocal)

!     !**** Heat diffusion coef calcul
!     do nsa = 1, nsamax
!       do nr = 1, nrmax
!         do np = 1, npmax
!           do nth = 1, nthmax
!             !**** Heat flux for FOW momentum not heat velocity [2022/3/22]
!             heatrw(nr,nsa) = heatrw(nr,nsa)&
!                            - (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
!                            * Drwlocal(nth,np,nr,nsa) &
!                            * JIR(nth,np,nr,nsa) &
!                            * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa)*1.d20 &
!                            / (AEE*1.D3)&  !** Unit convert [J] to [keV]
!                            * delp(ns)*delthm(nth,np,nr,nsa)
!             heatrwav(nr,nsa) = heatrwav(nr,nsa) &
!                              - (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
!                              * Drwba(nth,np,nr,nsa) &
!                              * JIR(nth,np,nr,nsa) &
!                              * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa)*1.d20 &
!                              / (AEE*1.D3)& !** Unit convert [J] to [keV]
!                              *delp(ns)* delthm(nth,np,nr,nsa)

! !            sumVI = sumVI - dfdrhom(nth,np,nr,nsa)&
! !                  * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa))/(AEE*1.D3)
! !*** Heat flux of heat velocity [2022/3/22]
! !            heatrw(nr,nsa) = heatrw(nr,nsa)&
! !                           + (vth(nr,nsa)*AMFP(nsa))**2/(2*AMFP(nsa)) &
! !                           * Drwlocal(nth,np,nr,nsa) &
! !                           * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
! !                           / (AEE*1.D3)  !**** Unit convert [J] to [keV]
! !            heatrwav(nr,nsa) = heatrwav(nr,nsa) &
! !                             + (vth(nr,nsa)*AMFP(nsa))**2/(2*AMFP(nsa)) &
! !                             * Drwba(nth,np,nr,nsa) &
! !                             * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
! !                             / (AEE*1.D3)  !****Unit convert [J] to [keV]
! !
!           end do
!         end do
! !        heatrw(nr, nsa)  = heatrw(nr, nsa)/sumVI
! !        heatrwav(nr,nsa) = heatrwav(nr,nsa)/sumVI
!       end do
!     end do

  ! end subroutine Heat_rw_baverage

  subroutine Rosenbluth_Hazeltine_Hinton_Neoclass(heatba,heatpla,chi_neo_ba,chi_neo_pla,Sr_ba,Sr_pla,Dnba,Dnpla)
  !---------------------------------------------------------------
  ! subroutine for neoclassical heat flux in banana, plateau region
  ! calculate heat flux and diffusion coefficient
  ! Plasma Transport in Toroidal Confinement Systems(1972)
  ! main calcul using[J]
  ! Calculation with [J]
  !--------------------------------------------------------------

    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: heatba, heatpla
    double precision,dimension(nrmax,nsamax),intent(out)  :: chi_neo_ba, chi_neo_pla
    ! double precision,dimension(nrmax,nsamax),intent(out)  :: inv_asp, pla_parm
    double precision,dimension(nrmax,nsamax),intent(out)  :: Dnba, Dnpla
    double precision,dimension(nrmax,nsamax),intent(out)  :: Sr_ba, Sr_pla
    double precision,dimension(nrmax,nsamax) :: Ta, dTadr, dndr
    double precision, dimension(nsamax) :: rho_a
    double precision fact, fact_s, tau_ele,tau_i, eps_t, Baxis, B_p
    double precision fact_ba,fact_pla,rho_e!,eps_t,Baxis,B_p
    integer nr, nsa

    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0) !** for SI unit collisional time
    Baxis = Bin_rg(1) ! approximation on B by Baxis

    !****temperature make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)
        !Ta(temperature)[J]
      end do
    end do

    !****first order derivation
    do nsa = 1, nsamax
      call first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
      call first_order_derivative(dndr(:,nsa), rnsl(:, nsa), rm)
    end do

    !****calculate flux and diffusion coefficient
    do nsa = 1, nsamax
      do nr = 1, nrmax

        eps_t = rm(nr)*RA/RR
        B_p = rm(nr)*RA*BB/(safety_factor(nr)*RR)

        tau_ele = fact*sqrt(AMFP(1))*((Ta(nr,1))**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**2*lnlam(nr,2,1)*AEE**2)
              ! write(*,'(e32.2)')tau_ele
        tau_i   = fact*sqrt(2.d0)*sqrt(AMFP(2))*((Ta(nr,2))**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2))

        rho_e = sqrt(2*Ta(nr,1)/AMFP(1))*AMFP(1)/(AEE*Baxis)
        ! rho_e = sqrt(Ta(nr,1)/AMFP(1))*AMFP(1)/(AEE*Baxis)
        rho_a(nsa) = sqrt(2*Ta(nr,nsa)/AMFP(nsa))*AMFP(nsa)/(AEFP(nsa)*Baxis)
        ! rho_a(nsa) = sqrt(Ta(nr,nsa)/AMFP(nsa))*AMFP(nsa)/(AEFP(nsa)*Baxis)
        ! inv_asp(nr,nsa) = eps_t**(1.5) !!!!! for plateau
        ! if (nsa == 1)then !!!! for plateau index number must 1.0 form
        !   pla_parm(nr,nsa) = (rm(nr)*RA*Baxis/(B_p*sqrt(3.d0*Ta(nr,nsa)/AMFP(1)))/tau_ele)!/inv_asp(nr,nsa)!**(1.0/3.0)
        ! else
        !   pla_parm(nr,nsa) = (rm(nr)*RA*Baxis/(B_p*sqrt(3.d0*Ta(nr,nsa)/AMFP(nsa)))/tau_i)!/inv_asp(nr,nsa)!**(1.0/3.0)
        ! end if

        !==== Calculation of Particle Fluxes and Coefficients ====
        fact_ba = (safety_factor(nr)**2)*1.d20*rnsl(nr,nsa) &
             *(rho_e**2)/((eps_t**1.5)*tau_ele)*RKEV ! modified by anzai[*RKEV]
        fact_pla = -sqrt(pi)*(eps_t**2)*Ta(nr,1)*RKEV &! modified by anzai [*RKEV]
             *rho_e*1.d20*rnsl(nr,nsa)/(2*AEE*B_p*rm(nr)*RA)

        Sr_ba(nr,nsa) = fact_ba*(-1.22d0*(1+Ta(nr,2)/Ta(nr,1)) &
                      * 1.d20 * dndr(nr, nsa) &
                      ! + 4.3d-1*dTadr(nr,1)/Ta(nr,1) &
                      ! + 1.9d-1*dTadr(nr,2)/Ta(nr,1) &
                      )
        !** neglect E_para term

        Sr_pla(nr,nsa) = fact_pla*((1+Ta(nr,2)/Ta(nr,1)) &
                       * 1.d20 * dndr(nr,nsa) &
                       / (1.d20*rnsl(nr,nsa)) &
                       + 1.5d0*dTadr(nr,1)/Ta(nr,1) + 1.5d0*dTadr(nr,2)/Ta(nr,1) )
        !** neglect B term

        Dnba(nr,nsa) = - Sr_ba(nr, nsa)/(1.d20*dndr(nr, nsa))
        Dnpla(nr,nsa) = - Sr_pla(nr, nsa)/(1.d20*dndr(nr, nsa))

        !==== Calculation of Heat Fluxes and Coeffisients ====
        if (nsa == 1) then
          fact_s = (safety_factor(nr)**2) &
             *(rho_a(1)**2)/((eps_t**1.5)*tau_ele)
          heatba(nr,nsa) = fact_s*rnsl(nr,nsa)*1.D20*Ta(nr,1)/RKEV &
            * ( 0.d0 &
              ! + 1.53d0*(1+Ta(nr,2)/Ta(nr,1))*dndr(nr, nsa))/rnsl(nr,nsa) &
              - 1.81d0*dTadr(nr,1)/Ta(nr,1) &
              !  - 0.27d0*dTadr(nr,2)/Ta(nr,1) &
              ) !test for [keV]
        !** neglect E_para term

          heatpla(nr,nsa) =  sqrt(PI)/4*eps_t**2*Ta(nr,nsa)/AEE/1.d3 &
               / (AEE*B_p)*rho_a(nsa)*rnsl(nr,nsa)*1.d20*Ta(nr,nsa)/(AEE*1.d3) &
               / (rm(nr)*RA)*( (1+Ta(nr,2)/Ta(nr,1))*dndr(nr,nsa)*1.d20 &
               / (rnsl(nr,nsa)*1.d20) &
               + 7.5d0*dTadr(nr,nsa)/Ta(nr,nsa) &
               + 1.5d0*dTadr(nr,2)/Ta(nr,nsa))*RKEV!test [keV]
        !** neglect E_para term

          ! chi_neo_ba(nr,nsa) = 1.81d0*(safety_factor(nr)**2)*rho_a(1)**2/(eps_t**1.5*tau_ele)
          chi_neo_ba(nr,nsa) = -heatba(nr,nsa)/(rnsl(nr,nsa)*1.d20)/dTadr(nr,nsa)*AEE*1.d3
          chi_neo_pla(nr,nsa)= -sqrt(PI)/4*eps_t**2*Ta(nr,nsa)/AEE/1.d3&
               / (AEE*B_p)*rho_a(nsa)&!*rnsl(nr,nsa)*1.d20 &
               / (rm(nr)*RA)*7.5d0*RKEV
        else
          fact_s = (safety_factor(nr)**2) &
              *(rho_a(nsa)**2)/((eps_t**1.5)*tau_i)
          heatba(nr,nsa) = - 0.68d0*fact_s*(1 + 0.48d0*sqrt(eps_t)) &
              *rnsl(nr,nsa)*1.d20*dTadr(nr,2)/RKEV !test for [keV]
          heatpla(nr,nsa) = - 1.5d0*sqrt(pi) &
                        * eps_t**2*Ta(nr,nsa)/AEE/1.d3/(AEE*B_p) &
                        * rho_a(nsa)/(rm(nr)*RA)*rnsl(nr,nsa)*1.d20 &
                        * dTadr(nr,nsa)/(AEE*1.d3)*RKEV !test for [keV]
          chi_neo_ba(nr,nsa) = 0.68d0*(1+0.48d0*sqrt(eps_t))*(safety_factor(nr)**2)*rho_a(nsa)**2/(eps_t**1.5*tau_i)
          chi_neo_pla(nr,nsa)= 1.5d0*sqrt(pi)*eps_t**2*Ta(nr,nsa)/AEE/1.d3*rho_a(nsa)&!*rnsl(nr,nsa)*1.d20 &
                             / (AEFP(nsa)*B_p*rm(nr)*RA)*RKEV
        end if
     end do
   end do

    ! !****calculate flux and diffusion coefficient
    ! do nsa = 1, nsamax
    !   do nr = 1, nrmax

    !     eps_t = rm(nr)*RA/RR
    !     tau_ele = fact*sqrt(AMFP(1))*((Ta(nr,1))**1.5d0)/ &
    !          (rnsl(nr,2)*1.d20*AEFP(2)**2*lnlam(nr,2,1)*AEE**2)
    !     tau_i   = fact*sqrt(2.d0)*sqrt(AMFP(2))*((Ta(nr,2))**1.5d0)/ &
    !          (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2))
    !     ! tau_ele = fact*sqrt(AMFP(1))*((Ta(nr,1))**1.5d0)/ &
    !     !      (rnsl(nr,2)*1.d20*AEFP(2)**2*AEE**2*lnlam(nr,2,1))
    !     ! tau_i   = fact*sqrt(2.d0)*sqrt(AMFP(2))*((Ta(nr,2))**1.5d0)/ &
    !     !      (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2))
    !     ! !**** *RKEV converts [keV] to [J]

    !     ! rho_a(nsa) = sqrt(Ta(nr,nsa)/AMFP(nsa))*AMFP(nsa)/(AEFP(nsa)*Baxis)
    !     !**** /RKEV converts [J] to [keV]
    !     !****calculate heat flux
    !     if (nsa == 1) then
    !       fact_s = (safety_factor(nr)**2) &
    !          *(rho_a(1)**2)/((eps_t**1.5)*tau_ele)
    !       heatba(nr,nsa) = fact_s*rnsl(nr,nsa)*1.D20*Ta(nr,1)/RKEV &
    !         * ((1.53d0*(1+Ta(nr,2)/Ta(nr,1))*dndr(nr, nsa))/rnsl(nr,nsa) &
    !           - 1.81d0*dTadr(nr,1)/Ta(nr,1) &
    !           - 0.27d0*dTadr(nr,2)/Ta(nr,1))!test for [keV]
    !     ! neglect E_para term

    !       heatpla(nr,nsa) =  sqrt(PI)/4*eps_t**2*Ta(nr,nsa)/AEE/1.d3 &
    !            / (AEE*B_p)*rho_a(nsa)*rnsl(nr,nsa)*1.d20*Ta(nr,nsa)/(AEE*1.d3) &
    !            / (rm(nr)*RA)*( (1+Ta(nr,2)/Ta(nr,1))*dndr(nr,nsa)*1.d20 &
    !            / (rnsl(nr,nsa)*1.d20) &
    !            + 7.5d0*dTadr(nr,nsa)/Ta(nr,nsa) &
    !            + 1.5d0*dTadr(nr,2)/Ta(nr,nsa))*RKEV!test [keV]
    !     ! neglect E_para term
    !       ! chi_neo_ba(nr,nsa) = 1.81d0*(safety_factor(nr)**2)*rho_a(1)**2/(eps_t**1.5*tau_ele)
    !       chi_neo_ba(nr,nsa) = heatba(nr,nsa)/(rnsl(nr,nsa)*1.d20)
    !       chi_neo_pla(nr,nsa)= -sqrt(PI)/4*eps_t**2*Ta(nr,nsa)/AEE/1.d3&
    !            / (AEE*B_p)*rho_a(nsa)&!*rnsl(nr,nsa)*1.d20 &
    !            / (rm(nr)*RA)*7.5d0*RKEV
    !     else
    !       fact_s = (safety_factor(nr)**2) &
    !           *(rho_a(nsa)**2)/((eps_t**1.5)*tau_i)
    !       heatba(nr,nsa) = - 0.68d0*fact_s*(1 + 0.48d0*sqrt(eps_t)) &
    !           *rnsl(nr,nsa)*1.d20*dTadr(nr,2)/RKEV !test for [keV]
    !       heatpla(nr,nsa) = - 1.5d0*sqrt(pi) &
    !                     * eps_t**2*Ta(nr,nsa)/AEE/1.d3/(AEE*B_p) &
    !                     * rho_a(nsa)/(rm(nr)*RA)*rnsl(nr,nsa)*1.d20 &
    !                     * dTadr(nr,nsa)/(AEE*1.d3)*RKEV !test for [keV]
    !       chi_neo_ba(nr,nsa) = 0.68d0*(1+0.48d0*sqrt(eps_t))*(safety_factor(nr)**2)*rho_a(nsa)**2/(eps_t**1.5*tau_i)
    !       chi_neo_pla(nr,nsa)= 1.5d0*sqrt(pi)*eps_t**2*Ta(nr,nsa)/AEE/1.d3*rho_a(nsa)&!*rnsl(nr,nsa)*1.d20 &
    !                          / (AEFP(nsa)*B_p*rm(nr)*RA)*RKEV
    !     end if
    !     ! chi_neo_ba      = -Sr_ba(nr,nsa)/(dTadr(nr,nsa)/RKEV)
    !     ! chi_neo_pla     = -Sr_pla(nr,nsa)/(dTadr(nr,nsa)/RKEV)
    !   end do
    ! end do

  end subroutine Rosenbluth_Hazeltine_Hinton_Neoclass

  subroutine Chang_Hinton_neoclass(Sr_ch, chi_ch)
  !---------------------------------------------------------------
  !subroutine for neoclassical heat flux by Chang and Hinton
  ! Effect of impurity particles on the finite aspect-ratio neoclassical
  ! ion thermal conductivity in a tokamak(1986)
  !--------------------------------------------------------------

    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: Sr_ch
    double precision,dimension(nrmax,nsamax),intent(out)  :: chi_ch
    double precision,dimension(nrmax,nsamax) :: Ta, dTadr
    double precision fact, tau_i, eps_t, Baxis, B_p
    double precision fact_a, fact_1, fact_2, f_1, f_2
    double precision sh_shift, alpha, rho_i, nu_i, mu_i
    integer nr, nsa

    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bin_rg(1) ! approximation on B by Baxis
    sh_shift = 0.0D0
    alpha = 0.0d0 !! n_I*Z_I^2/(n_i*Z_i^2) for impurity ion

    !****temperature make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)
        !Ta(temperature)[J]
      end do
    end do

    !****first order derivation
    do nsa = 1, nsamax
      call first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
    end do
            !  write(*,*)aefp(1)
            !  write(*,*)tau_ee
    !****calculate flux and diffusion coefficient
    do nsa = 1, nsamax
      do nr = 1, nrmax
        ! B_p = rm(nr)*RA*BB/(safety_factor(nr)*RR)
        eps_t = rm(nr)*RA/RR
        tau_i   = fact*sqrt(2.d0)*sqrt(AMFP(2))*((Ta(nr,2))**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2))
            !  write(*,'(e32.2)')tau_i
        rho_i = sqrt(2.d0*Ta(nr,2)/AMFP(2))*AMFP(2)/(AEFP(2)*Baxis)
        ! rho_i = sqrt(Ta(nr,2)/AMFP(2))*AMFP(2)/(AEFP(2)*Baxis)
        nu_i = RR*safety_factor(nr)*eps_t**(-1.5d0)*sqrt(2.d0*Ta(nr,2)*(AMFP(2)))*tau_i**(-1.d0)
        mu_i = nu_i*(1.d0 + 1.54d0 * alpha)

        fact_a = eps_t**(-1.5d0)*(rho_i**(2.d0))*tau_i**(-1.d0)*safety_factor(nr)**(2.d0)
        fact_1 = (0.66d0*(1+1.54d0*alpha) + (1.88d0*eps_t**(0.5) - 1.54d0*eps_t)*(1+3.75d0*alpha)) * &
                 (1.d0 + 1.03d0*mu_i**(0.5d0) + 0.31d0*mu_i)**(-1.d0)
        fact_2 = 0.59d0*mu_i*eps_t*(1.d0+0.74*mu_i*eps_t**(1.5))**(-1.d0) * &
                 ( 1.d0 + (1.33d0*alpha*(1.d0 + 0.6d0*alpha))*(1.d0+1.79d0*alpha)**(-1.d0))

        f_1 = (1.d0 + 1.5d0*(eps_t**(2.d0) + eps_t * sh_shift) + 3.d0*(8.d0)**(-1.d0)*eps_t**(3.d0)*sh_shift) * &
              (1 + 0.5 * eps_t * sh_shift)**(-1.d0)
        f_2 = sqrt(1.d0 - eps_t**(2.d0))*(1.d0 + eps_t*sh_shift*2**(-1.d0)) * &
              (1.d0 + sh_shift*(sqrt(1.d0 - eps_t**(2.d0)) - 1.D0)*eps_t**(-1.d0))**(-1.d0)

        if (nsa == 1) then !** Reference Value
          chi_ch(nr,nsa) = fact_a*(fact_1*f_1 + fact_2*(f_1 - f_2))*sqrt(AMFP(1)/AMFP(2))
          Sr_ch(nr,nsa) = -rnsl(nr,nsa)*1.d20*chi_ch(nr,nsa)*dTadr(nr,nsa)/AEE/1.d3
        else 
          chi_ch(nr,nsa) = fact_a*(fact_1*f_1 + fact_2*(f_1 - f_2))
          Sr_ch(nr,nsa) = -rnsl(nr,nsa)*1.d20*chi_ch(nr,nsa)*dTadr(nr,nsa)/AEE/1.d3
        end if

     end do
   end do

  end subroutine Chang_Hinton_neoclass

  subroutine Hinton_Hazeltine_neoclass(Sr_hh, Dr_hh, heat_hh, chi_hh)
  !---------------------------------------------------------------
  ! subroutine for neoclassical flux by Hinton and Hazeltine
  ! Theory of plasma transport in troidal confinement systems(1976) 
  ! Reviews of Modern Physics 48, 2, part 1
  ! neglect E_parallel term
  !--------------------------------------------------------------

    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: Sr_hh, heat_hh
    double precision,dimension(nrmax,nsamax),intent(out)  :: Dr_hh, chi_hh
    double precision,dimension(nrmax,nsamax) :: Ta, dTadr, Pa_, dPadr, dNadr
    double precision,dimension(3,3,3) :: K_, a_, b_, c_
    double precision,dimension(3,3,3) :: K_eff
    !**** dimension K is (species, m, n) in the paper p298(1976)
    double precision fact, tau_e, tau_i, eps_t, Baxis, B_p
    double precision fact_i
    double precision rho_ep, rho_ip 
    double precision inner_product, A_force
    double precision nu_i, nu_e
    integer nr, nsa, i, j, k

    !**** initialization
    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bin_rg(1) ! approximation on B by Baxis
    K_eff(:,:,:)=0.d0
    K_(:,:,:) = 0.000001d0 !** avoid floating error
    a_(:,:,:) = 0.000001d0 !** avoid floating error
    b_(:,:,:) = 0.000001d0 !** avoid floating error
    c_(:,:,:) = 0.000001d0 !** avoid floating error

    call make_HHCoef(K_, a_, b_, c_)

    !****temperature make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)
        !**Ta(temperature)[J]
        Pa_(nr,nsa) = rnsl(nr,nsa)*1.d20*Ta(nr,nsa) 
        !**Pressure[Pa]
      end do
    end do

    !****first order derivation
    do nsa = 1, nsamax
      call first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
      call first_order_derivative(dPadr(:,nsa), Pa_(:,nsa), rm)
      call first_order_derivative(dNadr(:,nsa), rnsl(:, nsa), rm)
    end do
            !  write(*,*)aefp(1)
            !  write(*,*)tau_ee
    !****calculate flux and diffusion coefficient
    do nsa = 1, nsamax
      do nr = 1, nrmax
        ! B_p = rm(nr)*RA*BB/(safety_factor(nr)*RR)
        B_p = rm(nr)*RA*Baxis/(safety_factor(nr)*RR)
        eps_t = rm(nr)*RA/RR
        tau_e = fact*sqrt(AMFP(1))*((Ta(nr,1))**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**2*lnlam(nr,2,1)*AEE**2) !** for SI unit
        ! tau_e = 3.d0/4.d0*(2.d0*pi)**(-0.5d0)*sqrt(AMFP(1))*((Ta(nr,1))**1.5d0)/ &
            !  (rnsl(nr,2)*1.d20*AEFP(2)**2*lnlam(nr,2,1)*AEE**2) !** for cgs unit
        tau_i   = fact*sqrt(2.d0)*sqrt(AMFP(2))*((Ta(nr,2))**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2)) !** for SI unit
        ! tau_i   = 3.d0/4.d0*sqrt(pi)**(-1)*sqrt(AMFP(2))*((Ta(nr,2))**1.5d0)/ &
            !  (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2))*EPS0**2.d0 !** for cgs unit
            !  write(*,'(e32.2)')tau_i
        ! nu_e = sqrt(2.d0)*rm(nr)*RA*Baxis/(B_p*sqrt(2.d0*Ta(nr,1)/AMFP(1))*tau_e*eps_t**(1.5d0))
        ! nu_i = sqrt(2.d0)*rm(nr)*RA*Baxis/(B_p*sqrt(2.d0*Ta(nr,2)/AMFP(2))*tau_i*eps_t**(1.5d0))
        nu_e = sqrt(2.d0)*rm(nr)*RA*Baxis/(B_p*sqrt(2.d0*Ta(nr,1)/AMFP(1))*tau_e*eps_t**(1.5d0))
        nu_i = sqrt(2.d0)*rm(nr)*RA*Baxis/(B_p*sqrt(2.d0*Ta(nr,2)/AMFP(2))*tau_i*eps_t**(1.5d0))
        inner_product = (1.17d0 - 0.35*nu_i**(0.5d0))*(1.d0 + 0.7*nu_i**(0.5d0))**(-1.d0)
        inner_product = (inner_product -2.1d0*nu_i**(2.d0)*eps_t**(3.d0))/(1.d0 + nu_i**(2.d0)*eps_t**(3.d0)) 
        rho_ep = sqrt(2.d0*Ta(nr,1)/AMFP(1))*AMFP(1)/(AEFP(1)*B_p)!*vc !** for cgs unit
        rho_ip = sqrt(2.d0*Ta(nr,2)/AMFP(2))*AMFP(2)/(AEFP(2)*B_p)!*vc !** for cgs unit
        ! rho_ep = sqrt(Ta(nr,1)/AMFP(1))*AMFP(1)/(AEFP(1)*B_p)!*vc !** for cgs unit
        ! rho_ip = sqrt(Ta(nr,2)/AMFP(2))*AMFP(2)/(AEFP(2)*B_p)!*vc !** for cgs unit

        fact_i = 0.66d0*( 1.d0*(1.d0 + 1.03d0 * nu_i**(0.5d0) + 0.31d0*nu_i)**(-1.d0) + &
                 (eps_t**(3.d0)*(0.74d0)**(2.d0)/(0.31d0)*nu_i) * &
                 (1.d0 + 0.74d0*nu_i*eps_t**(1.5d0))**(-1.d0) )
            !  write(*,'(e32.2)')tau_i

        !**** Making coefficient
        do k = 1, 3
          do j = 1, 3
            do i = 1, 3
              if (k .ne. 3) then
                K_eff(i,j,k)=K_(i,j,k)*(1.d0/(1.d0 + a_(i,j,k)*nu_e**(0.5d0) + b_(i,j,k)*nu_e) + &
                             (eps_t**(3.d0)*c_(i,j,k)**(2.d0)*b_(i,j,k)**(-1.d0)*nu_e)/ &
                             (1.d0 + c_(i,j,k)*nu_e*eps_t**(1.5d0)) )
              else if (k .eq. 3) then
                K_eff(i,j,k)= K_(i,j,k)/(1.d0 + a_(i,j,k)*nu_e**(0.5d0) + b_(i,j,k)*nu_e) / &
                              (1.d0+c_(i,j,k)*nu_e*eps_t**(1.5d0))
              end if
            end do!i
          end do !j
        end do !k
            !  write(*,'(e32.2)')tau_i

        !**** Making Fluxes
        A_force = &
                  dPadr(nr,1)*(Pa_(nr,nsa)**(-1.d0)) &
                  -2.5d0*dTadr(nr,1)*(Ta(nr,1)**(-1.d0))+ &
                  (Ta(nr,2)*(AEFP(2)**(-1.d0)*AEE)*Ta(nr,1)) * &
                  ( &
                  dPadr(nr,2)*(Pa_(nr,2)**(-1.d0))&
                   - inner_product*(1.d0 + nu_e**(2.d0)*eps_t**(2.d0))**(-1.d0) * &
                  dTadr(nr,2) * (Ta(nr,2)**(-1.d0)) )
            !  write(*,'(e32.2)')nu_e!EPS0
        !** Z_i = 1 selected --> nsa = 2
        Sr_hh(nr,nsa) = - rnsl(nr,1)*1.d20 * eps_t**(0.5d0)*rho_ep**(2.d0)*(tau_e)**(-1.d0)* &
                        (K_eff(1,1,1)*A_force &
                        + K_eff(1,1,2)*dTadr(nr,1)*(Ta(nr,1)**(-1.d0)) &
                        )

        if (nsa .eq. 1) then
          heat_hh(nr,nsa) = (-2.5d0*Ta(nr,nsa)*Sr_hh(nr,nsa) &
                            -rnsl(nr,nsa)*1.d20 * &
                            Ta(nr,nsa) * eps_t**(0.5d0)*(rho_ep**(2.d0)*(tau_e)**(-1.d0)) * &
                            ( & 
                            K_eff(1,1,2)*A_force &
                            + K_eff(1,2,2)*dTadr(nr,nsa)*(Ta(nr,nsa)**(-1.d0))))/AEE/1.d3
        else
          heat_hh(nr,nsa) = (- fact_i * rnsl(nr,nsa)*1.d20*eps_t**(0.5d0)*rho_ip**(2.d0)*(tau_i)**(-1.d0)* &
                            dTadr(nr,nsa) &
                             - inner_product* Ta(nr,nsa) * &
                            Sr_hh(nr,nsa)/AEFP(2)*AEE*(1.d0 + nu_e**(2.d0)*eps_t**(3.d0))**(-1.d0) &
                            )/AEE/1.d3
        end if

      Dr_hh(nr,nsa) = - Sr_hh(nr,nsa)/(dNadr(nr,nsa)*1.d20)
      chi_hh(nr,nsa) = - heat_hh(nr,nsa)/(rnsl(nr,nsa)*1.d20)/dTadr(nr,nsa)*AEE*1.d3
     end do
   end do

  end subroutine Hinton_Hazeltine_neoclass

  subroutine make_HHCoef(K_, a_, b_, c_)
  !------------------------------------------
  ! Subroutine for making coefficient in
  ! Theory of plasma transport(1976)
  ! ref: Hinton_Hazeltine_neoclass
  !------------------------------------------
    implicit none
    double precision,dimension(3,3,3),intent(inout) :: K_, a_, b_, c_

    !**** input values of K, a_, b_, c_ in p298
    !** Z_i = 1
    K_(1,1,1)=1.04d0
    K_(1,1,2)=1.20d0
    K_(1,1,3)=2.30d0
    K_(1,2,2)=2.55d0
    K_(1,2,3)=4.19d0  
    K_(1,3,3)=1.83d0

    a_(1,1,1)=2.01d0
    a_(1,1,2)=0.76d0
    a_(1,1,3)=1.02d0
    a_(1,2,2)=0.45d0
    a_(1,2,3)=0.57d0
    a_(1,3,3)=0.68d0

    b_(1,1,1)=1.53d0
    b_(1,1,2)=0.67d0
    b_(1,1,3)=1.07d0
    b_(1,2,2)=0.43d0
    b_(1,2,3)=0.61d0
    b_(1,3,3)=0.32d0

    c_(1,1,1)=0.89d0
    c_(1,1,2)=0.56d0
    c_(1,1,3)=1.07d0
    c_(1,2,2)=0.43d0
    c_(1,2,3)=0.61d0
    c_(1,3,3)=0.66d0

    !** Z_i = 2
    K_(2,1,1)=0.86d0
    K_(2,1,2)=0.95d0
    K_(2,1,3)=1.87d0
    K_(2,2,2)=1.99d0
    K_(2,2,3)=3.72d0
    K_(2,3,3)=1.56d0

    a_(2,1,1)=2.18d0
    a_(2,1,2)=0.78d0
    a_(2,1,3)=0.89d0
    a_(2,2,2)=0.46d0
    a_(2,2,3)=0.52d0
    a_(2,3,3)=0.56d0

    b_(2,1,1)=1.17d0
    b_(2,1,2)=0.50d0
    b_(2,1,3)=0.62d0
    b_(2,2,2)=0.26d0
    b_(2,2,3)=0.34d0
    b_(2,3,3)=0.25d0

    c_(2,1,1)=0.79d0
    c_(2,1,2)=0.51d0
    c_(2,1,3)=0.69d0
    c_(2,2,2)=0.34d0
    c_(2,2,3)=0.38d0
    c_(2,3,3)=0.58d0

    !** Z_i = 4
    K_(3,1,1)=0.76d0
    K_(3,1,2)=0.83d0
    K_(3,1,3)=1.65d0
    K_(3,2,2)=1.71d0
    K_(3,2,3)=3.54d0
    K_(3,3,3)=1.42d0

    a_(3,1,1)=2.30d0
    a_(3,1,2)=0.80d0
    a_(3,1,3)=0.79d0
    a_(3,2,2)=0.46d0
    a_(3,2,3)=0.48d0
    a_(3,3,3)=0.47d0

    b_(3,1,1)=0.98d0
    b_(3,1,2)=0.42d0
    b_(3,1,3)=0.56d0
    b_(3,2,2)=0.22d0
    b_(3,2,3)=0.33d0
    b_(3,3,3)=0.20d0

    c_(3,1,1)=0.74d0
    c_(3,1,2)=0.48d0
    c_(3,1,3)=0.51d0
    c_(3,2,2)=0.30d0
    c_(3,2,3)=0.28d0
    c_(3,3,3)=0.51d0
    
    end subroutine make_HHCoef
!===========================================
! for moment vector quantities calculation
!===========================================
  subroutine make_moment_vector(moment_vector)
  !---------------------------------------------------
  ! maybe using interpole module is good [2022/11/7]
  ! pm(np,nsa) --> momentvector
  ! real moment is pm(np,nsa)*pdfp0(nsa)
  !---------------------------------------------------

    use fpcomm
    use fowcomm
    ! use foworbit
    ! use libspl3d

    implicit none
    ! type(orbit):: orbit_m
    real(rkind),dimension(3,3,max_stp) :: dIdu
    real(rkind),dimension(3, nthmax, npmax+1, nrmax, nsamax) :: moment_vector_obp !by anzai
    real(rkind),dimension(nthmax, npmax+1, nrmax,nthpmax,nsamax) :: pg_tmp !by anzai
    real(rkind),intent(out),dimension(3, nthmax, npmax, nrmax, nsamax) :: moment_vector !by anzai
    !** momentvector(1,:,:,:,:) = pm(:,:)        moment direction element
    !** momentvector(2,:,:,:,:) = pm(:,:)*dthdp  pitch angle direction element
    !** momentvector(3,:,:,:,:) = pm(:,:)*drdp   radial direction element
    real(rkind),allocatable:: U_moment(:,:,:,:,:,:,:,:)
    real(rkind) :: moment_obp_temp, momentl_p, momentl_th, momentl_r
    real(rkind) :: cpitch_ob, thetap_ob, psip_ob
    real(rkind) :: dt
    integer :: nth, np, nr, nsa, nthp, mode(3), nstp, nstpmax, ierr = 0

    !**** initialize
    moment_vector_obp(:,:,:,:,:) = 0.0
    moment_vector(:,:,:,:,:) = 0.0
    moment_obp_temp = 0.0
    momentl_p = 0.0
    momentl_th = 0.0
    momentl_r = 0.0
    dIdu(:,:,:) = 0.0
    ! !** for moment grid
    ! allocate(U_moment(4,4,4,nthmax,nrmax,nthpmax,npmax+1,nsamax))
    ! do nsa = 1, NSAMAX
    !   do nthp = 1, nthpmax
    !     do nr = 1, nrmax
    !       do np = 1, npmax
    !         do nth = 1, nthmax
    !           pg_tmp(nth,np,nr,nthp,nsa) = pg(np,nsa)
    !         end do !** nth
    !       end do !** np
    !     end do !** nr
    !   end do !** nthp
    ! end do !** nsa
    ! do nsa = 1, nsamax
    !   do np = 1, npmax+1
    !   if ( pg(np,nsa) >= 1.d-70 )   call make_U_Dxy(U_moment, pg_tmp, 'p', nsa)
    !   end do !** np
    ! end do !** nsa

    ! !**** Bounce average and making moment vector at moment grid
    ! do nsa = 1, nsamax
    !   do nr = 1, nrmax
    !     do np = 1, npmax+1
    !       do nth = 1, nthmax

    !         if ( np == 1 ) then
    !           cycle
    !         end if

    !         nstpmax = orbit_p(nth,np,nr,nsa)%nstp_max

    !         mode = [0,1,0]
    !         call transformation_matrix(dIdu, orbit_p(nth,np,nr,nsa), nth, np, nr, nsa, mode)

    !         !**** time-integral over an orbit then divide poloidal period
    !         do nstp = 2, nstpmax
    !           ! cpitch_ob = orbit_p(nth,np,nr,nsa)%costh(nstp)
    !           ! psip_ob   = orbit_p(nth,np,nr,nsa)%psip(nstp)
    !           ! thetap_ob = orbit_p(nth,np,nr,nsa)%thetap(nstp)
              
    !           !**** Interpole moment along orbit
    !           ! call interpolate_D_unlessZero(moment_obp_temp, U_moment(:,:,:,:,:,:,np,nsa), &
    !           !                               pg(np,nsa), cpitch_ob, psip_ob, thetap_ob)

    !           dt = orbit_p(nth,np,nr,nsa)%time(nstp)-orbit_p(nth,np,nr,nsa)%time(nstp-1)

    !           !**** Making vector element
    !           moment_vector_obp(1,nth,np,nr,nsa) = moment_vector_obp(1,nth,np,nr,nsa)&
    !                                 + abs(pg(np,nsa)) * dIdu(1,1,nstp) * dt
    !                                 ! + abs(moment_obp_temp) * dIdu(1,1,nstp) * dt
    !           moment_vector_obp(2,nth,np,nr,nsa) = moment_vector_obp(2,nth,np,nr,nsa)&
    !                                 + abs(pg(np,nsa)) * dIdu(1,2,nstp) * dt
    !                                 ! + abs(moment_obp_temp) * dIdu(1,2,nstp) * dt
    !           moment_vector_obp(3,nth,np,nr,nsa) = moment_vector_obp(3,nth,np,nr,nsa)&
    !                                 + abs(pg(np,nsa)) * dIdu(1,3,nstp) * dt
    !                                 ! + abs(moment_obp_temp) * dIdu(1,3,nstp) * dt
    !         end do !** nstp

    !         moment_vector_obp(1,nth,np,nr,nsa) = moment_vector_obp(1,nth,np,nr,nsa) &
    !                                            / orbit_p(nth,np,nr,nsa)%time(nstpmax)
    !         moment_vector_obp(2,nth,np,nr,nsa) = moment_vector_obp(2,nth,np,nr,nsa) &
    !                                            / orbit_p(nth,np,nr,nsa)%time(nstpmax)
    !         moment_vector_obp(3,nth,np,nr,nsa) = moment_vector_obp(3,nth,np,nr,nsa) &
    !                                            / orbit_p(nth,np,nr,nsa)%time(nstpmax)
    !       end do !** nth
    !     end do !** np
    !   end do !** nr
    ! end do !** nsa

    ! do nsa = 1, nsamax
    !   do nr = 1, nrmax
    !     do np = 1, npmax
    !       do nth = 1, nthmax
    !         if (np == npmax + 1) then
    !           momentl_p  = moment_vector_obp(1,nth,np-1,nr,nsa)
    !           momentl_th = moment_vector_obp(2,nth,np-1,nr,nsa)
    !           momentl_r  = moment_vector_obp(3,nth,np-1,nr,nsa)
    !         else if ( np == 1 ) then
    !           momentl_p  = moment_vector_obp(1,nth,np,nr,nsa)
    !           momentl_th = moment_vector_obp(2,nth,np,nr,nsa)
    !           momentl_r  = moment_vector_obp(3,nth,np,nr,nsa)
    !         else
    !           momentl_p  = (moment_vector_obp(1,nth,np,nr,nsa) + moment_vector_obp(1,nth,np+1,nr,nsa)) * 0.50
    !           momentl_th = (moment_vector_obp(2,nth,np,nr,nsa) + moment_vector_obp(2,nth,np+1,nr,nsa)) * 0.50
    !           momentl_r  = (moment_vector_obp(3,nth,np,nr,nsa) + moment_vector_obp(3,nth,np+1,nr,nsa)) * 0.50
    !         end if
 
    !         moment_vector(1,nth,np,nr,nsa) = moment_vector(1,nth,np,nr,nsa) + momentl_p
    !         moment_vector(2,nth,np,nr,nsa) = moment_vector(2,nth,np,nr,nsa) + momentl_th
    !         moment_vector(3,nth,np,nr,nsa) = moment_vector(3,nth,np,nr,nsa) + momentl_r

    !       end do
    !     end do
    !   end do
    ! end do

    !**** Bounce average moment vector
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
         do nth = 1, nthmax
          
            nstpmax = orbit_m(nth,np,nr,nsa)%nstp_max
            call transformation_matrix_for_vector(dIdu, orbit_m(nth,np,nr,nsa), &
                                                   nth, np, nr, nsa)

            !**** for nstp = 1
            !moment_vector(1,nth,np,nr,nsa) = moment_vector(1,nth,np,nr,nsa) + abs(pm(np,nsa))*dIdu(1,1,1)
            !moment_vector(2,nth,np,nr,nsa) = moment_vector(2,nth,np,nr,nsa) + abs(pm(np,nsa))*dIdu(1,2,1)
            !moment_vector(3,nth,np,nr,nsa) = moment_vector(3,nth,np,nr,nsa) + abs(pm(np,nsa))*dIdu(1,3,1)

            do nstp = 2, nstpmax

              dt = orbit_m(nth,np,nr,nsa)%time(nstp)-orbit_m(nth,np,nr,nsa)%time(nstp-1)
              moment_vector(1,nth,np,nr,nsa) = moment_vector(1,nth,np,nr,nsa) & 
                                             + (pm(np,nsa))*dIdu(1,1,nstp)*dt
              moment_vector(2,nth,np,nr,nsa) = moment_vector(2,nth,np,nr,nsa) &
                                             + (pm(np,nsa))*dIdu(1,2,nstp)*dt
              moment_vector(3,nth,np,nr,nsa) = moment_vector(3,nth,np,nr,nsa) &
                                             + (pm(np,nsa))*dIdu(1,3,nstp)*dt
            end do !** nstp

            moment_vector(1,nth,np,nr,nsa) = moment_vector(1,nth,np,nr,nsa)/orbit_m(nth,np,nr,nsa)%time(nstpmax)
            moment_vector(2,nth,np,nr,nsa) = moment_vector(2,nth,np,nr,nsa)/orbit_m(nth,np,nr,nsa)%time(nstpmax)
            moment_vector(3,nth,np,nr,nsa) = moment_vector(3,nth,np,nr,nsa)/orbit_m(nth,np,nr,nsa)%time(nstpmax)

          end do
        end do
      end do
    end do

  end subroutine make_moment_vector

  subroutine transformation_matrix(dIdu, ob, nth, np, nr, nsa, mode)
  !----------------------------------------------------------------
  ! Calcuration of transform matrix  from Energy, momentum space
  ! to momentum, pitch angle, the largest magnetic flux
  !----------------------------------------------------------------
    use fpcomm
    use fowcomm
    use foworbit

    implicit none

    real(rkind),dimension(3,3,max_stp),intent(out) :: dIdu
    type(orbit),intent(in) :: ob
    integer,intent(in) :: nth, np, nr, nsa, mode(3)

    ! elements of transformation matrix, dIdu
    real(rkind) :: dpdp,  dthmdp,  drmdp,&
                   dpdth, dthmdth, drmdth,&
                   dpdr,  dthmdr,  drmdr
    real(rkind) :: pl, Bml, dBmdrl, Fml, dFmdrl, cthm, sthm, dpsimdrl
    real(rkind) :: Bob, Fob, cthob, sthob, dBobdr, dpsiobdr, dFobdr
    real(rkind) :: A(2,2), b(2), detA
    integer :: nstp, nstpmax

    select case(mode(1))
    case(0)
      if ( mode(2) == 0 .and. mode(3) == 0 ) cthm = COS( thetam(nth,np,nr,nsa) )
      if ( mode(2) == 1 .and. mode(3) == 0 ) cthm = COS( thetam_pg(nth,np,nr,nsa) )
      if ( mode(2) == 0 .and. mode(3) == 1 ) cthm = COS( thetam_rg(nth,np,nr,nsa) )
    case(1)
      cthm = COS( thetam_tg(nth,np,nr,nsa) )
    end select
    sthm = SQRT( 1.d0-cthm**2 )

    select case(mode(2))
    case(0)
      pl = pm(np,nsa)*ptfp0(nsa)
    case(1)
      pl = pg(np,nsa)*ptfp0(nsa)
    end select

    select case(mode(3))
    case(0)
      if ( cthm*aefp(nsa) >= 0.d0 ) then
        Bml = Bout(nr)
        dBmdrl = dBoutdr(nr)
      else
        Bml = Bin(nr)
        dBmdrl = dBindr(nr)
      end if
      dpsimdrl = dpsimdr(nr)
      Fml = Fpsi(nr)
      dFmdrl = dFdr(nr)
    case(1)
      if ( cthm*aefp(nsa) >= 0.d0 ) then
        Bml = Bout_rg(nr)
        dBmdrl = dBoutgdr(nr)
      else
        Bml = Bin_rg(nr)
        dBmdrl = dBingdr(nr)
      end if
      dpsimdrl = dpsimgdr(nr)
      Fml = Fpsi_rg(nr)
      dFmdrl = dFgdr(nr)
    end select

    A(1,1) = 2.d0*sthm*cthm/Bml
    A(1,2) = -1.d0*sthm**2*dBmdrl/Bml**2
    A(2,1) = Fml/Bml*pl*sthm
    A(2,2) = aefp(nsa)*dpsimdrl-(dFmdrl*Bml-Fml*dBmdrl)/Bml**2*pl*cthm
    detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)

    nstpmax = ob%nstp_max

    if ( detA /= 0.d0 ) then
      do nstp = 1, nstpmax
        Bob = ob%Babs(nstp)
        Fob = ob%F(nstp)
        cthob = ob%costh(nstp)
        sthob = ob%sinth(nstp)
        dBobdr = ob%dBdr(nstp)
        dpsiobdr = ob%dpsipdr(nstp)
        dFobdr = ob%dFdr(nstp)

        ! dX/dp
        b(1) = 0.d0
        b(2) = Fml/Bml*cthm-Fob/Bob*cthob
        dpdp   = 1.d0
        dthmdp = (A(2,2)*b(1)-A(1,2)*b(2))/detA
        drmdp  = (A(1,1)*b(2)-A(2,1)*b(1))/detA

        ! dX/dtheta
        b(1) = 2.d0*sthob*cthob/Bob
        b(2) = Fob/Bob*pl*sthob
        dpdth   = 0.d0
        dthmdth = (A(2,2)*b(1)-A(1,2)*b(2))/detA
        drmdth  = (A(1,1)*b(2)-A(2,1)*b(1))/detA

        ! dX/drho
        b(1) = -1.d0*sthob**2/Bob**2*dBobdr
        b(2) = aefp(nsa)*dpsiobdr - ( dFobdr*Bob-Fob*dBobdr )/Bob**2*pl*cthob
        dpdr   = 0.d0
        dthmdr = (A(2,2)*b(1)-A(1,2)*b(2))/detA
        drmdr  = (A(1,1)*b(2)-A(2,1)*b(1))/detA

        dIdu(1,1,nstp) = dpdp
        dIdu(1,2,nstp) = dthmdp*ptfp0(nsa)
        dIdu(1,3,nstp) = drmdp*ptfp0(nsa)
        dIdu(2,1,nstp) = dpdth
        dIdu(2,2,nstp) = dthmdth
        dIdu(2,3,nstp) = drmdth
        dIdu(3,1,nstp) = dpdr
        dIdu(3,2,nstp) = dthmdr
        dIdu(3,3,nstp) = drmdr

      end do

    else
      do nstp = 1, nstpmax
        dIdu(1,1,nstp) = 1.d0
        dIdu(1,2,nstp) = 0.d0
        dIdu(1,3,nstp) = 0.d0
        dIdu(2,1,nstp) = 0.d0
        dIdu(2,2,nstp) = 1.d0
        dIdu(2,3,nstp) = 0.d0
        dIdu(3,1,nstp) = 0.d0
        dIdu(3,2,nstp) = 0.d0
        dIdu(3,3,nstp) = 1.d0
      end do

    end if

    ! do nstp = 1, nstpmax
    !   write(6,'(I4)')nstp
    !   write(6,'(3ES12.4)')dIdu(1,1,nstp),dIdu(1,2,nstp),dIdu(1,3,nstp)
    !   write(6,'(3ES12.4)')dIdu(2,1,nstp),dIdu(2,2,nstp),dIdu(2,3,nstp)
    !   write(6,'(3ES12.4)')dIdu(3,1,nstp),dIdu(3,2,nstp),dIdu(3,3,nstp)
    ! end do

  end subroutine transformation_matrix

  subroutine transformation_matrix_for_vector(dIdu, ob, nth, np, nr, nsa)
  !----------------------------------------------------------------
  ! Calcuration of transform matrix  from Energy, momentum space
  ! to momentum, pitch angle, the largest magnetic flux
  !----------------------------------------------------------------
    use fpcomm
    use fowcomm
    ! use foworbit

    implicit none

    real(rkind),dimension(3,3,max_stp),intent(out) :: dIdu
    type(orbit),intent(in) :: ob
    integer,intent(in) :: nth, np, nr, nsa!, mode(3)

    ! elements of transformation matrix, dIdu
    real(rkind) :: dpdp,  dthmdp,  drmdp,&
                   dpdth, dthmdth, drmdth,&
                   dpdr,  dthmdr,  drmdr
    real(rkind) :: pl, Bml, dBmdrl, Fml, dFmdrl, cthm, sthm, dpsimdrl
    real(rkind) :: Bob, Fob, cthob, sthob, dBobdr, dpsiobdr, dFobdr
    real(rkind) :: A(2,2), b(2), detA
    integer :: nstp, nstpmax

    ! select case(mode(1))
    ! case(0)
    !   if ( mode(2) == 0 .and. mode(3) == 0 ) cthm = COS( thetam(nth,np,nr,nsa) )
    !   if ( mode(2) == 1 .and. mode(3) == 0 ) cthm = COS( thetam_pg(nth,np,nr,nsa) )
    !   if ( mode(2) == 0 .and. mode(3) == 1 ) cthm = COS( thetam_rg(nth,np,nr,nsa) )
    ! case(1)
    !   cthm = COS( thetam_tg(nth,np,nr,nsa) )
    ! end select

    !**** for pitch angle
    cthm = COS ( thetam(nth,np,nr,nsa))
    sthm = SQRT( 1.d0-cthm**2 )
    
    ! select case(mode(2))
    ! case(0)
    !   pl = pm(np,nsa)*ptfp0(nsa)
    ! case(1)
    !   pl = pg(np,nsa)*ptfp0(nsa)
    ! end select

    !**** for moment
      pl = pm(np,nsa)*ptfp0(nsa)

    ! select case(mode(3))
    ! case(0)
    !   if ( cthm*aefp(nsa) >= 0.d0 ) then
    !     Bml = Bout(nr)
    !     dBmdrl = dBoutdr(nr)
    !   else
    !     Bml = Bin(nr)
    !     dBmdrl = dBindr(nr)
    !   end if
    !   dpsimdrl = dpsimdr(nr)
    !   Fml = Fpsi(nr)
    !   dFmdrl = dFdr(nr)
    ! case(1)
    !   if ( cthm*aefp(nsa) >= 0.d0 ) then
    !     Bml = Bout_rg(nr)
    !     dBmdrl = dBoutgdr(nr)
    !   else
    !     Bml = Bin_rg(nr)
    !     dBmdrl = dBingdr(nr)
    !   end if
    !   dpsimdrl = dpsimgdr(nr)
    !   Fml = Fpsi_rg(nr)
    !   dFmdrl = dFgdr(nr)
    ! end select

    !**** for radial part
    if ( cthm*aefp(nsa) >= 0.d0 ) then
      Bml = Bout(nr)
      dBmdrl = dBoutdr(nr)
    else
      Bml = Bin(nr)
      dBmdrl = dBindr(nr)
    end if
      dpsimdrl = dpsimdr(nr)
      Fml = Fpsi(nr)
      dFmdrl = dFdr(nr)

    A(1,1) = 2.d0*sthm*cthm/Bml
    A(1,2) = -1.d0*sthm**2*dBmdrl/Bml**2
    A(2,1) = Fml/Bml*pl*sthm
    A(2,2) = aefp(nsa)*dpsimdrl-(dFmdrl*Bml-Fml*dBmdrl)/Bml**2*pl*cthm
    detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)

    nstpmax = ob%nstp_max

    if ( detA /= 0.d0 ) then
      do nstp = 1, nstpmax
        Bob = ob%Babs(nstp)
        Fob = ob%F(nstp)
        cthob = ob%costh(nstp)
        sthob = ob%sinth(nstp)
        dBobdr = ob%dBdr(nstp)
        dpsiobdr = ob%dpsipdr(nstp)
        dFobdr = ob%dFdr(nstp)

        ! dX/dp
        b(1) = 0.d0
        b(2) = Fml/Bml*cthm-Fob/Bob*cthob
        dpdp   = 1.d0
        dthmdp = (A(2,2)*b(1)-A(1,2)*b(2))/detA
        drmdp  = (A(1,1)*b(2)-A(2,1)*b(1))/detA

        ! dX/dtheta
        b(1) = 2.d0*sthob*cthob/Bob
        b(2) = Fob/Bob*pl*sthob
        dpdth   = 0.d0
        dthmdth = (A(2,2)*b(1)-A(1,2)*b(2))/detA
        drmdth  = (A(1,1)*b(2)-A(2,1)*b(1))/detA

        ! dX/drho
        b(1) = -1.d0*sthob**2/Bob**2*dBobdr
        b(2) = aefp(nsa)*dpsiobdr - ( dFobdr*Bob-Fob*dBobdr )/Bob**2*pl*cthob
        dpdr   = 0.d0
        dthmdr = (A(2,2)*b(1)-A(1,2)*b(2))/detA
        drmdr  = (A(1,1)*b(2)-A(2,1)*b(1))/detA

        dIdu(1,1,nstp) = dpdp
        dIdu(1,2,nstp) = dthmdp*ptfp0(nsa)
        dIdu(1,3,nstp) = drmdp*ptfp0(nsa)
        dIdu(2,1,nstp) = dpdth
        dIdu(2,2,nstp) = dthmdth
        dIdu(2,3,nstp) = drmdth
        dIdu(3,1,nstp) = dpdr
        dIdu(3,2,nstp) = dthmdr
        dIdu(3,3,nstp) = drmdr

      end do

    else
      do nstp = 1, nstpmax
        dIdu(1,1,nstp) = 1.d0
        dIdu(1,2,nstp) = 0.d0
        dIdu(1,3,nstp) = 0.d0
        dIdu(2,1,nstp) = 0.d0
        dIdu(2,2,nstp) = 1.d0
        dIdu(2,3,nstp) = 0.d0
        dIdu(3,1,nstp) = 0.d0
        dIdu(3,2,nstp) = 0.d0
        dIdu(3,3,nstp) = 1.d0
      end do

    end if
  end subroutine transformation_matrix_for_vector

  subroutine make_U_Dxy(U_Dxy, Dxyl, x, nsa)
    use fpcomm
    use fowcomm
    USE libspl3d

    implicit none

    real(rkind),intent(out) :: U_Dxy(:,:,:,:,:,:,:,:)
    real(rkind),intent(in) :: Dxyl(:,:,:,:,:)
    character(*),intent(in) :: x
    intent(in) :: nsa
    real(rkind),allocatable :: Dxyl_tmp(:,:,:), U_Dxy_tmp(:,:,:,:,:,:) &
                              , FX(:,:,:), FY(:,:,:), FZ(:,:,:), FXY(:,:,:), FYZ(:,:,:) &
                              , FZX(:,:,:), FXYZ(:,:,:), Xtmp(:)

    integer :: nth, np, nr, nsa, nthp, i, j ,k 
    integer :: nxmax, nymax, nzmax, p, t, ierr = 0

    if ( x == 'p' ) then
      p = 1
      t = 0

      allocate(Xtmp(nthmax))
      do nth = 1, nthmax
        Xtmp(nth) = cosm(nth)
      end do
        
    else if ( x == 't' ) then
      p = 0
      t = 1

      allocate(Xtmp(nthmax+1))
      do nth = 1, nthmax+1
        Xtmp(nth) = cosg(nth)
      end do

    else if ( x == 'm' ) then
      p = 0
      t = 0

      allocate(Xtmp(nthmax))
      do nth = 1, nthmax
        Xtmp(nth) = cosm(nth)
      end do

    end if

    allocate( Dxyl_tmp(nthmax+t,nrmax,nthpmax), U_Dxy_tmp(4,4,4,nthmax+t,nrmax,nthpmax) )
    allocate( FX(nthmax+t,nrmax,nthpmax), FY(nthmax+t,nrmax,nthpmax), FZ(nthmax+t,nrmax,nthpmax), FXY(nthmax+t,nrmax,nthpmax) )
    allocate( FZX(nthmax+t,nrmax,nthpmax), FYZ(nthmax+t,nrmax,nthpmax), FXYZ(nthmax+t,nrmax,nthpmax) )

    do np = 1, npmax+p

      do nthp = 1, nthpmax
        do nr = 1, nrmax
          do nth = 1, nthmax+t
            Dxyl_tmp(nth,nr,nthp) = Dxyl(nth,np,nr,nthp,nsa)
          end do
        end do
      end do

      call SPL3D(Xtmp,psim,theta_p,Dxyl_tmp,FX,FY,FZ,FXY,FYZ,FZX,FXYZ,&
                U_Dxy_tmp,nthmax+t,nrmax,nthmax+t,nrmax,nthpmax,0,0,0,IERR)

      do nthp = 1, nthpmax
        do nr = 1, nrmax
          do nth = 1, nthmax+t

            do k = 1, 4
              do j = 1, 4
                do i = 1, 4
                  U_Dxy(i,j,k,nth,nr,nthp,np,nsa) = U_Dxy_tmp(i,j,k,nth,nr,nthp)
                end do
              end do
            end do

          end do
        end do
      end do

    end do

  end subroutine make_U_Dxy

  subroutine interpolate_D_unlessZero(C_out, U, check0, cpitch_in, psip_in, thetap_in)
    use fpcomm
    use fowcomm
    USE libspl3d
    implicit none

    real(rkind),intent(out) :: C_out
    real(rkind),intent(in) :: U(:,:,:,:,:,:), check0, cpitch_in, psip_in, thetap_in
    integer :: nxmax, nymax, nzmax, ierr
    ierr = 0

    nxmax = size(U,4)
    nymax = size(U,5)
    nzmax = size(U,6)

    if ( check0 < 1.d-70 ) then
      C_out = 0.d0
      return
    end if

    if ( nxmax == nthmax .and. nymax == nrmax ) then
      call SPL3DF(cpitch_in,psip_in,thetap_in,C_out,cosm,psim,theta_p&
                  ,U,nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)

    else if (nxmax == nthmax+1 .and. nymax == nrmax ) then
      call SPL3DF(cpitch_in,psip_in,thetap_in,C_out,cosg,psim,theta_p&
                  ,U,nthmax+1,nrmax,nthmax+1,nrmax,nthpmax,IERR)

    else if (nxmax == nthmax .and. nymax == nrmax+1 ) then
      call SPL3DF(cpitch_in,psip_in,thetap_in,C_out,cosm,psim_rg,theta_p&
                  ,U,nthmax,nrmax+1,nthmax,nrmax+1,nthpmax,IERR)
                                        
    end if


  end subroutine interpolate_D_unlessZero
!===============================================
! Subroutine for factor check
!===============================================
  subroutine check_factorneo(eps_t,cyclo_rho, dTadr,Ta,tau_ele,tau_i, spVI,spVI2, sTI)
  !-------------------------------------
  ! calculate main factors in neoclassical theory
  !--------------------------------------

    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: cyclo_rho, dTadr
    double precision,dimension(nrmax,nsamax),intent(out) :: Ta, spVI,spVI2, sTI
    double precision,dimension(nrmax),intent(out) :: eps_t, tau_ele, tau_i
    double precision Baxis, fact
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom,dfdrhom2, fnsp_l,fnsp_l2
    integer nth,np,nr,nsa,ns

    Baxis = Bin_rg(1)
    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    !****temperature make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
        !Ta(temperature)[keV] rwsl[1/m^3*MJ]
        cyclo_rho(nr,nsa) = sqrt(Ta(nr,nsa)*RKEV/AMFP(nsa)) &
                          * AMFP(nsa)/(AEFP(nsa)*Baxis)
        !****Unit [J]
      end do
    end do

    !****first order derivation
    do nsa = 1, nsamax
      call first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
    !  call first_order_derivative(dndr(:,nsa), rnsl(:, nsa), rm)
    end do


    !**** Calculation of fators
    do nr = 1, nrmax
        eps_t(nr) = rm(nr)*RA/RR
        tau_ele(nr) = fact*sqrt(AMFP(1))*((Ta(nr,1)*RKEV)**1.5d0) &
                    / (rnsl(nr,2)*1.d20*AEFP(2)**2*AEE**2*lnlam(nr,2,1))
        tau_i(nr)   = fact*sqrt(2.d0)*sqrt(AMFP(2)) &
                    * ((Ta(nr,2)*RKEV)**1.5d0) &
                    / (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2))
        !**** *RKEV converts [keV] to [J]
    end do

    !**** initialization ****
    spVI(:,:) = 0.d0

    !********** for new dfdrhom module [2022/3/4]
    fnsp_l(:,:,:,:)=0.d0
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr= 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit
          do nth = 1, nthmax
          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa) !* JIR(nth,np,nr,nsa)
          fnsp_l2(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa)
          end do
        end do
      end do
    end do

    !**** first order derivative
    do nsa = 1, nsamax
      do np = 1, npmax
        do nth = 1, nthmax
          call first_order_derivative(dfdrhom(nth,np,:,nsa), &
                 fnsp_l(nth,np,:,nsa), rm)
          call first_order_derivative(dfdrhom2(nth,np,:,nsa), &
                 fnsp_l2(nth,np,:,nsa), rm)
        end do
      end do
    end do
    !!**** end of new dfdrhom module

    !****Integration over moment(np) and pitch angle(nth)
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
          do nth = 1, nthmax
            !**** dfdrhom is made of fnsp*dVI
            spVI(nr,nsa) = spVI(nr,nsa) &
                          + 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
!                          * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                          /AEE/1.d3 &![keV]
                          *JIR(nth,np,nr,nsa)&
                          !** unit [J][2022/5/29]
                          * delp(ns) * delthm(nth,np,nr,nsa)
            spVI2(nr,nsa) = spVI2(nr,nsa) &
                          + 2.d0/3.d0*dfdrhom2(nth,np,nr,nsa) &
                          * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                          /AEE/1.d3 &![keV]
                          *JIR(nth,np,nr,nsa)&
                          !** unit [J][2022/5/29]
                          * delp(ns) * delthm(nth,np,nr,nsa)
            sTI(nr,nsa) = sTI(nr,nsa) &
                          + 2.d0/3.d0*fnsp_l2(nth,np,nr,nsa) &
                          * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                          /AEE/1.d3 &![keV]
                          *JIR(nth,np,nr,nsa)&
                          !** unit [J][2022/5/29]
                          * delp(ns) * delthm(nth,np,nr,nsa)
          end do
        end do
!        heatfow_out(nr,nsa) = heatfow_out(nr,nsa)/sumVI(nr,nsa)
      end do
    end do

  end subroutine check_factorneo

  subroutine ch_intHD(hfowout_r,hfowout_p,hfowout_t,hfowout_f,&
                                      chi_Dr,chi_Dp, chi_Dt, chi_Fr)
  !--------------------------------------------
  !Subroutine for factor check of FOW heat currents
  !calculate only heat diffusion coefficient elements
  !--------------------------------------------

    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out) :: hfowout_p, &
                                      hfowout_r,hfowout_t, hfowout_f
    double precision,dimension(nrmax,nsamax),intent(out) :: chi_Dp, &
                                      chi_Dt, chi_Dr, chi_Fr
    double precision,dimension(nrmax,nsamax) :: heatfow_l,sumVI_l
    double precision,dimension(nrmax,nsamax) :: Ta, dTadr
    double precision,dimension(nthmax,npmax,nrmax,nsamax) ::fnsp_l, dfdrhom
    double precision,dimension(nthmax,npmax,nrmax,nsamax) ::dfdthm, dfdp
    integer nth,np,nr,nsa,ns

    !**** initialization ****
    hfowout_p(:,:) = 0.d0
    hfowout_r(:,:) = 0.d0
    hfowout_t(:,:) = 0.d0
    hfowout_f(:,:) = 0.d0
    chi_Dp(:,:) = 0.d0
    chi_Dr(:,:) = 0.d0
    chi_Dt(:,:) = 0.d0
    chi_Fr(:,:) = 0.d0
    !****temperature make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3![keV]
        !Ta(temperature)[keV] rwsl[1/m^3*MJ]
      end do
    end do
    !****first order derivation
    do nsa = 1, nsamax
      call first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
    end do

    !********** for new dfdrhom module [2022/3/4]
    fnsp_l(:,:,:,:)=0.d0
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr= 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit
          do nth = 1, nthmax
!          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa) * JIR(nth,np,nr,nsa)
          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa) * 1.d20
          end do
        end do
      end do
    end do

    !**** first order derivative
    do nsa = 1, nsamax
      do np = 1, npmax
        do nth = 1, nthmax
          call first_order_derivative(dfdrhom(nth,np,:,nsa), &
                 fnsp_l(nth,np,:,nsa), rm)
        end do
      end do
    end do
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          call first_order_derivative(dfdthm(:,np,nr,nsa), &
                 fnsp_l(:,np,nr,nsa), thetam(:,np,nr,nsa))
        end do
      end do
    end do
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do nth = 1, nthmax
          call first_order_derivative(dfdp(nth,:,nr,nsa), &
                 fnsp_l(nth,:,nr,nsa), pm(:,nsa))
        end do
      end do
    end do
    !!**** end of new dfdrhom module

    !****Integration over moment(np) and pitch angle(nth)
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
          do nth = 1, nthmax

           !**** for dfdrhom is made of fnsp*dVI
           hfowout_r(nr,nsa) = hfowout_r(nr,nsa) &
                               - (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa)) &
                               /(AEE*1.D3) &
                               * Drrfow(nth,np,nr,nsa) &
                               * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
                               !* JIR(nth,np,nr,nsa) &
                               !* JIR(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
                               !** unit [keV]
           hfowout_p(nr,nsa) = hfowout_p(nr,nsa) &
                               - (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa)) &
                               /(AEE*1.D3) &
                               * Drpfow(nth,np,nr,nsa) &
                               * 2.d0/3.d0*dfdp(nth,np,nr,nsa) &
                               !* JIR(nth,np,nr,nsa) &
                               !* JIR(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
           hfowout_t(nr,nsa) = hfowout_t(nr,nsa) &
                               - (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa)) &
                               /(AEE*1.D3) &
                               * Drtfow(nth,np,nr,nsa) &
                               !* JIR(nth,np,nr,nsa)&
                               !* JIR(nth,np,nr,nsa) &
                               * 2.d0/3.d0*dfdthm(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
           hfowout_f(nr,nsa) = hfowout_f(nr,nsa) &
                               + (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa) )&
                               /(AEE*1.D3) &
                               ! unit converter [J] to [keV] &
                               * Frrfow(nth,np,nr,nsa) &
                               * 2.d0/3.d0*fnsp_l(nth,np,nr,nsa) &
                               !* JIR(nth,np,nr,nsa) &
                               !* JIR(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
          end do
        end do
        chi_Dr(nr,nsa) = hfowout_r(nr,nsa)/dTadr(nr,nsa)
        chi_Dp(nr,nsa) = hfowout_p(nr,nsa)/dTadr(nr,nsa)
        chi_Dt(nr,nsa) = hfowout_t(nr,nsa)/dTadr(nr,nsa)
        chi_Fr(nr,nsa) = hfowout_f(nr,nsa)/dTadr(nr,nsa)
      end do
    end do

  end subroutine ch_intHD

end module check_neoclass