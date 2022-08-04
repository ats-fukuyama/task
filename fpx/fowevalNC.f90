! fowevalNC.f90

MODULE fowevalNC

  PRIVATE

  PUBLIC fow_fout_NC

CONTAINS

  SUBROUTINE fow_gout_NC

    USE fpcomm
    USE fowcomm
    IMPLICIT NONE

    REAL(rkind),DIMENSION(nrmax,nsamax):: &
         Dn_banana,Dn_plateu,chi_banana,chi_plateua, &
         Dn_rw,Dn_rw_bav,chi_rw,chi_rw_bav, &
         Dn_fow,chi_fow
    INTEGER:: nth, np, nr, nsa

    CALL eval_analytic(Dn_banana,Dn_plateu,chi_banana,chi_plateua)
  
    CALL eval_rw(Dn_rw,Dn_rw_bav,chi_rw,chi_rw_bav)

    CALL eval_fow(Dn_fow,chi_fow)

    CALL PAGES
    CALL grd1d(1,rm,Dn_banana(:,1),nrmax,nrmax,1,'@Dn_banana_e vs r@')
    CALL grd1d(2,rm,Dn_banana(:,2),nrmax,nrmax,1,'@Dn_banana_i vs r@')
    CALL grd1d(3,rm,Dn_plateau(:,1),nrmax,nrmax,1,'@Dn_plateau_e vs r@')
    CALL grd1d(4,rm,Dn_plateau(:,2),nrmax,nrmax,1,'@Dn_plateau_i vs r@')
    CALL PAGEE
  
    CALL PAGES
    CALL grd1d(1,rm,Dn_rw(:,1),nrmax,nrmax,1,'@Dn_rw_e vs r@')
    CALL grd1d(2,rm,Dn_rw(:,2),nrmax,nrmax,1,'@Dn_rw_i vs r@')
    CALL grd1d(3,rm,Dn_rw_bav(:,1),nrmax,nrmax,1,'@Dn_rw_bav_e vs r@')
    CALL grd1d(4,rm,Dn_rw_bav(:,2),nrmax,nrmax,1,'@Dn_rw_bav_i vs r@')
    CALL PAGEE
  
    CALL PAGES
    CALL grd1d(1,rm,Dn_fow(:,1),nrmax,nrmax,1,'@Dn_fow_e vs r@')
    CALL grd1d(2,rm,Dn_fow(:,2),nrmax,nrmax,1,'@Dn_fow_i vs r@')
    CALL PAGEE

  END SUBROUTINE fow_gout_NC
  
  SUBROUTINE fow_fout_NC
  !----------------------------------------
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
    double precision,dimension(nrmax,nsamax) :: heatfow, heatrw, heatrwav
    double precision,dimension(nrmax,nsamax) :: chi_a
    double precision,dimension(nrmax,nsamax) :: hfowout_r, hfowout_p
    double precision,dimension(nrmax,nsamax) :: hfowout_t, hfowout_f
    double precision,dimension(nrmax,nsamax) :: chi_Dp,chi_Dt
    double precision,dimension(nrmax,nsamax) :: chi_Dr, chi_Fr
    double precision,dimension(nrmax,nsamax) :: heatba,heatpla
    double precision,dimension(nrmax,nsamax) :: chi_neo_ba,chi_neo_pla
    double precision,dimension(nrmax,nsamax) :: cyclo_rho, dTadr
    double precision,dimension(nrmax,nsamax) :: Ta, spVI, spVI2,sumdf, sTI
    double precision,dimension(nrmax,nsamax) :: nu_eff
    double precision,dimension(nrmax) :: eps_t,tau_ele,tau_i
    integer nth, np, nr, nsa, nsb, mode(3)

    !**** make data folder
    call system('mkdir -p dat')

    !**** calculation of physical values
!    call nu_effective(nu_eff)
    call integral_ParticleDif(Sr_part,Dr)
    call integral_Heatdiff(heatfow,chi_a) ![2022/6/6] editted by anzai
!    call D_random_walk_baverage(Drwav, Drw, Drweff, Drweffav) ![2022/3/7] editted by anzai
    call Particleneoclass(Sr_ba,Sr_pla,Dnba,Dnpla) ![2022/3/7] editted by anzai
    call Heatneoclass(heatba,heatpla,chi_neo_ba,chi_neo_pla) ![2022/3/14] editted by anzai
    call Heat_rw_baverage(heatrwav, heatrw) ![2022/2/23] editted by anzai
    call check_factorneo(eps_t, cyclo_rho, dTadr,Ta,tau_ele,tau_i, spVI,spVI2,sTI) ![2022/2/27] editted by anzai
    call ch_intHD(hfowout_r,hfowout_p,hfowout_t,hfowout_f,chi_Dr,chi_Dp,chi_Dt,chi_Fr) ![2022/6/6]

    !**** write txt file
    call fptxt2D(Dr,"dat/Dr.txt")
    call fptxt2D(Sr_part,"dat/Sr_part.txt")
!    call fptxt2D(Drw,"dat/Drw.txt")
!    call fptxt2D(Drwav,"dat/Drwav.txt")
    call fptxt2D(Dnba,"dat/Dnba.txt")
    call fptxt2D(Dnpla,"dat/Dnpla.txt")
    call fptxt2D(Sr_ba, "dat/Sr_ba.txt") ![2022/7/25] edited by anzai
    call fptxt2D(Sr_pla, "dat/Sr_pla.txt") ![2022/7/25] editted by anzai
    call fptxt2D(heatba, "dat/heatba.txt")![2022/3/4] editted by anzai
    call fptxt2D(heatpla, "dat/heatpla.txt")![2022/3/4] editted by anzai
!    call fptxt2D(Drweff,"dat/Drweff.txt") ![2022/3/2] editted by anzai
!    call fptxt2D(Drweffav,"dat/Drweffav.txt") ![2022/3/2] editted by anzai
    call fptxt2D(heatfow, "dat/heatfow.txt")![2022/2/19] editted by anzai
    call fptxt2D(heatba, "dat/heatba.txt")![2022/3/4] editted by anzai
    call fptxt2D(heatpla, "dat/heatpla.txt")![2022/3/4] editted by anzai
    call fptxt2D(heatrw, "dat/heatrw.txt")![2022/2/23] editted by anzai
    call fptxt2D(heatrwav, "dat/heatrwav.txt")![2022/2/23] editted by anzai

    !***** For factor check
    call fptxt1D(eps_t, "dat/eps_t.txt") ![2022/2/27] editted by anzai
    call fptxt2D(cyclo_rho, "dat/cyclo_rho.txt") ![2022/2/27] editted by anzai
    call fptxt2D(spVI, "dat/spVI.txt") ![2022/3/4] editted by anzai
    call fptxt2D(spVI2, "dat/spVI2.txt") ![2022/3/4] editted by anzai
    call fptxt2D(sumdf, "dat/sumdf.txt") ![2022/3/4] editted by anzai
    call fptxt2D(dTadr, "dat/dTadr.txt") ![2022/2/27] editted by anzai
    call fptxt2D(Ta, "dat/Ta.txt") ![2022/2/27] editted by anzai
    call fptxt2D(sTI, "dat/sTI.txt") ![2022/6/10] editted by anzai
    call fptxt1D(tau_ele, "dat/tau_ele.txt") ![2022/2/27] editted by anzai
    call fptxt1D(tau_i, "dat/tau_i.txt") ![2022/2/27] editted by anzai
    call fptxt2D(nu_eff, "dat/nu_eff.txt") ![2022/3/24] editted by anzai
    call fptxt2D(hfowout_r, "dat/hfow_r.txt")![2022/6/6] editted by anzai
    call fptxt2D(hfowout_p, "dat/hfow_p.txt")![2022/6/6] editted by anzai
    call fptxt2D(hfowout_t, "dat/hfow_t.txt")![2022/6/6] editted by anzai
    call fptxt2D(hfowout_f, "dat/hfow_f.txt")![2022/6/6] editted by anzai
    call fptxt2D(chi_a, "dat/chi_a.txt")![2022/2/19] editted by anzai
    call fptxt2D(chi_Dr, "dat/chi_Dr.txt")![2022/2/19] editted by anzai
    call fptxt2D(chi_Dp, "dat/chi_Dp.txt")![2022/2/19] editted by anzai
    call fptxt2D(chi_Dt, "dat/chi_Dt.txt")![2022/2/19] editted by anzai
    call fptxt2D(chi_Fr, "dat/chi_Fr.txt")![2022/2/19] editted by anzai
    call fptxt2D(chi_neo_ba, "dat/chi_neo_ba.txt")![2022/2/19] editted by anzai
    call fptxt2D(chi_neo_pla, "dat/chi_neo_pla.txt")![2022/2/19] editted by anzai

  END SUBROUTINE fow_fout_NC
  

  subroutine Particleneoclass(Sr_ba,Sr_pla,Dnba,Dnpla)
    !------------------------------------------------------------
    ! Subroutine for neoclassical particle flux in banana region
    ! Calculate particle flux and diffusion coefficient
    ! Formulae from "Tokamaks fourth edition"
    !------------------------------------------------------------

    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: Dnba, Dnpla
    double precision,dimension(nrmax,nsamax),intent(out)  :: Sr_ba, Sr_pla
        double precision,dimension(nrmax,nsamax)  :: Ta,dTadr,dndr
    double precision fact,fact_ba,fact_pla,tau_ele,rho_e,eps_t,Baxis,B_p
    integer nr, nsa

    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bing(1) ! approximation on B by Baxis

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
        tau_ele = fact*sqrt(AMFP(1))*((Ta(nr,1))**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**2*lnlam(nr,2,1)*AEE**2)
        rho_e = sqrt(2*Ta(nr,1)/AMFP(1))*AMFP(1)/(AEE*Baxis)
        eps_t = rm(nr)*RA/RR
        B_p = rm(nr)*RA*BB/(safety_factor(nr)*RR)
        fact_ba = (safety_factor(nr)**2)*1.d20*rnsl(nr,nsa) &
             *(rho_e**2)/((eps_t**1.5)*tau_ele)*RKEV ! modified by anzai[*RKEV]
        fact_pla = -sqrt(pi)*(eps_t**2)*Ta(nr,1)*RKEV &! modified by anzai [*RKEV]
             *rho_e*1.d20*rnsl(nr,nsa)/(2*AEE*B_p*rm(nr)*RA)

        Sr_ba(nr,nsa) = fact_ba*(-1.22d0*(1+Ta(nr,2)/Ta(nr,1)) &
                      * 1.d20 * dndr(nr, nsa) &
                      + 4.3d-1*dTadr(nr,1)/Ta(nr,1) &
                      + 1.9d-1*dTadr(nr,2)/Ta(nr,1) )
        !**** neglect E_para term
        Sr_pla(nr,nsa) = fact_pla*((1+Ta(nr,2)/Ta(nr,1)) &
                       * 1.d20 * dndr(nr,nsa) &
                       / (1.d20*rnsl(nr,nsa)) &
                       + 1.5d0*dTadr(nr,1)/Ta(nr,1) + 1.5d0*dTadr(nr,2)/Ta(nr,1) )
        !**** neglect B term

        Dnba(nr,nsa) = - Sr_ba(nr, nsa)/(1.d20*dndr(nr, nsa))
        Dnpla(nr,nsa) = - Sr_pla(nr, nsa)/(1.d20*dndr(nr, nsa))
     end do
   end do

  end subroutine Particleneoclass

  subroutine Heatneoclass(heatba,heatpla,chi_neo_ba,chi_neo_pla)
  !---------------------------------------------------------------
  !subroutine for neoclassical heat flux in banana, plateau region
  !calculate heat flux and diffusion coefficient
  ! main calcul using[J]
  ! Formulae from "Tokamaks fourth edition"
  !--------------------------------------------------------------

    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: heatba, heatpla
    double precision,dimension(nrmax,nsamax),intent(out)  :: chi_neo_ba, chi_neo_pla
    double precision,dimension(nrmax,nsamax) :: Sr_ba, Sr_pla
    double precision,dimension(nrmax,nsamax) :: Ta, dTadr, dndr
    double precision, dimension(nsamax) :: rho_a
    double precision fact, fact_s, tau_ele,tau_i, eps_t, Baxis, B_p
    integer nr, nsa

    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bing(1) ! approximation on B by Baxis

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

        tau_ele = fact*sqrt(AMFP(1))*((Ta(nr,1))**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**2*AEE**2*lnlam(nr,2,1))
        tau_i   = fact*sqrt(2.d0)*sqrt(AMFP(2))*((Ta(nr,2))**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2))
        !**** *RKEV converts [keV] to [J]

        rho_a(nsa) = sqrt(Ta(nr,nsa)/AMFP(nsa))*AMFP(nsa)/(AEFP(nsa)*Baxis)
        !**** *RKEV converts [J] to [keV]
        eps_t = rm(nr)*RA/RR
        B_p = rm(nr)*RA*BB/(safety_factor(nr)*RR)
        !****calculate heat flux
        if (nsa == 1) then
         fact_s = (safety_factor(nr)**2) &
             *(rho_a(nsa)**2)/((eps_t**1.5)*tau_ele)
         Sr_ba(nr,nsa) = fact_s*rnsl(nr,nsa)*1.D20*Ta(nr,1) &
            * ((1.53d0*(1+Ta(nr,2)/Ta(nr,1))*dndr(nr, nsa))/rnsl(nr,nsa) &
              - 1.81d0*dTadr(nr,1)/Ta(nr,1) &
              - 0.27d0*dTadr(nr,2)/Ta(nr,1))/RKEV!test
        ! neglect E_para term

         Sr_pla(nr,nsa) = - sqrt(PI)/4*eps_t**2*Ta(nr,nsa) &
               / (AEE*B_p)*rho_a(nsa)*rnsl(nr,nsa)*1.d20*Ta(nr,nsa) &
               / rm(nr)*( (1+Ta(nr,2)/Ta(nr,1))*dndr(nr,nsa) &
               / (rnsl(nr,nsa)*1.d20) &
               + 7.5d0*dTadr(nr,nsa)/Ta(nr,nsa) &
               + 1.5d0*dTadr(nr,2)/Ta(nr,nsa))/RKEV!test
        ! neglect E_para term
        else
         fact_s = (safety_factor(nr)**2) &
              *(rho_a(nsa)**2)/((eps_t**1.5)*tau_i)
         Sr_ba(nr,nsa) = - 0.68d0*fact_s*(1 + 0.48d0*sqrt(eps_t)) &
              *rnsl(nr,nsa)*1.d20*dTadr(nr,2)/RKEV !test
         Sr_pla(nr,nsa) = - 1.5d0*sqrt(pi) &
                        * eps_t**2*Ta(nr,nsa)/(AEE*B_p) &
                        * rho_a(nsa)/rm(nr)*rnsl(nr,nsa) &
                        * 1.d20*dTadr(nr,nsa)/RKEV !test
        end if

        heatba(nr,nsa)  = Sr_ba(nr, nsa)
        heatpla(nr,nsa) = Sr_pla(nr, nsa)
        chi_neo_ba      = -Sr_ba(nr,nsa)/(dTadr(nr,nsa)/RKEV)
        chi_neo_pla     = -Sr_pla(nr,nsa)/(dTadr(nr,nsa)/RKEV)
      end do
    end do

  end subroutine Heatneoclass


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

  end subroutine bounce_average_for_Drw

  subroutine D_random_walk_baverage(Drwav, Drw, Drweff, Drweffav)
  !-----------------------------------------------------
  ! For particle diffusion in banana region
  ! modified and added by anzai [2022/3/2]
  !-----------------------------------------------------
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: Drw, Drwav
    double precision,dimension(nrmax,nsamax),intent(out)  :: Drweff, Drweffav
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: fnsp_l
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drweff_l, Drweff_ba
    double precision,dimension(npmax,nrmax,nsamax) :: nu
    double precision,dimension(nrmax,nsamax) :: Ta
    double precision,dimension(nrmax,nsamax) :: nu_eff
    double precision step_len, step_time, step_time_eff, eps_t
    double precision rho_theta, Baxis, pv, B_theta, rho
    double precision dVI, sumVI
    integer nth, np, nr, nsa, ns

!     !**** initialization ****
!     Baxis = Bing(1)
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
!             dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JI(nth,np,nr,nsa)
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

  end subroutine D_random_walk_baverage

  subroutine Heat_rw_baverage(heatrwav, heatrw)
  !---------------------------------------------
  ! heat diffusion coef calcul module
  ! This module corresponds to integral_Heatdiff
  !---------------------------------------------

    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out) :: heatrw, heatrwav
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom, fnsp_l
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
    double precision,dimension(nrmax,nsamax) :: nu, Ta, Drw, vth
    double precision step_len, step_time, eps_t, rho_theta, Baxis, pv, rho
    double precision dVI, sumVI
    integer nth, np, nr, nsa, ns

    !**** initialization ****
    Baxis = Bing(1)
    heatrw(:,:) = 0.d0
    heatrwav(:,:) = 0.d0
    sumVI = 0.d0
 !**** for new dfdrhom module [2022/3/4]
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
           call first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp_l(nth,np,:    ,nsa), rm)
         end do
       end do
     end do
 !**** end of new dfdrhom module

    !**** Temperature and flecency make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)
        !**** Ta(temperature)[J]
        ! /RKEV = 1/(AEE*1.D3) unit converter [keV]to [J]
      end do
    end do

   call nu_effective(nu)

    !**** Calculation of particle diff coef
    do nsa = 1, nsamax
      do nr = 1, nrmax
        eps_t = rm(nr)*RA/RR
        do np = 1, npmax
          step_time = 1.d0/nu(nr,nsa)
          !** thermal rho
!          vth(nr,nsa) = SQRT( 2.d0*Ta(nr,nsa)/AMFP(nsa)*RKEV)
!          rho = vth(nr,nsa)*AMFP(nsa)/(aefp(nsa)*Baxis)
          !** original rho
          rho = (pm(np,nsa)*ptfp0(nsa) ) &!/sqrt(AEE*1.d3)) &
              / (aefp(nsa)*Baxis)
          !****Unit [J]

          do nth = 1, nthmax
            step_len = rho * safety_factor(nr)/sqrt(eps_t)
            Drwlocal(nth,np,nr,nsa) = step_len**2/step_time*sqrt(eps_t)
          end do
        end do
      end do
    end do

    call bounce_average_for_Drw(Drwba, Drwlocal)

    !**** Heat diffusion coef calcul
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            !**** Heat flux for FOW momentum not heat velocity [2022/3/22]
            heatrw(nr,nsa) = heatrw(nr,nsa)&
                           - (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                           * Drwlocal(nth,np,nr,nsa) &
                           * JI(nth,np,nr,nsa) &
                           * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa)*1.d20 &
                           / (AEE*1.D3)&  !** Unit convert [J] to [keV]
                           * delp(ns)*delthm(nth,np,nr,nsa)
            heatrwav(nr,nsa) = heatrwav(nr,nsa) &
                             - (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                             * Drwba(nth,np,nr,nsa) &
                             * JI(nth,np,nr,nsa) &
                             * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa)*1.d20 &
                             / (AEE*1.D3)& !** Unit convert [J] to [keV]
                             *delp(ns)* delthm(nth,np,nr,nsa)

!            sumVI = sumVI - dfdrhom(nth,np,nr,nsa)&
!                  * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa))/(AEE*1.D3)
!*** Heat flux of heat velocity [2022/3/22]
!            heatrw(nr,nsa) = heatrw(nr,nsa)&
!                           + (vth(nr,nsa)*AMFP(nsa))**2/(2*AMFP(nsa)) &
!                           * Drwlocal(nth,np,nr,nsa) &
!                           * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
!                           / (AEE*1.D3)  !**** Unit convert [J] to [keV]
!            heatrwav(nr,nsa) = heatrwav(nr,nsa) &
!                             + (vth(nr,nsa)*AMFP(nsa))**2/(2*AMFP(nsa)) &
!                             * Drwba(nth,np,nr,nsa) &
!                             * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
!                             / (AEE*1.D3)  !****Unit convert [J] to [keV]
!
          end do
        end do
!        heatrw(nr, nsa)  = heatrw(nr, nsa)/sumVI
!        heatrwav(nr,nsa) = heatrwav(nr,nsa)/sumVI
      end do
    end do

  end subroutine Heat_rw_baverage

!===================================================
! Subroutines for particle diffusion coef
!===================================================
  subroutine integral_ParticleDif(Sr_part,Dr)
  !---------------------------------------------
  ! calculation of FOW diff coef
  !---------------------------------------------
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out) :: Sr_part, Dr
    double precision,dimension(nrmax,nsamax) :: Na, dNadr
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom,dfdp
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdthm
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: fnsp_l
    double precision Drr_, dVI
    integer nth,np,nr,nsa,ns

    !**** initialization ****
    Sr_part(:,:) = 0.d0
    Dr(:,:) = 0.d0
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
          call first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp_l(nth,np,:,nsa), rm)
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
    !**** end of new dfdrhom module

    !**** Source calculation
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit
          do nth = 1, nthmax
            Sr_part(nr,nsa) = Sr_part(nr,nsa) &
                            - 1.d0&
                            * (Drrfow(nth,np,nr,nsa)&
                            * dfdrhom(nth,np,nr,nsa)*1.d20 &
                            + Drpfow(nth,np,nr,nsa) &
                            * dfdp(nth,np,nr,nsa)*1.d20 &
                            + Drtfow(nth,np,nr,nsa) &
                            * dfdthm(nth,np,nr,nsa)*1.d20&
                            ) &
                            * JI(nth,np,nr,nsa) &
                            !* JI(nth,np,nr,nsa) &
                            * delp(ns) * delthm(nth,np,nr,nsa) &

                          + 1.d0&
                          * Frrfow(nth,np,nr,nsa) &
                          * fnsp_l(nth,np,nr,nsa)*1.d20 &
                          * JI(nth,np,nr,nsa) &
                          !* JI(nth,np,nr,nsa) &
                          * delp(ns) * delthm(nth,np,nr,nsa)
          end do
        end do
        Dr(nr,nsa) = -Sr_part(nr,nsa)/dNadr(nr,nsa)
      end do
    end do

  end subroutine integral_ParticleDif

!======================================================================
!subroutines for heat fluxes
!======================================================================

  subroutine integral_Heatdiff(heatfow_out,chi_a)
  !--------------------------------------------
  !Subroutine for FOW heat diffuion coefficient
  !calculate only heat diffusion coefficient
  !--------------------------------------------

    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out) :: heatfow_out
    double precision,dimension(nrmax,nsamax),intent(out) :: chi_a
    double precision,dimension(nrmax,nsamax) :: heatfow_l,sumVI_l,sumVI
    double precision,dimension(nrmax,nsamax) :: Ta, dTadr
    double precision,dimension(nthmax,npmax,nrmax,nsamax) ::fnsp_l, dfdrhom
    double precision,dimension(nthmax,npmax,nrmax,nsamax) ::dfdthm, dfdp
    integer nth,np,nr,nsa,ns

    !**** initialization ****
    heatfow_out(:,:) = 0.d0
    sumVI_l(:,:) = 0.d0
    fnsp_l(:,:,:,:)=0.d0

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

    !****Integration over moment(np) and pitch angle(nth)
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
          do nth = 1, nthmax

           !**** for dfdrhom is made of fnsp*dVI
           heatfow_out(nr,nsa) = heatfow_out(nr,nsa) &
                               - (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa)) &
                               / (AEE*1.D3)*2.d0/3.d0 * 1.d0&
                               *(Drrfow(nth,np,nr,nsa)&
                               * dfdrhom(nth,np,nr,nsa)*1.d20 &
                               + Drpfow(nth,np,nr,nsa) &
                               * dfdp(nth,np,nr,nsa)*1.d20 &
                               + Drtfow(nth,np,nr,nsa) &
                               * dfdthm(nth,np,nr,nsa)*1.d20&
                               ) &
                               * JI(nth,np,nr,nsa) &
                               !* JI(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa) &

                               + (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa) )&
                               / (AEE*1.D3)*2.d0/3.d0 *1.d0&
                               ! unit converter [J] to [keV] &
                               * Frrfow(nth,np,nr,nsa) &
                               * fnsp_l(nth,np,nr,nsa)*1.d20 &
                               * JI(nth,np,nr,nsa) &
                               !* JI(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
                               !** unit [keV] [22/6/6]
          end do
        end do
        chi_a(nr,nsa) = -heatfow_out(nr,nsa)/dTadr(nr,nsa)
      end do
    end do

  end subroutine integral_Heatdiff

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

    Baxis = Bing(1)
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
          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa) !* JI(nth,np,nr,nsa)
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
                          *JI(nth,np,nr,nsa)&
                          !** unit [J][2022/5/29]
                          * delp(ns) * delthm(nth,np,nr,nsa)
            spVI2(nr,nsa) = spVI2(nr,nsa) &
                          + 2.d0/3.d0*dfdrhom2(nth,np,nr,nsa) &
                          * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                          /AEE/1.d3 &![keV]
                          *JI(nth,np,nr,nsa)&
                          !** unit [J][2022/5/29]
                          * delp(ns) * delthm(nth,np,nr,nsa)
            sTI(nr,nsa) = sTI(nr,nsa) &
                          + 2.d0/3.d0*fnsp_l2(nth,np,nr,nsa) &
                          * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                          /AEE/1.d3 &![keV]
                          *JI(nth,np,nr,nsa)&
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
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
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
!          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa) * JI(nth,np,nr,nsa)
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
                               * JI(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
                               !** unit [keV]
           hfowout_p(nr,nsa) = hfowout_p(nr,nsa) &
                               - (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa)) &
                               /(AEE*1.D3) &
                               * Drpfow(nth,np,nr,nsa) &
                               * 2.d0/3.d0*dfdp(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
           hfowout_t(nr,nsa) = hfowout_t(nr,nsa) &
                               - (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa)) &
                               /(AEE*1.D3) &
                               * Drtfow(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa)&
                               * JI(nth,np,nr,nsa) &
                               * 2.d0/3.d0*dfdthm(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
           hfowout_f(nr,nsa) = hfowout_f(nr,nsa) &
                               + (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa) )&
                               /(AEE*1.D3) &
                               ! unit converter [J] to [keV] &
                               * Frrfow(nth,np,nr,nsa) &
                               * 2.d0/3.d0*fnsp_l(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
          end do
        end do
        chi_Dr(nr,nsa) = hfowout_r(nr,nsa)!/dTadr(nr,nsa)
        chi_Dp(nr,nsa) = hfowout_p(nr,nsa)!/dTadr(nr,nsa)
        chi_Dt(nr,nsa) = hfowout_t(nr,nsa)!/dTadr(nr,nsa)
        chi_Fr(nr,nsa) = hfowout_f(nr,nsa)!/dTadr(nr,nsa)
      end do
    end do

  end subroutine ch_intHD

  ! *** Modules of collision flecuencies ***

  SUBROUTINE cal_nu_ei(nu)
    USE fpcomm
    USE fowcomm
    IMPLICIT NONE
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

  END SUBROUTINE cal_nu_ei

  SUBROUTINE nu_effective(nu_eff)
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

  ! *** Utilities ***

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

end module fowevalNC
