! check_neoclass.f90
! [2022/3/7]
! ***************************************
!  Calculation of diffusion coefficients
! ***************************************
! made by ota / modified by anzai
!

module check_neoclass
  private

  public :: output_neoclass

contains

!===================================================
! COMMON modules
!===================================================
  subroutine output_neoclass
  !----------------------------------------
  ! module for file out put as *.txt
  !----------------------------------------
    use fpcomm
    use fowcomm
    use fpwrite

    implicit none
    double precision,dimension(nrmax,nsamax) :: Drhmrhm
    double precision,dimension(nrmax,nsamax) :: Drw,Drwav,Drweff,Drweffav
    double precision,dimension(nrmax,nsamax) :: Dnba,Dnpla
    double precision,dimension(nrmax,nsamax) :: heatfow, heatrw, heatrwav
    double precision,dimension(nrmax,nsamax) :: heatba,heatpla
    double precision,dimension(nrmax,nsamax) :: cyclo_rho, dTadr
    double precision,dimension(nrmax,nsamax) :: Ta, sumVI, sumdf
    double precision,dimension(nrmax) :: eps_t,tau_ele,tau_i
    integer nth, np, nr, nsa, nsb, mode(3)
   
    !**** make data folder
    call system('mkdir -p dat')

    !**** calculation of physical values
    call integral_Drr(Drhmrhm)
    call integral_Heatdiff(heatfow,sumVI) ![2022/3/4] editted by anzai
    call D_random_walk_baverage(Drwav, Drw, Drweff, Drweffav) ![2022/3/7] editted by anzai
    call Particleneoclass(Dnba,Dnpla) ![2022/3/7] editted by anzai
    call Heatneoclass(heatba,heatpla) ![2022/3/14] editted by anzai
    call Heat_rw_baverage(heatrwav, heatrw) ![2022/2/23] editted by anzai
    call check_factorneo(eps_t, cyclo_rho, dTadr,Ta,tau_ele,tau_i) ![2022/2/27] editted by anzai

    !**** write txt file
    call fptxt2D(Drhmrhm,"dat/Drhmrhm.txt")
    call fptxt2D(Drw,"dat/Drw.txt")
    call fptxt2D(Drwav,"dat/Drwav.txt")
    call fptxt2D(Dnba, "dat/Dnba.txt") ![2022/1/31] edited by anzai
    call fptxt2D(Dnpla, "dat/Dnpla.txt") ![2022/1/31] editted by anzai
    call fptxt2D(Drweff,"dat/Drweff.txt") ![2022/3/2] editted by anzai
    call fptxt2D(Drweffav,"dat/Drweffav.txt") ![2022/3/2] editted by anzai
    call fptxt2D(heatfow, "dat/heatfow.txt")![2022/2/19] editted by anzai
    call fptxt2D(heatba, "dat/heatba.txt")![2022/3/4] editted by anzai
    call fptxt2D(heatpla, "dat/heatpla.txt")![2022/3/4] editted by anzai
    call fptxt2D(heatrw, "dat/heatrw.txt")![2022/2/23] editted by anzai
    call fptxt2D(heatrwav, "dat/heatrwav.txt")![2022/2/23] editted by anzai

    !***** For factor check
    call fptxt1D(eps_t, "dat/eps_t.txt") ![2022/2/27] editted by anzai
    call fptxt2D(cyclo_rho, "dat/cyclo_rho.txt") ![2022/2/27] editted by anzai
    call fptxt2D(sumVI, "dat/sumVI.txt") ![2022/3/4] editted by anzai
    call fptxt2D(sumdf, "dat/sumdf.txt") ![2022/3/4] editted by anzai
    call fptxt2D(dTadr, "dat/dTadr.txt") ![2022/2/27] editted by anzai
    call fptxt2D(Ta, "dat/Ta.txt") ![2022/2/27] editted by anzai
    call fptxt1D(tau_ele, "dat/tau_ele.txt") ![2022/2/27] editted by anzai
    call fptxt1D(tau_i, "dat/tau_i.txt") ![2022/2/27] editted by anzai

  end subroutine output_neoclass

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
        nu(nr,nsa) = 1.d0/tau
!        WRITE(6,'(A,I6,6ES12.4)') 'nu:',nr,Te,lnlam(nr,2,1),tau,nu(nr,nsa)
      end do
    end do
    
  end subroutine

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

!            if() then
!             
!            else if
              nudb = fact*nuBra*(fx-gx)/x**3

            nud(np,nr,nsa) = nud(np,nr,nsa)+nudb

          end do
        end do
      end do
    end do

  end subroutine nu_deflection

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


  end subroutine

!===================================================
! Subroutines for particle diffusion coef
!===================================================
  subroutine integral_Drr(Drr_out)
  !---------------------------------------------
  ! calculation of FOW diff coef
  !---------------------------------------------
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out) :: Drr_out
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom
    double precision Drr_, dVI, sumVI, JIl_p, JIl_m
    integer nth,np,nr,nsa,ns
    
    !**** initialization ****
    Drr_out(:,:) = 0.d0
    sumVI = 0.d0
  
    !**** first order derivatives
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

  subroutine D_random_walk_baverage(Drwav, Drw, Drweff, Drweffav)
  !-----------------------------------------------------
  ! For particle diffusion in banana region
  !modified and added by anzai [2022/3/2]
  !-----------------------------------------------------
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: Drw, Drwav
    double precision,dimension(nrmax,nsamax),intent(out)  :: Drweff, Drweffav
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drweff_l, Drweff_ba
    double precision,dimension(npmax,nrmax,nsamax) :: nu
    double precision,dimension(nrmax,nsamax) :: Ta
    double precision,dimension(nrmax,nsamax) :: nu_eff
    double precision step_len, step_time, step_time_eff, eps_t, rho_theta, Baxis, pv, B_theta, rho
    double precision dVI, sumVI
    integer nth, np, nr, nsa

    !**** initialization ****
    Baxis = Bing(1)
    Drw(:,:) = 0.d0
    Drwav(:,:) = 0.d0
    sumVI = 0.d0

    call nu_deflection(nu)
 !   call ca_nu_ei(nu)

    do nsa = 1, nsamax
      do np = 1, npmax
        do nth = 1, nthmax
           call first_order_derivative(dfdrhom(nth,np,:,nsa), &
                fnsp(nth,np,:,nsa), rm)
        end do
      end do
    end do
    
    !************************by anzai[2022/2/15]
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/( 1.5d0*rnsl(nr,nsa)*1.d20)  / AEE/1.D3
        !** Ta(temperature)[keV]
        nu_eff(nr,nsa) = (1.D0/6.4D14)*rnsl(nr,nsa)*1.D20/(Ta(nr,nsa)**1.5D0)
        nu_eff(nr,nsa) = nu_eff(nr, nsa)/(rm(nr)*RA/RR)*sqrt(AMFP(1)/AMFP(nsa)) 
        !** Effective nu
      end do
    end do
    
    do nsa = 1, nsamax 
      do nr = 1, nrmax
        eps_t = rm(nr)*RA/RR
        ! B_theta = safety_factor(nr)/RR*dpsimdr(nr)/RA
        do np = 1, npmax
          
          step_time = 1.d0/nu(np, nr, nsa)
          step_time_eff = 1.d0/nu_eff(nr, nsa)
          rho = pm(np,nsa)*ptfp0(nsa)/(aefp(nsa)*Baxis) 
          !** [keV]verocity's length

          do nth = 1, nthmax
            !step_len = rho
            step_len = rho * safety_factor(nr)/sqrt(eps_t)
            !Drwlocal(nth,np,nr,nsa) = step_len**2/step_time*1.46d0*sqrt(eps_t)
            Drwlocal(nth,np,nr,nsa) = step_len**2/step_time*sqrt(eps_t)
            Drweff_l(nth,np,nr,nsa) = step_len**2/step_time_eff*sqrt(eps_t)
          end do
        end do
      end do
    end do

    call bounce_average_for_Drw(Drwba, Drwlocal)
    call bounce_average_for_Drw(Drweff_ba, Drweff_l)

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)

            Drw(nr,nsa) = Drw(nr,nsa)&
                          + Drwlocal(nth,np,nr,nsa)*dfdrhom(nth,np,nr,nsa)*dVI
            Drwav(nr,nsa) = Drwav(nr,nsa)&
                          + Drwba(nth,np,nr,nsa)*dfdrhom(nth,np,nr,nsa)*dVI
            Drweff(nr,nsa) = Drweff(nr,nsa)&
                          + Drweff_l(nth,np,nr,nsa)*dfdrhom(nth,np,nr,nsa)*dVI
            Drweffav(nr,nsa) = Drweffav(nr,nsa)&
                          + Drweff_ba(nth,np,nr,nsa)*dfdrhom(nth,np,nr,nsa)*dVI
            sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
          end do
        end do
        Drw(nr,nsa) = Drw(nr,nsa)/sumVI
        Drwav(nr,nsa) = Drwav(nr,nsa)/sumVI
        Drweff(nr,nsa) = Drweff(nr,nsa)/sumVI
        Drweffav(nr,nsa) = Drweffav(nr,nsa)/sumVI
      end do
    end do

  end subroutine D_random_walk_baverage

  subroutine Particleneoclass(Dnba,Dnpla)
    !----------------------------------
    !subroutine for neoclassical particle flux in banana region
    !calculate particle flux and diffusion coefficient
    !----------------------------------
    
    use fpcomm
    use fowcomm
    
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: Dnba, Dnpla 
    double precision,dimension(nrmax,nsamax) :: Sr_ba, Sr_pla,Ta,dTadr,dndr
    double precision fact,fact_ba,fact_pla,tau_ele,rho_e,eps_t,Baxis,B_p
    integer nr, nsa

    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bing(1) ! approximation on B by Baxis
 
    !****temperature make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6*RKEV/(1.5d0*rnsl(nr,nsa)*1.d20)
        !**** Ta [keV]
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
        !B_p = dpsimdr(nr)*RA/(safety_factor(nr)*RR) ! careful
        B_p = rm(nr)*RA*BB/(safety_factor(nr)*RR) 
        fact_ba = (safety_factor(nr)**2)*1.d20*rnsl(nr,nsa) &
             *(rho_e**2)/((eps_t**1.5)*tau_ele)
        fact_pla = -sqrt(pi)*(eps_t**2)*Ta(nr,1) &
             *rho_e*1.d20*rnsl(nr,nsa)/(2*AEE*B_p*rm(nr)*RA)        
       
        Sr_ba(nr,nsa) = fact_ba*(-1.22d0*(1+Ta(nr,2)/Ta(nr,1)) &
                      * 1.d20 * dndr(nr, nsa)) 
            ! + 4.3d-1*dTadr(nr,1)/Ta(nr,1) &
            ! + 1.9d-1*dTadr(nr,2)/Ta(nr,1) ) 
        !**** neglect E_para term
        Sr_pla(nr,nsa) = fact_pla*((1+Ta(nr,2)/Ta(nr,1)) & 
                       * 1.d20 * dndr(nr,nsa) &
                       / (1.d20*rnsl(nr,nsa)))! &
            ! + 1.5d0*dTadr(nr,1)/Ta(nr,1) + 1.5d0*dTadr(nr,2)/Ta(nr,1) )
        !**** neglect B term
   
        Dnba(nr,nsa) = - Sr_ba(nr, nsa)/(1.d20*dndr(nr, nsa))
        Dnpla(nr,nsa) = - Sr_pla(nr, nsa)/(1.d20*dndr(nr, nsa))
     end do
   end do

  end subroutine Particleneoclass
  
!======================================================================
!subroutines for heat fluxes
!======================================================================

  subroutine integral_Heatdiff(heatfow_out,sumVI)
  !--------------------------------------------
  !Subroutine for FOW heat diffuion coefficient
  !calculate only heat diffusion coefficient
  !--------------------------------------------
    
    use fpcomm
    use fowcomm
    
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out) :: heatfow_out, sumVI
    double precision,dimension(nrmax,nsamax) :: heatfow_l,sumVI_l, sumdf
    double precision,dimension(nthmax,npmax,nrmax,nsamax) ::fnsp_l, dfdrhom
    double precision Drr_, dVI, JIl_p, JIl_m
!    double precision sumVI
    integer nth,np,nr,nsa,ns
  
    !**** initialization ****
    heatfow_out(:,:) = 0.d0
    sumVI_l(:,:) = 0.d0
    sumdf(:,:) = 0.d0

!********** for new dfdrhom module [2022/3/4]
    fnsp_l(:,:,:,:)=0.d0
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr= 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit
          do nth = 1, nthmax
          dVI = delp(ns)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)
          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa) * dVI
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
!!**** end of new dfdrhom module


!************ original derivative
!    !****first order derivative
!    do nsa = 1, nsamax
!      do np = 1, npmax
!        do nth = 1, nthmax
!          call first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
!       end do
!      end do
!    end do

    !****Integration over moment(np) and pitch angle(nth)
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit ! fact_balk=5
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
            Drr_ = ( Drrfow(nth,np,nr+1,nsa)/JIl_p & 
                 + Drrfow(nth,np,nr,nsa)/JIl_m )*0.5d0
           
!******** original heatfow                        
!            heatfow_out(nr,nsa) = heatfow_out(nr,nsa) &
!                                + (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &  
!                                * Drr_*dfdrhom(nth,np,nr,nsa)*dVI &
!                                / (AEE*1.D3)  !****unit convert [J] to [keV]
!******* for dfdrhom is made of fnsp*dVI 
           heatfow_out(nr,nsa) = heatfow_out(nr,nsa) &
                                + (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &  
                                * Drr_* 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
                                / (AEE*1.D3)  !****unit convert [J] to [keV]

!**** fnsp sumVI_l
!            sumVI_l(nr,nsa) = sumVI_l(nr,nsa) +2.d0/3.d0* fnsp(nth,np,nr,nsa)*dVI &
!                 *(pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)*AEE*1.D3) ! AF220220

!**** dfdrhom is made of fnsp*dVI
            sumVI(nr,nsa) = sumVI(nr,nsa) + 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
                 *(pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)*AEE*1.D3) ! anzai[2022/3/4]

!**** original sumVI
!            sumVI(nr,nsa) = sumVI(nr,nsa) + dfdrhom(nth,np,nr,nsa)*dVI &
!                 *(pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)*AEE*1.D3) ! AF220220
         sumdf(nr,nsa)=sumdf(nr,nsa) + dfdrhom(nth,np,nr,nsa)
          end do
        end do
      end do
    end do

    !****first order derivative
    do nsa = 1, nsamax
!      call first_order_derivative(sumVI(:,nsa), sumVI_l(:,nsa), rm)
!      call first_order_derivative(heatfow_out(:,nsa), heatfow_l(:,nsa), rm)
      do nr = 1, nrmax
        heatfow_out(nr,nsa) = heatfow_out(nr,nsa)/sumVI(nr,nsa)
      end do
    end do

  end subroutine integral_Heatdiff
  
  subroutine Heat_rw_baverage(heatrwav, heatrw)
  !---------------------------------------------
  ! heat diffusion coef calcul module
  ! This module corresponds to integral_Heatdiff
  !--------------------------------------------- 
    
    use fpcomm
    use fowcomm
    
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: heatrw, heatrwav
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom, fnsp_l
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
    double precision,dimension(nrmax,nsamax) :: nu, Ta, Drw
    double precision step_len, step_time, eps_t, rho_theta, Baxis, pv, rho
    double precision dVI, sumVI
    integer nth, np, nr, nsa, ns

    !**** initialization ****
    Baxis = Bing(1)
    heatrw(:,:) = 0.d0
    heatrwav(:,:) = 0.d0
    sumVI = 0.d0
 !********** for new dfdrhom module [2022/3/4]
     fnsp_l(:,:,:,:)=0.d0
     do nsa = 1, nsamax
       ns = ns_nsa(nsa)
       do nr= 1, nrmax
         do np = 1, npmax
           if ( pm(np,ns) > fact_bulk ) exit
           do nth = 1, nthmax
           dVI = delp(ns)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)
           fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa) * dVI
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
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/( 1.5d0*rnsl(nr,nsa)*1.d20) / AEE/1.D3
        !**** Ta(temperature)[keV] RKEV = AEE*1.D3

        nu(nr,nsa) = (1.D0/6.4D14)*rnsl(nr,nsa)*1.D20/(Ta(nr,nsa)**1.5D0)
        nu(nr,nsa) = nu(nr, nsa)/(rm(nr)*RA/RR)*sqrt(AMFP(1)/AMFP(nsa)) 
        !**** Effective collisional flecuency nu
      end do
    end do
   
    !**** Calculation of particle diff coef
    do nsa = 1, nsamax
      do nr = 1, nrmax
        eps_t = rm(nr)*RA/RR
        do np = 1, npmax
          step_time = 1.d0/nu(nr,nsa)

          rho = pm(np,nsa)*ptfp0(nsa)/(aefp(nsa)*Baxis)
          !****Unit [keV]

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
            dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)

            heatrw(nr,nsa) = heatrw(nr,nsa)&
                           + (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                           * Drwlocal(nth,np,nr,nsa) & 
                           * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
                           / (AEE*1.D3)  !**** Unit convert [J] to [keV]
            heatrwav(nr,nsa) = heatrwav(nr,nsa) &
                             + (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                             * Drwba(nth,np,nr,nsa) & 
                             * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
                             / (AEE*1.D3)  !****Unit convert [J] to [keV]

            sumVI = sumVI + dfdrhom(nth,np,nr,nsa)&
                  * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)*AEE*1.D3)
          end do
        end do
        heatrw(nr, nsa)  = heatrw(nr, nsa)/sumVI
        heatrwav(nr,nsa) = heatrwav(nr,nsa)/sumVI
      end do
    end do

  end subroutine Heat_rw_baverage
  
  subroutine Heatneoclass(heatba,heatpla)
    !---------------------------------------------------------------
    !subroutine for neoclassical heat flux in banana, plateau region
    !calculate heat flux and diffusion coefficient
    ! main calcul using[J]
    !--------------------------------------------------------------

    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: heatba, heatpla
    double precision,dimension(nrmax,nsamax) :: Sr_ba, Sr_pla
    double precision,dimension(nrmax,nsamax) :: Ta, dTadr, dndr
    !double precision,dimension(nrmax) :: B_r
    double precision, dimension(nsamax) :: rho_a
    double precision fact, fact_s, tau_ele,tau_i, eps_t, Baxis, B_p
    integer nr, nsa

    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bing(1) ! approximation on B by Baxis

    !****temperature make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
        !Ta(temperature)[keV]
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

        tau_ele = fact*sqrt(AMFP(1))*((Ta(nr,1)*RKEV)**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**2*AEE**2*lnlam(nr,2,1))
        tau_i   = fact*sqrt(2.d0)*sqrt(AMFP(2))*((Ta(nr,2)*RKEV)**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2))
        !**** *RKEV converts [keV] to [J] 
      
        rho_a(nsa) = sqrt(Ta(nr,nsa)*RKEV/AMFP(nsa))*AMFP(nsa)/(AEFP(nsa)*Baxis)
        !**** *RKEV converts [keV] to [J] 
        eps_t = rm(nr)*RA/RR
        B_p = rm(nr)*RA*BB/(safety_factor(nr)*RR)
        !****calculate heat flux
        if (nsa == 1) then
         fact_s = (safety_factor(nr)**2) &
             *(rho_a(nsa)**2)/((eps_t**1.5)*tau_ele)
         Sr_ba(nr,nsa) = fact_s*rnsl(nr,nsa)*1.D20*Ta(nr,1)*RKEV & 
            * ((1.53d0*(1+Ta(nr,2)/Ta(nr,1))*dndr(nr, nsa))/rnsl(nr,nsa) &
              - 1.81d0*dTadr(nr,1)/Ta(nr,1) &
              - 0.27d0*dTadr(nr,2)/Ta(nr,1))
        ! neglect E_para term
        
         Sr_pla(nr,nsa) = - RKEV*sqrt(PI)/4*eps_t**2*Ta(nr,nsa) & 
               / (AEE*B_p)*rho_a(nsa)*rnsl(nr,nsa)*1.d20*Ta(nr,nsa) & 
               / rm(nr)*( (1+Ta(nr,2)/Ta(nr,1))*dndr(nr,nsa) &
               / (rnsl(nr,nsa)*1.d20) &
               + 7.5d0*dTadr(nr,nsa)/Ta(nr,nsa) & 
               + 1.5d0*dTadr(nr,2)/Ta(nr,nsa))
        ! neglect E_para term
        else
         fact_s = (safety_factor(nr)**2) &
              *(rho_a(nsa)**2)/((eps_t**1.5)*tau_i)
         Sr_ba(nr,nsa) = - 0.68d0*fact_s*(1 + 0.48d0*sqrt(eps_t)) &
              *rnsl(nr,nsa)*1.d20*dTadr(nr,2)*RKEV
         Sr_pla(nr,nsa) = - RKEV*1.5d0*sqrt(pi) &
                        * eps_t**2*Ta(nr,nsa)/(AEE*B_p) &
                        * rho_a(nsa)/rm(nr)*rnsl(nr,nsa) & 
                        * 1.d20*dTadr(nr,nsa) 
        end if  

        heatba(nr,nsa) = - Sr_ba(nr, nsa) &
             /(rnsl(nr,nsa)*1.d20*dTadr(nr, nsa)*RKEV)
        heatpla(nr,nsa) = - Sr_pla(nr, nsa) &
             /(rnsl(nr,nsa)*1.d20*dTadr(nr, nsa)*RKEV)
      end do
    end do

  end subroutine Heatneoclass

!===============================================
! Subroutine for factor check
!===============================================
  subroutine check_factorneo(eps_t,cyclo_rho, dTadr,Ta,tau_ele,tau_i)
  !-------------------------------------
  ! calculate main factors in neoclassical theory
  !--------------------------------------
  
    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: cyclo_rho, dTadr
    double precision,dimension(nrmax,nsamax),intent(out) :: Ta
    double precision,dimension(nrmax),intent(out) :: eps_t, tau_ele, tau_i
    double precision Baxis, fact
    integer nsa, nr

    Baxis = Bing(1)
    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    !****temperature make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
        !Ta(temperature)[keV]
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
        tau_ele(nr) = fact*sqrt(AMFP(1))*((Ta(nr,1)*RKEV)**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**2*AEE**2*lnlam(nr,2,1))
        tau_i(nr)   = fact*sqrt(2.d0)*sqrt(AMFP(2))*((Ta(nr,2)*RKEV)**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2))
        !**** *RKEV converts [keV] to [J]
    end do

  end subroutine check_factorneo 

end module check_neoclass
