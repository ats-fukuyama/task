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
                                                 
    double precision,dimension(nrmax,nsamax) :: Drw, Drwav, nud_mono, D_mono, Dbanana
    double precision,dimension(nrmax,nsamax) :: Dnewba, Dnewpla, jaceffect
    integer nth, np, nr, nsa, nsb, mode(3)

    call system('mkdir -p dat')

    call integral_Drr(Drhmrhm)
    call D_random_walk_baverage(Drwav, Drw)
    call nu_deflection_monoenergy(nud_mono)
!    call D_banana(Dbanana)
    call newneoclass_ba(Dnewba) ![2022/1/31] editted by anzai
    call newneoclass_pla(Dnewpla) ![2022/1/31] editted by anzai
    call check_jeffect(jaceffect) ![2022/2/5] editted by anzai

    call fptxt2D(Drhmrhm,"dat/Drhmrhm.txt")

    call fptxt2D(Drw,"dat/Drw.txt")
    call fptxt2D(Drwav,"dat/Drwav.txt")
!    call fptxt2D(Dbanana,"dat/Dbanana.txt")
    call fptxt2D(Dnewba, "dat/Dnewba.txt") ![2022/1/31] edited by anzai
    call fptxt2D(Dnewpla, "dat/Dnewpla.txt") ![2022/1/31] editted by anzai
    call fptxt2D(jaceffect, "dat/jaceffect.txt")![2022/2/5] editted by anzai
    call fptxt2D(nud_mono,"dat/nud_mono.txt")

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

!  subroutine D_banana(Dbanana)
!    use fpcomm
!    use fowcomm
!    implicit none
!    double precision,dimension(nrmax,nsamax),intent(out)  :: Dbanana
!    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom
!    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
!    double precision,dimension(nrmax,nsamax) :: nud
!    double precision step_len, step_time, eps_t, rho, Baxis, pv, Ta, vtha
!    double precision dVI, sumVI
!    integer nth, np, nr, nsa
!
!    Baxis = Bing(1)
!
!    call cal_nu_ei(nud)
!
!    do nsa = 1, nsamax
!      do np = 1, npmax
!        do nth = 1, nthmax
!          call first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
!        end do
!      end do
!    end do
!
!    do nsa = 1, nsamax
!      do nr = 1, nrmax
!        eps_t = rm(nr)*RA/RR
!        Ta = rwsl(nr,nsa)*1.d6/( 1.5d0*rnsl(nr,nsa)*1.d20 )
!        vtha = SQRT( 2.d0*Ta/amfp(nsa) )
!        rho = amfp(nsa)*vtha/(ABS(aefp(nsa))*Baxis)
!        Dbanana(nr,nsa) = safety_factor(nr)**2*rho**2*nud(nr,nsa)/eps_t**1.5d0
!        WRITE(6,'(A,I4,6ES12.4)') 'nr:',nr,eps_t,Ta,vtha,rho,nud(nr,nsa)
!      end do
!    end do
!
!  end subroutine

  subroutine cal_nu_ei(nu)
  !**************************************
  ![2022/2/8] Modified by anzai
  !**************************************
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

  subroutine D_random_walk_baverage(Drwav, Drw)
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: Drw, Drwav
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
!    double precision,dimension(npmax,nrmax,nsamax) :: nud
    double precision,dimension(nrmax,nsamax) :: nu, Ta
    double precision step_len, step_time, eps_t, rho_theta, Baxis, pv, B_theta, rho
    double precision dVI, sumVI
    integer nth, np, nr, nsa

    Baxis = Bing(1)

    !call nu_deflection(nud)
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
        !Ta(temperature)[keV]
        nu(nr,nsa) = (1.D0/6.4D14)*rnsl(nr,nsa)*1.D20/(Ta(nr,nsa)**1.5D0)
        nu(nr,nsa) = nu(nr, nsa)/(rm(nr)*RA/RR)*(AMFP(1)/AMFP(nsa)) !Effective nu
      end do
    end do
    
    do nsa = 1, nsamax 
      do nr = 1, nrmax
        eps_t = rm(nr)*RA/RR
        ! B_theta = safety_factor(nr)/RR*dpsimdr(nr)/RA
        do np = 1, npmax
          !step_time = 1.d0/nud(np,nr,nsa)
          step_time = 1.d0/nu(nr,nsa)
          rho = pm(np,nsa)*ptfp0(nsa)/(aefp(nsa)*Baxis)
          do nth = 1, nthmax
            !step_len = rho
            step_len = rho * safety_factor(nr)/sqrt(eps_t)
            !Drwlocal(nth,np,nr,nsa) = step_len**2/step_time*1.46d0*sqrt(eps_t)
            Drwlocal(nth,np,nr,nsa) = step_len**2/step_time*sqrt(eps_t)
          end do
        end do
      end do
    end do

    call bounce_average_for_Drw(Drwba, Drwlocal)

    do nsa = 1, nsamax
      do nr = 1, nrmax
        Drw(nr,nsa) = 0.d0
        Drwav(nr,nsa) = 0.d0
        sumVI = 0.d0
        do np = 1, npmax
          do nth = 1, nthmax
            dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)

            Drw(nr,nsa) = Drw(nr,nsa)&
                          + Drwlocal(nth,np,nr,nsa)*dfdrhom(nth,np,nr,nsa)*dVI
            Drwav(nr,nsa) = Drwav(nr,nsa)&
                          + Drwba(nth,np,nr,nsa)*dfdrhom(nth,np,nr,nsa)*dVI
            sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
          end do
        end do
        Drw(nr,nsa) = Drw(nr,nsa)/sumVI
        Drwav(nr,nsa) = Drwav(nr,nsa)/sumVI
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

            numerator  = rnsl(nr,nsb)*1.d20 &
                 *aefp(nsa)**2*aefp(nsb)**2 * lnlam(nr,nsb,nsa)
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

  subroutine nu_deflection_monoenergy(nud_mono)
    use fpcomm
    use fowcomm
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out) :: nud_mono
    double precision x, fx, gx
    double precision nudb, fact, vthb, Tb, vtha, Ta
    double precision nuBra, factBra, numerator, denominator
    integer nr,nsa,nsb,nssa,nssb

    fact = 3.d0/2.d0*SQRT( pi/2.d0 )
    factBra = 1.d0/( 12.d0*pi**1.5d0*eps0**2 )

    do nsa = 1, nsamax
      nssa = ns_nsa(nsa)
      do nr = 1, nrmax
        nud_mono(nr,nsa) = 0.d0
        Ta = rwsl(nr,nsa)*1.d6/( 1.5d0*rnsl(nr,nsa)*1.d20 )
        vtha = SQRT( 2.d0*Ta/amfp(nsa) )
        do nsb = 1, nsbmax
          nssb = ns_nsb(nsa)
          Tb = rwsl(nr,nsb)*1.d6/( 1.5d0*rnsl(nr,nsb)*1.d20 )
          vthb = SQRT( 2.d0*Tb/amfp(nsb) )

          numerator  = rnsl(nr,nsb)*1.d20 * ( aefp(nsa)*aefp(nsb) )**2 * lnlam(nr,nsb,nsa)
          denominator= SQRT( amfp(nsb) ) * Tb**1.5d0 * ( amfp(nsa)/amfp(nsb) )**2
          nuBra = factBra*numerator/denominator

          x = vtha/vthb
          call gosakn(x,fx,gx)

          nudb = fact*nuBra*(fx-gx)/x**3

          nud_mono(nr,nsa) = nud_mono(nr,nsa)+nudb

        end do

      end do
    end do

  end subroutine nu_deflection_monoenergy

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
!**********************[2022/1/31] added by anzai***************
  subroutine newneoclass_ba(Dnewba)
    !----------------------------------
    !subroutine for neoclassical particle flux in banana region
    !calculate particle flux and diffusion coefficient
    !----------------------------------
    
    use fpcomm
    use fowcomm
    
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: Dnewba
    double precision,dimension(nrmax,nsamax) :: Sr_nba, Ta, dTadr, dndr
    !double precision,dimension(nrmax) :: B_r
    double precision fact, fact_s, tau_ele, rho_e, eps_t, Baxis, B_p
    integer nr, nsa

    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bing(1) ! approximation on B by Baxis
 
    !****temperature make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)
      end do
    end do

    !****first order derivation
    do nsa = 1, nsamax
      call first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
      call first_order_derivative(dndr(:,nsa), rnsl(:, nsa)*1.d20, rm)
    end do
    
    !****calculate flux and diffusion coefficient
    do nsa = 1, nsamax
      do nr = 1, nrmax  
        tau_ele = fact*sqrt(AMFP(1))*((Ta(nr,1))**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**2*lnlam(nr,2,1)*AEE**2)
        rho_e = sqrt(2*Ta(nr,1)/AMFP(1))*AMFP(1)/(AEE*Baxis)
        eps_t = rm(nr)*RA/RR
        !B_p = dpsimdr(nr)*RA/(safety_factor(nr)*RR) ! careful
        fact_s = (safety_factor(nr)**2) &
             *(rho_e**2)/((eps_t**1.5)*tau_ele)
        Sr_nba(nr,nsa) = fact_s*(-1.22d0*(1+Ta(nr,2)/Ta(nr,1))*dndr(nr, nsa))! &
            ! + 4.3d-1*rnsl(nr, nsa)*1.D20*dTadr(nr,1)/Ta(nr,1) &
            ! + 1.9d-1*rnsl(nr, nsa)*1.D20*dTadr(nr,2)/Ta(nr,1) ) 

        ! neglect E_para term
        
        Dnewba(nr,nsa) = - Sr_nba(nr, nsa)/dndr(nr, nsa)
!        Dnewba(nr,nsa) = fact_s*1.22d0*(1+Ta(nr,2)/Ta(nr,1))

       ! WRITE(6,'(A,I4,6ES12.4)') 'Dn:',nr, &
            ! TA(nr,1)/AEE,safety_factor(nr),tau_ele,rho_e, &
            ! rho_e**2/tau_ele,Dnewba(nr,nsa)
      end do 
    end do

  end subroutine newneoclass_ba	

  subroutine newneoclass_pla(Dnewpla)
    !---------------------------------
    !subroutine for neoclassical particle flux in plateau region
    !calculate particle flux and diffusion coefficient
    !---------------------------------
    
    use fpcomm
    use fowcomm
    
    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: Dnewpla
    double precision,dimension(nrmax,nsamax) :: Sr_npla, Ta, dTadr, dndr
    !double precision,dimension(nrmax) :: B_r
    double precision fact, fact_s, tau_ele, rho_e, eps_t, Baxis, B_p
    integer nr, nsa

    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bing(1) ! approximation on B by Baxis

    !****temperature make
    do nsa = 1, nsamax
      do nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)
      end do
    end do

    !****first order derivation
    do nsa = 1, nsamax
      call first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
      call first_order_derivative(dndr(:,nsa), rnsl(:, nsa)*1.d20, rm)
    end do

    !****calculate flux and diffusion coefficient
    do nsa = 1, nsamax
      do nr = 1, nrmax
        tau_ele = fact*sqrt(AMFP(1))*Ta(nr,1)**1.5d0/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**2*lnlam(nr,2,1)*AEE**2)
        rho_e = sqrt(2*Ta(nr,1)/AMFP(1))*AMFP(1)/(AEE*Baxis)
        eps_t = rm(nr)*RA/RR
        B_p = dpsimdr(nr)*RA/(safety_factor(nr)*RR) ! careful
        fact_s = -sqrt(pi)*(eps_t**2)*Ta(nr,1) &
             *rho_e*(rnsl(nr,nsa)*1.d20)/(2*AEE*B_p*rm(nr))        
        
        Sr_npla(nr,nsa) = fact_s*((1+Ta(nr,2)/Ta(nr,1))*dndr(nr,nsa) &
             /(rnsl(nr,nsa)*1.d20))! &
            ! + 1.5d0*dTadr(nr,1)/Ta(nr,1) + 1.5d0*dTadr(nr,2)/Ta(nr,1) )
        !neglect B term
   
        Dnewpla(nr,nsa) = - Sr_npla(nr, nsa)/dndr(nr, nsa)
      end do 
    end do

  end subroutine newneoclass_pla

  subroutine check_jeffect(jaceffect)
  !-------------------------------------
  !Integrrate jacobian over theta_m and p
  !--------------------------------------
  
    use fpcomm
    use fowcomm

    implicit none
    double precision,dimension(nrmax,nsamax),intent(out)  :: jaceffect
     double precision,dimension(nthmax,npmax,nrmax,nsamax) :: dfdrhom
    double precision dVI, sumVI, JIl_p, JIl_m
    integer nth,np,nr,nsa,ns

     !****make dfdrhom
     do nsa = 1, nsamax
      do np = 1, npmax
        do nth = 1, nthmax
          call first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
        end do
      end do
    end do
   
    !****integrate
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        jaceffect(nr,nsa) = 0.d0
        sumVI = 0.d0
        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit
          do nth = 1, nthmax
          !-------------------------------------
          !p means plus, m means minus
          !-------------------------------------
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

            jaceffect(nr,nsa) = jaceffect(nr,nsa) + dfdrhom(nth,np,nr,nsa)*dVI
            sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
            
           ! WRITE(*,*) "jacobian:", jaceffect(nr,nsa), sumVI, dfdrhom(nth,np,nr,nsa)
          end do
        end do
        jaceffect(nr,nsa) = jaceffect(nr,nsa)/sumVI
      end do
    end do

  end subroutine check_jeffect


  ! subroutine banana_width_and_omega_bounce(w_b, omega_b)
  !   use fpcomm
  !   use fowcomm
  !   implicit none
  !   double precision,dimension(nthmax,npmax,nrmax,nsamax),intent(out) :: w_b, omega_b
  !   double precision,allocatable :: drdt(:), U(:,:)
  !   type(orbit) ob
  !   integer nth, np, nr, nsa, nstp, nstpmax, ierr
  !   double precision r_min

  !   ierr = 0

  !   do nsa = 1, nsamax
  !     do nr = 1, nrmax
  !       do np = 1, npmax
  !         do nth = 1, nthmax
  !           ob = orbit_m(nth,np,nr,nsa)
  !           nstpmax = ob%nstp_max
  !           omega_b(nth,np,nr,nsa) = 2.d0*pi/ob%time(nstpmax)
  !           allocate(drdt(nstpmax))
  !           allocate(U(4, nstpmax))
  !           call spl1D(ob%time, ob%r, drdt, U, nstpmax, 0, ierr)
  !           call spl1DF(ob%time(nstpmax)*0.5d0, r_min, ob%time, U ,nstpmax, ierr)
  !           w_b(nth,np,nr,nsa) = ( rm(nr)-r_min )*RA
  !           deallocate(drdt)
  !           deallocate(U)
  !         end do
  !       end do
  !     end do
  !   end do

  ! end subroutine banana_width_and_omega_bounce

  ! subroutine banana_width_and_omega_bounce_model(w_b, omega_b)
  !   use fpcomm
  !   use fowcomm
  !   implicit none
  !   double precision,dimension(nthmax,npmax,nrmax,nsamax),intent(out) :: w_b, omega_b
  !   double precision :: rho, eps,sign_prev,sign,Baxis,vl,pv,B_theta,thetapl
  !   type(orbit) :: ob
  !   logical :: isTrapped
  !   integer :: nth,np,nr,nsa,ns,nstp,nstpmax

  !   Baxis = Boutg(1)

  !   do nsa = 1, nsamax
  !     ns = ns_nsa(nsa)
  !     do nr = 1, nrmax
  !       eps = rm(nr)*RA/RR
  !       do np = 1, npmax
  !         pv = SQRT(1.d0+theta0(nsa)*pm(np,nsa)**2)
  !         vl = vc*SQRT(1.d0-1.d0/pv**2)
  !         do nth = 1, nthmax
  !           ob = orbit_m(nth,np,nr,nsa)
  !           nstpmax = ob%nstp_max

  !           isTrapped = .false.
  !           if ( ob%costh(1) < 0.d0 ) then
  !             sign_prev = -1.d0
  !           else if ( ob%costh(1) > 0.d0 ) then
  !             sign_prev = 1.d0
  !           else
  !             sign_prev = 0.d0
  !           end if

  !           do nstp = 2, nstpmax
  !             if ( ob%costh(nstp) /= 0.d0 ) then
  !               sign = ob%costh(nstp) / ABS( ob%costh(nstp) )
  !             else
  !               sign = 0.d0
  !             end if

  !             if ( sign_prev*sign > 0.d0 ) then
  !               cycle
  !             else if ( sign_prev*sign < 0.d0 ) then
  !               isTrapped = .True.
  !             else
  !               sign_prev = sign
  !             end if
  !           end do

  !           if ( isTrapped ) then
  !             rho = pm(np,ns)*ptfp0(ns)/ABS( pz(ns)*aee*Baxis )
  !             w_b(nth,np,nr,nsa) = safety_factor(nr)*rho/SQRT(eps)
  !             omega_b(nth,np,nr,nsa) = SQRT(eps)*vl/( safety_factor(nr)*RR )
  !           else
  !             if ( pz(ns)*COS( thetam(nth,np,nr,nsa) ) > 0.d0 ) then
  !               thetapl = 0.d0
  !             else
  !               thetapl = pi
  !             end if
  !             B_theta = safety_factor(nr)/( RR+RA*rm(nr)*COS( thetapl ) )*dpsimdr(nr)
  !             w_b(nth,np,nr,nsa) = 0.d0
  !             omega_b(nth,np,nr,nsa) = B_theta*vl*COS( thetam(nth,np,nr,nsa) )&
  !                                     /( Baxis*rm(nr)*RA )
  !             omega_b(nth,np,nr,nsa) = ABS( omega_b(nth,np,nr,nsa) )
  !             ! omega_b(nth,np,nr,nsa) = 0.d0
  !           end if

  !         end do
  !       end do
  !     end do
  !   end do

  ! end subroutine banana_width_and_omega_bounce_model

  ! subroutine mean_transmatrix(dIdu,mode)
  !   use fpcomm
  !   use fowcomm
  !   use fowcoef
  !   implicit none
  !   double precision,dimension(:,:,:,:,:,:),intent(out) :: dIdu
  !   integer,dimension(3),intent(in) :: mode
  !   double precision,dimension(3,3,max_stp) :: dIdul
  !   double precision :: dt, taup
  !   integer :: nth,np,nr,nsa,nstp,nstpmax,i,j
  !   type(orbit) :: ob

  !   do nsa = 1, nsamax
  !     do nr = 1, nrmax+mode(3)
  !       do np = 1, npmax+mode(2)
  !         do nth = 1, nthmax+mode(1)
  !           do j = 1, 3
  !             do i = 1, 3
  !               dIdu(i,j,nth,np,nr,nsa) = 0.d0
  !             end do
  !           end do
  !         end do
  !       end do
  !     end do
  !   end do

  !   do nsa = 1, nsamax
  !     do nr = 1+mode(3), nrmax+mode(3)
  !       do np = 1+mode(2), npmax+mode(2)
  !         do nth = 1, nthmax+mode(1)

  !           if ( mode(1) == 0 .and. mode(2) == 0 .and. mode(3) == 0 ) then
  !             ob = orbit_m(nth,np,nr,nsa)
  !           else if ( mode(1) == 1 .and. mode(2) == 0 .and. mode(3) == 0 ) then
  !             ob = orbit_th(nth,np,nr,nsa)
  !             if ( nth == nth_stg(nsa) ) cycle
  !           else if ( mode(1) == 0 .and. mode(2) == 1 .and. mode(3) == 0 ) then
  !             ob = orbit_p(nth,np,nr,nsa)
  !           else if ( mode(1) == 0 .and. mode(2) == 0 .and. mode(3) == 1 ) then
  !             ob = orbit_r(nth,np,nr,nsa)
  !           end if

  !           nstpmax = ob%nstp_max
  !           taup = ob%time(nstpmax)
  !           call transformation_matrix(dIdul, ob, nth, np, nr, nsa, mode)

  !           do nstp = 2, nstpmax
  !             dt = ( ob%time(nstp)-ob%time(nstp-1) )/taup
  !             do j = 1, 3
  !               do i = 1, 3
  !                 dIdu(i,j,nth,np,nr,nsa) = dIdu(i,j,nth,np,nr,nsa) + dIdul(i,j,nstp)*dt
  !               end do
  !             end do

  !           end do

  !         end do
  !       end do
  !     end do
  !   end do

  ! end subroutine  mean_transmatrix

  ! subroutine D_random_walk_baverage_old(Drw)
  !   use fpcomm
  !   use fowcomm
  !   implicit none
  !   double precision,dimension(nrmax,nsamax),intent(out)  :: Drw
  !   double precision,dimension(nthmax,npmax,nrmax,nsamax) :: w_b, omega_b, dfdrhom
  !   double precision,dimension(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
  !   double precision,dimension(npmax,nrmax,nsamax) :: nud
  !   double precision step_len, step_time
  !   double precision dVI, sumVI
  !   integer nth, np, nr, nsa

  !   call banana_width_and_omega_bounce(w_b, omega_b)
  !   call nu_deflection(nud)

  !   do nsa = 1, nsamax
  !     do np = 1, npmax
  !       do nth = 1, nthmax
  !         call first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
  !       end do
  !     end do
  !   end do

  !   do nsa = 1, nsamax
  !     do nr = 1, nrmax
  !       do np = 1, npmax
  !         step_time = 1.d0/nud(np,nr,nsa)
  !         do nth = 1, nthmax
  !           step_len = w_b(nth,np,nr,nsa)
  !           Drwlocal(nth,np,nr,nsa) = step_len**2/step_time
  !         end do
  !       end do
  !     end do
  !   end do

  !   call bounce_average_for_Drw(Drwba, Drwlocal)

  !   do nsa = 1, nsamax
  !     do nr = 1, nrmax
  !       Drw(nr,nsa) = 0.d0
  !       sumVI = 0.d0
  !       do np = 1, npmax
  !         do nth = 1, nthmax
  !           dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)

  !           Drw(nr,nsa) = Drw(nr,nsa)&
  !                         + Drwba(nth,np,nr,nsa)*dfdrhom(nth,np,nr,nsa)*dVI
  !           sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
  !         end do
  !       end do
  !       Drw(nr,nsa) = Drw(nr,nsa)/sumVI
  !     end do
  !   end do

  ! end subroutine D_random_walk_baverage_old

  ! subroutine ionHeatFlux(q_ion)
  !   use fpcomm
  !   use fowcomm
  !   implicit none
  !   double precision,dimension(nrmax),intent(out) :: q_ion
  !   double precision,dimension(nthmax,npmax,nrmax) :: dfdrhom
  !   double precision :: Drr_, dVI, K, pv
  !   integer :: nth,np,nr,nsa,ns

  !   do np = 1, npmax
  !     do nth = 1, nthmax
  !       call first_order_derivative(dfdrhom(nth,np,:), fnsp(nth,np,:,2), rm)
  !       do nr = 1, nrmax
  !         dfdrhom(nth,np,nr) = dfdrhom(nth,np,nr)/RA
  !       end do
  !     end do
  !   end do

  !   do nr = 1, nrmax
  !     q_ion(nr) = 0.d0
  !     do np = 1, npmax
  !       if ( pm(np,2) > fact_bulk ) exit
  !       do nth = 1, nthmax
  !         dVI = delp(2)*delthm(nth,np,nr,2)*JIR(nth,np,nr,2)
  !         Drr_ = ( Drrfow(nth,np,nr+1,2)+Drrfow(nth,np,nr,2) )*0.5d0/JI(nth,np,nr,2)
  !         pv = SQRT(1.d0+theta0(2)*pm(np,2)**2)
  !         K = amfp(2)*vc**2*(pv-1.d0)

  !         q_ion(nr) = q_ion(nr) - K*Drr_*dfdrhom(nth,np,nr)*dVI*1.d20
  !       end do
  !     end do
  !     write(*,*)"ion",q_ion(nr)
  !   end do

  ! end subroutine ionHeatFlux

  ! subroutine D_neoclass()
  !   !delete D_banana, D_plateau, nu_ei, nu_p, nu_b as variables [2022/1/17]
  !   use fowcomm
  !   use fpcomm

  !   implicit none
  !   double precision,dimension(nrmax),intent(out) :: D_banana, D_plateau, nu_ei, nu_p, nu_b
  !   integer :: nr, nsb, nssb
  !   double precision :: dens ,Te, Vt, omega_e, Baxis, rho_e, epst
  !   double precision :: factSpitz

  !   Baxis = Bing(1)
  !   omega_e = aee*Baxis/ame
  !   factSpitz = aee**4/( 93.d0*eps0**2*SQRT(ame) )
  !   do nr = 1, nrmax
  !     Te = rwsl(nr,1)*1.d6/( 1.5d0*rnsl(nr,1)*1.d20 ) ! [J]
  !     nu_ei(nr) = 0.d0
  !     do nsb = 1, nsbmax
  !       nssb = ns_nsb(nsb)
  !       if ( nssb == 1 ) cycle
  !       dens = rnsl(nr,nsb)*1.d20
  !       nu_ei(nr) = nu_ei(nr) + pz(nssb)**2*dens*lnlam(nr,nsb,1)
  !     end do
  !     nu_ei(nr) = nu_ei(nr)*factSpitz/Te**1.5d0
  !     ! nu_ei = rnud(nr,2,1)

  !     Vt = SQRT( Te/ame )
  !     rho_e = Vt/omega_e
  !     epst = rm(nr)*RA/RR
  !     nu_p(nr) = Vt/( safety_factor(nr)*RR )
  !     nu_b(nr) = epst**1.5d0*nu_p(nr)

  !     !D_banana (nr) = ( safety_factor(nr)*rho_e )**2 * nu_ei(nr) / epst**1.5d0
  !     D_plateau(nr) = ( safety_factor(nr)*rho_e )**2 * nu_p(nr)
  !   end do

  ! end subroutine D_neoclass

  ! subroutine D_random_walk(Drw)
  !   use fpcomm
  !   use fowcomm
  !   implicit none
  !   double precision,dimension(nrmax,nsamax),intent(out)  :: Drw
  !   double precision,dimension(nthmax,npmax,nrmax,nsamax) :: w_b, omega_b, dfdrhom
  !   double precision,dimension(npmax,nrmax,nsamax)        :: nud
  !   double precision step_len, step_time, Drwlocal
  !   double precision dVI, sumVI
  !   integer nth, np, nr, nsa

  !   call banana_width_and_omega_bounce(w_b, omega_b)
  !   call nu_deflection(nud)

  !   do nsa = 1, nsamax
  !     do np = 1, npmax
  !       do nth = 1, nthmax
  !         call first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
  !       end do
  !     end do
  !   end do

  !   do nsa = 1, nsamax
  !     do nr = 1, nrmax
  !       Drw(nr,nsa) = 0.d0
  !       sumVI = 0.d0

  !       do np = 1, npmax
  !         step_time = 1.d0/nud(np,nr,nsa)
  !         do nth = 1, nthmax
  !           step_len = w_b(nth,np,nr,nsa)/RA
  !           Drwlocal = step_len**2/step_time
            
  !           dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)
  !           sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
  !           Drw(nr,nsa) = Drw(nr,nsa) + Drwlocal*dfdrhom(nth,np,nr,nsa)*dVI

  !         end do
  !       end do

  !       Drw(nr,nsa) = Drw(nr,nsa)/sumVI

  !     end do
  !   end do

  ! end subroutine D_random_walk

end module check_neoclass
