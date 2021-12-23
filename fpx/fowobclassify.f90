! fowobclassify.f90

MODULE fowobclassify

  USE fpcomm,ONLY: rkind
  PRIVATE

  PUBLIC :: fow_ob_classify

  REAL(rkind),ALLOCATABLE,dimension(:) :: theta_m,xi
  REAL(rkind),ALLOCATABLE,dimension(:,:,:) :: beta_D,beta_stag
  REAL(rkind),ALLOCATABLE,dimension(:,:,:) :: beta_pinch,theta_pinch
  REAL(rkind),ALLOCATABLE,dimension(:,:) :: xi_Xtype_boundary
  PUBLIC:: theta_m,xi,beta_D,beta_stag,beta_pinch,theta_pinch,xi_Xtype_boundary

CONTAINS

  SUBROUTINE fow_ob_classify
    USE fpcomm
    USE fowcomm
    USE fpwrite
    IMPLICIT NONE

    call prep_orbit_classify
    call pinch_orbit
    call D_orbit
    call stagnation_orbit
    call stagnation_type

  end subroutine fow_ob_classify

  SUBROUTINE prep_orbit_classify
    USE fpcomm
    USE fowcomm
    IMPLICIT NONE
    INTEGER,SAVE:: nthmax_ob_save=0
    INTEGER,SAVE:: nrmax_save=0
    INTEGER,SAVE:: nsamax_save=0
    INTEGER:: nth

    IF(nthmax_ob_save.NE.nthmax_ob.OR. &
       nrmax.NE.nrmax_save.OR. &
       nsamax.NE.nsamax_save) THEN
       IF(ALLOCATED(theta_m)) THEN
          DEALLOCATE(theta_m,xi,beta_D,beta_stag,beta_pinch,theta_pinch)
          DEALLOCATE(xi_Xtype_boundary)
       END IF

       ALLOCATE(theta_m(nthmax_ob))
       ALLOCATE(xi(nthmax_ob))
       ALLOCATE(beta_D(nthmax_ob,nrmax,nsamax))
       ALLOCATE(beta_stag(nthmax_ob,nrmax,nsamax))
       ALLOCATE(beta_pinch(nthmax_ob,nrmax,nsamax))
       ALLOCATE(theta_pinch(nthmax_ob,nrmax,nsamax))
       ALLOCATE(xi_Xtype_boundary(nrmax,2))

       nthmax_ob_save=nthmax_ob
       nrmax_save=nrmax
       nsamax_save=nsamax
    END IF

    theta_m(1) = 0.d0
    xi(1) = 1.d0
    DO nth = 2, nthmax_ob
       theta_m(nth) = theta_m(nth-1)+pi/dble(nthmax_ob-1)
       xi(nth) = COS( theta_m(nth) )
    END do

  END SUBROUTINE prep_orbit_classify

  subroutine pinch_orbit
    use fpcomm
    use fowcomm
    implicit none
    REAL(rkind) :: psipin
    integer :: nr,nsa,ierr,nth

    ierr = 0

    do nsa = 1, nsamax
      do nr = 1, nrmax
          do nth = 1, nthmax_ob
            psipin = ( psim(nr)-1.d-8 )/(nthmax_ob+1)*dble(nth+1)
            call get_pinch_point(beta_pinch(nth,nr,nsa), &
                 theta_pinch(nth,nr,nsa), psipin, psim(nr), nsa)
          end do
      end do
    end do

  end subroutine pinch_orbit

  subroutine get_pinch_point(beta_pich,theta_pncp,psip_in, psim_in, nsa_in)
    ! input psip_in that is pinch point of pinch orbit with psi_m = psim_in.
    ! psi_p must be less than psim(nr_in)
    ! return momentum of pinch orbit with psi_m = psim(nr_in), psi_pnc=psip_in
    use fowcomm
    use fpcomm
    USE fowlib
    implicit none
    REAL(rkind),INTENT(OUT):: beta_pich,theta_pncp
    real(rkind),intent(in):: psip_in, psim_in
    integer,intent(in) :: nsa_in

    real(rkind) :: F_pncp, Bin_pncp, BFFB, ps_ratio
    REAL(rkind):: dFdpsi_pncp, dBdpsi_pncp, G_m, C(3)
    REAL(rkind):: w, FB_prime, xi_pncp, xi2
    real(rkind) :: F_m, B_m, p_ret
    COMPLEX(rkind) :: z(2)
    real(rkind),allocatable :: dFdpsi(:), dBdpsi(:)

    allocate(dFdpsi(nrmax+1), dBdpsi(nrmax+1))

    call first_order_derivative(dFdpsi, Fpsig, psimg)
    call first_order_derivative(dBdpsi, Bing, psimg)

    call fow_cal_spl(F_pncp, psip_in, Fpsig, psimg)
    call fow_cal_spl(dFdpsi_pncp, psip_in, dFdpsi, psimg)
    call fow_cal_spl(F_m, psim_in, Fpsig, psimg)
    call fow_cal_spl(B_m, psim_in, Boutg, psimg)    
    
    if ( psip_in == 0.d0 ) then
      Bin_pncp = Bing(1)
      dBdpsi_pncp = dBdpsi(1)
    else
      call fow_cal_spl(Bin_pncp, psip_in, Bing, psimg)
      call fow_cal_spl(dBdpsi_pncp, psip_in, dBdpsi, psimg)  
    end if

    ! dFdpsi_pncp = dFdpsi_pncp*psip_in
    ! dBdpsi_pncp = dBdpsi_pncp*psip_in
    dFdpsi_pncp = dFdpsi_pncp*psim_in
    dBdpsi_pncp = dBdpsi_pncp*psim_in

    G_m = aefp(nsa_in)*B_m*psim_in/(amfp(nsa_in)*vc*F_m)

    ps_ratio = 1.d0-psip_in/psim_in

    BFFB = 2.d0*Bin_pncp*dFdpsi_pncp - F_pncp*dBdpsi_pncp

    C(3) = -4.d0*B_m*Bin_pncp**3*F_m**2&
          +B_m**2*(ps_ratio*BFFB+2.d0*Bin_pncp*F_pncp)**2

    C(2) = -4.d0*(Bin_pncp-B_m)*Bin_pncp**3*F_m**2&
           -2.d0*B_m**2*dBdpsi_pncp*F_pncp*ps_ratio&
           *(ps_ratio*BFFB+2.d0*Bin_pncp*F_pncp)

    C(1) = (B_m*dBdpsi_pncp*F_pncp*ps_ratio)**2

    call solve_quadratic_equation(z, C)

    w = real(z(1))

    xi2 = 1.d0-(1.d0-w)*B_m/Bin_pncp
    ! write(*,*)"pnc",xi2,w
    if ( xi2 >= 1.d0 ) xi2 = 1.d0
    if ( xi2 <= 0.d0 ) then
      if( ABS(psip_in-psim_in) >= ABS(psip_in-0.d0) )p_ret = pmax(nsa_in)
      if( ABS(psip_in-psim_in) < ABS(psip_in-0.d0) ) p_ret = 0.d0
      return
    end if

    if ( aefp(nsa_in) >= 0.d0 ) then
      xi_pncp = sqrt(xi2)
    else 
      xi_pncp = -sqrt(xi2)
    end if

    FB_prime = (dFdpsi_pncp*Bin_pncp-F_pncp*dBdpsi_pncp)/Bin_pncp**2*B_m/F_m
    p_ret = G_m*sqrt(w) &
         /(w*FB_prime-0.5d0*(1.d0-xi_pncp**2)*F_pncp*dBdpsi_pncp/F_m/Bin_pncp)
                ! LHS = gamma*beta
    p_ret = vc*sqrt(p_ret**2/(1.d0+p_ret**2))
                ! LHS = velocity of stagnation orbit
    beta_pinch = MIN( p_ret, vc ) / vc
    theta_pncp = ACOS(xi_pncp)

  end subroutine get_pinch_point
  
  subroutine D_orbit
    ! Used for visualization
    ! beta_D is maximum momentum of trapped particles
    ! for given psi_m, xi and particle species
    use fowcomm
    use fpcomm
    implicit none
    integer :: nth, nr, nsa, ir
    real(rkind) :: dnr_bn, v_D_orbit
    real(rkind),dimension(nthmax_ob,nrmax,nsamax) :: Gm
    real(rkind),dimension(nthmax_ob,nrmax) :: psin, B_m, Bn

    do nsa = 1, nsamax

      do nr = 1, nrmax
        do nth = 1, nthmax_ob
          if ( aefp(nsa)*xi(nth) < 0.d0 ) then
            B_m(nth,nr) = Bin(nr)
          else if ( aefp(nsa)*xi(nth) > 0.d0 ) then
            B_m(nth,nr) = Bout(nr)
          end if
        end do
      end do
  
      do nr = 1, nrmax
        do nth = 1, nthmax_ob
           Gm(nth,nr,nsa) = aefp(nsa)*B_m(nth,nr)*psim(nr) &
                /(amfp(nsa)*vc*Fpsi(nr))
        end do
      end do
  
      do nr = 1, nrmax
        do nth = 1, nthmax_ob
          if ( xi(nth) == 0.d0 ) continue

          if ( ABS(xi(nth)) == 1.d0 ) then
            Bn(nth,nr) = B_m(nth,nr)/(1.d0-xi(nth)**2+1.d-8)
          else
            Bn(nth,nr) = B_m(nth,nr)/(1.d0-xi(nth)**2)
          end if
    
          if( Bn(nth,nr)<=Boutg(1) ) then
            do ir = 1, nrmax
              if( Boutg(ir+1)<=Bn(nth,nr) ) then
                dnr_bn = (Bn(nth,nr)-Boutg(ir))/(Boutg(ir+1)-Boutg(ir))
                psin(nth,nr) = psimg(ir)+(psimg(ir+1)-psimg(ir))*dnr_bn
                exit
              else if ( Bn(nth,nr)<Boutg(nrmax+1) ) then
                 dnr_bn = (Bn(nth,nr)-Boutg(nrmax+1)) &
                      /(Boutg(nrmax+1)-Boutg(nrmax))
                 psin(nth,nr) = psimg(nrmax+1) &
                      +(psimg(nrmax+1)-psimg(nrmax))*dnr_bn
                exit
              end if
            end do
          else
            do ir = 1, nrmax
              if( Bn(nth,nr)<=Bing(ir+1) ) then
                dnr_bn = (Bn(nth,nr)-Bing(ir))/(Bing(ir+1)-Bing(ir))
                psin(nth,nr) = psimg(ir)+(psimg(ir+1)-psimg(ir))*dnr_bn
                exit
              else if ( Bn(nth,nr)>Bing(nrmax+1) ) then
                dnr_bn = (Bn(nth,nr)-Bing(nrmax+1))/(Bing(nrmax+1)-Bing(nrmax))
                psin(nth,nr) = psimg(nrmax+1) &
                     +(psimg(nrmax+1)-psimg(nrmax))*dnr_bn
                exit
              end if
            end do
          end if
      
        end do
      end do
  
      do nr = 1, nrmax
        do nth = 1, nthmax_ob
          if ( xi(nth) == 0.d0 ) then
            beta_D(nth,nr,nsa) = 0.d0
            cycle
          end if
          if ( aefp(nsa)*xi(nth) < 0.d0 ) then
            beta_D(nth,nr,nsa) = 0.d0
          end if
          
          if( Boutg(nrmax+1)<=Bn(nth,nr) .and. &
               Bn(nth,nr)<=Bing(nrmax+1) .and. psin(nth,nr)<=psim(nr) ) then
             v_D_orbit = vc/sqrt(1.d0 &
                  +xi(nth)**2/(Gm(nth,nr,nsa)*(1.d0-psin(nth,nr)/psim(nr)))**2)
             ! beta_D(nth,nr,nsa) &
             !     = amfp(nsa)*v_D_orbit/sqrt(1.d0-v_D_orbit**2/vc**2)
             ! beta_D(nth,nr,nsa) &
             !     = beta_D(nth,nr,nsa)/ptfp0(nsa) ! normalize
             beta_D(nth,nr,nsa) = v_D_orbit/vc
          else
             beta_D(nth,nr,nsa) = 0.d0
          end if
        end do
      end do

    end do
  
  end subroutine D_orbit

  subroutine stagnation_orbit
    ! Used for visualization
    ! beta_stag is maximum momentum of not-forbitten particles &
    !   for given psi_m, xi and particle species
    use fowcomm
    use fpcomm
    USE fowlib
    implicit none
    integer :: nth,nr,nsa
    real(rkind) :: v_stagnation_orbit, F_p, B_p
    real(rkind),dimension(nthmax_ob,nrmax,nsamax) :: Gm
    real(rkind),dimension(nthmax_ob,nrmax) :: B_m, dBmdpsi
    real(rkind),dimension(nrmax) :: dFdpsi

    do nsa = 1,nsamax

      do nr = 1,nrmax
        do nth = 1, nthmax_ob
          if ( aefp(nsa)*xi(nth) <= 0.d0 ) then
            B_m(nth,nr) = Bin(nr)
          else
            B_m(nth,nr) = Bout(nr)
          end if
        end do
      end do
  
        do nr = 1,nrmax
          do nth = 1, nthmax_ob
             Gm(nth,nr,nsa) = aefp(nsa)*B_m(nth,nr)*psim(nr) &
                  /(amfp(nsa)*vc*Fpsi(nr))
          end do
        end do

      call first_order_derivative(dFdpsi, Fpsi, psim)

      do nth = 1, nthmax_ob
        call first_order_derivative(dBmdpsi(nth,:), B_m(nth,:), psim)
      end do
    
      do nr = 1,nrmax
        do nth = 1, nthmax_ob
          if ( xi(nth) == 0.d0 ) then
            beta_stag(nth,nr,nsa) = 0.d0
            cycle
          end if

          F_p = dFdpsi(nr)*psim(nr)/Fpsi(nr)
          B_p = dBmdpsi(nth,nr)*psim(nr)/B_m(nth,nr)

          v_stagnation_orbit &
               = Gm(nth,nr,nsa)*xi(nth) &
               /(xi(nth)**2*F_p-0.5d0*(1.d0+xi(nth)**2)*B_p) ! LHS=gamma*beta 
          v_stagnation_orbit &
               = vc*sqrt(v_stagnation_orbit**2/(1.d0+v_stagnation_orbit**2))
                                         ! LHS = velocity of stagnation orbit
  
          ! beta_stag(nth,nr,nsa) &
          !  = amfp(nsa)*v_stagnation_orbit &
          !    /(1.d0-v_stagnation_orbit**2/vc**2)
          !                          momentum of stagnation orbit
          ! beta_stag(nth,nr,nsa) &
          !  = beta_stag(nth,nr,nsa)/ptfp0(nsa) &
          !                           normalize
          
          beta_stag(nth,nr,nsa) = v_stagnation_orbit/vc
        end do
      end do
    end do
  
  end subroutine stagnation_orbit

  subroutine stagnation_type
    use fpcomm
    use fowcomm
    use fpwrite
    USE fowlib
    implicit none
    ! F_p : ( dF / d(psi_p/psi_m) )/F 
    ! F_pp : ( d^2F / d(psi_p/psi_m)^2 )/F
    integer :: nr
    real(rkind) :: C(3), F_p, F_pp, B_p, B_pp
    COMPLEX(rkind) :: z(2)
    real(rkind),dimension(nrmax,2) :: dBmdpsi, d2Bmdpsi, B_, B__
    real(rkind),dimension(nrmax) :: dFdpsi, d2Fdpsi, F_, F__

    call first_order_derivative(dFdpsi, Fpsi, psim)
    call second_order_derivative(d2Fdpsi, Fpsi, psim)
    
    call first_order_derivative(dBmdpsi(:,1), Bout, psim)
    call first_order_derivative(dBmdpsi(:,2), Bin, psim)

    call second_order_derivative(d2Bmdpsi(:,1), Bout, psim)
    call second_order_derivative(d2Bmdpsi(:,2), Bin, psim)

    do nr = 1,nrmax
      F_(nr)=dFdpsi(nr)*psim(nr)/Fpsi(nr)
      F__(nr)=d2Fdpsi(nr)*psim(nr)**2/Fpsi(nr)
      B_(nr,1) = dBmdpsi(nr,1)*psim(nr)/Bout(nr)
      B__(nr,1) = d2Bmdpsi(nr,1)*psim(nr)**2/Bout(nr)
      B_(nr,2) = dBmdpsi(nr,2)*psim(nr)/Bin(nr)
      B__(nr,2) = d2Bmdpsi(nr,2)*psim(nr)**2/Bin(nr)
    end do

    do nr = 1,nrmax
      F_p = dFdpsi(nr)*psim(nr)/Fpsi(nr)
      F_pp = d2Fdpsi(nr)*psim(nr)**2/Fpsi(nr)

      ! calculate for xi > 0
      B_p = dBmdpsi(nr,1)*psim(nr)/Bout(nr)
      B_pp = d2Bmdpsi(nr,1)*psim(nr)**2/Bout(nr)

      C(3) = 4.d0*F_pp-4.d0*B_p*F_p-2.d0*B_pp+3.d0*B_p**2
      C(2) = 6.d0*B_p**2-2.d0*B_pp-4.d0*B_p*F_p
      C(1) = -1.d0*B_p**2
      ! write(*,*)"+nr",nr,F_p,F_pp,B_p,B_pp

      call solve_quadratic_equation(z, C)

      if ( aimag(z(2)) == 0.d0 ) then
        if ( real(z(2)) >= 1.d0 ) then
          xi_Xtype_boundary(nr,1) = 1.d0
        else if ( real(z(2)) <= 0.d0 ) then
          xi_Xtype_boundary(nr,1) = 0.d0
        else
          xi_Xtype_boundary(nr,1) = sqrt(real(z(2)))
        end if
      else
        xi_Xtype_boundary(nr,1) = 0.d0
      end if

      ! calculate for xi < 0
      B_p = dBmdpsi(nr,2)*psim(nr)/Bin(nr)
      B_pp = d2Bmdpsi(nr,2)*psim(nr)**2/Bin(nr)

      C(3) = 4.d0*F_pp-4.d0*B_p*F_p-2.d0*B_pp+3.d0*B_p**2
      C(2) = 6.d0*B_p**2-2.d0*B_pp-4.d0*B_p*F_p
      C(1) = -1.d0*B_p**2

      call solve_quadratic_equation(z, C)

      if ( aimag(z(2)) == 0.d0 ) then
        if ( real(z(2)) >= 1.d0 ) then
          xi_Xtype_boundary(nr,2) = 1.d0
        else if ( real(z(2)) <= 0.d0 ) then
          xi_Xtype_boundary(nr,2) = 0.d0
        else
          xi_Xtype_boundary(nr,2) = -sqrt(real(z(2)))
        end if
      else
        xi_Xtype_boundary(nr,2) = 0.d0
      end if

    end do

  end subroutine stagnation_type

  recursive function func_kaijou(n) result(m)
    implicit none
    integer,intent(in) :: n
    integer :: m

    if(n == 1) then
      m = 1
    else
      m = n*func_kaijou(n-1)
    end if

  end function func_kaijou

end module fowobclassify
