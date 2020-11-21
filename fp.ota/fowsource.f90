module fowsource
  private
  public :: fow_calculate_source

  type beam_quantities
    integer :: id, np_near
    double precision :: p, psip, psipw, delps, E, thetap, &
                        costh, sinth, total, B, F
    double precision :: mu, Pzeta, psimw, thetamw, psim, thetam
  end type beam_quantities

  double precision,allocatable :: dBdpsip(:,:), dBdthp(:,:), dpsipdr(:)
  double precision,allocatable :: dBmdpsim(:,:), dFdpsip(:)
  double precision,allocatable :: thetap(:)

  double precision,allocatable,dimension(:,:,:,:) :: S_beam_total

contains

  subroutine fow_calculate_source

    use fowcomm
    use fpcomm
    use fowprep
    implicit none
    integer :: nth,np,nr,nsa

    allocate(S_beam_total(nthmax,npmax,nrmax,nsamax))

    call fow_beam_source

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            SPPB(nth,np,nr,nsa) = S_beam_total(nth,np,nr,nsa)
          end do
        end do
      end do
    end do

  end subroutine fow_calculate_source

  subroutine fow_beam_source

    use fowcomm
    use fpcomm
    use fowprep,only:fow_cal_spl, first_order_derivative

    use fpwrite

    implicit none
    real(rkind),allocatable :: S_beam(:,:,:,:,:)
    real(rkind) :: fact_thm, fact_psm, sum_S_beam, normalize, dV
    integer :: nth,np,nr,nsa, nthp, nb
    type(beam_quantities) :: beam

    allocate(S_beam(nthmax,npmax,nrmax,nsamax,nbeammax))
    allocate(dBdpsip(nrmax+1,nthpmax),dBdthp(nrmax+1,nthpmax),dpsipdr(nrmax))
    allocate(dFdpsip(nrmax),dBmdpsim(2,nrmax))
    allocate(thetap(nthpmax))

    do nthp = 1, nthpmax
      thetap(nthp) = dble(nthp-1)/dble(nthpmax)*2.d0*pi
    end do 

    call first_order_derivative(dFdpsip,Fpsi,psim)
    call first_order_derivative(dBmdpsim(1,:),Bout,psim)
    call first_order_derivative(dBmdpsim(2,:),Bin,psim)
    call first_order_derivative(dpsipdr,psim,rm)
    do nr = 1, nrmax+1
      call first_order_derivative(dBdthp(nr,:),Babs(nr,:),thetap)
    end do
    do nthp = 1, nthpmax
      call first_order_derivative(dBdpsip(:,nthp),Babs(:,nthp),psimg)
    end do
        
    do nb = 1, nbeammax

      ! ion term ( if NSSPB(nb) = 1 then electron term )
      beam = construct_beam(NSSPB(nb), SPBTOT(nb), SPBR0(nb), &
                            SPBRW(nb), SPBENG(nb), SPBANG(nb), SPBPANG(nb))

      call calculate_beam_COMs(beam)
      write(*,*)"beam_COM",beam%psim,beam%thetam

      call calculate_beam_COMs_width(beam)
      write(*,*)"beam_COM",beam%psim,beam%thetam

      sum_S_beam = 0.d0
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            if ( np == beam%np_near ) then
              fact_thm = -(thetam(nth,np,nr,beam%id)-beam%thetam)**2/beam%thetamw**2/2.d0
              fact_thm = EXP( fact_thm ) / SQRT( TWOPI*beam%thetamw**2 )

              fact_psm = -(psim(nr)-beam%psim)**2/beam%psimw**2/2.d0
              fact_psm = EXP( fact_psm ) / SQRT( TWOPI*beam%psimw**2 )

              S_beam(nth,np,nr,beam%id,nb) = fact_psm * fact_thm

              dV = delthm(nth,np,nr,beam%id)*delp(beam%id)*(psimg(nr+1)-psimg(nr))
              sum_S_beam = sum_S_beam + Jacobian_I(nth,np,nr,beam%id)*S_beam(nth,np,nr,beam%id,nb)*dV
            else
              S_beam(nth,np,nr,beam%id,nb) = 0.d0
            end if
          end do
        end do
      end do
      
      normalize = beam%total / sum_S_beam
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            S_beam(nth,np,nr,beam%id,nb) = S_beam(nth,np,nr,beam%id,nb) * normalize
          end do
        end do
      end do

      call fpcsv2D(S_beam(:,beam%np_near,:,beam%id,nb),"./csv/beam.csv")
      call fpcsv2D(thetam(:,beam%np_near,:,beam%id),"./csv/thetam_r.csv")


      ! electron term
      if ( NSSPB(nb) /= 1 ) then
        beam = construct_beam(1, SPBTOT(nb), SPBR0(nb), &
                SPBRW(nb), SPBENG(nb)*amfp(1)/amfp(NSSPB(nb)), SPBANG(nb), SPBPANG(nb))

        call calculate_beam_COMs(beam)

        call calculate_beam_COMs_width(beam)

        sum_S_beam = 0.d0
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              if ( np == beam%np_near ) then
                fact_thm = -(thetam(nth,np,nr,beam%id)-beam%thetam)**2/beam%thetamw**2/2.d0
                fact_thm = EXP( fact_thm ) / SQRT( TWOPI*beam%thetamw**2 )
  
                fact_psm = -(psim(nr)-beam%psim)**2/beam%psimw**2/2.d0
                fact_psm = EXP( fact_psm ) / SQRT( TWOPI*beam%psimw**2 )
  
                S_beam(nth,np,nr,beam%id,nb) = fact_psm * fact_thm
  
                dV = delthm(nth,np,nr,beam%id)*delp(beam%id)*(psimg(nr+1)-psimg(nr))
                sum_S_beam = sum_S_beam + Jacobian_I(nth,np,nr,beam%id)*S_beam(nth,np,nr,beam%id,nb)*dV
              else
                S_beam(nth,np,nr,beam%id,nb) = 0.d0
              end if
            end do
          end do
        end do
        
        normalize = beam%total / sum_S_beam
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              S_beam(nth,np,nr,beam%id,nb) = S_beam(nth,np,nr,beam%id,nb) * normalize
            end do
          end do
        end do
  
    
      end if

    end do

    deallocate(dBdpsip,dBdthp,dFdpsip)
    deallocate(dBmdpsim,thetap)

  end subroutine fow_beam_source

  type(beam_quantities) function construct_beam(NS, TOT, R0, RW, ENG, ANG, PANG)
    use fowcomm
    use fpcomm
    use fowprep,only:fow_cal_spl

    implicit none

    integer,intent(in) :: NS
    real(rkind),intent(in) :: TOT, R0, RW, ENG, ANG, PANG
    real(rkind) :: dpsipdrl, PV

    construct_beam%id = NS

    PV = 1.d0 + ENG*aee/(amfp(NS)*vc**2)

    construct_beam%p = amfp(ns)*vc*SQRT( PV**2-1.d0 )

    ! construct_beam%p  = SQRT( 2.d0*amfp(NS)*ENG*aee )
    construct_beam%E  = ENG*aee
    construct_beam%thetap = PANG*pi/180.d0
    construct_beam%costh = cos(ANG*pi/180.d0)
    construct_beam%sinth = sin(ANG*pi/180.d0)

    call fow_cal_spl(construct_beam%psip, R0, psimg, rg)
    call fow_cal_spl(dpsipdrl, R0, dpsipdr, rm)
    construct_beam%psipw = dpsipdrl*RW

    call fow_cal_spl2D(construct_beam%B, construct_beam%psip, construct_beam%thetap, Babs, psimg, thetap)
    call fow_cal_spl(construct_beam%F, construct_beam%psip, Fpsig, psimg)

    construct_beam%mu = construct_beam%E*construct_beam%sinth**2/construct_beam%B
    construct_beam%Pzeta = construct_beam%F/construct_beam%B*construct_beam%p*construct_beam%costh&
                          -aefp(NS)*construct_beam%psip

    construct_beam%total = TOT
  end function

  subroutine calculate_beam_COMs(beam)

    use plinit
    use obcomm
    use obinit
    use obparm
    use obprep
    use obcalc

    use fowcomm,only:psi0fp => psi0
    use fpcomm,only:pmfp => pm, ipmax =>npmax, ptfp0, ithmax => nthmax
    use fowprep,only:fow_cal_spl

    implicit none
    type(beam_quantities),intent(inout) :: beam
    integer :: nstp, nmax, i, im, ip, ierr = 0
    real(rkind),allocatable :: theta_ob_tmp(:), psip_ob_tmp(:), vpara_ob_tmp(:), vperp_ob_tmp(:)
    real(rkind) :: psimax, n_psim, cthm, vparam, vperpm &
                  ,diff_thetap_min, diff_thetap, thetap_eq(6)

    call ob_init
    call ob_parm(1,'../fp.ota/fpparm',ierr)
    call ob_prep(ierr)
    ns_ob = beam%id
    nobt_max = 1
    call ob_allocate
    pcangle_in(1) = beam%costh
    penergy_in(1) = beam%E/aee/1.0d3
    psipn_in(1)   = beam%psip/psi0fp
    theta_in(1)   = beam%thetap*180.d0/pi
    zeta_in(1)    = 0.d0

    call ob_calc(ierr)

    psimax = -1.d0
    do nstp = 0, nstp_max_nobt(1)
      if ( psip_ob(nstp,1) >= psimax ) then
        psimax = psip_ob(nstp,1)
        nmax = nstp
      end if
    end do

    if ( ABS( theta_ob(nmax,1) ) <= pi/dble(2*ithmax) ) then
      beam%psim = psip_ob(nmax,1)
      cthm = vpara_ob(nmax,1)/SQRT( vpara_ob(nmax,1)**2 + vperp_ob(nmax,1)**2 )
      beam%thetam = ACOS( cthm )

    else
      allocate(theta_ob_tmp(0:nmax), psip_ob_tmp(0:nmax), vpara_ob_tmp(0:nmax), vperp_ob_tmp(0:nmax))
      do nstp = 0, nmax
        theta_ob_tmp(nstp) = theta_ob(nstp,1)
        psip_ob_tmp(nstp)  = psip_ob(nstp,1)
        vpara_ob_tmp(nstp) = vpara_ob(nstp,1)
        vpara_ob_tmp(nstp) = vpara_ob(nstp,1)
      end do

      ! search poloidal angle where psip = psim.  (psim occurs on the equater.)
      thetap_eq = [-2.d0*pi, -1.d0*pi, 0.d0, pi, 2.d0*pi, 3.d0*pi]
      diff_thetap_min = 1.d8
      do i = 1, 6
        diff_thetap = ABS( thetap_eq(i)-theta_ob(nmax,1) )
        if ( diff_thetap < diff_thetap_min ) then
          diff_thetap_min = diff_thetap
          im = i
        end if
      end do

      call fow_cal_spl(beam%psim, thetap_eq(im), psip_ob_tmp, theta_ob_tmp)
      call fow_cal_spl(vparam, thetap_eq(im), vpara_ob_tmp, theta_ob_tmp)
      call fow_cal_spl(vperpm, thetap_eq(im), vperp_ob_tmp, theta_ob_tmp)

      cthm = vparam / SQRT( vparam**2+vperpm**2 )
      beam%thetam = ACOS( cthm )

    end if

    ! search nearest np to the momentum of beam particle
    do ip = 1, ipmax
      if ( pmfp(ip,beam%id)*ptfp0(beam%id) >= beam%p ) then
        beam%np_near = ip
        exit
      end if
    end do

  end subroutine calculate_beam_COMs

  subroutine calculate_beam_COMs_width(beam)
    ! calculate derivative of COM in term of poloidal flux and poloidal angle
    use fowcomm
    use fpcomm
    use fowprep,only:fow_cal_spl

    implicit none
    type(beam_quantities),intent(inout) :: beam
    integer :: nth ,np ,nr ,nsa, ns_beam, np_near
    real(rkind) :: dthmdpsip, dthmdthp, dpsimdpsip, dpsimdthp
    real(rkind) :: dBdpsipl, dBdthpl, dFdpsipl
    real(rkind) :: cthm, sthm, pl, Bm_, Fm_, dBm, dFm
    real(rkind) :: C(2,2), v_ps(2), v_th(2), det


    call fow_cal_spl(dFdpsipl,beam%psip,dFdpsip,psim)
    call fow_cal_spl2D(dBdpsipl, beam%psip, beam%thetap, dBdpsip, psimg, thetap)
    call fow_cal_spl2D(dBdthpl, beam%psip, beam%thetap, dBdthp, psimg, thetap)

    ns_beam  = beam%id
    np_near  = beam%np_near

    pl = pm(np_near,ns_beam)*ptfp0(ns_beam)
    cthm = cos( beam%thetam )
    sthm = sin( beam%thetam )


    if ( cthm*aefp(ns_beam) >= 0.d0 ) then
      call fow_cal_spl(Bm_, beam%psim, Bout, psim)
      call fow_cal_spl(dBm, beam%psim, dBmdpsim(1,:), psim)
    else
      call fow_cal_spl(Bm_, beam%psim, Bin, psim)
      call fow_cal_spl(dBm, beam%psim, dBmdpsim(2,:), psim)
    end if
    call fow_cal_spl(Fm_, beam%psim, Fpsi, psim)
    call fow_cal_spl(dFm, beam%psim, dFdpsip, psim)

    C(1,1) = 2.d0*sthm*cthm
    C(1,2) = -1.d0*sthm**2*dBm/Bm_**2
    C(2,1) = -1.d0*pl*sthm*Fm_/Bm_
    C(2,2) = pl*cthm*(dFm*Bm_-Fm_*dBm)/Bm_**2 - aefp(ns_beam)

    det = C(1,1)*C(2,2)-C(1,2)*C(2,1)

    if ( det == 0.d0 ) then
      dthmdpsip  = 0.d0
      dpsimdpsip = 0.d0
      ! dthmdthp   = 0.d0
      ! dpsimdthp) = 0.d0

    else
      v_ps(1) = -1.d0*beam%sinth**2/beam%B**2*dBdpsipl
      v_ps(2) = beam%p*beam%costh*( &
        dFdpsipl*beam%B-beam%F*dBdpsipl &
      )/beam%B**2-aefp(ns_beam)

      dthmdpsip  = (C(2,2)*v_ps(1)-C(1,2)*v_ps(2))/det
      dpsimdpsip = (-1.d0*C(2,1)*v_ps(1)+C(1,1)*v_ps(2))/det

      ! v_th(1) = -1.d0*sin0**2*dBdthpl/B0**2
      ! v_th(2) = -1.d0*F0*p0*cos0*dBdthpl/B0**2

      ! dthmdthp  = (C(2,2)*v_th(1)-C(1,2)*v_th(2))/det
      ! dpsimdthp = (-1.d0*C(2,1)*v_th(1)+C(1,1)*v_th(2))/det  
      
    end if

    beam%psimw   = ABS( dpsimdpsip * beam%psipw ) ! + dpsimdthp*beam%thetapw
    beam%thetamw = ABS( dthmdpsip  * beam%psipw ) ! + dthmdthp*beam%thetapw

  end subroutine calculate_beam_COMs_width

  subroutine fow_cal_spl2D(f_out, x_in, y_in, f, x, y)
    ! 2D-version of fow_cal_spl
    use fpcomm,only:rkind

    implicit none
    real(rkind),intent(out) :: f_out
    real(rkind),intent(in) :: x_in, y_in, f(:,:), x(:), y(:)
    integer :: i, j, imax, jmax, ierr = 0
    real(rkind),allocatable :: U(:,:,:,:), fx(:,:), fy(:,:) , fxy(:,:)

    imax = size(x)
    jmax = size(y)

    allocate(U(4,4,imax,jmax),fx(imax,jmax),fy(imax,jmax),fxy(imax,jmax))

    if ( size(f,1) /= imax &
        .or. size(f,2) /= jmax ) then
      write(*,*)"error at fow_cal_spl2D"
      STOP
    end if

    call SPL2D(x,y,f,fx,fy,fxy,U,imax,imax,jmax,0,0,IERR)

    call SPL2DF(x_in,y_in,f_out,x,y,U,imax,imax,jmax,ierr)

  end subroutine

end module fowsource