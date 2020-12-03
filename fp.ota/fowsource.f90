module fowsource
  private
  public :: fow_calculate_source

  type beam_quantities
    logical :: isLoss
    integer :: id, np_near
    double precision :: p, psip, psipw, delps, E, thetap, &
                        costh, sinth, total, B, F
    double precision :: psimw, thetamw, psim, thetam
  end type beam_quantities

  double precision,allocatable :: dBdpsip(:,:), dBdthp(:,:), dpsipdr(:)
  double precision,allocatable :: dBmdpsim(:,:), dFdpsip(:), psim0(:), psimg0(:)

  double precision,allocatable,dimension(:,:,:,:) :: S_beam_total

contains

  subroutine fow_calculate_source

    use fowcomm
    use fpcomm
    implicit none
    integer :: nth,np,nr,nsa

    allocate(S_beam_total(nthmax,npmax,nrmax,nsamax))

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            S_beam_total(nth,np,nr,nsa) = 0.d0
          end do
        end do
      end do
    end do

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

    use fpwrite

    implicit none
    real(rkind),allocatable :: S_beam(:,:,:,:,:)
    real(rkind) :: fact_thm, fact_psm, sum_S_beam, normalize, dV
    integer :: nth,np,nr,nsa, nthp, nb
    type(beam_quantities) :: beam

    allocate(S_beam(nthmax,npmax,nrmax,nsamax,nbeammax))
    allocate(dBdpsip(nrmax+1,nthpmax),dBdthp(nrmax+1,nthpmax),dpsipdr(nrmax))
    allocate(dFdpsip(nrmax),dBmdpsim(nrmax,2))
    allocate(psim0(nrmax),psimg0(nrmax+1))

    do nr = 1, nrmax
      psim0(nr) = psim(nr) * psi0
      psimg0(nr) = psimg(nr) * psi0
    end do
    psimg0(nrmax+1) = psimg(nrmax+1) * psi0


    call first_order_derivative(dFdpsip,Fpsi,psim0)
    call first_order_derivative(dBmdpsim(:,1),Bout,psim0)
    call first_order_derivative(dBmdpsim(:,2),Bin,psim0)
    call first_order_derivative(dpsipdr,psim,rm)

    do nr = 1, nrmax+1
      call first_order_derivative(dBdthp(nr,:),Babs(nr,:),theta_p)
    end do
    do nthp = 1, nthpmax
      call first_order_derivative(dBdpsip(:,nthp),Babs(:,nthp),psimg0)
    end do
        
    do nb = 1, nbeammax

      ! ion term ( if NSSPB(nb) = 1 then electron term )
      beam = construct_beam(NSSPB(nb), SPBTOT(nb), SPBR0(nb), &
                            SPBRW(nb), SPBENG(nb), SPBANG(nb), SPBPANG(nb))

      call calculate_beam_COMs(beam)

      if ( PMAX(beam%id) < beam%p ) then
        write(*,*)"Beam's momentum > PMAX  "
        cycle
      end if
  
      if ( beam%isLoss ) then
        write(*,*)"Beam's COM is in the loss region. "
        write(*,*)"psim",beam%psim
        write(*,*)"thetam",beam%thetam
        cycle
      end if

      call calculate_beam_COMs_width(beam)
      write(*,*)"psim",beam%psim,beam%psimw
      write(*,*)"thetam",beam%thetam,beam%thetamw
      write(*,*)"ip",beam%np_near

      sum_S_beam = 0.d0
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            if ( np == beam%np_near ) then

              fact_psm = -(psim(nr)-beam%psim)**2/beam%psimw**2
              fact_psm = EXP( fact_psm ) / SQRT( pi*beam%psimw**2 )

              if ( thetamg(nth,np,nr,beam%id) <= beam%thetam &
                  .and. beam%thetam <= thetamg(nth+1,np,nr,beam%id) ) then
                fact_thm = 1.d0
              else
                fact_thm = -(thetam(nth,np,nr,beam%id)-beam%thetam)**2/beam%thetamw**2
                fact_thm = EXP( fact_thm ) / SQRT( pi*beam%thetamw**2 )  
              end if
              S_beam(nth,np,nr,beam%id,nb) = fact_psm * fact_thm

              dV = delthm(nth,np,nr,beam%id)*delp(beam%id)*delps(nr)*Jacobian_I(nth,np,nr,beam%id)
              sum_S_beam = sum_S_beam + S_beam(nth,np,nr,beam%id,nb)*dV

            else
              S_beam(nth,np,nr,beam%id,nb) = 0.d0
            end if
          end do
        end do
      end do
      
      normalize = beam%total / sum_S_beam * ptfp0(beam%id)
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

        if ( PMAX(beam%id) < beam%p ) then
          write(*,*)"Beam's momentum > PMAX  "
          cycle
        end if
    
        if ( beam%isLoss ) then
          write(*,*)"Beam's COM is in the loss region. "
          write(*,*)"psim",beam%psim
          write(*,*)"thetam",beam%thetam
          cycle
        end if 

        call calculate_beam_COMs_width(beam)
        write(*,*)"psim",beam%psim,beam%psimw
        write(*,*)"thetam",beam%thetam,beam%thetamw
        write(*,*)"ip",beam%np_near

        sum_S_beam = 0.d0
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              if ( np == beam%np_near ) then
  
                fact_psm = -(psim(nr)-beam%psim)**2/beam%psimw**2
                fact_psm = EXP( fact_psm ) / SQRT( pi*beam%psimw**2 )
  
                if ( thetamg(nth,np,nr,beam%id) <= beam%thetam &
                    .and. beam%thetam <= thetamg(nth+1,np,nr,beam%id) ) then
                  fact_thm = 1.d0
                else
                  fact_thm = -(thetam(nth,np,nr,beam%id)-beam%thetam)**2/beam%thetamw**2
                  fact_thm = EXP( fact_thm ) / SQRT( pi*beam%thetamw**2 )  
                end if
                S_beam(nth,np,nr,beam%id,nb) = fact_psm * fact_thm

                dV = delthm(nth,np,nr,beam%id)*delp(beam%id)*delps(nr)*Jacobian_I(nth,np,nr,beam%id)
                sum_S_beam = sum_S_beam + S_beam(nth,np,nr,beam%id,nb)*dV
                                          
                write(*,*)"beam",sum_S_beam,S_beam(nth,np,nr,beam%id,nb)
              else
                S_beam(nth,np,nr,beam%id,nb) = 0.d0
              end if
            end do
          end do
        end do

        normalize = beam%total / (sum_S_beam * ptfp0(beam%id))
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
    deallocate(dBmdpsim)

  end subroutine fow_beam_source

  type(beam_quantities) function construct_beam(NS, TOT, R0, RW, ENG, ANG, PANG)
    use fowcomm
    use fpcomm

    implicit none

    integer,intent(in) :: NS
    real(rkind),intent(in) :: TOT, R0, RW, ENG, ANG, PANG
    real(rkind) :: dpsipdrl, PV

    construct_beam%id = NS

    PV = 1.d0 + ENG*aee/(amfp(NS)*vc**2)

    construct_beam%p = amfp(ns)*vc*SQRT( PV**2-1.d0 )/ptfp0(NS)

    ! construct_beam%p  = SQRT( 2.d0*amfp(NS)*ENG*aee )
    construct_beam%E  = ENG*aee
    construct_beam%thetap = PANG*pi/180.d0
    construct_beam%costh = cos(ANG*pi/180.d0)
    construct_beam%sinth = sin(ANG*pi/180.d0)

    call fow_cal_spl(construct_beam%psip, R0, psimg, rg)
    call fow_cal_spl(dpsipdrl, R0, dpsipdr, rm)
    construct_beam%psipw = dpsipdrl*RW

    call fow_cal_spl2D(construct_beam%B, construct_beam%psip, construct_beam%thetap, Babs, psimg, theta_p)
    call fow_cal_spl(construct_beam%F, construct_beam%psip, Fpsig, psimg)

    construct_beam%total = TOT

  end function

  subroutine calculate_beam_COMs(beam)
    use fpcomm
    use fowcomm
    use foworbit

    implicit none
    type(beam_quantities),intent(inout) :: beam
    real(rkind) :: thetaml, psiml, tau_loss
    integer :: ip, ierr


    call fow_set_obparm(ierr)
    call fow_cal_local_COMs(thetaml, psiml, tau_loss, beam%p*ptfp0(beam%id), ACOS( beam%costh ), beam%thetap, beam%psip, beam%id)
    beam%thetam = thetaml
    beam%psim   = psiml

    if ( tau_loss == 0 ) then
      beam%isLoss = .false.
    else
      beam%isLoss = .true.
    end if

    ! search nearest np to the momentum of beam particle
    do ip = 1, npmax
      if ( pg(ip,beam%id) <= beam%p .and. beam%p < pg(ip+1,beam%id) ) then
        beam%np_near = ip
        exit
      end if
    end do

  end subroutine calculate_beam_COMs

  subroutine calculate_beam_COMs_width(beam)
    ! calculate derivative of COM in term of poloidal flux and poloidal angle
    use fowcomm
    use fpcomm

    implicit none
    type(beam_quantities),intent(inout) :: beam
    integer :: nth ,np ,nr ,nsa, ns_beam, np_near
    real(rkind) :: dthmdpsip, dthmdthp, dpsimdpsip, dpsimdthp
    real(rkind) :: dBdpsipl, dBdthpl, dFdpsipl
    real(rkind) :: cthm, sthm, pl, Bm_, Fm_, dBm, dFm
    real(rkind) :: C(2,2), v_ps(2), v_th(2), det

    
    call fow_cal_spl(dFdpsipl,beam%psip*psi0,dFdpsip,psim0)
    call fow_cal_spl2D(dBdpsipl, beam%psip*psi0, beam%thetap, dBdpsip, psimg0, theta_p)
    call fow_cal_spl2D(dBdthpl, beam%psip*psi0, beam%thetap, dBdthp, psimg0, theta_p)

    ns_beam  = beam%id
    np_near  = beam%np_near

    pl = pm(np_near,ns_beam)*ptfp0(ns_beam)
    cthm = cos( beam%thetam )
    sthm = sin( beam%thetam )


    if ( cthm*aefp(ns_beam) >= 0.d0 ) then
      call fow_cal_spl(Bm_, beam%psim*psi0, Bout, psim0)
      call fow_cal_spl(dBm, beam%psim*psi0, dBmdpsim(:,1), psim0)
    else
      call fow_cal_spl(Bm_, beam%psim*psi0, Bin, psim0)
      call fow_cal_spl(dBm, beam%psim*psi0, dBmdpsim(:,2), psim0)
    end if
    call fow_cal_spl(Fm_, beam%psim*psi0, Fpsi, psim0)
    call fow_cal_spl(dFm, beam%psim*psi0, dFdpsip, psim0)

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
      v_ps(2) = beam%p*ptfp0(ns_beam)*beam%costh*( &
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
    beam%thetamw = ABS( dthmdpsip  * beam%psipw * psi0 ) ! + dthmdthp*beam%thetapw

  end subroutine calculate_beam_COMs_width

end module fowsource