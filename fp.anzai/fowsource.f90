module fowsource
  private
  public :: fow_calculate_source

  type beam_quantities
    logical :: isLoss
    integer :: ns, np_near
    double precision :: p, rhomw, thetamw, rhom, thetam, total
  end type beam_quantities

contains

  subroutine fow_calculate_source(nt)
    use fowcomm
    use fpcomm
    use foworbit

    implicit none
    integer,intent(in) :: nt
    integer :: nth,np,nr,nsa,nb,ierr
    type(beam_quantities) :: beam
    real(rkind),allocatable :: sppb_nb(:,:,:,:,:)

    ierr = 0

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            SPPB(nth,np,nr,nsa) = 0.d0
          end do
        end do
      end do
    end do

    ! Beam
    allocate(sppb_nb(nthmax,npmax,nrmax,nsamax,nbeammax))
    if( nt == 0 ) call fow_set_obparm(ierr)

    do nb = 1, nbeammax

      ! ion term
      beam = construct_beam(NSSPB(nb), SPBTOT(nb), SPBR0(nb), &
                            SPBRW(nb), SPBENG(nb), SPBANG(nb), SPBPANG(nb))
      call beam_source(sppb_nb(:,:,:,NSSPB(nb),nb), beam)

      ! electron term
      if ( NSSPB(nb) /= 1.d0 ) then
        beam = construct_beam(1, SPBTOT(nb), SPBR0(nb), &
                              SPBRW(nb), SPBENG(nb)*amfp(1)/amfp(NSSPB(nb)), SPBANG(nb), SPBPANG(nb))
        call beam_source(sppb_nb(:,:,:,1,nb), beam)
      end if

      do nsa = 1, nsamax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              sppb(nth,np,nr,nsa) = sppb_nb(nth,np,nr,nsa,nb)
            end do
          end do
        end do
      end do

    end do

  end subroutine fow_calculate_source

  subroutine beam_source(s_beam, beam)
    use fowcomm
    use fpcomm

    implicit none
    real(rkind),dimension(nthmax,npmax,nrmax),intent(out) :: S_beam
    type(beam_quantities),intent(in) :: beam
    real(rkind) :: ex_thm, ex_rm, fact_thm, fact_rm
    real(rkind) :: n_total, normalize, dV
    integer :: nth, np, nr, nsa


    if ( PMAX(beam%ns) < beam%p ) then
      write(*,*)"Beam's momentum > PMAX. "
      return
    end if

    if ( beam%isLoss ) then
      write(*,*)"Beam's COM is in the loss region. "
      write(*,*)"psim",beam%rhom
      write(*,*)"thetam",beam%thetam
      return
    end if

    write(*,*)"rhom",beam%rhom,beam%rhomw
    write(*,*)"thetam",beam%thetam,beam%thetamw
    write(*,*)"ip",beam%np_near

    n_total = 0.d0
    do nr = 1, nrmax
      do np = 1, npmax
        do nth = 1, nthmax
          if ( np == beam%np_near ) then

            ex_rm = -(rm(nr)-beam%rhom)**2/beam%rhomw**2
            fact_rm = EXP( ex_rm ) / SQRT( pi*beam%rhomw**2 )

            if ( thetamg(nth,np,nr,beam%ns) <= beam%thetam &
                .and. beam%thetam <= thetamg(nth+1,np,nr,beam%ns) ) then
              fact_thm = 1.d0
            else
              ex_thm = -(thetam(nth,np,nr,beam%ns)-beam%thetam)**2/beam%thetamw**2
              fact_thm = EXP( ex_thm ) / SQRT( pi*beam%thetamw**2 )
            end if

            S_beam(nth,np,nr) = fact_rm * fact_thm

            dV = delthm(nth,np,nr,beam%ns)*delp(beam%ns)*JI(nth,np,nr,beam%ns)

            n_total = n_total + S_beam(nth,np,nr)*dV*volr(nr)

          else
            S_beam(nth,np,nr) = 0.d0
          end if
          
        end do
      end do
    end do

    if ( n_total > 0.d0 ) then
      normalize = beam%total/n_total
    else
      normalize = 0.d0
    end if

    do nr = 1, nrmax
      do np = 1, npmax
        do nth = 1, nthmax
          S_beam(nth,np,nr) = S_beam(nth,np,nr)*normalize
        end do
      end do
    end do

  end subroutine beam_source

  type(beam_quantities) function construct_beam(NS, TOT, R0, RW, ENG, ANG, PANG)! result(beam)
    use fowcomm
    use fpcomm
    use foworbit

    implicit none
    integer,intent(in) :: NS
    real(rkind),intent(in) :: TOT, R0, RW, ENG, ANG, PANG
    real(rkind) :: energy_b, thetap_b, costh_b, sinth_b, PV, theta_b, tau_loss
    real(rkind) :: dBdr_b, dBdthp_b, dFdr_b, dpsipdr_b, B_b, F_b, psip_b
    real(rkind) :: thetam_b, rhom_b, cthm, sthm, dBmdrl, dFmdrl, Bml, Fml, dpsimdrl, pl
    real(rkind) :: dthmdr, dthmdthp, drmdr, drmdthp
    real(rkind) :: A(2,2), b(2), detA, THPW
    integer :: np, ierr

    ierr = 0

    construct_beam%ns = NS

    PV = 1.d0 + ENG*aee/(amfp(NS)*vc**2)

    construct_beam%p = amfp(ns)*vc*SQRT( PV**2-1.d0 )/ptfp0(NS)
    pl = construct_beam%p*ptfp0(NS)

    energy_b = ENG*aee
    thetap_b = PANG*pi/180.d0
    theta_b  = ANG*pi/180.d0
    costh_b  = COS( theta_b )
    sinth_b  = SIN( theta_b )

    construct_beam%total = TOT

    call fow_cal_spl(psip_b, R0, psim, rm)
    call fow_cal_spl(dpsipdr_b, R0, dpsimdr, rm)
    call fow_cal_spl(F_b, R0, Fpsi, rm)
    call fow_cal_spl(dFdr_b, R0, dFdr, rm)
    call fow_cal_spl2D(B_b, R0, thetap_b, Babs, rm, theta_p)
    call fow_cal_spl2D(dBdr_b, R0, thetap_b, dBdr, rm, theta_p)
    call fow_cal_spl2D(dBdthp_b, R0, thetap_b, dBdthp, rm, theta_p)

    ! calculate COM of beam
    call fow_cal_local_COMs(thetam_b, rhom_b, tau_loss, pl, theta_b, thetap_b, psip_b, NS)

    construct_beam%thetam = thetam_b
    construct_beam%rhom   = rhom_b
    cthm = COS( thetam_b )
    sthm = SIN( thetam_b )

    if ( tau_loss <= 0.d0 ) then
      construct_beam%isLoss = .false.
    else
      construct_beam%isLoss = .true.
    end if

    ! search nearest np to the momentum of beam particle
    do np = 1, npmax
      if ( pg(np,ns) <= construct_beam%p .and. construct_beam%p < pg(np+1,ns) ) then
        construct_beam%np_near = np
        exit
      end if
    end do

    ! calculate beam width in COM space
    call fow_cal_spl(Fml, rhom_b, Fpsi, rm)
    call fow_cal_spl(dFmdrl, rhom_b, dFdr, rm)
    call fow_cal_spl(dpsimdrl, rhom_b, dpsimdr, rm)
    if ( cthm*aefp(ns) >= 0.d0 ) then
      call fow_cal_spl(Bml, rhom_b, Bout, rm)
      call fow_cal_spl2D(dBmdrl, rhom_b, 0.d0, dBdr, rm, theta_p)
    else
      call fow_cal_spl(Bml, rhom_b, Bin, rm)
      call fow_cal_spl2D(dBmdrl, rhom_b, pi, dBdr, rm, theta_p)
    end if

    A(1,1) = 2.d0*sthm*cthm/Bml
    A(1,2) = -1.d0*sthm**2*dBmdrl/Bml**2
    A(2,1) = Fml/Bml*pl*sthm
    A(2,2) = aefp(ns)*dpsimdrl-(dFmdrl*Bml-Fml*dBmdrl)/Bml**2*pl*cthm
    detA   = A(1,1)*A(2,2)-A(1,2)*A(2,1)

    if ( detA /= 0.d0 ) then
      b(1) = -1.d0*sinth_b**2/B_b**2*dBdr_b
      b(2) = aefp(ns)*dpsipdr_b&
            - ( dFdr_b*B_b-F_b*dBdr_b )/B_b**2*pl*costh_b
      dthmdr = (A(2,2)*b(1)-A(1,2)*b(2))/detA
      drmdr  = (A(1,1)*b(2)-A(2,1)*b(1))/detA

      ! b(1) = -1.d0*sinth_b**2/B_b**2*dBdr_b
      ! b(2) = F_b/B_b**2*pl*costh_b*dBdthp_b
      ! dthmdthp = (A(2,2)*b(1)-A(1,2)*b(2))/detA
      ! drmdthp  = (A(1,1)*b(2)-A(2,1)*b(1))/detA

    else
      dthmdr   = 0.d0 
      drmdr    = 1.d0
      dthmdthp = 0.d0
      drmdthp  = 0.d0

    end if

    construct_beam%rhomw   = drmdr *RW !+ drmdthp *THPW
    construct_beam%thetamw = dthmdr*RW !+ dthmdthp*THPW

  end function

end module fowsource