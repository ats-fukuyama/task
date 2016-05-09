module tx_ripple
  implicit none
  public

  integer(4) :: N_RPIN
  real(8), dimension(:), allocatable :: DltRP_rim, theta_rim
  real(8), dimension(:), allocatable :: RHOL, DltRP_rimL, theta_rimL, rip_ratL, DltRPL, &
       &                                DltRP_midL

contains

!**********************************************************************************
!
!   Toroidal ripple effects
!
!     Input  (real*8) : dQdrho (0:NRMAX) : rho-derivative of the safety factor
!     Output (real*8) : rNubrp1(0:NRMAX) : defraction frequency
!                       rNubrp2(0:NRMAX) : defraction frequency
!                       Ubrp   (0:NRMAX) : convective velocity due to grad B drift
!                       Dbrp   (0:NRMAX) : diffusivity
!
!**********************************************************************************

  subroutine ripple_effect(dQdrho)

    use libell, only : ellfc, ellec
    use tx_commons
    use tx_interface, only : fgaussian

    real(8), dimension(0:NRMAX), intent(in) :: dQdrho

    integer(4) :: NR, IER, i, imax, irip, nr_potato
    real(8) :: thetab, sinthb, rhol, ellE, ellK, EpsL, theta1, theta2, dlt, &
         &     width0, width1, dltwidth, ARC, diff_min, theta_min, &
         &     rhobl, rNueff, rNubnc, DRP, Dltcr, DltR, Vdrift, Rpotato, dQdrL
!!    real(8) :: sum_rp, DltRP_ave, DCB, Dlteff
    real(8) :: AITKEN2P
    real(8), dimension(0:NRMAX) :: th1, th2
!!rp_conv         &                         ,PNbrpL, DERIV
!!rp_conv    real(8), dimension(1:4,0:NRMAX) :: U

    if(allocated(DltRP_rim) .EQV. .FALSE.) allocate(DltRP_rim(0:NRMAX),theta_rim(0:NRMAX))

    ! Ripple amplitude
    thetab = 0.5D0 * PI ! pitch angle of a typical banana particle
    sinthb = sin(thetab)
    if(IRPIN == 0) then
       DO NR = 0, NRMAX
          rhol = rho(NR) * (1.D0 + (kappa - 1.D0) * sin(thetab))
          DltRP(NR) = ripple(rhol,thetab,FSRP) ! Ripple amplitude at banana tip point
                                             ! Mainly use for estimation of diffusive processes
       END DO
    end if

    ellK = ELLFC(sin(0.5d0*thetab),IER) ! first kind of complete elliptic function 
    ellE = ELLEC(sin(0.5d0*thetab),IER) ! second kind of complete elliptic function 

!!rp_conv    PNbrpLV(0:NRMAX) = 0.D0

        !     ***** Ripple loss transport *****

    rip_rat(0:NRMAX) = 0.D0
    IF(ABS(FSRP) > 0.D0) THEN
       ! +++ Convective loss +++

       ! *** Physical aspect of trapped particles in local ripple wells ***
       ! As pointed out in Goldston and Towner (1981) in section 2.5, banana particles
       ! which have small perpendicular component of their velocity, i.e. those
       ! residing near the banana tip point, are apt to be trapped in ripple wells,
       ! provided that they lie in ripple well region with alpha < 1.
       ! In this code, however, we cannot calculate banana tip point for each particle.
       ! Then we assume that all the ripple trapped particles are generated from the
       ! banana particles at the rim of the ripple well region, i.e. alpha = 1.
       ! We readily find two rims in a upper-half plane of a flux surface: LFS and HFS.
       ! Since the contribution from the HFS is usually negligible of their smallness,
       ! we only take account of the rim at the LFS.
       ! The banana particles have an opportunity to be trapped in local ripple wells
       ! when their banana tips lie in ripple well region. As will be discussed in 
       ! "Diffusive loss", we assume the poloidal angle of banana tip is 90 degrees on
       ! behalf of all the banana particles. However, on most of magnetic surfaces,
       ! the banana tips of 90 degrees do not lie in ripple well region, hence we cannot
       ! assume it for ripple well trapping. In this case we assume that the banana
       ! particles have an maxwellian distribution over the poloidal plane whose peak
       ! is taken at theta=90 degree, which is taken to be consistent with the assumption
       ! for diffusive loss case. Therefore, when we consider the scattering at the rim
       ! of ripple well region, we estimate the banana particle density which are going
       ! to be trapped as {sqrt(delta(theta_rim)) * fgaussian(theta_rim)} * nb.

       ! Ripple well region, i.e. area where "width0 - width1 > 0".
       DO NR = 1, NRMAX
          EpsL = epst(NR)

          ! ===== Calculation of ripple field and ripple well region =====
          if(IRPIN == 0) then
             ! For LFS
             theta1  = 0.d0
             i = 0
             imax = 101
             dlt = 1.d0 / (imax - 1)
             irip = 0
             do
                i = i + 1
                irip = irip + 1
                if(i == imax) then
!                   write(6,'(A,I3)') "LFS rim of ripple well region not detected at NR = ",NR
                   theta1 = PI
                   exit
                end if
                theta1 = theta1 + PI * dlt
                rho = rho(NR) * (1.D0 + (kappa - 1.D0) * sin(theta1))
                width0 = ripple(rhol,theta1,FSRP)
                width1 = EpsL * sin(theta1) / (NTCOIL * Q(NR))
                ! Poloidal angle at which the difference between width0 and width1 is minimized.
                dltwidth = abs(width0 - width1)
                if(i == 1 .or. dltwidth < diff_min) then
                   diff_min  = dltwidth
                   theta_min = theta1
                end if
                ! Convergence: Rim of ripple well region detected
                if(dltwidth < 1.d-6) exit
                ! Overreached a rim of ripple well. Go back and use finer step size.
                if(width0 < width1) then
                   theta1  = theta1 - PI * dlt
                   dlt = 0.1d0 * dlt
                   i = 0
                   irip = irip - 1
                   cycle
                end if
                ! Seek ripple amplitude just inside the rim of the ripple well region
                DltRP_rim(nr) = width0
                theta_rim(nr) = theta1
             end do
             ARC = 2.d0 * theta1
             th1(nr) = theta_min ! save for graphics

             ! For HFS
             theta2 = PI
             i = 0
             imax = 101
             dlt = 1.d0 / (imax - 1)
             irip = 0
             do 
                i = i + 1
                irip = irip + 1
                if(i == imax) then
!                   write(6,'(A,I3)') "HFS rim of ripple well region not detected at NR = ",NR
                   theta2 = PI
                   exit
                end if
                theta2 = theta2 - PI * dlt
                rhol = rho(NR) * (1.D0 + (kappa - 1.D0) * sin(theta2))
                width0 = ripple(rhol,theta2,FSRP)
                width1 = EpsL * sin(theta2) / (NTCOIL * Q(NR))
                dltwidth = abs(width0 - width1)
                if(i == 1 .or. dltwidth < diff_min) then
                   diff_min  = dltwidth
                   theta_min = theta2
                end if
                if(dltwidth < 1.d-6) exit
                if(width0 < width1) then
                   theta2  = theta2 + PI * dlt
                   dlt = 0.1d0 * dlt
                   i = 0
                   irip = irip - 1
                   cycle
                end if
             end do
             ARC = ARC + 2.d0 * (PI - theta2)
             ! Ratio of ripple well region in a certain flux surface
             rip_rat(NR) = ARC / (2.d0 * PI)
             th2(nr) = theta_min ! save for graphics
          end if

!!rpl_ave          sum_rp = 0.d0
!!rpl_ave          imax = 51
!!rpl_ave          dlt = 1.d0 / (imax - 1)
!!rpl_ave          do i = 1, imax
!!rpl_ave             theta = (i - 1) * PI * dlt
!!rpl_ave             rhol = rho(NR) * (1.D0 + (kappa - 1.D0) * sin(theta))
!!rpl_ave             sum_rp = sum_rp + ripple(rhol,theta,FSRP)**2
!!rpl_ave          end do
!!rpl_ave          DltRP_ave = sqrt(sum_rp/imax)

!!$          ! alpha_l : ripple well parameter
!!$          alpha_l = EpsL * sin(thetab) / (NTCOIL * Q(NR) * DltRP(NR))
!!$          ! Dlteff : effective depth of well along the magnetic field line
!!$          Dlteff = 2.D0*DltRP(NR)*(SQRT(1.D0-alpha_l**2)-alpha_l*acos(alpha_l))

          ! ===== Coefficients calculation (Part I) =====

          ! effective time of detrapping
          ! (Yushmanov NF (1982), Stringer NF (1972) 689, Takamura (5.31))
          if(DltRP_rim(nr) == 0.d0) then
             rNubrp1(NR) = 0.d0
          else
             rNubrp1(NR) = rNuD(NR) / DltRP_rim(nr)
          end if
          ! See the description of "Convective loss"
          rNubrp2(NR) = rNubrp1(NR) * SQRT(DltRP_rim(nr)) * fgaussian(theta_rim(nr),0.5D0*PI,0.85D0)
!!rpl_ave          rNubrp2(NR) = rNubrp1(NR) * SQRT(DltRP_ave)

          ! Convectitve loss (vertical grad B drift velocity)
          Vdrift = 0.5D0 * amb * amqp * Vb**2 / (achgb * RR * SQRT(BphV(NR)**2 + BthV(NR)**2))
!          RUbrp(NR)=(NTCOIL*Q(NR)*RR*DltRP(NR))*Vdrift
!          if(nr/=0) Ubrp(NR)=(NTCOIL*Q(NR)*RR*DltRP(NR))/R(NR)*Vdrift
          RUbrp(NR)=(NTCOIL*Q(NR)*RR*DltRP_rim(nr))*Vdrift
          if(nr/=0) Ubrp(NR)=(NTCOIL*Q(NR)*RR*DltRP_rim(nr))/R(NR)*Vdrift
!!$          Ubrp(NR) = 0.5D0 * Vdrift
!!$            &  * (theta1*sin(theta1) + (PI - theta2)*sin(theta2)) / (PI + theta1 - theta2))
!!$          IF(NR == NRMAX) THEN
!!$             rNubL(NR) = rNubL(NR-1)
!!$          ELSE
!!$             rNubL(NR) = Ubrp(NR) / SQRT(RB**2 - R(NR)**2)
!!$          END IF
!!$          rNubL(NR) = Ubrp(NR) / (R(NR) * sin(theta1))
!          if(nr >=5) stop
       END DO
!!$       RV0 = AITKEN2P(R(0),r(1)*(pi-th2(1)),r(1)*th1(1),r(2)*th1(2),-R(1),R(1),R(2))
!!$       rNubL(0) = Ubrp(0) / RV0
       Ubrp(0) = AITKEN2P(R(0),Ubrp(1),Ubrp(2),Ubrp(3),R(1),R(2),R(3))
       ! On the axis ripple amplitude is uniquely defined because of no poloidal variation.
       rNubrp1(0) = rNuD(0) / DltRP(0)
       rNubrp2(0) = rNubrp1(0) * SQRT(DltRP(0))

       ! Save for graphic
       if(IRPIN == 0) then
          thrp(1:nrmax) = th2(nrmax:1:-1)
          thrp(nrmax+1:2*nrmax) = th1(1:nrmax)
       end if

!!rp_conv       CALL SPL1D(R,PNbrpV,DERIV,U,NRMAX+1,0,IER)
!!rp_conv       do nr = 0, nrmax
!!rp_conv          if(rho(NR) <= rhob * cos(th1(nr))) then
!!rp_conv             tmp = r(nr)*cos(th1(nr))
!!rp_conv             call wherenr(r,tmp,nrl,Left,Right)
!!rp_conv             CALL SPL1DF(tmp,PNbrpL(NR),R,U,NRMAX+1,IER)!
!!rp_conv             PNbrpLV(NRL-1) = PNbrpLV(NRL-1) + Left  * PNbrpL(NR)
!!rp_conv             PNbrpLV(NRL)   = PNbrpLV(NRL)   + Right * PNbrpL(NR)
!!rp_conv             write(6,*) nrl,real(tmp),real(r(nrl-1)),real(r(nrl)),real(PNbrpL(NR)),real(Left  * PNbrpL(NR)),real(Right * PNbrpL(NR)),real(PNbrpLV(NRL-1)),real(PNbrpLV(NRL))
!!rp_conv          end if
!!rp_conv       end do

       !  +++ Diffusive loss +++

       ! *** Physical aspect of banana particles suffered from diffusive loss ***
       ! Banana particles are not affected by local ripple wells all the way to
       ! the bounce motion except banana tip points, even if they lie in a region
       ! with ripple wells. All the diffusion processes for them occur at the
       ! banana tip point "thetab", hence we only consider the ripple amplitude at
       ! the banana tip point, DltRP(NR).

       ! ===== Coefficients calculation (Part II) =====

       !  -- Collisional diffusion of trapped fast particles --
       IF(PNBH == 0.D0) THEN
          Dbrp(0:NRMAX) = 0.D0
       ELSE
       do nr = 1, nrmax
          EpsL = epst(NR)
          dQdrL = dQdrho(NR) / RA

          ! rhobl : Larmor radius of beam ions
          rhobl = amb * amqp * Vb / (achgb * SQRT(BphV(NR)**2 + BthV(NR)**2))
          ! DltR : Step size of banana particles
!          DltR = SQRT(PI * NTCOIL * (Q(NR) / EpsL)**3 * DltRP(NR)**2 * rhobl**2)
          DltR = rhobl*DltRP(NR)*SQRT(PI*NTCOIL*Q(nr)**3)*ABS(sinthb) &
               & /((NTCOIL*Q(NR)*DltRP(NR))**1.5D0+(EpsL*ABS(sinthb))**1.5D0)

          ! effective collisional frequency (G. Park et al., PoP 10 (2003) 4004, eq.(6))
!          rNueff = 1.82d0 * Q(NR)**2 * NTCOIL**2 / EpsL * rNuD(NR) ! for thetab=0.5*PI
          rNueff = 8.D0 * Q(NR)**2 * NTCOIL**2 * rNuD(NR) / (EpsL * sinthb**2) &
               & * (ellE/ellK - cos(0.5d0*thetab)**2)
          ! rNubnc : bounce frequency of beam ions (Helander and Sigmar, p132 eq.(7.27))
!          rNubnc = SQRT(EpsL) * Vb / (10.5D0 * Q(NR) * (RR + R(NR))) ! for thetab=0.5*PI
          rNubnc = SQRT(EpsL) * Vb / (4.D0 * SQRT(2.D0) * ellK * Q(NR) * RR)
!!$          ! DCB : confined banana diffusion coefficient
!!$          ! (V. Ya Goloborod'ko, et al., Physica Scripta T16 (1987) 46)
!!$          DCB = NTCOIL**2.25D0*Q(NR)**3.25D0*RR*rhobl*DltRP(NR)**1.5d0*rNuD(NR)/EpsL**2.5D0
          ! DRP : ripple-plateau diffusion coefficient
          DRP = rNubnc * DltR**2

          ! Collisional ripple well diffusion
!!$          if (rNueff < rNubnc) then
!!$             Dbrp(NR) = DCB
!!$          else
!!$             Dbrp(NR) = DRP
!!$          end if
!!$!          Dbrp(NR) = DCB * DRP / (DCB + DRP)
          Dbrp(NR) = DRP/SQRT(1.D0+(rNubnc/rNueff*DltR*NTCOIL*dQdrL)**2)

          ! Dltcr : GWB criterion of stochastic diffusion at banana tip point
!          Dltcr = (EpsL / (PI * NTCOIL * Q(NR)))**1.5D0 / (rhobl * dQdrL)
          Dltcr = DltRP(NR) / (DltR * NTCOIL * (2.D0 * thetab * dQdrL &
          &     + 2.D0 * Q(NR) / R(NR) * cos(thetab) / sinthb))
!!$          ! Fraction of stochastic region occupied in a flux surface
!!$          theta1 = 0.d0
!!$          i = 0
!!$          dlt = 1.d0 / (imax - 1)
!!$          ist = -1
!!$          do
!!$             i = i + 1
!!$             rhol = rho(NR) * (1.D0 + (kappa - 1.D0) * sin(theta1))
!!$             if(ripple(rhol,theta1,fsrp) > Dltcr) ist = ist + 1
!!$             theta1 = theta1 + PI * dlt
!!$             if(i == imax) exit
!!$          end do
!!$          ! Collisionless stochastic (ergodic) diffusion (whose value is 
!!$          ! equivalent to that of ripple-plateau diffusion)
!!$          if (DltRP(NR) > Dltcr) then
!!$             ! facST : Fraction of stochastic region in a flux surface
!!$             facST = dble(ist) / dble(imax - 1)
!!$             Dbrp(NR) = DRP * facST + Dbrp(NR) * (1.d0 - facST)
!!$          end if
          if (DltRP(NR) > Dltcr) Dbrp(NR) = DRP

!!$          ! Old version for display and comparison
!!$          DRP = rNubnc * PI * NTCOIL * (Q(NR) / EpsL)**3 * DltRP(NR)**2 * rhobl**2
!!$          DCB = NTCOIL**2.25D0*Q(NR)**3.25D0*RR*rhobl*DltRP(NR)**1.5d0*rNuD(NR)/EpsL**2.5D0
!!$          Dbrp(NR) = DCB * DRP / (DCB + DRP)
!!$          if(DltRP(NR) > Dltcr) Dbrp(NR) = DRP
!!$          write(6,*) r(nr)/ra,ft(NR)*rip_rat(NR)*Dbrp(NR)
       end do
       Dbrp(0) = AITKEN2P(R(0),Dbrp(1),Dbrp(2),Dbrp(3),R(1),R(2),R(3))

       ! Potato orbit effect
       !   Approaching to the magnetic axis, we reach the point where the banana
       !   particles changes themselves to the potato particles.
       !   Inside the point, the particle orbit is no longer changed.
       do nr = 0, nrmax
          ! potato width (Helander and Sigmar, p133)
          Rpotato = (Q(NR)**2*(amb * amqp * Vb / (achgb * BphV(NR)))**2*RR)**(1.D0/3.D0)
          if(r(nr) > Rpotato) then
             nr_potato = nr - 1
             exit
          end if
       end do
       do nr = 0, nr_potato
          Dbrp(nr)  = Dbrp(nr_potato+1)
       end do
       END IF

    ELSE
       rNubrp1(0:NRMAX) = 0.D0
       rNubrp2(0:NRMAX) = 0.D0
       Ubrp(0:NRMAX) = 0.D0
       Dbrp(0:NRMAX) = 0.D0
    END IF

  end subroutine ripple_effect

!***************************************************************
!
!   Ripple amplitude function
!     (Yushmanov review pp.122, 123)
!
!***************************************************************

  real(8) function ripple(rhol,theta,FSRP) result(f)
    use libbes, only : BESIN
    use tx_commons, only : RR, NTCOIL, DltRPn, RA
    real(8), intent(in) :: rhol, theta, FSRP
    real(8) :: a, L0, rl, Rmag0 = 2.4D0 ! specific value for JT-60U

    if(FSRP /= 0.D0) then
       rl = rhol * RA
       L0 = RR - Rmag0
       a = sqrt((RL**2+L0**2+2.D0*RL*L0*cos(theta))*(RR-L0)/(RR+RL*cos(theta)))
       f = DltRPn * BESIN(0,NTCOIL/(RR-L0)*a)
    else
       f = 0.d0
    end if

  end function ripple

!!rp_conv  ! Search minimum radial number NR satisfying R(NR) > X.
!!rp_conv
!!rp_conv  subroutine wherenr(R,X,NR,Left,Right)
!!rp_conv    real(8), dimension(0:NRMAX), intent(in) :: R
!!rp_conv    real(8), intent(in) :: X
!!rp_conv    integer, intent(out) :: NR
!!rp_conv    real(8), intent(out) :: Left, Right
!!rp_conv    integer :: NRL
!!rp_conv
!!rp_conv    if(X == 0.d0) then
!!rp_conv       NR = 1
!!rp_conv       Left  = 0.d0
!!rp_conv       Right = 0.d0
!!rp_conv       return
!!rp_conv    end if
!!rp_conv
!!rp_conv    do nrl = 1, nrmax
!!rp_conv       if(r(nrl) > x) then
!!rp_conv          NR = nrl
!!rp_conv          Right = (x - r(nr-1)) / (r(nr) - r(nr-1))
!!rp_conv          Left  = 1.d0 - Right
!!rp_conv          exit
!!rp_conv       end if
!!rp_conv    end do
!!rp_conv
!!rp_conv  end subroutine wherenr

!***************************************************************
!
!   Read ripple data file from TASK/EQ
!
!***************************************************************

  subroutine ripple_input(ier)

    integer(4) :: ier, nrpin, ist, nrpmax, n, i, j
    character(len=130) :: kline

    !  *** Read ripple data from file ***

    nrpin = 23
    call FROPEN(nrpin,'tx_ripple.dat',1,0,'RIPPLE FROM EQ',ier)
    if(ier /= 0) return

    read(nrpin,'(I4)',iostat=ist) nrpmax
    if(ist /= 0) return
    read(nrpin,'(A130)',iostat=ist) kline
    if(ist /= 0) return

    N = nrpmax
    allocate(RHOL(1:N), DltRP_rimL(1:N), theta_rimL(1:N), rip_ratL(1:N), DltRPL(1:N), &
         &   DltRP_midL(1:N))

    do i = 1, N
       read(nrpin,'(I4,6E15.7)',iostat=ist) j, RHOL(i), DltRP_rimL(i), theta_rimL(i), &
            &                                  rip_ratL(i), DltRPL(i), DltRP_midL(i)
       if(ist /= 0) return
!!$       write(6,'(I4,1P6E15.7)') i,RHOL(i), DltRP_rimL(i), theta_rimL(i), &
!!$            &                                  rip_ratL(i), DltRPL(i), DltRP_midL(i)
    end do

    close(nrpin)

    N_RPIN = nrpmax

  end subroutine ripple_input

!***************************************************************
!
!   Spline interpolated ripple data
!
!***************************************************************

  subroutine ripple_spl

    use tx_commons, only : NRA, Rho, rip_rat, DltRP, DltRP_mid, DltRPn, NRMAX

    integer(4) :: nr, ier, N
    real(8), dimension(:), allocatable :: DERIV
    real(8), dimension(:,:), allocatable :: U
    real(8) :: deriv4

    N = N_RPIN

    allocate(DERIV(1:N),U(1:4,1:N)) 
    if(allocated(DltRP_rim) .EQV. .FALSE.) allocate(DltRP_rim(0:NRMAX),theta_rim(0:NRMAX))

    !   *** Interpolation ***

    DERIV(1) = deriv4(1,RHOL,DltRP_rimL,N,1)
    DERIV(N) = deriv4(N,RHOL,DltRP_rimL,N,1)
    call spl1d(RHOL,DltRP_rimL,DERIV,U,N,3,ier)
    if(ier /= 0) stop 'Error at spl1d in ripple_input: DltRP_rim'
    do nr = 0, nra
       call spl1df(Rho(nr),DltRP_rim(nr),RHOL,U,N,ier)
       if(ier /= 0) stop 'Error at spl1df in ripple_input: DltRP_rim'
    end do
    where(DltRP_rim < 0.d0) DltRP_rim = 0.d0

    DERIV(1) = deriv4(1,RHOL,theta_rimL,N,1)
    DERIV(N) = deriv4(N,RHOL,theta_rimL,N,1)
    call spl1d(RHOL,theta_rimL,DERIV,U,N,3,ier)
    if(ier /= 0) stop 'Error at spl1d in ripple_input: theta_rim'
    do nr = 0, nra
       call spl1df(Rho(nr),theta_rim(nr),RHOL,U,N,ier)
       if(ier /= 0) stop 'Error at spl1df in ripple_input: theta_rim'
    end do
    where(theta_rim < 0.d0) theta_rim = 0.d0
    
    DERIV(1) = deriv4(1,RHOL,rip_ratL,N,1)
    DERIV(N) = deriv4(N,RHOL,rip_ratL,N,1)
    call spl1d(RHOL,rip_ratL,DERIV,U,N,3,ier)
    if(ier /= 0) stop 'Error at spl1d in ripple_input: rip_rat'
    do nr = 0, nra
       call spl1df(Rho(nr),rip_rat(nr),RHOL,U,N,ier)
       if(ier /= 0) stop 'Error at spl1df in ripple_input: rip_rat'
    end do
    where(rip_rat < 0.d0) rip_rat = 0.d0

    DERIV(1) = deriv4(1,RHOL,DltRPL,N,1)
    DERIV(N) = deriv4(N,RHOL,DltRPL,N,1)
    call spl1d(RHOL,DltRPL,DERIV,U,N,3,ier)
    if(ier /= 0) stop 'Error at spl1d in ripple_input: DltRPL'
    do nr = 0, nra
       call spl1df(Rho(nr),DltRP(nr),RHOL,U,N,ier)
       if(ier /= 0) stop 'Error at spl1df in ripple_input: DltRPL'
    end do

    DERIV(1) = deriv4(1,RHOL,DltRP_midL,N,1)
    DERIV(N) = deriv4(N,RHOL,DltRP_midL,N,1)
    call spl1d(RHOL,DltRP_midL,DERIV,U,N,3,ier)
    if(ier /= 0) stop 'Error at spl1d in ripple_input: DltRP_midL'
    do nr = 0, nra
       call spl1df(Rho(nr),DltRP_mid(nr),RHOL,U,N,ier)
       if(ier /= 0) stop 'Error at spl1df in ripple_input: DltRP_midL'
    end do

    DltRPn = DltRP_midL(1)

    !   *** Linear extrapolation ***

    do nr = nra+1, nrmax
       call aitken(Rho(nr),DltRP_rim(nr),Rho,DltRP_rim,1,nr)
       call aitken(Rho(nr),theta_rim(nr),Rho,theta_rim,1,nr)
       call aitken(Rho(nr),rip_rat(nr)  ,Rho,rip_rat  ,1,nr)
       call aitken(Rho(nr),DltRP(nr)    ,Rho,DltRP    ,1,nr)
       call aitken(Rho(nr),DltRP_mid(nr),Rho,DltRP_mid,1,nr)
    end do

!!$    do nr = 0, nrmax
!!$       write(6,'(I4,1P5E15.7)') nr,Rho(nr), DltRP_rim(nr), theta_rim(nr), &
!!$            &                                  rip_rat(nr), DltRP(nr)       
!!$    end do
!!$    stop
    
    deallocate(DERIV,U)

  end subroutine ripple_spl

!***************************************************************
!
!   Deallocate ripple data
!
!***************************************************************

  subroutine deallocate_ripple

    if(allocated(RHOL)) deallocate(RHOL, DltRP_rimL, theta_rimL, rip_ratL, DltRPL, DltRP_midL)
    if(allocated(DltRP_rim)) deallocate(DltRP_rim,theta_rim)

  end subroutine deallocate_ripple

end module tx_ripple
