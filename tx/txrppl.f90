module tx_ripple
  implicit none
  public

contains

!**********************************************************************************
!
!   Toroidal ripple effects
!
!     Inputs (real*8) : dQdr   (0:NRMAX) : r-derivative of the safety factor
!     Output (real*8) : rNubrp1(0:NRMAX) : defraction frequency
!                       rNubrp2(0:NRMAX) : defraction frequence
!                       Ubrp   (0:NRMAX) : convective velocity due to grad B drift
!                       Dbrp   (0:NRMAX) : diffusivity
!
!**********************************************************************************

  subroutine ripple_effect(dQdr)

    use tx_commons
    use tx_interface, only : fgaussian

    real(8), dimension(0:NRMAX), intent(in) :: dQdr

    integer(4) :: NR, IER, i, imax, irip, nr_potato
    real(8) :: thetab, sinthb, RL, ellE, ellK, EpsL, theta1, theta2, dlt, width0, width1, &
         &     ARC, DltRP_rim, theta_rim, diff_min, theta_min, sum_rp, DltRP_ave, &
         &     rhob, rNueff, rNubnc, DCB, DRP, Dltcr, Dlteff, DltR, Vdrift, Rpotato
    real(8) :: ELLFC, ELLEC, AITKEN2P
    real(8), dimension(0:NRMAX) :: th1, th2
!!rp_conv         &                         ,PNbrpL, DERIV
!!rp_conv    real(8), dimension(1:4,0:NRMAX) :: U

    ! Ripple amplitude
    thetab = 0.5D0 * PI ! pitch angle of a typical banana particle
    sinthb = sin(thetab)
    DO NR = 0, NRMAX
       RL = R(NR) * (1.D0 + (kappa - 1.D0) * sin(thetab))
       DltRP(NR) = ripple(RL,thetab,FSRP) ! Ripple amplitude at banana tip point
                                          ! Mainly use for estimation of diffusive processes
    END DO

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

       ! Ripple well region
       DO NR = 1, NRMAX
          RL = R(NR)
          EpsL = RL / RR
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
!                write(6,'(A,I3)') "LFS rim of ripple well region not detected at NR = ",NR
                theta1 = PI
                exit
             end if
             theta1 = theta1 + PI * dlt
             RL = R(NR) * (1.D0 + (kappa - 1.D0) * sin(theta1))
             width0 = ripple(RL,theta1,FSRP)
             width1 = EpsL * sin(theta1) / (NTCOIL * Q(NR))
             ! Poloidal angle at which the difference between width0 and width1 is minimized.
             if(i == 1) then
                theta_min = theta1
                diff_min = abs(width0 - width1)
             else
                if(abs(width0 - width1) < diff_min) then
                   diff_min = abs(width0 - width1)
                   theta_min = theta1
                end if
             end if
             ! Rim of ripple well region detected
             if(abs(width0 - width1) < 1.d-6) exit
             ! Overreached a rim of ripple well. Go back and use finer step size.
             if(width0 < width1) then
                theta1  = theta1 - PI * dlt
                dlt = 0.1d0 * dlt
                i = 0
                irip = irip - 1
                cycle
             end if
             ! Ripple amplitude at the rim of the ripple well region
             DltRP_rim = width0
             theta_rim = theta1
          end do
          ARC = 2.d0 * theta1
          th1(nr) = theta_min

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
!                write(6,'(A,I3)') "HFS rim of ripple well region not detected at NR = ",NR
                theta2 = PI
                exit
             end if
             theta2 = theta2 - PI * dlt
             RL = R(NR) * (1.D0 + (kappa - 1.D0) * sin(theta2))
             width0 = ripple(RL,theta2,FSRP)
             width1 = EpsL * sin(theta2) / (NTCOIL * Q(NR))
             if(i == 1) then
                theta_min = theta2
                diff_min = abs(width0 - width1)
             else
                if(abs(width0 - width1) < diff_min) then
                   diff_min = abs(width0 - width1)
                   theta_min = theta2
                end if
             end if
             if(abs(width0 - width1) < 1.d-6) exit
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
          th2(nr) = theta_min

!!rpl_ave          sum_rp = 0.d0
!!rpl_ave          imax = 51
!!rpl_ave          dlt = 1.d0 / (imax - 1)
!!rpl_ave          do i = 1, imax
!!rpl_ave             theta = (i - 1) * PI * dlt
!!rpl_ave             RL = R(NR) * (1.D0 + (kappa - 1.D0) * sin(theta))
!!rpl_ave             sum_rp = sum_rp + ripple(RL,theta,FSRP)**2
!!rpl_ave          end do
!!rpl_ave          DltRP_ave = sqrt(sum_rp/imax)

!!$          ! alpha_l : ripple well parameter
!!$          alpha_l = EpsL * sin(thetab) / (NTCOIL * Q(NR) * DltRP(NR))
!!$          ! Dlteff : effective depth of well along the magnetic field line
!!$          Dlteff = 2.D0*DltRP(NR)*(SQRT(1.D0-alpha_l**2)-alpha_l*acos(alpha_l))

          ! effective time of detrapping
          ! (Yushmanov NF (1982), Stringer NF (1972) 689, Takamura (5.31))
          rNubrp1(NR) = rNuD(NR) / DltRP_rim
          ! See the description of "Convective loss"
          rNubrp2(NR) = rNubrp1(NR) * SQRT(DltRP_rim) * fgaussian(theta_rim,0.5D0*PI,0.85D0)
!!rpl_ave          rNubrp2(NR) = rNubrp1(NR) * SQRT(DltRP_ave)

          ! Convectitve loss (vertical grad B drift velocity)
          Vdrift = 0.5D0 * AMb * Vb**2 / (PZ * AEE * RR * SQRT(BphV(NR)**2 + BthV(NR)**2))
!          RUbrp(NR)=(NTCOIL*Q(NR)*RR*DltRP(NR))*Vdrift
!          if(nr/=0) Ubrp(NR)=(NTCOIL*Q(NR)*RR*DltRP(NR))/R(NR)*Vdrift
          RUbrp(NR)=(NTCOIL*Q(NR)*RR*DltRP_rim)*Vdrift
          if(nr/=0) Ubrp(NR)=(NTCOIL*Q(NR)*RR*DltRP_rim)/R(NR)*Vdrift
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
       thrp(1:nrmax) = th2(nrmax:1:-1)
       thrp(nrmax+1:2*nrmax) = th1(1:nrmax)

!!rp_conv       CALL SPL1D(R,PNbrpV,DERIV,U,NRMAX+1,0,IER)
!!rp_conv       do nr = 0, nrmax
!!rp_conv          if(R(NR) <= RB * cos(th1(nr))) then
!!rp_conv             tmp = r(nr)*cos(th1(nr))
!!rp_conv             call wherenr(r,tmp,nrl,Left,Right)
!!rp_conv             CALL SPL1DF(tmp,PNbrpL(NR),R,U,NRMAX+1,IER)!
!!rp_conv             PNbrpLV(NRL-1) = PNbrpLV(NRL-1) + Left  * PNbrpL(NR)
!!rp_conv             PNbrpLV(NRL)   = PNbrpLV(NRL)   + Right * PNbrpL(NR)
!!rp_conv             write(6,*) nrl,sngl(tmp),sngl(r(nrl-1)),sngl(r(nrl)),sngl(PNbrpL(NR)),sngl(Left  * PNbrpL(NR)),sngl(Right * PNbrpL(NR)),sngl(PNbrpLV(NRL-1)),sngl(PNbrpLV(NRL))
!!rp_conv          end if
!!rp_conv       end do

       !  +++ Diffusive loss +++

       ! *** Physical aspect of banana particles suffered from diffusive loss ***
       ! Banana particles are not affected by local ripple wells all the way to
       ! the bounce motion except banana tip points, even if they lie in a region
       ! with ripple wells. All the diffusion processes for them occur at the
       ! banana tip point "thetab", hence we only consider the ripple amplitude at
       ! the banana tip point, DltRP(NR).

       !  -- Collisional diffusion of trapped fast particles --
       IF(PNBH == 0.D0) THEN
          Dbrp(0:NRMAX) = 0.D0
       ELSE
       do nr = 1, nrmax
          RL = R(NR)
          EpsL = RL / RR

          ! rhob : Larmor radius of beam ions
          rhob = AMb * Vb / (PZ * AEE * SQRT(BphV(NR)**2 + BthV(NR)**2))
          ! DltR : Step size of banana particles
!          DltR = SQRT(PI * NTCOIL * (Q(NR) / EpsL)**3 * DltRP(NR)**2 * rhob**2)
          DltR = rhob*DltRP(NR)*SQRT(PI*NTCOIL*Q(nr)**3)*ABS(sinthb) &
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
!!$          DCB = NTCOIL**2.25D0*Q(NR)**3.25D0*RR*rhob*DltRP(NR)**1.5d0*rNuD(NR)/EpsL**2.5D0
          ! DRP : ripple-plateau diffusion coefficient
          DRP = rNubnc * DltR**2

          ! Collisional ripple well diffusion
!!$          if (rNueff < rNubnc) then
!!$             Dbrp(NR) = DCB
!!$          else
!!$             Dbrp(NR) = DRP
!!$          end if
!!$!          Dbrp(NR) = DCB * DRP / (DCB + DRP)
          Dbrp(NR) = DRP/SQRT(1.D0+(rNubnc/rNueff*DltR*NTCOIL*dQdr(NR))**2)

          ! Dltcr : GWB criterion of stochastic diffusion at banana tip point
!          Dltcr = (EpsL / (PI * NTCOIL * Q(NR)))**1.5D0 / (rhob * dQdr(NR))
          Dltcr = DltRP(NR) / (DltR * NTCOIL * (2.D0 * thetab * dQdr(NR) &
          &     + 2.D0 * Q(NR) / R(NR) * cos(thetab) / sinthb))
!!$          ! Fraction of stochastic region occupied in a flux surface
!!$          theta1 = 0.d0
!!$          i = 0
!!$          dlt = 1.d0 / (imax - 1)
!!$          ist = -1
!!$          do
!!$             i = i + 1
!!$             RL = R(NR) * (1.D0 + (kappa - 1.D0) * sin(theta1))
!!$             if(ripple(RL,theta1,fsrp) > Dltcr) ist = ist + 1
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
!!$          DRP = rNubnc * PI * NTCOIL * (Q(NR) / EpsL)**3 * DltRP(NR)**2 * rhob**2
!!$          DCB = NTCOIL**2.25D0*Q(NR)**3.25D0*RR*rhob*DltRP(NR)**1.5d0*rNuD(NR)/EpsL**2.5D0
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
          Rpotato = (Q(NR)**2*(AMb * Vb / (PZ * AEE * BphV(NR)))**2*RR)**(1.D0/3.D0)
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

  real(8) function ripple(RL,theta,FSRP) result(f)
    use tx_commons, only : RR, NTCOIL, DltRPn, RA
    real(8), intent(in) :: RL, theta, FSRP
    real(8) :: a, L0, Rmag0 = 2.4D0 ! specific value for JT-60U
    real(8) :: BESIN

    if(FSRP /= 0.D0) then
       L0 = RR - Rmag0
       a = sqrt((RL**2+L0**2+2.D0*RL*L0*cos(theta))*(RR-L0)/(RR+RL*cos(theta)))
       f = DltRPn * BESIN(0,NTCOIL/(RR-L0)*a)
    else
       f = 0.d0
    end if

  end function ripple

!!$  real(8) function ripple(NR,theta,FSRP) result(f)
!!$    use tx_commons, only : RR, R, RA, NTCOIL
!!$    integer(4), intent(in) :: NR
!!$    real(8), intent(in) :: theta, FSRP
!!$    real(8) :: DIN = 0.2D0, DltRP0 = 0.015D0
!!$
!!$    if(FSRP /= 0.D0) then
!!$       f = DltRP0 * (       ((RR + R(NR) * cos(theta)) / (RR + RA))**(NTCOIL-1) &
!!$            &        + DIN *((RR - RA) / (RR + R(NR) * cos(theta)))**(NTCOIL+1))
!!$    else
!!$       f = 0.D0
!!$    end if
!!$       
!!$  end function ripple

end module tx_ripple
