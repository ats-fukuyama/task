!***************************************************************
!
!   Metrics associated with equilibrium
!
!***************************************************************

subroutine txequ
  use tx_commons, only : ieqread, epst, aat, rrt, ckt, suft, sst, vro, vlt, art, NRMAX, &
       & rr, Pisq, ra, rho, bb, d_rrr, bit, bbrt, Rax, Zax, surflcfs, Pi, elip, trig
  use eqread_mod, only : intequ, txmesh

  implicit none
  integer(4) :: nr

  ! surflcfs ; surface area of a torus enclosed by Last Closed Flux Surface
  ! epst = inverse aspect ratio
  ! aat  = <R^-2>
  ! rrt  = <R^2>
  ! ckt  = <B_p^2>(dV/dpsi)^2=<|nabla V|^2/R^2>
  ! suft = <|nabla V|>
  ! sst  = <|nabla V|^2>=<R^2B_p^2>(dV/dpsi)^2
  ! vro  = dV/drho
  ! vlt  = volume
  ! art  = area (cross section area)
  ! d_rrr = <1/R> ;  <1/R> = 2 Pi dS/dV
  ! bit  = <B^-2> ; for NCLASS ; This variable does not affect the estimate of friction and 
  !                              viscosity coefficients because it is used solely for the
  !                              estimate of the fluxes.
  ! bbrt = <B> ; for NCLASS
  ! elip = elongation
  ! trig = triangularity

  select case(ieqread)
  case(0) ! Large aspect ratio limit, circular cross section
     Rax = rr
     Zax = 0.d0
     surflcfs = 4.d0 * Pisq * RR * RA
     do NR = 0, NRMAX
        epst(NR)  = (rho(NR) * ra) / rr
        aat(NR)   = 1.d0 / rr**2
        rrt(NR)   = rr**2
        ckt(NR)   = 16.d0 * Pisq**2 * rr**2 * epst(NR)**2
        suft(NR)  = 4.d0 * Pisq * rr**2 * epst(NR)
        sst(NR)   = suft(NR)**2
        vro(NR)   = 4.d0 * Pisq * rr * ra**2 * rho(NR)
        vlt(NR)   = 2.d0 * Pisq * rr * (ra * rho(NR))**2
        art(NR)   = Pi * ra**2 * rho(NR)**2
        d_rrr(NR) = 1.d0 / rr
        bit(NR)   = 1.d0 / bb**2
        bbrt(NR)  = bb
        elip(NR)  = 1.d0
        trig(NR)  = 0.d0
     end do

     call txmesh
  case(1) ! Large aspect ratio approximation, circular cross section
     Rax = rr
     Zax = 0.d0
     surflcfs = 4.d0 * Pisq * RR * RA
     do NR = 0, NRMAX
        epst(NR)  = (rho(NR) * ra) / rr
        aat(NR)   = 1.d0 / ( rr**2 * sqrt(1.d0 - epst(NR)**2))
        rrt(NR)   = (1.d0 + 1.5d0 * epst(NR)**2) * rr**2
        ckt(NR)   = 16.d0 * Pisq**2 * rr**2 * epst(NR)**2 / sqrt(1.d0 - epst(NR)**2)
        suft(NR)  = 4.d0 * Pisq * rr**2 * epst(NR)
        sst(NR)   = suft(NR)**2
        vro(NR)   = 4.d0 * Pisq * rr * ra**2 * rho(NR)
        vlt(NR)   = 2.d0 * Pisq * rr * (ra * rho(NR))**2
        art(NR)   = Pi * ra**2 * rho(NR)**2
        d_rrr(NR) = 1.d0 / rr
        bit(NR)   = (1.d0 + 1.5d0 * epst(NR)**2) / bb**2
        bbrt(NR)  = bb
        elip(NR)  = 1.d0
        trig(NR)  = 0.d0
     end do

     call txmesh
  case(2) ! read eqdata once
     call intequ
  case default
     stop 'Under construction!'
  end select

end subroutine txequ
