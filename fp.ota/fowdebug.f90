module fowdebug
  
contains
  subroutine check_Jacobian
    use fowcomm
    use fpcomm

    implicit none 
    integer :: np, nth, nr, nsa
    real(rkind) :: V_v, V_x, J_sum, deltaps, deltath, deltap, pl, PVmax, Vmax

    do nsa = 1, nsamax
      PVmax = SQRT(1.D0+THETA0(NSA)*PM(NPMAX,NSA)**2)
      Vmax = vc * sqrt(1.d0-1.d0/PVmax**2)

      V_v = 4.d0/3.d0*pi*Vmax**3
      
      V_x = (pi*ra**2)*(2.d0*pi*rr)

      J_sum = 0.d0
      do nr = 1, nrmax
        deltaps = psimg(nr+1)-psimg(nr)
        do np = 1, npmax
          deltap = (pg(np+1,nsa)-pg(np,nsa))*ptfp0(nsa)
          do nth = 1, nthmax
            deltath = thetamg(nth+1,np,nr,nsa)-thetamg(nth,np,nr,nsa)
            J_sum = J_sum + deltap*deltath*deltaps*Jacobian_I(nth,np,nr,nsa)
          end do
        end do
      end do
      write(*,*)"nsa = ",nsa
      write(*,*)"Jacobian,V_xv,Vmax",J_sum,V_x*V_v, Vmax,PVmax

    end do

  end subroutine

end module fowdebug