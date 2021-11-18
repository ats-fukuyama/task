module fowdistribution
  implicit none
  private
  public :: fI_Maxwellian
  public :: convert_fI_to_fu
  public :: moment_0th_order_COM, moment_2nd_order_COM
  public :: particle_flux, particle_flux_element
  public :: total_N, effective_diffusion_coefficient

contains

  subroutine fI_Maxwellian(fI)
    use fpcomm
    use fpsub
    use fowcomm
    use foworbit

    implicit none
    real(rkind),intent(out) :: &
         fI(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend)

    real(rkind) :: rtfd0l, ptfd0l, ex
    real(rkind) :: rnfd0l, rnfdl, rtfdl, normalize
    integer :: nth, np, nr, ns, nsa
    real(rkind) :: sumI, sumu, r0, psip0, costh0, sinth0
    real(rkind) :: B0, F0, dBdr0, dFdr0, dpsipdr0

    do nsa = nsastart,nsaend
       ns = ns_nsa(nsa)

       sumI = 0.d0
       sumu = 0.d0
       do nr = nrstartw, nrendwm
          do np = npstartw, npendwm
             do nth = 1, nthmax
  
                call mean_ra_quantities(orbit_m(nth,np,nr,nsa), &
                     r0, psip0, costh0, sinth0, B0, F0, dBdr0, dFdr0, dpsipdr0)
  
                rtfd0l = ( ptpr(ns)+2.d0*ptpp(ns) )/3.d0
                ptfd0l = SQRT( rtfd0l*1.d3*aee*amfp(ns) )
                rnfd0l = pn(ns)
              
                ! rnfdl=rn_temp(nr,ns)/rnfd0l
                ! rtfdl=rt_temp(nr,ns)
                if ( nr == 1 ) then
                   rnfdl = rn_temp(nr,ns)/rnfd0l &
                         +(rn_temp(nr+1,ns)/rnfd0l-rn_temp(nr,ns) &
                         /rnfd0l)/delr*(r0-rm(nr))
                   rtfdl = rt_temp(nr,ns) &
                         +(rt_temp(nr+1,ns)-rt_temp(nr,ns))/delr*(r0-rm(nr))
                else
                   rnfdl = rn_temp(nr,ns)/rnfd0l &
                         +(rn_temp(nr,ns)/rnfd0l-rn_temp(nr-1,ns) &
                         /rnfd0l)/delr*(r0-rm(nr))
                   rtfdl = rt_temp(nr,ns) &
                         +(rt_temp(nr,ns)-rt_temp(nr-1,ns))/delr*(r0-rm(nr))
                end if
            
                ex = pm(np,nsa)**2/(2.d0*rtfdl/rtfd0l)
                fI(nth,np,nr,nsa) = EXP( -1.d0*ex )! / JI(nth,np,nr,nsa)
        
             end do
          end do
       end do
  
       do nr = nrstartw, nrendwm
          sumI = 0.d0
          do np = npstartw, npendwm
             do nth = 1, nthmax
                sumI = sumI &
                     + fI(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)*delp(nsa) &
                       *delthm(nth,np,nr,nsa)
             end do
          end do
  
          normalize = rn_temp(nr,ns)/(rnfd0l*sumI)
          do np = npstartw, npendwm
             do nth = 1, nthmax
                fI(nth,np,nr,nsa) = fI(nth,np,nr,nsa)*normalize
             end do
          end do
  
       end do
  
    end do
    
  end subroutine fI_Maxwellian

  subroutine convert_fI_to_fu(fu_l, fI)

    use fpcomm
    use fowcomm
    USE libmpi
    implicit none
    real(rkind),intent(in) :: &
         fI(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend)
    real(rkind),intent(out) :: &
         fu_l(nthmax,npstartw:npendwm,nrstartw:nrendwm, &
              nthpmax,nsastart:nsaend)

    integer :: nth,np,nr,nthp,nsa
    integer :: ir, ith, irr, irl, ithr, ithl
    real(rkind) :: rhoml, thetaml, suml, sumI
    real(rkind) :: dxm, dxp, dym, dyp, f11, f12, f21, f22, dd, vr
    REAL(rkind) :: suml_tot,sumI_tot
    INTEGER :: suml_loc,sumI_loc

    do nsa = nsastart,nsaend

       suml = 0.d0
       do nthp = 1,nthpmax
          do nr = nrstartw, nrendwm
             do np = npstartw, npendwm
                do nth = 1, nthmax
                   if ( time_loss(nth,np,nr,nthp,nsa) /= 0 ) then ! loss region
                      fu_l(nth,np,nr,nthp,nsa) = 0.d0

                   else
                      thetaml = thetam_local(nth,np,nr,nthp,nsa)
                      rhoml   = rhom_local(nth,np,nr,nthp,nsa)
  
                      ! search the cell include thetaml, psiml
                      do ir = nrstartw, nrendwm-1
                         if ( rm(ir) <= rhoml .and. rhoml <= rm(ir+1) ) then
                            exit
                         end if
                      end do
  
                      do ith = 1, nthmax-1
                         if ( thetam(ith,np,ir,nsa) <= thetaml &
                              .and. thetaml <= thetam(ith+1,np,ir,nsa) ) then
                            exit
                         end if
                      end do

                      ! use bilinear intepolation
                      ithr = min( ith+1, nthmax )
                      ithl = ithr-1
                      irr  = min( ir+1, nrmax )
                      irl  = irr-1

                      f11 = fI(ithl,np,irl,nsa)*JI(ithl,np,irl,nsa) &
                           *delthm(ithl,np,irl,nsa)*delp(nsa)*delr

                      f12 = fI(ithl,np,irr,nsa)*JI(ithl,np,irr,nsa) &
                           *delthm(ithl,np,irr,nsa)*delp(nsa)*delr

                      f21 = fI(ithr,np,irl,nsa)*JI(ithr,np,irl,nsa) &
                           *delthm(ithr,np,irl,nsa)*delp(nsa)*delr

                      f22 = fI(ithr,np,irr,nsa)*JI(ithr,np,irr,nsa) &
                           *delthm(ithr,np,irr,nsa)*delp(nsa)*delr

                      dxm = thetaml - thetam(ithl,np,irl,nsa)
                      dxp = thetam(ithr,np,irl,nsa) - thetaml 
                      dym = rhoml - rm(irl)
                      dyp = rm(irr) - rhoml
                      dd = (delr+delr)/2.d0 &
                           *(delthm(ithl,np,irl,nsa)+delthm(ithr,np,irl,nsa)) &
                           /2.d0

                      fu_l(nth,np,nr,nthp,nsa) &
                           = ( dyp*dxp*f11 + dyp*dxm*f21 &
                             + dym*dxp*f12 + dym*dxm*f22 )/dd

                      vr = twopi*rr*(rg(nr+1)-rg(nr))*rm(nr) &
                           *twopi/dble(nthpmax)

                      fu_l(nth,np,nr,nthp,nsa) = fu_l(nth,np,nr,nthp,nsa) &
                           / (vr*volp(nth,np,nsa)) / RNFP0(NSA)*1.D20

                      suml = suml &
                           + vr*volp(nth,np,nsa)*fu_l(nth,np,nr,nthp,nsa)
  
                   end if
              
                end do
             end do
          end do
       end do

       sumI = 0.d0
       do nr = nrstartw, nrendwm
          do np = npstartw, npendwm
             do nth = 1, nthmax
                sumI = sumI + delr*delp(nsa)*delthm(nth,np,nr,nsa) &
                              *JI(nth,np,nr,nsa)*fI(nth,np,nr,nsa)
             end do
          end do
       end do

       CALL mtx_set_communicator(comm_nrnp)
       CALL mtx_allreduce1_real8(suml,3,suml_tot,suml_loc)
       CALL mtx_allreduce1_real8(sumI,3,sumI_tot,sumI_loc)
       CALL mtx_reset_communicator
       suml=suml_tot
       sumI=sumI_tot

       do nthp = 1, nthpmax
          do nr = nrstartw, nrendwm
             do np = npstartw, npendwm
                do nth = 1, nthmax
                   vr = twopi*rr*(rg(nr+1)-rg(nr))*rm(nr)*twopi/dble(nthpmax)
                   fu_l(nth,np,nr,nthp,nsa) &
                        = fu_l(nth,np,nr,nthp,nsa) * sumI / suml 
                end do
             end do
          end do
       end do

    end do

  end subroutine convert_fI_to_fu

  subroutine moment_0th_order_COM(M0, fI_in)
    use fpcomm
    use fowcomm
    use foworbit
    USE libmpi
    implicit none
    real(rkind),intent(in) :: &
         fI_in(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend)
    real(rkind),intent(out) :: M0(nrstartw:nrendwm,nsastart:nsaend)
    integer :: nth, np, nr, nsa
    real(rkind) :: dVI,sum,sum_tot
    INTEGER:: sum_loc

    CALL mtx_set_communicator(comm_np)
    do nsa = nsastart, nsaend
       do nr = nrstartw, nrendwm
          sum = 0.d0
          do np = npstartw, npendwm
             do nth = 1, nthmax
                dVI = delthm(nth,np,nr,nsa)*delp(nsa)
                sum = sum + fI_in(nth,np,nr,nsa) * dVI * JIR(nth,np,nr,nsa)
             end do
          end do
          CALL mtx_allreduce1_real8(sum,3,sum_tot,sum_loc)
          M0(nr,nsa) = sum_tot*rnfp0(nsa)
       end do
    end do
    CALL mtx_reset_communicator

  end subroutine moment_0th_order_COM

  subroutine moment_2nd_order_COM(M2, fI_in)
    use fpcomm
    use fowcomm
    use foworbit
    USE libmpi
    implicit none
    real(rkind),intent(in) :: &
         fI_in(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend)
    real(rkind),intent(out) :: M2(nrstartw:nrendwm,nsastart:nsaend)
    integer :: nth, np, nr, nsa
    real(rkind) :: dVI, K, PV, sum, sum_tot
    INTEGER:: sum_loc

    CALL mtx_set_communicator(comm_np)
    do nsa = nsastart, nsaend
       do nr = nrstartw, nrendwm
          sum = 0.d0
          do np = npstartw, npendwm
             PV = sqrt(1.d0+theta0(nsa)*pm(np,nsa)**2)
             do nth = 1, nthmax
                K = (PV-1.d0)*amfp(nsa)*vc**2
                dVI = delthm(nth,np,nr,nsa)*delp(nsa)
                M2(nr,nsa) = M2(nr,nsa) &
                     + K*fI_in(nth,np,nr,nsa) * dVI * JIR(nth,np,nr,nsa)
             end do
          end do
          CALL mtx_allreduce1_real8(sum,3,sum_tot,sum_loc)
          M2(nr,nsa) = sum_tot*rnfp0(nsa)*1.d20*1.d-6 ! [*1.e-20 MJ]

       end do
    end do
    CALL mtx_reset_communicator

  end subroutine

  subroutine particle_flux(Sr)
    use fpcomm
    use fowcomm
    USE fowlib
    USE libmpi
    implicit none
    real(rkind),intent(out) :: &
         Sr(nrstartw:nrendwm,nsastart:nsaend)
    real(rkind), &
         dimension(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend) &
         :: dfdp, dfdth, dfdrm
    real(rkind) :: sum, dVI, gammaI, Drx(3), Fr
    integer :: nth, np, nr, nsa
    REAL(rkind):: sum_tot
    INTEGER:: sum_loc

    CALL mtx_set_communicator(comm_np)

    do nsa = nsastart, nsaend
       do nr = nrstartw, nrendwm
          do nth = 1, nthmax
             call first_order_derivative( &
                  dfdp(nth,:,nr,nsa), fnsp(nth,:,nr,nsa), pm(:,nsa))
          end do

          do np = npstartw, npendwm
             call first_order_derivative( &
                  dfdth(:,np,nr,nsa), fnsp(:,np,nr,nsa), thetam(:,np,nr,nsa))
          end do
       end do

       do np = npstartw, npendwm
          do nth = 1, nthmax
             call first_order_derivative( &
                  dfdrm(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
          end do
       end do

       do nr = nrstartw, nrendwm
          sum = 0.d0
          do np = npstartw, npendwm
             do nth = 1, nthmax
                dVI = JIR(nth,np,nr,nsa) * delp(nsa) * delthm(nth,np,nr,nsa)
                Drx(1) = ( Drpfow(nth,np,nr+1,nsa) &
                         + Drpfow(nth,np,nr,nsa) )*0.5d0
                Drx(2) = ( Drtfow(nth,np,nr+1,nsa) &
                         + Drtfow(nth,np,nr,nsa) )*0.5d0
                Drx(3) = ( Drrfow(nth,np,nr+1,nsa) &
                         + Drrfow(nth,np,nr,nsa) )*0.5d0
                Fr     = ( Frrfow(nth,np,nr+1,nsa) &
                         + Frrfow(nth,np,nr,nsa) )*0.5d0

                gammaI = -1.d0*( &
                     Drx(1) * dfdp(nth,np,nr,nsa) &
                     + Drx(2) * dfdth(nth,np,nr,nsa)/pm(np,nsa) &
                     + Drx(3) * dfdrm(nth,np,nr,nsa) &
                     - Fr     * fnsp(nth,np,nr,nsa) &
                     )/JI(nth,np,nr,nsa)

                sum = sum + gammaI*dVI
             end do
          end do
          CALL mtx_allreduce1_real8(sum,3,sum_tot,sum_loc)
          Sr(nr,nsa) = rnfp0(nsa)*sum_tot
       end do
    end do
    CALL mtx_reset_communicator

  end subroutine particle_flux

  subroutine particle_flux_element(Sr, Sr_Dp, Sr_Dth, Sr_Dr, Sr_F)
    use fpcomm
    use fowcomm
    USE fowlib
    USE libmpi
    implicit none
    real(rkind),dimension(nrstartw:nrendwm,nsastart:nsaend),intent(out) :: &
         Sr, Sr_Dp, Sr_Dth, Sr_Dr, Sr_F
    real(rkind), &
         dimension(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend)::&
         dfdp, dfdth, dfdrm
    real(rkind) :: dVI, Drx(3), Fr, sum_tot
    integer :: nth, np, nr, nsa, sum_loc

    CALL mtx_set_communicator(comm_np)

    do nsa = nsastart,nsaend
       do nr = nrstartw, nrendwm
          do nth = 1, nthmax
             call first_order_derivative( &
                  dfdp(nth,npstartw:npendwm,nr,nsa), &
                  fnsp(nth,npstartw:npendwm,nr,nsa), &
                  pm(npstartw:npendwm,nsa))
          end do
       END do

       do nr = nrstartw, nrendwm
          do np = npstartw, npendwm
             call first_order_derivative( &
                  dfdth(:,np,nr,nsa), fnsp(:,np,nr,nsa), thetam(:,np,nr,nsa))
          end do
       END do

       do np = npstartw, npendwm
          do nth = 1, nthmax
             call first_order_derivative( &
                  dfdrm(nth,np,nrstartw:nrendwm,nsa), &
                  fnsp(nth,np,nrstartw:nrendwm,nsa), &
                  rm(nrstartw:nrendwm))
          end do
       end do

       do nr = nrstartw, nrendwm
          Sr_Dp (nr,nsa) = 0.d0
          Sr_Dth(nr,nsa) = 0.d0
          Sr_Dr (nr,nsa) = 0.d0
          Sr_F  (nr,nsa) = 0.d0
          do np = npstartw, npendwm
             do nth = 1, nthmax
                dVI = JIR(nth,np,nr,nsa) * delp(nsa) * delthm(nth,np,nr,nsa)
                Drx(1) = ( Drpfow(nth,np,nr+1,nsa) &
                         + Drpfow(nth,np,nr,nsa) )*0.5d0
                Drx(2) = ( Drtfow(nth,np,nr+1,nsa) &
                         + Drtfow(nth,np,nr,nsa) )*0.5d0
                Drx(3) = ( Drrfow(nth,np,nr+1,nsa) &
                         + Drrfow(nth,np,nr,nsa) )*0.5d0
                Fr     = ( Frrfow(nth,np,nr+1,nsa) &
                         + Frrfow(nth,np,nr,nsa) )*0.5d0

                Sr_Dp (nr,nsa) = Sr_Dp (nr,nsa) &
                     - rnfp0(nsa)*Drx(1)*dfdp (nth,np,nr,nsa) &
                       /JI(nth,np,nr,nsa)*dVI
                Sr_Dth(nr,nsa) = Sr_Dth(nr,nsa) &
                     - rnfp0(nsa)*Drx(2)*dfdth(nth,np,nr,nsa) &
                       /JI(nth,np,nr,nsa)*dVI
                Sr_Dr (nr,nsa) = Sr_Dr (nr,nsa) &
                     - rnfp0(nsa)*Drx(3)*dfdrm(nth,np,nr,nsa) &
                       /JI(nth,np,nr,nsa)*dVI
                Sr_F  (nr,nsa) = Sr_F  (nr,nsa) &
                     + rnfp0(nsa)*Fr    *fnsp (nth,np,nr,nsa) &
                       /JI(nth,np,nr,nsa)*dVI
             end do
          end do
          CALL mtx_allreduce1_real8(Sr_Dp (nr,nsa),3,sum_tot,sum_loc)
          Sr_Dp (nr,nsa)=sum_tot
          CALL mtx_allreduce1_real8(Sr_Dth(nr,nsa),3,sum_tot,sum_loc)
          Sr_Dth(nr,nsa)=sum_tot
          CALL mtx_allreduce1_real8(Sr_Dr (nr,nsa),3,sum_tot,sum_loc)
          Sr_Dr (nr,nsa)=sum_tot
          CALL mtx_allreduce1_real8(Sr_F  (nr,nsa),3,sum_tot,sum_loc)
          Sr_F  (nr,nsa)=sum_tot
          Sr(nr,nsa) = Sr_Dp(nr,nsa)+Sr_Dth(nr,nsa)+Sr_Dr(nr,nsa)+Sr_F(nr,nsa)
       end do
    end do

    CALL mtx_reset_communicator

  end subroutine particle_flux_element

  subroutine effective_diffusion_coefficient(Deff)
    use fpcomm
    use fowcomm
    USE fowlib
    implicit none
    real(rkind),dimension(nrmax,nsamax),intent(out) :: Deff
    real(rkind),dimension(nrmax,nsamax) :: Sr, dndr
    integer :: nr, nsa

    call particle_flux(Sr)
    do nsa = nsastart, nsaend
       call first_order_derivative(dndr(nrstartw:nrendwm,nsa), &
                                   rnsl(nrstartw:nrendwm,nsa), &
                                   rm(nrstartw:nrendwm))
    end do

    do nsa = nsastart, nsaend
       do nr = nrstartw, nrendwm
          Deff(nr,nsa) = -1.d0*Sr(nr,nsa)/dndr(nr,nsa)
       end do
    end do

  end subroutine effective_diffusion_coefficient

  subroutine total_N(N,fI,nsa)
    use fpcomm
    use fowcomm
    USE libmpi
    implicit none
    real(rkind),intent(in) :: fi(:,:,:,:)
    integer,intent(in) :: nsa
    real(rkind),intent(out) :: N
    integer :: nth,np,nr
    real(rkind) :: sumI, dVI, sum_tot
    INTEGER:: sum_loc

    
    sumI = 0.d0
    do nr = nrstartw, nrendwm
       do np = npstartw, npendwm
          do nth = 1, nthmax
             dVI = delp(nsa)*delthm(nth,np,nr,nsa)*delr
             sumI = sumI+FI(nth,np,nr,nsa)*dVI*JI(nth,np,nr,nsa)
          end do
       end do
    end do

    CALL mtx_set_communicator(comm_nrnp)
    CALL mtx_allreduce1_real8(sumI,3,sum_tot,sum_loc)
    N=sum_tot
    
  end subroutine total_N

end module fowdistribution
