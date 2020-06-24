module fpcalr

! ****************************************
!     Radial transport
! ****************************************

contains

  subroutine fp_calr

      use fpcomm
      use plprof
      use cdbmfp_mod
      use libmpi,only: mtx_set_communicator, mtx_reset_communicator
      use fpmpi,only: p_theta_integration
      implicit none
      integer:: nsa, nsba, ns, nr, nth, np, ng, nsb, nr1, nr2
      real(kind8):: rhon, rtfpl, factr, factp, sv, ppp, ppm
      real(kind8):: psib, pcos, x, etal, sumd, sumf, delh
      real(kind8):: dndr, nedge, fact, dint_d, dint_f, wrl
      real(kind8):: sum1, temp1, srhodp, srhodm, srhofp, srhofm
      real(kind8):: wrh, dfdr_d, dfdr_f, f_r2, dfdr_r2, f_r1, dfdr_r1
      real(kind8):: shear,pnel,rhoni,dpdr,dvexbdr,calf,ckap,cexb,fsz,fez
      real(kind8),dimension(1:nrmax+1):: chi_cdbm
      type(pl_plf_type),dimension(nsmax):: plf
      double precision:: densm, densp, rgama
      integer:: isw_d

!---- calculation of cdbm diffusion coefficient ----

      if(modeld_rdep.eq.2) then
         if(nrank.eq.0) then
            nr1=1
            nr2=nrmax+1
         else
            nr1=nrstart
            nr2=nrendwg
         end if
         do nr=nr1,nr2
            rhon=rg(nr)

            if(nr.eq.1) then        ! magnetic shear s=(r/q)(dq/dr)
               shear=0.d0
            else if(nr.eq.nrmax+1) then  !qlm(nrmax+1)=qlg(rhon=1)
               shear=(rhon/qlm(nr))*(qlm(nr)-qlm(nr-1)) &
                     *0.5d0/(rg(nr)-rm(nr-1))
            else
               shear=(rhon/(rm(nr)-rm(nr-1)))*(qlm(nr)-qlm(nr-1)) &
                     *2.d0/(qlm(nr)+qlm(nr-1))
            end if

            pnel=0.d0               ! electron density
            do ns=1,nsmax
               if(id_ns(ns).eq.-1) then
                  if(nr.eq.1) then
                     pnel=pnel+rn_temp(nr,ns)*1.d20
                  else if(nr.eq.nrmax+1) then
                     pnel=pnel+rn_temp(nr-1,ns)*1.d20
                  else
                     pnel=pnel+0.5d0*(rn_temp(nr-1,ns)+rn_temp(nr,ns))*1.d20
                  end if
               end if
            end do

            rhoni=0.d0              ! ion mass density
            do ns=1,nsmax
               if(id_ns(ns).eq.1) then
                  if(nr.eq.1) then
                     rhoni=rhoni+pa(ns)*amp*rn_temp(nr,ns)*1.d20
                  else if(nr.eq.nrmax+1) then
                     rhoni=rhoni+pa(ns)*amp*rn_temp(nr-1,ns)*1.d20
                  else
                     rhoni=rhoni+pa(ns)*amp &
                                *0.5d0*(rn_temp(nr-1,ns)+rn_temp(nr,ns))*1.d20
                  end if
               end if
            end do

            dpdr=0.d0               ! pressure gradient
            ppp=0.d0
            ppm=0.d0
            do ns=1,nsmax
               if(id_ns(ns).eq.1.or.id_ns(ns).eq.-1) then
                  if(nr.eq.1) then
                     ppp=ppp+rhoni+rn_temp(nr,ns)*1.d20*rt_temp(nr,ns)*rkev
                     ppm=ppm+rhoni+rn_temp(nr,ns)*1.d20*rt_temp(nr,ns)*rkev
                  else if(nr.eq.nrmax) then
                     ppp=ppp+rhoni+rn_temp(nr  ,ns)*1.d20*rt_temp(nr  ,ns)*rkev
                     ppm=ppm+rhoni+rn_temp(nr-1,ns)*1.d20*rt_temp(nr-1,ns)*rkev
                  else if(nr.eq.nrmax+1) then
                     ppp=ppp+rhoni+rn_temp(nr-1,ns)*1.d20*rt_temp(nr-1,ns)*rkev
                     ppm=ppm+rhoni+rn_temp(nr-2,ns)*1.d20*rt_temp(nr-2,ns)*rkev
                  else
                     ppp=ppp+rhoni+rn_temp(nr+1,ns)*1.d20*rt_temp(nr+1,ns)*rkev
                     ppm=ppm+rhoni+rn_temp(nr-1,ns)*1.d20*rt_temp(nr-1,ns)*rkev
                  end if
               end if
            end do
            if(nr.eq.1) then
               dpdr=0.d0
            else if(nr.eq.nrmax) then
               dpdr=(ppp-ppm)/(ra*(rm(nr  )-rm(nr-1)))
            else if(nr.eq.nrmax+1) then
               dpdr=(ppp-ppm)/(ra*(rm(nr-1)-rm(nr-2)))
            else
               dpdr=(ppp-ppm)/(ra*(rm(nr+1)-rm(nr-1)))
            end if

            dvexbdr=0.d0
            calf=1.d0
            ckap=1.d0
            cexb=0.d0
!            fsz=1.d0     ! option
!            curvz=0.d0   ! option
!            fez=0.d0     ! option

            call cdbmfp(bb,rr,ra*rhon,rkap,qlm(nr),shear,pnel,rhoni,dpdr, &
                      dvexbdr,calf,ckap,cexb,modeld_cdbm,chi_cdbm(nr))
         end do
      end if

      do nsa=nsastart,nsaend
         ns=ns_nsa(nsa)
         nsba=nsb_nsa(nsa)

      do nr=nrstart,nrendwg
         rhon=rg(nr)
         if(nr.eq.1) then
            nr1=nr
            nr2=nr
         else if(nr.eq.nrmax+1) then
            nr1=nr-1
            nr2=nr-1
         else
            nr1=nr-1
            nr2=nr
         endif

         if(nr.eq.1) then
            rtfpl= rtfp_g(nr,nsa)/rtfp0(nsa)
         else
            rtfpl=rtfps(nsa)/rtfp0(nsa)
         end if

!------------- set r dependence
         select case(modeld_rdep)
         case(0)
            factr=(drr0-drrs)*(1.d0-rhon**2)+drrs
         case(1)
            factr=pi*rr*qlm(nr)*deltab_b**2 *ptfp0(ns)/amfp(ns)
         case(2)
            factr=factor_cdbm*chi_cdbm(nr)
         case default
            write(6,*) 'xx fpcalr: undefined modeld_rdep'
            stop
         end select

!------------- set edge value
         select case(modeld_edge)
         case(1)
            if(rhon.gt.rho_edge) factr=drr_edge
         case(2)
            if(rhon.gt.rho_edge) factr=factor_drr_edge*factr
         end select

         do np=npstart,npend
            do nth=1,nthmax

!------------- set p dependence
               select case(modeld_pdep)
               case(0) ! no p dependence
                  factp=1.d0
               case(1) ! proportional to 1/sqrt{p_perp}
                  factp=(rtfpl/(rtfpl+pm(np,ns)**2*sinm(nth)**2))**0.25d0
               case(2) ! proportional to 1/p_perp
                  factp=sqrt(rtfpl/(rtfpl+pm(np,ns)**2*sinm(nth)**2))
               case(3) ! proportional to 1/p_perp^2
                  factp=rtfpl/(rtfpl+pm(np,ns)**2*sinm(nth)**2)
               case(4) ! stochastic delta b /b; relativistic
                  rgama=sqrt(1.d0+theta0(ns)*pm(np,ns)**2)
                  factp=pm(np,ns)*abs(cosm(nth))/rgama
               case(5) ! proportional to  p^pdep_exp
                  factp=sqrt(pm(np,ns)**2)**pdep_exp &
                        /sqrt(1.d0+pm(np,ns)**2*sinm(nth)**2)
               ! case(6) ! proportional to p_perp^pdep_exp
               !    factp=sqrt(pm(np,ns)**2*sinm(nth)**2)**pdep_exp &
               !          /sqrt(1.d0+pm(np,ns)**2*sinm(nth)**2)
               ! case(7) ! proportional to p^pdep_exp
               !    factp=sqrt(pm(np,ns)**2)**pdep_exp &
               !          /sqrt(1.d0+pm(np,ns)**2)
               end select
               drr(nth,np,nr,nsa)= factr*factp/ra**2*rlamda_rg(nth,nr)   ! normalization for rhon
            enddo
         enddo

! ------ set pinch term

         select case(modeld_pinch)
         case(0) ! no pinch
            fact=0.d0
         case(1) ! no radial particle transport (particle flux = 0)
            dint_d=0.d0
            dint_f=0.d0
            do np=npstart,npend
               do nth=1,nthmax
                  if(timefp.eq.0.and.n_impl.eq.0)then
                     if(nr.eq.1)then
                        wrl=0.25d0 ! not necessary
                     else
                        wrl=(4.d0*rg(nr)+delr)/(8.d0*rg(nr))
                     end if
                  else
                     wrl=weighr(nth,np,nr,nsa)
                  end if

                  if(nr.eq.1)then
                     dfdr_r1 = ( fnsp(nth,np,nr,nsa)-fs0(nth,np,nsa) ) &
                               / delr *2.d0*0
                     f_r1 = fs0(nth,np,nsa)
                     srhodm=dfdr_r1 * drr(nth,np,nr,nsa)
                     srhofm=f_r1    * drr(nth,np,nr,nsa)
                  elseif(nr.eq.nrmax+1)then ! fs2 = f at rho=1+delr/2
                     dfdr_r1 = ( fs2(nth,np,nsa)-fnsp(nth,np,nr-1,nsa) ) / delr
                     if(modeld_boundary.eq.0)then
                        f_r1 = ( (1.d0-wrl)*fs2(nth,np,nsa) &
                             + wrl*fnsp(nth,np,nr-1,nsa) )
                     elseif(modeld_boundary.eq.1)then
                        f_r1 = fs1(nth,np,nsa)
                     end if
                     srhodm=dfdr_r1 * drr(nth,np,nr,nsa)
                     srhofm=f_r1    * drr(nth,np,nr,nsa)
                  else
                     dfdr_r1 = ( fnsp(nth,np,nr,nsa)-fnsp(nth,np,nr-1,nsa) ) &
                               / delr
                     f_r1 = ( (1.d0-wrl)*fnsp(nth,np,nr,nsa) &
                            + wrl*fnsp(nth,np,nr-1,nsa) )
                     srhodm=dfdr_r1 * drr(nth,np,nr,nsa)
                     srhofm=f_r1    * drr(nth,np,nr,nsa)
                  end if
                  dint_d = dint_d + volp(nth,np,ns)*srhodm
                  dint_f = dint_f + volp(nth,np,ns)*srhofm
               end do
            end do
! integration
            call mtx_set_communicator(comm_np)
            call p_theta_integration(dint_d)
            call p_theta_integration(dint_f)
            call mtx_reset_communicator

!         write(*,'(a,2i3,2e14.6)') "dint=", nsa,nr,dint_d, dint_f

            fact=dint_d/dint_f

         case(2) ! pinch for gaussian profile
            fact=-2.d0*factor_pinch*rhon/ra
         end select

         do np=npstart,npend
            do nth=1,nthmax
               frr(nth,np,nr,nsa) = fact * drr(nth,np,nr,nsa)
            end do
         end do

      end do ! nr
      end do ! nsa

      return
      end subroutine fp_calr

    end module fpcalr
