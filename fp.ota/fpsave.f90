!     $Id: fpsave.f90,v 1.41 2013/02/08 07:36:24 nuga Exp $
!
! *************************
!     SAVE DATA ROUTINE
! *************************
!
      module fpsave

      use fpcomm
      use plprof
      use fpexec
      use libbes,only: BESEKNX,besjnv

      contains

!----------------------------------

      subroutine fpssub
!
      use libmtx
      use fpmpi
      implicit none
      integer:: nr, nsa, ns,j=0

      if(isave.ne.0) return

      do nsa=nsastart,nsaend
         rndrs(nsa) =0.d0
         rpdrs(nsa) =0.d0
         do nr=nrstart, nrend
            rpdrl(nr,nsa)=0.d0
            rndrl(nr,nsa)=0.d0
         end do
      end do

      call mtx_reset_communicator

      call define_bulk_np
      call moment_0th_order(fnsp,rnsl)
      call moment_0th_order(fnsm,rnsl_prev)
      call moment_1st_order(fnsp,rjsl)
      call moment_2nd_order(fnsp,rwsl)
      call moment_2nd_order(fnsm,rwsl_prev)
      call power_from_diffusion_coef
      call power_from_source_term
      call mtx_reset_communicator
      do nr=nrstart,nrend
         do nsa=nsastart,nsaend
            call bulk_temperature(np_bulk(nr,nsa),nr,nsa)
         enddo ! nsa
      enddo ! nr
      call power_from_radial_transport
      call total_particle_source
      call mtx_reset_communicator
      call fpsavecomm
      isave=1
      return
      end subroutine fpssub
!
! *************************
!     save global data
! *************************
!
      subroutine fpsglb
!
        use fpsub
      implicit none
      integer:: nsa, nsb, nr
      real(8):: rtemp, rtemp2

      ntg1=ntg1+1
      call fp_adjust_ntg1

      ptg(ntg1)=timefp

      do nsa=1,nsamax
         pnt(nsa,ntg1)=0.d0
         pit(nsa,ntg1)=0.d0
         pwt(nsa,ntg1)=0.d0
         ppct(nsa,ntg1)=0.d0
         ppwt(nsa,ntg1)=0.d0
         ppet(nsa,ntg1)=0.d0
         plht(nsa,ntg1)=0.d0
         pfwt(nsa,ntg1)=0.d0
         pect(nsa,ntg1)=0.d0
         pwrt(nsa,ntg1)=0.d0
         pwmt(nsa,ntg1)=0.d0
         pwt2(nsa,ntg1)=0.d0
         pspbt(nsa,ntg1)=0.d0
         pspft(nsa,ntg1)=0.d0
         pspst(nsa,ntg1)=0.d0
         psplt(nsa,ntg1)=0.d0
         pwtd(nsa,ntg1)=0.d0
         do nsb=1,nsbmax
            ppct2(nsb,nsa,ntg1)= 0.d0
         end do
         pdr(nsa,ntg1)=0.d0
         pndr(nsa,ntg1)=0.d0
         ptt_bulk(nsa,ntg1)=0.d0
         ppst(nsa,ntg1)=0.d0
         pplt(nsa,ntg1)=0.d0
      enddo

      do nsa=1,nsamax
         do nr=1,nrmax
            pnt(nsa,ntg1) =pnt(nsa,ntg1) +rns(nr,nsa)*volr(nr)
            pit(nsa,ntg1) =pit(nsa,ntg1) +rjs(nr,nsa)*volr(nr)
            pwt(nsa,ntg1) =pwt(nsa,ntg1) +rws(nr,nsa)*volr(nr)
            pwtd(nsa,ntg1)=pwtd(nsa,ntg1)+rws123(nr,nsa)*volr(nr)
            ppct(nsa,ntg1)=ppct(nsa,ntg1)+rpcs(nr,nsa)*volr(nr)
            ppwt(nsa,ntg1)=ppwt(nsa,ntg1)+rpws(nr,nsa)*volr(nr)
            ppet(nsa,ntg1)=ppet(nsa,ntg1)+rpes(nr,nsa)*volr(nr)
            plht(nsa,ntg1)=plht(nsa,ntg1)+rlhs(nr,nsa)*volr(nr)
            pfwt(nsa,ntg1)=pfwt(nsa,ntg1)+rfws(nr,nsa)*volr(nr)
            pect(nsa,ntg1)=pect(nsa,ntg1)+recs(nr,nsa)*volr(nr)
            pwrt(nsa,ntg1)=pwrt(nsa,ntg1)+rwrs(nr,nsa)*volr(nr)
            pwmt(nsa,ntg1)=pwmt(nsa,ntg1)+rwms(nr,nsa)*volr(nr)
            pspbt(nsa,ntg1)=pspbt(nsa,ntg1)+rspb(nr,nsa)*volr(nr)
            pspft(nsa,ntg1)=pspft(nsa,ntg1)+rspf(nr,nsa)*volr(nr)
            pspst(nsa,ntg1)=pspst(nsa,ntg1)+rsps(nr,nsa)*volr(nr)
            psplt(nsa,ntg1)=psplt(nsa,ntg1)+rspl(nr,nsa)*volr(nr)
            pdr(nsa,ntg1) = pdr(nsa,ntg1) + rpdr(nr,nsa)*volr(nr)
            pndr(nsa,ntg1) = pndr(nsa,ntg1) + rndr(nr,nsa)*volr(nr)
            ptt_bulk(nsa,ntg1)=ptt_bulk(nsa,ntg1) &
                              +rt_bulk(nr,nsa)*volr(nr)*rns(nr,nsa)
            ppst(nsa,ntg1)=ppst(nsa,ntg1)+rpss(nr,nsa)*volr(nr)
            pplt(nsa,ntg1)=pplt(nsa,ntg1)+rpls(nr,nsa)*volr(nr)

            if(modelr.eq.1) then
               call fpnewton(nr,nsa,rns(nr,nsa),rws(nr,nsa),rtemp)
            else
               rtemp=0.d0
               rtemp2=0.d0
            endif
            pwt2(nsa,ntg1) =pwt2(nsa,ntg1) +rtemp*volr(nr)/1.d6 &
                                  *(1.5d0*rns(nr,nsa)*1.d20*aee*1.d3)
            do nsb=1,nsbmax
               ppct2(nsb,nsa,ntg1)=ppct2(nsb,nsa,ntg1)       &
                                  +rpcs2(nr,nsb,nsa)*volr(nr)
            end do
         enddo

         pit(nsa,ntg1) =pit(nsa,ntg1)/(2.d0*pi*rr)
         pspt(nsa,ntg1)=pspbt(nsa,ntg1)+pspft(nsa,ntg1) &
                       +pspst(nsa,ntg1)+psplt(nsa,ntg1)

         if(pnt(nsa,ntg1).ne.0.d0) then
            ptt(nsa,ntg1) =pwt(nsa,ntg1)*1.d6   &
                 /(1.5d0*pnt(nsa,ntg1)*1.d20*aee*1.d3)
            ptt2(nsa,ntg1) =pwt2(nsa,ntg1)*1.d6 &
                 /(1.5d0*pnt(nsa,ntg1)*1.d20*aee*1.d3)
            ptt_bulk(nsa,ntg1) = ptt_bulk(nsa,ntg1)/pnt(nsa,ntg1)
         else
            ptt(nsa,ntg1)=0.d0
            ptt2(nsa,ntg1)=0.d0
            ptt_bulk(nsa,ntg1) = 0.d0
         endif
         pnt(nsa,ntg1) =pnt(nsa,ntg1)/tvolr
      enddo

      return
      end subroutine fpsglb
!
! *************************
!     save profile data
! *************************
!
      subroutine fpsprf
!
        use fpsub
      implicit none
      integer:: nr, nsa, nsb
      real(8):: rs, rtemp

      ntg2=ntg2+1
      call fp_adjust_ntg2

      rtg(ntg2)=timefp

      do nsa=1,nsamax
         do nr=1,nrmax
            rnt(nr,nsa,ntg2) = rns(nr,nsa)
            if(model_disrupt.eq.1) rate_runaway(nr,ntg2)=rfp(nr)
            rjt(nr,nsa,ntg2) = rjs(nr,nsa)
            rwt(nr,nsa,ntg2) = rws(nr,nsa)
            rpct(nr,nsa,ntg2)= rpcs(nr,nsa)
            rpwt(nr,nsa,ntg2)= rpws(nr,nsa)
            rpet(nr,nsa,ntg2)= rpes(nr,nsa)
            rlht(nr,nsa,ntg2)= rlhs(nr,nsa)
            rfwt(nr,nsa,ntg2)= rfws(nr,nsa)
            rect(nr,nsa,ntg2)= recs(nr,nsa)
            rwrt(nr,nsa,ntg2)= rwrs(nr,nsa)
            rwmt(nr,nsa,ntg2)= rwms(nr,nsa)
            do nsb=1,nsbmax
               rpct2(nr,nsb,nsa,ntg2)= rpcs2(nr,nsb,nsa)
            end do
            rspbt(nr,nsa,ntg2)= rspb(nr,nsa)
            rspft(nr,nsa,ntg2)= rspf(nr,nsa)
            rspst(nr,nsa,ntg2)= rsps(nr,nsa)
            rsplt(nr,nsa,ntg2)= rspl(nr,nsa)
            rpdrt(nr,nsa,ntg2)= rpdr(nr,nsa)
            rndrt(nr,nsa,ntg2)= rndr(nr,nsa)
            rtt_bulk(nr,nsa,ntg2) = rt_bulk(nr,nsa)
!
            if(rns(nr,nsa).ne.0.d0) then
               if(modelr.eq.0)then
                  rtt(nr,nsa,ntg2) = rws(nr,nsa)*1.d6 &
                       /(1.5d0*rns(nr,nsa)*1.d20*aee*1.d3)
               elseif(modelr.eq.1)then
                  call fpnewton(nr,nsa,rns(nr,nsa),rws(nr,nsa),rtemp)
                  rtt(nr,nsa,ntg2) = rtemp
               end if
            else
               rtt(nr,nsa,ntg2) = 0.d0
            endif
            ret(nr,ntg2) = e1(nr)
            rs=rsrhon(rm(nr))
            rqt(nr,ntg2) = rs*bb*2.d0/(rr*(bp(nr)+bp(nr+1)))
!            rt_t(nr,nsa) = rtt(nr,nsa,ntg2)
         enddo
      enddo
      do nr=1,nrmax
         eptr(nr,ntg2)=e1(nr)
      end do

      return
      end subroutine fpsprf
!
! ***********************************************************
!
!                         result
!
! ***********************************************************
!
      subroutine fpwrtglb
!
      implicit none
      integer:: nsa, nsb
      real(8):: rtotalpw, rtotalpc,rtotalsp,rtotalpc2
      real(8):: rtotaldr,rtotallh,rtotalfw,rtotalec,rtotalwr,rtotalwm,rtotalip
      character:: fmt0*50
!
      write(6,*)"--------------------------------------------"
      write(6,*)"-----global data"
      write(6,101) ptg(ntg1)*1000

      do nsa=1,nsamax
         if(modelr.eq.0)then
            write(6,112) nsa,ns_nsa(nsa), &
              pnt(nsa,ntg1),ptt(nsa,ntg1),pwt(nsa,ntg1),pit(nsa,ntg1),rndrs(nsa),ptt_bulk(nsa,ntg1)
         else
            write(6,112) nsa,ns_nsa(nsa), &
              pnt(nsa,ntg1),ptt2(nsa,ntg1),pwt(nsa,ntg1),pit(nsa,ntg1),rndrs(nsa),ptt_bulk(nsa,ntg1)
         end if
      enddo

      rtotalpw=0.d0
      rtotalpc=0.d0
      rtotalsp=0.d0
      rtotalpc2=0.d0
      rtotallh=0.d0
      rtotalfw=0.d0
      rtotalec=0.d0
      rtotalwr=0.d0
      rtotalwm=0.d0
      rtotalip=0.d0

      write(fmt0,'(a48)') &
           '(8x,2i2," pc,pw,pe,pdr,pplt=",6x,1p5e12.4)'
      do nsa=1,nsamax
         write(6,fmt0) nsa,ns_nsa(nsa), &
              ppct(nsa,ntg1),ppwt(nsa,ntg1),ppet(nsa,ntg1),rpdrs(nsa),pplt(nsa,ntg1)
         rtotalpw=rtotalpw + ppwt(nsa,ntg1)
         rtotalpc=rtotalpc + ppct(nsa,ntg1)
         rtotalsp=rtotalsp + pspt(nsa,ntg1)
         rtotalpc2 = rtotalpc2 +ppct(nsa,ntg1)-ppct2(nsa,nsa,ntg1)
         rtotallh=rtotallh + plht(nsa,ntg1)
         rtotalfw=rtotalfw + pfwt(nsa,ntg1)
         rtotalec=rtotalec + pect(nsa,ntg1)
         rtotalwr=rtotalwr + pwrt(nsa,ntg1)
         rtotalwm=rtotalwm + pwmt(nsa,ntg1)
         rtotalip=rtotalip + pit(nsa,ntg1)
      end do
      do nsa=1,nsamax
         if(nsbmax.gt.1) then
            write(6,104) nsa,ns_nsa(nsa), &
                 (ppct2(nsb,nsa,ntg1),nsb=1,nsbmax)
         endif
      end do

      do nsa=1,nsamax
         write(6,108) nsa,ns_nsa(nsa),pspbt(nsa,ntg1),pspft(nsa,ntg1), &
                                      pspst(nsa,ntg1),psplt(nsa,ntg1), &
                                      pspst(nsa,ntg1)+psplt(nsa,ntg1)
      end do
      write(6,105) rtotalpw,rtotalwr,rtotalwm
      write(6,115) rtotallh,rtotalfw,rtotalec
      write(6,107) rtotalpc
      write(6,109) rtotalsp
      write(6,110) rtotalpc2
      write(6,'("total plasma current   [ma]",1pe12.4)') rtotalip

      return

  101 format(' time=',f12.3,' ms')
  112 format(' nsa,ns=',2i2,' n,t,w,i,dn,t2=',1pe11.4,1p6e12.4)
  104 format('        ',2i2,' pcab    =',10x,1p14e12.4)
  105 format('total absorption power [mw]', 1pe12.4,'    wr:',1pe12.4,'    wm:',1pe12.4)
  115 format('   absorption power [mw] lh', 1pe12.4,'    fw:',1pe12.4,'    ec:',1pe12.4)
 106  format(f12.4, 8e12.4)
 107  format('total collision power  [mw]', 1pe12.4)
 108  format('        ',2i2,' pspb/f/s/l/s+l=',4x,1p5e12.4)
 109  format('total source power     [mw]', 1pe12.4)
 110  format('collision balance      [mw]', 1pe12.4)

      end subroutine fpwrtglb

! ***********************************************************

      subroutine fpwrtprf
!
        use fpsub
      implicit none
      integer:: nsa, nr, ns
      real(8):: rtemp
      character:: fmt0*50
!
!      write(fmt0,'(a15)') '(2i3,1p20e13.4)'
!      write(fmt0,'(a44)') '(2i3,1p8e16.8,1p5e12.3e3,1pe12.4,1p5e12.3e3)'

      write(6,*)"-----radial profile data"
      write(6,'(a,f12.3)') " time=", timefp*1000

      if(model_disrupt.eq.0)then
         write(fmt0,'(a44)') '(2i3,1p8e12.4,1p5e12.3e3,1pe12.4,1p5e12.3e3)'
         write(6,106)
      else
         write(fmt0,'(a44)') '(2i3,1p8e12.4,1p5e12.3e3,1pe12.4,1p5e12.3e3)'
         write(6,108)
      end if

      do nsa=1,nsamax
         ns=ns_nsa(nsa)
         do nr=1,nrmax
            if(model_disrupt.ne.0)then
               write(6,fmt0) nsa,ns_nsa(nsa),&
                    rm(nr),rnt(nr,nsa,ntg2),rtt(nr,nsa,ntg2), &
                    post_tau_ta(nr,nsa),     &
                    rsps(nr,nsa), &
                    e1(nr), conduct_sp(nr), rn_disrupt(nr),   &
                    rn_runaway(nr), rn_runaway(nr)-rn_drei(nr), rj_ohm(nr), rj_runaway(nr), &
                    rj_bs(nr), &
                    rt_quench(nr), &
                    rfp(nr), rconnor(nr), rfp_ava(nr), &
                    er_crit(nr)
            else
               write(6,fmt0) nsa,ns_nsa(nsa),&
                    rm(nr),rnt(nr,nsa,ntg2),rtt(nr,nsa,ntg2), &
                    rjt(nr,nsa,ntg2),rpct(nr,nsa,ntg2),       &
                    rpet(nr,nsa,ntg2), &
!                    rpwt(nr,nsa,ntg2),&
                    rns_delf(nr,ns), &
                    rspbt(nr,nsa,ntg2), &
                    rsps_cx(nr,nsa), &
                    rt_bulk(nr,nsa)!, &
            end if
         enddo
      enddo
      return
  106 format( &
           'nsa/ns',5x,'rm',10x,' n',8x,' t    ',6x, &
           ' j    ',5x,'pc     ',5x,'pe     ',5x,   &
           'n_b',9x,'pnbi',8x,'pdrp',8x,'tbulk' )
  107 format( &
           'nsa/ns',5x,'rm',10x,' n',8x,' t    ',6x, &
           ' j     ',5x,'pc     ',5x,'pc12     ',5x,  &
           'pc11   ',5x,'pe     ',5x,'e_ind    ',5x,  &
           'pe_ind ',4x,'psip   ',5x,'dpsip ',5x,'sigma' )
  108 format( &
           'nsa/ns',5x,'rm',10x,' n',8x,' t    ',6x, &
           ' nbi   ',5x,'j_fp   ',5x,'e1    ',5x,  &
           'sigma  ',5x,'n_bulk ',5x,'n_run  ',5x,'n_second',4x,  &
           'j_ohm  ',5x,'j_run  ',5x,'j_bs   ',5x,'t      ',5x, &
           'r_rate ',5x,'rconnor',5x,'rate_a ',5x,'e_c')

      end subroutine fpwrtprf
! ***********************************************************
!
      subroutine fpwrtsnap
!
      implicit none
      integer:: nr, nsa, nsb, ns
      real(8):: rnute, resist
      real(8):: taue_col, sigma_sp, fact
      real(8):: taue_col2, sigma_sp2

!-----check of conductivity--------
      if(ntg1.ne.1.and.nrank.eq.0)then
         nr=1
         nsa=1
         nsb=2
         ns=ns_nsa(nsa)
!         do nsa=1,nsamax
!         do nsb=1,nsbmax
            rnute=rnud(nr,nsa,nsa)*sqrt(2.d0)*rnfp(nr,nsa)/rnfp0(nsa)     &
                 *(ptfp0(nsa)/ptfd(1,nsa))**2
            resist=rnute*amfp(nsa)/rnfp(1,nsa)/aefp(nsa)**2/1.d20
!----------
            fact=aefp(nsa)**2*aefd(nsb)**2*lnlam(nr,nsb,nsa)/(eps0**2)
            taue_col=3.d0*sqrt((2.d0*pi)**3)/fact*sqrt(amfp(1)*(aee*rtt(nr,nsa,ntg2)*1.d3)**3)/rns(nr,nsb)*1.d-20
            taue_col2=3.d0*sqrt((2.d0*pi)**3)/fact*sqrt(amfp(1)*(aee*ptt_bulk(nsa,ntg1)*1.d3)**3)/rns(nr,nsb)*1.d-20
            sigma_sp=1.96d0*rns(nr,nsa)*1.d20*aefp(nsa)**2*taue_col/amfp(nsa) ! p. 174
            sigma_sp2=1.96d0*rns(nr,nsa)*1.d20*aefp(nsa)**2*taue_col2/amfp(nsa) ! p. 174
!            if(modela.eq.0)then
!               sigma_sp=sigma_sp*zeff*(1.d0+0.27d0*zeff-0.27d0)/(1.d0+0.47d0*zeff-0.47d0)
!            end if

            if(e0.ne.0)then
               write(6,'(a,1pe16.8)') &
                    " theta0= ", theta0(ns)
               write(6,'(a,1pe16.8,a,1pe16.8,a,1pe16.8,a,1pe16.8)') &
                    " whole space sigma_norm= ",(rjs(1,1)+rjs(1,2))/e1(1)*1.d6*resist, &
                    " sigma=j/e*1.d6=", (rjs(1,nsa)+rjs(1,2))/e1(1)*1.d6
               write(6,'(a,1pe16.8,a,1pe16.8)') &
                    " sigma_spitzer= ",sigma_sp, &
                    " sigma_spitzer_bulk=", sigma_sp2
               write(6,'(a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
                    " j_norm=j/q_e/n(t)/v_ta0=", rjs(1,nsa)/(abs(aefp(1))*rns(1,nsa)*1.d20*vtfp0(1))*1.d6, &
                    " j_norm=sigma_norm*e0=",(rjs(1,1))/e1(1)*1.d6*resist*e0, &
                    " j_norm=(j_e+j_i)=", (rjs(1,1)+rjs(1,2))/(abs(aefp(1))*rns(1,nsa)*1.d20*vtfp0(1))*1.d6
               write(6,'(a,1pe14.6,a,1pe14.6)') &
                    " sigma_bulk_karney =    ", &
                    ( rjs(1,nsa)/(abs(aefp(1))*rns(1,nsa)*1.d20*vtfp0(1))*1.d6 &
                    -0.5d0*(rate_runaway(1,ntg2)/e0 )*pg(npmax+1,1)**2 )/ &
                    e0/resist, &
                    " j_bulk_karney=", &
                    rjs(1,nsa)/(abs(aefp(1))*rns(1,nsa)*1.d20*vtfp0(1))*1.d6 &
                    -0.5d0*(rate_runaway(1,ntg2)/e0 )*pg(npmax+1,1)**2

            end if
!         end do
!         end do
      end if
!----end of conductivity check---------

      return

      end subroutine fpwrtsnap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fpsavecomm

      use libmtx
      use fpmpi
      implicit none
      integer:: nr, nsa, nsb, np, nsw, n,i
      double precision,dimension(nrstart:nrend,nsamax):: work
      double precision,dimension(nrmax,nsamax):: workg
      double precision,dimension(nsamax):: temp_nsanr
      double precision,dimension(npmax,nrstart:nrend):: temp_npnr1
      double precision,dimension(npmax,nrmax):: temp_npnr2
      double precision,dimension(npmax,nrmax,nsastart:nsaend):: temp_npnr3
      integer,dimension(nsamax):: vloc

      call mtx_set_communicator(comm_nsanr)
      nsw=nsaend-nsastart+1
      do n=1,nsw
         nsa=n+nsastart-1
         call fp_gatherv_real8_sav(rnsl,savlen(nrank+1),rns,n,nsa)
         call fp_gatherv_real8_sav(rnsl_prev,savlen(nrank+1),rns_prev,n,nsa)
         call fp_gatherv_real8_sav(rjsl,savlen(nrank+1),rjs,n,nsa)
         call fp_gatherv_real8_sav(rwsl,savlen(nrank+1),rws,n,nsa)
         call fp_gatherv_real8_sav(rwsl_prev,savlen(nrank+1),rws_prev,n,nsa)
         call fp_gatherv_real8_sav(rws123l,savlen(nrank+1),rws123,n,nsa)
         call fp_gatherv_real8_sav(rpcsl,savlen(nrank+1),rpcs,n,nsa)
         call fp_gatherv_real8_sav(rpwsl,savlen(nrank+1),rpws,n,nsa)
         call fp_gatherv_real8_sav(rpesl,savlen(nrank+1),rpes,n,nsa)
         call fp_gatherv_real8_sav(rlhsl,savlen(nrank+1),rlhs,n,nsa)
         call fp_gatherv_real8_sav(rfwsl,savlen(nrank+1),rfws,n,nsa)
         call fp_gatherv_real8_sav(recsl,savlen(nrank+1),recs,n,nsa)
         call fp_gatherv_real8_sav(rwrsl,savlen(nrank+1),rwrs,n,nsa)
         call fp_gatherv_real8_sav(rwmsl,savlen(nrank+1),rwms,n,nsa)
         call fp_gatherv_real8_sav(rspbl,savlen(nrank+1),rspb,n,nsa)
         call fp_gatherv_real8_sav(rspfl,savlen(nrank+1),rspf,n,nsa)
         call fp_gatherv_real8_sav(rspsl,savlen(nrank+1),rsps,n,nsa)
         call fp_gatherv_real8_sav(rspll,savlen(nrank+1),rspl,n,nsa)
         call fp_gatherv_real8_sav(rpdrl,savlen(nrank+1),rpdr,n,nsa)
         call fp_gatherv_real8_sav(rndrl,savlen(nrank+1),rndr,n,nsa)
         call fp_gatherv_real8_sav(rtl_bulk,savlen(nrank+1),rt_bulk,n,nsa)
!         call fp_gatherv_real8_sav(rdidtl,savlen(nrank+1),rdidt,n,nsa)
         call fp_gatherv_real8_sav(rpssl,savlen(nrank+1),rpss,n,nsa)
         call fp_gatherv_real8_sav(rplsl,savlen(nrank+1),rpls,n,nsa)
         call fp_gatherv_real8_sav(rspsl_cx,savlen(nrank+1),rsps_cx,n,nsa)
         call fp_gatherv_real8_sav(tpsl,savlen(nrank+1),tps,n,nsa)


      end do

      do nsa=1,nsamax
        if(rns(1,nsa).lt.0)then
          do i=1,nrmax
            rns(1,nsa)=rns(i,nsa)
            if(rns(1,nsa).ge.0) goto 6003
          end do
        end if
6003    continue
        if(rns(nrmax,nsa).lt.0)then
          do i=1,nrmax
            rns(nrmax,nsa)=rns(nrmax-i+1,nsa)
            if(rns(nrmax,nsa).ge.0) cycle
          end do
        end if
      end do

        do nsa=1,nsamax
          do nr=2,nrmax-1
            if(rns(nr,nsa).lt.0)then
              do i=1,nrmax-nr
                rns(nr,nsa)=(rt_bulk(nr+i,nsa)+rns(nr-1,nsa))/2.d0
                if(rns(nr,nsa).ge.0) cycle
              end do
            end if
          end do
        end do

      work(:,:)=0.d0
      workg(:,:)=0.d0
      do nsb=1,nsbmax
         do nsa=nsastart,nsaend
            do nr=nrstart,nrend
               work(nr,nsa) = rpcs2l(nr,nsb,nsa)
            end do
         end do
         do n=1,nsw
            nsa=n+nsastart-1
            call fp_gatherv_real8_sav(work,savlen(nrank+1),workg,n,nsa)
         enddo

         do nsa=1,nsamax
            do nr=1,nrmax
               rpcs2(nr,nsb,nsa) = workg(nr,nsa)
            end do
         end do
      enddo

      work(:,:)=0.d0
      workg(:,:)=0.d0
      do nsb=1,nsbmax
         do nsa=nsastart,nsaend
            do nr=nrstart,nrend
               work(nr,nsa) = rpcs2l_del(nr,nsb,nsa)
            end do
         end do
         do n=1,nsw
            nsa=n+nsastart-1
            call fp_gatherv_real8_sav(work,savlen(nrank+1),workg,n,nsa)
         enddo

         do nsa=1,nsamax
            do nr=1,nrmax
               rpcs2_del(nr,nsb,nsa) = workg(nr,nsa)
            end do
         end do
      enddo
      call mtx_reset_communicator

      call mtx_broadcast_real8(rns,nrmax*nsamax)
      call mtx_broadcast_real8(rws,nrmax*nsamax)
      call mtx_broadcast_real8(rt_bulk,nrmax*nsamax)
!      call mtx_broadcast_real8(rfp,nrmax)

      if(modeld.ne.0)then
         temp_nsanr(:)=0.d0
         call mtx_set_communicator(comm_nsanr)
         call mtx_reduce_real8(rpdrs,nsamax,3,temp_nsanr,vloc)
         rpdrs(:)=temp_nsanr
         call mtx_reduce_real8(rndrs,nsamax,3,temp_nsanr,vloc)
         rndrs(:)=temp_nsanr
         call mtx_reset_communicator
      end if

      call mtx_set_communicator(comm_nr)
      do nsa=nsastart, nsaend
         do nr=nrstart, nrend
            do np=1,npmax
               temp_npnr1(np,nr) = rpl_bulk(np,nr,nsa)
            end do
         end do
         call mtx_allgather_real8(temp_npnr1,npmax*(nrend-nrstart+1),temp_npnr2)
         do nr=1, nrmax
            do np=1, npmax
               temp_npnr3(np,nr,nsa) = temp_npnr2(np,nr)
            end do
         end do
      end do
      call mtx_set_communicator(comm_nsa)
      call mtx_allgather_real8(temp_npnr3,npmax*nrmax*(nsaend-nsastart+1),rp_bulk)
      call mtx_reset_communicator

      end subroutine fpsavecomm
!^------------------------------
!==============================================================
!     update rns_delf, rws_para, rws_perp
      subroutine count_beam_density

      use libmtx
      use fpmpi
      implicit none
      integer:: nth, np, nr, nsa, ns, n, nsw
      double precision:: fact, rsum3_para, rsum3_perp

      rns_delf_nsa(:,:)=0.d0
      rws_delf_para(:,:)=0.d0
      rws_delf_perp(:,:)=0.d0
      rnsl_delf(:,:)=0.d0
      rwsl_para(:,:)=0.d0
      rwsl_perp(:,:)=0.d0

      call moment_0th_order(fnsp_del,rnsl_delf)
      call mtx_set_communicator(comm_np)

      do nsa=nsastart, nsaend
         ns=ns_nsa(nsa)
         do nr=nrstart, nrend
            rsum3_para=0.d0
            rsum3_perp=0.d0
!pressure
            if(modela.eq.0) then
               do np=npstart,npend
                  do nth=1,nthmax
                     rsum3_para = rsum3_para                   &
                          +volp(nth,np,ns)*fnsp_del(nth,np,nr,nsa) &
                          *(pm(np,ns)*cosm(nth))**2
                     rsum3_perp = rsum3_perp                   &
                          +volp(nth,np,ns)*fnsp_del(nth,np,nr,nsa) &
                          *0.5d0*(pm(np,ns)*sinm(nth))**2
                  end do
               enddo
            else ! modela=1
               do np=npstart,npend
                  do nth=1,nthmax
                     rsum3_para = rsum3_para                   &
                          +volp(nth,np,ns)*fnsp_del(nth,np,nr,nsa) &
                          *(pm(np,ns)*cosm(nth))**2*rlamda(nth,nr)*rfsadg(nr)
                     rsum3_perp = rsum3_perp                   &
                          +volp(nth,np,ns)*fnsp_del(nth,np,nr,nsa) &
                          *0.5d0*(pm(np,ns)*sinm(nth))**2*rlamda(nth,nr)*rfsadg(nr)
                  end do
               enddo
            end if

            call p_theta_integration(rsum3_para)
            call p_theta_integration(rsum3_perp)

            fact=rnfp0(nsa)*1.d20*ptfp0(nsa)**2/amfp(nsa)
            rwsl_para(nr,nsa) = rsum3_para*fact
            rwsl_perp(nr,nsa) = rsum3_perp*fact

         end do
      end do
      call mtx_reset_communicator

      call mtx_set_communicator(comm_nsanr)
      nsw=nsaend-nsastart+1
      do n=1,nsw
         nsa=n+nsastart-1
         call fp_gatherv_real8_sav(rnsl_delf,savlen(nrank+1),rns_delf_nsa,n,nsa)
         call fp_gatherv_real8_sav(rwsl_para,savlen(nrank+1),rws_delf_para,n,nsa)
         call fp_gatherv_real8_sav(rwsl_perp,savlen(nrank+1),rws_delf_perp,n,nsa)
      end do
      call mtx_reset_communicator

      rns_delf(:,:)=0.d0
      do nsa=1, nsamax
         ns=ns_nsa(nsa)
         do nr=1, nrmax
            rns_delf(nr,ns)=rns_delf_nsa(nr,nsa)
            rws_para(nr,ns)=rws_delf_para(nr,nsa)
            rws_perp(nr,ns)=rws_delf_perp(nr,nsa)
         end do
      end do

      call mtx_broadcast_real8(rns_delf,nrmax*nsmax)
      call mtx_broadcast_real8(rws_para,nrmax*nsmax)
      call mtx_broadcast_real8(rws_perp,nrmax*nsmax)

      end subroutine count_beam_density
!==============================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fpsavecomm2

      use libmtx
      use fpmpi
      implicit none
      integer:: nsa, nsw, n,nr

      call mtx_set_communicator(comm_nsanr)
      nsw=nsaend-nsastart+1
      do n=1,nsw
         nsa=n+nsastart-1
         call fp_gatherv_real8_sav(rnsl,savlen(nrank+1),rns,n,nsa)
         call fp_gatherv_real8_sav(rwsl,savlen(nrank+1),rws,n,nsa)
         call fp_gatherv_real8_sav(rtl_bulk,savlen(nrank+1),rt_bulk,n,nsa)
      end do
      call mtx_reset_communicator


      call mtx_broadcast_real8(rns,nrmax*nsamax)
      call mtx_broadcast_real8(rws,nrmax*nsamax)
      call mtx_broadcast_real8(rt_bulk,nrmax*nsamax)

      end subroutine fpsavecomm2
!==============================================================
      subroutine update_rn_rt(recv)

      use libmtx
      use fpmpi
      use eg_read
      use plprof
      implicit none
      integer:: nr, nsa, ns,i
      real(8):: rhon
      type(pl_plf_type),dimension(nsmax):: plf
      double precision,dimension(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend),intent(in):: recv

      call mtx_reset_communicator
      do nsa=nsastart,nsaend
         ns=ns_nsa(nsa)
         do nr=nrstart,nrend
            rhon=rm(nr)
            call pl_prof(rhon,plf)
            rn_temp(nr,ns)=plf(ns)%rn
            rt_temp(nr,ns)=(plf(ns)%rtpr+2.d0*plf(ns)%rtpp)/3.d0
         enddo
      end do

      call define_bulk_np
      call moment_0th_order(recv,rnsl)
      call moment_2nd_order(recv,rwsl)
      do nr=nrstart,nrend
         do nsa=nsastart,nsaend
            call bulk_temperature(np_bulk(nr,nsa),nr,nsa)
         enddo ! nsa
      enddo ! nr
      call mtx_reset_communicator
      call fpsavecomm2

      if(model_ex_read_tn.eq.0)then
         do ns=1,nsmax
            do nr=1,nrmax
               rhon=rm(nr)
               call pl_prof(rhon,plf)
               rn_temp(nr,ns)=plf(ns)%rn
               rt_temp(nr,ns)=(plf(ns)%rtpr+2.d0*plf(ns)%rtpp)/3.d0
            enddo
         end do
         do nsa=1,nsamax
            ns=ns_nsa(nsa)
            do nr=1,nrmax
               rn_temp(nr,ns) = rns(nr,nsa)
               rt_temp(nr,ns) = rt_bulk(nr,nsa)
            end do
         enddo
         do nsa=1,nsamax
           if(rt_temp(1,nsa).lt.0)then
             do i=1,nrmax
               rt_temp(1,nsa)=rt_temp(i,nsa)
               if(rt_temp(1,nsa).ge.0) goto 5001
             end do
           end if
5001       continue
           if(rt_temp(nrmax,nsa).lt.0)then
             do i=1,nrmax
               rt_temp(nrmax,nsa)=rt_temp(nrmax-i+1,nsa)
               if(rt_temp(nrmax,nsa).ge.0) cycle
             end do
           end if
         end do

           do nsa=1,nsamax
             do nr=2,nrmax-1
               if(rt_temp(nr,nsa).lt.0)then
                 do i=1,nrmax-nr
                   rt_temp(nr,nsa)=(rt_temp(nr+i,nsa)+rt_temp(nr-1,nsa))/2.d0
                   if(rt_temp(nr,nsa).ge.0) cycle
                 end do
               end if
             end do
           end do

           do nsa=1,nsamax
             if(rn_temp(1,nsa).lt.0)then
               do i=1,nrmax
                 rn_temp(1,nsa)=rn_temp(i,nsa)
                 if(rn_temp(1,nsa).ge.0) goto 5002
               end do
             end if
5002         continue
             if(rn_temp(nrmax,nsa).lt.0)then
               do i=1,nrmax
                 rn_temp(nrmax,nsa)=rn_temp(nrmax-i+1,nsa)
                 if(rn_temp(nrmax,nsa).ge.0) cycle
               end do
             end if
           end do

             do nsa=1,nsamax
               do nr=2,nrmax-1
                 if(rn_temp(nr,nsa).lt.0)then
                   do i=1,nrmax-nr
                     rn_temp(nr,nsa)=(rn_temp(nr+i,nsa)+rn_temp(nr-1,nsa))/2.d0
                     if(rn_temp(nr,nsa).ge.0) cycle
                   end do
                 end if
               end do
             end do

      elseif(model_ex_read_tn.ne.0)then
         call make_exp_prof(timefp+delt)
      end if

      end subroutine update_rn_rt
!==============================================================
      subroutine bulk_temperature(npb,nr,nsa)

        use fpsub
      use fpmpi
      implicit none
      integer,intent(in):: npb, nr, nsa
!      real(8),intent(out):: rtl
      integer:: isw_bulk, np, nth, ns
      real(8):: rsum_t, rsum_v, pv, dfdp, wpl, ffp, rsumn, rsumw, fact
      real(8),dimension(nthmax,npmax):: t_bulk
      real(8),dimension(npstart:npend):: rpl_bulk_send
      real(8),dimension(npmax):: rpl_bulk_recv
      real(8):: rnl_bulk, rwl_bulk, rtemp

!      isw_bulk=0 ! sometimes dfdp becomes 0 and
!                   then it makes density nan. (for fact_bulk < 4)
      isw_bulk=1  ! requires higher fact_bulk (fact_bulk >= 4)
!                   to obtain accurate rt_bulk

      ns=ns_nsa(nsa)

!      rtl_bulk(:,:)=0.d0
!      rpl_bulk(:,:,:)=0.d0

      call mtx_set_communicator(comm_np)
      if(model_ex_read_tn.eq.0)then
      if(isw_bulk.eq.0)then ! bulk t calculation using dfdp
         rsum_t=0.d0
         rsum_v=0.d0
         do np=npstart,npend
            if(np.ge.2.and.np.le.npb)then
               pv=sqrt(1.d0+theta0(ns)*pg(np,ns)**2)
               do nth=1,nthmax
                  if(fnsp(nth,np,nr,nsa).gt.0.d0.and. &
                     fnsp(nth,np-1,nr,nsa).gt.0.d0)then
                     dfdp=delp(ns) &
                          /( log(fnsp(nth,np  ,nr,nsa)) &
                            -log(fnsp(nth,np-1,nr,nsa)) )
                  else
                     wpl=weighp(nth  ,np,nr,nsa)
                     ffp=   ( (1.d0-wpl)*fnsp(nth  ,np  ,nr,nsa)  &
                                   +wpl *fnsp(nth  ,np-1,nr,nsa) )
                     dfdp=delp(ns)*ffp/(                         &
                          fnsp(nth,np,nr,nsa)-fnsp(nth,np-1,nr,nsa) )
                  end if
                  if(model_disrupt.ne.1)then
                     if(dfdp.ne.dfdp)then
                                   ! nan never equals to any other variables.
!                              "dfdp is nan in fpsave. timefp= ", &
!                               timefp, np, nth, nsa, &
!                               fnsp(nth,np,nr,nsa),fnsp(nth,np-1,nr,nsa), &
!                               fnsp(nth,np,nr,nsa)-fnsp(nth,np-1,nr,nsa)
                        dfdp=0.d0
                     end if
                     if(dfdp.gt.0.d0)then
!                        write(16,'(a,e14.6,3i4,3e14.6)') &
!                             "dfdp is positive in fpsave. timefp= ", &
!                             timefp, np, nth, nsa, &
!                             fnsp(nth,np,nr,nsa),fnsp(nth,np-1,nr,nsa), &
!                             fnsp(nth,np,nr,nsa)-fnsp(nth,np-1,nr,nsa)
                        dfdp=0.d0
                     end if
                  end if

                  t_bulk(nth,np)=-pg(np,ns)*ptfp0(nsa)*dfdp &
                       /aee/1.d3*vtfp0(nsa)/pv
                  rsum_t = rsum_t + t_bulk(nth,np)*volp(nth,np,ns)
                  rsum_v = rsum_v + volp(nth,np,ns)
               end do
               rpl_bulk_send(np) = rsum_t/rsum_v
            else
               rpl_bulk_send(np) = 0.d0
            end if
         end do
         call mtx_allgather_real8(rpl_bulk_send,npend-npstart+1,rpl_bulk_recv)
         do np=1, npmax
            rpl_bulk(np,nr,nsa) = rpl_bulk_recv(np)
         end do
         call p_theta_integration(rsum_t)
         call p_theta_integration(rsum_v)

         rtl_bulk(nr,nsa)=rsum_t/rsum_v
      else ! isw_bulk=1  bulk t using t=w/n
         rsumn=0.d0
         rsumw=0.d0

         if(modela.eq.0)then
            do np=npstart,npend
               if(np.le.npb)then
                  do nth=1,nthmax
                     rsumn = rsumn+volp(nth,np,ns)*fnsp(nth,np,nr,nsa)
                  end do
               end if
            enddo
         else
            do np=npstart,npend
               if(np.le.npb)then
                  do nth=1,nthmax
                     rsumn = rsumn+volp(nth,np,ns)*fnsp(nth,np,nr,nsa) &
                          *rlamda(nth,nr)*rfsadg(nr)
                  end do
               end if
            enddo
         end if
         if(modela.eq.0) then
            if(modelr.eq.0) then
               do np=npstart,npend
                  if(np.le.npb)then
                     do nth=1,nthmax
                        rsumw = rsumw                       &
                             +volp(nth,np,ns)*fnsp(nth,np,nr,nsa) &
                             *0.5d0*pm(np,ns)**2
                     end do
                  end if
               enddo
            else
               do np=npstart,npend
                  if(np.le.npb)then
                     pv=sqrt(1.d0+theta0(ns)*pm(np,ns)**2)
                     do nth=1,nthmax
                        rsumw = rsumw                       &
                             +volp(nth,np,ns)*fnsp(nth,np,nr,nsa) &
                             *(pv-1.d0)/theta0(ns)
                     end do
                  end if
               end do
            endif
         else ! modela=1
            if(modelr.eq.0) then
               do np=npstart,npend
                  if(np.le.npb)then
                     do nth=1,nthmax
                        rsumw = rsumw                        &
                             +volp(nth,np,ns)*fnsp(nth,np,nr,nsa)  &
                             *0.5d0*pm(np,ns)**2*rlamda(nth,nr)*rfsadg(nr)
                     end do
                  end if
               enddo
            else
               do np=npstart,npend
                  if(np.le.npb)then
                     pv=sqrt(1.d0+theta0(ns)*pm(np,ns)**2)
                     do nth=1,nthmax
                        rsumw = rsumw                        &
                             +volp(nth,np,ns)*fnsp(nth,np,nr,nsa)  &
                             *(pv-1.d0)/theta0(ns)*rlamda(nth,nr)*rfsadg(nr)
                     end do
                  end if
               end do
            endif
         end if
         call p_theta_integration(rsumn)
         call p_theta_integration(rsumw)

         fact=rnfp0(nsa)*1.d20
         rnl_bulk = rsumn*fact*1.d-20
         fact=rnfp0(nsa)*1.d20*ptfp0(nsa)**2/amfp(nsa)
         rwl_bulk = rsumw*fact*1.d-6

         if(rnl_bulk.ne.0.d0) then
            if(modelr.eq.0)then
               rtl_bulk(nr,nsa) = rwl_bulk*1.d6 &
                    /(1.5d0*rnl_bulk*1.d20*aee*1.d3)
            elseif(modelr.eq.1)then
               call fpnewton(nr,nsa,rnl_bulk,rwl_bulk,rtemp)
               rtl_bulk(nr,nsa) = rtemp
            end if
         end if

         do np=1, npmax
            rpl_bulk(np,nr,nsa) = 0.d0
         end do
      end if
      else ! model_ex_read_tn!=0
         rtl_bulk(nr,nsa)=rt_read(nr,ns)
      end if
      call mtx_reset_communicator

      end subroutine bulk_temperature
!==============================================================
      subroutine define_bulk_np

      use fpmpi
      use libmtx
      implicit none
      integer:: np, nr, nsa, ns,i
      double precision:: pmax_bulk, p_bulk_r, rhon, rtfpl
      type(pl_plf_type),dimension(nsmax):: plf

      do nsa=1,nsamax
        if(rt_bulk(1,nsa).lt.0)then
          do i=1,nrmax
            rt_bulk(1,nsa)=rt_bulk(i,nsa)
            if(rt_bulk(1,nsa).gt.0) goto 5003
          end do
        end if
5003    continue
        if(rt_bulk(nrmax,nsa).lt.0)then
          do i=1,nrmax
            rt_bulk(nrmax,nsa)=rt_bulk(nrmax-i+1,nsa)
            if(rt_bulk(nrmax,nsa).gt.0) cycle
          end do
        end if
      end do

        do nsa=1,nsamax
          do nr=2,nrmax-1
            if(rt_bulk(nr,nsa).lt.0)then
              do i=1,nrmax-nr
                rt_bulk(nr,nsa)=(rt_bulk(nr+i,nsa)+rt_bulk(nr-1,nsa))/2.d0
                if(rt_bulk(nr,nsa).gt.0) cycle
              end do
            end if
          end do
        end do


!     define bulk momentum range: 0 < np < np_bulk(nr,nsa)
      if(model_bulk_const.eq.0)then
         do nsa=nsastart,nsaend
            ns=ns_nsa(nsa)
            do nr=nrstart,nrend
               if(nt_init.eq.0)then
                  pmax_bulk = fact_bulk
                  p_bulk_r = pmax_bulk* rt_temp(nr,ns)/rtfp0(nsa)
               else
                  rhon=rm(nr)
                  call pl_prof(rhon,plf)
                  rtfpl=(plf(ns)%rtpr+2.d0*plf(ns)%rtpp)/3.d0
                  pmax_bulk = sqrt( rt_bulk(nr,nsa)/rtfpl )*fact_bulk
                  p_bulk_r = pmax_bulk* rt_bulk(nr,nsa)/rtfp0(nsa)
               end if
               np_bulk(nr,nsa) = npmax
               do np=npmax, 1, -1
                  if(p_bulk_r.lt.pg(np,ns)) np_bulk(nr,nsa)=np
               end do
            end do
         end do
      else
         do nsa=nsastart, nsaend
            ns=ns_nsa(nsa)
            do nr=nrstart, nrend
               pmax_bulk = fact_bulk
               if(model_ex_read_tn.eq.0)then
                  p_bulk_r = pmax_bulk* rtfp(nr,nsa)/rtfp0(nsa)
               else
                  p_bulk_r = pmax_bulk* rt_temp(nr,ns)/rtfp0(nsa)
               end if
               np_bulk(nr,nsa) = npmax
               do np=npmax, 1, -1
                  if(p_bulk_r.lt.pg(np,ns)) np_bulk(nr,nsa)=np
               end do
            end do
         end do
      end if

      end subroutine define_bulk_np
!==============================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine moment_0th_order(send,recv)
!     send=fnsp, fnsp_del,    recv=rnsl, rnsl_del

      use libmtx
      use fpmpi
      implicit none
      integer:: np, nth, nr, nsa, ns
      double precision:: rsum1, fact
      double precision,dimension(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend),intent(in)::send
      double precision,dimension(nrstart:nrend,nsamax),intent(out):: recv

      call mtx_set_communicator(comm_np)

      do nsa=nsastart,nsaend
         ns=ns_nsa(nsa)
         do nr=nrstart,nrend
            rsum1=0.d0
            if(modela.eq.0)then
               do np=npstart,npend
                  do nth=1,nthmax
                     rsum1 = rsum1+volp(nth,np,ns)*send(nth,np,nr,nsa)
                  end do
               enddo
            else
               do np=npstart,npend
                  do nth=1,nthmax
                     rsum1 = rsum1+volp(nth,np,ns)*send(nth,np,nr,nsa) &
                          *rlamda(nth,nr)*rfsadg(nr)
                  end do
               enddo
            end if

            call p_theta_integration(rsum1)

            fact=rnfp0(nsa)*1.d20
            recv(nr,nsa) = rsum1*fact*1.d-20
!            rnsl(nr,nsa) = rsum1*fact*1.d-20
         end do
      end do

      call mtx_reset_communicator

      end subroutine moment_0th_order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine moment_1st_order(send,recv)

      use libmtx
      use fpmpi
      implicit none
      integer:: np, nth, nr, nsa, ns
      double precision:: fact, rsum2, pv
      double precision,dimension(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend),intent(in)::send
      double precision,dimension(nrstart:nrend,nsamax),intent(out):: recv

      call mtx_set_communicator(comm_np)

      do nsa=nsastart,nsaend
         ns=ns_nsa(nsa)
         do nr=nrstart,nrend
            rsum2=0.d0

            if(modela.eq.0) then
               if(modelr.eq.0) then
                  do np=npstart,npend
                     do nth=1,nthmax
                        rsum2 = rsum2                       &
                             +volp(nth,np,ns)*send(nth,np,nr,nsa) &
                             *pm(np,ns)*cosm(nth)
                     end do
                  enddo
               else
                  do np=npstart,npend
                     pv=sqrt(1.d0+theta0(ns)*pm(np,ns)**2)
                     do nth=1,nthmax
                        rsum2 = rsum2                       &
                             +volp(nth,np,ns)*send(nth,np,nr,nsa) &
                             *pm(np,ns)*cosm(nth)/pv
                     end do
                  end do
               endif
            else ! modela=1
               if(modelr.eq.0) then
                  do np=npstart,npend
                     do nth=1,nthmax
                        rsum2 = rsum2                        &
                             +volp(nth,np,ns)*send(nth,np,nr,nsa)  &
                             *pm(np,ns)*cosm(nth)*rlamda(nth,nr)*rfsadg(nr)
                     end do
                  enddo
               else
                  do np=npstart,npend
                     pv=sqrt(1.d0+theta0(ns)*pm(np,ns)**2)
                     do nth=1, itl(nr)-1
                        rsum2 = rsum2                        &
                             +volp(nth,np,ns)*send(nth,np,nr,nsa)  &
                             *pm(np,ns)*cosm(nth)/pv &
                             *(1.d0+epsrm2(nr))
                     end do
                     do nth=itu(nr)+1, nthmax
                        rsum2 = rsum2                        &
                             +volp(nth,np,ns)*send(nth,np,nr,nsa)  &
                             *pm(np,ns)*cosm(nth)/pv &
                             *(1.d0+epsrm2(nr))
                     end do
                  end do
               endif
            end if

            call p_theta_integration(rsum2)

            fact=rnfp0(nsa)*1.d20
            recv(nr,nsa) = rsum2*fact*aefp(nsa)*ptfp0(nsa) &
                           /amfp(nsa)*1.d-6
         end do
      end do


      call mtx_reset_communicator
      end subroutine moment_1st_order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine moment_2nd_order(send,recv)

      use libmtx
      use fpmpi
      implicit none
      integer:: np, nth, nr, nsa, ns
      double precision:: fact, rsum3, pv
      double precision,dimension(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend),intent(in)::send
      double precision,dimension(nrstart:nrend,nsamax),intent(out):: recv

      call mtx_set_communicator(comm_np)

      do nsa=nsastart,nsaend
         ns=ns_nsa(nsa)
         do nr=nrstart,nrend
            rsum3=0.d0

            if(modela.eq.0) then
               if(modelr.eq.0) then
                  do np=npstart,npend
                     do nth=1,nthmax
                        rsum3 = rsum3                       &
                             +volp(nth,np,ns)*send(nth,np,nr,nsa) &
                             *0.5d0*pm(np,ns)**2
                     end do
                  enddo
               else
                  do np=npstart,npend
                     pv=sqrt(1.d0+theta0(ns)*pm(np,ns)**2)
                     do nth=1,nthmax
                        rsum3 = rsum3                       &
                             +volp(nth,np,ns)*send(nth,np,nr,nsa) &
                             *(pv-1.d0)/theta0(ns)
                     end do
                  end do
               endif
            else ! modela=1
               if(modelr.eq.0) then
                  do np=npstart,npend
                     do nth=1,nthmax
                        rsum3 = rsum3                        &
                             +volp(nth,np,ns)*send(nth,np,nr,nsa)  &
                             *0.5d0*pm(np,ns)**2*rlamda(nth,nr)*rfsadg(nr)
                     end do
                  enddo
               else
                  do np=npstart,npend
                     pv=sqrt(1.d0+theta0(ns)*pm(np,ns)**2)
                     do nth=1,nthmax
                        rsum3 = rsum3                        &
                             +volp(nth,np,ns)*send(nth,np,nr,nsa)  &
                             *(pv-1.d0)/theta0(ns)*rlamda(nth,nr)*rfsadg(nr)
                     end do
                  end do
               endif
            end if

            call p_theta_integration(rsum3)

            fact=rnfp0(nsa)*1.d20*ptfp0(nsa)**2/amfp(nsa)
            recv(nr,nsa) = rsum3*fact               *1.d-6
         end do
      end do

      call mtx_reset_communicator
      end subroutine moment_2nd_order
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine power_from_diffusion_coef

      use libmtx
      use fpmpi
      implicit none
      integer:: np, nth, nr, nsa, ns, nps, nsb
      double precision:: fact, dfp, dft, ffp
      double precision,dimension(nthmax,npstart:npendwg,nrstart:nrendwm,nsamax):: no_coef
      double precision:: rsum1, rsum2, rsum3, rsum4, rsum5, rsum6, rsum8, rsum9, rsum10
      double precision,dimension(nsbmax):: rsum11, rsum12
      double precision:: rsum_wr, rsum_wm

      no_coef(:,:,:,:)=0.d0
      call mtx_set_communicator(comm_np)

      if(npstart.eq.1)then
         nps=2
      else
         nps=npstart
      end if
      do nsa=nsastart,nsaend
         ns=ns_nsa(nsa)
         do nr=nrstart,nrend
            rsum1=0.d0
            rsum2=0.d0
            rsum3=0.d0
            rsum4=0.d0
            rsum5=0.d0
            rsum6=0.d0
            rsum8=0.d0
            rsum9=0.d0
            rsum10=0.d0
            rsum11(:)=0.d0
            rsum12(:)=0.d0
            rsum_wr=0.d0
            rsum_wm=0.d0
            do np=nps,npend
               do nth=1,nthmax
                  call phase_space_derivative_f(nth,np,nr,nsa,fnsp,dfp,dft,ffp)

                  call integrand_in_m2od(nth,np,nr,nsa,dfp,dft,ffp,dcpp,dcpt,fcpp,rsum1)
                  call integrand_in_m2od(nth,np,nr,nsa,dfp,dft,ffp,dwpp,dwpt,no_coef,rsum2)
                  call integrand_in_m2od(nth,np,nr,nsa,dfp,dft,ffp,no_coef,no_coef,fepp,rsum3)
                  call integrand_in_m2od(nth,np,nr,nsa,dfp,dft,ffp,dwlhpp,dwlhpt,no_coef,rsum4)
                  call integrand_in_m2od(nth,np,nr,nsa,dfp,dft,ffp,dwfwpp,dwfwpt,no_coef,rsum5)
                  call integrand_in_m2od(nth,np,nr,nsa,dfp,dft,ffp,dwecpp,dwecpt,no_coef,rsum6)
                  call integrand_in_m2od(nth,np,nr,nsa,dfp,dft,ffp,no_coef,no_coef,fspp,rsum8)
                  call integrand_in_m2od(nth,np,nr,nsa,dfp,dft,ffp,dlpp,no_coef,flpp,rsum9)
                  call integrand_in_m2od(nth,np,nr,nsa,dfp,dft,ffp,dpp,dpt,fpp,rsum10)
                  call integrand_in_m2od(nth,np,nr,nsa,dfp,dft,ffp,dwwrpp,dwwrpt,no_coef,rsum_wr)
                  call integrand_in_m2od(nth,np,nr,nsa,dfp,dft,ffp,dwwmpp,dwwmpt,no_coef,rsum_wm)
                  do nsb=1,nsbmax
                     call integrand_in_m2od_nsb(nth,np,nr,nsa,nsb,dfp,dft,ffp,dcpp2,dcpt2,fcpp2,rsum11)
                  end do
               end do
            end do
            call p_theta_integration(rsum1)
            call p_theta_integration(rsum2)
            call p_theta_integration(rsum3)
            call p_theta_integration(rsum4)
            call p_theta_integration(rsum5)
            call p_theta_integration(rsum6)
            call p_theta_integration(rsum8)
            call p_theta_integration(rsum9)
            call p_theta_integration(rsum10)
            call p_theta_integration(rsum_wr)
            call p_theta_integration(rsum_wm)
            do nsb=1,nsbmax
               call p_theta_integration(rsum11(nsb))
            end do

            fact=rnfp0(nsa)*1.d20*ptfp0(nsa)**2/amfp(nsa)*rfsadg(nr)

            rpcsl(nr,nsa)=-rsum1*fact*2.d0*pi*delp(ns)*delth *1.d-6
            rpwsl(nr,nsa)=-rsum2*fact*2.d0*pi*delp(ns)*delth *1.d-6
            rpesl(nr,nsa)=-rsum3*fact*2.d0*pi*delp(ns)*delth *1.d-6
            rlhsl(nr,nsa)=-rsum4*fact*2.d0*pi*delp(ns)*delth *1.d-6
            rfwsl(nr,nsa)=-rsum5*fact*2.d0*pi*delp(ns)*delth *1.d-6
            recsl(nr,nsa)=-rsum6*fact*2.d0*pi*delp(ns)*delth *1.d-6
            rpssl(nr,nsa)=-rsum8*fact*2.d0*pi*delp(ns)*delth *1.d-6
            rplsl(nr,nsa)=-rsum9*fact*2.d0*pi*delp(ns)*delth *1.d-6
            rws123l(nr,nsa) =-rsum10*fact*2.d0*pi*delp(ns)*delth *1.d-6
            rwrsl(nr,nsa)=-rsum_wr*fact*2.d0*pi*delp(ns)*delth *1.d-6
            rwmsl(nr,nsa)=-rsum_wm*fact*2.d0*pi*delp(ns)*delth *1.d-6
            do nsb=1,nsbmax
               rpcs2l(nr,nsb,nsa)=-rsum11(nsb) &
                    *fact*2.d0*pi*delp(ns)*delth *1.d-6
            end do

!      for delta f
            if(model_delta_f(ns).eq.1)then
!      collisional power transfer delta f
               do np=nps,npend
                  do nth=1,nthmax
                     call phase_space_derivative_f(nth,np,nr,nsa,fnsp_del,dfp,dft,ffp)
                     do nsb=1,nsbmax
                        call integrand_in_m2od_nsb(nth,np,nr,nsa,nsb,dfp,dft,ffp,dcpp2,dcpt2,fcpp2,rsum12)
                     end do
                  end do
               end do
               do nsb=1,nsbmax
                  call p_theta_integration(rsum12(nsb))
               end do
               do nsb=1,nsbmax
                  rpcs2l_del(nr,nsb,nsa)=-rsum12(nsb) &
                       *fact*2.d0*pi*delp(ns)*delth *1.d-6
               end do
            end if

         end do ! nr
      end do ! nsa

      call mtx_reset_communicator

      end subroutine power_from_diffusion_coef
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine phase_space_derivative_f(nth,np,nr,nsa,send,dfp,dft,ffp)

      implicit none
      integer,intent(in):: nth,np,nr,nsa
      integer:: ns
      double precision,dimension(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend),intent(in)::send
      double precision,intent(out):: dfp, dft, ffp
      double precision:: wpp, wpm, wpl

      ns=ns_nsa(nsa)

      wpl=weighp(nth  ,np,nr,nsa)
      if(nth.eq.1) then
         wpm=0.d0
      else
         wpm=weighp(nth-1,np,nr,nsa)
      endif
      if(nth.eq.nthmax) then
         wpp=0.d0
      else
         wpp=weighp(nth+1,np,nr,nsa)
      endif
      dfp=    pg(np,ns) &
           /delp(ns)*(send(nth,np,nr,nsa)-send(nth,np-1,nr,nsa))
      if(nth.eq.1) then

         dft=1.d0/delth                             &
              *(                                     &
              ((1.d0-wpp)*send(nth+1,np  ,nr,nsa)   &
              +wpp *send(nth+1,np-1,nr,nsa)) &
              -                                    &
              ((1.d0-wpm)*send(nth,np  ,nr,nsa)     &
              +wpm *send(nth,np-1,nr,nsa)) &
              )

      else if(nth.eq.nthmax) then
         dft=    1.d0/delth                         &
              *(-                                    &
              ((1.d0-wpm)*send(nth-1,np  ,nr,nsa)   &
              +wpm *send(nth-1,np-1,nr,nsa)) &
              +                                     &
              ((1.d0-wpp)*send(nth,np  ,nr,nsa)     &
              +wpp *send(nth,np-1,nr,nsa)) &
              )
      else
         dft=    1.d0/(2.d0*delth)                  &
              *(                                     &
              ((1.d0-wpp)*send(nth+1,np  ,nr,nsa)   &
              +wpp *send(nth+1,np-1,nr,nsa)) &
              -                                    &
              ((1.d0-wpm)*send(nth-1,np  ,nr,nsa)   &
              +wpm *send(nth-1,np-1,nr,nsa)) &
              )
      endif
      ffp=    pg(np,ns)                           &
           *((1.d0-wpl)*send(nth  ,np  ,nr,nsa)  &
           +wpl *send(nth  ,np-1,nr,nsa))

      end subroutine phase_space_derivative_f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine integrand_in_m2od(nth,np,nr,nsa,dfp,dft,ffp,difpp,difpt,fricp,rsum)

      implicit none
      integer,intent(in):: nth,np,nr,nsa
      integer:: ns
      double precision,intent(in):: dfp,dft,ffp
      double precision,dimension(nthmax,npstart :npendwg,nrstart:nrendwm,nsamax),intent(in)::difpp, difpt, fricp
      double precision,intent(inout):: rsum
      double precision:: pv

      ns=ns_nsa(nsa)
      pv=sqrt(1.d0+theta0(ns)*pg(np,ns)**2)
      rsum = rsum+pg(np,ns)**2*sinm(nth)/pv   &
           *(difpp(nth,np,nr,nsa)*dfp           &
           +difpt(nth,np,nr,nsa)*dft           &
           -fricp(nth,np,nr,nsa)*ffp           &
           )

      end subroutine integrand_in_m2od
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine integrand_in_m2od_nsb(nth,np,nr,nsa,nsb,dfp,dft,ffp,difpp,difpt,fricp,rsum)

      implicit none
      integer,intent(in):: nth,np,nr,nsa,nsb
      integer:: ns
      double precision,intent(in):: dfp,dft,ffp
      double precision,dimension(nthmax,npstart :npendwg,nrstart:nrendwm,nsbmax,nsamax),intent(in)::difpp, difpt, fricp
      double precision,dimension(nsbmax),intent(inout):: rsum
      double precision:: pv

      ns=ns_nsa(nsa)
      pv=sqrt(1.d0+theta0(ns)*pg(np,ns)**2)
      rsum(nsb) = rsum(nsb)+pg(np,ns)**2*sinm(nth)/pv   &
           *(difpp(nth,np,nr,nsb,nsa)*dfp           &
           +difpt(nth,np,nr,nsb,nsa)*dft           &
           -fricp(nth,np,nr,nsb,nsa)*ffp           &
           )

      end subroutine integrand_in_m2od_nsb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine power_from_source_term

      use fpmpi
      use libmpi
      use libmtx
      implicit none
      integer:: nth, np, nr, nsa, ns,k
      double precision:: fact, rsum11b, rsum11f, rsum11s, rsum11l, rsum11s_cx, pv
      call mtx_set_communicator(comm_np)
      do nr=nrstart,nrend
         do nsa=nsastart,nsaend
            ns=ns_nsa(nsa)
            rsum11b=0.d0
            rsum11f=0.d0
            rsum11s=0.d0
            rsum11l=0.d0
            rsum11s_cx=0.d0

            if(modelr.eq.1)then
               do np=npstart,npend
                  pv=sqrt(1.d0+theta0(ns)*pm(np,ns)**2)
                  do nth=1,nthmax
                     rsum11b = rsum11b + pm(np,ns)**2*sinm(nth) &
                          *(pv-1.d0)/theta0(ns)*sppb(nth,np,nr,nsa)*rlamda(nth,nr)!*rfsadg(nr)
                     rsum11f = rsum11f + pm(np,ns)**2*sinm(nth) &
                          *(pv-1.d0)/theta0(ns)*sppf(nth,np,nr,nsa)*rlamda(nth,nr)!*rfsadg(nr)
                     rsum11s = rsum11s + pm(np,ns)**2*sinm(nth) &
                          *(pv-1.d0)/theta0(ns)*sppl(nth,np,nr,nsa)*rlamda(nth,nr)!*rfsadg(nr)
                     if(model_delta_f(ns).eq.0)then
                        rsum11l = rsum11l + pm(np,ns)**2*sinm(nth) &
                             *(pv-1.d0)/theta0(ns)* ppl(nth,np,nr,nsa)&
                             *fnsp(nth,np,nr,nsa)!*rlamda(nth,nr)!*rfsadg(nr) ppl includes rlamda
                     elseif(model_delta_f(ns).eq.1)then
                        rsum11l = rsum11l + pm(np,ns)**2*sinm(nth) &
                             *(pv-1.d0)/theta0(ns)* ppl(nth,np,nr,nsa)&
                             *fnsp_del(nth,np,nr,nsa)!*rlamda(nth,nr)!*rfsadg(nr) ppl includes rlamda
                     end if
                     rsum11s_cx = rsum11s_cx + pm(np,ns)**2*sinm(nth) &
                          *(pv-1.d0)/theta0(ns)*sppl_cx(nth,np,nr,nsa)*rlamda(nth,nr)!*rfsadg(nr)
                  end do
               end do
            else ! modelr=0
               do np=npstart,npend
                  do nth=1,nthmax
                     rsum11b = rsum11b + pm(np,ns)**2*sinm(nth) &
                          *0.5d0*pm(np,ns)**2*sppb(nth,np,nr,nsa)*rlamda(nth,nr)!*rfsadg(nr)
                     rsum11f = rsum11f + pm(np,ns)**2*sinm(nth) &
                          *0.5d0*pm(np,ns)**2*sppf(nth,np,nr,nsa)*rlamda(nth,nr)!*rfsadg(nr)
                     rsum11s = rsum11s + pm(np,ns)**2*sinm(nth) &
                          *0.5d0*pm(np,ns)**2*sppl(nth,np,nr,nsa)*rlamda(nth,nr)!*rfsadg(nr)
                     if(model_delta_f(ns).eq.0)then
                        rsum11l = rsum11l + pm(np,ns)**2*sinm(nth) &
                             *0.5d0*pm(np,ns)**2*ppl(nth,np,nr,nsa) &
                             *fnsp(nth,np,nr,nsa)!*rlamda(nth,nr)!*rfsadg(nr)
                     elseif(model_delta_f(ns).eq.1)then
                        rsum11l = rsum11l + pm(np,ns)**2*sinm(nth) &
                             *0.5d0*pm(np,ns)**2*ppl(nth,np,nr,nsa) &
                             *fnsp_del(nth,np,nr,nsa)!*rlamda(nth,nr)!*rfsadg(nr)
                     end if
                     rsum11s_cx = rsum11s_cx + pm(np,ns)**2*sinm(nth) &
                          *0.5d0*pm(np,ns)**2*sppl_cx(nth,np,nr,nsa)*rlamda(nth,nr)!*rfsadg(nr)
                  end do
               end do
            end if
            call p_theta_integration(rsum11b) ! beam power
            call p_theta_integration(rsum11f) ! fusion power
            call p_theta_integration(rsum11s) ! tloss power
            call p_theta_integration(rsum11l) ! ppl power including tloss and thermalization
            call p_theta_integration(rsum11s_cx) ! cx loss power

            fact=rnfp0(nsa)*1.d20*ptfp0(nsa)**2/amfp(nsa)*rfsadg(nr)
            rspbl(nr,nsa)= rsum11b*fact*2.d0*pi*delp(ns)*delth*1.d-6
            rspfl(nr,nsa)= rsum11f*fact*2.d0*pi*delp(ns)*delth*1.d-6
            rspsl(nr,nsa)= rsum11s*fact*2.d0*pi*delp(ns)*delth*1.d-6
            rspll(nr,nsa)= rsum11l*fact*2.d0*pi*delp(ns)*delth*1.d-6
            rspsl_cx(nr,nsa)= rsum11s_cx*fact*2.d0*pi*delp(ns)*delth*1.d-6
         end do
      end do
      end subroutine power_from_source_term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine power_from_radial_transport

      use fpmpi
      use libmpi
      implicit none
      integer:: nth, np, nr, nsa, ns
      real(8):: wrl, wrh, dint_dfdt_r1, dint_dfdt_r2
      real(8):: dfdr_r1, dfdr_r2, dfdt_r1, dfdt_r2, dint_dr, rsum_dr,rgama,f_r1,f_r2, rsumn_dr
      real(8):: srhor1, srhor2, rsum_drs, rsumn_drs
      real(8)::drrm, drrp, rl

      if(modeld.ne.0)then
         call mtx_set_communicator(comm_np)
         do nr=nrstart,nrend
            do nsa=nsastart,nsaend
               ns=ns_nsa(nsa)

               rsum_dr=0.d0
               rsum_drs=0.d0
               rsumn_dr=0.d0
               rsumn_drs=0.d0
               dint_dfdt_r1=0.d0
               dint_dfdt_r2=0.d0
               srhor1 = 0.d0
               srhor2 = 0.d0
               if(modela.eq.0)then
                  drrm=rg(nr)
                  drrp=rg(nr+1)
                  rl=rm(nr)
               else
                  drrm=1.d0
                  drrp=1.d0
                  rl=1.d0
               end if
               do np=npstart,npend
                  rgama=sqrt(1.d0+theta0(ns)*pm(np,ns)**2)
                  do nth=1,nthmax
                     wrl=weighr(nth,np,nr,nsa)
                     wrh=weighr(nth,np,nr+1,nsa)
                     if(nr.ne.1.and.nr.ne.nrmax)then
                        dfdr_r1 = ( fnsp(nth,np,nr,nsa)-fnsp(nth,np,nr-1,nsa) ) / delr
                        f_r1 = ( (1.d0-wrl)*fnsp(nth,np,nr,nsa) + wrl*fnsp(nth,np,nr-1,nsa) )

                        dfdr_r2 = ( fnsp(nth,np,nr+1,nsa)-fnsp(nth,np,nr,nsa) ) / delr
                        f_r2 = ( (1.d0-wrh)*fnsp(nth,np,nr+1,nsa) + wrh*fnsp(nth,np,nr,nsa) )

                        dfdt_r1 = ( drr(nth,np,nr,nsa)*dfdr_r1 - frr(nth,np,nr,nsa)*f_r1 )*drrm
                        dfdt_r2 = ( drr(nth,np,nr+1,nsa)*dfdr_r2 - frr(nth,np,nr+1,nsa)*f_r2)*drrp
                     elseif(nr.eq.1)then
                        f_r2 = ( (1.d0-wrh)*fnsp(nth,np,nr+1,nsa) + wrh*fnsp(nth,np,nr,nsa) )
                        dfdr_r2 = ( fnsp(nth,np,nr+1,nsa)-fnsp(nth,np,nr,nsa) ) / delr

                        dfdt_r1 = 0.d0
                        dfdt_r2 = ( drr(nth,np,nr+1,nsa)*dfdr_r2 - frr(nth,np,nr+1,nsa)*f_r2)*drrp
                     elseif(nr.eq.nrmax)then
                        f_r1 = ( (1.d0-wrl)*fnsp(nth,np,nr,nsa) + wrl*fnsp(nth,np,nr-1,nsa) )
                        dfdr_r1 = ( fnsp(nth,np,nr,nsa)-fnsp(nth,np,nr-1,nsa) ) / delr

                        f_r2 = ( (1.d0-wrh)*fs2(nth,np,nsa) + wrh*fnsp(nth,np,nr,nsa) )
                        dfdr_r2 = ( fs2(nth,np,nsa)-fnsp(nth,np,nr,nsa) ) / delr

                        dfdt_r1 = ( drr(nth,np,nr,nsa)*dfdr_r1 - frr(nth,np,nr,nsa)*f_r1 )*drrm
                        dfdt_r2 = ( drr(nth,np,nr+1,nsa)*dfdr_r2 - frr(nth,np,nr+1,nsa)*f_r2)*drrp
                     end if

                     dint_dr = ( dfdt_r2 - dfdt_r1 )/delr/rl*rfsadg(nr)

                     rsumn_dr=rsumn_dr +             dint_dr*volp(nth,np,ns)
                     if(modelr.eq.1)then
                        rsum_dr=rsum_dr+(rgama-1.d0)/theta0(ns)*dint_dr*volp(nth,np,ns)
                     elseif(modelr.eq.0)then
                        rsum_dr=rsum_dr + 0.5d0*pg(np,ns)**2*dint_dr*volp(nth,np,ns)
                     end if
!
                     if(nr.eq.nrmax)then
                        rsumn_drs=rsumn_drs +      dfdt_r2*volp(nth,np,ns)/delr *rfsadg(nr)/rl
                        if(modelr.eq.1)then
                           rsum_drs=rsum_drs+(rgama-1.d0)/theta0(ns)*dfdt_r2*volp(nth,np,ns)/delr *rfsadg(nr)
                        else
                           rsum_drs=rsum_drs+0.5d0*pg(np,ns)**2*dfdt_r2*volp(nth,np,ns)/delr *rfsadg(nr)
                        end if
                     end if
                  end do
               end do
               !     reduce rsum
               call p_theta_integration(rsumn_dr)
               call p_theta_integration(rsum_dr)

               rndrl(nr,nsa)=rnfp0(nsa)*rsumn_dr/(delr*rm(nr))!*rfsadg(nr)
               rpdrl(nr,nsa)=rnfp0(nsa)*1.d20*ptfp0(nsa)**2/amfp(nsa)*1.d-6 &
                    *rsum_dr/(delr*rm(nr))!*rfsadg(nr)
               if(nr.eq.nrmax)then
                  call p_theta_integration(rsumn_drs)
                  call p_theta_integration(rsum_drs)
                  rndrs(nsa) = rnfp0(nsa)*rsumn_drs/(delr*rm(nr))!* 2.d0
                  rpdrs(nsa) = rnfp0(nsa)*1.d20*ptfp0(nsa)**2/amfp(nsa)*1.d-6 &
                       *rsum_drs/(delr*rm(nr))
               end if
! --------- end of radial transport
            enddo ! nsa
         enddo ! nr
         call mtx_reset_communicator
      end if ! modeld

      end subroutine power_from_radial_transport
!------------------------------------------------------------------------------
subroutine total_particle_source
  use fpcomm
  use libmpi
  use fpmpi

  implicit none

  real(rkind)::psum
  integer::nsa,nr,np,nth

  tpsl=0.d0
  call mtx_set_communicator(comm_np)

  do nsa=nsastart,nsaend
    do nr=nrstart,nrend
      psum=0.d0
      do np=npstart,npend
        do nth=1,nthmax
          psum=psum+spp(nth,np,nr,nsa)*volp(nth,np,nsa)*rlamda(nth,nr)
        end do
      end do
      call p_theta_integration(psum)
      tpsl(nr,nsa)=psum*rnfp0(nsa)
    end do
  end do
  call mtx_reset_communicator
end subroutine total_particle_source
!==============================================================


      end module fpsave
