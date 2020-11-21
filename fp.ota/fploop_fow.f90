!     $Id: fploop.f90,v 1.40 2013/02/08 07:36:24 nuga Exp $

! *****************
!     MAIN LOOP
! *****************

      MODULE fowloop

      USE fpcomm
      use fpexec
      use fpsave
      use libmpi
      use fpmpi
      use fpdisrupt
      use eg_read
      use fpoutdata

      contains

!-----------------------------

      subroutine fow_loop

      use fpcoef
      use fpsub
      use libmtx
      use plprof
      use fpmpi
      use fpprep, only: coulomb_log 
      use fpnfrr
      use fpcaltp
      implicit none
      real(kind8):: deps,ip_all_fp,deps_e2

      integer:: nt, nr, np, nth, nsa, ns, ierr, nsb
      real(4):: gut_exe1, gut_exe2, gut_coef1, gut_coef2, gut_coef3
      real(4):: gut_loop1, gut_loop2, gut1, gut2, gut_conv3
      real(4):: sum_gut_ex, sum_gut_coef, gut_1step, sum_gut_conv
      logical:: flag

!     +++++ time loop +++++

    do nt=1,ntmax

      n_impl=0
      deps=1.d0
      do nsa=nsastart,nsaend
        ns=ns_nsa(nsa)
        do nr=nrstartw,nrendwm ! local
          do np=npstartw,npendwm
            do nth=1,nthmax
              fnsm(nth,np,nr,nsa)=fnsp(nth,np,nr,nsa) ! minus step: invariant during n_impl 
            end do
          end do
        end do
      end do
 
      sum_gut_ex = 0.0
      sum_gut_conv=0.0
      sum_gut_coef= 0.0
      gut_1step= 0.0
      call gutime(gut_loop1)

      call gutime(gut_coef3)

         do while(n_impl.le.lmaxfp) ! start do while

            call gutime(gut_exe1)
            call solve_matrix_update_fns0(ierr)
            call gutime(gut_exe2)
            sum_gut_ex = sum_gut_ex + (gut_exe2-gut_exe1)

            call implicit_convergence_update_fnsp(nt,deps)
            call gutime(gut_conv3)
            sum_gut_conv = sum_gut_conv + (gut_conv3-gut_exe2)

            deps_e2=0.d0

            if(nrank.eq.0.and.deps.le.epsfp.and.deps_e2.le.epsfp)then
               n_impl=1+lmaxfp ! exit dowhile
            endif
            call mtx_broadcast1_integer(n_impl)

!-- updating diffusion coef
            call gutime(gut_coef1)
            call update_rn_rt(fnsp) ! update rn_temp, rt_temp for coulomb ln, fpcalcnr, and fpmxwl_exp
            if(model_lnl.eq.0) call coulomb_log ! update coulomb log

            call fusion_source_init
!           update fnsb (fnsb is required by nl collsion and nf reaction)
            flag=.false.
            do nsb=1,nsbmax
               ns=ns_nsb(nsb)
               if(modelc(ns).ge.4) flag=.true.
            end do
            if(flag.or.models.ge.2)then
               call mtx_set_communicator(comm_nsa)
               call update_fnsb_maxwell
               call update_fnsb
               call mtx_reset_communicator
            end if
!           end of update fnsb

! if model_disrupt=1, fp_coef is already called in top_of_time_loop_disrupt
            if(model_disrupt.eq.0)then 
               if (mod(nt,ntstep_coef).eq.0) call fp_coef(nt)
            end if
!           sum up sppf
            if(models.ne.0)then
               call mtx_set_communicator(comm_nsa) 
               call source_allreduce(sppf)
               call mtx_reset_communicator
            end if
!           end of sum up sppf
            call gutime(gut_coef2)
            sum_gut_coef = sum_gut_coef + gut_coef2 - gut_coef1
!-- end of updating diffusion coef

            write(6,'(a,1p3e12.4)') '---deps,deps2=', &
                 deps,deps_e2,epsfp


         end do ! end of dowhile
         sum_gut_coef= sum_gut_coef + gut_coef3 - gut_loop1

         call mtx_reset_communicator
         timefp=timefp+delt

!        bulk f is replaced by maxwellian
         if(model_bulk_const.eq.1)then
            if(model_ex_read_tn.eq.0)then
               call bulk_const1_for_non_exp
            else
               call bulk_const1_for_exp
            end if
         end if

         if(model_disrupt.eq.1)then
            call file_output_disrupt(nt,ip_all_fp) 
         end if

         if(models.ge.2)then
            call allreduce_nf_rate
            call prof_of_nf_reaction_rate(1)
         end if

         do nsa=nsastart,nsaend
            do nr=nrstart,nrend
               do np=npstartw,npendwm
                  do nth=1,nthmax
                     f(nth,np,nr)=fnsp(nth,np,nr,nsa) ! used at fpweight only!
                  end do
               end do
            end do
            call fpweight(nsa,ierr) ! set fpweight for fpsave
         end do

         call gutime(gut_loop2)
         gut_1step = gut_loop2-gut_loop1

         if (mod(nt,ntg1step).eq.0) then
         if(nrank.eq.0) &
              write(6,'(a,e12.4, a,e12.4, a,e12.4, a,e12.4, a,e12.4)') &
              " gut_exec= ", sum_gut_ex,   " gut_conv= ",sum_gut_conv, &
              " gut_coef= ", sum_gut_coef, " gut_1step= ", gut_1step, &
              " exec_ratio = ", sum_gut_ex/gut_1step
         end if
!     +++++ calculate and save global data +++++

         call gutime(gut1)

         isave=0
         call fpssub
         if (mod(nt,ntg1step).eq.0) then
            if(nrank.eq.0) then
               call fpsglb
               call fpwrtglb
               call fp_caltp
            endif
         endif
         if (mod(nt,ntg2step).eq.0) then
            if(nrank.eq.0) then
               call fpsprf
               call fpwrtprf
            endif
         endif
!         call mtx_broadcast_real8(rt_t,nrmax*nsamax)
         call mtx_broadcast_real8(rns,nrmax*nsamax)
         call mtx_broadcast1_integer(ntg1)
         call mtx_broadcast1_integer(ntg2)
         call gutime(gut2)
         if (mod(nt,ntg1step).eq.0) then
            if(nrank.eq.0) write(*,'(a,e14.6)') "--------save_time=",gut2-gut1
         end if

         do nsa=1,nsamax
            do nr=1,nrmax
               if(rns(nr,nsa).lt.0)then
                  ierr_g = ierr_g + 1
                  if(nrank.eq.0)then
                     write(*,'(a,3i5,e14.6)') "negative dens. at nr= ", nr, nsa, ierr_g, rns(nr,nsa)
                  end if
               end if
            end do
         end do

         if(ierr.ne.0)then
            call mtx_abort(ierr)
         end if
         if(ierr_g.ne.0)then
            call mtx_abort(ierr_g)
         end if
         call mtx_reset_communicator 
         if(nrank.eq.0)then
            if(output_txt_beam_dens.eq.1) call number_of_nonthermal_ions
            if(output_txt_heat_prof.eq.1) call out_txt_heat_prof
         end if

      enddo ! end of nt loop

!     +++++ end of time loop +++++

      call gutime(gut1)
      call update_fns
      call gutime(gut2)
!      if(model_disrupt.eq.1) call fluxs_pth
      if(nrank.eq.0) write(6,'(a,e14.6)') "---------time update fns =",gut2-gut1

!  txt format output
      if(nrank.eq.0)then
         if(output_txt_f1.eq.1) call out_txt_f1
         if(output_txt_beam_width.eq.1) call out_txt_beam_width
         if(output_txt_delta_f.eq.1) call out_txt_fns_del
      end if

      return
      end subroutine fow_loop
!------------------------------------
!*************************************************************
!     included in do while
      subroutine solve_matrix_update_fns0(ierr)

      implicit none
      integer:: nsa, nr, np, nth, ns, its
      integer,intent(out):: ierr
      double precision,dimension(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend)::&
           send

      n_impl=n_impl+1
      do nsa=nsastart,nsaend 
         ns=ns_nsa(nsa)
         do nr=nrstart,nrend
            do np=npstartw,npendwm
               do nth=1,nthmax
                  f(nth,np,nr)=fnsp(nth,np,nr,nsa) ! used at fpweight only!
               end do
            end do
         end do
                  
         if(model_connor_fp.eq.1.or.model_disrupt.eq.0)then ! connor model doesn't require f evolution
            call fp_exec(nsa,ierr,its) ! f1 and fns0 are changed
         end if
         ierr=0
         nt_init=1
      end do


      end subroutine solve_matrix_update_fns0
!------------------------------------
      subroutine implicit_convergence_update_fnsp(nt,deps)

      use libmtx
      use fpmpi
      implicit none
      integer:: nsa, nr, np, nth, iloc1, nsw
      integer,intent(in):: nt
      integer,dimension(nsastart:nsaend):: ilocl
      integer,dimension(nsamax):: iloc
      real(kind8),intent(out):: deps
      real(kind8),dimension(nsamax)::rsumf,rsumf0,rsum_ss
      real(kind8),dimension(nsastart:nsaend):: deps_maxvl, depsv
      real(kind8),dimension(nsamax):: deps_maxv
      real(kind8):: rsumf_, rsumf0_, deps_max, deps1
      character:: fmt*40

      nsw = nsaend-nsastart+1      
      do nsa=nsastart,nsaend 
         rsumf(nsa)=0.d0
         rsumf0(nsa)=0.d0
         rsum_ss(nsa)=0.d0
         do nr=nrstart,nrend
            do np=npstart,npend
               do nth=1,nthmax
!                  if(nrank.eq.0) write(6,'(a,4i5,1p4e12.4)') &
!                       'nsa,nr,np,nth,fnsp,fns0,diff,rsumf=', &
!                       nsa,nr,np,nth,fnsp(nth,np,nr,nsa),fns0(nth,np,nr,nsa),&
!                       fnsp(nth,np,nr,nsa)-fns0(nth,np,nr,nsa),rsumf(nsa)
                  rsumf(nsa)=(fnsp(nth,np,nr,nsa)-fns0(nth,np,nr,nsa) )**2 &
                       + rsumf(nsa)
                  rsumf0(nsa)=(fnsp(nth,np,nr,nsa))**2 + rsumf0(nsa)
               enddo
            enddo
         enddo
         do nr=nrstartw,nrendwm
            do np=npstartw,npendwm
               do nth=1,nthmax
                  fnsp(nth,np,nr,nsa)=fns0(nth,np,nr,nsa)
               enddo
            enddo
         enddo
         if(modeld_boundary.eq.1.and.nrend.eq.nrmax)then
            call update_radial_f_boundary(nsa)
         end if
      enddo ! end of nsa
!!---------- convergence criterion
      call mtx_set_communicator(comm_np) 
      do nsa=nsastart,nsaend
         rsumf_=rsumf(nsa)
         rsumf0_=rsumf0(nsa)
         call p_theta_integration(rsumf_) 
         call p_theta_integration(rsumf0_) 
         rsumf(nsa)=rsumf_
         rsumf0(nsa)=rsumf0_
      end do

      deps=0.d0
      do nsa=nsastart,nsaend
         depsv(nsa) = rsumf(nsa)/rsumf0(nsa)
         deps1 = depsv(nsa)
         deps=max(deps,deps1)
      end do
      deps_max=0.d0
      call mtx_set_communicator(comm_nsa) 
      call mtx_reduce1_real8(deps,1,deps_max,iloc1) ! max deps among nsa
      deps = deps_max
      call mtx_set_communicator(comm_nr) 
      call mtx_reduce1_real8(deps,1,deps_max,iloc1) ! max deps among nr
      call mtx_reset_communicator
      deps = deps_max

      call mtx_set_communicator(comm_nr) !3d
      call mtx_allreduce_real8(depsv,nsw,4,deps_maxvl,ilocl) ! the peak depsv for each nsa

      call mtx_set_communicator(comm_nsa) !3d
      call mtx_gather_real8(deps_maxvl,nsw,deps_maxv) 
      call mtx_gather_integer(ilocl,nsw,iloc) 
      call mtx_reset_communicator

      if (mod(nt,ntg1step).eq.0) then
         if(nrank.eq.0) then
            write(fmt,'(a16,i1,a6,i1,a3)') &
                 '(a,1pe12.4,i4,1p',nsamax,'e12.4,',nsamax,'i4)'
            write(6,fmt) 'deps',&
                 deps,n_impl,(deps_maxv(nsa),nsa=1,nsamax) &
                 ,(iloc(nsa),nsa=1,nsamax)
         endif
      end if
      
      end subroutine implicit_convergence_update_fnsp
!------------------------------------
!*************************************************************
      subroutine bulk_const1_for_non_exp ! bulk is const

      use plprof
      implicit none
      integer:: nth, np, nr, nsa, ns
      real(8),dimension(nrmax,nsmax):: tempt, tempn
      type(pl_plf_type),dimension(nsmax):: plf
      real(8):: rhon, fl

!     bulk f is replaced by initial maxwellian
      call define_bulk_np
      do ns=1, nsmax
         do nr=1, nrmax
            tempn(nr,ns)=rn_temp(nr,ns) ! rns
            tempt(nr,ns)=rt_temp(nr,ns) ! rt_bulk
         end do
      end do
      
      do nsa=nsastart, nsaend
         ns=ns_nsa(nsa)
         do nr=nrstartw, nrendwm
            rhon=rm(nr)
            call pl_prof(rhon,plf) ! bulk values are fixed to initial values
            rn_temp(nr,ns)=plf(ns)%rn
            rt_temp(nr,ns)=(plf(ns)%rtpr+2.d0*plf(ns)%rtpp)/3.d0
            if(model_delta_f(ns).eq.0)then
               do np=npstartw, npendwm
                  if(np.le.np_bulk(nr,nsa))then
                     do nth=1, nthmax
                        fl=fpmxwl_exp(pm(np,ns),nr,ns)
                        fns0(nth,np,nr,nsa)=fl
                        fnsp(nth,np,nr,nsa)=fl
                     end do
                  end if
               end do
            else ! delta_f mode: update f_m and delta_f
               if(np.le.np_bulk(nr,nsa))then! eliminate delta f in bulk region
                  do np=npstartw, npendwm
                     do nth=1, nthmax
                        fnsp_del(nth,np,nr,nsa)=0.d0
                     end do
                  end do
               end if
               do nth=1, nthmax
                  fl=fpmxwl_exp(pm(np,ns),nr,ns)
                  fnsp_mxwl(nth,np,nr,nsa)=fl
               end do
               do nth=1, nthmax
                  fns0(nth,np,nr,nsa)=fnsp_del(nth,np,nr,nsa)+fnsp_mxwl(nth,np,nr,nsa)
                  fnsp(nth,np,nr,nsa)=fnsp_del(nth,np,nr,nsa)+fnsp_mxwl(nth,np,nr,nsa)
               end do
            end if
         end do
      end do
      
      do ns=1, nsmax
         do nr=1, nrmax
            rn_temp(nr,ns)=tempn(nr,ns)
            rt_temp(nr,ns)=tempt(nr,ns)
         end do
      end do

      end subroutine bulk_const1_for_non_exp
!-------------------------------------------------------------
      subroutine bulk_const1_for_exp ! bulk is updated by exp

      implicit none
      integer:: nth, np, nr, nsa, ns
      real(8):: fl

!     bulk f is replaced by maxwellian
      call define_bulk_np

      do nsa=nsastart, nsaend
         ns=ns_nsa(nsa)
         if(model_delta_f(ns).eq.0)then
            do nr=nrstart, nrend
               do np=npstartw, npendwm
                  if(np.le.np_bulk(nr,nsa))then
                     do nth=1, nthmax
                        fl=fpmxwl_exp(pm(np,ns),nr,ns)
                        fnsp(nth,np,nr,nsa)=fl
                     end do
                  end if
               end do
            end do
         elseif(model_delta_f(ns).eq.1)then
            do nr=nrstart, nrend
               do np=npstartw, npendwm
                  if(np.le.np_bulk(nr,nsa))then
                     do nth=1, nthmax
                        fnsp_del(nth,np,nr,nsa)=0.d0
                     end do
                  end if
                  do nth=1, nthmax
                     fnsp_mxwl(nth,np,nr,nsa)=fpmxwl_exp(pm(np,ns),nr,ns)
                  end do
                  do nth=1, nthmax
                     fns0(nth,np,nr,nsa)=fnsp_del(nth,np,nr,nsa)+fnsp_mxwl(nth,np,nr,nsa)
                     fnsp(nth,np,nr,nsa)=fnsp_del(nth,np,nr,nsa)+fnsp_mxwl(nth,np,nr,nsa)  
                  end do
               end do
            end do
         end if
      end do

      end subroutine bulk_const1_for_exp
!*************************************************************
!-------------------------------------------------------------
      subroutine update_radial_f_boundary(nsa)

      implicit none
      integer,intent(in):: nsa
      integer:: nth, np

      do np=npstartw, npendwm
         do nth=1, nthmax
            fs2(nth,np,nsa) = 2.d0*fs1(nth,np,nsa) - fnsp(nth,np,nrmax,nsa)
         end do
      end do

      end subroutine update_radial_f_boundary
!------------------------------------------------------------
! ************************************
!     prediction of electric field
! ************************************

!      subroutine fpnewe
!
!      implicit none
!      integer:: nr
!
!      do nr=2,nrmax
!         bp(nr)=bp(nr)+(e1(nr)-e1(nr-1))*delt/(ra*delr)
!      enddo
!
!      do nr=nrstart,nrend
!         rj2(nr)=(rg(nr+1)*bp(nr+1)-rg(nr)*bp(nr)) &
!                 /(rmu0*rm(nr)*delr*ra)
!      enddo
!
!      do nr=nrstart,nrend
!         e2(nr)=rj2(nr)*e1(nr)/rj1(nr)
!      enddo
!
!      return
!      end subroutine fpnewe

!------------------------------------
      end module fowloop
