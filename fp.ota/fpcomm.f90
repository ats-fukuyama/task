module fpcomm_parm

      use bpsd_kinds
      use bpsd_constants
      use commpi
      use plcomm
      implicit none

      public

!     --- input parameters ---

      integer,parameter:: kind8=rkind
      integer,parameter:: nbeamm=20
      real(rkind),parameter:: rkev=aee*1.d3

      integer:: nsamax,nsbmax,ns_nsa(nsm),ns_nsb(nsm)
      integer:: nsa_f1,nth_f1,nr_f1
      integer:: lmax_wr,ncmin(nsm),ncmax(nsm)
      integer:: nbeammax,nsspb(nbeamm),nsspf
      integer:: npmax,nthmax,nrmax,navmax,np2max
      integer:: ntmax,ntstep_coef,ntstep_coll
      integer:: ntg1step,ntg1min,ntg1max
      integer:: ntg2step,ntg2min,ntg2max
      integer:: modele,modela,modelc(nsm),modelr,modeld,models,modelw(nsm)
      integer:: model_delta_f(nsm)
      integer:: model_ne_d,modeld_rdep,modeld_pdep,modeld_edge,modeld_pinch
      integer:: modeld_boundary,modeld_cdbm
      integer:: model_loss,model_synch,model_nbi,model_wave
      integer:: imtx,model_ksp,model_pc,lmaxfp,lmaxe
      integer:: ngline,ngraph,llmax,llmax_nf,idbgfp
      integer:: model_disrupt,model_connor_fp,model_bs,model_jfp,model_lnl
      integer:: model_re_pmax,modeld_n_re,model_impurity,model_sink,n_impu
      integer:: model_ex_read_tn,model_ex_read_dh_ratio
      integer:: model_bulk_const,model_cx_loss
      integer:: n_partition_r,n_partition_s,n_partition_p
      integer:: output_txt_delta_f, output_txt_f1, output_txt_beam_width, output_txt_heat_prof
      integer:: output_txt_beam_dens

      real(rkind):: pmax(nsm),pmax_bb(nsm),emax(nsm)
      real(rkind):: r1,delr1,rmin,rmax
      real(rkind):: e0,zeff
      real(rkind):: pabs_ec,pabs_lh,pabs_fw,pabs_wr,pabs_wm
      real(rkind):: fact_wr,fact_wm,delnpr_wr,delnpr_wm
      real(rkind):: rf_wm,eps_wr,dely_wr,dely_wm,y0_wm
      real(rkind):: dec,pec1,pec2,pec3,pec4,rfec,delyec
      real(rkind):: dlh,plh1,plh2,rlh
      real(rkind):: dfw,pfw1,pfw2,rfw
      complex(rkind):: cewr,cewth,cewph
      real(rkind):: rkwr,rkwth,rkwph
      real(rkind),dimension(nbeamm):: spbtot,spbr0,spbrw,spbeng,spbang,spbpang
      real(rkind):: spftot,spfr0,spfrw,spfeng
      real(rkind):: drr0,drrs,factor_cdbm,drr_edge,rho_edge,factor_drr_edge
      real(rkind):: factor_pinch,deltab_b
      real(rkind),dimension(nsm):: tloss
      real(rkind):: delt,rimpl,epsfp,epsm,epse,epsde,h0de
      real(rkind):: pgmax,rgmax,rgmin
      real(rkind):: t0_quench,tau_quench,tau_mgi
      real(rkind):: time_quench_start,rjprof1,rjprof2
      real(rkind):: v_re,target_zeff,spitot,fact_bulk
      real(rkind):: rn_neu0, rn_neus ! temporal
      real(rkind):: ni_ratio(nsm)
      real(rkind):: pdep_exp

!     for read experiment data
      character(len=80):: eg_name_tms, eg_name_cx, eg_name_ha3
      real(rkind),dimension(:,:),pointer:: read_tms_double, read_cx_double !!    containers of profile data
      integer,dimension(:,:),pointer:: read_tms_int
      real(rkind),dimension(:),pointer:: cte_fit, cti_fit !!    fitting param
      real(rkind),dimension(:),pointer:: cne_fit !!    fitting param
      integer:: nend_tms, nend_cx !!    # time step
      real(rkind),dimension(:),pointer:: rne_exp, rte_exp, rti_exp !! container of exp. profile data at timefp
      real(rkind):: time_exp_offset, rne_exp_edge, rte_exp_edge, rti_exp_edge

!     for read fit3d result
      character(len=80):: sv_file_name_h, sv_file_name_d
      double precision,dimension(:),pointer:: time_grid_fit_h, time_grid_fit_d
      integer:: ntmax_fit_h, ntmax_fit_d
      integer,dimension(:,:),pointer:: i_fit_h, i_fit_d
      double precision,dimension(:),pointer:: d_fit_h, d_fit_d
      integer,dimension(:,:),pointer:: number_of_lines_fit_h, number_of_lines_fit_d

end module fpcomm_parm

module fpcomm
!
      use fpcomm_parm
      use libmpi,only: mtx_mpi_type
      implicit none

      public

!      --- internal variables ---

!         +++ mpi and petsc variables +++
      type(mtx_mpi_type):: comm_nr,comm_nsa,comm_np,&
           comm_nrnp,comm_nsanp,comm_nsanr
      integer:: imtxsize,imtxwidth,imtxstart,imtxend
      integer:: ntg1m,ntg2m
      integer:: nrstart,nrend,nmstart,nmend
      integer:: nsastart,nsaend,npstart,npend
      integer:: npendwm,npendwg,npstartw,nrstartw,nrendwm,nrendwg
      integer:: modeld_temp
      integer,dimension(:),pointer:: mtxlen,mtxpos
      integer,dimension(:),pointer:: savlen
      integer,dimension(:,:),pointer:: savpos
      integer,dimension(:,:),pointer:: rank_partition_data

      integer::isave
      integer,dimension(nsm):: nsb_nsa,nsa_nsb
      real(rkind):: delr, delth
      real(rkind):: timefp
      real(rkind),dimension(:),pointer :: delp
      real(rkind),dimension(:),pointer :: &
           rnfp0,rnfps,rtfp0,rtfps,amfp,aefp,ptfp0,vtfp0, &
           aefd,amfd,ptfd0,vtfd0,theta0,rnfd0,rtfd0,rtfds, &
           rn0_mgi
      integer:: ntg1,ntg2
      real(rkind):: tvolr
      real(rkind):: px
      integer:: nrx,nthx
      integer,dimension(6):: nsa1_nf,nsa2_nf,nsb1_nf,nsb2_nf
      real(rkind),dimension(6):: eng1_nf,eng2_nf

      real(rkind),dimension(:,:,:),pointer :: & ! (nthm,npm,nrm)
           f,f1
      integer,dimension(:),pointer :: & ! (nrm)
           itl,itu, itl_judge, itlg_judge
      integer,dimension(:),pointer :: & ! (nrm)
           itlg,itug,itlg_rg,itug_rg
      real(rkind),dimension(:),pointer :: & ! (nrm)
           volr
      real(rkind),dimension(:,:),pointer :: & ! (nrm,nsam)
           rlamdag,etamg,etamg_rg,etagg_rg,rlamdag_rg
      real(rkind),dimension(:),pointer :: & ! (nrm)
           rg,rm
      real(rkind),dimension(:,:),pointer :: & ! (npm:nsam)
           pg,pm
      real(rkind),dimension(:),pointer :: & ! (npm)
           pgb,pmb
      real(rkind),dimension(:),pointer :: & ! (nthm)
           thg,thm
      real(rkind),dimension(:),pointer :: & ! (nthm)
           re_pitch, rlamda_nrmaxp1, etam_nrmaxp1, etag_nrmaxp1
      real(rkind),dimension(:),pointer :: & ! (nrm)
           bp,qr,rj1,e1,rj2,e2,bpg,bpm,qlm,qlg, &
           ep,em,em_w, &
           rn_disrupt, rn_runaway, rj_ohm, rj_runaway, conduct_sp, &
           sigma_spp, sigma_spm, er_drei, er_crit, rconnor, lnl_gl, rfp_ava, &
           rfpl,rfp,rp_crit,rt_quench,rt_quench_f,previous_rate, rj_bs, &
           previous_rate_p, rn_drei, rj_bsm, rn_runaway_m, r_djdt, &
           previous_rate_g, previous_rate_p_g
      real(rkind),dimension(:),pointer :: & ! (nrm)
           epsrm,epsrg,epsrm2,epsrg2
      real(rkind),dimension(:),pointer :: & ! (nrm)
           epsrmx,epsrgx
      real(rkind),dimension(:,:,:),pointer :: & ! (nthm,npm,nsbm)
           volp
      real(rkind),dimension(:,:),pointer :: & ! (nthm,nrmp)
           etag,etam,rlamda,rlamdc,etam_rg,etag_rg,rlamda_rg,rlamdc_g
      real(rkind),dimension(:),pointer:: & !(nr)
           rfsad,rfsadg, rfsadg_rg, ratio_navmax, a_chi0, line_element

      real(rkind),dimension(:),pointer :: & ! (nthm)
           sing,cosg,sinm,cosm

      real(rkind),dimension(:,:,:,:),pointer :: & ! (nthm,npm,nrm,nsbm)
           fns
      real(rkind),dimension(:,:,:,:),pointer :: & ! (nthm,npm,nrs:nre,nsam)
           fns0,fnsp,fnsm,fnsp_del,fnsp_mxwl
      real(rkind),dimension(:,:,:,:),pointer :: & ! (nthm,npm,nrm,nsbm)
           fns_l
      real(rkind),dimension(:,:,:,:),pointer :: & ! (nthm,npm,nrm,nsbmax)
           fnsb

      real(rkind),dimension(:,:),pointer :: & ! (nrm,nsam)
           rnfp,rtfp,ptfp,vtfp,theta,dkbsr, post_tau_ta,rnfp_g,rtfp_g
      real(rkind),dimension(:,:),pointer :: & ! (nrm,nsbm)
           rnfd,rtfd,ptfd,vtfd, rn_mgi, rn_mgi_g
      real(rkind),dimension(:,:,:),pointer :: & ! (nrm,nsbm,nsam)
           rnuf,rnud,lnlam,post_lnlam_f,post_lnlam
      real(rkind),dimension(:,:,:),pointer :: & ! (nthm,npm,nsam)
           fs0,fs2,fs1
      real(rkind),dimension(:,:,:,:),pointer :: & ! (nthm,npm,nrm,nsam)
           weighp,weight,weighr,weighr_g
      real(rkind),dimension(:,:,:,:),pointer :: & ! (nthm,npm,nrm,nsam)
           dpp,dpt,dtp,dtt,fpp,fth,drr,frr,spp,ppl, &
           fepp,feth,dcpp,dcpt,dctp,dctt,fcpp,fcth, &
           dwpp,dwpt,dwtp,dwtt, &
           dwecpp,dwecpt,dwlhpp,dwlhpt,dwfwpp,dwfwpt, &
           dwwrpp,dwwrpt,dwwmpp,dwwmpt, &
           sppb,sppf,spps,sppi, &
           dwpp_p,dwpt_p,dwtp_p,dwtt_p, &
           dwwrpp_p,dwwrpt_p,dwwmpp_p,dwwmpt_p, &
           dwectp,dwectt,dwwrtp,dwwrtt,dwwmtp,dwwmtt, &
           dcppb,dcptb,fcppb, &
           fspp,fsth,dlpp,flpp,sppl,sppl_cx
      real(rkind),dimension(:,:,:),pointer :: sppd
      real(rkind),dimension(:,:,:,:,:),pointer :: &
           dcpp2,dcpt2,dctp2,dctt2,fcpp2,fcth2  !(nthm,npm,nrm,nsam,nsbm)
      real(rkind),dimension(:,:,:,:,:),pointer :: &
           dcpp2b,dcpt2b,fcpp2b  !(nthm,npm,nrm,nsam,nsbm)
      real(rkind),dimension(:,:),pointer :: & ! (nrm,nsam)
           rnsl,rjsl,rwsl,rpcsl,rpwsl,rpesl,rlhsl,rfwsl,recsl,rws123l, &
           rspbl,rspfl,rspsl,rspll,rpdr,rndr, rtl_bulk, rt_bulk, &
           rwrsl,rwmsl,tpsl, &
           rdidtl, rjsrl, rpssl, rplsl, rspsl_cx
      real(rkind),dimension(:,:),pointer :: & ! (nrm,nsam)
           rnsl_delf, rwsl_para, rwsl_perp, rnsl_prev,rwsl_prev
      real(rkind),dimension(:,:,:),pointer :: & ! (npm,nrm,nsam)
           rp_bulk,rpl_bulk
      real(rkind),dimension(:,:,:),pointer :: & ! (nrm,nsam,nsbm)
           rpcs2l, rpcs2l_del

      real(rkind),dimension(:,:),pointer :: & ! (nrm,nsam)
           rpwec_l,rpwlh_l,rpwfw_l,rpwwr_l,rpwwm_l

      real(rkind),dimension(:,:),pointer :: & ! (nrm,nsam)
           rns,rjs,rws,rpcs,rpws,rpes,rlhs,rfws,recs,rws123, &
           rspb,rspf,rsps,rspl,rpdrl,rndrl,rwrs,rwms,rdidt,&
           rjsr,rpss,rpls,rjs_m,rsps_cx, rns_delf,rns_prev,rws_prev,tps
      real(rkind),dimension(:,:),pointer :: & ! (nrm,nsam)
           rns_delf_nsa, rws_delf_para, rws_delf_perp, rws_para, rws_perp
      real(rkind),dimension(:),pointer :: & ! (nsam)
           rns_s2
      real(rkind),dimension(:,:,:),pointer :: & ! (nrm,nsam,nsbm)
           rpcs2, rpcs2_del

      real(rkind),dimension(:),pointer :: & ! (ntg1m)
           ptg,pet,pqt,q_eng
      real(rkind),dimension(:,:),pointer :: & ! (ntg1m,nsam)
           pnt,pwt,ptt,pit,ppct,ppwt,ppet,plht,pfwt,pect,ptt3, &
           pitt,pwtt,pwrt,pwmt,pirt,ppst,pplt
      real(rkind),dimension(:,:),pointer :: & ! (ntg1m,nsam)
           pspt,pspbt,pspft,psplt,pspst
      real(rkind),dimension(:,:),pointer :: & ! (ntg1m,nsam)
           pnt2,pwt2,ptt2,pit2,pwtd,pdr,pndr,ptt_bulk
      real(rkind),dimension(:,:,:),pointer :: & ! (ntg1m,nsam,nsbm)
           ppct2
      real(rkind),dimension(:,:),pointer :: & ! (nrm,nsm)
           ripp
      real(rkind),dimension(:),pointer :: & ! (ntg2m)
           rtg
      real(rkind),dimension(:,:),pointer :: & ! (nrm,ntg2m)
           ret,rqt,rate_runaway
      real(rkind),dimension(:,:,:),pointer :: & ! (nrm,ntg2m,nsam)
           rnt,rwt,rtt,rjt,rpct,rpwt,rpet,rlht,rfwt,rect, &
           rspbt,rspft,rsplt,rspst,rpdrt,rndrt,rtt_bulk,rwrt,rwmt,&
           rate_runaway2,rjrt
      real(rkind),dimension(:,:,:,:),pointer :: & ! (nrm,ntg2m,nsam,nsbm)
           rpct2
      real(rkind),dimension(:,:),pointer:: & !(nrm, nsam)
           rn_temp, rt_temp, rn_read, rt_read
      integer:: nmmax,nlmaxm
      integer,dimension(:,:,:),pointer :: & ! (nthm,npm,nrm)
           nma
      integer,dimension(:),pointer :: & ! (nmm)
           nlmax
      integer,dimension(:,:),pointer :: & ! (nmm,nlm)
           ll
      real(rkind),dimension(:),pointer :: & ! (nmm)
           dl,bm
      real(rkind),dimension(:,:),pointer :: & ! (nmm,nlm)
           al
      real(rkind),dimension(:),pointer :: & ! (nrm*nthm*npm)
           fm, fm_shadow_m, fm_shadow_p!,bmtot

      real(rkind),dimension(:,:,:,:,:),pointer :: &
           sigmav_nf ! (nthmax+1,npmax+1,nthmax+1,npmax+1,6)
      real(rkind),dimension(:,:,:,:,:),pointer :: &
           sigmav_lg ! (0:llmax_nf,npstart:npend,0:llmax_nf,npstart:npend,6)
      real(rkind),dimension(:,:),pointer :: &
           pl_nf ! (0:llmax_nf,nthmax)
      real(rkind),dimension(:,:),pointer :: &
           rate_nf, rate_nf_bb ! (nrstart:nrend,6)
      real(rkind),dimension(:,:,:,:),pointer :: &
           rate_nf_d1, rate_nf_d2 ! (nthmax,npmax,nrstart:nrend,6)
      real(rkind),dimension(:,:),pointer :: & ! (npm:nsam)
           pg2,pm2
      real(rkind),dimension(:),pointer :: & ! (nsam)
           deps_ss, rpdrs, rndrs,tau_ta0,e_drei0,e_crit0,post_tau_ta0_f
      integer:: n_impl
      real,dimension(10):: gut_comm
      real(rkind),dimension(:,:),pointer:: eptr
      real(rkind):: e_edgem, sigp_e, rn_e, rt_e, rlnrl_e
      real(rkind):: pc_runaway, rf_wr
      real(rkind):: zeff_imp, taue0_nb_e, taue0_nb_i, taue0_nb
      integer:: npc_runaway
      integer:: nt_init, n_f1
      integer:: ierr_g
      integer,dimension(:,:),pointer :: & ! (nrm,nsm)
           np_bulk
      real(rkind):: ebeam0, ebeam1
      real(rkind),dimension(:),allocatable::taup,taue
      contains

        subroutine fp_allocate
          implicit none
          integer,save:: nrstart_save=0,nrend_save=0,nrmax_save=0
          integer,save:: nthmax_save=0,npmax_save=0
          integer,save:: nsamax_save=0,nsbmax_save=0
          integer,save:: init=0

          if(init.eq.0) then
             init=1
          else
             if((npmax.eq.npmax_save).and. &
                (nthmax.eq.nthmax_save).and. &
                (nrmax.eq.nrmax_save).and. &
                (nrstart.eq.nrstart_save).and. &
                (nrend.eq.nrend_save).and. &
                (nsamax.eq.nsamax_save).and. &
                (nsbmax.eq.nsbmax_save)) return

             call fp_deallocate
          endif

          allocate( f(nthmax,npstartw:npendwm,nrstart:nrend))
          allocate(f1(nthmax,npstartw:npendwm,nrstart:nrend))

          allocate(rg(nrmax+1),rm(nrmax),volr(nrmax))
          allocate(rlamdag(nthmax,nrmax+1),rlamdag_rg(nthmax,nrmax+1))
          allocate(etamg(nthmax,nrmax+1),etamg_rg(nthmax,nrmax+1))
          allocate(bp(nrmax+1),qr(nrmax))
          allocate(bpg(nrmax+1),bpm(nrmax+1))
          allocate(qlg(nrmax+1),qlm(nrmax+1))
          allocate(rj1(nrmax),e1(nrmax))
          allocate(rj2(nrmax),e2(nrmax))
          allocate(ep(nrstart:nrend), em(nrstart:nrend), em_w(nrstart-1:nrend+1))

          allocate(epsrm(nrmax+1),epsrg(nrmax+1))
          allocate(epsrm2(nrmax+1),epsrg2(nrmax+1))
          allocate(epsrmx(nrmax+1),epsrgx(nrmax+1))
          allocate(itl(nrmax+1),itu(nrmax+1))
          allocate(itlg(nrmax+1),itug(nrmax+1))
          allocate(itlg_rg(nrmax+1),itug_rg(nrmax+1))
          allocate(itl_judge(nrmax+1),itlg_judge(nrmax+1))

          allocate(pg(npmax+1,nsmax),pm(npmax,nsmax))
          allocate(thg(nthmax+1),thm(nthmax))
          allocate(delp(nsmax))
          allocate(pg2(np2max+1,nsmax),pm2(np2max,nsmax))

          allocate(volp(nthmax,npstart:npend,nsmax))
          allocate(etag(nthmax+1,nrstart:nrend+1),etam(nthmax,nrstart:nrend+1))
          allocate(etag_rg(nthmax+1,nrstart:nrend+1),etam_rg(nthmax,nrstart:nrend+1))
          allocate(rlamda(nthmax,nrstart:nrend),rlamdc(nthmax+1,nrstart:nrend))
          allocate(rlamda_nrmaxp1(nthmax),etam_nrmaxp1(nthmax),etag_nrmaxp1(nthmax+1))

          allocate(rfsad(nrstart:nrend),rfsadg(nrmax+1))
          allocate(rfsadg_rg(nrmax+1))
          allocate(ratio_navmax(nrstart:nrend))
          allocate(a_chi0(nrstart:nrend))
          allocate(line_element(nrmax+1))

          allocate(rlamda_rg(nthmax,nrstart:nrendwg),rlamdc_g(nthmax+1,nrstart:nrend))
          allocate(sing(nthmax+1),cosg(nthmax+1))
          allocate(sinm(nthmax),cosm(nthmax))

!          allocate(fns(nthmax,npmax,nrmax,nsbmax))
          allocate(fns(nthmax,npmax,nrmax,nsamax))

          allocate(fns0(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend))
          allocate(fnsp(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend))
          allocate(fnsm(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend))
          allocate(fnsp_del(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend))
          allocate(fnsp_mxwl(nthmax,npstartw:npendwm,nrstartw:nrendwm,nsastart:nsaend))
          allocate(fnsb(nthmax,npstart:npend,nrstart:nrend,nsbmax)) ! backgroud f

          allocate(fs0(nthmax,npstartw:npendwm,nsamax))
          allocate(fs2(nthmax,npstartw:npendwm,nsamax))
          allocate(fs1(nthmax,npstartw:npendwm,nsamax))

          allocate(rnfp0(nsamax),rnfps(nsamax))
          allocate(rtfp0(nsamax),rtfps(nsamax))
          allocate(rtfd0(nsbmax),rtfds(nsbmax))
          allocate(amfp(nsamax),aefp(nsamax))
          allocate(ptfp0(nsamax),vtfp0(nsamax))
          allocate(rnfd0(nsbmax))
          allocate(aefd(nsbmax),amfd(nsbmax))
          allocate(ptfd0(nsbmax),vtfd0(nsbmax))
          allocate(theta0(nsmax))

          allocate(rnfp(nrstart:nrend+1,nsamax),rtfp(nrstart:nrend+1,nsamax))
          allocate(rnfp_g(nrstart:nrendwg,nsamax),rtfp_g(nrstart:nrendwg,nsamax))
          allocate(ptfp(nrstart:nrend+1,nsamax),vtfp(nrstart:nrend+1,nsamax))
          allocate(theta(nrstart:nrend,nsmax),dkbsr(nrstart:nrend,nsamax))
          allocate(weighp(nthmax  ,npstart:npendwg,nrstart:nrend+1,nsamax))
          allocate(weight(nthmax+1,npstartw:npendwm,nrstart:nrend+1,nsamax))
          allocate(weighr(nthmax,npstart:npend,nrstart:nrend+1,nsamax))
          if(modeld.ne.0)then
             allocate(weighr_g(nthmax,npmax,nrmax+1,nsamax))
          end if
          allocate(rnfd(nrstart:nrendwm,nsbmax),rtfd(nrstart:nrendwm,nsbmax))
          allocate(rn0_mgi(nsbmax))
          allocate(rn_mgi(nrstart:nrend,nsbmax))
          allocate(rn_mgi_g(nrmax,nsbmax))
          allocate(ptfd(nrstart:nrend+1,nsbmax),vtfd(nrstart:nrend+1,nsbmax))
          allocate(rnuf(nrmax,nsbmax,nsbmax))
          allocate(rnud(nrmax,nsbmax,nsbmax))
          allocate(lnlam(nrstart:nrend+1,nsbmax,nsbmax))

          allocate(dpp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dpt(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dtp(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          allocate(dtt(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          allocate(fpp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(fth(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          allocate(drr(nthmax  ,npstart:npend,nrstart:nrendwg,nsamax))
          allocate(frr(nthmax  ,npstart:npend,nrstart:nrendwg,nsamax))

          allocate(fepp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(feth(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))

          allocate(dcpp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dcpt(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dctp(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          allocate(dctt(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          allocate(fcpp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(fcth(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))

          allocate(dlpp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(flpp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))

          allocate(fspp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(fsth(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))

          allocate(dcpp2(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsbmax,nsamax))
          allocate(dcpt2(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsbmax,nsamax))
          allocate(dctp2(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsbmax,nsamax))
          allocate(dctt2(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsbmax,nsamax))
          allocate(fcpp2(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsbmax,nsamax))
          allocate(fcth2(nthmax+1,npstartw:npendwg,nrstart:nrendwm,nsbmax,nsamax))

!          allocate(dcppb(4,npmax+1,nrstart:nrend+1,nsamax))
!          allocate(dcptb(4,npmax+1,nrstart:nrend+1,nsamax))
!          allocate(fcppb(4,npmax+1,nrstart:nrend+1,nsamax))

!          allocate(dcpp2b(4,npmax+1,nrstart:nrend+1,nsbmax,nsamax))
!          allocate(dcpt2b(4,npmax+1,nrstart:nrend+1,nsbmax,nsamax))
!          allocate(fcpp2b(4,npmax+1,nrstart:nrend+1,nsbmax,nsamax))

          allocate(dwpp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dwpt(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dwtp(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          allocate(dwtt(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
!          allocate(dwpp(nthmax  ,npmax+1,nrstart:nrend+1,nsamax))
!          allocate(dwpt(nthmax  ,npmax+1,nrstart:nrend+1,nsamax))
!          allocate(dwtp(nthmax+1,npmax  ,nrstart:nrend+1,nsamax))
!          allocate(dwtt(nthmax+1,npmax  ,nrstart:nrend+1,nsamax))
          allocate(dwlhpp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dwlhpt(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dwfwpp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dwfwpt(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dwecpp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dwecpt(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dwectp(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          allocate(dwectt(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          allocate(dwwrpp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dwwrpt(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dwwrtp(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          allocate(dwwrtt(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          allocate(dwwmpp(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dwwmpt(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
          allocate(dwwmtp(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          allocate(dwwmtt(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          if(model_wave.ne.0)then
             allocate(dwpp_p(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
             allocate(dwpt_p(nthmax  ,npstart :npendwg,nrstart:nrendwm,nsamax))
             allocate(dwtp_p(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
             allocate(dwtt_p(nthmax+1,npstartw:npendwm,nrstart:nrendwm,nsamax))
          end if
          if(model_disrupt.ne.0)then
             allocate(er_drei(nrmax),er_crit(nrmax),rconnor(nrmax), lnl_gl(nrmax),rp_crit(nrmax))
             allocate(rn_disrupt(nrmax),rn_runaway(nrmax), rn_drei(nrmax),rn_runaway_m(nrmax))
             allocate(rj_ohm(nrmax),rj_runaway(nrmax),rj_bs(nrmax),r_djdt(nrmax))
             allocate(rj_bsm(nrstart:nrend))
             allocate(previous_rate(nrstart:nrend),previous_rate_p(nrstart:nrend))
             allocate(previous_rate_g(nrmax),previous_rate_p_g(nrmax))
             allocate(sigma_spp(nrstart:nrend),sigma_spm(nrstart:nrend))
             allocate(re_pitch(nthmax))
             allocate(rt_quench(nrmax),rt_quench_f(nrmax))
             allocate(post_lnlam_f(nrstart:nrend+1,nsbmax,nsbmax))
             allocate(post_lnlam(nrstart:nrend+1,nsbmax,nsbmax))
             allocate(post_tau_ta0_f(nsamax))
             allocate(post_tau_ta(nrmax,nsamax))
             allocate(e_drei0(nsamax),e_crit0(nsamax))
          end if
          allocate(conduct_sp(nrmax))

          allocate(spp (nthmax,npstart:npend,nrstart:nrend,nsamax))
          allocate(ppl (nthmax,npstart:npend,nrstart:nrend,nsamax))
          allocate(sppb(nthmax,npstart:npend,nrstart:nrend,nsamax))
          allocate(sppf(nthmax,npstart:npend,nrstart:nrend,nsamax))
          allocate(spps(nthmax,npstart:npend,nrstart:nrend,nsamax))
          allocate(sppi(nthmax,npstart:npend,nrstart:nrend,nsamax))
          allocate(sppd(nthmax,npstart:npend,nsamax))
          allocate(sppl(nthmax,npstart:npend,nrstart:nrend,nsamax))
          allocate(sppl_cx(nthmax,npstart:npend,nrstart:nrend,nsamax))

          allocate(rn_temp(nrmax,nsmax),rt_temp(nrmax,nsmax) )
          if(model_ex_read_tn.ne.0)then
             allocate(rn_read(nrmax,nsmax),rt_read(nrmax,nsmax) )
             allocate(rne_exp(nrmax), rte_exp(nrmax), rti_exp(nrmax))
             allocate(cte_fit(5),cne_fit(6),cti_fit(5))
          end if
          allocate(rnsl(nrstart:nrend,nsamax),rjsl(nrstart:nrend,nsamax))
          allocate(rnsl_prev(nrstart:nrend,nsamax),rjsl(nrstart:nrend,nsamax))
          allocate(rnsl_delf(nrstart:nrend,nsamax))
          allocate(rwsl_para(nrstart:nrend,nsamax),rwsl_perp(nrstart:nrend,nsamax))
          allocate(rfpl(nrstart:nrend),rjsrl(nrstart:nrend,nsamax))
          allocate(rwsl(nrstart:nrend,nsamax),rws123l(nrstart:nrend,nsamax))
          allocate(rspbl(nrstart:nrend,nsamax),rspfl(nrstart:nrend,nsamax))
          allocate(rspsl(nrstart:nrend,nsamax),rspll(nrstart:nrend,nsamax))
          allocate(rspsl_cx(nrstart:nrend,nsamax))
          allocate(rpcsl(nrstart:nrend,nsamax),rpesl(nrstart:nrend,nsamax))
          allocate(rpwsl(nrstart:nrend,nsamax),rlhsl(nrstart:nrend,nsamax))
          allocate(rfwsl(nrstart:nrend,nsamax),recsl(nrstart:nrend,nsamax))
          allocate(rwrsl(nrstart:nrend,nsamax),rwmsl(nrstart:nrend,nsamax))
          allocate(rpcs2l(nrstart:nrend,nsbmax,nsamax))
          allocate(rpcs2l_del(nrstart:nrend,nsbmax,nsamax))
          allocate(rdidtl(nrstart:nrend,nsamax))
          allocate(rdidt(nrmax,nsamax))
          allocate(rpssl(nrstart:nrend,nsamax))
          allocate(rplsl(nrstart:nrend,nsamax))
          allocate(rwsl_prev(nrstart:nrend,nsamax))
          allocate(tpsl(nrstart:nrend,nsamax))


          allocate(rns(nrmax,nsamax),rjs(nrmax,nsamax),rjs_m(nrmax,nsamax))
          allocate(rns_delf_nsa(nrmax,nsamax))
          allocate(rws_delf_para(nrmax,nsamax),rws_delf_perp(nrmax,nsamax))
          allocate(rns_delf(nrmax,nsmax))
          allocate(rns_prev(nrmax,nsamax))
          allocate(rws_para(nrmax,nsmax),rws_perp(nrmax,nsmax))
          allocate(rws_prev(nrmax,nsamax))
          allocate(tps(nrmax,nsamax))

          allocate(rfp(nrmax),rjsr(nrmax,nsamax))
          allocate(rfp_ava(nrmax))
          allocate(rns_s2(nsamax))
          allocate(rws(nrmax,nsamax),rws123(nrmax,nsamax))
          allocate(rspb(nrmax,nsamax),rspf(nrmax,nsamax))
          allocate(rspl(nrmax,nsamax),rsps(nrmax,nsamax))
          allocate(rsps_cx(nrmax,nsamax))
          allocate(rpcs(nrmax,nsamax),rpes(nrmax,nsamax))
          allocate(rpws(nrmax,nsamax),rlhs(nrmax,nsamax))
          allocate(rfws(nrmax,nsamax),recs(nrmax,nsamax))
          allocate(rwrs(nrmax,nsamax),rwms(nrmax,nsamax))
          allocate(rpcs2(nrmax,nsbmax,nsamax))
          allocate(rpcs2_del(nrmax,nsbmax,nsamax))
          allocate(rpss(nrmax,nsamax))
          allocate(rpls(nrmax,nsamax))

          allocate(rpdr(nrmax,nsamax),rndr(nrmax,nsamax))
          allocate(rpdrs(nsamax),rndrs(nsamax))
          allocate(tau_ta0(nsamax))
          allocate(rpdrl(nrstart:nrend,nsamax),rndrl(nrstart:nrend,nsamax))
          allocate(rt_bulk(nrmax,nsamax))
          allocate(rtl_bulk(nrstart:nrend,nsamax))
          allocate(rp_bulk(npmax,nrmax,nsamax))
          allocate(rpl_bulk(npmax,nrstart:nrend,nsastart:nsaend))

          allocate(rpwlh_l(nrstart:nrend,nsamax))
          allocate(rpwfw_l(nrstart:nrend,nsamax))
          allocate(rpwec_l(nrstart:nrend,nsamax))
          allocate(rpwwr_l(nrstart:nrend,nsamax))
          allocate(rpwwm_l(nrstart:nrend,nsamax))

          allocate(ripp(nrmax,nsamax))
!         nlmaxm= 8   ! this is for analysis without bounce average
!         nlmaxm=11   ! this is for analysis without radial transport
          nlmaxm=15   ! this is for analysis with a simple radial transport
          nmmax=nrmax*nthmax*npmax

          allocate(nma(nthmax,npstartw:npendwm,nrstartw:nrendwm))
          allocate(nlmax(nmstart:nmend))
          allocate(ll(nmstart:nmend,nlmaxm))
          allocate(dl(nmstart:nmend),bm(nmstart:nmend))
          allocate(al(nmstart:nmend,nlmaxm))
          allocate(fm(nmstart-nthmax:nmend+nthmax))
          allocate(fm_shadow_m(nmstart-nthmax-nthmax*npmax:nmend+nthmax-nthmax*npmax))
          allocate(fm_shadow_p(nmstart-nthmax+nthmax*npmax:nmend+nthmax+nthmax*npmax))
!          allocate(fm(1+nthmax*(npstartw-1)+nthmax*npmax*(nrstartw-1):npmax*nthmax*nrendwm) )
!          allocate(bmtot(nmmax))

          if(models.eq.2.or.models.eq.3)then
!             if(models.eq.2) allocate(sigmav_nf(nthmax,npmax,nthmax,npmax,6))
             if(models.eq.2) allocate(sigmav_nf(nthmax,npstart:npend,nthmax,npmax,6))
             if(models.eq.3) allocate( &
                  sigmav_lg(0:llmax_nf,npstart:npend,0:llmax_nf,npmax,6), &
                  pl_nf(0:llmax_nf,nthmax))
             allocate(rate_nf(nrstart:nrend,6),rate_nf_bb(nrstart:nrend,6))
!             allocate(rate_nf_d1(nthmax,npmax,nrstart:nrend,6))
             allocate(rate_nf_d2(nthmax,npmax,nrstart:nrend,6))
             allocate(rate_nf_d1(nthmax,npstart:npend,nrstart:nrend,6))
!             allocate(rate_nf_d2(nthmax,npstart:npend,nrstart:nrend,6))
          endif
          allocate(np_bulk(nrmax,nsamax))
          allocate(taup(nsamax),taue(nsamax))
          allocate(deps_ss(nsamax))
          npmax_save=npmax
          nthmax_save=nthmax
          nrmax_save=nrmax
          nrstart_save=nrstart
          nrend_save=nrend
          nsamax_save=nsamax
          nsbmax_save=nsbmax
          return
        end subroutine fp_allocate

        subroutine fp_deallocate
          implicit none

          deallocate(mtxlen,mtxpos,savlen)

          deallocate(f)
          deallocate(f1)

          deallocate(rg,rm,volr)
          deallocate(rlamdag,rlamdag_rg)
          deallocate(etamg,etamg_rg)
          deallocate(bp,qr)
          deallocate(bpg,bpm)
          deallocate(qlg,qlm)
          deallocate(rj1,e1,rj2,e2)
          deallocate(ep,em,em_w)

          deallocate(epsrm,epsrg)
          deallocate(epsrm2,epsrg2)
          deallocate(epsrmx,epsrgx)
          deallocate(itl,itu)
          deallocate(itlg,itug)
          deallocate(itlg_rg,itug_rg)
          deallocate(itl_judge,itlg_judge)

          deallocate(pg,pm)
          deallocate(thg,thm)
          deallocate(delp)
          deallocate(pg2,pm2)

          deallocate(volp)
          deallocate(etag,etam)
          deallocate(etag_rg,etam_rg)
          deallocate(rlamda,rlamdc,rlamda_nrmaxp1)
          deallocate(etam_nrmaxp1, etag_nrmaxp1)
          deallocate(rfsad,rfsadg)
          deallocate(rfsadg_rg)
          deallocate(ratio_navmax, a_chi0)
          deallocate(line_element)
          deallocate(rlamda_rg,rlamdc_g)
          deallocate(sing,cosg,sinm,cosm)

          deallocate(fns)
          deallocate(fns0)
          deallocate(fnsp)
          deallocate(fnsm)
          deallocate(fnsp_del,fnsp_mxwl)
          deallocate(fnsb)
          deallocate(fs0,fs2,fs1)

          deallocate(rnfp0,rnfps)
          deallocate(rtfp0,rtfps,rtfd0,rtfds)
          deallocate(amfp,aefp)
          deallocate(ptfp0,vtfp0)
          deallocate(rnfd0)
          deallocate(aefd,amfd)
          deallocate(ptfd0,vtfd0)
          deallocate(theta0)

          deallocate(rnfp,rtfp)
          deallocate(rnfp_g,rtfp_g)
          deallocate(ptfp,vtfp)
          deallocate(theta,dkbsr)
          deallocate(weighp,weight)
          deallocate(weighr)
          if(modeld.ne.0)then
             deallocate(weighr)
          end if

          deallocate(rnfd,rtfd,ptfd,vtfd)
          deallocate(rn_mgi, rn_mgi_g)
          deallocate(rn0_mgi)
          deallocate(rnuf,rnud,lnlam)
          deallocate(dpp,dpt)
          deallocate(dtp,dtt)
          deallocate(fpp,fth)
          deallocate(drr,frr)
          deallocate(spp,ppl)
          deallocate(fepp,feth)

          deallocate(fspp,fsth)

          deallocate(dlpp,flpp)

          deallocate(dcpp,dcpt)
          deallocate(dctp,dctt)
          deallocate(fcpp,fcth)

          if(model_wave.ne.0)then
             deallocate(dwpp,dwpt)
             deallocate(dwtp,dwtt)
             deallocate(dwpp_p,dwpt_p)
             deallocate(dwtp_p,dwtt_p)

             deallocate(dwlhpp,dwlhpt)
             deallocate(dwfwpp,dwfwpt)
             deallocate(dwecpp,dwecpt,dwectp,dwectt)
             deallocate(dwwrpp,dwwrpt,dwwrtp,dwwrtt)
             deallocate(dwwmpp,dwwmpt,dwwmtp,dwwmtt)
          end if
          if(model_disrupt.ne.0)then
             deallocate(er_drei, er_crit,rp_crit)
             deallocate(previous_rate, previous_rate_p)
             deallocate(previous_rate_g, previous_rate_p_g)
             deallocate(rt_quench, rt_quench_f, conduct_sp)
             deallocate(re_pitch)
             deallocate(rn_disrupt, rn_runaway, rj_ohm, rj_runaway, rn_drei,rn_runaway_m)
             deallocate(rj_bs, rj_bsm, r_djdt)
             deallocate(sigma_spp,sigma_spm)
             deallocate(post_lnlam_f,post_lnlam)
             deallocate(e_drei0,e_crit0)
             deallocate(post_tau_ta0_f)
             deallocate(post_tau_ta)
          end if
          deallocate(tau_ta0)

          deallocate(sppb,sppf,spps,sppd,sppi,sppl)
          deallocate(sppl_cx)


          deallocate(dcpp2,dcpt2)
          deallocate(dctp2,dctt2)
          deallocate(fcpp2,fcth2)

!          deallocate(dcppb,dcptb,fcppb)
!          deallocate(dcpp2b,dcpt2b,fcpp2b)

          deallocate(rn_temp,rt_temp)
          if(model_ex_read_tn.ne.0)then
             deallocate(rn_read, rt_read)
             deallocate(rne_exp, rte_exp, rti_exp)
             deallocate(cte_fit, cne_fit, cti_fit)
          end if
          deallocate(rnsl,rjsl,rwsl,rws123l,rfpl,rjsrl,tpsl,rwsl_prev)
          deallocate(rspbl,rspfl,rspsl,rspll,rpcsl,rpesl,rspsl_cx)
          deallocate(rlhsl,rfwsl,recsl,rwrsl,rwmsl,rpcs2l,rpcs2l_del)
          deallocate(rdidt, rdidtl)
          deallocate(rpssl, rplsl)
          deallocate(rnsl_delf,rwsl_para,rwsl_perp)

          deallocate(rns,rjs,rjs_m,rfp,rjsr,rconnor,rfp_ava,tps,rws_prev)
          deallocate(rns_delf_nsa,rws_delf_para,rws_delf_perp)
          deallocate(rns_delf,rws_para,rws_perp)
          deallocate(rns_s2)
          deallocate(rws,rws123)
          deallocate(rspb,rspf)
          deallocate(rspl,rsps,rsps_cx)
          deallocate(rpcs,rpes)
          deallocate(rpws,rlhs)
          deallocate(rfws,recs,rwrs,rwms,rpcs2,rpcs2_del)
          deallocate(rpss, rpls)

          deallocate(rpdr,rndr,rpdrs,rndrs)
          deallocate(rpdrl,rndrl,rt_bulk,rtl_bulk)
          deallocate(rp_bulk)

          deallocate(rpwlh_l,rpwfw_l,rpwec_l,rpwwr_l,rpwwm_l)
          deallocate(ripp)

          deallocate(nma)
          deallocate(nlmax,ll)
          deallocate(dl,bm)
          deallocate(al)
          deallocate(fm)
          deallocate(fm_shadow_m, fm_shadow_p)
!          deallocate(bmtot)
          deallocate(np_bulk)

          if(models.eq.2.or.models.eq.3)then
             if(models.eq.2) deallocate(sigmav_nf)
             if(models.eq.3) deallocate(sigmav_lg,pl_nf)
             deallocate(rate_nf,rate_nf_d1,rate_nf_d2,rate_nf_bb)
          end if
          deallocate(deps_ss)
          return

        end subroutine fp_deallocate

        subroutine fp_allocate_ntg1
          implicit none
          integer,save:: nsamax_save=0,nsbmax_save=0,ntg1m_save=0
          integer,save:: init=0

          if((nsamax.eq.nsamax_save).and. &
             (nsbmax.eq.nsbmax_save).and. &
             (ntg1m.eq.ntg1min)) return

          if(init.eq.0) then
             init=1
          else
             call fp_deallocate_ntg1
          endif

          ntg1m=ntg1min

          allocate(ptg(ntg1m))
          allocate(pet(ntg1m))
          allocate(pqt(ntg1m))
          allocate(q_eng(ntg1m))
          allocate(pnt(nsamax,ntg1m))
          allocate(pwt(nsamax,ntg1m))
          allocate(ptt(nsamax,ntg1m))
          allocate(pit(nsamax,ntg1m))
          allocate(pirt(nsamax,ntg1m))
          allocate(ppct(nsamax,ntg1m))
          allocate(ppst(nsamax,ntg1m))
          allocate(pplt(nsamax,ntg1m))
          allocate(ppwt(nsamax,ntg1m))
          allocate(ppet(nsamax,ntg1m))
          allocate(plht(nsamax,ntg1m))
          allocate(pfwt(nsamax,ntg1m))
          allocate(pect(nsamax,ntg1m))
          allocate(pwrt(nsamax,ntg1m))
          allocate(pwmt(nsamax,ntg1m))
          allocate(ptt3(nsamax,ntg1m))
          allocate(pitt(nsamax,ntg1m))
          allocate(pwtt(nsamax,ntg1m))
          allocate(pnt2(nsamax,ntg1m))
          allocate(pwt2(nsamax,ntg1m))
          allocate(pwtd(nsamax,ntg1m))
          allocate(ptt2(nsamax,ntg1m))
          allocate(pit2(nsamax,ntg1m))
          allocate(pspt(nsamax,ntg1m))
          allocate(pspbt(nsamax,ntg1m))
          allocate(pspft(nsamax,ntg1m))
          allocate(psplt(nsamax,ntg1m))
          allocate(pspst(nsamax,ntg1m))
          allocate(ppct2(nsbmax,nsamax,ntg1m))
          allocate(pdr(nsamax,ntg1m),pndr(nsamax,ntg1m))
          allocate(ptt_bulk(nsamax,ntg1m))

          nsamax_save=nsamax
          nsbmax_save=nsbmax
          return
        end subroutine fp_allocate_ntg1

        subroutine fp_deallocate_ntg1
          implicit none

          deallocate(ptg)
          deallocate(pet)
          deallocate(pqt)
          deallocate(q_eng)
          deallocate(pnt)
          deallocate(pwt)
          deallocate(ptt)
          deallocate(pit)
          deallocate(pirt)
          deallocate(ppct)
          deallocate(ppst)
          deallocate(pplt)
          deallocate(ppwt)
          deallocate(ppet)
          deallocate(plht)
          deallocate(pfwt)
          deallocate(pect)
          deallocate(pwrt)
          deallocate(pwmt)
          deallocate(ptt3)
          deallocate(pitt)
          deallocate(pwtt)
          deallocate(pnt2)
          deallocate(pwt2)
          deallocate(pwtd)
          deallocate(ptt2)
          deallocate(pit2)
          deallocate(pspt,pspbt,pspft,psplt,pspst)
          deallocate(ppct2)
          deallocate(pdr,pndr)
          deallocate(ptt_bulk)
          deallocate(taup,taue)

        end subroutine fp_deallocate_ntg1

        subroutine fp_adjust_ntg1
          implicit none
          real(rkind),dimension(:),pointer:: tempa
          real(rkind),dimension(:,:),pointer:: tempb
          real(rkind),dimension(:,:,:),pointer:: tempc
          integer:: ntg,nsa,nsb,ntg1m_new

          if(ntg1.gt.ntg1m) then
             write(6,*) '# fp_adjust_ntg1:',&
                        ntg1,size(ptg),ntg1m,ntg1min,ntg1max
             if(ntg1m.ge.ntg1max) then
                ntg1=(ntg1-1)/2
                do ntg=1,ntg1
                   ptg(ntg)=ptg(2*ntg-1)
                   pet(ntg)=pet(2*ntg-1)
                   pqt(ntg)=pqt(2*ntg-1)
                   q_eng(ntg)=q_eng(2*ntg-1)
                   do nsa=1,nsamax
                      pnt(nsa,ntg)=pnt(nsa,2*ntg-1)
                      pwt(nsa,ntg)=pwt(nsa,2*ntg-1)
                      pwt2(nsa,ntg)=pwt2(nsa,2*ntg-1)
                      ptt(nsa,ntg)=ptt(nsa,2*ntg-1)
                      pit(nsa,ntg)=pit(nsa,2*ntg-1)
                      pirt(nsa,ntg)=pirt(nsa,2*ntg-1)
                      ppct(nsa,ntg)=ppct(nsa,2*ntg-1)
                      ppst(nsa,ntg)=ppst(nsa,2*ntg-1)
                      pplt(nsa,ntg)=pplt(nsa,2*ntg-1)
                      ppwt(nsa,ntg)=ppwt(nsa,2*ntg-1)
                      ppet(nsa,ntg)=ppet(nsa,2*ntg-1)
                      plht(nsa,ntg)=plht(nsa,2*ntg-1)
                      pfwt(nsa,ntg)=pfwt(nsa,2*ntg-1)
                      pect(nsa,ntg)=pect(nsa,2*ntg-1)
                      pwrt(nsa,ntg)=pwrt(nsa,2*ntg-1)
                      pwmt(nsa,ntg)=pwmt(nsa,2*ntg-1)
                      ptt3(nsa,ntg)=ptt3(nsa,2*ntg-1)
                      pitt(nsa,ntg)=pitt(nsa,2*ntg-1)
                      pwtt(nsa,ntg)=pwtt(nsa,2*ntg-1)
                      pnt2(nsa,ntg)=pnt2(nsa,2*ntg-1)
                      pwt2(nsa,ntg)=pwt2(nsa,2*ntg-1)
                      pwtd(nsa,ntg)=pwtd(nsa,2*ntg-1)
                      ptt2(nsa,ntg)=ptt2(nsa,2*ntg-1)
                      pit2(nsa,ntg)=pit2(nsa,2*ntg-1)
                      pspt(nsa,ntg)=pspt(nsa,2*ntg-1)
                      pspbt(nsa,ntg)=pspbt(nsa,2*ntg-1)
                      pspft(nsa,ntg)=pspft(nsa,2*ntg-1)
                      psplt(nsa,ntg)=psplt(nsa,2*ntg-1)
                      pdr(nsa,ntg)  =pdr(nsa,2*ntg-1)
                      pndr(nsa,ntg) =pndr(nsa,2*ntg-1)
                      ptt_bulk(nsa,ntg)=ptt_bulk(nsa,2*ntg-1)
                      do nsb=1,nsbmax
                         ppct2(nsb,nsa,ntg)=ppct2(nsb,nsa,2*ntg-1)
                      end do
                   end do
                end do
                ntg1=ntg1+1
             else
                ntg1m_new=2*ntg1m
                if(ntg1m_new.gt.ntg1max) ntg1m_new=ntg1max
                allocate(tempa(ntg1m))
                call fp_adjust_ntg1_a(ptg,tempa,ntg1m_new)
                call fp_adjust_ntg1_a(pet,tempa,ntg1m_new)
                call fp_adjust_ntg1_a(pqt,tempa,ntg1m_new)
                call fp_adjust_ntg1_a(q_eng,tempa,ntg1m_new)
                deallocate(tempa)
                allocate(tempb(nsamax,ntg1m))
                call fp_adjust_ntg1_b(pnt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pwt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pwt2,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(ptt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pit,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pirt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(ppct,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(ppst,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pplt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(ppwt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(ppet,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(plht,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pfwt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pect,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pwrt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pwmt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(ptt3,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pitt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pwtt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pnt2,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pwt2,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pwtd,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(ptt2,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pit2,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pspt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pspbt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pspft,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(psplt,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pspst,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pdr,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(pndr,tempb,ntg1m_new)
                call fp_adjust_ntg1_b(ptt_bulk,tempb,ntg1m_new)
                deallocate(tempb)
                allocate(tempc(nsbmax,nsamax,ntg1m))
                call fp_adjust_ntg1_c(ppct2,tempc,ntg1m_new)
                deallocate(tempc)
                ntg1m=ntg1m_new
             endif
          endif
        end subroutine fp_adjust_ntg1
!------
        subroutine fp_adjust_ntg1_a(data,temp,ntg1m_new)
          implicit none
          real(rkind),dimension(:),pointer:: data,temp
          integer,intent(in):: ntg1m_new
          integer ntg1,ntg

          do ntg=1,ntg1m
             temp(ntg)=data(ntg)
          enddo
          deallocate(data)
          allocate(data(ntg1m_new))
          do ntg=1,ntg1m
             data(ntg)=temp(ntg)
          enddo
        end subroutine fp_adjust_ntg1_a
!------
        subroutine fp_adjust_ntg1_b(data,temp,ntg1m_new)
          implicit none
          real(rkind),dimension(:,:),pointer:: data,temp
          integer,intent(in):: ntg1m_new
          integer ntg1,ntg,nr,nsa

          do ntg=1,ntg1m
             do nsa=1,nsamax
                temp(nsa,ntg)=data(nsa,ntg)
             enddo
          enddo
          deallocate(data)
          allocate(data(nsamax,ntg1m_new))
          do ntg=1,ntg1m
             do nsa=1,nsamax
                data(nsa,ntg)=temp(nsa,ntg)
             enddo
          enddo
        end subroutine fp_adjust_ntg1_b
!-------
        subroutine fp_adjust_ntg1_c(data,temp,ntg1m_new)
          implicit none
          real(rkind),dimension(:,:,:),pointer:: data,temp
          integer,intent(in):: ntg1m_new
          integer ntg1,ntg,nr,nsa,nsb

          do ntg=1,ntg1m
             do nsb=1,nsbmax
             do nsa=1,nsamax
                temp(nsb,nsa,ntg)=data(nsb,nsa,ntg)
             enddo
             enddo
          enddo
          deallocate(data)
          allocate(data(nsbmax,nsamax,ntg1m_new))
          do ntg=1,ntg1m
             do nsb=1,nsbmax
             do nsa=1,nsamax
                data(nsb,nsa,ntg)=temp(nsb,nsa,ntg)
             enddo
             enddo
          enddo
        end subroutine fp_adjust_ntg1_c
!------
        subroutine fp_allocate_ntg2
          implicit none
          integer,save:: nrmax_save=0,nsamax_save=0,nsbmax_save=0
          integer,save:: init=0

          if((nrmax.eq.nrmax_save).and. &
             (nsamax.eq.nsamax_save).and. &
             (nsbmax.eq.nsbmax_save).and. &
             (ntg2m.eq.ntg2min)) return

          if(init.eq.0) then
             init=1
          else
             call fp_deallocate_ntg2
          endif

          ntg2m=ntg2min

          allocate(rtg(ntg2m))
          allocate(ret(nrmax,ntg2m))
          allocate(rqt(nrmax,ntg2m))
          allocate(eptr(nrmax,ntg2m))
          allocate(rate_runaway(nrmax,ntg2m))
          allocate(rnt(nrmax,nsamax,ntg2m))
          allocate(rwt(nrmax,nsamax,ntg2m))
          allocate(rtt(nrmax,nsamax,ntg2m))
          allocate(rjt(nrmax,nsamax,ntg2m))
          allocate(rjrt(nrmax,nsamax,ntg2m))
          allocate(rpct(nrmax,nsamax,ntg2m))
          allocate(rpwt(nrmax,nsamax,ntg2m))
          allocate(rpet(nrmax,nsamax,ntg2m))
          allocate(rlht(nrmax,nsamax,ntg2m))
          allocate(rfwt(nrmax,nsamax,ntg2m))
          allocate(rect(nrmax,nsamax,ntg2m))
          allocate(rwrt(nrmax,nsamax,ntg2m))
          allocate(rwmt(nrmax,nsamax,ntg2m))
          allocate(rspbt(nrmax,nsamax,ntg2m))
          allocate(rspft(nrmax,nsamax,ntg2m))
          allocate(rsplt(nrmax,nsamax,ntg2m))
          allocate(rspst(nrmax,nsamax,ntg2m))
          allocate(rpdrt(nrmax,nsamax,ntg2m))
          allocate(rndrt(nrmax,nsamax,ntg2m))
          allocate(rtt_bulk(nrmax,nsamax,ntg2m))
          allocate(rate_runaway2(nrmax,nsamax,ntg2m))
          allocate(rpct2(nrmax,nsbmax,nsamax,ntg2m))

          nrmax_save=nrmax
          nsamax_save=nsamax
          nsbmax_save=nsbmax
          return
        end subroutine fp_allocate_ntg2
!------
        subroutine fp_deallocate_ntg2
          implicit none

          deallocate(rtg)
          deallocate(ret)
          deallocate(rqt)
          deallocate(eptr)
          deallocate(rate_runaway)
          deallocate(rnt)
          deallocate(rwt)
          deallocate(rtt)
          deallocate(rjt)
          deallocate(rjrt)
          deallocate(rpct)
          deallocate(rpwt)
          deallocate(rpet)
          deallocate(rlht)
          deallocate(rfwt)
          deallocate(rect)
          deallocate(rwrt,rwmt)
          deallocate(rspbt,rspft,rsplt,rspst)
          deallocate(rpdrt,rndrt)
          deallocate(rtt_bulk)
          deallocate(rate_runaway2)
          deallocate(rpct2)

          return
        end subroutine fp_deallocate_ntg2
!------
        subroutine fp_adjust_ntg2
          implicit none
          real(rkind),dimension(:),pointer:: temp0
          real(rkind),dimension(:,:),pointer:: tempa
          real(rkind),dimension(:,:,:),pointer:: tempb
          real(rkind),dimension(:,:,:,:),pointer:: tempc
          integer:: ntg,nr,nsa,nsb,ntg2m_new

          if(ntg2.gt.ntg2m) then
             write(6,*) '# fp_adjust_ntg2:', &
                        ntg2,size(rtg),ntg2m,ntg2min,ntg2max
             if(ntg2m.ge.ntg2max) then
                ntg2=(ntg2-1)/2
                do ntg=1,ntg2
                   rtg(ntg)=rtg(2*ntg-1)
                   do nr=nrstart,nrend
                      ret(nr,ntg)=ret(nr,2*ntg-1)
                      rqt(nr,ntg)=rqt(nr,2*ntg-1)
                      eptr(nr,ntg)=eptr(nr,2*ntg-1)
                      rate_runaway(nr,ntg)=rate_runaway(nr,2*ntg-1)
                   enddo
                   do nsa=1,nsamax
                      do nr=nrstart,nrend
                         rnt(nr,nsa,ntg)=rnt(nr,nsa,2*ntg-1)
                         rwt(nr,nsa,ntg)=rwt(nr,nsa,2*ntg-1)
                         rtt(nr,nsa,ntg)=rtt(nr,nsa,2*ntg-1)
                         rjt(nr,nsa,ntg)=rjt(nr,nsa,2*ntg-1)
                         rjrt(nr,nsa,ntg)=rjrt(nr,nsa,2*ntg-1)
                         rpct(nr,nsa,ntg)=rpct(nr,nsa,2*ntg-1)
                         rpwt(nr,nsa,ntg)=rpwt(nr,nsa,2*ntg-1)
                         rpet(nr,nsa,ntg)=rpet(nr,nsa,2*ntg-1)
                         rlht(nr,nsa,ntg)=rlht(nr,nsa,2*ntg-1)
                         rfwt(nr,nsa,ntg)=rfwt(nr,nsa,2*ntg-1)
                         rect(nr,nsa,ntg)=rect(nr,nsa,2*ntg-1)
                         rwrt(nr,nsa,ntg)=rwrt(nr,nsa,2*ntg-1)
                         rwmt(nr,nsa,ntg)=rwmt(nr,nsa,2*ntg-1)
                         rspbt(nr,nsa,ntg)=rspbt(nr,nsa,2*ntg-1)
                         rspft(nr,nsa,ntg)=rspft(nr,nsa,2*ntg-1)
                         rsplt(nr,nsa,ntg)=rsplt(nr,nsa,2*ntg-1)
                         rspst(nr,nsa,ntg)=rspst(nr,nsa,2*ntg-1)
                         rpdrt(nr,nsa,ntg)=rpdrt(nr,nsa,2*ntg-1)
                         rndrt(nr,nsa,ntg)=rndrt(nr,nsa,2*ntg-1)
                         rtt_bulk(nr,nsa,ntg)=rtt_bulk(nr,nsa,2*ntg-1)
                      end do
                      do nsb=1,nsbmax
                         do nr=nrstart,nrend
                            rpct2(nr,nsb,nsa,ntg)=rpct2(nr,nsb,nsa,2*ntg-1)
                         end do
                      end do
                   end do
                end do
                ntg2=ntg2+1
             else
                ntg2m_new=2*ntg2m
                if(ntg2m_new.gt.ntg2max) ntg2m_new=ntg2max
                allocate(temp0(ntg2m))
                call fp_adjust_ntg2_0(rtg,temp0,ntg2m_new)
                deallocate(temp0)
                allocate(tempa(nrmax,ntg2m))
                call fp_adjust_ntg2_a(ret,tempa,ntg2m_new)
                call fp_adjust_ntg2_a(rqt,tempa,ntg2m_new)
                call fp_adjust_ntg2_a(eptr,tempa,ntg2m_new)
                call fp_adjust_ntg2_a(rate_runaway,tempa,ntg2m_new)
                deallocate(tempa)
                allocate(tempb(nrmax,nsamax,ntg2m))
                call fp_adjust_ntg2_b(rnt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rwt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rtt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rjt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rjrt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rpct,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rpwt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rpet,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rlht,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rfwt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rect,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rwrt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rwmt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rspbt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rspft,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rsplt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rspst,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rpdrt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rndrt,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rtt_bulk,tempb,ntg2m_new)
                call fp_adjust_ntg2_b(rate_runaway2,tempb,ntg2m_new)
                deallocate(tempb)
                allocate(tempc(nrmax,nsbmax,nsamax,ntg2m))
                call fp_adjust_ntg2_c(rpct2,tempc,ntg2m_new)
                deallocate(tempc)
                ntg2m=ntg2m_new
             endif
          endif
        end subroutine fp_adjust_ntg2
!-----
        subroutine fp_adjust_ntg2_0(data,temp,ntg2m_new)
          implicit none
          real(rkind),dimension(:),pointer:: data,temp
          integer,intent(in):: ntg2m_new
          integer ntg2,ntg

          do ntg=1,ntg2m
             temp(ntg)=data(ntg)
          enddo
          deallocate(data)
          allocate(data(ntg2m_new))
          do ntg=1,ntg2m
             data(ntg)=temp(ntg)
          enddo
        end subroutine fp_adjust_ntg2_0
!-----
        subroutine fp_adjust_ntg2_a(data,temp,ntg2m_new)
          implicit none
          real(rkind),dimension(:,:),pointer:: data,temp
          integer,intent(in):: ntg2m_new
          integer ntg2,ntg,nr

          do ntg=1,ntg2m
             do nr=nrstart,nrend
                temp(nr,ntg)=data(nr,ntg)
             enddo
          enddo
          deallocate(data)
          allocate(data(nrmax,ntg2m_new))
          do ntg=1,ntg2m
             do nr=nrstart,nrend
                data(nr,ntg)=temp(nr,ntg)
             enddo
          enddo
        end subroutine fp_adjust_ntg2_a
!------
        subroutine fp_adjust_ntg2_b(data,temp,ntg2m_new)
          implicit none
          real(rkind),dimension(:,:,:),pointer:: data,temp
          integer,intent(in):: ntg2m_new
          integer ntg2,ntg,nr,nsa

          do ntg=1,ntg2m
             do nsa=1,nsamax
             do nr=nrstart,nrend
                temp(nr,nsa,ntg)=data(nr,nsa,ntg)
             enddo
             enddo
          enddo
          deallocate(data)
          allocate(data(nrmax,nsamax,ntg2m_new))
          do ntg=1,ntg2m
             do nsa=1,nsamax
             do nr=nrstart,nrend
                data(nr,nsa,ntg)=temp(nr,nsa,ntg)
             enddo
             enddo
          enddo
        end subroutine fp_adjust_ntg2_b
!--------
        subroutine fp_adjust_ntg2_c(data,temp,ntg2m_new)
          implicit none
          real(rkind),dimension(:,:,:,:),pointer:: data,temp
          integer,intent(in):: ntg2m_new
          integer ntg2,ntg,nr,nsa,nsb

          do ntg=1,ntg2m
             do nsb=1,nsbmax
             do nsa=1,nsamax
             do nr=nrstart,nrend
                temp(nr,nsb,nsa,ntg)=data(nr,nsb,nsa,ntg)
             enddo
             enddo
             enddo
          enddo
          deallocate(data)
          allocate(data(nrmax,nsbmax,nsamax,ntg2m_new))
          do ntg=1,ntg2m
             do nsb=1,nsbmax
             do nsa=1,nsamax
             do nr=nrstart,nrend
                data(nr,nsb,nsa,ntg)=temp(nr,nsb,nsa,ntg)
             enddo
             enddo
             enddo
          enddo
        end subroutine fp_adjust_ntg2_c
!-----
     end module fpcomm
