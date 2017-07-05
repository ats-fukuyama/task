
MODULE fpcomm_parm

      use bpsd_kinds
      use bpsd_constants
      use commpi
      use plcomm
      implicit none

      public

!     --- input parameters ---

      integer,parameter:: kind8=rkind
      integer,parameter:: NBEAMM=20
      real(rkind),parameter:: rkev=aee*1.D3

      integer:: NSAMAX,NSBMAX,NS_NSA(NSM),NS_NSB(NSM), NS_F1, NTH_F1, NR_F1
      integer:: LMAXNWR,NCMIN(NSM),NCMAX(NSM),NBEAMMAX,NSSPB(NBEAMM),NSSPF
      integer:: NPMAX,NTHMAX,NRMAX,NAVMAX,NP2MAX
      integer:: NTMAX,NTSTEP_COEF,NTSTEP_COLL
      integer:: NTG1STEP,NTG1MIN,NTG1MAX
      integer:: NTG2STEP,NTG2MIN,NTG2MAX
      integer:: MODELE,MODELA,MODELC,MODELR,MODELD,MODELS,MODELW(NSM),MODEL_DELTA_F(NSM)
      integer:: MODEL_ne_D,MODELD_RDEP,MODELD_PDEP,MODELD_EDGE,MODELD_PINCH
      integer:: MODELD_BOUNDARY,MODELD_CDBM
      integer:: MODEL_LOSS,MODEL_SYNCH,MODEL_NBI,MODEL_WAVE
      integer:: IMTX,MODEL_KSP,MODEL_PC,LMAXFP,LMAXE
      integer:: NGLINE,NGRAPH,LLMAX,LLMAX_NF,IDBGFP
      integer:: MODEL_DISRUPT,MODEL_Connor_fp,MODEL_BS,MODEL_jfp,MODEL_LNL
      integer:: MODEL_RE_pmax,MODELD_n_RE,MODEL_IMPURITY,MODEL_SINK,N_IMPU
      integer:: MODEL_EX_READ, MODEL_BULK_CONST, MODEL_CX_LOSS
      integer:: N_partition_r,N_partition_s,N_partition_p

      real(rkind):: PMAX(NSM),PMAX_BB(NSM),EMAX(NSM)
      real(rkind):: R1,DELR1,RMIN,RMAX
      real(rkind):: E0,ZEFF
      real(rkind):: PWAVE,RFDW,DELNPR,EPSNWR,REWY,DREWY,FACTWM
      real(rkind):: DEC,PEC1,PEC2,PEC3,PEC4,RFEC,DELYEC
      real(rkind):: DLH,PLH1,PLH2,RLH
      real(rkind):: DFW,PFW1,PFW2,RFW
      complex(rkind):: CEWR,CEWTH,CEWPH
      real(rkind):: RKWR,RKWTH,RKWPH
      real(rkind),dimension(NBEAMM):: SPBTOT,SPBR0,SPBRW,SPBENG,SPBANG,SPBPANG
      real(rkind):: SPFTOT,SPFR0,SPFRW,SPFENG
      real(rkind):: DRR0,DRRS,FACTOR_CDBM,DRR_EDGE,RHO_EDGE,FACTOR_DRR_EDGE
      real(rkind):: FACTOR_PINCH,deltaB_B
      real(rkind),dimension(NSM):: TLOSS
      real(rkind):: DELT,RIMPL,EPSFP,EPSM,EPSE,EPSDE,H0DE
      real(rkind):: PGMAX,RGMAX,RGMIN
      real(rkind):: T0_quench,tau_quench,tau_mgi
      real(rkind):: time_quench_start,RJPROF1,RJPROF2
      real(rkind):: v_RE,target_zeff,SPITOT,FACT_BULK
      real(rkind):: RN_NEU0, RN_NEUS ! temporal 

!     for read experiment data
      CHARACTER(len=80):: EG_NAME_TMS, EG_NAME_CX
      real(rkind),dimension(:,:),pointer:: read_tms_double, read_cx_double !!    containers of profile data
      integer,dimension(:,:),pointer:: read_tms_int
      real(rkind),dimension(:),pointer:: cte_fit, cti_fit !!    fitting param
      real(rkind),dimension(:),pointer:: cne_fit !!    fitting param
      integer:: nend_tms, nend_cx !!    # time step
      real(rkind),dimension(:),pointer:: RNE_EXP, RTE_EXP, RTI_EXP !! container of exp. profile data at TIMEFP
      real(rkind):: time_exp_offset, RNE_EXP_EDGE, RTE_EXP_EDGE, RTI_EXP_EDGE

!     for read FIT3D result
      CHARACTER(len=80):: SV_FILE_NAME_H, SV_FILE_NAME_D
      double precision,dimension(:),pointer:: time_grid_fit_H, time_grid_fit_D
      integer:: ntmax_fit_H, ntmax_fit_D
      integer,dimension(:,:),pointer:: I_FIT_H, I_FIT_D
      double precision,dimension(:),pointer:: D_FIT_H, D_FIT_D
      integer,dimension(:,:),pointer:: number_of_lines_fit_H, number_of_lines_fit_D

END module fpcomm_parm

module fpcomm
!
      use fpcomm_parm
      use libmpi,ONLY: mtx_mpi_type
      implicit none

      public

!      --- internal variables ---

!         +++ mpi and petsc variables +++
      TYPE(mtx_mpi_type):: comm_nr,comm_nsa,comm_np,&
           comm_nrnp,comm_nsanp,comm_nsanr
      integer:: imtxsize,imtxwidth,imtxstart,imtxend
      integer:: NTG1M,NTG2M
      integer:: nrstart,nrend,nmstart,nmend
      integer:: NSASTART,NSAEND,NPSTART,NPEND
      integer:: NPENDWM,NPENDWG,NPSTARTW,NRSTARTW,NRENDWM,NRENDWG
      integer:: MODELD_temp
      integer,dimension(:),POINTER:: mtxlen,mtxpos
      integer,dimension(:),POINTER:: savlen
      integer,dimension(:,:),POINTER:: savpos
      integer,dimension(:,:),POINTER:: Rank_Partition_Data

      integer::ISAVE
      integer,dimension(NSM):: nsb_nsa,nsa_nsb
      real(rkind):: DELR, DELTH
      real(rkind):: TIMEFP
      real(rkind),dimension(:),POINTER :: DELP
      real(rkind),dimension(:),POINTER :: &
           RNFP0,RNFPS,RTFP0,RTFPS,AMFP,AEFP,PTFP0,VTFP0, &
           AEFD,AMFD,PTFD0,VTFD0,THETA0,RNFD0,RTFD0,RTFDS, &
           RN0_MGI
      integer:: NTG1,NTG2
      real(rkind):: TVOLR
      real(rkind):: PX
      integer:: NRX,NTHX
      integer,dimension(6):: NSA1_NF,NSA2_NF,NSB1_NF,NSB2_NF
      real(rkind),dimension(6):: ENG1_NF,ENG2_NF

      real(rkind),dimension(:,:,:),POINTER :: & ! (NTHM,NPM,NRM)
           F,F1
      integer,dimension(:),POINTER :: & ! (NRM)
           ITL,ITU, ITL_judge, ITLG_judge
      integer,dimension(:),POINTER :: & ! (NRM)
           ITLG,ITUG,ITLG_RG,ITUG_RG
      real(rkind),dimension(:),POINTER :: & ! (NRM)
           volr
      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSAM)
           rlamdag,ETAMG,ETAMG_RG,ETAGG_RG,RLAMDAG_RG
      real(rkind),dimension(:),POINTER :: & ! (NRM)
           RG,RM
      real(rkind),dimension(:,:),POINTER :: & ! (NPM:NSAM)
           PG,PM
      real(rkind),dimension(:),POINTER :: & ! (NPM)
           PGB,PMB
      real(rkind),dimension(:),POINTER :: & ! (NTHM)
           THG,THM
      real(rkind),dimension(:),POINTER :: & ! (NTHM)
           RE_PITCH, RLAMDA_NRMAXP1, ETAM_NRMAXP1, ETAG_NRMAXP1

      real(rkind),dimension(:),POINTER :: & ! (NRM)
           BP,QR,RJ1,E1,RJ2,E2,BPG,BPM,QLM,QLG, &
           EP,EM,EM_W, &
           RN_disrupt, RN_runaway, Rj_ohm, RJ_runaway, conduct_sp, &
           SIGMA_SPP, SIGMA_SPM, ER_drei, ER_crit, Rconnor, LNL_GL, RFP_ava, &
           RFPL, RFP, RP_crit, RT_quench,RT_quench_f,previous_rate, RJ_bs, &
           previous_rate_p, rn_drei, RJ_bsm, RN_runaway_M, R_djdt, &
           previous_rate_G, previous_rate_p_G
      real(rkind),dimension(:),POINTER :: & ! (NRM)
           EPSRM,EPSRG,EPSRM2,EPSRG2
      real(rkind),dimension(:),POINTER :: & ! (NRM)
           EPSRMX,EPSRGX
      real(rkind),dimension(:,:,:),POINTER :: & ! (NTHM,NPM,NSBM)
           VOLP
      real(rkind),dimension(:,:),POINTER :: & ! (NTHM,NRMP)
           ETAG,ETAM,RLAMDA,RLAMDC,ETAM_RG,ETAG_RG,RLAMDA_RG,RLAMDC_G
      real(rkind),dimension(:),POINTER:: & !(NR)
           RFSAD,RFSADG, RFSADG_RG, RATIO_NAVMAX, A_chi0, Line_Element

      real(rkind),dimension(:),POINTER :: & ! (NTHM)
           SING,COSG,SINM,COSM

      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NTHM,NPM,NRM,NSBM)
           FNS
      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NTHM,NPM,NRS:NRE,NSAM)
           FNS0,FNSP,FNSM,FNSP_DEL,FNSP_MXWL
      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NTHM,NPM,NRM,NSBM)
           FNS_L
      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NTHM,NPM,NRM,NSBMAX)
           FNSB

      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSAM)
           RNFP,RTFP,PTFP,VTFP,THETA,DKBSR, POST_tau_ta,RNFP_G,RTFP_G,RT_T
      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSBM)
           RNFD,RTFD,PTFD,VTFD, RN_MGI, RN_MGI_G
      real(rkind),dimension(:,:,:),POINTER :: & ! (NRM,NSBM,NSAM)
           RNUF,RNUD,LNLAM,POST_LNLAM_f,POST_LNLAM
      real(rkind),dimension(:,:,:),POINTER :: & ! (NTHM,NPM,NSAM)
           FS0,FS2,FS1
      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NTHM,NPM,NRM,NSAM)
           WEIGHP,WEIGHT,WEIGHR,WEIGHR_G
      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NTHM,NPM,NRM,NSAM)
           DPP,DPT,DTP,DTT,FPP,FTH,DRR,FRR,SPP,PPL, &
           FEPP,FETH,DCPP,DCPT,DCTP,DCTT,FCPP,FCTH, &
           DWPP,DWPT,DWTP,DWTT,DWLHPP,DWLHPT,DWFWPP,DWFWPT, &
           DWECPP,DWECPT,SPPB,SPPF,SPPS,SPPI,DWICPP,DWICPT, &
           DWPP_P, DWPT_P, DWTP_P, DWTT_P, &
           DWICPP_P, DWICPT_P, DWECPP_P, DWECPT_P, &
           DWECTP, DWECTT, DCPPB, DCPTB, FCPPB, &
           FSPP, FSTH, DLPP, FLPP, SPPL, SPPL_CX
      real(rkind),dimension(:,:,:),POINTER :: SPPD
      real(rkind),dimension(:,:,:,:,:),POINTER :: &
           DCPP2,DCPT2,DCTP2,DCTT2,FCPP2,FCTH2  !(NTHM,NPM,NRM,NSAM,NSBM)
      real(rkind),dimension(:,:,:,:,:),POINTER :: &
           DCPP2B,DCPT2B,FCPP2B  !(NTHM,NPM,NRM,NSAM,NSBM)
      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSAM)
           RNSL,RJSL,RWSL,RPCSL,RPWSL,RPESL,RLHSL,RFWSL,RECSL,RWS123L, &
           RSPBL,RSPFL,RSPSL,RSPLL,RPDR,RNDR, RTL_BULK, RT_BULK, RICSL,&
           RDIDTL, RJSRL, RPSSL, RPLSL, RSPSL_CX
      real(rkind),dimension(:,:,:),POINTER :: & ! (NPM,NRM,NSAM)
           RP_BULK,RPL_BULK
      real(rkind),dimension(:,:,:),POINTER :: & ! (NRM,NSAM,NSBM)
           RPCS2L, RPCS2L_DEL

      real(rkind),dimension(:,:,:),POINTER :: & ! (NRM,NSAM,NSBM)
           RPW_IMPL, RPWEC_IMPL, RPWIC_IMPL
      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSAM,NSBM)
           RPW_INIT, RPWEC_INIT, RPWIC_INIT

      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSAM)
           RNS,RJS,RWS,RPCS,RPWS,RPES,RLHS,RFWS,RECS,RWS123, &
           RSPB,RSPF,RSPS,RSPL, RPDRL,RNDRL, RICS, RDIDT,&
           RJSR, RPSS, RPLS, RJS_M,RSPS_CX
      real(rkind),dimension(:),POINTER :: & ! (NSAM)
           RNS_S2
      real(rkind),dimension(:,:,:),POINTER :: & ! (NRM,NSAM,NSBM)
           RPCS2, RPCS2_DEL

      real(rkind),dimension(:),POINTER :: & ! (NTG1M)
           PTG,PET,PQT,Q_ENG
      real(rkind),dimension(:,:),POINTER :: & ! (NTG1M,NSAM)
           PNT,PWT,PTT,PIT,PPCT,PPWT,PPET,PLHT,PFWT,PECT,PTT3, &
           PITT,PWTT,PICT,PIRT, PPST, PPLT
      real(rkind),dimension(:,:),POINTER :: & ! (NTG1M,NSAM)
           PSPT,PSPBT,PSPFT,PSPLT,PSPST
      real(rkind),dimension(:,:),POINTER :: & ! (NTG1M,NSAM)
           PNT2,PWT2,PTT2,PIT2,PWTD,PDR,PNDR,PTT_BULK
      real(rkind),dimension(:,:,:),POINTER :: & ! (NTG1M,NSAM,NSBM)
           PPCT2
      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSM)
           RIPP

      real(rkind),dimension(:),POINTER :: & ! (NTG2M)
           RTG
      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NTG2M)
           RET,RQT,RATE_RUNAWAY
      real(rkind),dimension(:,:,:),POINTER :: & ! (NRM,NTG2M,NSAM)
           RNT,RWT,RTT,RJT,RPCT,RPWT,RPET,RLHT,RFWT,RECT, &
           RSPBT,RSPFT,RSPLT,RSPST,RPDRT,RNDRT,RTT_BULK,RICT,&
           RATE_RUNAWAY2,RJRT
      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NRM,NTG2M,NSAM,NSBM)
           RPCT2
      real(rkind),dimension(:,:),POINTER:: & !(NRM, NSAM)
           RN_TEMP, RT_TEMP, RN_READ, RT_READ
      integer:: NMMAX,NLMAXM
      integer,dimension(:,:,:),POINTER :: & ! (NTHM,NPM,NRM)
           NMA
      integer,dimension(:),POINTER :: & ! (NMM)
           NLMAX
      integer,dimension(:,:),POINTER :: & ! (NMM,NLM)
           LL
      real(rkind),dimension(:),POINTER :: & ! (NMM)
           DL,BM
      real(rkind),dimension(:,:),POINTER :: & ! (NMM,NLM)
           AL
      real(rkind),dimension(:),POINTER :: & ! (NRM*NTHM*NPM)
           FM, FM_shadow_m, FM_shadow_p!,BMTOT

      real(rkind),dimension(:,:,:,:,:),POINTER :: & 
           SIGMAV_NF ! (NTHMAX+1,NPMAX+1,NTHMAX+1,NPMAX+1,6)
      real(rkind),dimension(:,:,:,:,:),POINTER :: &
           SIGMAV_LG ! (0:LLMAX_NF,NPSTART:NPEND,0:LLMAX_NF,NPSTART:NPEND,6)
      real(rkind),dimension(:,:),POINTER :: &
           PL_NF ! (0:LLMAX_NF,NTHMAX)
      real(rkind),dimension(:,:),POINTER :: & 
           RATE_NF, RATE_NF_BB ! (NRSTART:NREND,6)
      real(rkind),dimension(:,:,:,:),POINTER :: & 
           RATE_NF_D1, RATE_NF_D2 ! (NTHMAX,NPMAX,NRSTART:NREND,6)
      real(rkind),dimension(:,:),POINTER :: & ! (NPM:NSAM)
           PG2,PM2
      real(rkind),dimension(:),POINTER :: & ! (NSAM)
           DEPS_SS, RPDRS, RNDRS,tau_ta0,E_drei0,E_crit0,POST_tau_ta0_f
      integer:: N_IMPL
      real,dimension(10):: gut_comm
      real(rkind),dimension(:,:),POINTER:: EPTR
      real(rkind):: E_EDGEM, SIGP_E, RN_E, RT_E, RLNRL_E
      real(rkind):: pc_runaway
      real(rkind):: Zeff_imp, tauE0_NB_e, tauE0_NB_i, tauE0_NB
      integer:: NPC_runaway
      integer:: nt_init, N_f1
      integer:: ierr_g
      integer,dimension(:,:),POINTER :: & ! (NRM,NSM)
           NP_BULK

      contains

        subroutine fp_allocate
          implicit none
          integer,save:: NRSTART_save=0,NREND_save=0,NRMAX_save=0
          integer,save:: NTHMAX_save=0,NPMAX_save=0
          integer,save:: NSAMAX_save=0,NSBMAX_save=0
          integer,save:: init=0

          if(init.eq.0) then
             init=1
          else
             if((NPMAX.eq.NPMAX_save).and. &
                (NTHMAX.eq.NTHMAX_save).and. &
                (NRMAX.eq.NRMAX_save).and. &
                (NRSTART.eq.NRSTART_save).and. &
                (NREND.eq.NREND_save).and. &
                (NSAMAX.eq.NSAMAX_save).and. &
                (NSBMAX.eq.NSBMAX_save)) return

             call fp_deallocate
          endif

          allocate( F(NTHMAX,NPSTARTW:NPENDWM,NRSTART:NREND))
          allocate(F1(NTHMAX,NPSTARTW:NPENDWM,NRSTART:NREND))

          allocate(RG(NRMAX+1),RM(NRMAX),VOLR(NRMAX))
          allocate(RLAMDAG(NTHMAX,NRMAX+1),RLAMDAG_RG(NTHMAX,NRMAX+1))
          allocate(ETAMG(NTHMAX,NRMAX+1),ETAMG_RG(NTHMAX,NRMAX+1))
          allocate(BP(NRMAX+1),QR(NRMAX))
          allocate(BPG(NRMAX+1),BPM(NRMAX+1))
          allocate(QLG(NRMAX+1),QLM(NRMAX+1))
          allocate(RJ1(NRMAX),E1(NRMAX))
          allocate(RJ2(NRMAX),E2(NRMAX))
          allocate(EP(NRSTART:NREND), EM(NRSTART:NREND), EM_W(NRSTART-1:NREND+1))

          allocate(EPSRM(NRMAX+1),EPSRG(NRMAX+1))
          allocate(EPSRM2(NRMAX+1),EPSRG2(NRMAX+1))
          allocate(EPSRMX(NRMAX+1),EPSRGX(NRMAX+1))
          allocate(ITL(NRMAX+1),ITU(NRMAX+1))
          allocate(ITLG(NRMAX+1),ITUG(NRMAX+1))
          allocate(ITLG_RG(NRMAX+1),ITUG_RG(NRMAX+1))
          allocate(ITL_judge(NRMAX+1),ITLG_judge(NRMAX+1))

          allocate(PG(NPMAX+1,NSBMAX),PM(NPMAX,NSBMAX))
          allocate(THG(NTHMAX+1),THM(NTHMAX))
          allocate(DELP(NSBMAX))
          allocate(PG2(NP2MAX+1,NSBMAX),PM2(NP2MAX,NSBMAX))

!          allocate(VOLP(NTHMAX,NPMAX,NSBMAX))
          allocate(VOLP(NTHMAX,NPSTART:NPEND,NSBMAX))
          allocate(ETAG(NTHMAX+1,NRSTART:NREND+1),ETAM(NTHMAX,NRSTART:NREND+1))
          allocate(ETAG_RG(NTHMAX+1,NRSTART:NREND+1),ETAM_RG(NTHMAX,NRSTART:NREND+1))
          allocate(RLAMDA(NTHMAX,NRSTART:NREND),RLAMDC(NTHMAX+1,NRSTART:NREND))
          allocate(RLAMDA_NRMAXP1(NTHMAX),ETAM_NRMAXP1(NTHMAX),ETAG_NRMAXP1(NTHMAX+1))

          allocate(RFSAD(NRSTART:NREND),RFSADG(NRMAX+1))
          allocate(RFSADG_RG(NRMAX+1))
          allocate(RATIO_NAVMAX(NRSTART:NREND))
          allocate(A_chi0(NRSTART:NREND))
          allocate(Line_Element(NRMAX+1))

          allocate(RLAMDA_RG(NTHMAX,NRSTART:NRENDWG),RLAMDC_G(NTHMAX+1,NRSTART:NREND))
          allocate(SING(NTHMAX+1),COSG(NTHMAX+1))
          allocate(SINM(NTHMAX),COSM(NTHMAX))

!          allocate(FNS(NTHMAX,NPMAX,NRMAX,NSBMAX))
          allocate(FNS(NTHMAX,NPMAX,NRMAX,NSAMAX))

          allocate(FNS0(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND))
          allocate(FNSP(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND)) 
          allocate(FNSM(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND))
          allocate(FNSP_DEL(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND)) 
          allocate(FNSP_MXWL(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND)) 
          allocate(FNSB(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSBMAX)) ! backgroud f

          allocate(FS0(NTHMAX,NPSTARTW:NPENDWM,NSAMAX))
          allocate(FS2(NTHMAX,NPSTARTW:NPENDWM,NSAMAX))
          allocate(FS1(NTHMAX,NPSTARTW:NPENDWM,NSAMAX))

          allocate(RNFP0(NSAMAX),RNFPS(NSAMAX))
          allocate(RTFP0(NSAMAX),RTFPS(NSAMAX))
          allocate(RTFD0(NSBMAX),RTFDS(NSBMAX))
          allocate(AMFP(NSAMAX),AEFP(NSAMAX))
          allocate(PTFP0(NSAMAX),VTFP0(NSAMAX))
          allocate(RNFD0(NSBMAX))
          allocate(AEFD(NSBMAX),AMFD(NSBMAX))
          allocate(PTFD0(NSBMAX),VTFD0(NSBMAX))
          allocate(THETA0(NSBMAX))

          allocate(RNFP(NRSTART:NREND+1,NSAMAX),RTFP(NRSTART:NREND+1,NSAMAX))
          allocate(RNFP_G(NRSTART:NRENDWG,NSAMAX),RTFP_G(NRSTART:NRENDWG,NSAMAX))
          allocate(RT_T(NRMAX,NSBMAX)) ! 
          allocate(PTFP(NRSTART:NREND+1,NSAMAX),VTFP(NRSTART:NREND+1,NSAMAX))
          allocate(THETA(NRSTART:NREND,NSBMAX),DKBSR(NRSTART:NREND,NSAMAX))
          allocate(WEIGHP(NTHMAX  ,NPSTART:NPENDWG,NRSTART:NREND+1,NSAMAX))
          allocate(WEIGHT(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND+1,NSAMAX))
          allocate(WEIGHR(NTHMAX,NPSTART:NPEND,NRSTART:NREND+1,NSAMAX))
          IF(MODELD.ne.0)THEN
             allocate(WEIGHR_G(NTHMAX,NPMAX,NRMAX+1,NSAMAX))
          END IF
          allocate(RNFD(NRSTART:NREND+1,NSBMAX),RTFD(NRSTART:NREND+1,NSBMAX))
          allocate(RN0_MGI(NSBMAX))
          allocate(RN_MGI(NRSTART:NREND,NSBMAX))
          allocate(RN_MGI_G(NRMAX,NSBMAX))
          allocate(PTFD(NRSTART:NREND+1,NSBMAX),VTFD(NRSTART:NREND+1,NSBMAX))
          allocate(RNUF(NRMAX,NSBMAX,NSBMAX))
          allocate(RNUD(NRMAX,NSBMAX,NSBMAX))
          allocate(LNLAM(NRSTART:NREND+1,NSBMAX,NSBMAX))

          allocate(DPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DPT(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DTP(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))
          allocate(DTT(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))
          allocate(FPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(FTH(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))
          allocate(DRR(NTHMAX  ,NPSTART:NPEND,NRSTART:NRENDWG,NSAMAX))
          allocate(FRR(NTHMAX  ,NPSTART:NPEND,NRSTART:NRENDWG,NSAMAX))

          allocate(FEPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(FETH(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))

          allocate(DCPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DCPT(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DCTP(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))
          allocate(DCTT(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))
          allocate(FCPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(FCTH(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))

          allocate(DLPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(FLPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))

          allocate(FSPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(FSTH(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))

          allocate(DCPP2(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSBMAX,NSAMAX))
          allocate(DCPT2(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSBMAX,NSAMAX))
          allocate(DCTP2(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSBMAX,NSAMAX))
          allocate(DCTT2(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSBMAX,NSAMAX))
          allocate(FCPP2(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSBMAX,NSAMAX))
          allocate(FCTH2(NTHMAX+1,NPSTARTW:NPENDWG,NRSTART:NRENDWM,NSBMAX,NSAMAX))

!          allocate(DCPPB(4,NPMAX+1,NRSTART:NREND+1,NSAMAX))
!          allocate(DCPTB(4,NPMAX+1,NRSTART:NREND+1,NSAMAX))
!          allocate(FCPPB(4,NPMAX+1,NRSTART:NREND+1,NSAMAX))

!          allocate(DCPP2B(4,NPMAX+1,NRSTART:NREND+1,NSBMAX,NSAMAX))
!          allocate(DCPT2B(4,NPMAX+1,NRSTART:NREND+1,NSBMAX,NSAMAX))
!          allocate(FCPP2B(4,NPMAX+1,NRSTART:NREND+1,NSBMAX,NSAMAX))

          allocate(DWPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DWPT(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DWTP(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))
          allocate(DWTT(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))
!          allocate(DWPP(NTHMAX  ,NPMAX+1,NRSTART:NREND+1,NSAMAX))
!          allocate(DWPT(NTHMAX  ,NPMAX+1,NRSTART:NREND+1,NSAMAX))
!          allocate(DWTP(NTHMAX+1,NPMAX  ,NRSTART:NREND+1,NSAMAX))
!          allocate(DWTT(NTHMAX+1,NPMAX  ,NRSTART:NREND+1,NSAMAX))
          allocate(DWLHPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DWLHPT(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DWFWPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DWFWPT(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DWECPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DWECPT(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DWECTP(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))
          allocate(DWECTT(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))
          allocate(DWICPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          allocate(DWICPT(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          IF(MODEL_WAVE.ne.0)THEN
             allocate(DWPP_P(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
             allocate(DWPT_P(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
             allocate(DWTP_P(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))
             allocate(DWTT_P(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX))

             allocate(DWECPP_P(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
             allocate(DWECPT_P(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
             allocate(DWICPP_P(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
             allocate(DWICPT_P(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX))
          END IF
          IF(MODEL_DISRUPT.ne.0)THEN
             allocate(ER_drei(NRMAX),ER_crit(NRMAX),Rconnor(NRMAX),lnl_gl(NRMAX),RP_crit(NRMAX))
             allocate(RN_disrupt(NRMAX),RN_runaway(NRMAX), RN_drei(NRMAX),RN_runaway_M(NRMAX))
             allocate(Rj_ohm(NRMAX),RJ_runaway(NRMAX),RJ_bs(NRMAX),R_djdt(NRMAX))
             allocate(RJ_bsm(NRSTART:NREND))
             allocate(previous_rate(nrstart:nrend),previous_rate_p(nrstart:nrend))
             allocate(previous_rate_G(nrmax),previous_rate_p_G(nrmax))
             allocate(SIGMA_SPP(NRSTART:NREND),SIGMA_SPM(NRSTART:NREND))
             allocate(RE_PITCH(NTHMAX))
             allocate(RT_quench(NRMAX),RT_quench_f(NRMAX))          
             allocate(POST_LNLAM_f(NRSTART:NREND+1,NSBMAX,NSBMAX))
             allocate(POST_LNLAM(NRSTART:NREND+1,NSBMAX,NSBMAX))
             allocate(POST_tau_ta0_f(NSAMAX))
             allocate(POST_tau_ta(NRMAX,NSAMAX))
             allocate(E_drei0(NSAMAX),E_crit0(NSAMAX))
          END IF
          allocate(conduct_sp(NRMAX))

          allocate(SPP (NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX))
          allocate(PPL (NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX))
          allocate(SPPB(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX))
          allocate(SPPF(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX))
          allocate(SPPS(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX))
          allocate(SPPI(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX))
          allocate(SPPD(NTHMAX,NPSTART:NPEND,NSAMAX))
          allocate(SPPL(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX))
          allocate(SPPL_CX(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX))

          allocate(RN_TEMP(NRMAX,NSMAX),RT_TEMP(NRMAX,NSMAX) )
          IF(MODEL_EX_READ.ne.0)THEN
             allocate(RN_READ(NRMAX,NSMAX),RT_READ(NRMAX,NSMAX) )
             allocate(RNE_EXP(NRMAX), RTE_EXP(NRMAX), RTI_EXP(NRMAX))
             allocate(cte_fit(5),cne_fit(6),cti_fit(5))
          END IF
          allocate(RNSL(NRSTART:NREND,NSAMAX),RJSL(NRSTART:NREND,NSAMAX))
          allocate(RFPL(NRSTART:NREND),RJSRL(NRSTART:NREND,NSAMAX))
          allocate(RWSL(NRSTART:NREND,NSAMAX),RWS123L(NRSTART:NREND,NSAMAX))
          allocate(RSPBL(NRSTART:NREND,NSAMAX),RSPFL(NRSTART:NREND,NSAMAX))
          allocate(RSPSL(NRSTART:NREND,NSAMAX),RSPLL(NRSTART:NREND,NSAMAX))
          allocate(RSPSL_CX(NRSTART:NREND,NSAMAX))
          allocate(RPCSL(NRSTART:NREND,NSAMAX),RPESL(NRSTART:NREND,NSAMAX))
          allocate(RPWSL(NRSTART:NREND,NSAMAX),RLHSL(NRSTART:NREND,NSAMAX))
          allocate(RFWSL(NRSTART:NREND,NSAMAX),RECSL(NRSTART:NREND,NSAMAX))
          allocate(RICSL(NRSTART:NREND,NSAMAX))
          allocate(RPCS2L(NRSTART:NREND,NSBMAX,NSAMAX))
          allocate(RPCS2L_DEL(NRSTART:NREND,NSBMAX,NSAMAX))
          allocate(RDIDTL(NRSTART:NREND,NSAMAX))
          allocate(RDIDT(NRMAX,NSAMAX))
          allocate(RPSSL(NRSTART:NREND,NSAMAX))
          allocate(RPLSL(NRSTART:NREND,NSAMAX))

          allocate(RNS(NRMAX,NSAMAX),RJS(NRMAX,NSAMAX),RJS_M(NRMAX,NSAMAX))
          allocate(RFP(NRMAX),RJSR(NRMAX,NSAMAX))
          allocate(RFP_ava(NRMAX))
          allocate(RNS_S2(NSAMAX))
          allocate(RWS(NRMAX,NSAMAX),RWS123(NRMAX,NSAMAX))
          allocate(RSPB(NRMAX,NSAMAX),RSPF(NRMAX,NSAMAX))
          allocate(RSPL(NRMAX,NSAMAX),RSPS(NRMAX,NSAMAX))
          allocate(RSPS_CX(NRMAX,NSAMAX))
          allocate(RPCS(NRMAX,NSAMAX),RPES(NRMAX,NSAMAX))
          allocate(RPWS(NRMAX,NSAMAX),RLHS(NRMAX,NSAMAX))
          allocate(RFWS(NRMAX,NSAMAX),RECS(NRMAX,NSAMAX))
          allocate(RICS(NRMAX,NSAMAX))
          allocate(RPCS2(NRMAX,NSBMAX,NSAMAX))
          allocate(RPCS2_DEL(NRMAX,NSBMAX,NSAMAX))
          allocate(RPSS(NRMAX,NSAMAX))
          allocate(RPLS(NRMAX,NSAMAX))

          allocate(RPDR(NRMAX,NSAMAX),RNDR(NRMAX,NSAMAX))
          allocate(RPDRS(NSAMAX),RNDRS(NSAMAX))
          allocate(tau_ta0(NSAMAX))
          allocate(RPDRL(NRSTART:NREND,NSAMAX),RNDRL(NRSTART:NREND,NSAMAX))
          allocate(RT_BULK(NRMAX,NSAMAX))
          allocate(RTL_BULK(NRSTART:NREND,NSAMAX))
          allocate(RP_BULK(NPMAX,NRMAX,NSAMAX))
          allocate(RPL_BULK(NPMAX,NRSTART:NREND,NSASTART:NSAEND))

          allocate(RPW_IMPL(NRSTART:NREND,NSAMAX,0:LMAXFP+1))
          allocate(RPWEC_IMPL(NRSTART:NREND,NSAMAX,0:LMAXFP+1))
          allocate(RPWIC_IMPL(NRSTART:NREND,NSAMAX,0:LMAXFP+1))
          allocate(RPW_INIT(NRSTART:NREND,NSAMAX))
          allocate(RPWEC_INIT(NRSTART:NREND,NSAMAX))
          allocate(RPWIC_INIT(NRSTART:NREND,NSAMAX))

          allocate(RIPP(NRMAX,NSAMAX))
!         NLMAXM= 8   ! this is for analysis without bounce average
!         NLMAXM=11   ! this is for analysis without radial transport
          NLMAXM=15   ! this is for analysis with a simple radial transport
          NMMAX=NRMAX*NTHMAX*NPMAX

          allocate(NMA(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM))
          allocate(NLMAX(NMSTART:NMEND))
          allocate(LL(NMSTART:NMEND,NLMAXM))
          allocate(DL(NMSTART:NMEND),BM(NMSTART:NMEND))
          allocate(AL(NMSTART:NMEND,NLMAXM))
          allocate(FM(NMSTART-NTHMAX:NMEND+NTHMAX))
          allocate(FM_shadow_m(NMSTART-NTHMAX-NTHMAX*NPMAX:NMEND+NTHMAX-NTHMAX*NPMAX))
          allocate(FM_shadow_p(NMSTART-NTHMAX+NTHMAX*NPMAX:NMEND+NTHMAX+NTHMAX*NPMAX))
!          allocate(FM(1+NTHMAX*(NPSTARTW-1)+NTHMAX*NPMAX*(NRSTARTW-1):NPMAX*NTHMAX*NRENDWM) )
!          allocate(BMTOT(NMMAX))

          IF(MODELS.eq.2.OR.MODELS.eq.3)THEN
!             IF(MODELS.eq.2) allocate(SIGMAV_NF(NTHMAX,NPMAX,NTHMAX,NPMAX,6))
             IF(MODELS.eq.2) allocate(SIGMAV_NF(NTHMAX,NPSTART:NPEND,NTHMAX,NPMAX,6))
             IF(MODELS.eq.3) allocate( &
                  SIGMAV_LG(0:LLMAX_NF,NPSTART:NPEND,0:LLMAX_NF,NPMAX,6), &
                  PL_NF(0:LLMAX_NF,NTHMAX))
             allocate(RATE_NF(NRSTART:NREND,6),RATE_NF_BB(NRSTART:NREND,6))
!             allocate(RATE_NF_D1(NTHMAX,NPMAX,NRSTART:NREND,6))
             allocate(RATE_NF_D2(NTHMAX,NPMAX,NRSTART:NREND,6))
             allocate(RATE_NF_D1(NTHMAX,NPSTART:NPEND,NRSTART:NREND,6))
!             allocate(RATE_NF_D2(NTHMAX,NPSTART:NPEND,NRSTART:NREND,6))
          ENDIF
          allocate(NP_BULK(NRMAX,NSBMAX))

          allocate(DEPS_SS(NSAMAX))
          NPMAX_save=NPMAX
          NTHMAX_save=NTHMAX
          NRMAX_save=NRMAX
          NRSTART_save=NRSTART
          NREND_save=NREND
          NSAMAX_save=NSAMAX
          NSBMAX_save=NSBMAX
          return
        end subroutine fp_allocate

        subroutine fp_deallocate
          implicit none

          deallocate(MTXLEN,MTXPOS,SAVLEN)

          deallocate(F)
          deallocate(F1)

          deallocate(RG,RM,VOLR)
          deallocate(RLAMDAG,RLAMDAG_RG)
          deallocate(ETAMG,ETAMG_RG)
          deallocate(BP,QR)
          deallocate(BPG,BPM)
          deallocate(QLG,QLM)
          deallocate(RJ1,E1,RJ2,E2)
          deallocate(EP,EM,EM_W)

          deallocate(EPSRM,EPSRG)
          deallocate(EPSRM2,EPSRG2)
          deallocate(EPSRMX,EPSRGX)
          deallocate(ITL,ITU)
          deallocate(ITLG,ITUG)
          deallocate(ITLG_RG,ITUG_RG)
          deallocate(ITL_judge,ITLG_judge)

          deallocate(PG,PM)
          deallocate(THG,THM)
          deallocate(DELP)
          deallocate(PG2,PM2)

          deallocate(VOLP)
          deallocate(ETAG,ETAM)
          deallocate(ETAG_RG,ETAM_RG)
          deallocate(RLAMDA,RLAMDC,RLAMDA_NRMAXP1)
          deallocate(ETAM_NRMAXP1, ETAG_NRMAXP1)
          deallocate(RFSAD,RFSADG)
          deallocate(RFSADG_RG)
          deallocate(RATIO_NAVMAX, A_chi0)
          deallocate(Line_Element)
          deallocate(RLAMDA_RG,RLAMDC_G)
          deallocate(SING,COSG,SINM,COSM)

          deallocate(FNS)
          deallocate(FNS0)
          deallocate(FNSP)
          deallocate(FNSM)
          deallocate(FNSP_DEL,FNSP_MXWL)
          deallocate(FNSB)
          deallocate(FS0,FS2,FS1)

          deallocate(RNFP0,RNFPS)
          deallocate(RTFP0,RTFPS,RTFD0,RTFDS)
          deallocate(AMFP,AEFP)
          deallocate(PTFP0,VTFP0)
          deallocate(RNFD0)
          deallocate(AEFD,AMFD)
          deallocate(PTFD0,VTFD0)
          deallocate(THETA0)

          deallocate(RNFP,RTFP)
          deallocate(RNFP_G,RTFP_G)
          deallocate(PTFP,VTFP)
          deallocate(THETA,DKBSR)
          deallocate(WEIGHP,WEIGHT)
          deallocate(WEIGHR)
          IF(MODELD.ne.0)THEN
             deallocate(WEIGHR)
          END IF

          deallocate(RNFD,RTFD,PTFD,VTFD)
          deallocate(RN_MGI, RN_MGI_G)
          deallocate(RN0_MGI)
          deallocate(RNUF,RNUD,LNLAM)
          deallocate(DPP,DPT)
          deallocate(DTP,DTT)
          deallocate(FPP,FTH)
          deallocate(DRR,FRR)
          deallocate(SPP,PPL)
          deallocate(FEPP,FETH)

          deallocate(FSPP,FSTH)

          deallocate(DLPP,FLPP)

          deallocate(DCPP,DCPT)
          deallocate(DCTP,DCTT)
          deallocate(FCPP,FCTH)

          IF(MODEL_WAVE.ne.0)THEN
             deallocate(DWPP,DWPT)
             deallocate(DWTP,DWTT)
             deallocate(DWPP_P,DWPT_P)
             deallocate(DWTP_P,DWTT_P)

             deallocate(DWLHPP,DWLHPT)
             deallocate(DWFWPP,DWFWPT)
             deallocate(DWECPP,DWECPT)
             deallocate(DWECTP,DWECTT)
             deallocate(DWICPP,DWICPT)
             deallocate(DWECPP_P,DWECPT_P)
             deallocate(DWICPP_P,DWICPT_P)
          END IF
          IF(MODEL_DISRUPT.ne.0)THEN
             deallocate(ER_drei, ER_crit,RP_crit)
             deallocate(previous_rate, previous_rate_p)
             deallocate(previous_rate_G, previous_rate_p_G)
             deallocate(RT_quench, RT_quench_f, conduct_sp)
             deallocate(RE_PITCH)
             deallocate(RN_disrupt, RN_runaway, Rj_ohm, RJ_runaway, RN_drei,RN_runaway_M)
             deallocate(RJ_bs, RJ_bsm, R_djdt)
             deallocate(SIGMA_SPP,SIGMA_SPM)
             deallocate(POST_LNLAM_f,POST_LNLAM)
             deallocate(E_drei0,E_crit0)
             deallocate(POST_tau_ta0_f)
             deallocate(POST_tau_ta)
          END IF
          deallocate(tau_ta0)

          deallocate(SPPB,SPPF,SPPS,SPPD,SPPI,SPPL)
          deallocate(SPPL_CX)
          

          deallocate(DCPP2,DCPT2)
          deallocate(DCTP2,DCTT2)
          deallocate(FCPP2,FCTH2)
          
!          deallocate(DCPPB,DCPTB,FCPPB)
!          deallocate(DCPP2B,DCPT2B,FCPP2B)

          deallocate(RN_TEMP,RT_TEMP,RT_T)
          IF(MODEL_EX_READ.ne.0)THEN
             deallocate(RN_READ, RT_READ)
             deallocate(RNE_EXP, RTE_EXP, RTI_EXP)
             deallocate(cte_fit, cne_fit, cti_fit)
          END IF
          deallocate(RNSL,RJSL,RWSL,RWS123L,RFPL,RJSRL)
          deallocate(RSPBL,RSPFL,RSPSL,RSPLL,RPCSL,RPESL,RSPSL_CX)
          deallocate(RLHSL,RFWSL,RECSL,RICSL,RPCS2L,RPCS2L_DEL)
          deallocate(RDIDT, RDIDTL)
          deallocate(RPSSL, RPLSL)

          deallocate(RNS,RJS,RJS_M,RFP,RJSR,Rconnor,RFP_ava)
          deallocate(RNS_S2)
          deallocate(RWS,RWS123)
          deallocate(RSPB,RSPF)
          deallocate(RSPL,RSPS,RSPS_CX)
          deallocate(RPCS,RPES)
          deallocate(RPWS,RLHS)
          deallocate(RFWS,RECS,RICS,RPCS2,RPCS2_DEL)
          deallocate(RPSS, RPLS)

          deallocate(RPDR,RNDR,RPDRS,RNDRS)
          deallocate(RPDRL,RNDRL,RT_BULK,RTL_BULK)
          deallocate(RP_BULK)
          
          deallocate(RPW_IMPL,RPWEC_IMPL,RPWIC_IMPL)
          deallocate(RPW_INIT,RPWEC_INIT,RPWIC_INIT)
          deallocate(RIPP)

          deallocate(NMA)
          deallocate(NLMAX,LL)
          deallocate(DL,BM)
          deallocate(AL)
          deallocate(FM)
          deallocate(FM_shadow_m, FM_shadow_p)
!          deallocate(BMTOT)
          deallocate(NP_BULK)

          IF(MODELS.eq.2.OR.MODELS.eq.3)THEN
             IF(MODELS.eq.2) deallocate(SIGMAV_NF)
             IF(MODELS.eq.3) deallocate(SIGMAV_LG,PL_NF)
             deallocate(RATE_NF,RATE_NF_D1,RATE_NF_D2,RATE_NF_BB)
          END IF
          deallocate(DEPS_SS)
          return

        end subroutine fp_deallocate

        subroutine fp_allocate_ntg1
          implicit none
          integer,save:: NSAMAX_save=0,NSBMAX_save=0,NTG1M_save=0
          integer,save:: init=0
          
          if((NSAMAX.eq.NSAMAX_save).and. &
             (NSBMAX.eq.NSBMAX_save).and. &
             (NTG1M.eq.NTG1MIN)) return

          if(init.eq.0) then
             init=1
          else
             call fp_deallocate_ntg1
          endif

          NTG1M=NTG1MIN

          allocate(PTG(NTG1M))
          allocate(PET(NTG1M))
          allocate(PQT(NTG1M))
          allocate(Q_ENG(NTG1M))
          allocate(PNT(NSAMAX,NTG1M))
          allocate(PWT(NSAMAX,NTG1M))
          allocate(PTT(NSAMAX,NTG1M))
          allocate(PIT(NSAMAX,NTG1M))
          allocate(PIRT(NSAMAX,NTG1M))
          allocate(PPCT(NSAMAX,NTG1M))
          allocate(PPST(NSAMAX,NTG1M))
          allocate(PPLT(NSAMAX,NTG1M))
          allocate(PPWT(NSAMAX,NTG1M))
          allocate(PPET(NSAMAX,NTG1M))
          allocate(PLHT(NSAMAX,NTG1M))
          allocate(PFWT(NSAMAX,NTG1M))
          allocate(PECT(NSAMAX,NTG1M))
          allocate(PICT(NSAMAX,NTG1M))
          allocate(PTT3(NSAMAX,NTG1M))
          allocate(PITT(NSAMAX,NTG1M))
          allocate(PWTT(NSAMAX,NTG1M))
          allocate(PNT2(NSAMAX,NTG1M))
          allocate(PWT2(NSAMAX,NTG1M))
          allocate(PWTD(NSAMAX,NTG1M))
          allocate(PTT2(NSAMAX,NTG1M))
          allocate(PIT2(NSAMAX,NTG1M))
          allocate(PSPT(NSAMAX,NTG1M))
          allocate(PSPBT(NSAMAX,NTG1M))
          allocate(PSPFT(NSAMAX,NTG1M))
          allocate(PSPLT(NSAMAX,NTG1M))
          allocate(PSPST(NSAMAX,NTG1M))
          allocate(PPCT2(NSBMAX,NSAMAX,NTG1M))
          allocate(PDR(NSAMAX,NTG1M),PNDR(NSAMAX,NTG1M))
          allocate(PTT_BULK(NSAMAX,NTG1M))

          NSAMAX_save=NSAMAX
          NSBMAX_save=NSBMAX
          return
        end subroutine fp_allocate_ntg1

        subroutine fp_deallocate_ntg1
          implicit none

          deallocate(PTG)
          deallocate(PET)
          deallocate(PQT)
          deallocate(Q_ENG)
          deallocate(PNT)
          deallocate(PWT)
          deallocate(PTT)
          deallocate(PIT)
          deallocate(PIRT)
          deallocate(PPCT)
          deallocate(PPST)
          deallocate(PPLT)
          deallocate(PPWT)
          deallocate(PPET)
          deallocate(PLHT)
          deallocate(PFWT)
          deallocate(PECT)
          deallocate(PICT)
          deallocate(PTT3)
          deallocate(PITT)
          deallocate(PWTT)
          deallocate(PNT2)
          deallocate(PWT2)
          deallocate(PWTD)
          deallocate(PTT2)
          deallocate(PIT2)
          deallocate(PSPT,PSPBT,PSPFT,PSPLT,PSPST)
          deallocate(PPCT2)
          deallocate(PDR,PNDR)
          deallocate(PTT_BULK)

        end subroutine fp_deallocate_ntg1

        subroutine fp_adjust_ntg1
          implicit none
          real(rkind),dimension(:),POINTER:: tempA
          real(rkind),dimension(:,:),POINTER:: tempB
          real(rkind),dimension(:,:,:),POINTER:: tempC
          integer:: NTG,NSA,NSB,NTG1M_NEW

          if(NTG1.GT.NTG1M) then
             write(6,*) '# fp_adjust_ntg1:',&
                        NTG1,SIZE(PTG),NTG1M,NTG1MIN,NTG1MAX
             if(NTG1M.GE.NTG1MAX) then
                NTG1=(NTG1-1)/2
                DO NTG=1,NTG1
                   PTG(NTG)=PTG(2*NTG-1)
                   PET(NTG)=PET(2*NTG-1)
                   PQT(NTG)=PQT(2*NTG-1)
                   Q_ENG(NTG)=Q_ENG(2*NTG-1)
                   DO NSA=1,NSAMAX
                      PNT(NSA,NTG)=PNT(NSA,2*NTG-1)
                      PWT(NSA,NTG)=PWT(NSA,2*NTG-1)
                      PWT2(NSA,NTG)=PWT2(NSA,2*NTG-1)
                      PTT(NSA,NTG)=PTT(NSA,2*NTG-1)
                      PIT(NSA,NTG)=PIT(NSA,2*NTG-1)
                      PIRT(NSA,NTG)=PIRT(NSA,2*NTG-1)
                      PPCT(NSA,NTG)=PPCT(NSA,2*NTG-1)
                      PPST(NSA,NTG)=PPST(NSA,2*NTG-1)
                      PPLT(NSA,NTG)=PPLT(NSA,2*NTG-1)
                      PPWT(NSA,NTG)=PPWT(NSA,2*NTG-1)
                      PPET(NSA,NTG)=PPET(NSA,2*NTG-1)
                      PLHT(NSA,NTG)=PLHT(NSA,2*NTG-1)
                      PFWT(NSA,NTG)=PFWT(NSA,2*NTG-1)
                      PECT(NSA,NTG)=PECT(NSA,2*NTG-1)
                      PICT(NSA,NTG)=PICT(NSA,2*NTG-1)
                      PTT3(NSA,NTG)=PTT3(NSA,2*NTG-1)
                      PITT(NSA,NTG)=PITT(NSA,2*NTG-1)
                      PWTT(NSA,NTG)=PWTT(NSA,2*NTG-1)
                      PNT2(NSA,NTG)=PNT2(NSA,2*NTG-1)
                      PWT2(NSA,NTG)=PWT2(NSA,2*NTG-1)
                      PWTD(NSA,NTG)=PWTD(NSA,2*NTG-1)
                      PTT2(NSA,NTG)=PTT2(NSA,2*NTG-1)
                      PIT2(NSA,NTG)=PIT2(NSA,2*NTG-1)
                      PSPT(NSA,NTG)=PSPT(NSA,2*NTG-1)
                      PSPBT(NSA,NTG)=PSPBT(NSA,2*NTG-1)
                      PSPFT(NSA,NTG)=PSPFT(NSA,2*NTG-1)
                      PSPLT(NSA,NTG)=PSPLT(NSA,2*NTG-1)
                      PDR(NSA,NTG)  =PDR(NSA,2*NTG-1)
                      PNDR(NSA,NTG) =PNDR(NSA,2*NTG-1)
                      PTT_BULK(NSA,NTG)=PTT_BULK(NSA,2*NTG-1)
                      DO NSB=1,NSBMAX
                         PPCT2(NSB,NSA,NTG)=PPCT2(NSB,NSA,2*NTG-1)
                      END DO
                   END DO
                END DO
                NTG1=NTG1+1
             else
                NTG1M_NEW=2*NTG1M
                IF(NTG1M_NEW.GT.NTG1MAX) NTG1M_NEW=NTG1MAX
                allocate(tempA(NTG1M))
                call fp_adjust_ntg1_A(PTG,tempA,NTG1M_NEW)
                call fp_adjust_ntg1_A(PET,tempA,NTG1M_NEW)
                call fp_adjust_ntg1_A(PQT,tempA,NTG1M_NEW)
                call fp_adjust_ntg1_A(Q_ENG,tempA,NTG1M_NEW)
                deallocate(tempA)
                allocate(tempB(NSAMAX,NTG1M))
                call fp_adjust_ntg1_B(PNT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PWT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PWT2,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PTT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PIT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PIRT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PPCT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PPST,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PPLT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PPWT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PPET,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PLHT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PFWT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PECT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PICT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PTT3,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PITT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PWTT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PNT2,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PWT2,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PWTD,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PTT2,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PIT2,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PSPT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PSPBT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PSPFT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PSPLT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PSPST,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PDR,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PNDR,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PTT_BULK,tempB,NTG1M_NEW)
                deallocate(tempB)
                allocate(tempC(NSBMAX,NSAMAX,NTG1M))
                call fp_adjust_ntg1_C(PPCT2,tempC,NTG1M_NEW)
                deallocate(tempC)
                NTG1M=NTG1M_NEW
             endif
          endif
        end subroutine fp_adjust_ntg1
!------
        subroutine fp_adjust_ntg1_A(data,temp,NTG1M_NEW)
          implicit none
          real(rkind),dimension(:),POINTER:: data,temp
          integer,intent(in):: NTG1M_NEW
          integer NTG1,NTG
         
          DO NTG=1,NTG1M
             temp(NTG)=data(NTG)
          ENDDO
          deallocate(data)
          allocate(data(NTG1M_NEW))
          DO NTG=1,NTG1M
             data(NTG)=temp(NTG)
          ENDDO
        end subroutine fp_adjust_ntg1_A
!------
        subroutine fp_adjust_ntg1_B(data,temp,NTG1M_NEW)
          implicit none
          real(rkind),dimension(:,:),POINTER:: data,temp
          integer,intent(in):: NTG1M_NEW
          integer NTG1,NTG,NR,NSA
         
          DO NTG=1,NTG1M
             DO NSA=1,NSAMAX
                temp(NSA,NTG)=data(NSA,NTG)
             ENDDO
          ENDDO
          deallocate(data)
          allocate(data(NSAMAX,NTG1M_NEW))
          DO NTG=1,NTG1M
             DO NSA=1,NSAMAX
                data(NSA,NTG)=temp(NSA,NTG)
             ENDDO
          ENDDO
        end subroutine fp_adjust_ntg1_B
!-------
        subroutine fp_adjust_ntg1_C(data,temp,NTG1M_NEW)
          implicit none
          real(rkind),dimension(:,:,:),POINTER:: data,temp
          integer,intent(in):: NTG1M_NEW
          integer NTG1,NTG,NR,NSA,NSB
         
          DO NTG=1,NTG1M
             DO NSB=1,NSBMAX
             DO NSA=1,NSAMAX
                temp(NSB,NSA,NTG)=data(NSB,NSA,NTG)
             ENDDO
             ENDDO
          ENDDO
          deallocate(data)
          allocate(data(NSBMAX,NSAMAX,NTG1M_NEW))
          DO NTG=1,NTG1M
             DO NSB=1,NSBMAX
             DO NSA=1,NSAMAX
                data(NSB,NSA,NTG)=temp(NSB,NSA,NTG)
             ENDDO
             ENDDO
          ENDDO
        end subroutine fp_adjust_ntg1_C
!------
        subroutine fp_allocate_ntg2
          implicit none
          integer,save:: NRMAX_save=0,NSAMAX_save=0,NSBMAX_save=0
          integer,save:: init=0
          
          if((NRMAX.eq.NRMAX_save).and. &
             (NSAMAX.eq.NSAMAX_save).and. &
             (NSBMAX.eq.NSBMAX_save).and. &
             (NTG2M.eq.NTG2MIN)) return

          if(init.eq.0) then
             init=1
          else
             call fp_deallocate_ntg2
          endif

          NTG2M=NTG2MIN

          allocate(RTG(NTG2M))
          allocate(RET(NRMAX,NTG2M))
          allocate(RQT(NRMAX,NTG2M))
          allocate(EPTR(NRMAX,NTG2M))
          allocate(RATE_RUNAWAY(NRMAX,NTG2M))
          allocate(RNT(NRMAX,NSAMAX,NTG2M))
          allocate(RWT(NRMAX,NSAMAX,NTG2M))
          allocate(RTT(NRMAX,NSAMAX,NTG2M))
          allocate(RJT(NRMAX,NSAMAX,NTG2M))
          allocate(RJRT(NRMAX,NSAMAX,NTG2M))
          allocate(RPCT(NRMAX,NSAMAX,NTG2M))
          allocate(RPWT(NRMAX,NSAMAX,NTG2M))
          allocate(RPET(NRMAX,NSAMAX,NTG2M))
          allocate(RLHT(NRMAX,NSAMAX,NTG2M))
          allocate(RFWT(NRMAX,NSAMAX,NTG2M))
          allocate(RECT(NRMAX,NSAMAX,NTG2M))
          allocate(RICT(NRMAX,NSAMAX,NTG2M))
          allocate(RSPBT(NRMAX,NSAMAX,NTG2M))
          allocate(RSPFT(NRMAX,NSAMAX,NTG2M))
          allocate(RSPLT(NRMAX,NSAMAX,NTG2M))
          allocate(RSPST(NRMAX,NSAMAX,NTG2M))
          allocate(RPDRT(NRMAX,NSAMAX,NTG2M))
          allocate(RNDRT(NRMAX,NSAMAX,NTG2M))
          allocate(RTT_BULK(NRMAX,NSAMAX,NTG2M))
          allocate(RATE_RUNAWAY2(NRMAX,NSAMAX,NTG2M))
          allocate(RPCT2(NRMAX,NSBMAX,NSAMAX,NTG2M))

          NRMAX_save=NRMAX
          NSAMAX_save=NSAMAX
          NSBMAX_save=NSBMAX
          return
        end subroutine fp_allocate_ntg2
!------
        subroutine fp_deallocate_ntg2
          implicit none

          deallocate(RTG)
          deallocate(RET)
          deallocate(RQT)
          deallocate(EPTR)
          deallocate(RATE_RUNAWAY)
          deallocate(RNT)
          deallocate(RWT)
          deallocate(RTT)
          deallocate(RJT)
          deallocate(RJRT)
          deallocate(RPCT)
          deallocate(RPWT)
          deallocate(RPET)
          deallocate(RLHT)
          deallocate(RFWT)
          deallocate(RECT)
          deallocate(RICT)
          deallocate(RSPBT,RSPFT,RSPLT,RSPST)
          deallocate(RPDRT,RNDRT)
          deallocate(RTT_BULK)
          deallocate(RATE_RUNAWAY2)
          deallocate(RPCT2)

          return
        end subroutine fp_deallocate_ntg2
!------
        subroutine fp_adjust_ntg2
          implicit none
          real(rkind),dimension(:),POINTER:: temp0
          real(rkind),dimension(:,:),POINTER:: tempA
          real(rkind),dimension(:,:,:),POINTER:: tempB
          real(rkind),dimension(:,:,:,:),POINTER:: tempC
          integer:: NTG,NR,NSA,NSB,NTG2M_NEW

          if(NTG2.GT.NTG2M) then
             write(6,*) '# fp_adjust_ntg2:', &
                        NTG2,SIZE(RTG),NTG2M,NTG2MIN,NTG2MAX
             if(NTG2M.GE.NTG2MAX) then
                NTG2=(NTG2-1)/2
                DO NTG=1,NTG2
                   RTG(NTG)=RTG(2*NTG-1)
                   DO NR=NRSTART,NREND
                      RET(NR,NTG)=RET(NR,2*NTG-1)
                      RQT(NR,NTG)=RQT(NR,2*NTG-1)
                      EPTR(NR,NTG)=EPTR(NR,2*NTG-1)
                      RATE_RUNAWAY(NR,NTG)=RATE_RUNAWAY(NR,2*NTG-1)
                   ENDDO
                   DO NSA=1,NSAMAX
                      DO NR=NRSTART,NREND
                         RNT(NR,NSA,NTG)=RNT(NR,NSA,2*NTG-1)
                         RWT(NR,NSA,NTG)=RWT(NR,NSA,2*NTG-1)
                         RTT(NR,NSA,NTG)=RTT(NR,NSA,2*NTG-1)
                         RJT(NR,NSA,NTG)=RJT(NR,NSA,2*NTG-1)
                         RJRT(NR,NSA,NTG)=RJRT(NR,NSA,2*NTG-1)
                         RPCT(NR,NSA,NTG)=RPCT(NR,NSA,2*NTG-1)
                         RPWT(NR,NSA,NTG)=RPWT(NR,NSA,2*NTG-1)
                         RPET(NR,NSA,NTG)=RPET(NR,NSA,2*NTG-1)
                         RLHT(NR,NSA,NTG)=RLHT(NR,NSA,2*NTG-1)
                         RFWT(NR,NSA,NTG)=RFWT(NR,NSA,2*NTG-1)
                         RECT(NR,NSA,NTG)=RECT(NR,NSA,2*NTG-1)
                         RICT(NR,NSA,NTG)=RICT(NR,NSA,2*NTG-1)
                         RSPBT(NR,NSA,NTG)=RSPBT(NR,NSA,2*NTG-1)
                         RSPFT(NR,NSA,NTG)=RSPFT(NR,NSA,2*NTG-1)
                         RSPLT(NR,NSA,NTG)=RSPLT(NR,NSA,2*NTG-1)
                         RSPST(NR,NSA,NTG)=RSPST(NR,NSA,2*NTG-1)
                         RPDRT(NR,NSA,NTG)=RPDRT(NR,NSA,2*NTG-1)
                         RNDRT(NR,NSA,NTG)=RNDRT(NR,NSA,2*NTG-1)
                         RTT_BULK(NR,NSA,NTG)=RTT_BULK(NR,NSA,2*NTG-1)
                      END DO
                      DO NSB=1,NSBMAX
                         DO NR=NRSTART,NREND
                            RPCT2(NR,NSB,NSA,NTG)=RPCT2(NR,NSB,NSA,2*NTG-1)
                         END DO
                      END DO
                   END DO
                END DO
                NTG2=NTG2+1
             else
                NTG2M_NEW=2*NTG2M
                IF(NTG2M_NEW.GT.NTG2MAX) NTG2M_NEW=NTG2MAX
                allocate(temp0(NTG2M))
                call fp_adjust_ntg2_0(RTG,temp0,NTG2M_NEW)
                deallocate(temp0)
                allocate(tempA(NRMAX,NTG2M))
                call fp_adjust_ntg2_A(RET,tempA,NTG2M_NEW)
                call fp_adjust_ntg2_A(RQT,tempA,NTG2M_NEW)
                call fp_adjust_ntg2_A(EPTR,tempA,NTG2M_NEW)
                call fp_adjust_ntg2_A(RATE_RUNAWAY,tempA,NTG2M_NEW)
                deallocate(tempA)
                allocate(tempB(NRMAX,NSAMAX,NTG2M))
                call fp_adjust_ntg2_B(RNT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RWT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RTT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RJT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RJRT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RPCT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RPWT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RPET,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RLHT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RFWT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RECT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RICT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RSPBT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RSPFT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RSPLT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RSPST,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RPDRT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RNDRT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RTT_BULK,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RATE_RUNAWAY2,tempB,NTG2M_NEW)
                deallocate(tempB)
                allocate(tempC(NRMAX,NSBMAX,NSAMAX,NTG2M))
                call fp_adjust_ntg2_C(RPCT2,tempC,NTG2M_NEW)
                deallocate(tempC)
                NTG2M=NTG2M_NEW
             endif
          endif
        end subroutine fp_adjust_ntg2
!-----
        subroutine fp_adjust_ntg2_0(data,temp,NTG2M_NEW)
          implicit none
          real(rkind),dimension(:),POINTER:: data,temp
          integer,intent(in):: NTG2M_NEW
          integer NTG2,NTG
         
          DO NTG=1,NTG2M
             temp(NTG)=data(NTG)
          ENDDO
          deallocate(data)
          allocate(data(NTG2M_NEW))
          DO NTG=1,NTG2M
             data(NTG)=temp(NTG)
          ENDDO
        end subroutine fp_adjust_ntg2_0
!-----
        subroutine fp_adjust_ntg2_A(data,temp,NTG2M_NEW)
          implicit none
          real(rkind),dimension(:,:),POINTER:: data,temp
          integer,intent(in):: NTG2M_NEW
          integer NTG2,NTG,NR
         
          DO NTG=1,NTG2M
             DO NR=NRSTART,NREND
                temp(NR,NTG)=data(NR,NTG)
             ENDDO
          ENDDO
          deallocate(data)
          allocate(data(NRMAX,NTG2M_NEW))
          DO NTG=1,NTG2M
             DO NR=NRSTART,NREND
                data(NR,NTG)=temp(NR,NTG)
             ENDDO
          ENDDO
        end subroutine fp_adjust_ntg2_A
!------
        subroutine fp_adjust_ntg2_B(data,temp,NTG2M_NEW)
          implicit none
          real(rkind),dimension(:,:,:),POINTER:: data,temp
          integer,intent(in):: NTG2M_NEW
          integer NTG2,NTG,NR,NSA
         
          DO NTG=1,NTG2M
             DO NSA=1,NSAMAX
             DO NR=NRSTART,NREND
                temp(NR,NSA,NTG)=data(NR,NSA,NTG)
             ENDDO
             ENDDO
          ENDDO
          deallocate(data)
          allocate(data(NRMAX,NSAMAX,NTG2M_NEW))
          DO NTG=1,NTG2M
             DO NSA=1,NSAMAX
             DO NR=NRSTART,NREND
                data(NR,NSA,NTG)=temp(NR,NSA,NTG)
             ENDDO
             ENDDO
          ENDDO
        end subroutine fp_adjust_ntg2_B
!--------
        subroutine fp_adjust_ntg2_C(data,temp,NTG2M_NEW)
          implicit none
          real(rkind),dimension(:,:,:,:),POINTER:: data,temp
          integer,intent(in):: NTG2M_NEW
          integer NTG2,NTG,NR,NSA,NSB
         
          DO NTG=1,NTG2M
             DO NSB=1,NSBMAX
             DO NSA=1,NSAMAX
             DO NR=NRSTART,NREND
                temp(NR,NSB,NSA,NTG)=data(NR,NSB,NSA,NTG)
             ENDDO
             ENDDO
             ENDDO
          ENDDO
          deallocate(data)
          allocate(data(NRMAX,NSBMAX,NSAMAX,NTG2M_NEW))
          DO NTG=1,NTG2M
             DO NSB=1,NSBMAX
             DO NSA=1,NSAMAX
             DO NR=NRSTART,NREND
                data(NR,NSB,NSA,NTG)=temp(NR,NSB,NSA,NTG)
             ENDDO
             ENDDO
             ENDDO
          ENDDO
        end subroutine fp_adjust_ntg2_C
!-----
     end module fpcomm
