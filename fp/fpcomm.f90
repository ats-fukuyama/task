!     $Id$

      module fpcomm
!
      use bpsd
      use plcomm
      implicit none

      public

!     --- input parameters ---

      integer:: NPMAX,NTHMAX,NRMAX,NAVMAX,IMTX
      integer:: NTMAX,NTCLSTEP,LMAXE,NGLINE,LMAXNWR
      integer:: MODELE,MODELA,MODELC,MODELR,MODELD,MODELS
      integer:: LLMAX,IDBGFP
      integer:: NTG1STEP,NTG1MIN,NTG1MAX
      integer:: NTG2STEP,NTG2MIN,NTG2MAX
      integer,dimension(NSM):: MODELW
      integer:: NTEST, NGRAPH
      integer,parameter:: NBEAMM=20
      integer:: NBEAMMAX,NP2MAX
      integer:: MODEL_KSP, MODEL_PC
      integer,parameter:: kind8=rkind
      real(rkind):: PGMAX, RGMAX, RGMIN
      real(rkind):: DRR0,E0,R1,DELR1,RMIN,RMAX,DRRS
      real(rkind):: DEC,PEC1,PEC2,PEC3,PEC4,RFEC,DELYEC
      real(rkind):: DLH,PLH1,PLH2,RLH
      real(rkind):: DFW,PFW1,PFW2,RFW
      real(rkind):: RFDW,DELNPR
      complex(rkind):: CEWR,CEWTH,CEWPH
      real(rkind):: RKWR,RKWTH,RKWPH,REWY,DREWY,FACTWM
      real(rkind):: ZEFF,DELT,RIMPL,EPSM,EPSE,EPSDE,H0DE
      real(rkind):: PWAVE,EPSNWR

      integer:: nsamax,nsbmax
      integer,dimension(NSM):: ns_nsa,ns_nsb
      real(rkind),dimension(NSM):: pmax,tloss
!      real(rkind),dimension(NSM) :: SPTOT,SPR0,SPRW,SPENG,SPANG
      integer:: NSSPF
      integer,dimension(NBEAMM) :: NSSPB
      real(rkind),dimension(NBEAMM) :: SPBTOT,SPBR0,SPBRW,SPBENG,SPBANG,SPBPANG
      real(rkind) :: SPFTOT,SPFR0,SPFRW,SPFENG,SPFANG

      integer::LMAXFP
      real(rkind):: epsfp
      integer,dimension(NSM):: NCMIN, NCMAX
      
!      --- internal variables ---

!         +++ mpi and petsc variables +++
      integer:: nrank,nprocs,imtxsize,imtxwidth,imtxstart,imtxend
      integer:: nrstart,nrend,nrendx,nmstart,nmend
      integer,dimension(:),POINTER:: mtxlen,mtxpos
      
      integer::ISAVE
      integer,dimension(NSM):: nsb_nsa,nsa_nsb
      real(rkind):: DELR, DELTH
      real(rkind):: TIMEFP
      real(rkind),dimension(:),POINTER :: DELP
      real(rkind),dimension(:),POINTER :: &
           RNFP0,RNFPS,RTFP0,RTFPS,AMFP,AEFP,PTFP0,VTFP0, &
           AEFD,AMFD,PTFD0,VTFD0,THETA0,RNFD0
      integer:: NTG1,NTG2,NTG1M,NTG2M
      real(rkind):: TVOLR
      real(rkind):: PX
      integer:: NRX,NTHX
      integer,dimension(6):: NSA1_NF,NSA2_NF,NSB1_NF,NSB2_NF
      real(rkind),dimension(6):: ENG1_NF,ENG2_NF

      real(rkind),dimension(:,:,:),POINTER :: & ! (NTHM,NPM,NRM)
           F,F1
      integer,dimension(:),POINTER :: & ! (NRM)
           ITL,ITU,ITL_G,ITU_G
      integer,dimension(:),POINTER :: & ! (NRM)
           ITLG,ITUG,ITLG_G,ITUG_G
      real(rkind),dimension(:),POINTER :: & ! (NRM,NSBM)
           RCOEF, RCOEF_G, RCOEFN, RCOEFN_G
      real(rkind),dimension(:),POINTER :: & ! (NSAM)
           RCOEF1,RCOEF2,RCOEF2_G
      real(rkind),dimension(:),POINTER :: & ! (NSAM)
           RCOEFG, RCOEF_GG, RCOEFNG, RCOEFN_GG
      real(rkind),dimension(:),POINTER :: & ! (NRM)
           volr
      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSAM)
           rlamdag,ETAMG,ETAM_GG,ETAG_GG,RLAMDA_GG,ETAG_G_GL
      real(rkind),dimension(:),POINTER :: & ! (NRM)
           RG,RM
      real(rkind),dimension(:,:),POINTER :: & ! (NPM:NSAM)
           PG,PM
      real(rkind),dimension(:),POINTER :: & ! (NPM)
           PGB,PMB
      real(rkind),dimension(:),POINTER :: & ! (NTHM)
           THG,THM

      real(rkind),dimension(:),POINTER :: & ! (NRM)
           BP,QR,RJ1,E1,RJ2,E2,BPG,BPM,QLM,QLG
      real(rkind),dimension(:),POINTER :: & ! (NRM)
           EPSRM,EPSRG,EPSRM2,EPSRG2
      real(rkind),dimension(:),POINTER :: & ! (NRM)
           EPSRMX,EPSRGX
      real(rkind),dimension(:,:,:),POINTER :: & ! (NTHM,NPM,NSBM)
           VOLP
      real(rkind),dimension(:,:),POINTER :: & ! (NTHM,NRMP)
           ETAG,ETAM,RLAMDA,RLAMDC,ETAM_G,ETAG_G,RLAMDA_G,RlAMDC_G
      real(rkind),dimension(:),POINTER:: & !(NR)
           RFSAD,RFSADG, RFSAD_G, RFSAD_GG

      real(rkind),dimension(:),POINTER :: & ! (NTHM)
           SING,COSG,SINM,COSM

      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NTHM,NPM,NRM,NSBM)
           FNS
      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NTHM,NPM,NRM,NSAM)
           FNS1,FNS2,FNS22
      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NTHM,NPM,NRM,NSBM)
           FNS_L


      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSAM)
           RNFP,RTFP,PTFP,VTFP,THETA,DKBSR
      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSBM)
           RNFD,RTFD,PTFD,VTFD
      real(rkind),dimension(:,:,:),POINTER :: & ! (NRM,NSBM,NSAM)
           RNUF,RNUD
      real(rkind),dimension(:,:,:),POINTER :: & ! (NTHM,NPM,NSAM)
           FS1,FS2,FS3
      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NTHM,NPM,NRM,NSAM)
           WEIGHP,WEIGHT,WEIGHR
      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NTHM,NPM,NRM,NSAM)
           DPP,DPT,DTP,DTT,FPP,FTH,DRR,FRR,SPP,PPL, &
           FEPP,FETH,DCPP,DCPT,DCTP,DCTT,FCPP,FCTH, &
           DWPP,DWPT,DWTP,DWTT,DWLHPP,DWLHPT,DWFWPP,DWFWPT, &
           DWECPP,DWECPT,SPPB,SPPF,SPPS,DWICPP,DWICPT, &
           DWPP_P, DWPT_P, DWTP_P, DWTT_P, &
           DWICPP_P, DWICPT_P, DWECPP_P, DWECPT_P, &
           DWECTP, DWECTT, DCPPB, DCPTB, FCPPB
      real(rkind),dimension(:,:,:,:,:),POINTER :: &
           DCPP2,DCPT2,DCTP2,DCTT2,FCPP2,FCTH2  !(NTHM,NPM,NRM,NSAM,NSBM)
      real(rkind),dimension(:,:,:,:,:),POINTER :: &
           DCPP2B,DCPT2B,FCPP2B  !(NTHM,NPM,NRM,NSAM,NSBM)

      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSAM)
           RNSL,RJSL,RWSL,RPCSL,RPWSL,RPESL,RLHSL,RFWSL,RECSL,RWS123L, &
           RSPBL,RSPFL,RSPSL,RSPLL,RPDR,RNDR, RTL_BULK, RT_BULK, RICSL!, RNSL_test
      real(rkind),dimension(:,:,:),POINTER :: & ! (NRM,NSAM,NSBM)
           RPCS2L

      real(rkind),dimension(:,:,:),POINTER :: & ! (NRM,NSAM,NSBM)
           RPW_IMPL, RPWEC_IMPL, RPWIC_IMPL
      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSAM,NSBM)
           RPW_INIT, RPWEC_INIT, RPWIC_INIT

      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NSAM)
           RNS,RJS,RWS,RPCS,RPWS,RPES,RLHS,RFWS,RECS,RWS123, &
           RSPB,RSPF,RSPS,RSPL, RPDRL,RNDRL, RICS!, RNS_test
      real(rkind),dimension(:),POINTER :: & ! (NSAM)
           RNS_S2
      real(rkind),dimension(:,:,:),POINTER :: & ! (NRM,NSAM,NSBM)
           RPCS2

      real(rkind),dimension(:),POINTER :: & ! (NTG1M)
           PTG,PET,PQT
      real(rkind),dimension(:,:),POINTER :: & ! (NTG1M,NSAM)
           PNT,PWT,PTT,PIT,PPCT,PPWT,PPET,PLHT,PFWT,PECT,PTT3,PITT,PWTT,PICT
      real(rkind),dimension(:,:),POINTER :: & ! (NTG1M,NSAM)
           PSPT,PSPBT,PSPFT,PSPLT,PSPST
      real(rkind),dimension(:,:),POINTER :: & ! (NTG1M,NSAM)
           PNT2,PWT2,PTT2,PIT2,PWTD,PDR,PNDR,PTT_BULK
      real(rkind),dimension(:,:,:),POINTER :: & ! (NTG1M,NSAM,NSBM)
           PPCT2

      real(rkind),dimension(:),POINTER :: & ! (NTG2M)
           RTG
      real(rkind),dimension(:,:),POINTER :: & ! (NRM,NTG2M)
           RET,RQT
      real(rkind),dimension(:,:,:),POINTER :: & ! (NRM,NTG2M,NSAM)
           RNT,RWT,RTT,RJT,RPCT,RPWT,RPET,RLHT,RFWT,RECT, &
           RSPBT,RSPFT,RSPLT,RSPST,RPDRT,RNDRT,RTT_BULK,RICT!, RNT_test
      real(rkind),dimension(:,:,:,:),POINTER :: & ! (NRM,NTG2M,NSAM,NSBM)
           RPCT2

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
           FM,BMTOT

      real(rkind),dimension(:,:,:,:,:),POINTER :: & 
           SIGMAV_NF ! (NTHMAX+1,NPMAX+1,NTHMAX+1,NPMAX+1,6)
      real(rkind),dimension(:,:),POINTER :: & 
           RATE_NF ! (NRSTART:NREND,6)
      real(rkind),dimension(:,:,:,:),POINTER :: & 
           RATE_NF_D1, RATE_NF_D2 ! (NTHMAX,NPMAX,NRSTART:NREND,6)
      real(rkind),dimension(:,:),POINTER :: & ! (NPM:NSAM)
           PG2,PM2
      real(rkind),dimension(:),POINTER :: & ! (NSAM)
           DEPS_SS, RPDRS, RNDRS
      integer:: N_IMPL, NCALCNR

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

          allocate(MTXLEN(nprocs),MTXPOS(nprocs))

          allocate(F(NTHMAX,NPMAX,NRSTART:NREND))
          allocate(F1(NTHMAX,NPMAX,NRSTART:NREND))

          allocate(RG(NRMAX+1),RM(NRMAX+1),VOLR(NRMAX))
          allocate(RLAMDAG(NTHMAX,NRMAX+1),RLAMDA_GG(NTHMAX,NRMAX+1))
          allocate(ETAMG(NTHMAX,NRMAX+1),ETAM_GG(NTHMAX,NRMAX+1))
          allocate(ETAG_G_GL(NTHMAX+1,NRMAX+1))
          allocate(BP(NRMAX+1),QR(NRMAX))
          allocate(BPG(NRMAX+1),BPM(NRMAX+1))
          allocate(QLG(NRMAX+1),QLM(NRMAX+1))
          allocate(RJ1(NRMAX),E1(NRMAX))
          allocate(RJ2(NRMAX),E2(NRMAX))
          allocate(EPSRM(NRMAX+1),EPSRG(NRMAX+1))
          allocate(EPSRM2(NRMAX+1),EPSRG2(NRMAX+1))
          allocate(EPSRMX(NRMAX+1),EPSRGX(NRMAX+1))
          allocate(ITL(NRMAX+1),ITU(NRMAX+1))
          allocate(ITLG(NRMAX+1),ITUG(NRMAX+1))
          allocate(ITL_G(NRMAX+1),ITU_G(NRMAX+1))
          allocate(ITLG_G(NRMAX+1),ITUG_G(NRMAX+1))

          allocate(RCOEF(NRSTART:NREND), RCOEF_G(NRSTART:NREND))
          allocate(RCOEFN(NRSTART:NREND), RCOEFN_G(NRSTART:NREND))
          allocate(RCOEF1(NSAMAX),RCOEF2(NSAMAX))
          allocate(RCOEF2_G(NSAMAX))
          allocate(RCOEFG(NRMAX+1), RCOEF_GG(NRMAX+1))
          allocate(RCOEFNG(NRMAX+1), RCOEFN_GG(NRMAX+1))
          allocate(PG(NPMAX+1,NSBMAX),PM(NPMAX,NSBMAX))
          allocate(THG(NTHMAX+1),THM(NTHMAX))
          allocate(DELP(NSBMAX))
          allocate(PG2(NP2MAX+1,NSBMAX),PM2(NP2MAX,NSBMAX))

          allocate(VOLP(NTHMAX,NPMAX,NSBMAX))
          allocate(ETAG(NTHMAX+1,NRSTART:NREND+1),ETAM(NTHMAX,NRSTART:NREND+1))
          allocate(ETAG_G(NTHMAX+1,NRSTART:NREND+1),ETAM_G(NTHMAX,NRSTART:NREND+1))
          allocate(RLAMDA(NTHMAX,NRSTART:NREND),RLAMDC(NTHMAX+1,NRSTART:NREND))

          allocate(RFSAD(NRSTART:NREND),RFSADG(NRMAX))
          allocate(RFSAD_G(NRSTART:NREND),RFSAD_GG(NRMAX+1))

          allocate(RLAMDA_G(NTHMAX,NRSTART:NREND),RLAMDC_G(NTHMAX+1,NRSTART:NREND))
          allocate(SING(NTHMAX+1),COSG(NTHMAX+1))
          allocate(SINM(NTHMAX),COSM(NTHMAX))

          allocate(FNS(NTHMAX+1,NPMAX+1,NRMAX+1,NSBMAX))

          allocate(FNS1(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(FNS2(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(FNS22(NTHMAX+1,NPMAX+1,NRMAX+1,NSAMAX))
          allocate(FS1(NTHMAX+1,NPMAX+1,NSAMAX))
          allocate(FS2(NTHMAX+1,NPMAX+1,NSAMAX))
          allocate(FS3(NTHMAX+1,NPMAX+1,NSAMAX))

          allocate(RNFP0(NSAMAX),RNFPS(NSAMAX))
          allocate(RTFP0(NSAMAX),RTFPS(NSAMAX))
          allocate(AMFP(NSAMAX),AEFP(NSAMAX))
          allocate(PTFP0(NSAMAX),VTFP0(NSAMAX))
          allocate(RNFD0(NSBMAX))
          allocate(AEFD(NSBMAX),AMFD(NSBMAX))
          allocate(PTFD0(NSBMAX),VTFD0(NSBMAX))
          allocate(THETA0(NSBMAX))

          allocate(RNFP(NRSTART:NREND+1,NSAMAX),RTFP(NRSTART:NREND+1,NSAMAX))
          allocate(PTFP(NRSTART:NREND+1,NSAMAX),VTFP(NRSTART:NREND+1,NSAMAX))
          allocate(THETA(NRSTART:NREND,NSAMAX),DKBSR(NRSTART:NREND,NSAMAX))
          allocate(WEIGHP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(WEIGHT(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(WEIGHR(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))

          allocate(RNFD(NRSTART:NREND+1,NSBMAX),RTFD(NRSTART:NREND+1,NSBMAX))
          allocate(PTFD(NRSTART:NREND+1,NSBMAX),VTFD(NRSTART:NREND+1,NSBMAX))
          allocate(RNUF(NRSTART:NREND+1,NSBMAX,NSBMAX))
          allocate(RNUD(NRSTART:NREND+1,NSBMAX,NSBMAX))

          allocate(DPP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DPT(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DTP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DTT(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(FPP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(FTH(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DRR(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(FRR(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(SPP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(PPL(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))

          allocate(FEPP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(FETH(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))

          allocate(DCPP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DCPT(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DCTP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DCTT(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(FCPP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(FCTH(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))

          allocate(DWPP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWPT(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWTP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWTT(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWPP_P(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWPT_P(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWTP_P(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWTT_P(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))

          allocate(DWLHPP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWLHPT(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWFWPP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWFWPT(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))

          allocate(DWECPP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWECPT(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWECTP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWECTT(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWICPP(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWICPT(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWECPP_P(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWECPT_P(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWICPP_P(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DWICPT_P(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))

          allocate(SPPB(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(SPPF(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(SPPS(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DCPP2(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSBMAX,NSAMAX))
          allocate(DCPT2(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSBMAX,NSAMAX))
          allocate(DCTP2(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSBMAX,NSAMAX))
          allocate(DCTT2(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSBMAX,NSAMAX))
          allocate(FCPP2(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSBMAX,NSAMAX))
          allocate(FCTH2(NTHMAX+1,NPMAX+1,NRSTART:NREND+1,NSBMAX,NSAMAX))

          allocate(DCPPB(4,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(DCPTB(4,NPMAX+1,NRSTART:NREND+1,NSAMAX))
          allocate(FCPPB(4,NPMAX+1,NRSTART:NREND+1,NSAMAX))

          allocate(DCPP2B(4,NPMAX+1,NRSTART:NREND+1,NSBMAX,NSAMAX))
          allocate(DCPT2B(4,NPMAX+1,NRSTART:NREND+1,NSBMAX,NSAMAX))
          allocate(FCPP2B(4,NPMAX+1,NRSTART:NREND+1,NSBMAX,NSAMAX))

          
          allocate(RNSL(NRSTART:NRENDX,NSAMAX),RJSL(NRSTART:NRENDX,NSAMAX))
!          allocate(RNSL_test(NRSTART:NRENDX,NSAMAX))
          allocate(RWSL(NRSTART:NRENDX,NSAMAX),RWS123L(NRSTART:NRENDX,NSAMAX))
          allocate(RSPBL(NRSTART:NRENDX,NSAMAX),RSPFL(NRSTART:NRENDX,NSAMAX))
          allocate(RSPSL(NRSTART:NRENDX,NSAMAX),RSPLL(NRSTART:NRENDX,NSAMAX))
          allocate(RPCSL(NRSTART:NRENDX,NSAMAX),RPESL(NRSTART:NRENDX,NSAMAX))
          allocate(RPWSL(NRSTART:NRENDX,NSAMAX),RLHSL(NRSTART:NRENDX,NSAMAX))
          allocate(RFWSL(NRSTART:NRENDX,NSAMAX),RECSL(NRSTART:NRENDX,NSAMAX))
          allocate(RICSL(NRSTART:NRENDX,NSAMAX))
          allocate(RPCS2L(NRSTART:NRENDX,NSBMAX,NSAMAX))

          allocate(RNS(NRMAX,NSAMAX),RJS(NRMAX,NSAMAX))
!          allocate(RNS_test(NRMAX,NSAMAX))
          allocate(RNS_S2(NSAMAX))
          allocate(RWS(NRMAX,NSAMAX),RWS123(NRMAX,NSAMAX))
          allocate(RSPB(NRMAX,NSAMAX),RSPF(NRMAX,NSAMAX))
          allocate(RSPL(NRMAX,NSAMAX),RSPS(NRMAX,NSAMAX))
          allocate(RPCS(NRMAX,NSAMAX),RPES(NRMAX,NSAMAX))
          allocate(RPWS(NRMAX,NSAMAX),RLHS(NRMAX,NSAMAX))
          allocate(RFWS(NRMAX,NSAMAX),RECS(NRMAX,NSAMAX))
          allocate(RICS(NRMAX,NSAMAX))
          allocate(RPCS2(NRMAX,NSBMAX,NSAMAX))

          allocate(RPDR(NRMAX,NSAMAX),RNDR(NRMAX,NSAMAX))
          allocate(RPDRS(NSAMAX),RNDRS(NSAMAX))
          allocate(RPDRL(NRSTART:NRENDX,NSAMAX),RNDRL(NRSTART:NRENDX,NSAMAX))
          allocate(RT_BULK(NRMAX,NSAMAX))
          allocate(RTL_BULK(NRSTART:NRENDX,NSAMAX))

          allocate(RPW_IMPL(NRSTART:NRENDX,NSAMAX,0:LMAXFP+1))
          allocate(RPWEC_IMPL(NRSTART:NRENDX,NSAMAX,0:LMAXFP+1))
          allocate(RPWIC_IMPL(NRSTART:NRENDX,NSAMAX,0:LMAXFP+1))
          allocate(RPW_INIT(NRSTART:NRENDX,NSAMAX))
          allocate(RPWEC_INIT(NRSTART:NRENDX,NSAMAX))
          allocate(RPWIC_INIT(NRSTART:NRENDX,NSAMAX))

!         NLMAXM= 8   ! this is for analysis without bounce average
!         NLMAXM=11   ! this is for analysis without radial transport
          NLMAXM=15   ! this is for analysis with a simple radial transport
          NMMAX=NRMAX*NTHMAX*NPMAX

          allocate(NMA(NTHMAX,NPMAX,NRMAX))
          allocate(NLMAX(NMSTART:NMEND))
          allocate(LL(NMSTART:NMEND,NLMAXM))
          allocate(DL(NMSTART:NMEND),BM(NMSTART:NMEND))
          allocate(AL(NMSTART:NMEND,NLMAXM))
          allocate(FM(NMMAX),BMTOT(NMMAX))

          allocate(SIGMAV_NF(NTHMAX,NPMAX,NTHMAX,NPMAX,6))
          allocate(RATE_NF(NRSTART:NREND,6))
          allocate(RATE_NF_D1(NTHMAX,NPMAX,NRSTART:NREND,6))
          allocate(RATE_NF_D2(NTHMAX,NPMAX,NRSTART:NREND,6))

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

          deallocate(F)
          deallocate(F1)

          deallocate(FNS)
          deallocate(FNS1)
          deallocate(FNS2)
          deallocate(FNS22)

          deallocate(RCOEF,RCOEF_G)
          deallocate(RCOEF1,RCOEF2,VOLR,RCOEF2_G)
          deallocate(RCOEFG,RCOEF_GG)
          deallocate(ITL,ITU)
          deallocate(ITLG,ITUG)
          deallocate(ITL_G,ITU_G)
          deallocate(ITLG_G,ITUG_G)
          deallocate(RG,RM)
          deallocate(PG,PM)
          deallocate(THG,THM)

          deallocate(RNFP0,RNFPS)
          deallocate(RTFP0,RTFPS)
          deallocate(AMFP,AEFP)
          deallocate(PTFP0,VTFP0)
          deallocate(RNFD0)
          deallocate(AEFD,AMFD)
          deallocate(PTFD0,VTFD0)
          deallocate(THETA0,DELP)

          deallocate(RNFP,RTFP)
          deallocate(PTFP,VTFP)
          deallocate(BP,QR,RJ1,E1,RJ2,E2)
          deallocate(BPG,BPM)
          deallocate(QLG,QLM)
          deallocate(EPSRM,EPSRG,EPSRMX,EPSRGX)
          deallocate(EPSRM2,EPSRG2)
          deallocate(THETA,DKBSR)

          deallocate(RNFD,RTFD,PTFD,VTFD)
          deallocate(RNUF,RNUD)
          deallocate(FS1,FS2,VOLP,FS3)
          deallocate(ETAG,ETAM)
          deallocate(RLAMDA,RLAMDC)
          deallocate(RFSAD,RFSADG)
          deallocate(RFSAD_G,RFSAD_GG)
          deallocate(SING,COSG,SINM,COSM)
          deallocate(WEIGHP,WEIGHT)
          deallocate(WEIGHR)

          deallocate(DPP,DPT)
          deallocate(DTP,DTT)
          deallocate(FPP,FTH)
          deallocate(DRR,FRR)
          deallocate(SPP,PPL)
          deallocate(FEPP,FETH)

          deallocate(DCPP,DCPT)
          deallocate(DCTP,DCTT)
          deallocate(FCPP,FCTH)

          deallocate(DWPP,DWPT)
          deallocate(DWTP,DWTT)
          deallocate(DWPP,DWPT)
          deallocate(DWTP,DWTT)
          deallocate(DWPP_P,DWPT_P)
          deallocate(DWTP_P,DWTT_P)
          deallocate(DWPP_P,DWPT_P)
          deallocate(DWTP_P,DWTT_P)

          deallocate(DWLHPP,DWLHPT)
          deallocate(DWFWPP,DWFWPT)
          deallocate(DWECPP,DWECPT)
          deallocate(DWECTP,DWECTT)
          deallocate(DWICPP,DWICPT)
          deallocate(DWECPP_P,DWECPT_P)
          deallocate(DWICPP_P,DWICPT_P)
          deallocate(SPPB,SPPF,SPPS)

          deallocate(DCPP2,DCPT2)
          deallocate(DCTP2,DCTT2)
          deallocate(FCPP2,FCTH2)
          
          deallocate(RNSL,RJSL,RWSL)
!          deallocate(RNSL_test)
          deallocate(RPCSL,RPESL)
          deallocate(RPWSL,RLHSL,RFWSL,RECSL,RICSL)
          deallocate(RWS123L,RPCS2L)
          deallocate(RSPBL,RSPFL,RSPLL,RSPSL)
          deallocate(RPDRL,RNDRL)
          deallocate(RPDR,RNDR)
          deallocate(RPDRS,RNDRS)
          deallocate(RT_BULK,RTL_BULK)
          deallocate(RPW_IMPL, RPWEC_IMPL, RPWIC_IMPL)
          deallocate(RPW_INIT, RPWEC_INIT, RPWIC_INIT)

          deallocate(RNS,RJS,RWS)
!          deallocate(RNS_test)
          deallocate(RNS_S2)
          deallocate(RPCS,RPES)
          deallocate(RPWS,RLHS,RFWS,RECS,RICS)
          deallocate(RWS123,RPCS2)
          deallocate(RSPB,RSPF,RSPL,RSPS)

          deallocate(NMA)
          deallocate(NLMAX,LL)
          deallocate(DL,BM)
          deallocate(AL)
          deallocate(FM,BMTOT)

          deallocate(SIGMAV_NF,RATE_NF)
          deallocate(RATE_NF_D1,RATE_NF_D2)

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
          allocate(PNT(NSAMAX,NTG1M))
          allocate(PWT(NSAMAX,NTG1M))
          allocate(PTT(NSAMAX,NTG1M))
          allocate(PIT(NSAMAX,NTG1M))
          allocate(PPCT(NSAMAX,NTG1M))
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
          deallocate(PNT)
          deallocate(PWT)
          deallocate(PTT)
          deallocate(PIT)
          deallocate(PPCT)
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
                   DO NSA=1,NSAMAX
                      PNT(NSA,NTG)=PNT(NSA,2*NTG-1)
                      PWT(NSA,NTG)=PWT(NSA,2*NTG-1)
                      PWT2(NSA,NTG)=PWT2(NSA,2*NTG-1)
                      PTT(NSA,NTG)=PTT(NSA,2*NTG-1)
                      PIT(NSA,NTG)=PIT(NSA,2*NTG-1)
                      PPCT(NSA,NTG)=PPCT(NSA,2*NTG-1)
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
                deallocate(tempA)
                allocate(tempB(NSAMAX,NTG1M))
                call fp_adjust_ntg1_B(PNT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PWT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PWT2,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PTT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PIT,tempB,NTG1M_NEW)
                call fp_adjust_ntg1_B(PPCT,tempB,NTG1M_NEW)
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
          allocate(RNT(NRMAX,NSAMAX,NTG2M))
!          allocate(RNT_test(NRMAX,NSAMAX,NTG2M))
          allocate(RWT(NRMAX,NSAMAX,NTG2M))
          allocate(RTT(NRMAX,NSAMAX,NTG2M))
          allocate(RJT(NRMAX,NSAMAX,NTG2M))
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
          allocate(RPCT2(NRMAX,NSBMAX,NSAMAX,NTG2M))
          allocate(RPDRT(NRMAX,NSAMAX,NTG2M),RNDRT(NRMAX,NSAMAX,NTG2M))
          allocate(RTT_BULK(NRMAX,NSAMAX,NTG2M))

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
          deallocate(RNT)
!          deallocate(RNT_test)
          deallocate(RWT)
          deallocate(RTT)
          deallocate(RJT)
          deallocate(RPCT)
          deallocate(RPWT)
          deallocate(RPET)
          deallocate(RLHT)
          deallocate(RFWT)
          deallocate(RECT)
          deallocate(RICT)
          deallocate(RSPBT,RSPFT,RSPLT,RSPST)
          deallocate(RPDRT,RNDRT)
          deallocate(RPCT2)
          deallocate(RTT_BULK)

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
                   ENDDO
                   DO NSA=1,NSAMAX
                      DO NR=NRSTART,NREND
                         RNT(NR,NSA,NTG)=RNT(NR,NSA,2*NTG-1)
                         RJT(NR,NSA,NTG)=RJT(NR,NSA,2*NTG-1)
                         RWT(NR,NSA,NTG)=RWT(NR,NSA,2*NTG-1)
                         RTT(NR,NSA,NTG)=RTT(NR,NSA,2*NTG-1)
                         RPCT(NR,NSA,NTG)=RPCT(NR,NSA,2*NTG-1)
                         RPWT(NR,NSA,NTG)=RPWT(NR,NSA,2*NTG-1)
                         RPET(NR,NSA,NTG)=RPET(NR,NSA,2*NTG-1)
                         RLHT(NR,NSA,NTG)=RLHT(NR,NSA,2*NTG-1)
                         RFWT(NR,NSA,NTG)=RFWT(NR,NSA,2*NTG-1)
                         RECT(NR,NSA,NTG)=RECT(NR,NSA,2*NTG-1)
                         RICT(NR,NSA,NTG)=RICT(NR,NSA,2*NTG-1)
                         RSPBT(NR,NSA,NTG)=RSPBT(NR,NSA,2+NTG-1)
                         RSPFT(NR,NSA,NTG)=RSPFT(NR,NSA,2+NTG-1)
                         RSPLT(NR,NSA,NTG)=RSPLT(NR,NSA,2+NTG-1)
                         RSPST(NR,NSA,NTG)=RSPST(NR,NSA,2+NTG-1)
                         RPDRT(NR,NSA,NTG)=RPDRT(NR,NSA,2+NTG-1)
                         RNDRT(NR,NSA,NTG)=RNDRT(NR,NSA,2+NTG-1)
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
                deallocate(tempA)
                allocate(tempB(NRMAX,NSAMAX,NTG2M))
                call fp_adjust_ntg2_B(RNT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RWT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RTT,tempB,NTG2M_NEW)
                call fp_adjust_ntg2_B(RJT,tempB,NTG2M_NEW)
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
                deallocate(tempB)
                allocate(tempC(NRMAX,NSAMAX,NSBMAX,NTG2M))
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
                temp(NR,NSA,NSB,NTG)=data(NR,NSA,NSB,NTG)
             ENDDO
             ENDDO
             ENDDO
          ENDDO
          deallocate(data)
          allocate(data(NRMAX,NSAMAX,NSBMAX,NTG2M_NEW))
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
