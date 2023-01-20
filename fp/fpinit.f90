!     $Id: fpinit.f90,v 1.18 2013/01/22 16:21:46 fukuyama Exp $

      module fpinit

      PRIVATE
      PUBLIC fp_init

      contains

!     ***************************
!         INITIAL PARAMETERS
!     ***************************

      SUBROUTINE fp_init

      use fpcomm_parm
      IMPLICIT NONE
      integer:: ns,nsa,nsb,nbeam

!-----PARTICLE SPECIES--------------------------------------------------
!     nsamax: number of test particle species
!     nsbmax: number of field particle species (0 for nsbmax=nsmax: default)
!     ns_nsa(nsa): mapping from NSA to NS in pl 
!     ns_nsb(nsb): mapping from NSB to NS in pl
!     pmax(nsb)  : maximum momentum (normailzed by central thermal momentum)
!     Emax(nsb)  : if Emax is not 0, p upper boundary is defined from not pmax but sqrt(2*m*Emax)/PTFP0

      nsamax = 1
      nsbmax = 1
      DO nsa=1,nsm
         ns_nsa(nsa)=nsa          ! default test particle species list
      ENDDO
      DO nsb=1,nsm
         ns_nsb(nsb)=nsb          ! default field particle species lise
         pmax(nsb)=7.d0           ! default pmax=7
         pmax_bb(nsb)=5.d0        ! default pmax=5 R_beam_beam
         pmax_tt(nsb)=3.d0        ! default pmax=3 R_thermal_thermal
         Emax(nsb)=0.D0
      ENDDO

!-----RADIAL MESH---------------------------------------------------------
!     R1    : radial position for NRMAX=1 (r/a)
!     DELR1 : radial spacing for NRMAX=1 (r/a)
!     RMIN  : minimum radius for NRMAX<>1 (r/a)
!     RMIN  : maximum radius for NRMAX<>1 (r/a)

      R1    = 0.1D0*RR
      DELR1 = 0.1D0
      RMIN  = 0.05D0
      RMAX  = 0.3D0

!-----E0----------------------------------------------------------------
!     E0    : toroidal electric field (V/m)
!     ZEFF  : effective ion charge for simple collision term

      E0    = 0.D0
      ZEFF  = 1.D0

!-----WM/WR-------------------------------------------------------------
!     PABS_LH : input power of LH [MW] (0 for given DLH) MODELW=0
!     PABS_FW : input power of FW [MW] (0 for given DFW) MODELW=0
!     PABS_EC : input power of EC [MW] (0 for given DEC) MODELW=0
!     PABS_WR : input power of WR [MW] (0 for given E) MODELW=1 or 2
!     PABS_WM : input power of WM [MW] (0 for given E) MODELW=3 or 4
!     RF_WM   : wave frequency [MHz] used for MODELW=3

!     FACT_WM: Numerical factor for wave amplitude for WR
!     FACT_WR: Numerical factor for wave amplitude for WM
!     DELNPR_WR: width of toroidal mode number for WR
!     DELNPR_WM: width of toroidal mode number for WM
!     LMAX_WR: max loop count in newton method to find ray position
!     EPS_WR: convergence criterion in newton method to find ray position
!     DELY_WR: vertical half-width of ray [r/a]
!     Y0_WM: vertical position of wave beam [r/a]
!     DELY_WM: vertical half-width of wave beam [r/a]
!     NCMIN(NS): minimum order of cyclotron harmonics for species NS
!     NCMAX(NS): maximum order of cyclotron harmonics for species NS

      PABS_LH = 0.0D0
      PABS_FW = 0.0D0
      PABS_EC = 0.0D0
      PABS_WR = 0.0D0
      PABS_WM = 0.0D0
      RF_WM   = 64.D0

      FACT_WR  = 1.D0
      FACT_WM  = 1.D0
      DELNPR_WR= 0.05D0
      DELNPR_WM= 0.05D0
      LMAX_WR = 100
      EPS_WR  = 1.D-6
      DELY_WR = 0.1D0
      Y0_WM   = 0.D0
      DELY_WM =10.D0
      DO NS=1,NSM
         NCMIN(NS) = -3
         NCMAX(NS) = 3
      ENDDO

!-----EC/LH/FW----------------------------------------------------------
!     DEC   : Electron Cyclotron wave diffusion coefficient (normalized)
!     PEC1  : EC wave spectrum N para center
!     PEC2  : EC wave spectrum N para width
!     PEC3  : EC wave horizontal minimum x [m]
!     PEC4  : EC wave horizontal decay length [m]
!     RFEC  : Electron cyclotron frequency at r = 0 / EC wave frequency
!     DELYEC: EC wave beam vertical  width (m)
!     DLH   : Lower Hybrid wave diffusion coefficient (normalized)
!     PLH1  : LH wave spectrum (velocity minimum/N para center)
!     PLH2  : LH wave spectrum (velocity maximum/N para width)
!     RLH   : LH wave minimum accessible minor radius (m)
!     DFW   : Fast Wave diffusion coefficient (normalized)
!     PFW1  : FW spectrum (velocity minimum/N para center)
!     PFW2  : FW spectrum (velocity maximum/N para width)
!     RFW   : FW mimum accessible minor radius (m)

      DEC   = 0.D0
      PEC1  = 0.5D0
      PEC2  = 0.1D0
      PEC3  =-3.D0
      PEC4  = 0.1D0
      RFEC  = 0.8D0
      DELYEC= 0.1D0
      DLH   = 0.D0
!      PLH1  = 4.D0
!      PLH2  = 7.D0
      PLH1  = 1.92D0
      PLH2  = 0.3D0
      RLH   = 0.D0
      DFW   = 0.D0
!      PFW1  = 4.D0
!      PFW2  = 7.D0
      PFW1  = 0.D0
      PFW2  = 1.1D0
      RFW   = 0.D0

!-----SIMPLE RF----------------------------------------------------------------
!     CEWR  : radial component of wave electric field [V/m]
!     CEWTH : poloidal component of wave electric field [V/m]
!     CEWPH : toroidal component of wave electric field [V/m]
!     RKWR  : radial component of wave number [1/m]
!     RKWTH : poloidal component of wave number [1/m]
!     RKWPH : toroidal component of wave number [1/m]

      CEWR  = (0.D0,0.D0)
      CEWTH = (0.D0,0.D0)
      CEWPH = (1.D3,0.D0)
      RKWR  = 0.D0
      RKWTH = 0.D0
      RKWPH = 1.D3

!-----NBI---------------------------------------------------------------
!     NBEAMMAX     : Number of NBI
!     NSSPB(nbeam) : NBI particle species
!     SPBTOT(nbeam): Particle source [1/s] 
!     SPBR0(nbeam) : Source radius [r/a]
!     SPBRW(nbeam) : Source width [r/a]
!     SPBENG(nbeam): Particle energy [eV]
!     SPBPANG(nbeam): Source poloidal angle [degree]
!     SPBANG(nbeam): Source pitch angle [degree] ! at chi=SPBPANG
!     for f1_1.dat
!     MODEL_CD     :=1 calculate time evolution of j including inductive current
!     MODEL_NBCD   : MODEL of electron cancel current
!                   =0: no cancel current
!                   =1: Simplest Ohkawa current

      DO NBEAM=1,NBEAMM
         NSSPB(NBEAM)=2
         SPBTOT(NBEAM)=0.d0
         SPBR0(NBEAM)=0.d0
         SPBRW(NBEAM)=0.2d0
         SPBENG(NBEAM)=1.D6
         SPBANG(NBEAM)=20.D0
         SPBPANG(NBEAM)=0.D0
      ENDDO
      NSA_F1=1
      NTH_F1=1
      NR_F1=1
      MODEL_NBCD=0
      MODEL_CD=0
!-----FUSION REACTION----------------------------------------------------
!     NSSPF  : Fusion product particle species
!     SPFTOT : Particle source [1/m^3 s]
!     SPFR0  :  Source radius [r/a]
!     SPFRW  :  Source width [r/a]
!     SPFENG : Particle energy [eV]

      NSSPF=4
      SPFTOT=0.d0
      SPFR0=0.d0
      SPFRW=0.2d0
      SPFENG=3.5D6

!-----CHARGE EXCHANGE REACTION-------------------------------------------
!     MODEL_CX_LOSS: 0=off, 1=on 
!     RN_NEU0      : neutral gas density on axis [m^-3] D or H gas is assumed temporally
!     RN_NEUS      : neutral gas density on edge [m^-3] D or H gas is assumed temporally
!                    [1.e20]
      MODEL_CX_LOSS=0
      RN_NEU0 = 1.D-6
      RN_NEUS = 1.D-5
!-----RADIAL DIFFUSION--------------------------------------------------
!     DRR0  : radial diffusion coefficient at magnetic axis [m^2/s]
!     DRRS  : radial diffusion coefficient at plasma surface [m^2/s]
!     FACTOR_CDBM : Multiplication factor for CDBM model
!     DRR_EDGE : DRR at plasma edge [m^2/s]
!     RHO_EDGE : Normalized radius of inner boundary of edge region
!     FACTOR_DRR_EDGE : Reduction factor of DRR at plasma edge
!     FACTOR_PINCH : Pinch factor for exp(-factor r^2/a^2)

      DRR0       = 0.D0
      DRRS       = 0.D0
      DELTAB_B   = 0.D0 
      FACTOR_CDBM= 0.D0
      DRR_EDGE   = 0.1D0
      RHO_EDGE   = 0.95D0
      FACTOR_DRR_EDGE=0.1D0
      FACTOR_PINCH=1.0D0

!-----LOSS--------------------------------------------------------------
!     TLOSS(ns): loss time [s] (0.D0 for no loss)

      DO nsa=1,nsm
         TLOSS(nsa)=0.d0          ! default no loss
         TLOSS_PARA(nsa)=0.d0          ! default no loss
         TLOSS_PERP(nsa)=0.d0          ! default no loss
      ENDDO

!-----NUMBER OF MESH----------------------------------------------------
!     NPMAX : momentum magnitude division number
!     NTHMAX: momentum angle division number
!     NRMAX : radial division number
!     NAVMAX: wave diuffusion bounce average division number
!     NP2MAX: minor mesh number for electron momentum
!     FACT_BULK: Definition of the bulk region, 0<p<FACT_BULK*p_th is the bulk region

      NPMAX = 50
      NTHMAX= 50
      NRMAX = 1
      NAVMAX= 100
      NP2MAX= 20
      FACT_BULK = 5.D0

!-----NUMBER OF TIME STEP-----------------------------------------------
!     NTMAX   : maximum time step count
!     NTSTEP_COEF: time step interval for recalculating transport coefficients
!     NTSTEP_COLL: time step interval for recalculating collision coefficients
!     NTG1STEP: time step interval for storing global data
!     NTG1MIN:  minimum number of NTG1 save (initial allocation)
!     NTG1MAX:  maximum number of NTG1 save (thin out if exceeded)
!     NTG2STEP: time step interval for storing radial profile  data
!     NTG2MIN:  minimum number of NTG2 save (initial allocation)
!     NTG2MAX:  maximum number of NTG2 save (thin out if exceeded)

      NTMAX    = 10
      NTSTEP_COEF = 1000
      NTSTEP_COLL = 1000
      NTG1STEP = 1
      NTG1MIN  = 101
      NTG1MAX  = 1001
      NTG2STEP = 1
      NTG2MIN  = 101
      NTG2MAX  = 1001

!-----MODEL SELECTION---------------------------------------------------
!     MODELE: 0 for fixed electric field
!             1 for predicting electric field
!     MODELA: 0 without bounce average
!             1 with bounce average
!     MODELC: 0 : non-relative background Maxwell
!             1 : isotropic background f
!             2 : isotropic background f, Temperature is updated
!             4 : nonlinear collision operator (require to satisfy NSAMAX=NSBMAX)
!             5 : linear coll. operator for different species (for debug)
!             6 : nonlinear coll. operator for different species (for debug)
!            -1 : linear collision operator for same with ion scattering
!            -2 : nonlinear collision operator for same with ion scattering
!     MODELR: 0 : without relativistic effect
!             1 : with relativistic effect
!     MODELS : 0 No fusion reaction
!              1 Constant isotropic fast ion source (For example fast alpha)
!              2 Self-consistent fusion reaction source and loss term
!              3 Using Legendre expansion in fusion reaction rate calculation
!             -2 Self-consistent fusion reaction rate (no source and loss term)
!                Only Fusion reaction rate is calculated.
!     MODELS_* : assumed condition MODELS=2 or -2
!              MODELS_full : 0 double integration is omitted
!                            1 double integration is calculated (heavy)
!              MODELS_bt   : 0 beam thermal fusion reaction rate is not calcularted
!                            1 Mikkelsen models is used for beam-thermal fusion reaction.
!                            2 beam-thermal fusion reaction is evaluated by pmax_bb.
!              MODELS_tt :   0 thermal-thermal component is not calculated.
!                            1 thermal-thermal component is calculated.
!     MODELW(ns): 0 for given diffusion coefficient model
!                 1 for wave E field calculated by WR(without beam radius)
!                 2 for wave E field calculated by WR(with beam radius)
!                 3 for given wave E field model
!                 4 for wave E field calculated by WM
!     MODELD: 0 : without radial transport
!             1 : with radial transport
!     MODELD_RDEP : 0 : fixed:    (DRR0-DRRS)*(1.D0-RHON**2)+DRRS 
!                   1 : magnetic: QLM(NR)*deltaB_B**2 
!                   2 : CDBM with FACTOR_CDBM
!     MODELD_PDEP : 0 : no p dependence
!                   1 : 1/sqrt(p_perp) dependence
!                   2 : 1/p_perp dependence
!                   3 : 1/p_perp^2 dependence
!                   4 : v_para
!     MODELD_EDGE : 0 : as it is near edge
!                   1 : reduced to DRR_edge for rho > rho_edge
!                   2 : reduced by FACTOR_DRR_edge for rho > rho_edge
!     MODELD_PINCH  0 : no pinch
!                   1 : no radial particle transport (particle flux = 0)
!                   2 : v = - 2 * FACTOR_PINCH * r * DRR / RA**2
!     MODELD_BOUNDARY : 0 fix f at rho=1+DELR/2, namely FS2
!                     : 1 fix f at rho = 1, namely FS1. FS2 is variable.
!     MODELD_CDBM     :  0: CDBM original
!                        1: CDBM05 including elongation
!                        2: CDBM original with weak ExB shear
!                        3: CDBM05 with weak ExB shear
!                        4: CDBM original with strong ExB shear
!                        5: CDBM05 with strong ExB shear
!     MODEL_LOSS      : 1 for LOSS TERM
!     MODEL_SYNCH     : 1 for synchlotron radiation
!     MODEL_NBI       : 1 for NBI calculation with simple model
!                     : 2 read FIT3D data (limited)
!     MODEL_WAVE      : 1 for wave calculation
!     MODEL_BULK_CONST: 0 ordinary
!                     : 1 If MODEL_DELTA_F=0, FNSP in bulk region is Maxwellian
!                         IF MODEL_DELTA_F=1, FNSP_DEL in bulk region replaced by FNSP_DEL=0 
!                     : 2 Strong artifical sink term reduce FNSP_DEL in bulk (DELTA_F ONLY)
!     MODEL_NF_CS     : 0 NF cross secion in NRL with relative energy
!                     : 1 NRL 
!                     : 2 Bosch 
!     MODEL_DELTA_F(NSA): 0 ordinary full f mode
!                       : 1 f is described as f=f_M + delta f
!                         Evolution of f_M is not solved.
!          MODEL_DELTA_F_NI_RATIO
!                              Requires MODEL_DELTA_F(NSA) = 1
!                     : 0 for prediction
!                     : 1 USE NI_RATIO. 
!                         NI_RATIO, which means bulk ion density ratios to the electron density, is const. 
!                         This option leads over-estimation of n_i. (n_e < n_i_bulk + n_i_EP)
!                     : 2 USE charge neutrality, zeff, DH_RATIO, HHe_RATIO, and DT_RATIO. 
!          MODEL_DELTA_F_CN    Does charge neutrarility consider beam density?
!                              Requires MODEL_DELTA_F_NI_RATIO = 2
!                     : 0 Does not satisfy charge neutrality : n_e = n_i^bulk, n_i = n_i^bulk + n_i^beam
!                     : 1 Satisfies charge neutrality        : n_e = n_i,      n_i = n_i^bulk + n_i^beam
!     NI_RATIO(NS)    : Density ratio for each species. Electron is unity. MODEL_EX_READ_Tn or MODEL_DELTA_F
!     MODEL_BULK_T(NS): How to calculate the bulk temperature
!                     : 0 calculate from f (need MODEL_DELTA_F(NSA)=0). 
!                       IF MODEL_DELTA_F(NSA)=0 & MODEL_BULK_T(NS)=0, MODEL_BULK_T(NS) is replaced by 1.
!                     : 1 keep initial value
!                     : 2 read external file (Require MODEL_EX_READ_Tn!=0)

      MODELE= 0
      MODELA= 0
      MODELC= 0
      MODELR= 0
      MODELS= 0
      DO NSA=1,NSM
         MODELW(NSA)=0
      END DO
      MODELS_full=1
      MODELS_bt=0
      MODELS_tt=1

      MODELD= 0
      MODELD_RDEP = 0
      MODELD_PDEP = 0
      MODELD_EDGE = 0
      MODELD_PINCH = 0
      MODELD_BOUNDARY= 0
      MODELD_CDBM= 0

      MODEL_LOSS=0
      MODEL_SYNCH=0
      MODEL_NBI=0
      MODEL_WAVE=0 ! 0=no wave calc., 1=wave calc.

      MODEL_BULK_CONST=0

      DO NS=1,NSM
         MODEL_DELTA_F(NS)=0
      END DO
      MODEL_NF_CS = 2
      NF_IDMAX = 6
      OUTPUT_NFID=1

      MODEL_DELTA_F_NI_RATIO=2
      MODEL_DELTA_F_CN=1
      DO NS=1,NSM
         NI_RATIO(NS)=1.D0
      END DO
      DH_RATIO=1.D0 ! pure D, nD/(nD+nH)
      HHe_RATIO=0.D0 ! no He, nH/(nH+nHe)
      DT_RATIO=1.D0 ! no T, nD/(nD+nT)
      given_zeff=1.D0 ! no impurity
      DO NS=1,NSM
         MODEL_BULK_T(NS)=0
      END DO

!-----READ EXP DATA------------------------------------------
!     MODEL_EX_READ_Tn: 0 Read no temperature and density
!                     : 1 Read electron temperature and density (no ion data case)
!                     : 2 Read electron and ion temperature and density 
!     time_exp_offset : TIMEFP + time_exp_offset = time in experiment
!     Ti_Te_ratio(NSMAX): If you assume T_e!=T_i, T_i=(Ti_Te_ratio(NS))*T_e
!                         premise MODEL_EX_READ_Tn=1
!
!     te_poly(1:4), ne_poly(1:5): LHD like profile fitting coef 
!     IF MODEL_PROF_POLY=1, these are used. IF MODEL_EX_READ_Tn!=0, MODEL_EX_READ_Tn has higher priority.
!     MODEL_Q_PROF      : 0 Default. q profile is defined according to MODELG
!                     : 1 Read VMEC output (kspdiag_data)

      MODEL_EX_READ_Tn=0
      time_exp_offset=0.D0
      DO NS=1,NSM
         Ti_Te_ratio(NS)=1.D0
      END DO

      MODEL_PROF_POLY=0
      te_poly(:)=0.D0
      te_poly(1)=3.D0
      ne_poly(:)=0.D0
      ne_poly(1)=3.D0
      MODEL_Q_PROF=0

!-----TXT TYPE OUTPUT COMMAND------------------------------------------
!     OUTPUT_TXT_DELTA_F: OUTPUT DELTA f
!     OUTPUT_TXT_F1: OUTPUT E-f1 on NR_F1 DIRECT TO NTH_F1

      OUTPUT_TXT_DELTA_F=0
      OUTPUT_TXT_F1=0
      OUTPUT_TXT_BEAM_WIDTH=0
      OUTPUT_TXT_HEAT_PROF=0
      OUTPUT_TXT_BEAM_DENS=0
      OUTPUT_BEAM_BIRTH_PROF=0

!-----COMPUTATION PARAMETERS------------------------------------------
!     DELT  : time step size (s)
!     RIMPL : implicit computation parameter
!     IMTX  : info level of KSP solver
!     MODEL_KSP : type of KSP solver
!     MODEL_PC  : type of pre-conditioning
!     EPSFP : convergence limit in nonlinear effects
!     LMAXFP: maximum loop count in nonlinear effects
!     EPSM  : convergence limit in matrix equation solver
!     EPSE  : convergence limit in electric field prediction
!     LMAXE : maximum loop count in electric field prediction
!     EPSDE : convergence limit in double-exponential integration
!     H0DE  : initial step size in double-exponential integration

      DELT  = 1.D-2
      RIMPL = 1.D0
      IMTX  = 1
      MODEL_KSP=5
      MODEL_PC =1
      EPSFP = 1.D-9
      LMAXFP = 10
      EPSM  = 1.D-8
      EPSE  = 1.D-4
      LMAXE = 10
      EPSDE = 1.D-8
      H0DE  = 0.25D0

!-----GRAPHI PARAMETERS------------------------------------------
!     PGMAX:  maximum p in graphics (if not 0)
!     RGMIN:  minimum rho in graphics (if not 0)
!     RGMAX:  maximum rho in graphics (if not 1)
!     NGLINE: maximum number of contour lines
!     NGRAPH: graphic mode: 0 for file out, 1 for contour plot

      PGMAX=10.D0
      RGMIN=0.D0
      RGMAX=1.D0
      NGLINE= 25
      NGRAPH= 1

!-----LEGENDRE EXPANSION-----------------------------------------------
!     LLMAX    : dimension of legendre polynomials's calculation
!     LLMAX_NF : dimension of legendre polynomials's calculation for NF rate

      LLMAX = 2
      LLMAX_NF = 2

!-----------------------------------------------------------------------
!     IDBGFP : debug graphic control parameter 
!          1 : Legendre polynomials
!          2 : fpl, M_l, N_l
!          4 : psy, phy and their derivatives
!          8 : dcpp, dctt, fcp

      IDBGFP=0

!-----------------------------------------------------------------------
!     Parameters relevant to DISRUPTION 
!
!     T0_quench        : temperature after thermal quench [keV] at r=0
!     tau_quench       : thermal quench time [sec]
!     tau_mgi          : MGI duration [sec]
!     MODEL_DISRUPT    : 0=no disruption, 1=disruption calc.
!     MODEL_Connor_FP  : runaway rate 0= Connor, 1=FP
!     MODEL_BS         : bootstrap current 0= off, 1=simple model
!     MODEL_jfp        : current evaluation 
!                            0: independent of f
!                            1: depend on f
!     MODEL_LNL        : Coulomb logarithm 
!                            0: variable w T 
!                            1: fixed initial value
!                            2: fixed disrupted value
!     MODEL_RE_pmax    : RE non-RE boundary 0=NPMAX, 1=NPC_runaway
!     MODELD_n_RE      : radial transport of RE density 0=off, 1=on
!     MODEL_IMPURITY   :     0: Default
!                            1: MGI: satisfy quasi-neutrality
!     MODEL_SINK       :     0: Default
!                            1: deltaB/B like sink term
!     time_quench_start: ?
!     RJPROF1          : ?
!     RJPROF2          : ?
!     v_RE             : RE velocity / VC
!     target_zeff      : ?
!     N_IMPU           : ?
!     SPITOT           : ?   m^-3 s^-1 on axis
      
      T0_quench=2.D-2
      tau_quench=1.D-3
      tau_mgi=5.D-3
      MODEL_DISRUPT=0
      MODEL_Connor_FP=0
      MODEL_BS=0
      MODEL_jfp=0
      MODEL_LNL=0
      MODEL_RE_pmax=0
      MODELD_n_RE=0
      MODEL_IMPURITY=0
      MODEL_SINK=0

      time_quench_start=0.D0
      RJPROF1=2.D0
      RJPROF2=2.D0
      v_RE=1.D0
      target_zeff=3.D0
      N_IMPU=3
      SPITOT=0.D0 ! m^-3 s^-1 on axis

!-----------------------------------------------------------------------
!     MPI Partition number : 
!          N_partition_s: the number of patition for NSA
!          N_partition_r: the number of patition for NR
!          N_partition_p: the number of patition for NP

      N_partition_s = 1
      N_partition_r = nsize
      N_partition_p = 1

      RETURN
      END SUBROUTINE fp_init

  END module fpinit
