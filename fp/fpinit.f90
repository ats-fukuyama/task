! fpinit.f90

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
      integer:: ns,nsa,nsb,nbeam,nray

!-----PARTICLE SPECIES--------------------------------------------------
!     nsamax: number of test particle species
!     nsbmax: number of field particle species
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
         Emax(nsb)=0.D0
      ENDDO

!-----RADIAL MESH---------------------------------------------------------
!     R1    : radial position for NRMAX=1 (r/a)
!     DELR1 : radial spacing for NRMAX=1 (r/a)
!     RMIN  : minimum radius for NRMAX<>1 (r/a)
!     RMIN  : maximum radius for NRMAX<>1 (r/a)

      R1    = 0.3D0
      DELR1 = 0.1D0
      RMIN  = 0.05D0
      RMAX  = 0.3D0

!-----E0----------------------------------------------------------------
!     E0    : toroidal electric field (V/m)
!     ZEFF  : effective ion charge for simple collision term

      E0    = 0.D0
      ZEFF  = 1.D0

!-----WM/WR-------------------------------------------------------------
!     PABS_LH : absorbed power of LH [MW] (0 for given DLH) MODELW=0
!     PABS_FW : absorbed power of FW [MW] (0 for given DFW) MODELW=0
!     PABS_EC : absorbed power of EC [MW] (0 for given DEC) MODELW=0
!     PIN_WR  : input power of WR [MW] (0 for given E) MODELW=1 or 2
!     PABS_WM : absorbed power of WM [MW] (0 for given E) MODELW=3 or 4
!     RF_WM   : wave frequency [MHz] used for MODELW=3

!     FACT_WM: Numerical factor for wave amplitude for WM
!     FACT_WR: Numerical factor for wave amplitude for WR
!     PIN_WR_NRAY(NRAYM): input power of each ray (for PIN_WR=0)
!     DELNPR_WR: width of toroidal mode number for WR
!     DELNPR_WM: width of toroidal mode number for WM
!     LMAX_WR: max loop count in newton method to find ray position
!     EPS_WR: convergence criterion in newton method to find ray position
!     DELY_WR: vertical half-width of ray [m]
!     Y0_WM: vertical position of wave beam [r/a]
!     DELY_WM: vertical half-width of wave beam [r/a]
!     NRAYS_WR: start of NRAY (IF 0, NRAYS_WR=1)
!     NRAYE_WR: end of NRAY (IF 0, NRAYE_WR=NRAYMAX)
!     NCMIN(NS): minimum order of cyclotron harmonics for species NS
!     NCMAX(NS): maximum order of cyclotron harmonics for species NS

      PABS_LH = 0.0D0
      PABS_FW = 0.0D0
      PABS_EC = 0.0D0
      PIN_WR  = 0.0D0
      PABS_WM = 0.0D0
      RF_WM   = 64.D0

      FACT_WR  = 1.D0
      FACT_WM  = 1.D0
      DO NRAY=1,NRAYM
         PIN_WR_NRAY(NRAY)=1.D0
      END DO
      DELNPR_WR= 0.05D0
      DELNPR_WM= 0.05D0
      LMAX_WR = 100
      EPS_WR  = 1.D-6
      DELY_WR = 0.1D0
      Y0_WM   = 0.D0
      DELY_WM =0.1D0
      NRAYS_WR=0
      NRAYE_WR=0
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

      NBEAMMAX=0
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
      FACTOR_CDBM= 1.D0
      DRR_EDGE   = 0.1D0
      RHO_EDGE   = 0.95D0
      FACTOR_DRR_EDGE=0.1D0
      FACTOR_PINCH=1.0D0

!-----LOSS--------------------------------------------------------------
!     TLOSS(ns): loss time [s] (0.D0 for no loss)

      DO nsa=1,nsm
         TLOSS(nsa)=0.d0          ! default no loss
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
!     MODELC(ns): 0 : non-relativistic background Maxwell
!                 1 : isotropic background f
!                 2 : isotropic background f, Temperature is updated
!                 4 : nonlinear collision operator
!                 5 : linear col. operator for different species (for debug)
!                 6 : nonlinear col. operator for different species (for debug)
!                -1 : linear collision operator for same with ion scattering
!                -2 : nonlinear collision operator for same with ion scattering
!     MODELR: 0 : without relativistic effect
!             1 : with relativistic effect
!     MODELS: 0 No fusion reaction
!             1 DT reaction source (NSSPF,SPFTOT,SPFR0,SPFRW,SPFENG)
!             2 DT reaction source (self-consistent reactioin rate)
!             3 Using Legendre expansion in fusion reaction rate calculation
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
!     MODEL_EX_READ_Tn: 0 Read temperature and density
!                     : 1 Read electron temperature and density (no ion data case)
!                     : 2 Read electron and ion temperature and density 
!     MODEL_EX_READ_DH_RATIO
!                     : 0 for prediction
!                     : 1 constant ratio. Use NI_RATIO
!                     : 2 read egdata 
!     time_exp_offset : TIMEFP + time_exp_offset = time in experiment
!
!     MODEL_DELTA_F(NSA): 0 ordinary mode
!                       : 1 f is described as f=f_M + delta f
!                         Evolution of f_M is not solved.

      MODELE= 0
      MODELA= 0
      MODELR= 0
      MODELS= 0
      DO NSA=1,NSM
         MODELW(NSA)=0
      END DO

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

      MODEL_EX_READ_Tn=0
      MODEL_EX_READ_DH_RATIO=0
      DO NS=1,NSM
         NI_RATIO(NS)=1.D0
      END DO

      MODEL_BULK_CONST=0
      time_exp_offset=0.D0

      DO NS=1,NSM
         MODEL_DELTA_F(NS)=0
      END DO
!-----TXT TYPE OUTPUT COMMAND------------------------------------------
!     OUTPUT_TXT_DELTA_F: OUTPUT DELTA f
!     OUTPUT_TXT_F1: OUTPUT E-f1 on NR_F1 DIRECT TO NTH_F1

      OUTPUT_TXT_DELTA_F=0
      OUTPUT_TXT_F1=0
      OUTPUT_TXT_BEAM_WIDTH=0
      OUTPUT_TXT_HEAT_PROF=0
      OUTPUT_TXT_BEAM_DENS=0

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
