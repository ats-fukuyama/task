!     $Id: fpinit.f90,v 1.18 2013/01/22 16:21:46 fukuyama Exp $

      module fpinit

      IMPLICIT NONE

      contains

!     ***************************
!         INITIAL PARAMETERS
!     ***************************

      SUBROUTINE fp_init

      use fpcomm_parm
      integer:: ns,nsa,nsb,nbeam

!-----PARTICLE SPECIES--------------------------------------------------
!     nsamax: number of test particle species
!     nsbmax: number of field particle species (0 for nsbmax=nsmax: default)
!     ns_nsa(nsa): mapping from NSA to NS in pl 
!     ns_nsb(nsb): mapping from NSB to NS in pl
!     pmax(nsb)  : maximum momentum (normailzed by central thermal momentum)

      nsamax = 1
      nsbmax = 1
      DO nsa=1,nsm
         ns_nsa(nsa)=nsa          ! default test particle species list
      ENDDO
      DO nsb=1,nsm
         ns_nsb(nsb)=nsb          ! default field particle species lise
         pmax(nsb)=7.d0           ! default pmax=7
         pmax_bb(nsb)=5.d0        ! default pmax=5 R_beam_beam
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
!     PWAVE : input power
!     RFDW  : wave frequency [MHz]
!     DELNPR: width of toroidal mode number
!     LMAXNWR: max loop count in newton method to find ray position
!     EPSNWR: convergence criterion in newton method to find ray position
!     REWY  : vertical position of ray [r/a]
!     DREWY : vertical half-width of ray [r/a]
!     FACTWM: Numerical factor for wave amplitude
!     NCMIN(NSA): minimum order of cyclotron harmonics
!     NCMAX(NSA): maximum order of cyclotron harmonics

      PWAVE = 1.0D0
      RFDW  = 170.D3
      DELNPR= 0.05D0
      LMAXNWR=100
      EPSNWR=1.D-6
      REWY  = 0.D0
      DREWY = 0.1D0
      FACTWM= 1.D0
      DO NSA=1,NSM
         NCMIN(NSA) = -3
         NCMAX(NSA) = 3
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
!     SPBTOT(nbeam): Particle source [1/m^3 s] ! [1/s]? 
!     SPBR0(nbeam) : Source radius [r/a]
!     SPBRW(nbeam) : Source width [r/a]
!     SPBENG(nbeam): Particle energy [eV]
!     SPBPANG(nbeam): Source poloidal angle [degree]
!     SPBANG(nbeam): Source pitch angle [degree] ! at phi=SPBPANG

      DO NBEAM=1,NBEAMM
         NSSPB(NBEAM)=2
         SPBTOT(NBEAM)=0.d0
         SPBR0(NBEAM)=0.d0
         SPBRW(NBEAM)=0.2d0
         SPBENG(NBEAM)=1.D6
         SPBANG(NBEAM)=20.D0
         SPBPANG(NBEAM)=20.D0
      ENDDO

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
      ENDDO

!-----NUMBER OF MESH----------------------------------------------------
!     NPMAX : momentum magnitude division number
!     NTHMAX: momentum angle division number
!     NRMAX : radial division number
!     NAVMAX: wave diuffusion bounce average division number
!     NP2MAX: minor mesh number for electron momentum

      NPMAX = 50
      NTHMAX= 50
      NRMAX = 1
      NAVMAX= 100
      NP2MAX= 20

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
!             4 : nonlinear collision operator
!             5 : linear coll. operator for different species (for debug)
!             6 : nonlinear coll. operator for different species (for debug)
!            -1 : linear collision operator for same with ion scattering
!            -2 : nonlinear collision operator for same with ion scattering
!     MODELR: 0 : without relativistic effect
!             1 : with relativistic effect
!     MODELS : 0 No fusion reaction
!              1 DT reaction source (NSSPF,SPFTOT,SPFR0,SPFRW,SPFENG)
!              2 DT reaction source (self-consistent reactioin rate)
!              3 Using Legendre expansion in fusion reaction rate calculation
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
!                   2 : v = - 2 * FACTOR_PINCH * r * DRR
!     MODELD_BOUNDARY : 0 fix f at rho=1+DELR/2, namely FS2
!                     : 1 fix f at rho = 1, namely FS1. FS2 is variable.
!     MODEL_LOSS      : 1 for LOSS TERM
!     MODEL_SYNCH     : 1 for synchlotron radiation
!     MODEL_NBI       : 1 for NBI calculation
!     MODEL_WAVE      : 1 for wave calculation

      MODELE= 0
      MODELA= 1
      MODELC= 0
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

      MODEL_LOSS=0
      MODEL_SYNCH=0
      MODEL_NBI=0
      MODEL_WAVE=0 ! 0=no wave calc., 1=wave calc.

!-----COMPUTATION PARAMETERS------------------------------------------
!     DELT  : time step size (s)
!     RIMPL : implicit computation parameter
!     imtx  : type of matrix solver 
!                    0: petsc ksp (GMRES ILU(0) no initial guess)
!                    1: petsc ksp (GMRES ILU(0)    initial guess) *default
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

!     ***********************
!          PARAMETER INPUT
!     ***********************

      SUBROUTINE fp_parm(mode,kin,ierr)

!     mode=0 : standard namelinst input
!     mode=1 : namelist file input
!     mode=2 : namelist line input

!     ierr=0 : normal end
!     ierr=1 : namelist standard input error
!     ierr=2 : namelist file does not exist
!     ierr=3 : namelist file open error
!     ierr=4 : namelist file read error
!     ierr=5 : namelist file abormal end of file
!     ierr=6 : namelist line input error
!     ierr=7 : unknown MODE
!     ierr=10X : input parameter out of range

      IMPLICIT NONE
      INTEGER,INTENT(IN):: mode
      CHARACTER(LEN=*),INTENT(IN)::  kin
      INTEGER,INTENT(OUT):: ierr

    1 CALL task_parm(mode,'FP',kin,fp_nlin,fp_plst,ierr)
      IF(ierr.NE.0) RETURN

      CALL fp_check(ierr)
      IF(mode.EQ.0.AND.ierr.NE.0) GO TO 1
      IF(ierr.NE.0) ierr=ierr+100

      RETURN
      END SUBROUTINE fp_parm

!     ****** INPUT NAMELIST ******

      SUBROUTINE fp_nlin(nid,ist,ierr)

      use fpcomm_parm
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: nid
      INTEGER,INTENT(OUT) :: ist,ierr

      NAMELIST /FP/ &
           NSMAX,MODELG,MODELN,MODELQ,IDEBUG,MODEFR,MODEFW, &
           RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           RHOMIN,QMIN,RHOEDG,RHOITB,RHOGMN,RHOGMX, &
           PA,PZ,PZ0,PN,PNS,PTPR,PTPP,PTS,PU,PUS, &
           PNITB,PTITB,PUITB,PZCL, &
           KNAMEQ,KNAMWR,KNAMFP,KNAMWM,KNAMPF, &
           KNAMFO,KNAMTR,KNAMEQ2,KID_NS,ID_NS, &
           NSAMAX,NSBMAX,NS_NSA,NS_NSB, &
           LMAXNWR,NCMIN,NCMAX,NBEAMMAX,NSSPB,NSSPF, &
           NPMAX,NTHMAX,NRMAX,NAVMAX,NP2MAX, &
           NTMAX,NTSTEP_COEF,NTSTEP_COLL, &
           NTG1STEP,NTG1MIN,NTG1MAX, &
           NTG2STEP,NTG2MIN,NTG2MAX, &
           MODELE,MODELA,MODELC,MODELR,MODELS,MODELW, &
           MODELD,MODELD_RDEP,MODELD_PDEP,MODELD_EDGE, &
           MODELD_PINCH,MODELD_BOUNDARY, &
           MODEL_LOSS,MODEL_SYNCH,MODEL_NBI,MODEL_WAVE, &
           IMTX,MODEL_KSP,MODEL_PC,LMAXFP,LMAXE, &
           NGLINE,NGRAPH,LLMAX,LLMAX_NF,IDBGFP, &
           MODEL_DISRUPT,MODEL_Connor_fp,MODEL_BS,MODEL_jfp, &
           MODEL_LNL,MODEL_RE_pmax,MODELD_n_RE,MODEL_IMPURITY, &
           MODEL_SINK,N_IMPU, &
           N_partition_r,N_partition_s,N_partition_p, &
           PMAX,PMAX_BB, &
           R1,DELR1,RMIN,RMAX,E0,ZEFF, &
           PWAVE,RFDW,DELNPR,EPSNWR,REWY,DREWY,FACTWM, &
           DEC,PEC1,PEC2,PEC3,PEC4,RFEC,DELYEC, &
           DLH,PLH1,PLH2,RLH,DFW,PFW1,PFW2,RFW, &
           CEWR,CEWTH,CEWPH,RKWR,RKWTH,RKWPH, &
           SPBTOT,SPBR0,SPBRW,SPBENG,SPBANG,SPBPANG, &
           SPFTOT,SPFR0,SPFRW,SPFENG, &
           DRR0,DRRS,FACTOR_CDBM,DRR_EDGE,RHO_EDGE, &
           FACTOR_DRR_EDGE,FACTOR_PINCH,deltaB_B,TLOSS, &
           DELT,RIMPL,EPSFP,EPSM,EPSE,EPSDE,H0DE, &
           PGMAX,RGMAX,RGMIN, &
           T0_quench,tau_quench,tau_mgi, &
           time_quench_start,RJPROF1,RJPROF2, &
           v_RE,target_zeff,SPITOT

      READ(nid,FP,IOSTAT=ist,ERR=9800,END=9900)

      ierr=0
      RETURN

 9800 ierr=8
      RETURN
 9900 ierr=9
      RETURN
      END SUBROUTINE fp_nlin

!     ***** INPUT PARAMETER LIST *****

      SUBROUTINE fp_plst

      WRITE(6,*) '&FP : NSMAX,MODELG,MODELN,MODELQ,IDEBUG,MODEFR,MODEFW,'
      WRITE(6,*) '      RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'
      WRITE(6,*) '      PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'
      WRITE(6,*) '      RHOMIN,QMIN,RHOEDG,RHOITB,RHOGMN,RHOGMX,'
      WRITE(6,*) '      PA,PZ,PZ0,PN,PNS,PTPR,PTPP,PTS,PU,PUS,'
      WRITE(6,*) '      PNITB,PTITB,PUITB,PZCL,'
      WRITE(6,*) '      KNAMEQ,KNAMWR,KNAMFP,KNAMWM,KNAMPF,'
      WRITE(6,*) '      KNAMFO,KNAMTR,KNAMEQ2,KID_NS,ID_NS,'
      WRITE(6,*) '      NSAMAX,NSBMAX,NS_NSA,NS_NSB,'
      WRITE(6,*) '      LMAXNWR,NCMIN,NCMAX,NBEAMMAX,NSSPB,NSSPF,'
      WRITE(6,*) '      NPMAX,NTHMAX,NRMAX,NAVMAX,NP2MAX,'
      WRITE(6,*) '      NTMAX,NTSTEP_COEF,NTSTEP_COLL,'
      WRITE(6,*) '      NTG1STEP,NTG1MIN,NTG1MAX,'
      WRITE(6,*) '      NTG2STEP,NTG2MIN,NTG2MAX,'
      WRITE(6,*) '      MODELE,MODELA,MODELC,MODELR,MODELS,MODELW,'
      WRITE(6,*) '      MODELD,MODELD_RDEP,MODELD_PDEP,MODELD_EDGE,'
      WRITE(6,*) '      MODELD_BOUNDARY,MODELD_PINCH,'
      WRITE(6,*) '      MODEL_LOSS,MODEL_SYNCH,MODEL_NBI,MODEL_WAVE,'
      WRITE(6,*) '      IMTX,MODEL_KSP,MODEL_PC,LMAXFP,LMAXE,'
      WRITE(6,*) '      NGLINE,NGRAPH,LLMAX,LLMAX_NF,IDBGFP,'
      WRITE(6,*) '      MODEL_DISRUPT,MODEL_Connor_fp,MODEL_BS,MODEL_jfp,'
      WRITE(6,*) '      MODEL_LNL,MODEL_RE_pmax,MODELD_n_RE,MODEL_IMPURITY,'
      WRITE(6,*) '      MODEL_SINK,N_IMPU,'
      WRITE(6,*) '      N_partition_r,N_partition_s,N_partition_p,'
      WRITE(6,*) '      PMAX,PMAX_BB,'
      WRITE(6,*) '      R1,DELR1,RMIN,RMAX,E0,ZEFF,'
      WRITE(6,*) '      PWAVE,RFDW,DELNPR,EPSNWR,REWY,DREWY,FACTWM,'
      WRITE(6,*) '      DEC,PEC1,PEC2,PEC3,PEC4,RFEC,DELYEC,'
      WRITE(6,*) '      DLH,PLH1,PLH2,RLH,DFW,PFW1,PFW2,RFW,'
      WRITE(6,*) '      CEWR,CEWTH,CEWPH,RKWR,RKWTH,RKWPH,'
      WRITE(6,*) '      SPBTOT,SPBR0,SPBRW,SPBENG,SPBANG,SPBPANG,'
      WRITE(6,*) '      SPFTOT,SPFR0,SPFRW,SPFENG,'
      WRITE(6,*) '      DRR0,DRRS,FACTOR_CDBM,DRR_EDGE,RHO_EDGE,'
      WRITE(6,*) '      FACTOR_DRR_EDGE,FACTOR_PINCH,deltaB_B,TLOSS,'
      WRITE(6,*) '      DELT,RIMPL,EPSFP,EPSM,EPSE,EPSDE,H0DE,'
      WRITE(6,*) '      PGMAX,RGMAX,RGMIN,'
      WRITE(6,*) '      T0_quench,tau_quench,tau_mgi,'
      WRITE(6,*) '      time_quench_start,RJPROF1,RJPROF2,'
      WRITE(6,*) '      v_RE,target_zeff,SPITOT,'

      RETURN
    END SUBROUTINE fp_plst

!     ***** CHECK INPUT PARAMETERS *****

      SUBROUTINE fp_check(ierr)

      USE fpcomm
      IMPLICIT NONE
      integer,intent(out):: ierr

      ierr=0

      RETURN
      END SUBROUTINE fp_check

!     ***** BROADCAST INPUT PARAMETERS *****

      SUBROUTINE fp_broadcast

      USE fpcomm_parm
      USE libmpi
      USE libmtx
      IMPLICIT NONE
      INTEGER,DIMENSION(99):: idata
      real(8),DIMENSION(99):: rdata
      complex(8),DIMENSION(3):: cdata
      INTEGER:: NS

!----- PL input parameters -----     

      idata( 1)=NSMAX
      idata( 2)=MODELG
      idata( 3)=MODELN
      idata( 4)=MODELQ
      idata( 5)=IDEBUG
      idata( 6)=MODEFR
      idata( 7)=MODEFW

      CALL mtx_broadcast_integer(idata,7)
      NSMAX =idata( 1)
      MODELG=idata( 2)
      MODELN=idata( 3)
      MODELQ=idata( 4)
      IDEBUG=idata( 5)
      MODEFR=idata( 6)
      MODEFW=idata( 7)

      rdata( 1)=RR
      rdata( 2)=RA
      rdata( 3)=RB
      rdata( 4)=RKAP
      rdata( 5)=RDLT
      rdata( 6)=BB
      rdata( 7)=Q0
      rdata( 8)=QA
      rdata( 9)=RIP
      rdata(10)=PROFJ
      rdata(11)=PROFN1
      rdata(12)=PROFN2
      rdata(13)=PROFT1
      rdata(14)=PROFT2
      rdata(15)=PROFU1
      rdata(16)=PROFU2
      rdata(17)=RHOMIN
      rdata(18)=QMIN
      rdata(19)=RHOEDG
      rdata(20)=RHOITB
      rdata(21)=RHOGMN
      rdata(22)=RHOGMX

      CALL mtx_broadcast_real8(rdata,22)
      RR    =rdata( 1)
      RA    =rdata( 2)
      RB    =rdata( 3)
      RKAP  =rdata( 4)
      RDLT  =rdata( 5)
      BB    =rdata( 6)
      Q0    =rdata( 7)
      QA    =rdata( 8)
      RIP   =rdata( 9)
      PROFJ =rdata(10)
      PROFN1=rdata(11)
      PROFN2=rdata(12)
      PROFT1=rdata(13)
      PROFT2=rdata(14)
      PROFU1=rdata(15)
      PROFU2=rdata(16)
      RHOMIN=rdata(17)
      QMIN  =rdata(18)
      RHOEDG=rdata(19)
      RHOITB=rdata(20)
      RHOGMN=rdata(21)
      RHOGMX=rdata(22)

      CALL mtx_broadcast_real8(PA,NSMAX)
      CALL mtx_broadcast_real8(PZ,NSMAX)
      CALL mtx_broadcast_real8(PZ0,NSMAX)
      CALL mtx_broadcast_integer(ID_NS,NSMAX)
      CALL mtx_broadcast_real8(PN,NSMAX)
      CALL mtx_broadcast_real8(PNS,NSMAX)
      CALL mtx_broadcast_real8(PTPR,NSMAX)
      CALL mtx_broadcast_real8(PTPP,NSMAX)
      CALL mtx_broadcast_real8(PTS,NSMAX)
      CALL mtx_broadcast_real8(PU,NSMAX)
      CALL mtx_broadcast_real8(PUS,NSMAX)
      CALL mtx_broadcast_real8(PNITB,NSMAX)
      CALL mtx_broadcast_real8(PTITB,NSMAX)
      CALL mtx_broadcast_real8(PUITB,NSMAX)
      CALL mtx_broadcast_real8(PZCL,NSMAX)

      CALL mtx_broadcast_character(KNAMEQ,80)
      CALL mtx_broadcast_character(KNAMWR,80)
      CALL mtx_broadcast_character(KNAMFP,80)
      CALL mtx_broadcast_character(KNAMWM,80)
      CALL mtx_broadcast_character(KNAMPF,80)
      CALL mtx_broadcast_character(KNAMFO,80)
      CALL mtx_broadcast_character(KNAMTR,80)
      CALL mtx_broadcast_character(KNAMEQ2,80)
      DO NS=1,NSMAX
         CALL mtx_broadcast_character(KID_NS(NS),2)
      END DO

!----- FP input parameters -----

      idata( 1)=NSAMAX
      idata( 2)=NSBMAX
      idata( 3)=LMAXNWR
      idata( 4)=NBEAMMAX
      idata( 5)=NSSPF
      idata( 6)=NPMAX
      idata( 7)=NTHMAX
      idata( 8)=NRMAX
      idata( 9)=NAVMAX
      idata(10)=NP2MAX
      idata(11)=NTMAX
      idata(12)=NTSTEP_COEF
      idata(13)=NTSTEP_COLL
      idata(14)=NTG1STEP
      idata(15)=NTG1MIN
      idata(16)=NTG1MAX
      idata(17)=NTG2STEP
      idata(18)=NTG2MIN
      idata(19)=NTG2MAX
      idata(20)=MODELE
      idata(21)=MODELA
      idata(22)=MODELC
      idata(23)=MODELR
      idata(24)=MODELS
      idata(25)=MODELD
      idata(26)=MODELD_RDEP
      idata(27)=MODELD_PDEP
      idata(28)=MODELD_EDGE
      idata(29)=MODELD_PINCH
      idata(30)=MODELD_BOUNDARY
      idata(31)=MODEL_LOSS
      idata(32)=MODEL_SYNCH
      idata(33)=MODEL_NBI
      idata(34)=MODEL_WAVE
      idata(35)=IMTX
      idata(36)=MODEL_KSP
      idata(37)=MODEL_PC
      idata(38)=LMAXFP
      idata(39)=LMAXE
      idata(40)=NGLINE
      idata(41)=NGRAPH
      idata(42)=LLMAX
      idata(43)=LLMAX_NF
      idata(44)=IDBGFP
      idata(45)=MODEL_DISRUPT
      idata(46)=MODEL_Connor_FP
      idata(47)=MODEL_BS
      idata(48)=MODEL_jfp
      idata(49)=MODEL_LNL
      idata(50)=MODEL_RE_pmax
      idata(51)=MODELD_n_RE
      idata(52)=MODEL_IMPURITY
      idata(53)=MODEL_SINK
      idata(54)=n_impu
      idata(55)=N_partition_r
      idata(56)=N_partition_s
      idata(57)=N_partition_p

      CALL mtx_broadcast_integer(idata,57)
      NSAMAX         =idata( 1)
      NSBMAX         =idata( 2)
      LMAXNWR        =idata( 3)
      NBEAMMAX       =idata( 4)
      NSSPF          =idata( 5)
      NPMAX          =idata( 6)
      NTHMAX         =idata( 7)
      NRMAX          =idata( 8)
      NAVMAX         =idata( 9)
      NP2MAX         =idata(10)
      NTMAX          =idata(11)
      NTSTEP_COEF    =idata(12)
      NTSTEP_COLL    =idata(13)
      NTG1STEP       =idata(14)
      NTG1MIN        =idata(15)
      NTG1MAX        =idata(16)
      NTG2STEP       =idata(17)
      NTG2MIN        =idata(18)
      NTG2MAX        =idata(19)
      MODELE         =idata(20)
      MODELA         =idata(21)
      MODELC         =idata(22)
      MODELR         =idata(23)
      MODELS         =idata(24)
      MODELD         =idata(25)
      MODELD_RDEP    =idata(26)
      MODELD_PDEP    =idata(27)
      MODELD_EDGE    =idata(28)
      MODELD_PINCH   =idata(29)
      MODELD_BOUNDARY=idata(30)
      MODEL_LOSS     =idata(31)
      MODEL_SYNCH    =idata(32)
      MODEL_NBI      =idata(33)
      MODEL_WAVE     =idata(34)
      IMTX           =idata(35)
      MODEL_KSP      =idata(36)
      MODEL_PC       =idata(37)
      LMAXFP         =idata(38)
      LMAXE          =idata(39)
      NGLINE         =idata(40)
      NGRAPH         =idata(41)
      LLMAX          =idata(42)
      LLMAX_NF       =idata(43)
      IDBGFP         =idata(44)
      MODEL_DISRUPT  =idata(45)
      MODEL_Connor_FP=idata(46)
      MODEL_BS       =idata(47)
      MODEL_jfp      =idata(48)
      MODEL_LNL      =idata(49)
      MODEL_RE_pmax  =idata(50)
      MODELD_n_RE    =idata(51)
      MODEL_IMPURITY =idata(52)
      MODEL_SINK     =idata(53)
      n_impu         =idata(54)
      N_partition_r  =idata(55)
      N_partition_s  =idata(56)
      N_partition_p  =idata(57)

      CALL mtx_broadcast_integer(NS_NSA,NSAMAX)
      CALL mtx_broadcast_integer(NS_NSB,NSBMAX)
      CALL mtx_broadcast_integer(NCMIN,NSAMAX)
      CALL mtx_broadcast_integer(NCMAX,NSAMAX)
      CALL mtx_broadcast_integer(NSSPB,NBEAMMAX)
      CALL mtx_broadcast_integer(MODELW,NSAMAX)

      rdata( 1)=R1
      rdata( 2)=DELR1
      rdata( 3)=RMIN
      rdata( 4)=RMAX
      rdata( 5)=E0
      rdata( 6)=ZEFF
      rdata( 7)=PWAVE
      rdata( 8)=RFDW
      rdata( 9)=DELNPR
      rdata(10)=EPSNWR
      rdata(11)=REWY
      rdata(12)=DREWY
      rdata(13)=FACTWM
      rdata(14)=DEC
      rdata(15)=PEC1
      rdata(16)=PEC2
      rdata(17)=PEC3
      rdata(18)=PEC4
      rdata(19)=RFEC
      rdata(20)=DELYEC
      rdata(21)=DLH
      rdata(22)=PLH1
      rdata(23)=PLH2
      rdata(24)=RLH
      rdata(25)=DFW
      rdata(26)=PFW1
      rdata(27)=PFW2
      rdata(28)=RFW
      rdata(29)=RKWR
      rdata(30)=RKWTH
      rdata(31)=RKWPH
      rdata(32)=SPFTOT
      rdata(33)=SPFR0
      rdata(34)=SPFRW
      rdata(35)=SPFENG
      rdata(36)=DRR0
      rdata(37)=DRRS
      rdata(38)=FACTOR_CDBM
      rdata(39)=DRR_EDGE
      rdata(40)=RHO_EDGE
      rdata(41)=FACTOR_DRR_EDGE
      rdata(42)=FACTOR_PINCH
      rdata(43)=DELTAB_B
      rdata(44)=DELT
      rdata(45)=RIMPL
      rdata(46)=EPSFP
      rdata(47)=EPSM
      rdata(48)=EPSE
      rdata(49)=EPSDE
      rdata(50)=H0DE
      rdata(51)=PGMAX
      rdata(52)=RGMAX
      rdata(53)=RGMIN
      rdata(54)=T0_quench
      rdata(55)=tau_quench
      rdata(56)=tau_mgi
      rdata(57)=time_quench_start
      rdata(58)=RJPROF1
      rdata(59)=RJPROF2
      rdata(60)=v_RE
      rdata(61)=target_zeff
      rdata(62)=SPITOT

      CALL mtx_broadcast_real8(rdata,62)
      R1               =rdata( 1)
      DELR1            =rdata( 2)
      RMIN             =rdata( 3)
      RMAX             =rdata( 4)
      E0               =rdata( 5)
      ZEFF             =rdata( 6)
      PWAVE            =rdata( 7)
      RFDW             =rdata( 8)
      DELNPR           =rdata( 9)
      EPSNWR           =rdata(10)
      REWY             =rdata(11)
      DREWY            =rdata(12)
      FACTWM           =rdata(13)
      DEC              =rdata(14)
      PEC1             =rdata(15)
      PEC2             =rdata(16)
      PEC3             =rdata(17)
      PEC4             =rdata(18)
      RFEC             =rdata(19)
      DELYEC           =rdata(20)
      DLH              =rdata(21)
      PLH1             =rdata(22)
      PLH2             =rdata(23)
      RLH              =rdata(24)
      DFW              =rdata(25)
      PFW1             =rdata(26)
      PFW2             =rdata(27)
      RFW              =rdata(28)
      RKWR             =rdata(29)
      RKWTH            =rdata(30)
      RKWPH            =rdata(31)
      SPFTOT           =rdata(32)
      SPFR0            =rdata(33)
      SPFRW            =rdata(34)
      SPFENG           =rdata(35)
      DRR0             =rdata(36)
      DRRS             =rdata(37)
      FACTOR_CDBM      =rdata(38)
      DRR_EDGE         =rdata(39)
      RHO_EDGE         =rdata(40)
      FACTOR_DRR_EDGE  =rdata(41)
      FACTOR_PINCH     =rdata(42)
      DELTAB_B         =rdata(43)
      DELT             =rdata(44)
      RIMPL            =rdata(45)
      EPSFP            =rdata(46)
      EPSM             =rdata(47)
      EPSE             =rdata(48)
      EPSDE            =rdata(49)
      H0DE             =rdata(50)
      PGMAX            =rdata(51)
      RGMAX            =rdata(52)
      RGMIN            =rdata(53)
      T0_quench        =rdata(54)
      tau_quench       =rdata(55)
      tau_mgi          =rdata(56)
      time_quench_start=rdata(57)
      RJPROF1          =rdata(58)
      RJPROF2          =rdata(59)
      v_RE             =rdata(60)
      target_zeff      =rdata(61)
      SPITOT           =rdata(62)

      CALL mtx_broadcast_real8(pmax,NSAMAX)
      CALL mtx_broadcast_real8(pmax_bb,NSAMAX)
      CALL mtx_broadcast_real8(TLOSS,NSMAX)

      CALL mtx_broadcast_real8(SPBTOT,NBEAMMAX)
      CALL mtx_broadcast_real8(SPBR0 ,NBEAMMAX)
      CALL mtx_broadcast_real8(SPBRW ,NBEAMMAX)
      CALL mtx_broadcast_real8(SPBENG,NBEAMMAX)
      CALL mtx_broadcast_real8(SPBANG,NBEAMMAX)
      CALL mtx_broadcast_real8(SPBPANG,NBEAMMAX)

      cdata (1)=CEWR
      cdata (2)=CEWTH
      cdata (3)=CEWPH
      CALL mtx_broadcast_complex8(cdata,3)
      CEWR =cdata(1)
      CEWTH=cdata(2)
      CEWPH=cdata(3)

      RETURN
      END SUBROUTINE fp_broadcast

!     ***********************
!          PARAMETER VIEW
!     ***********************

      SUBROUTINE fp_view

      use fpcomm_parm
      IMPLICIT NONE
      integer:: nsa,nsb,ns,NBEAM

      WRITE(6,600) 'E0      ',E0      ,'ZEFF    ',ZEFF
      IF(NRMAX.EQ.1) THEN
         WRITE(6,600) &
              'R1      ',R1      ,'DELR1   ',DELR1   ,'DELT    ',DELT 
      ELSE
         WRITE(6,600) &
              'RMIN    ',RMIN    ,'RMAX    ',RMAX    ,'DELT    ',DELT 
      ENDIF

      WRITE(*,*) "----- ANALYZED SPECIES ----- "
      WRITE(*,*) "-----    TEST SPECIES    ----- "
      DO nsa=1,nsamax
         WRITE(6,603) 'ns_nsa  ',ns_nsa(nsa)
      END DO
      WRITE(*,*) "-----    BACKGROUND SPECIES    ----- "
      DO nsb=1,nsbmax
         WRITE(6,603) 'ns_nsb  ',ns_nsb(nsb)
      END DO
      WRITE(*,*) "----- Maximum normalized momentum and loss time ----- "
      DO nsb=1,nsbmax
         WRITE(6,600) 'pmax    ',pmax(nsb), 'pmax_bb ',pmax_bb(nsb), &
                      'tloss   ',tloss(nsb)
      END DO

      WRITE(*,*) "----- PARAMETERS OF HEATINGS -----"
      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         WRITE(6,'(A,I2,A,I2)') 'NSA = ',NSA,'  NS  = ',NS
         
         IF(MODELW(NS).EQ.0) THEN
            
            WRITE(6,600) 'DEC     ',DEC     ,'RFEC    ',RFEC    , &
                         'DELYEC  ',DELYEC
            WRITE(6,600) 'PEC1    ',PEC1    ,'PEC2    ',PEC2
            WRITE(6,600) 'PEC3    ',PEC3    ,'PEC4    ',PEC4
            WRITE(6,600) 'DLH     ',DLH     ,'RLH     ',RLH
            WRITE(6,600) 'PLH1    ',PLH1    ,'PLH2    ',PLH2
            WRITE(6,600) 'DFW     ',DFW     ,'RFW     ',RFW
            WRITE(6,600) 'PFW1    ',PFW1    ,'PFW2    ',PFW2
            
         ELSEIF(MODELW(NS).EQ.1) THEN
            WRITE(6,600) 'RFDW    ',RFDW    ,'DELNPR  ',DELNPR
            WRITE(6,600) 'PWAVE   ',PWAVE   ,'DELYEC  ',DELYEC
            WRITE(6,602) 'EPSNWR  ',EPSNWR  ,'LMAXNWR ',LMAXNWR
            
         ELSEIF(MODELW(NS).EQ.2) THEN
            WRITE(6,600) 'RFDW    ',RFDW    ,'DELNPR  ',DELNPR
            WRITE(6,600) 'PWAVE   ',PWAVE   ,'DELYEC  ',DELYEC
            WRITE(6,602) 'EPSNWR  ',EPSNWR  ,'LMAXNWR ',LMAXNWR
            
         ELSEIF(MODELW(NS).EQ.3) THEN
            WRITE(6,600) 'RFDW    ',RFDW    ,'DELNPR  ',DELNPR
            WRITE(6,600) 'CEWR/R  ',DBLE(CEWR) ,'CEWR/I  ',DIMAG(CEWR)
            WRITE(6,600) 'CEWTH/R ',DBLE(CEWTH),'CEWTH/I ',DIMAG(CEWTH)
            WRITE(6,600) 'CEWPH/R ',DBLE(CEWPH),'CEWPH/I ',DIMAG(CEWPH)
            WRITE(6,600) 'REWY    ',REWY    ,'DREWY   ',DREWY
            WRITE(6,600) 'RKWR    ',RKWR    ,'RKWTH   ',RKWTH   , &
                         'RKWPH   ',RKWPH
            
         ELSEIF(MODELW(NS).EQ.4) THEN
            WRITE(6,600) 'RFDW    ',RFDW    ,'DELNPR  ',DELNPR
            WRITE(6,600) 'FACTWM  ',FACTWM
         ENDIF
         
         IF(TLOSS(NS).NE.0.D0) THEN
            WRITE(6,600) 'TLOSS   ',TLOSS(NS)
         ENDIF

         IF(MODELW(NS).EQ.0)THEN
            WRITE(6,*) 'GIVEN WAVE DIFFUSION COEFFICIENTS'
         ELSE IF(MODELW(NS).EQ.1)THEN
            WRITE(6,*) 'RAY TRACING WAVE DATA'
         ELSE IF(MODELW(NS).EQ.2)THEN
            WRITE(6,*) 'BEAM TRACING WAVE DATA'
         ELSE IF(MODELW(NS).EQ.3)THEN
            WRITE(6,*) 'GIVEN WAVE AMPLITUDE'
         ELSE IF(MODELW(NS).EQ.4)THEN
            WRITE(6,*) 'FULL WAVE DATA'
         ELSE
            WRITE(6,*) 'XX UNKNOWN MODELW: MODELW =',MODELW(NSA)
         END IF
      
      END DO

      DO NBEAM=1,NBEAMMAX
         IF(SPBTOT(NBEAM).NE.0.D0) THEN
            WRITE(6,*) "NBI NUMBER NBEAM & NSSPB=",NBEAM,NSSPB(NBEAM)
            WRITE(6,600) 'SPBTOT   ',SPBTOT(NBEAM), &
                         'SPBENG   ',SPBENG(NBEAM), &
                         'SPBANG   ',SPBANG(NBEAM)
            WRITE(6,600) 'SPBPANG   ',SPBPANG(NBEAM), &
                         'SPBR0    ',SPBR0(NBEAM), &
                         'SPBRW    ',SPBRW(NBEAM)
         ENDIF
      END DO

      IF(MODELS.eq.0)THEN
         WRITE(6,*) 'No fusion reaction'
      ELSEIF(MODELS.eq.1)THEN
         IF(SPFTOT.NE.0.D0) THEN
            WRITE(6,*) 'Given profile fusion reaction'
            WRITE(6,*) "NSSPF=",NSSPF
            WRITE(6,600) 'SPFTOT   ',SPFTOT, &
                         'SPFENG   ',SPFENG
            WRITE(6,600) 'SPFR0    ',SPFR0, &
                         'SPFRW    ',SPFRW
         ENDIF
      ELSE
         WRITE(6,*) 'Self-consistent fusion reaction'
      END IF

      WRITE(6,604)'DRR0            ',DRR0            , &
                  'DRRS            ',DRRS
      WRITE(6,604)'FACTOR_CDBM     ',FACTOR_CDBM     , &
                  'DELTAB_B        ',DELTAB_B
      WRITE(6,604)'DRR_EDGE        ',DRR_EDGE        , &
                  'RHO_EDGE        ',RHO_EDGE
      WRITE(6,604)'FACTOR_DRR_EDGE ',FACTOR_DRR_EDGE , &
                  'FACTOR_PINCH    ',FACTOR_PINCH


      WRITE(6,*) "----- OTHER PARAMETERS -----"
      WRITE(6,600) 'RIMPL   ',RIMPL   ,'EPSM    ',EPSM    ,'EPSFP   ',EPSFP
      WRITE(6,600) 'EPSE    ',EPSE    ,'EPSDE   ',EPSDE   ,'H0DE    ',H0DE    
      WRITE(6,600) 'PGMAX   ',PGMAX   ,'RGMAX   ',RGMAX   ,'RGMIN   ',RGMIN
      WRITE(6,603) 'LMAXE   ',LMAXE   ,'LLMAX   ',LLMAX   ,'NGLINE  ',NGLINE
      WRITE(6,603) 'IDBGFP  ',IDBGFP  ,'NGRAPH  ',NGRAPH  ,'LLMAX_NF',LLMAX_NF

      WRITE(6,603) 'NPMAX   ',NPMAX   ,'NTHMAX  ',NTHMAX  ,'NRMAX   ',NRMAX
      WRITE(6,603) 'NAVMAX  ',NAVMAX  ,'NP2MAX  ',NP2MAX  ,'NTMAX   ',NTMAX   
      WRITE(6,603) 'NTG1STEP',NTG1STEP,'NTG1MIN ',NTG1MIN ,'NTG1MAX ',NTG1MAX
      WRITE(6,603) 'NTG2STEP',NTG2STEP,'NTG2MIN ',NTG2MIN ,'NTG2MAX ',NTG2MAX
      WRITE(6,606) 'NTSTEP_COEF     ',NTSTEP_COEF, &
                   'NTSTEP_COLL     ',NTSTEP_COLL
      WRITE(6,603) 'MODELE  ',MODELE  ,'MODELA  ',MODELA  ,'MODELC  ',MODELC
      WRITE(6,603) 'MODELR  ',MODELR  ,'MODELS  ',MODELS  ,'MODELD  ',MODELD
      WRITE(6,606) 'MODELD_RDEP     ',MODELD_RDEP    , &
                   'MODELD_PDEP     ',MODELD_PDEP
      WRITE(6,606) 'MODELD_EDGE     ',MODELD_EDGE    , &
                   'MODELD_PINCH    ',MODELD_PINCH
      WRITE(6,606) 'MODELD_BOUNDARY ',MODELD_BOUNDARY
      WRITE(6,606) 'MODEL_LOSS      ',MODEL_LOSS     , &
                   'MODEL_SYNCH     ',MODEL_SYNCH
      WRITE(6,606) 'MODEL_NBI       ',MODEL_NBI      , &
                   'MODEL_WAVE      ',MODEL_WAVE
      WRITE(6,603) 'IMTX    ',IMTX    , &
                   'LMAXFP  ',LMAXFP  , &
                   'LMAXE   ',LMAXE
      WRITE(6,606) 'MODEL_KSP       ',MODEL_KSP      , &
                   'MODEL_PC        ',MODEL_PC
      WRITE(6,603) 'NGLINE  ',NGLINE  , &
                   'NGRAPH  ',NGRAPH
      WRITE(6,603) 'LLMAX   ',LLMAX   , &
                   'LLMAX_NF',LLMAX_NF, &
                   'IDBGFP  ',IDBGFP
      WRITE(6,606) 'MODEL_DISRUPT   ',MODEL_DISRUPT   , &
                   'MODEL_Connor_fp ',MODEL_Connor_fp
      WRITE(6,606) 'MODEL_BS        ',MODEL_BS        , &
                   'MODEL_jfp       ',MODEL_jfp
      WRITE(6,606) 'MODEL_LNL       ',MODEL_LNL       , &
                   'MODEL_RE_pmax   ',MODEL_RE_pmax
      WRITE(6,606) 'MODELD_n_RE     ',MODELD_n_RE     , &
                   'MODEL_IMPURITY  ',MODEL_IMPURITY
      WRITE(6,606) 'MODEL_LNL       ',MODEL_LNL       , &
                   'MODEL_RE_pmax   ',MODEL_RE_pmax
      WRITE(6,606) 'MODEL_SINK      ',MODEL_SINK      , &
                   'N_IMPU          ',N_IMPU
      WRITE(6,604) 'T0_quench       ',T0_quench       , &
                   'tau_quench      ',tau_quench
      WRITE(6,604) 'tau_mgi         ',tau_mgi         , &
                   'time_quench_star',time_quench_start
      WRITE(6,600) 'RJPROF1 ',RJPROF1 , &
                   'RJPROF2 ',RJPROF2 , &
                   'v_RE    ',v_RE
      WRITE(6,604) 'target_zeff     ',target_zeff     , &
                   'SPITOT          ',SPITOT

      WRITE(6,*) "-------- PLASMA MODELS --------"

      IF(MODELE.EQ.0)THEN
         WRITE(6,*) 'FIXED ELECTRIC FIELD'
      ELSE IF(MODELE.EQ.1)THEN
         WRITE(6,*) 'CONSISTENT ELECTRIC FIELD'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELE: MODELE =',MODELE
      ENDIF

      IF(MODELR.EQ.0)THEN
         WRITE(6,*) 'NONRELATIVISTIC'
      ELSE IF(MODELR.EQ.1)THEN
         WRITE(6,*) 'RELATIVISTIC'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELR: MODELR =',MODELR
      ENDIF

      IF(MODELD.EQ.0)THEN
         WRITE(6,*) 'WITHOUT RADIAL TRANPORT'
      ELSE IF(MODELD.EQ.1)THEN
         WRITE(6,*) 'WITH RADIAL TRANSPORT (const. for r,p,th, without pinch)'
      ELSE IF(MODELD.EQ.2)THEN
         WRITE(6,*) 'WITH RADIAL TRANSPORT (const. for r,p,th, with pinch)'
      ELSE IF(MODELD.EQ.3)THEN
         WRITE(6,*) 'WITH RADIAL TRANSPORT (p dependence without pinch)'
      ELSE IF(MODELD.EQ.4)THEN
         WRITE(6,*) 'WITH RADIAL TRANSPORT (p dependence with pinch)'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELD: MODELD =',MODELD
      ENDIF

      IF(MODELA.EQ.0)THEN
         WRITE(6,*) 'NOT BOUNCE AVERAGED'
      ELSE IF(MODELA.EQ.1)THEN
         WRITE(6,*) 'BOUNCE AVERAGED'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELA: MODELA =',MODELA
      ENDIF

      IF(MODELC.EQ.0)THEN
         WRITE(6,*) 'LINEAR COLLISION OPERATOR & CONST. T'
      ELSE IF(MODELC.eq.1)THEN
         WRITE(6,*) 'LINEAR COLLISION OPERATOR & VARIABLE. T'
      ELSE IF(MODELC.EQ.2)THEN
         WRITE(6,*) &
              'NONLINEAR COLLISION OPERATOR FOR LIKE PARTILCES & CONST. T'
      ELSE IF(MODELC.eq.3)THEN
         WRITE(6,*) &
              'NONLINEAR COLLISION OPERATOR FOR LIKE PARTILCES & VARIABLE T'
      ELSE IF(MODELC.eq.4)THEN
         WRITE(6,*) 'NONLINEAR COLLISION OPERATOR'
      ELSE IF(MODELC.EQ.5)THEN
         WRITE(6,*) 'NONLINEAR COLLISION OPERATOR'
      ELSE IF(MODELC.EQ.-1)THEN
         WRITE(6,*) 'LINEAR COLLISION OPERATOR WITH ION SCATTERING'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELC: MODELC =',MODELC
      END IF

      IF(MODELG.EQ.2)THEN
         WRITE(6,*) 'GIVEN PLASMA GEOMETRY'
      ELSE IF(MODELG.EQ.3)THEN
         WRITE(6,*) 'MHD EQUILIBRIUM FROM TASK/EQ'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELG: MODELG =',MODELG
      END IF

      WRITE(6,*) "-------- MPI CONFIGURATION --------"
      WRITE(6,'(A,I4)') "THE NUMBER MPI PROCESS   =", nsize
      WRITE(6,'(A,I4)') "PARTITION NUMBER FOR NSA =", N_partition_s
      WRITE(6,'(A,I4)') "PARTITION NUMBER FOR NR  =", N_partition_r
      WRITE(6,'(A,I4)') "PARTITION NUMBER FOR NP  =", N_partition_p

      RETURN

  600 FORMAT(' ',A8,'=',1PE12.4:3X,A8,'=',1PE12.4:3X,A8,'=',1PE12.4)
  601 FORMAT(' ',A8,'=',1PE12.4:3X,A8,'=',1PE12.4:3X,A8,'=',I8)
  602 FORMAT(' ',A8,'=',1PE12.4:3X,A8,'=',I8,4X  :3X,A8,'=',I8)
  603 FORMAT(' ',A8,'=',I8,4X  :3X,A8,'=',I8,4X  :3X,A8,'=',I8)
  604 FORMAT(' ',A16,'=',1PE12.4:3X,A16,'=',1PE12.4)
  605 FORMAT(' ',A16,'=',1PE12.4:3X,A16,'=',I8)
  606 FORMAT(' ',A16,'=',I8,4X  :3X,A16,'=',I8)
    END SUBROUTINE fp_view

  END module fpinit
