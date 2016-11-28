!     $Id: fpinit.f90,v 1.18 2013/01/22 16:21:46 fukuyama Exp $

      module fpinit

      IMPLICIT NONE

      contains

!     ***************************
!         INITIAL PARAMETERS
!     ***************************

      SUBROUTINE fp_init

      use fpcomm
      integer:: ns,nsa,nsb,nbeam
!-----------------------------------------------------------------------
!     nsamax: number of test particle species
!     nsbmax: number of field particle species (0 for nsbmax=nsmax: default)
!     ns_nsa(nsa): mapping from NSA to NS in pl 
!     ns_nsb(nsb): mapping from NSB to NS in pl
!     pmax(nsb)  : maximum momentum (normailzed by central thermal momentum)
!     tloss(nsb) : loss time [s] (0.D0 for no loss)
!     zeff  : effective ion charge for simple collision term
!     imtx  : type of matrix solver 
!                    0: petsc ksp (GMRES ILU(0) no initial guess)
!                    1: petsc ksp (GMRES ILU(0)    initial guess) *default

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

      zeff  = 1.D0

      imtx=1

!-----------------------------------------------------------------------
!     DRR0  : radial diffusion coefficient at magnetic axis (m^2/s)
!     DRRS  : radial diffusion coefficient at plasma surface (m^2/s)
!     E0    : toroidal electric field (V/m)
!     R1    : radial position for NRMAX=1 (r/a)
!     DELR1 : radial spacing for NRMAX=1 (r/a)
!     RMIN  : minimum radius for NRMAX<>1 (r/a)
!     RMIN  : maximum radius for NRMAX<>1 (r/a)

      DRR0  = 0.D0
      DRRS  = 0.D0
      E0    = 0.D0
      R1    = 0.1D0*RR
      DELR1 = 0.1D0
      RMIN  = 0.05D0
      RMAX  = 0.3D0

!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
!     RFDW  : wave frequency [MHz]
!     DELNPR: width of toroidal mode number
!     CEWR  : radial component of wave electric field [V/m]
!     CEWTH : poloidal component of wave electric field [V/m]
!     CEWPH : toroidal component of wave electric field [V/m]
!     RKWR  : radial component of wave number [1/m]
!     RKWTH : poloidal component of wave number [1/m]
!     RKWPH : toroidal component of wave number [1/m]
!     REWY  : vertical position of ray [r/a]
!     DREWY : vertical half-width of ray [r/a]
!     FACTWM: Numerical factor for wave amplitude
!     NCMIN(NS): minimum order of cyclotron harmonics
!     NCMAX(NS): maximum order of cyclotron harmonics

      RFDW  = 170.D3
      DELNPR= 0.05D0
      CEWR  = (0.D0,0.D0)
      CEWTH = (0.D0,0.D0)
      CEWPH = (1.D3,0.D0)
      RKWR  = 0.D0
      RKWTH = 0.D0
      RKWPH = 1.D3
      REWY  = 0.D0
      DREWY = 0.1D0
      FACTWM= 1.D0
      DO NS=1,NSM
         NCMIN(NS) = -3
         NCMAX(NS) = 3
      ENDDO

!-----------------------------------------------------------------------
!     TLOSS(ns): loss time [s] (0.D0 for no loss)

!     NSSPB(nbeam) : NBI particle species
!     SPBTOT(nbeam): Particle source [1/m^3 s] ! [1/s]? 
!     SPBR0(nbeam) : Source radius [r/a]
!     SPBRW(nbeam) : Source width [r/a]
!     SPBENG(nbeam): Particle energy [eV]
!     SPBPANG(nbeam): Source poloidal angle [degree]
!     SPBANG(nbeam): Source pitch angle [degree] ! at chi=SPBPANG

!     NSSPF  : Fusion product particle species
!     SPFTOT : Particle source [1/m^3 s]
!     SPFR0  :  Source radius [r/a]
!     SPFRW  :  Source width [r/a]
!     SPFENG : Particle energy [eV]

      DO ns=1,nsm
         TLOSS(ns)=0.d0          ! default no loss
      ENDDO

      DO NBEAM=1,NBEAMM
         NSSPB(NBEAM)=2
         SPBTOT(NBEAM)=0.d0
         SPBR0(NBEAM)=0.d0
         SPBRW(NBEAM)=0.2d0
         SPBENG(NBEAM)=1.D6
         SPBANG(NBEAM)=20.D0
         SPBPANG(NBEAM)=0.D0
      ENDDO

      NSSPF=4
      SPFTOT=0.d0
      SPFR0=0.d0
      SPFRW=0.2d0
      SPFENG=3.5D6

!-----------------------------------------------------------------------
!     DELT  : time step size (s)
!     RIMPL : implicit computation parameter
!     EPSM  : convergence limit in matrix equation solver
!     EPSE  : convergence limit in electric field prediction
!     LMAXE : maximum loop count in electric field prediction
!     EPSDE : convergence limit in double-exponential integration
!     H0DE  : initial step size in double-exponential integration
!     PGMAX:  maximum p in graphics (if not 0)
!     RGMIN:  minimum rho in graphics (if not 0)
!     RGMAX:  maximum rho in graphics (if not 1)
!     NGLINE: maximum number of contour lines
!     NGRAPH: graphic mode: 0 for file out, 1 for contour plot

      DELT  = 1.D-2
      RIMPL = 1.D0
      EPSM  = 1.D-8
      EPSE  = 1.D-4
      LMAXE = 10
      EPSDE = 1.D-8
      H0DE  = 0.25D0
      PGMAX=10.D0
      RGMIN=0.D0
      RGMAX=1.D0
      NGLINE= 25
      NGRAPH= 1

!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
!     NTMAX   : maximum time step count
!     NTG1STEP: time step interval for storing global data
!     NTG1MIN:  minimum number of NTG1 save (initial allocation)
!     NTG1MAX:  maximum number of NTG1 save (thin out if exceeded)
!     NTG2STEP: time step interval for storing radial profile  data
!     NTG2MIN:  minimum number of NTG2 save (initial allocation)
!     NTG2MAX:  maximum number of NTG2 save (thin out if exceeded)
!     NTCLSTEP: time step interval for recalculating coefficients

      NTMAX    = 10
      NTG1STEP = 1
      NTG1MIN  = 101
      NTG1MAX  = 1001
      NTG2STEP = 1
      NTG2MIN  = 101
      NTG2MAX  = 1001
      NTCLSTEP = 1000

!-----------------------------------------------------------------------
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
!     MODELD: 0 : without radial transport
!           > 0 : with radial transport
!           > 10: with radial transport and pinch effect
! mod(MODELD,10)=1 : no p dependence
!                2 : 1/p dependence
!                3 : 1/sqrt(p) dependence
!                4 : 1/p^2 dependence
!                5 : deltaB/B stchastic diffusion
!     MODELD_boundary : 0 fix f at rho=1+DELR/2, namely FS2
!                    : 1 fix f at rho = 1, namely FS1. FS2 is variable.
!
!     MODELS : 0 No fusion reaction
!              1 DT reaction source (NSSPF,SPFTOT,SPFR0,SPFRW,SPFENG)
!              2 DT reaction source (self-consistent reactioin rate)
!              3 Using Legendre expansion in fusion reaction rate calculation
!     MODELW(ns): 0 for given diffusion coefficient model
!                 1 for wave E field calculated by WR(without beam radius)
!                 2 for wave E field calculated by WR(with beam radius)
!                 3 for given wave E field model
!                 4 for wave E field calculated by WM

      MODELE= 0
      MODELA= 0
      MODELC= 0
      MODELR= 0
      MODELD= 0
      MODELD_boundary= 0
      MODELS= 0
      DO NS=1,NSM
         MODELW(NS)=0
      END DO

      MODEL_KSP=5
      MODEL_PC =1

      MODEL_LOSS=0
      MODEL_synch=0
      MODEL_NBI=0
      MODEL_WAVE=0 ! 0=no wave calc., 1=wave calc.
!-----------------------------------------------------------------------
!     LLMAX : dimension of legendre polynomials's calculation
!     LLMAX_NF : dimension of legendre polynomials's calculation for nuclear fusion rate

      LLMAX = 2
      LLMAX_NF = 2
      EPSFP = 1.D-9
      LMAXFP = 10

!-----------------------------------------------------------------------
!     PWAVE : input power
!     LMAXNWR: max loop count in newton method to find ray position
!     EPSNWR: convergence criterion in newton method to find ray position

      PWAVE = 1.0D0
      LMAXNWR=100
      EPSNWR=1.D-6

!-----------------------------------------------------------------------
!     IDBGFP : debug graphic control parameter 
!          1 : Legendre polynomials
!          2 : fpl, M_l, N_l
!          4 : psy, phy and their derivatives
!          8 : dcpp, dctt, fcp

      IDBGFP=0

      NTG1M=0
      NTG2M=0

!-----------------------------------------------------------------------
!     Parameters relevant with DISRUPTION 
!
!     T0_quench : temperature after thermal quench [keV] at r=0
!     tau_quench: thermal quench time [sec]
!     tau_mgi   : MGI duration [sec]

      MODEL_DISRUPT=0 ! 0=no disruption, 1=disruption calc.
      MODEL_Conner_FP=0 ! runaway rate 0= Conner, 1=FP
      MODEL_BS=0 ! bootstrap current 0= off, 1=simple model
      MODEL_jfp=0 ! current evaluation 0= independent on f, 1=depend on f
      MODEL_LNL=0 ! Coulomb logarithm 0= variable w T , 1=fixed initial value, 2=fixed disrupted value
      MODEL_RE_pmax=0 ! RE non-RE boundary 0=NPMAX, 1=NPC_runaway
      MODELD_n_RE=0 ! radial transport of RE density 0=off, 1=on
      MODEL_IMPURITY=0 ! =1 MGI: satisfy quasi-neutrarity
      MODEL_SINK=0 ! 1= deltaB/B like sink term

      time_quench_start=0.D0
      RJPROF1=2.D0
      RJPROF2=2.D0
      v_RE=1.D0 ! RE velocity / VC

      T0_quench=2.D-2
      tau_quench=1.D-3
      tau_mgi=5.D-3
      
! IF MODEL_IMPURITY=1 
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

      use fpcomm, only: &
           RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
           NSMAX,PA,PZ,PZ0,PN,PNS, &
           PTPR,PTPP,PTS,PU,PUS, &
           PNITB,PTITB,PUITB, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           RHOMIN,QMIN,RHOITB,RHOEDG, &
           MODELG,MODELN,MODELQ,RHOGMN,RHOGMX, &
           KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF,KNAMEQ2, &
           MODEFR,MODEFW,IDEBUG, &
           NPMAX,NTHMAX,NRMAX,NAVMAX,NP2MAX,IMTX, &
           NTMAX,NTCLSTEP,LMAXE,NGLINE,NGRAPH,LMAXNWR, &
           MODELE,MODELA,MODELC,MODELR,MODELD,LLMAX,LLMAX_NF,IDBGFP, &
           NTG1STEP,NTG1MIN,NTG1MAX, &
           NTG2STEP,NTG2MIN,NTG2MAX, &
           DRR0,E0,R1,DELR1,RMIN,RMAX, &
           DEC,PEC1,PEC2,PEC3,PEC4,RFEC,DELYEC,DLH,PLH1,PLH2,RLH, &
           DFW,PFW1,PFW2,RFW,RFDW,DELNPR,CEWR,CEWTH,CEWPH, &
           RKWR,RKWTH,RKWPH,REWY,DREWY,FACTWM,PWAVE,EPSNWR, &
           ZEFF,DELT,RIMPL,EPSM,EPSE,EPSDE,H0DE, &
           nsamax,nsbmax, &
           ns_nsa,ns_nsb, &
           pmax,tloss,MODELW,MODELS,NBEAMMAX, &
           NSSPB,SPBTOT,SPBR0,SPBRW,SPBENG,SPBANG,SPBPANG,&
           NSSPF,SPFTOT,SPFR0,SPFRW,SPFENG,&
           LMAXFP, EPSFP,NCMIN,NCMAX, DRRS, MODEL_KSP, MODEL_PC, &
           N_partition_s, N_partition_r, N_partition_p, MODEL_DISRUPT, &
           MODEL_synch, MODEL_LOSS, MODEL_SINK, T0_quench, tau_quench, deltaB_B, &
           MODEL_NBI, MODEL_WAVE, MODEL_IMPURITY, MODEL_Conner_FP, MODEL_BS, MODEL_jfp, &
           MODEL_LNL, time_quench_start, MODEL_RE_pmax, RJPROF1, RJPROF2, MODELD_n_RE, &
           pmax_bb, v_RE, target_zeff, n_impu, SPITOT, MODELD_boundary, tau_mgi

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: nid
      INTEGER,INTENT(OUT) :: ist,ierr

      NAMELIST /FP/ &
           RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
           NSMAX,PA,PZ,PZ0,PN,PNS,PTPR,PTPP,PTS,PU,PUS, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           RHOMIN,QMIN,RHOITB,RHOEDG, &
           MODELG,MODELN,MODELQ,RHOGMN,RHOGMX, &
           KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF,KNAMEQ2, &
           MODEFR,MODEFW,IDEBUG, &
           NPMAX,NTHMAX,NRMAX,NAVMAX,NP2MAX,IMTX, &
           NTMAX,NTCLSTEP,LMAXE,NGLINE,NGRAPH,LMAXNWR, &
           MODELE,MODELA,MODELC,MODELW,MODELR,MODELD,LLMAX,LLMAX_NF,IDBGFP, &
           NTG1STEP,NTG1MIN,NTG1MAX, &
           NTG2STEP,NTG2MIN,NTG2MAX, &
           DRR0,E0,R1,DELR1,RMIN,RMAX, &
           DEC,PEC1,PEC2,PEC3,PEC4,RFEC,DELYEC,DLH,PLH1,PLH2,RLH, &
           DFW,PFW1,PFW2,RFW,RFDW,DELNPR,CEWR,CEWTH,CEWPH, &
           RKWR,RKWTH,RKWPH,REWY,DREWY,FACTWM,PWAVE,EPSNWR, &
           ZEFF,DELT,RIMPL,EPSM,EPSE,EPSDE,H0DE, &
           NSSPB,SPBTOT,SPBR0,SPBRW,SPBENG,SPBANG,SPBPANG, &
           NSSPF,SPFTOT,SPFR0,SPFRW,SPFENG, &
           pmax,tloss,LMAXFP,EPSFP,MODELS,NBEAMMAX, &
           nsamax,nsbmax,ns_nsa,ns_nsb,NCMIN,NCMAX,DRRS, MODEL_KSP, MODEL_PC, &
           N_partition_s, N_partition_r, N_partition_p, MODEL_DISRUPT, &
           MODEL_synch, MODEL_LOSS, MODEL_SINK, T0_quench, tau_quench, deltaB_B, &
           MODEL_NBI, MODEL_WAVE, MODEL_IMPURITY, MODEL_Conner_FP, MODEL_BS, MODEL_jfp, &
           MODEL_LNL, time_quench_start, MODEL_RE_pmax, RJPROF1, RJPROF2, MODELD_n_RE, &
           pmax_bb, v_RE, target_zeff, n_impu, SPITOT, MODELD_boundary, tau_mgi


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

      WRITE(6,*) '&FP : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'
      WRITE(6,*) '      NSMAX,PA,PZ,PZ0,PN,PNS,PTPR,PTPP,PTS,PU,PUS,'
      WRITE(6,*) '      PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'
      WRITE(6,*) '      RHOMIN,QMIN,RHOITB,RHOEDG,'
      WRITE(6,*) '      MODELG,MODELN,MODELQ,RHOGMN,RHOGMX,'
      WRITE(6,*) '      KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF,KNAMEQ2,'
      WRITE(6,*) '      MODEFR,MODEFW,IDEBUG,'
      WRITE(6,*) '      NPMAX,NTHMAX,NRMAX,NAVMAX,NP2MAX,IMTX,'
      WRITE(6,*) '      NTMAX,NTCLSTEP,LMAXE,NGLINE,NGRAPH,LMAXNWR,'
      WRITE(6,*) '      MODELE,MODELA,MODELC,MODELW,MODELR,MODELD,'
      WRITE(6,*) '      NTG1STEP,NTG1MIN,NTG1MAX,LLMAX,LLMAX_NF,IDBGFP,'
      WRITE(6,*) '      NTG2STEP,NTG2MIN,NTG2MAX,'
      WRITE(6,*) '      DRR0,E0,R1,DELR1,RMIN,RMAX,'
      WRITE(6,*) '      DEC,PEC1,PEC2,PEC3,PEC4,RFEC,DELYEC,DLH,PLH1,PLH2,RLH,'
      WRITE(6,*) '      DFW,PFW1,PFW2,RFW,RFDW,DELNPR,CEWR,CEWTH,CEWPH,'
      WRITE(6,*) '      RKWR,RKWTH,RKWPH,REWY,DREWY,FACTWM,PWAVE,EPSNWR,'
      WRITE(6,*) '      NSSPB,SPBTOT,SPBR0,SPBRW,SPBENG,SPBANG,SPBPANG,'
      WRITE(6,*) '      NSSPF,SPFTOT,SPFR0,SPFRW,SPFENG,'
      WRITE(6,*) '      ZEFF,DELT,RIMPL,EPSM,EPSE,EPSDE,H0DE,'
      WRITE(6,*) '      nsamax,nsbmax,ns_nsa,ns_nsb,pmax,tloss,'
      WRITE(6,*) '      MODELS,NBEAMMAX,DRRS,MODEL_KSP,MODEL_PC'
      WRITE(6,*) '      N_partition_s, N_partition_r, N_partition_p, MODEL_DISRUPT'
      WRITE(6,*) '      MODEL_synch, MODEL_loss, MODEL_SINK, T0_quench, tau_quench, deltaB_B,'
      WRITE(6,*) '      MODEL_NBI, MODEL_WAVE, MODEL_IMPURITY, MODEL_Conner_FP, MODEL_BS, MODEL_jfp'
      WRITE(6,*) '      MODEL_LNL, time_quench_start, MODEL_RE_pmax, RJPROF1, RJPROF2, MODELD_n_RE'
      WRITE(6,*) '      pmax_bb, v_RE, target_zeff, n_impu, SPITOT, MODELD_boundary, tau_mgi'


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

      USE fpcomm
      USE libmpi
      USE libmtx
      IMPLICIT NONE
      INTEGER,DIMENSION(99):: idata
      real(8),DIMENSION(99):: rdata
      complex(8),DIMENSION(3):: cdata

!----- PL input parameters -----     

      idata( 1)=NSMAX
      idata( 2)=MODELG
      idata( 3)=MODELN
      idata( 4)=MODELQ
      idata( 5)=MODEFR
      idata( 6)=MODEFW
      idata( 7)=IDEBUG

      CALL mtx_broadcast_integer(idata,7)
      NSMAX =idata( 1)
      MODELG=idata( 2)
      MODELN=idata( 3)
      MODELQ=idata( 4)
      MODEFR=idata( 5)
      MODEFW=idata( 6)
      IDEBUG=idata( 7)

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
      rdata(23)=RJPROF1
      rdata(24)=RJPROF2
      CALL mtx_broadcast_real8(rdata,24)
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
      RJPROF1=rdata(23)
      RJPROF2=rdata(24)

      CALL mtx_broadcast_real8(PA,NSMAX)
      CALL mtx_broadcast_real8(PZ,NSMAX)
      CALL mtx_broadcast_real8(PZ0,NSMAX)
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

      CALL mtx_broadcast_character(KNAMEQ,80)
      CALL mtx_broadcast_character(KNAMWR,80)
      CALL mtx_broadcast_character(KNAMFP,80)
      CALL mtx_broadcast_character(KNAMWM,80)
      CALL mtx_broadcast_character(KNAMPF,80)
      CALL mtx_broadcast_character(KNAMFO,80)
      CALL mtx_broadcast_character(KNAMTR,80)
      CALL mtx_broadcast_character(KNAMEQ2,80)

!----- FP input parameters -----

      idata( 1)=NPMAX
      idata( 2)=NTHMAX
      idata( 3)=NRMAX
      idata( 4)=NTMAX
      idata( 5)=NTG1STEP
      idata( 6)=NTG1MIN
      idata( 7)=NTG1MAX
      idata( 8)=NTG2STEP
      idata( 9)=NTG2MIN
      idata(10)=NTG2MAX
      idata(11)=NTCLSTEP
      idata(12)=NSAMAX
      idata(13)=NSBMAX

      idata(14)=MODELE
      idata(15)=MODELA
      idata(16)=MODELC
      idata(17)=MODELR
      idata(18)=MODELD

      idata(19)=LLMAX
      idata(20)=NAVMAX
      idata(21)=LMAXFP
      idata(22)=LMAXE
      idata(23)=LMAXNWR
      idata(24)=NGRAPH
      idata(25)=NGLINE
      idata(26)=IMTX
      idata(27)=IDBGFP
      idata(28)=NTEST
      idata(29)=NP2MAX
      idata(30)=NBEAMMAX
      idata(31)=NSSPF
      idata(32)=MODELS
      idata(33)=MODEL_KSP
      idata(34)=MODEL_PC
      idata(35)=N_partition_s
      idata(36)=N_partition_r
      idata(37)=N_partition_p
      idata(38)=MODEL_DISRUPT
      idata(39)=MODEL_synch
      idata(40)=MODEL_loss
      idata(41)=MODEL_NBI
      idata(42)=MODEL_IMPURITY
      idata(43)=MODEL_Conner_FP
      idata(44)=MODEL_BS
      idata(45)=MODEL_jfp
      idata(46)=MODEL_LNL
      idata(47)=MODEL_RE_pmax
      idata(48)=MODEL_WAVE
      idata(49)=MODEL_SINK
      idata(50)=MODELD_n_RE
      idata(51)=LLMAX_NF
      idata(52)=n_impu
      idata(53)=MODELD_boundary

      CALL mtx_broadcast_integer(idata,53)
      NPMAX   =idata( 1)
      NTHMAX  =idata( 2)
      NRMAX   =idata( 3)
      NTMAX   =idata( 4)
      NTG1STEP=idata( 5)
      NTG1MIN =idata( 6)
      NTG1MAX =idata( 7)
      NTG2STEP=idata( 8)
      NTG2MIN =idata( 9)
      NTG2MAX =idata(10)
      NTCLSTEP=idata(11)
      NSAMAX  =idata(12)
      NSBMAX  =idata(13)

      MODELE  =idata(14)
      MODELA  =idata(15)
      MODELC  =idata(16)
      MODELR  =idata(17)
      MODELD  =idata(18)

      LLMAX   =idata(19)
      NAVMAX  =idata(20)
      LMAXFP  =idata(21)
      LMAXE   =idata(22)
      LMAXNWR =idata(23)
      NGRAPH  =idata(24)
      NGLINE  =idata(25)
      IMTX    =idata(26)
      IDBGFP  =idata(27)
      NTEST   =idata(28)
      NP2MAX  =idata(29)
      NBEAMMAX=idata(30)
      NSSPF   =idata(31)
      MODELS  =idata(32)
      MODEL_KSP=idata(33)
      MODEL_PC =idata(34)
      N_partition_s = idata(35)
      N_partition_r = idata(36)
      N_partition_p = idata(37)
      MODEL_DISRUPT = idata(38)
      MODEL_synch = idata(39)
      MODEL_loss  = idata(40)
      MODEL_NBI  = idata(41)
      MODEL_IMPURITY  = idata(42)
      MODEL_Conner_FP  = idata(43)
      MODEL_BS  = idata(44)
      MODEL_jfp  = idata(45)
      MODEL_LNL  = idata(46)
      MODEL_RE_pmax  = idata(47)
      MODEL_WAVE = idata(48)
      MODEL_SINK = idata(49)
      MODELD_n_RE= idata(50)
      LLMAX_NF   = idata(51)
      n_impu = idata(52)
      MODELD_boundary = idata(53)

      CALL mtx_broadcast_integer(NS_NSA,NSAMAX)
      CALL mtx_broadcast_integer(NS_NSB,NSBMAX)
      CALL mtx_broadcast_integer(MODELW,NSMAX)
      CALL mtx_broadcast_integer(NCMIN,NSMAX)
      CALL mtx_broadcast_integer(NCMAX,NSMAX)
      CALL mtx_broadcast_integer(NSSPB,NBEAMMAX)

      rdata( 1)=DELT
      rdata( 2)=RMIN
      rdata( 3)=RMAX
      rdata( 4)=R1
      rdata( 5)=DELR1
      rdata( 6)=DRR0
      rdata( 7)=E0
      rdata( 8)=DEC
      rdata( 9)=PEC1
      rdata(10)=PEC2
      rdata(11)=RFEC
      rdata(12)=DELYEC
      rdata(13)=DLH
      rdata(14)=PLH1
      rdata(15)=PLH2
      rdata(16)=RLH
      rdata(17)=DFW
      rdata(18)=PFW1
      rdata(19)=PFW2
      rdata(20)=RFW
      rdata(21)=RFDW
      rdata(22)=DELNPR
      rdata(23)=RKWTH
      rdata(24)=RKWPH
      rdata(25)=REWY
      rdata(26)=DREWY
      rdata(27)=FACTWM
      rdata(28)=RIMPL
      rdata(29)=EPSFP
      rdata(30)=EPSE
      rdata(31)=EPSDE
      rdata(32)=H0DE
      rdata(33)=EPSNWR
      rdata(34)=EPSM
      rdata(35)=PGMAX
      rdata(36)=RGMAX
      rdata(37)=RGMIN
      rdata(38)=SPFTOT
      rdata(39)=SPFR0
      rdata(40)=SPFRW
      rdata(41)=SPFENG
      rdata(42)=DRRS
      rdata(43)=PEC3
      rdata(44)=PEC4
      rdata(45)=T0_quench
      rdata(46)=tau_quench
      rdata(47)=time_quench_start
      rdata(48)=deltaB_B
      rdata(49)=v_RE
      rdata(50)=target_zeff
      rdata(51)=SPITOT
      rdata(52)=tau_mgi

      CALL mtx_broadcast_real8(rdata,52)
      DELT  =rdata( 1)
      RMIN  =rdata( 2)
      RMAX  =rdata( 3)
      R1    =rdata( 4)
      DELR1 =rdata( 5)
      DRR0  =rdata( 6)
      E0    =rdata( 7)
      DEC   =rdata( 8)
      PEC1  =rdata( 9)
      PEC2  =rdata(10)
      RFEC  =rdata(11)
      DELYEC=rdata(12)
      DLH   =rdata(13)
      PLH1  =rdata(14)
      PLH2  =rdata(15)
      RLH   =rdata(16)
      DFW   =rdata(17)
      PFW1  =rdata(18)
      PFW2  =rdata(19)
      RFW   =rdata(20)
      RFDW  =rdata(21)
      DELNPR=rdata(22)
      RKWTH =rdata(23)
      RKWPH =rdata(24)
      REWY  =rdata(25)
      DREWY =rdata(26)
      FACTWM=rdata(27)
      RIMPL =rdata(28)
      EPSFP =rdata(29)
      EPSE  =rdata(30)
      EPSDE =rdata(31)
      H0DE  =rdata(32)
      EPSNWR=rdata(33)
      EPSM  =rdata(34)
      PGMAX =rdata(35)
      RGMAX =rdata(36)
      RGMIN =rdata(37)
      SPFTOT=rdata(38)
      SPFR0 =rdata(39)
      SPFRW =rdata(40)
      SPFENG=rdata(41)
      DRRS  =rdata(42)
      PEC3  =rdata(43)
      PEC4  =rdata(44)
      T0_quench=rdata(45)
      tau_quench=rdata(46)
      time_quench_start=rdata(47)
      deltaB_B = rdata(48)
      v_RE= rdata(49)
      target_zeff= rdata(50)
      SPITOT= rdata(51)
      tau_mgi= rdata(52)

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

      use fpcomm, only: &
           RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
           NSMAX,PA,PZ,PZ0,PN,PNS,PTPR,PTPP,PTS,PU,PUS, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           RHOMIN,QMIN,RHOITB,RHOEDG, &
           MODELG,MODELN,MODELQ,RHOGMN,RHOGMX, &
           KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF,KNAMEQ2, &
           MODEFR,MODEFW,IDEBUG, &
           NPMAX,NTHMAX,NRMAX,NAVMAX,NP2MAX,IMTX, &
           NTMAX,NTCLSTEP,LMAXE,NGLINE,NGRAPH,LMAXNWR, &
           MODELE,MODELA,MODELC,MODELW,MODELR,MODELD,LLMAX,LLMAX_NF,IDBGFP, &
           NTG1STEP,NTG1MIN,NTG1MAX, &
           NTG2STEP,NTG2MIN,NTG2MAX, &
           DRR0,E0,R1,DELR1,RMIN,RMAX, &
           DEC,PEC1,PEC2,PEC3,PEC4,RFEC,DELYEC,DLH,PLH1,PLH2,RLH, &
           DFW,PFW1,PFW2,RFW,RFDW,DELNPR,CEWR,CEWTH,CEWPH, &
           RKWR,RKWTH,RKWPH,REWY,DREWY,FACTWM,PWAVE,EPSNWR, &
           NSSPB,SPBTOT,SPBR0,SPBRW,SPBENG,SPBANG,SPBPANG,&
           NSSPF,SPFTOT,SPFR0,SPFRW,SPFENG,&
           ZEFF,DELT,RIMPL,EPSM,EPSE,EPSDE,H0DE, &
           nsamax,nsbmax,ns_nsa,ns_nsb,pmax,tloss,MODELS,NCMIN,NCMAX, &
           nbeammax,DRRS,MODEL_KSP,MODEL_PC,N_partition_s,N_partition_r,N_partition_p, &
           nsize, MODEL_DISRUPT, MODEL_synch, MODEL_LOSS, MODEL_SINK, T0_quench, tau_quench, deltaB_B, &
           MODEL_NBI, MODEL_WAVE, MODEL_IMPURITY, MODEL_Conner_FP, MODEL_BS, MODEL_jfp, MODEL_LNL, &
           time_quench_start, MODEL_RE_pmax, RJPROF1, RJPROF2, MODELD_n_RE, pmax_bb, v_RE, &
           target_zeff, n_impu, SPITOT, MODELD_boundary, tau_mgi

      IMPLICIT NONE
      integer:: nsa,nsb,ns,NBEAM

      WRITE(6,600) 'E0      ',E0      ,'DRR0    ',DRR0    ,'ZEFF    ',ZEFF
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
         WRITE(6,600) 'pmax    ',pmax(nsb),'tloss   ',tloss(nsb)
      END DO

      WRITE(*,*) "----- PARAMETERS OF HEATINGS -----"
      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         WRITE(6,'(A,I2,A,I2)') 'NSA = ',NSA,'  NS  = ',NS
         
         IF(MODELW(NS).EQ.0) THEN
            
            WRITE(6,600) 'DEC     ',DEC     ,'RFEC    ',RFEC  
            WRITE(6,600) 'PEC1    ',PEC1    ,'PEC2    ',PEC2    ,'DELYEC  ',DELYEC
            WRITE(6,600) 'PEC3    ',PEC3    ,'PEC4    ',PEC4
            WRITE(6,600) 'DLH     ',DLH     ,'RLH     ',RLH
            WRITE(6,600) 'PLH1    ',PLH1    ,'PLH2    ',PLH2
            WRITE(6,600) 'DFW     ',DFW     ,'RFW     ',RFW
            WRITE(6,600) 'PFW1    ',PFW1    ,'PFW2    ',PFW2
            
         ELSEIF(MODELW(NS).EQ.1) THEN
            WRITE(6,600) 'RFDW    ',RFDW    ,'DELNPR  ',DELNPR
            WRITE(6,600) 'PWAVE   ',PWAVE   ,'DELYEC  ',DELYEC
            WRITE(6,602) 'EPSNWR  ',EPSNWR  ,'LMAXNW  ',LMAXNWR
            
         ELSEIF(MODELW(NS).EQ.2) THEN
            WRITE(6,600) 'RFDW    ',RFDW    ,'DELNPR  ',DELNPR
            WRITE(6,600) 'PWAVE   ',PWAVE   ,'DELYEC  ',DELYEC
            WRITE(6,602) 'EPSNWR  ',EPSNWR  ,'LMAXNW  ',LMAXNWR
            
         ELSEIF(MODELW(NS).EQ.3) THEN
            WRITE(6,600) 'RFDW    ',RFDW    ,'DELNPR  ',DELNPR
            WRITE(6,600) 'CEWR/R  ',DBLE(CEWR) ,'CEWR/I  ',DIMAG(CEWR)
            WRITE(6,600) 'CEWTH/R ',DBLE(CEWTH),'CEWTH/I ',DIMAG(CEWTH)
            WRITE(6,600) 'CEWPH/R ',DBLE(CEWPH),'CEWPH/I ',DIMAG(CEWPH)
            WRITE(6,600) 'REWY    ',REWY    ,'DREWY   ',DREWY
            WRITE(6,600) 'RKWR    ',RKWR    ,'RKWTH   ',RKWTH   ,'RKWPH   ',RKWPH
            
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
                         'SPBANG   ',SPBANG(NBEAM), &
                         'SPBPANG   ',SPBPANG(NBEAM)
            WRITE(6,600) 'SPBR0    ',SPBR0(NBEAM), &
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

      WRITE(6,*) "----- THE OTHER PARAMETERS -----"
      WRITE(6,600) 'RIMPL   ',RIMPL   ,'EPSM    ',EPSM
      WRITE(6,600) 'EPSDE   ',EPSDE   ,'H0DE    ',H0DE    ,'EPSE    ',EPSE
      WRITE(6,603) 'LMAXE   ',LMAXE   ,'LLMAX   ',LLMAX   ,'NGLINE  ',NGLINE
      WRITE(6,603) 'IDBGFP  ',IDBGFP  ,'NGRAPH  ',NGRAPH  ,'LLMAX_NF',LLMAX_NF

      WRITE(6,603) 'NPMAX   ',NPMAX   ,'NTHMAX  ',NTHMAX  ,'NRMAX   ',NRMAX
      WRITE(6,603) 'NAVMAX  ',NAVMAX  ,'NP2MAX  ',NP2MAX  ,'NTMAX   ',NTMAX   
      WRITE(6,603) 'NTG1STEP',NTG1STEP,'NTG1MIN ',NTG1MIN ,'NTG1MAX ',NTG1MAX
      WRITE(6,603) 'NTG2STEP',NTG2STEP,'NTG2MIN ',NTG2MIN ,'NTG2MAX ',NTG2MAX
      WRITE(6,603) 'MODELE  ',MODELE  ,'MODELA  ',MODELA  ,'MODELC  ',MODELC
      WRITE(6,603) 'MODEFR  ',MODEFR  ,'MODELD  ',MODELD  ,'MODEFW  ',MODEFW
      WRITE(6,603) 'IMTX    ',IMTX    ,'MODEL_KSP',MODEL_KSP,'MODEL_PC ',MODEL_PC

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

  600 FORMAT(1H ,A8,'=',1PE12.4:3X,A8,'=',1PE12.4:3X,A8,'=',1PE12.4)
  601 FORMAT(1H ,A8,'=',1PE12.4:3X,A8,'=',1PE12.4:3X,A8,'=',I8)
  602 FORMAT(1H ,A8,'=',1PE12.4:3X,A8,'=',I8,4X  :3X,A8,'=',I8)
  603 FORMAT(1H ,A8,'=',I8,4X  :3X,A8,'=',I8,4X  :3X,A8,'=',I8)
    END SUBROUTINE fp_view

  END module fpinit
