MODULE trinit

  PUBLIC tr_init  ! initialize input parameters
  PUBLIC tr_term  ! terminate procedure
  PUBLIC tr_parm  ! parameter input through namelist
  PUBLIC tr_view  ! show values of input parameters

  PRIVATE

  INTERFACE 
     SUBROUTINE task_klin(line,kid,mode,xxparm)
       CHARACTER(LEN=80),INTENT(OUT):: line
       CHARACTER(LEN=1),INTENT(OUT):: kid
       INTEGER(4),INTENT(OUT):: mode
       EXTERNAL:: xxparm
     END SUBROUTINE task_klin
     SUBROUTINE task_parm(mode,kwd,kin,xxnlin,xxplist,ierr)
       INTEGER(4),INTENT(IN):: mode
       CHARACTER(LEN=*),INTENT(IN):: kwd
       CHARACTER(LEN=*),INTENT(IN):: kin
       EXTERNAL:: xxnlin
       EXTERNAL:: xxplist
       INTEGER(4),INTENT(OUT):: ierr
     END SUBROUTINE task_parm
  END INTERFACE

CONTAINS

! ***** initialize  input parameters

  SUBROUTINE tr_init

    USE plcomm
    USE trcomm, ONLY: &
           nrmax,ntmax,dt,rg_fixed,nsamax,ns_nsa, &
           phia,pa_mion,pz_mion,pa_mimp,pz_mimp,      &
           lmaxtr,epsltr,mdltr_nc,mdltr_tb,mdltr_prv, &
           mdluf,mdlxp,mdlni,mdler,modelg,nteqit, &
           time_slc,time_snap,mdlugt, &
           dtr0,dtr1,ltcr,ph0,phs,dprv1,dprv2,rhog_prv, &
           cdtrn,cdtru,cdtrt, &
           ntstep,ngtmax,ngtstp,rips,ripe,profj1,profj2, &
           pnb_tot,pnb_eng,pnb_rw,pnb_r0, &
           ufid_bin,kuf_dir,kuf_dev,kuf_dcg, &
           mdleqb,mdleqn,mdlequ,mdleqt, &
           mdlijq,mdlgmt,mdlsrc,mdlglb, &
           mdlwrt,nwrstp,kwpnam,kwtnam

    USE plinit
    IMPLICIT NONE
    INTEGER(ikind):: nsa,ns

    OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
    CALL pl_init

!     ======( DEVICE PARAMETERS )======
!       ***  variables declared in PL ***
!        RR    : Plasma major radius                             (m)
!        RA    : Plasma minor radius                             (m)
!        RB    : Wall minor radius                               (m)
!        RKAP  : Plasma shape elongation
!        RDLT  : Plasma shape triangularity
!        BB    : Magnetic field at center                        (T)
!        Q0    : Safety factor at center
!        QA    : Safety factor on plasma surface
!        RIP   : Plasma current                                 (MA)
!        PROFJ : Curren density profile parameter (power of (1 - rho^2))  

!       *** variables declared in TR ***
!        PHIA  : total toroidal flux enclosed by the plasma (Wb)
!        RIPS  : toroidal current at the beginning [MA]
!        RIPE  : toroidal current at the end       [MA]
!

      RR    = 3.D0
      RA    = 1.D0
      RB    = 1.2D0
      RKAP  = 1.D0
      RDLT  = 0.D0
      BB    = 3.D0
      Q0    = 1.D0
      QA    = 3.D0
      RIP   = 3.D0
      PROFJ = 2.D0

      PHIA = 0.D0
      RIPS = 1.D0
      RIPE = 1.D0

!     ======( PLASMA PARAMETERS )======

!        NSMAX : Number of particle species
!        PA    : Mass number
!        PZ    : Charge number
!        PZ0   : Atomic number (-1 for electron)
!        PN    : Density at center                     (1.0E20/m**3)
!        PNS   : Density on plasma surface             (1.0E20/m**3)
!        PTPR  : Parallel temperature at center                (keV)
!        PTPP  : Perpendicular temperature at center           (keV)
!        PTS   : Temperature on surface                        (keV)
!        PU    : Toroidal rotation velocity at center          (m/s)
!        PUS   : Toroidal rotation velocity on surface         (m/s)
!        PNITB : Density increment at ITB              (1.0E20/Mm*3)
!        PTITB : Temperature increment at ITB                  (keV)
!        PUITB : Toroidal rotation velocity increment at ITB   (m/s)

!        KIDNS : index of particle species
!        IDION :  1 = fast ion particle
!                 0 = else                                                     

    ! *** Species to be condidered are defined in PL     ***
    ! *** This part is only for initilization of profile ***

    NS = 1 ! electron
    PN(NS)   = 0.5D0
    PNS(NS)  = 0.2D0
    PTPR(NS) = 1.50D0
    PTPP(NS) = 1.50D0
    PTS(NS)  = 0.05D0
    PU(NS)   = 0.D0
    PUS(NS)  = 0.D0
    PNITB(NS)= 0.D0
    PTITB(NS)= 0.D0
    PUITB(NS)= 0.D0

    NS = 2 ! D (bulk)
    PN(NS)   = 0.5D0
    PNS(NS)  = 0.2D0
    PTPR(NS) = 1.50D0
    PTPP(NS) = 1.50D0
    PTS(NS)  = 0.05D0
    PU(NS)   = 0.D0
    PUS(NS)  = 0.D0
    PNITB(NS)= 0.D0
    PTITB(NS)= 0.D0
    PUITB(NS)= 0.D0
    
    NS = 3 ! D (fast)
    PN(NS)   = 0.D0
    PNS(NS)  = 0.D0
    PTPR(NS) = 5.0D0
    PTPP(NS) = 5.0D0
    PTS(NS)  = 0.05D0
    PU(NS)   = 0.D0
    PUS(NS)  = 0.D0
    PNITB(NS)= 0.D0
    PTITB(NS)= 0.D0
    PUITB(NS)= 0.D0

    NS = 4 ! C (impurity)
    PN(NS)   = 0.0D0
    PNS(NS)  = 0.0D0
    PTPR(NS) = 5.0D0
    PTPP(NS) = 5.0D0
    PTS(NS)  = 0.05D0
    PU(NS)   = 0.D0
    PUS(NS)  = 0.D0
    PNITB(NS)= 0.D0
    PTITB(NS)= 0.D0
    PUITB(NS)= 0.D0


!     ======( PROFILE PARAMETERS )======

!        PROFN1: Density profile parameter (power of rho)
!        PROFN2: Density profile parameter (power of (1 - rho^PROFN1))
!        PROFT1: Temperature profile parameter (power of rho)
!        PROFT2: Temperature profile parameter (power of (1 - rho^PROFN1))
!        PROFU1: Rotation profile parameter (power of rho)
!        PROFU2: Rotation profile parameter (power of (1 - rho^PROFN1))
!        PROFJ1: Current density profile parameter (power of rho)
!        PROFJ2: Current density profile parameter (power of (1 - rho^PROFJ1))

    PROFN1 = 2.D0
    PROFN1 = 1.5D0
    PROFU1 = 2.D0
    PROFU1 = 1.D0
    PROFT1 = 2.D0
    PROFT2 = 1.D0
    PROFJ1 = 2.D0
    PROFJ2 = 1.D0

!     ======( MODEL PARAMETERS )======

!        MODELG: Control plasma geometry model
!                   0: Slab geometry
!                   1: Cylindrical geometry
!                   2: Toroidal geometry
!                   3: TASK/EQ output geometry
!                   4: VMEC output geometry
!                   5: EQDSK output geometry
!                   6: Boozer output geometry
!                   7: new VMEC output geometry
!                   8: call TOPICS/EQU
!                   9: call TASK/EQ
!
!        NTEQIT: Step interval of EQ calulation
!                  0 : Initial equilibrium only


!        MODELN: Control plasma profile
!                   0: Calculated from PN,PNS,PTPR,PTPP,PTS,PU,PUS; 0 in SOL
!                   1: Calculated from PN,PNS,PTPR,PTPP,PTS,PU,PUS; PNS in SOL
!                   7: Read from file by means of WMDPRF routine (DIII-D)
!                   8: Read from file by means of WMXPRF routine (JT-60)
!                   9: Read from file KNAMTR (TASK/TR)
!        MODELQ: Control safety factor profile (for MODELG=0,1,2)
!                   0: Parabolic q profile (Q0,QA,RHOMIN,RHOITB)
!                   1: Given current profile (RIP,PROFJ)
!        RHOMIN: rho at minimum q (0 for positive shear)
!        QMIN  : q minimum for reversed shear
!        RHOITB: rho at ITB (0 for no ITB)
!        RHOEDG: rho at EDGE for smoothing (1 for no smooth)
!

    MODELG = 0 ! For now, 'MDLGMT' is used instead.
    NTEQIT = 0

!     ======( GRAPHIC PARAMETERS )======

!        RHOGMN: minimum rho in radial profile
!        RHOGMX: maximum rho in radial profile

!     ======( IO FILE NAMES )======

!        KNAMEQ: Filename of equilibrium data
!        KNAMWR: Filename of ray tracing data
!        KNAMWM: Filename of full wave data
!        KNAMFP: Filename of Fokker-Planck data
!        KNAMFO: Filename of File output
!        KNAMPF: Filename of profile data

!     ======( FILE IO MODES )======

!        MODEFR: File name interaction at reading file
!                 0 : WITHOUT PROMPT
!                 1 : WITH FILE NAME INPUT
!        MODEFW: File name interaction at writing file
!                 0 : WITHOUT PROMPT, ALWAYS OVERWRITE
!                 1 : WITHOUT PROMPT, CONFIRM, IF FILE EXISTS
!                 2 : WITHOUT PROMPT, ASK NEW NAME, IF FILE EXISTS
!                 3 : WITHOUT PROMPT, ERROR, IF FILE EXISTS
!                 4 : WITH FILE NAME INPUT, ALWAYS OVERWRITE
!                 5 : WITH FILE NAME INPUT, CONFIRM, IF FILE EXISTS
!                 6 : WITH FILE NAME INPUT, ASK NEW NAME, IF FILE EXISTS
!                 7 : WITH FILE NAME INPUT, ERROR, IF FILE EXISTS


! =========================================================================
!     Variables above is declared in TASK/PL
! =========================================================================
! =========================================================================


!     ==== TR PARAMETERS ====

!     ==== Eqs. Selection Parameter ====
    MDLEQB = 1  ! 0/1 for B_theta
    MDLEQN = 0  ! 0/1 for density
    MDLEQU = 0  ! 0/1 for rotation
    MDLEQT = 1  ! 0/1 for heat

!        NRMAX  : NUMBER OF RADIAL MESH POINTS
!        NTMAX  : NUMBER OF TIME STEP
!        DT     : SIZE OF TIME STEP

!        nsamax   : number of active particle species
!        rg_fixed(3,nsm) : minimum radius of fixed profile
!        nitmax   : maximum number of iterations
!        epsltr   : tolerance of iteration
!
!        pa_mion  : atomic number of main hydrogenic ion
!        pz_mion  : charge number of main hydrogenic ion
!        pa_mimp  : atomic number of main impurity ion
!        pz_mion  : charge number of main impurity ion

    nrmax    = 50
    ntmax    = 100
    dt       = 0.01D0
    rg_fixed(1:3,0:nsm) = rb

    nsamax   = 3

    pa_mion = 2.d0  ! Deuterium
    pz_mion = 1.d0 
    pa_mimp = 6.d0  ! Carbon
    pz_mimp = 12.d0

!     ==== Convergence Parameter ====
!        epsltr : convergence criterion of iteration
!        lmaxtr : maximum count of iteration

    epsltr = 1.D-6
    lmaxtr = 100


!       mdltr_nc  = 0 : no neoclassical transport
!                   1 : NCLASS trasport model
!       mdltr_tb  = 0 : no turbulent transport
!                   1 : constant diffusion
!                   2 : STIFF MODEL (Pereverzev)
!       dtr0        : lower diffusion coefficient for simple turbulent model
!       dtr1        : upper diffusion coefficnent for simple turbulent model
!       ltcr      : critical scale length [m]
!
!       cdtrn     : factor for particle diffusivity
!       cdtru     : factor for toroidal viscosity
!       cdtrt     : factor for thermal diffusivity
!
!     ==== TR PARAMETERS for stiff modeling by Ikari  ====
!
!       mdltr_prv = 0 : no Pereverzev method
!                   1 : Pereverzev method applied : Denh=dprv1
!                   2 : Pereverzev method applied : Denh=dprv2*Dorg
!                   3 : Pereverzev method applied : Denh=dprv2*Dorg+dprv1
!                   4 : Pereverzev method applied : Denh=MIN(dprv2*Dorg,dprv1)
!                   5 : Pereverzev method applied : Denh=MAX(dprv2*Dorg,dprv1)
!       dprv1     : enhanced diffusion coefficient for Pereverzev mothod
!       dprv2     : diffusion enhancement factor for Pereverzev mothod
!       rhog_prv  : enhanced diffusion region (rhog(nr) > rhog_prv)

    mdltr_nc  = 1
    mdltr_tb  = 1
    mdltr_prv = 0
    dtr0  = 0.1D0
    dtr1  = 0.5D0
    ltcr  = 1.D0

    cdtrn = 1.D0
    cdtru = 1.D0
    cdtrt = 1.D0

    dprv1     = 1.D0
    dprv2     = 1.D0
    rhog_prv  = 0.D0

!   === Auxiliary heating ===
!       ph0       : heating power density [MW/m^3] at r = 0
!       phs       : heating power density [MW/m^3] at r = a

    ph0   = 0.1D0
    phs   = 0.1D0

    pnb_tot = 0.d0  ! [MW]
    pnb_r0  = 0.3d0
    pnb_rw  = 0.5d0
    pnb_eng = 80.d0 ! [keV]

!   === Model switches ===
!
!       MDLIJQ : Control how to create d psi/d rho profile
!            1 : create from jtot using RIP as boundary condition
!            2 : create from qp using RIP as boundary condition
!            3 : create from jtot not using RIP
!            4 : create from qp not using RIP
!        * MDLIJQ = 3, 4 is especially for using exp. data
!
!        MDLER : the type of radial electric field
!            0 : pressure gradient only (nabla p)
!            1 : nabla p + toroidal rotation (V_tor)
!            2 : nabla p + V_tor + poloidal rotation (V_pol)

    MDLIJQ = 1
    MDLER  = 2

!     ==== Input from experimental data ==== 
!     MDLUF :
!     0 : not used
!     1 : read exp. data only for initial condition (default t=0.d0: time_slc)
!     2 : read exp. data successively in time evolution
!     3 : compared with TOPICS
!
!        MDLXP : Experimental data format
!            0 : UFILEs
!         else : MDSplus
!
!     UFID_BIN : Parameter which determines how to handle UFILEs.
!            0 : Binary files are loaded if available, or ASCII files
!                 are loaded and aftermath binary files are created.
!            1 : Only binary files are loaded.
!            2 : Only ASCII files are loaded and binary files are NOT
!                  created.
!
!        MDLNI : Switch how to determine main ion density, impurity density
!                 and effective charge number                             
!            1 : complete n_i and n_imp  from Zeff, n_e (and n_bulk)
!            2 : complete n_imp and Zeff from n_e, n_i (and n_bulk)
!            3 : complete n_i and Zeff   from n_e, n_imp (and n_bulk)    
!
!       MDLGMT :
!            0 : Slab geometry
!            1 : Cylindrical geometry
!            2 : Toroidal geometry
!            3 : TASK/EQ output geometry
!            4 : none
!            5 : none
!            6 : read experimental data only for initial profile (mdluf=1)
!            7 : read experimental data                          (mdluf=2)
!            8 : call TASK/EQ only for initial profile
!            9 : call TASK/EQ
!
!       MDLSRC :
!            1 : simple source model
!            6 : read experimental data only for initial profile (mdluf=1)
!            7 : read experimental data                          (mdluf=2)
!
!       MDLGLB :
!            1 : setup global variables by trinit
!            6 : read experimental data only for initial profile (mdluf=1)
!            7 : read experimental data                          (mdluf=2)
!
      mdluf     = 0

      mdlxp     = 0
      ufid_bin  = 0
      mdlni     = 1

      mdlgmt = 1
      mdlsrc = 1
      mdlglb = 1

!     ==== graphic output of experimental data ====
!      MDLUGT : set the time of snap shot
!           0 : --- from standard input every time graphic pages are opened.
!           1 : --- by 'time_snap' in namelist input (trparm)
!           2 : --- to the lastest time of the data
!       * This switch is valid only in the case of MDLUF = 2, 3

      mdlugt    = 0
      time_snap = -1.d0 ! not determined

!     ==== DEVICE NAME AND SHOT NUMBER IN UFILE DATA ====
!        KUF_DIR : UFILE database directory
!        KUF_DEV : Device name
!        KUF_DCG : Discharge number

!      kuf_dir = '../../../profiledb/profile_data'
      kuf_dir = '../../../profiledb/itpa_pr08'
      kuf_dev = 'd3d'
      kuf_dcg = '103818'

!     ==== TR PARAMETERS for data saving  ====
!        ntstep    : number of time step for status report
!        ngtmax    : maximum number of saved data
!        ngtstp    : number of time step for data save

    ntstep =    10
    ngtmax = 10001
    ngtstp =     1

!     ==== output raw data in csv format ===
!      MDLWRT : switch for the data output in csv format
!           0 : none
!        else : csv output
!      KWPNAM : csv file name for radial profiles
!      KWTNAM : csv file name for time evolution data
!      NWRSTP : number of time step for csv output

      mdlwrt = 0
      kwpnam = 'tr_rad.csv'
      kwtnam = 'tr_time.csv'

      nwrstp = 1

    RETURN
  END SUBROUTINE tr_init


!     ***** Terminalte procedure *****

  SUBROUTINE tr_term

    CLOSE(7)
    RETURN
  END SUBROUTINE tr_term


!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE tr_plist

    IMPLICIT NONE
    WRITE(6,601)
    RETURN

  601 FORMAT(' ','# &TR : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'/ &
             9X,'NSMAX,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,'/ &
             9X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'/ &
             9X,'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,'/ &
             9X,'MODELG,MODELN,MODELQ,RHOGMN,RHOGMX,'/ &
             9X,'KNAMEQ,KNAMWR,KNAMFP,KNAMFO,IDEBUG'/ &
             9X,'MODEFW,MODEFR'/ &
             9X,'nrmax,ntmax,dt,rg_fixed,nsamax,ns_nsa'/ &
             9X,'lmaxtr,epsltr,mdltr_nc,mdltr_tb,mdltr_prv,'/ &
             9X,'d0,d1,ltcr,ph0,phs,dprv1,dprv2,cdtrn,cdtru,cdtrt,'/ &
             9X,'ngtmax,ngtstep'/ &
             9X,'rips,ripe')
  END SUBROUTINE tr_plist

! ***** namelist input *****

  SUBROUTINE tr_nlin(nid,ist,ierr)

    USE plcomm
    USE trcomm, ONLY: &
           nrmax,ntmax,dt,rg_fixed,nsamax,ns_nsa, &
           ntstep,ngtmax,ngtstp, &
           phia,pa_mion,pz_mion,pa_mimp,pz_mimp,       &
           lmaxtr,epsltr,mdltr_nc,mdltr_tb,mdltr_prv,mdler,&
           dtr0,dtr1,ltcr,ph0,phs,dprv1,dprv2,rhog_prv,cdtrn,cdtru,cdtrt, &
           profj1,profj2,rips,ripe,nteqit,time_slc,time_snap,mdlugt,mdlni, &
           pnb_tot,pnb_eng,pnb_rw,pnb_r0, &
           ufid_bin,mdluf,mdlxp,kuf_dir,kuf_dev,kuf_dcg, &
           mdleqb,mdleqn,mdlequ,mdleqt, &
           mdlijq,mdlgmt,mdlsrc,mdlglb, &
           mdlwrt,nwrstp,kwpnam,kwtnam
    IMPLICIT NONE
    INTEGER(ikind),INTENT(IN) :: nid
    INTEGER(ikind),INTENT(OUT):: ist
    INTEGER(ikind),INTENT(OUT):: ierr

    NAMELIST /TR/ &
         RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ1,PROFJ2, &
         NSMAX,PA,PZ,PZ0,PN,PNS,PTPR,PTPP,PTS,PU,PUS, &
         PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
         RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG, &
         MODELG,MODELN,MODELQ,RHOGMN,RHOGMX, &
         KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF, &
         MODEFR,MODEFW,IDEBUG, &
         ufid_bin,mdluf,mdlxp,mdlugt,mdlni, &
         kuf_dir,kuf_dev,kuf_dcg,time_slc,time_snap, &
         mdleqb,mdleqn,mdlequ,mdleqt, &
         mdlijq,mdlgmt,mdlsrc,mdlglb, &
         nrmax,ntmax,dt,rg_fixed,nsamax,ns_nsa, &
         phia,pa_mion,pz_mion,pa_mimp,pz_mimp,      &
         lmaxtr,epsltr,mdltr_nc,mdltr_tb,mdltr_prv, &
         mdler,nteqit, &
         dtr0,dtr1,ltcr,ph0,phs,dprv1,dprv2,rhog_prv, &
         cdtrn,cdtru,cdtrt, &
         ntstep,ngtmax,ngtstp, &
         rips,ripe, &
         pnb_tot,pnb_eng,pnb_rw,pnb_r0, &
         mdlwrt,nwrstp,kwpnam,kwtnam

    READ(nid,TR,IOSTAT=ist,ERR=9800,END=9900)
    IERR=0
    RETURN

9800 IERR=8
    RETURN
9900 IERR=9
    RETURN
  END SUBROUTINE tr_nlin

! ***** parameter input through namelist *****

  SUBROUTINE tr_parm(mode,kin,ierr)

!     MODE=0 : standard namelinst input
!     MODE=1 : namelist file input
!     MODE=2 : namelist line input

!     IERR=0 : normal end
!     IERR=1 : namelist standard input error
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE
!     IERR=10X : input parameter out of range

    USE plcomm
    IMPLICIT NONE
    INTEGER(ikind),INTENT(IN) :: mode
    CHARACTER(LEN=*),INTENT(IN)::  kin
    INTEGER(ikind),INTENT(OUT):: ierr

1   CALL TASK_PARM(mode,'TR',kin,tr_nlin,tr_plist,ierr)
    IF(ierr /= 0) RETURN

    CALL tr_check_parm(ierr)

    IF(mode == 0 .AND. IERR /= 0) GOTO 1
    IF(IERR.NE.0) IERR=IERR+100

    RETURN
  END SUBROUTINE TR_PARM

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE tr_check_parm(ierr)

    USE trcomm, ONLY : ikind,nrum,nrmax,nsamax,mdltr_nc,mdltr_tb, &
         mdluf,mdlgmt,mdlsrc,mdlglb
    IMPLICIT NONE
    INTEGER(ikind), INTENT(OUT):: IERR
    CHARACTER(LEN=32) :: fmt1

    IERR=0

    ! ---------------------------------------------------------------------
    ! 
    IF(nrmax < 1) THEN
       WRITE(6,*) 'XX tr_check_parm: input error : illegal nrmax'
       WRITE(6,*) '                  nrmax =',nrmax
       IERR=1
    ENDIF

    IF(nsamax < 2) THEN
       WRITE(6,*) 'XX tr_check_parm: input error : illegal nsamax'
       WRITE(6,*) '                  nsamax =',nsamax
       IERR=1
    ENDIF

    IF(mdltr_nc==0 .AND. mdltr_tb==0) THEN
       WRITE(6,*) 'XX tr_check_parm: input error : no trasport'
       WRITE(6,*) '                mdltr_nc =',mdltr_nc
       WRITE(6,*) '                mdltr_tb =',mdltr_tb
       IERR=1
    ENDIF

    ! ---------------------------------------------------------------------
    ! check switches for experimental data
    SELECT CASE(mdluf)
    CASE(0)
       IF(mdlgmt==6 .OR. mdlgmt==7 .OR. mdlsrc==6 .OR. mdlsrc==7 .OR. &
          mdlglb==6 .OR. mdlglb==7)THEN
          WRITE(6,*) 'XX tr_check_parm: input error: experimental data are not read.'
          IERR=2
       END IF

    CASE(1)
       IF(mdlgmt==7 .OR. mdlsrc==7 .OR. mdlglb==7)THEN
          WRITE(6,*) 'XX tr_check_parm: input error: time evolution experimental data are not read.'
          IERR=2
       END IF
    END SELECT

    IF(IERR==2)THEN
       fmt1='(1X,3(A10,I2))'
       WRITE(6,fmt1) 'mdlglb= ',mdlglb,'mdlgmt= ',mdlgmt,'mdlsrc= ',mdlsrc
    END IF

    IF(mdluf > 0)THEN
       IF(nrmax+1 > nrum)THEN
          WRITE(6,*) 'XX tr_check_parm: input error: radial grid number must be less than "NRUM".'
          fmt1='(1X,2(A10,I4))'
          WRITE(6,*) 'NRMAX+1 = ',nrmax+1, 'NRUM = ',nrum
       END IF
    END IF

    ! ---------------------------------------------------------------------

    RETURN
  END SUBROUTINE tr_check_parm

! ***** show input parameters *****

  SUBROUTINE tr_view

    USE plcomm
    USE plinit
    USE trcomm, ONLY: &
           nrmax,ntmax,dt,rg_fixed,nsamax,ns_nsa, &
           lmaxtr,epsltr,nitmax,mdltr_nc,mdltr_tb,mdltr_prv, &
           dtr0,dtr1,ltcr,ph0,phs,dprv1,dprv2,cdtrn,cdtru,cdtrt, &
           ntstep,ngtmax,ngtstp
    IMPLICIT NONE
    INTEGER(ikind):: nsa

    CALL pl_view

    WRITE(6,*) ! --------------------------------------------------
    WRITE(6,*) '** TRANSPORT **'
    WRITE(6,603) 'nrmax  ',nrmax ,'ntmax  ',ntmax , &
                 'nsamax ',nsamax,'lmaxtr ',lmaxtr
    WRITE(6,603) 'ntstep ',ntstep,'ngtstp ',ngtstp ,&
                 'ngtmax ',ngtmax,'nitmax ',nitmax
    WRITE(6,'(1X,A,15I4)') 'NSA:',(NSA,NSA=1,NSAMAX)
    WRITE(6,'(1X,A,15I4)') 'NS :',(NS_NSA(NSA),NSA=1,NSMAX)

    WRITE(6,*) ! --------------------------------------------------
    WRITE(6,601) 'dt        ',dt,       'epsltr    ',epsltr
    DO nsa=1,nsamax
       IF(MIN(rg_fixed(1,nsa),rg_fixed(2,nsa),rg_fixed(3,nsa)) < rb) THEN
          WRITE(6,'(A,I5,1P3E12.4)') 'rg_fixed(n,u,t): nsa=', &
               nsa,rg_fixed(1,nsa),rg_fixed(2,nsa),rg_fixed(3,nsa)
       END IF
    END DO
    WRITE(6,602) 'mdltr_nc  ',mdltr_nc, 'mdltr_tb  ',mdltr_tb, &
                 'mdltr_prv ',mdltr_prv
    WRITE(6,601) 'dtr0      ',dtr0,     'dtr1      ',dtr1, &
                 'ltcr      ',ltcr
    WRITE(6,601) 'drpv1     ',dprv1,    'dprv2     ',dprv2
    WRITE(6,601) 'cdtrn     ',cdtrn,    'cdtru     ',cdtru     , &
                 'cdtrt     ',cdtrt

    WRITE(6,*) ! --------------------------------------------------
    WRITE(6,601) 'ph0       ',ph0,      'phs       ',phs
    RETURN

  601 FORMAT(' ',A10,'=',1PE12.4 :2X,A10,'=',1PE12.4: &
     &        2X,A10,'=',1PE12.4)
  602 FORMAT(' ',A10,'=',I12 :2X,A10,'=',I12: &
     &        2X,A10,'=',I12)
  603 FORMAT(' ',A7, '=',I7,3X   :2X,A7, '=',I7,3X  : &
     &        2X,A7, '=',I7,3X   :2X,A7, '=',I7)

  END SUBROUTINE tr_view

END MODULE trinit
