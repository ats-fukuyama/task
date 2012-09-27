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
           lmaxtr,epsltr,mdltr_nc,mdltr_tb,mdltr_prv, &
           mdluf,mdlxp,mdler,modelg,mdleqn,nteqit,time_slc, &
           dtr0,dtr1,ltcr,ph0,phs,dprv1,dprv2,cdtrn,cdtru,cdtrt, &
           ntstep,ngtmax,ngtstp,rips,ripe,profj1,profj2, &
           pnb_tot,pnb_eng,pnb_rw,pnb_r0, &
           kuf_dir,kuf_dev,kuf_dcg
    USE plinit
    IMPLICIT NONE
    INTEGER(ikind):: nsa

    OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
    CALL pl_init

!     ======( DEVICE PARAMETERS )======

!        RR    : Plasma major radius                             (m)
!        RA    : Plasma minor radius                             (m)
!        RB    : Wall minor radius                               (m)
!        RKAP  : Plasma shape elongation
!        RDEL  : Plasma shape triangularity *
!        BB    : Magnetic field at center                        (T)
!        Q0    : Safety factor at center
!        QA    : Safety factor on plasma surface
!        RIP   : Plasma current                                 (MA)
!        PROFJ : Curren density profile parameter (power of (1 - rho^2))  

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


!     ======( PROFILE PARAMETERS )======

!        PROFN1: Density profile parameter (power of rho)
!        PROFN2: Density profile parameter (power of (1 - rho^PROFN1))
!        PROFT1: Temperature profile parameter (power of rho)
!        PROFT2: Temperature profile parameter (power of (1 - rho^PROFN1))
!        PROFU1: Rotation profile parameter (power of rho)
!        PROFU2: Rotation profile parameter (power of (1 - rho^PROFN1))


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

    modelg = 0
    nteqit = 0


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

!     ==== Parameters above are initialized in plx/plinit ====

    CALL pl_init

!     ==== TR PARAMETERS ====

!        NRMAX  : NUMBER OF RADIAL MESH POINTS
!        NTMAX  : NUMBER OF TIME STEP
!        DT     : SIZE OF TIME STEP

!        rg_fixed(3,nsm) : minimum radius of fixed profile
!        nsamax   : number of active particle species
!        nitmax   : maximum number of iterations
!        epsit    : tolerance of iteration

    nrmax    = 50
    ntmax    = 100
    dt       = 0.01D0
    rg_fixed(1:3,0:nsm) = rb
    nsamax   = 2
    DO nsa = 1, nsamax
       ns_nsa(nsa) = nsa
    END DO

!     ==== Convergence Parameter ====
!        epsltr : convergence criterion of iteration
!        lmaxtr : maximum count of iteration

    epsltr = 1.D-6
!    epsltr = 1.D99
    lmaxtr = 100
        
!     ==== TR PARAMETERS for stiff modeling by Ikari  ====

!       mdltr_nc  = 0 : no neoclassical transport
!                   1 : NCLASS trasport model
!       mdltr_tb  = 0 : no turbulent transport
!                   1 : constant diffusion
!                   2 : STIFF MODEL (Pereverzev)
!       mdltr_prv = 0 : no Pereverzev method
!                   1 : Pereverzev method applied : Denh=dprv1
!                   2 : Pereverzev method applied : Denh=dprv2*Dorg
!                   3 : Pereverzev method applied : Denh=dprv2*Dorg+dprv1
!                   4 : Pereverzev method applied : Denh=MIN(dprv2*Dorg,dprv1)
!                   5 : Pereverzev method applied : Denh=MAX(dprv2*Dorg,dprv1)
!       dtr0        : lower diffusion coefficient
!       dtr1        : upper diffusion coefficnet
!       ltcr      : critical scale length [m]
!       ph0       : heating power density [MW/m^3] at r = 0
!       phs       : heating power density [MW/m^3] at r = a
!       dprv1     : enhanced diffusion coefficient
!       dprv2     : diffusion enhancement factor
!       cdtrn     : factor for particle diffusivity
!       cdtru     : factor for toroidal viscosity
!       cdtrt     : factor for thermal diffusivity

    mdltr_nc  = 1
    mdltr_tb  = 1
    mdltr_prv = 0
    dtr0    = 0.1D0
    dtr1    = 0.5D0
    ltcr  = 1.D0
    ph0   = 0.1D0
    phs   = 0.1D0
    dprv1 = 1.D0
    dprv2 = 1.D0
    cdtrn = 1.D0
    cdtru = 1.D0
    cdtrt = 1.D0

!!$    !     ==== Eqs. Selection Parameter ====
!!$
!!$    MDLEQB=1  ! 0/1 for B_theta
    MDLEQN=0  ! 0/1 for density
!!$    MDLEQT=1  ! 0/1 for heat
!!$    MDLEQU=0  ! 0/1 for rotation
!!$    MDLEQZ=0  ! 0/1 for impurity
!!$    MDLEQ0=0  ! 0/1 for neutral
!!$    MDLEQE=0  ! 0/1/2 for electron density
!!$    !               ! 0: electron only, 1: both, 2: ion only
!!$    MDLEOI=0  ! 0/1/2 for electron only or bulk ion only if NSMAX=1
!!$    !               ! 0: both, 1: electron, 2: ion

!     ==== TR PARAMETERS for data saving  ====
!        ntstep    : number of time step for status report
!        ngtmax    : maximum number of saved data
!        ngtstp    : number of time step for data save

    ntstep =    10
    ngtmax = 10001
    ngtstp =     1

!     ==== Input from experimental data ==== 
!        MDLXP :
!            0 : from UFILEs
!         else : MDSplus

!        MDLUF :
!            0 : not used
!            1 : steady state
!            2 : time evolution
!            3 : compared with TOPICS

      mdlxp = 0
      mdluf = 0

      time_slc = -1.d0 ! not determined

!     ==== DEVICE NAME AND SHOT NUMBER IN UFILE DATA ====
!        KUF_DIR : UFILE directory
!        KUF_DEV : Device name
!        KUF_DCG : Discharge number

      kuf_dir = '../../../profiledb/profile_data'
      kuf_dev = 'd3d'
      kuf_dcg = '103818'

!     ==== IMPURITY TREATMENT ====
!        MDNI  :
!            0 : NSMAX=2, ne=ni
!            1 : Use exp. ZEFFR profile if available
!            2 : Use exp. NM1 (or NM2) profile if available
!            3 : Use exp. NIMP profile if available
!      MDNI=0



! ==== for the time being ===
!        PROFJ1: Current density profile parameter (power of rho)
!        PROFJ2: Current density profile parameter (power of (1 - rho^PROFJ1))

!        MDLER: the type of radial electric field

    profj1 = 1.5d0
    profj2 = 1.5d0

    rips = 1.d0
    ripe = 1.d0

    mdler  = 2

!   ==== NBI heating ===
    pnb_tot = 1.d0  ! [MW]
    pnb_r0  = 0.3d0
    pnb_rw  = 0.5d0
    pnb_eng = 80.d0 ! [keV]          

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
           lmaxtr,epsltr,mdltr_nc,mdltr_tb,mdltr_prv,mdler,&
           dtr0,dtr1,ltcr,ph0,phs,dprv1,dprv2,cdtrn,cdtru,cdtrt, &
           profj1,profj2,rips,ripe,nteqit,time_slc, &
           pnb_tot,pnb_eng,pnb_rw,pnb_r0, &
           mdluf,mdlxp,kuf_dir,kuf_dev,kuf_dcg
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
         mdluf,mdlxp, &
         kuf_dir,kuf_dev,kuf_dcg,time_slc, &
         nrmax,ntmax,dt,rg_fixed,nsamax,ns_nsa, &
         lmaxtr,epsltr,mdltr_nc,mdltr_tb,mdltr_prv, &
         mdler,nteqit, &
         dtr0,dtr1,ltcr,ph0,phs,dprv1,dprv2,cdtrn,cdtru,cdtrt, &
         ntstep,ngtmax,ngtstp, &
         rips,ripe, &
         pnb_tot,pnb_eng,pnb_rw,pnb_r0

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

    CALL tr_input_check(ierr)

    IF(mode == 0.AND. IERR /= 0) GOTO 1
    IF(IERR.NE.0) IERR=IERR+100

    RETURN
  END SUBROUTINE TR_PARM

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE tr_input_check(ierr)

    USE trcomm, ONLY : ikind,nrmax,nsamax,mdltr_nc,mdltr_tb
    IMPLICIT NONE
    INTEGER(ikind), INTENT(OUT):: IERR

    IERR=0

    IF(nrmax < 1) THEN
       WRITE(6,*) 'XXX tr_check: input error : illegal nrmax'
       WRITE(6,*) '                  nrmax =',nrmax
       IERR=1
    ENDIF

    IF(nsamax < 2) THEN ! nsamax > nsmax ??
       WRITE(6,*) 'XXX tr_check: input error : illegal nsamax'
       WRITE(6,*) '                  nsamax =',nsamax
       IERR=1
    ENDIF

    IF(mdltr_nc==0 .AND. mdltr_tb==0) THEN
       WRITE(6,*) 'XXX tr_check: input error : no trasport'
       WRITE(6,*) '                mdltr_nc =',mdltr_nc
       WRITE(6,*) '                mdltr_tb =',mdltr_tb
       IERR=1
    ENDIF

    RETURN
  END SUBROUTINE tr_input_check

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
