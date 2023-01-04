!---------------------------------------------------------------------
!
!         MODULE FOR SETUP DEFAULT INPUT PARAMETER 
!
!                   LAST UPDATE 2014-06-05 H.Seto
!
!---------------------------------------------------------------------
MODULE T2INIT
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC T2_INIT !C INITIALIZE INPUT PARAMETERS
  
CONTAINS
  
  !
  !
  !
  SUBROUTINE T2_INIT
    
    USE T2CNST, ONLY: i0lmaxm,i0spcsm,Ame,Amp
    
    USE T2COMM,ONLY:&
         c10rname,i0fnum,&
         NNMAX,NQMAX,NDMAX,NSMAX,NPMIN,NLMAX,&
         !
         i1mlvl,i1rdn2,d1rec,d0rw,RR,RA,&
         i0nm,i0nn,i0tm,i0tn,d0qc,d0qs,d0bc,&
         d1nc,d1ns,d1nw,d1tc,d1ts,d1tw,Pa,Pz,&
         ! AF
         dt,time_init,eps_conv, &
         ntmax,ntstep,nt0dmax,nt0dstep,nt2dmax,nt2dstep,nconvmax, &
         idfile,idprint,idplot,idmode,idebug,&
         !
         UsePotentialDescription,UseNormalization,&
         UseSUPGFEM,             UseCoefficientCheck,&
         UseAnomalousTransportFT,UseAnomalousTransportGT,&
         !
         SolveElectron,SolveIons,&
         SolveBp,SolveBt,SolveEt,SolveEp,SolveEr,&
         SolveNn,SolveFr,SolveFb,SolveFt,SolveFp,&
         SolvePp,SolveQr,SolveQb,SolveQt,SolveQp,&
         !
         LockBpOnAxis,LockBtOnAxis,LockEtOnAxis,LockEpOnAxis,&
         LockErOnAxis,&
         LockNnOnAxis,LockFrOnAxis,LockFbOnAxis,LockFtOnAxis,&
         LockFpOnAxis,&
         LockPpOnAxis,LockQrOnAxis,LockQbOnAxis,LockQtOnAxis,&
         LockQpOnAxis,&
         !
         LockBpOnWall,LockBpOnWall,LockEtOnWall,LockEpOnWall,&
         LockErOnWall,&
         LockNnOnWall,LockFrOnWall,LockFbOnWall,LockFtOnWall,&
         LockFpOnWall,&
         LockPpOnWall,LockQrOnWall,LockQbOnWall,LockQtOnWall,&
         LockQpOnWall,&
         !
         CoordinateSwitch
    
    ! ----------------------------------------------------------------

    c10rname = 'TEST'
    i0fnum  = 10
    !
    CoordinateSwitch = 1

    ! settings for computational domain
    NNMAX =  4
    NQMAX = 32
    NDMAX =  2
    NSMAX =  2
    NPMIN = 10
    NLMAX =  1
    
    RR   =  1.30D0
    RA   =  0.30D0
    d0rw =  1.00D0
    d1rec(0:i0lmaxm) = 0.D0 ! least radial point in a level
    d1rec(1) = d0rw
    

    i1mlvl(0:i0lmaxm+1) = 0
    i1mlvl(1)=1        ! 10 x 2^0
    
    i1rdn2(-1:i0lmaxm) = 0  
    i1rdn2(1) = 11   ! number of radial nodes in a level

    ! plasma parameters
    
    d0qc     =  1.D0
    d0qs     =  3.D0
    d0bc     =  1.00D0

    i0nm = 2
    i0nn = 2
    i0tm = 2
    i0tn = 2
    
    !Electron 
    Pa(1) =  Ame/Amp
    Pz(1) = -1.D0
    
    d1nc(1) = 0.30D0
    d1ns(1) = 0.03D0
    d1nw(1) = 0.03D0
    d1tc(1) = 0.50D0
    d1ts(1) = 0.05D0
    d1tw(1) = 0.05D0
    !Deuterium
    Pa(2) = 2.D0
    Pz(2) = 1.D0
    
    d1nc(2) = 0.30D0
    d1ns(2) = 0.03D0
    d1nw(2) = 0.03D0
    d1tc(2) = 0.50D0
    d1ts(2) = 0.05D0
    d1tw(2) = 0.05D0
    
    !Other plasma speces
    d1nc(3:i0spcsm) = 0.D0
    d1ns(3:i0spcsm) = 0.D0 
    d1tc(3:i0spcsm) = 0.D0
    d1ts(3:i0spcsm) = 0.D0
    Pa(3:i0spcsm)   = 0.D0
    Pz(3:i0spcsm)   = 0.D0
    
    dt        = 1.D-4   ! time step [s]
    time_init = 0.D0    ! initial time [s]
    ntmax     = 1       ! number of time steps to go
    ntstep    = 1       ! time step to print snap shot of global data
    nt0dmax   = 1       ! maximim number of global data to be saved
    nt0dstep  = 1       ! time step to save global data
    nt2dmax   = 1       ! maximum number of profile data to be saved
    nt2dstep  = 1       ! time step to save profile data

    nconvmax  = 255     ! maximum number of convergence steps for implicit loop
    eps_conv  = 1.D-4   ! relative convergence criterion for implicit loop

    idfile    = 0       ! control id for file output: 0 for none, 9 for all
    idprint   = 9       ! control id for print output: 0 for none, 9 for all
    idplot    = 1       ! control id for plot type: 1 for contour
    idmode    = 1       ! control id for equation level: 
                        !     1 for electron density   
                        !     2 for electron and ion
                        !     3 for electron energy
                        !     4 for electron and ion
                        !     5 for Ampere, Faraday
                        !     6 for Gauss's law
                        !     7 with flux surface average
                        !     8 with equilibrium
    idebug    = 0       ! control id for debug mode
    
    ! 
    UsePotentialDescription = .FALSE.
    UseNormalization        = .TRUE.
    UseSUPGFEM              = .FALSE.
    UseCoefficientCheck     = .FALSE.
    UseAnomalousTransportFT = .FALSE.
    UseAnomalousTransportGT = .FALSE.
    
    ! set equations to be solved
    SolveElectron = .TRUE.
    SolveIons     = .TRUE.
    SolveBp       = .TRUE.
    SolveBt       = .TRUE.
    SolveEt       = .TRUE.
    SolveEp       = .TRUE.
    SolveEr       = .TRUE.
    SolveNn       = .TRUE.
    SolveFr       = .TRUE.
    SolveFb       = .TRUE.
    SolveFt       = .TRUE.
    SolveFp       = .TRUE.
    SolvePp       = .TRUE.
    SolveQr       = .TRUE.
    SolveQb       = .TRUE.
    SolveQt       = .TRUE.
    SolveQp       = .TRUE.
    
    ! set dirichlet boundary condition on magnetic axis
    LockBpOnAxis  = .FALSE.
    LockBtOnAxis  = .FALSE.
    LockEtOnAxis  = .FALSE.
    LockEpOnAxis  = .FALSE.
    LockErOnAxis  = .FALSE.
    LockNnOnAxis  = .FALSE.
    LockFrOnAxis  = .FALSE.
    LockFbOnAxis  = .FALSE.
    LockFtOnAxis  = .FALSE.
    LockFpOnAxis  = .FALSE.
    LockPpOnAxis  = .FALSE.
    LockQrOnAxis  = .FALSE.
    LockQbOnAxis  = .FALSE.
    LockQtOnAxis  = .FALSE.
    LockQpOnAxis  = .FALSE.
    
    ! set dirichlet boundary condition on first wall

    LockBpOnWall  = .TRUE.
    LockBpOnWall  = .TRUE.
    LockEtOnWall  = .TRUE.
    LockEpOnWall  = .TRUE.
    LockErOnWall  = .TRUE.
    LockNnOnWall  = .FALSE.
    LockFrOnWall  = .TRUE.
    LockFbOnWall  = .TRUE.
    LockFtOnWall  = .TRUE.
    LockFpOnWall  = .FALSE.
    LockPpOnWall  = .FALSE.
    LockQrOnWall  = .TRUE.
    LockQbOnWall  = .TRUE.
    LockQtOnWall  = .TRUE.
    LockQpOnWall  = .FALSE.

    RETURN
    
  END SUBROUTINE T2_INIT
END MODULE T2INIT
    
