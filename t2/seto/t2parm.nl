 &t2
 ntmax=1
 ntstep=1
 dt=1.D-10
 idebug=0
 NLMAX = 1
 NPMIN = 60
 i1mlvl(1) = 1
 i1rdn2(1) = 66
 eps_conv = 1.D-2

    ! 
    UsePotentialDescription = .FALSE.
    UseNormalization        = .TRUE.
    UseSUPGFEM              = .FALSE.
    UseCoefficientCheck     = .FALSE.
    UseAnomalousTransportFT = .TRUE.
    UseAnomalousTransportGT = .TRUE.
    
    ! set equations to be solved
    SolveElectron = .TRUE.
    SolveIons     = .TRUE.

    SolveBp       = .FALSE.
    SolveBt       = .FALSE.
    SolveEt       = .FALSE.
    SolveEp       = .FALSE.
    SolveEr       = .FALSE.
    SolveNn       = .TRUE.
    SolveFr       = .TRUE.
    SolveFb       = .TRUE.
    SolveFt       = .TRUE.
    SolveFp       = .TRUE.
    SolvePp       = .FALSE.
    SolveQr       = .FALSE.
    SolveQb       = .FALSE.
    SolveQt       = .FALSE.
    SolveQp       = .FALSE.
    
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

    LockBpOnWall  = .FALSE.
    LockBpOnWall  = .FALSE.
    LockEtOnWall  = .FALSE.
    LockEpOnWall  = .FALSE.
    LockErOnWall  = .FALSE.
    LockNnOnWall  = .FALSE.
    LockFrOnWall  = .TRUE.
    LockFbOnWall  = .TRUE.
    LockFtOnWall  = .TRUE.
    LockFpOnWall  = .FALSE.
    LockPpOnWall  = .FALSE.
    LockQrOnWall  = .FALSE.
    LockQbOnWall  = .FALSE.
    LockQtOnWall  = .FALSE.
    LockQpOnWall  = .FALSE.

 &end
