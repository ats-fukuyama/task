 &t2
 ntmax=1
 ntstep=1
 dt=1.D-10
 idebug=0
 NLMAX = 1
 NPMIN = 10
 i1mlvl(1) = 1
 i1rdn2(1) = 22
 eps_conv = 1.D-4

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
    SolveNn       = .FALSE.
    SolveFr       = .FALSE.
    SolveFb       = .TRUE.
    SolveFt       = .TRUE.
    SolveFp       = .TRUE.
    SolvePp       = .FALSE.
    SolveQr       = .FALSE.
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

    LockBpOnWall  = .FALSE.
    LockBpOnWall  = .FALSE.
    LockEtOnWall  = .FALSE.
    LockEpOnWall  = .FALSE.
    LockErOnWall  = .FALSE.
    LockNnOnWall  = .FALSE.
    LockFrOnWall  = .FALSE.
    LockFbOnWall  = .FALSE.
    LockFtOnWall  = .TRUE.
    LockFpOnWall  = .TRUE.
    LockPpOnWall  = .FALSE.
    LockQrOnWall  = .FALSE.
    LockQbOnWall  = .FALSE.
    LockQtOnWall  = .TRUE.
    LockQpOnWall  = .TRUE.

 &end
