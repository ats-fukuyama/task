 &t2
 ntmax = 100
 ntstep= 100
 dt    = 1.D-3
 idebug= 0
 NLMAX = 1
 NPMIN = 300
 NSMAX = 0
 i1mlvl(1) = 1
 i1rdn2(1) = 299
 eps_conv  = 1.D-1
 
 CoordinateSwitch = 2
 TestCase         = 3

    TestMS  = .TURE.
    TestAV  = .TRUE.
    TestAT  = .FALSE.
    TestDT  = .FALSE.
    TestGV  = .FALSE.
    TestGT  = .FALSE.
    TestES  = .FALSE.
    TestEV  = .FALSE.
    TestET  = .FALSE.
    TestSS  = .FALSE.
    TestLEQ = .FALSE.
    TestLAX = .TRUE.
    TestLWL = .TRUE.
    
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

    SolveBp       = .TRUE.
    SolveBt       = .FALSE.
    SolveEt       = .FALSE.
    SolveEp       = .FALSE.
    SolveEr       = .FALSE.
    SolveNn       = .FALSE.
    SolveFr       = .FALSE.
    SolveFb       = .FALSE.
    SolveFt       = .FALSE.
    SolveFp       = .FALSE.
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

    LockBpOnWall  = .TRUE.
    LockBpOnWall  = .FALSE.
    LockEtOnWall  = .TRUE.
    LockEpOnWall  = .FALSE.
    LockErOnWall  = .FALSE.
    LockNnOnWall  = .FALSE.
    LockFrOnWall  = .FALSE.
    LockFbOnWall  = .FALSE.
    LockFtOnWall  = .FALSE.
    LockFpOnWall  = .FALSE.
    LockPpOnWall  = .FALSE.
    LockQrOnWall  = .FALSE.
    LockQbOnWall  = .FALSE.
    LockQtOnWall  = .FALSE.
    LockQpOnWall  = .FALSE.

 &end
