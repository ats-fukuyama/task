 &t2
 ntmax = 100
 ntstep= 100
 dt    = 1.D-7
 idebug= 0
 NLMAX = 1
 NPMIN = 20
! NSMAX = 0
 i1mlvl(1) = 1
 i1rdn2(1) = 10
 eps_conv  = 1.D-2
 
 CoordinateSwitch = 1
 TestCase         = 3
 EqSet		  = 2

    TestMS  = .FALSE.
    TestAV  = .FALSE.
    TestAT  = .FALSE.
    TestDT  = .FALSE.
    TestGV  = .FALSE.
    TestGT  = .FALSE.
    TestES  = .FALSE.
    TestEV  = .FALSE.
    TestET  = .FALSE.
    TestSS  = .FALSE.
    TestLEQ = .FALSE.
    TestLAX = .FALSE.
    TestLWL = .FALSE.
    
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
    SolveBt       = .TRUE.
    SolveEt       = .TRUE.
    SolveEp       = .FALSE.
    SolveEr       = .FALSE.
    SolveNn       = .TRUE.
    SolveFr       = .TRUE.
    SolveFb       = .FALSE.
    SolveFt       = .TRUE.
    SolveFp       = .TRUE.
    SolvePp       = .FALSE.
    SolveQr       = .FALSE.
    SolveQb       = .FALSE.
    SolveQt       = .FALSE.
    SolveQp       = .FALSE.
    
    ! set dirichlet boundary condition on magnetic axis
    LockBpOnAxis  = .TRUE.
    LockBtOnAxis  = .TRUE.
    LockEtOnAxis  = .TRUE.
    LockEpOnAxis  = .FALSE.
    LockErOnAxis  = .FALSE.
    LockNnOnAxis  = .FALSE.
    LockFrOnAxis  = .TRUE.
    LockFbOnAxis  = .FALSE.
    LockFtOnAxis  = .FALSE.
    LockFpOnAxis  = .FALSE.
    LockPpOnAxis  = .FALSE.
    LockQrOnAxis  = .TRUE.
    LockQbOnAxis  = .FALSE.
    LockQtOnAxis  = .FALSE.
    LockQpOnAxis  = .FALSE.
    
    ! set dirichlet boundary condition on first wall

    LockBpOnWall  = .FALSE.
    LockBtOnWall  = .FALSE.
    LockEtOnWall  = .FALSE.
    LockEpOnWall  = .FALSE.
    LockErOnWall  = .TRUE.
    LockNnOnWall  = .FALSE.
    LockFrOnWall  = .TRUE.
    LockFbOnWall  = .FALSE.
    LockFtOnWall  = .TRUE.
    LockFpOnWall  = .FALSE.
    LockPpOnWall  = .FALSE.
    LockQrOnWall  = .TRUE.
    LockQbOnWall  = .FALSE.
    LockQtOnWall  = .TRUE.
    LockQpOnWall  = .TRUE.

 &end
