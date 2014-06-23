 &t2	
 ntmax = 1
 ntstep= 1	 
 dt    = 1.D-5
 idebug= 0
 NLMAX = 1
 NPMIN = 10
! NSMAX = 0
 i1mlvl(1) = 1
 i1rdn2(1) = 10
 eps_conv  = 1.D-3
 
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
    !UseNormalization        = .FALSE.
    UseSUPGFEM              = .FALSE.
    UseCoefficientCheck     = .FALSE.
    UseAnomalousTransportFT = .TRUE.
    UseAnomalousTransportGT = .TRUE.
    
    ! set equations to be solved
    SolveElectron = .TRUE.
    SolveIons     = .FALSE.

    SolveBp       = .TRUE.
    SolveBt       = .TRUE.
    SolveEt       = .TRUE.
    SolveEp       = .FALSE.
    SolveEr       = .FALSE.
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
    LockBpOnAxis  = .TRUE.
    LockBtOnAxis  = .TRUE.
    LockEtOnAxis  = .TRUE.
    LockEpOnAxis  = .TRUE.
    LockErOnAxis  = .TRUE.
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
    LockErOnWall  = .FALSE.
    LockNnOnWall  = .FALSE.
    LockFrOnWall  = .FALSE.
    LockFbOnWall  = .FALSE.
    LockFtOnWall  = .TRUE.
    LockFpOnWall  = .FALSE.
    LockPpOnWall  = .FALSE.
    LockQrOnWall  = .FALSE.
    LockQbOnWall  = .FALSE.
    LockQtOnWall  = .TRUE.
    LockQpOnWall  = .FALSE.

 &end
