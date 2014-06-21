 &t2	
 ntmax = 1
 ntstep= 1	 
 dt    = 1.D-5
 idebug= 0
 NLMAX = 1
 NPMIN = 50
! NSMAX = 0
 i1mlvl(1) = 1
 i1rdn2(1) = 50
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
    UseSUPGFEM              = .FALSE.
    UseCoefficientCheck     = .FALSE.
    UseAnomalousTransportFT = .TRUE.
    UseAnomalousTransportGT = .TRUE.
    
    ! set equations to be solved
    SolveElectron = .FALSE.
    SolveIons     = .FALSE.

    SolveBp       = .TRUE.
    SolveBt       = .TRUE.
    SolveEt       = .TRUE.
    SolveEp       = .TRUE.
    SolveEr       = .FALSE.
    SolveNn       = .FALSE.
    SolveFr       = .TRUE.
    SolveFb       = .FALSE.
    SolveFt       = .FALSE.
    SolveFp       = .FALSE.
    SolvePp       = .FALSE.
    SolveQr       = .FALSE.
    SolveQb       = .FALSE.
    SolveQt       = .FALSE.
    SolveQp       = .FALSE.
    
    ! set dirichlet boundary condition on magnetic axis
    LockBpOnAxis  = .TRUE.
    LockBtOnAxis  = .TRUE.
    LockEtOnAxis  = .TRUE.
    LockEpOnAxis  = .TRUE.
    LockErOnAxis  = .TRUE.
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
    LockQpOnWall  = .TRUE.

 &end
