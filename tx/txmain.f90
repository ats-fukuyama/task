!     $Id$
!
!***********************************************************
!
!             TOKAMAK TRANSPORT SIMULATION CODE
!    INCLUDING RADIAL ELECTRIC FIELD AND PLASMA ROTATION
!
!      FEM scheme with Linear interpolation function
!
!     DEVELOPED BY M. HONDA, Y. FUJI AND A. FUKUYAMA
!
!             DEPARTMENT OF NUCLEAR ENGINEERING
!              GRADUATE SCHOOL OF ENGINEERING
!                     KYOTO UNIVERSITY
!                      KYOTO 606-8501
!
!                JAPAN ATOMIC ENERGY AGENCY
!                   NAKA FUSION INSTITUTE
!                  TOKAMAK ANALYSIS GROUP
!
!     Reference:
!       M. Honda and A. Fukuyama, "Dynamic transport simulation
!       code including plasma rotation and radial electric field",
!       Journal of Computational Physics, 227 (2008) 2808-2844
!
!***********************************************************
!
!   Variables                                         Boundary Cond.
!                                                      axis   edge
!   X1  = phi                                           N      0
!   X2  = r * Atheta'                                   0      D 
!   X3  = Aphi'                                         N      D
!   X4  = Aphi                                          *      *
!   X5  = r * Atheta                                    *      *
!   X6  = Ne                 X11  = Ni                  N(*)   N(*)
!   X7  = r * Ne * Uer       X12  = r * Ni * Uir        0      *(0)
!   X8  = r * Ne * UeTheta   X13  = r * Ni * UiTheta    0      0
!   X9  = Ne * UePhi         X14  = Ni * UiPhi          N      0
!   X10 = Ne * Te            X15  = Ni * Ti             N      N
!   X16 = Nb                                            *(N)   *(N)
!   X17 = r * Nb * UbTheta                              *(0)   *(0)
!   X18 = Nb * UbPhi                                    *(N)   *(0)
!   X19 = Slow neutral                                  N      D
!   X20 = Fast neutral                                  N      0
!   X21 = Nb ripple                                     N      N
!
!   In case of MDFIXT /=0, X10 = Te, X15 = Ti
!
!   G16 = Total n0                        
!   G20 = Q                               
!   G21 = Jphi                            
!   G22 = JePhi     G23 = JiPhi           
!   G24 = JbPhi                           
!   G25 = UePerp    G26 = UePara          
!   G27 = UiPerp    G28 = UiPara          
!   G29 = D    G30 = G1h2    G31 = MuI    
!   G32 = S    G33 = Alpha   G34 = rKappa 
!   G35 = Slow N0   G36 = Fast N0
  

PROGRAM TASK_TX
  
  use tx_commons, only : SLID
  use tx_parameter_control, only : TXPARF, TXPARM_CHECK
  implicit none
  character(len=80) :: KPNAME

  CALL GSOPEN
  OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')

  !     ***** Version ID *****
  !     SLID is used to identify data file.
  SLID = 'tx454.0'
  WRITE(6,*) '######## TASK/TX V4.54.10 09/01/22 ########'

  CALL TXINIT
  KPNAME='txparm'
  CALL TXPARF(KPNAME)
  CALL TXPARM_CHECK

  CALL TXMENU

  CLOSE(7)
  CALL GSCLOS
  STOP
END PROGRAM TASK_TX
