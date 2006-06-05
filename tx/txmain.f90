!     $Id$
!
!***********************************************************
!
!               TOKAMAK TRANSPORT SIMULATION
!    INCLUDING RADIAL ELECTRIC FIELD AND PLASMA ROTATION
!
!     DEVELOPED BY Y. FUJI, M. HONDA AND A. FUKUYAMA
!
!             DEPARTMENT OF NUCLEAR ENGINEERING
!              GRADUATE SCHOOL OF ENGINEERING
!                      KYOTO UNIVERSITY
!                      KYOTO 606-8501
!
!***********************************************************
!
!   Varibles                                  Boundary Cond.
!                                              axis   edge
!   X1  = Er                                            0      N
!   X2  = r * Etheta                                    0      N 
!   X3  = Ephi                                          N      N
!   X4  = 1 / r d/dr(r * Btheta)                        N      N
!   X5  = Bphi                                          N      BB
!   X6  = Ne                 X11  = Ni                  N      N
!   X7  = Ne * Uer           X12  = Ni * Uir            0      N
!   X8  = Ne * UeTheta / r   X13  = Ni * UiTheta / r    N      N
!   X9  = Ne * UePhi         X14  = Ni * UiPhi          N      N
!   X10 = Ne * Te            X15  = Ni * Ti             N      N
!   X16 = Nb                                            N      N
!   X17 = Nb * UbTheta / r                              N      N
!   X18 = Nb * UbPhi                                    N      N
!   X19 = Slow n01                                      N      N
!   X20 = Fast n02                                      N      N
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
  
  use menu, only : TXMENU
  use init_prof, only : TXINIT
  use commons, only : SLID
  use parameter_control, only : TXPARF
  implicit none
  character(len=80) :: KPNAME

  CALL GSOPEN
  OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')

  !     ***** Version ID *****
  !     SLID is used to identify data file.
  SLID = 'tx300.0'
  WRITE(6,*) '######## TASK/TX V3.00.00 06/06/05 ########'

  CALL TXINIT
  KPNAME='txparm'
  CALL TXPARF(KPNAME)

  CALL TXMENU

  CLOSE(7)
  CALL GSCLOS
  STOP
END PROGRAM TASK_TX
