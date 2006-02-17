!     $Id$
!
!***********************************************************
!
!               TOKAMAK TRANSPORT SIMULATION
!    INCLUDING RADIAL ELECTRIC FIELD AND PLASMA ROTATION
!
!     DEVELOPED BY Y. FUJI AND A. FUKUYAMA
!
!             DEPARTMENT OF NUCLEAR ENGINEERING
!              GRADUATE SCHOOL OF ENGINEERING
!                      KYOTO UNIVERSITY
!                      KYOTO 606-8501
!
!***********************************************************
!
!   X1  = Er                                       integer
!   X2  = Etheta                                   integer
!   X3  = Ephi                                half integer
!   X4  = Btheta                                   integer
!   X5  = Bphi                                half integer
!   X6  = Ne             X11  = Ni            half integer
!   X7  = Ne * Uer       X12  = Ni * Uir           integer
!   X8  = Ne * UeTheta   X13  = Ni * UiTheta       integer
!   X9  = Ne * UePhi     X14  = Ni * UiPhi    half integer
!   X10 = Ne * Te        X15  = Ni * Ti       half integer
!   X16 = Nb                                  half integer
!   X17 = Nb * UbTheta                             integer
!   X18 = Nb * UbPhi                          half integer
!   X19 = Slow n01                            half integer
!   X20 = Fast n02                            half integer
!
!   G16 = Total n0                            half integer
!   G20 = Q                                        integer
!   G21 = Jphi                                     integer
!   G22 = JePhi     G23 = JiPhi                    integer
!   G24 = JbPhi                                    integer
!   G25 = UePerp    G26 = UePara                   integer
!   G27 = UiPerp    G28 = UiPara                   integer
!   G29 = D    G30 = G1h2    G31 = MuI        half integer
!   G32 = S    G33 = Alpha   G34 = rKappa     half integer
!   G35 = Slow N0   G36 = Fast N0             half integer

PROGRAM TASK_TX
  
  use menu
  use init_prof

  INCLUDE 'txcomm.inc'
  CHARACTER(80) :: KPNAME

  CALL GSOPEN
  OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')

  !     ***** Version ID *****
  !     SLID is used to identify data file.
  SLID = 'tx300.0'
  WRITE(6,*) '######## TASK/TX V3.00.00 04/11/03 ########'

  CALL TXINIT
  KPNAME='txparm'
  CALL TXPARF(KPNAME)

  CALL TXMENU

  CLOSE(7)
  CALL GSCLOS
  STOP
END PROGRAM TASK_TX
