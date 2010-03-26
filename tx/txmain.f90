!     $Id$
!
!***********************************************************
!
!             TOKAMAK TRANSPORT SIMULATION CODE
!    INCLUDING RADIAL ELECTRIC FIELD AND PLASMA ROTATION
!                        "TASK/TX"
!
!      FEM scheme with Linear interpolation function
!      in the r^2 coordinates
!
!     DEVELOPED BY M. HONDA, Y. FUJI AND A. FUKUYAMA
!
!                JAPAN ATOMIC ENERGY AGENCY
!                   NAKA FUSION INSTITUTE
!                 JT-60 PLASMA DESIGN GROUP
!
!             DEPARTMENT OF NUCLEAR ENGINEERING
!              GRADUATE SCHOOL OF ENGINEERING
!                     KYOTO UNIVERSITY
!
!     References:
!       M. Honda and A. Fukuyama
!         "Dynamic transport simulation code including plasma rotation
!            and radial electric field",
!         Journal of Computational Physics, 227 (2008) 2808-2844.
!       M. Honda, T. Takizuka, A. Fukuyama, M. Yoshida and T. Ozeki
!         "Numerical analysis of the effect of fast-ion losses
!            on plasma rotation in a tokamak with toroidal field ripple"
!         Nuclear Fusion, 48 (2008) 085003 (12pp)
!       M. Honda, T. Takizuka, A. Fukuyama, M. Yoshida and T. Ozeki
!         "Self-consistent simulation of torque generation
!            by radial current due to fast particles"
!         Nuclear Fusion, 49 (2009) 035009 (10pp)
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
!   X20 = Thermal neutral                               N      0
!   X21 = Halo Neutral                                  N      0
!   X22 = Nb ripple                                     N      N
!
!   In case of MDFIXT /=0, X10 = Te, X15 = Ti
!
!   Nomenclature:
!     N : Neumann condition, 0 : Dirichlet condition with zero
!     D : Neumann or Dirichelt condition with finite value
!     * : No boundary condition imposed
  

PROGRAM TASK_TX
  
  use tx_commons, only : SLID
  use tx_parameter_control, only : TXPARF, TXPARM_CHECK
  implicit none
  character(len=80) :: KPNAME

  CALL GSOPEN
  OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')

  !     ***** Version ID *****
  !     SLID is used to identify data file.
  SLID = 'tx458.0'
  WRITE(6,*) '######## TASK/TX V4.58.00 10/03/26 ########'

  CALL TXINIT
  KPNAME='txparm'
  CALL TXPARF(KPNAME)
  CALL TXPARM_CHECK

  CALL TXMENU

  CLOSE(7)
  CALL GSCLOS
  STOP
END PROGRAM TASK_TX
