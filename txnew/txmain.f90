!     $Id: txmain.f90,v 1.38 2011/04/13 05:50:25 honda Exp $
!
!***********************************************************************
!
!             TOKAMAK TRANSPORT SIMULATION CODE
!    INCLUDING RADIAL ELECTRIC FIELD AND PLASMA ROTATION
!                        "TASK/TX"
!
!      FEM scheme with Linear interpolation function
!      in the axisymmetric flux coordinates
!
!     DEVELOPED BY Mitsuru HONDA AND Atsushi FUKUYAMA
!
!               ADVANCED PLASMA MODELING GROUP
!                   NAKA FUSION INSTITUTE
!                JAPAN ATOMIC ENERGY AGENCY
!
!             DEPARTMENT OF NUCLEAR ENGINEERING
!              GRADUATE SCHOOL OF ENGINEERING
!                     KYOTO UNIVERSITY
!
!   References:
!     (Original)
!     M. Honda and A. Fukuyama
!       "Dynamic transport simulation code including plasma rotation
!          and radial electric field",
!       Journal of Computational Physics, 227 (2008) 2808-2844.
!     (With OFMC)
!     M. Honda, T. Takizuka, A. Fukuyama, M. Yoshida and T. Ozeki
!       "Self-consistent simulation of torque generation
!          by radial current due to fast particles"
!       Nuclear Fusion, 49 (2009) 035009 (10pp)
!     (Neutral model)
!     M. Honda, T. Takizuka, A. Fukuyama and K. Shimizu
!       "Modeling of neutral transport for self-consistent transport
!          simulations in tokamaks"
!       Journal of Plasma Fusion Research SERIES, 9 (2010) 529-534
!     (Radial electric field, jxB torque)
!     M. Honda, A. Fukuyama and N. Nakajima
!       "On the Neoclassical Relationship between the Radial Electric
!          Field and Radial Current in Tokamak Plasmas"
!       Journal of the Physical Society of Japan, 80 (2011) 114502 (14pp)
!
!   Obsolete references:
!     (Ripple model)
!     M. Honda, T. Takizuka, A. Fukuyama, M. Yoshida and T. Ozeki
!       "Numerical analysis of the effect of fast-ion losses
!          on plasma rotation in a tokamak with toroidal field ripple"
!       Nuclear Fusion, 48 (2008) 085003 (12pp)
!     (Anomalous particle transport model)
!     M. Honda, A. Fukuyama, T. Takizuka and K. Shimizu
!       "Modelling of anomalous particle transport
!          for dynamic transport simulations"
!       Nuclear Fusion, 50 (2010) 095012 (14pp)
!     (Neoclassical transport model)
!     M. Honda, A. Fukuyama and N. Nakajima
!       "Neoclassical Transport Modeling Compatible with a Two-Fluid
!          Transport Equation System"
!       Plasma and Fusion Research, 6 (2011) 1403008 (11pp)
!
!************************************************************************
!
!   Variables                                         Boundary Cond.
!                                                      axis   edge
!   X1  = phi                                           N      0
!   X2  = psit'                                         0      D 
!   X3  = psi'                                          N      D
!   X4  = psi                                           *      *
!   X5  = psit                                          *      *
!   X6  = Ne                 X13  = Ni                  N      N
!   X7  = Ne * UerV          X14  = Ni * UirV           0      *
!   X8  = <B Uepara>         X15  = <B Uipara>          N(*)   N(*)
!   X9  = Ne * <R UePhi>     X16  = Ni * <R UiPhi>      N      N
!   X10 = Ne * Te            X17  = Ni * Ti             N      N
!   X11 = <B qhatepara>      X18  = <B qhatipara>       *      *
!   X12 = Ne <UePhi/R>       X19 = Ni <UiPhi/R>         *      *
!   X20 = Nb                                            *      *
!   X21 = Nb * UbrV                                     0      *
!   X22 = <B Ubpara>                                    *      *
!   X23 = Nb * <R UbPhi>                                *      *
!   X24 = Nb <UbPhi/R>                                  *      *
!   X25 = Slow neutral                                  N      D
!   X26 = Thermal neutral                               N      0
!   X27 = Halo Neutral                                  N      0
!   X28 = Nb ripple                                     N      N
!
!   In case of MDFIXT /= 0, X10 = Te and X17 = Ti.
!   In case of MDBEAM == 0, X21 vanishes.
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
  SLID = 'tx500.0'
  WRITE(6,*) '######## TASK/TX V5.00.00 16/01/21 ########'

  CALL TXINIT
  KPNAME='txparm'
  CALL TXPARF(KPNAME)
  CALL TXPARM_CHECK

  CALL TXMENU

  CLOSE(7)
  CALL GSCLOS
END PROGRAM TASK_TX
