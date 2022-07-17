!***********************************************************************
!
!             TOKAMAK TRANSPORT SIMULATION CODE
!    INCLUDING RADIAL ELECTRIC FIELD AND PLASMA ROTATION
!                        "TASK/TX"
!           (TX : acronym of Transport Extended)
!
!      FEM scheme with Linear interpolation function
!      in the axisymmetric flux coordinates
!
!     DEVELOPED BY HONDA Mitsuru and FUKUYAMA Atsushi
!
!          ENGINEERING EDUCATION RESEARCH CENTER
!              GRADUATE SCHOOL OF ENGINEERING
!                     KYOTO UNIVERSITY
!
!               ADVANCED PLASMA MODELING GROUP
!                   NAKA FUSION INSTITUTE
!      NATIONAL INSTITUTES FOR QUANTUM AND RADIOLOGICAL
!                  SCIENCE AND TECHNOLOGY
!
!             DEPARTMENT OF NUCLEAR ENGINEERING
!              GRADUATE SCHOOL OF ENGINEERING
!                     KYOTO UNIVERSITY
!
!   References:
!     (Current version)
! [1] M. Honda and A. Fukuyama
!       "Development of the fluid-type transport code
!          on the flux coordinates in a tokamak",
!       Computer Physics Communications, 208 (2016) 117-134.
!     (Original)
! [2] M. Honda and A. Fukuyama
!       "Dynamic transport simulation code including plasma rotation
!          and radial electric field",
!       Journal of Computational Physics, 227 (2008) 2808-2844.
!     (With OFMC)
! [3] M. Honda, T. Takizuka, A. Fukuyama, M. Yoshida and T. Ozeki
!       "Self-consistent simulation of torque generation
!          by radial current due to fast particles"
!       Nuclear Fusion, 49 (2009) 035009 (10pp)
!     (Neutral model)
! [4] M. Honda, T. Takizuka, A. Fukuyama and K. Shimizu
!       "Modeling of neutral transport for self-consistent transport
!          simulations in tokamaks"
!       Journal of Plasma Fusion Research SERIES, 9 (2010) 529-534
!     (Radial electric field, jxB torque)
! [5] M. Honda, A. Fukuyama and N. Nakajima
!       "On the Neoclassical Relationship between the Radial Electric
!          Field and Radial Current in Tokamak Plasmas"
!       Journal of the Physical Society of Japan, 80 (2011) 114502 (14pp)
!     (Superstaged impurity)
! [6] A. Matsuyama
!       "Reduction of reaction rate equations of impurities" in Japanese
!       Private Communication, (2022)
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
!   Variables                                                            Boundary Cond.
!                                                                         axis   edge
!   X1  = phi                                                              N      0
!   X2  = psit'                                                            0      D 
!   X3  = psi'                                                             N      D
!   X4  = psi                                                              *      *
!   X5  = psit                                                             *      *
!   X6  = Ne               X14  = Ni               X22  = Nz               N      N
!   X7  = Ne * UerV        X15  = Ni * UirV        X23  = Nz * UzrV        0      *
!   X8  = <B Uepara>       X16  = <B Uipara>       X24  = <B Uzpara>       N(*)   N(*)
!   X9  = Ne * <R UePhi>   X17  = Ni * <R UiPhi>   X25  = Nz * <R UzPhi>   N      N
!   X10 = Ne * Te          X18  = Ni * Ti          X26  = Nz * Tz          N      N
!   X11 = <B qhatepara>    X19  = <B qhatipara>    X27  = <B qhatzpara>    *      *
!   X12 = Ne <UePhi/R>     X20 = Ni <UiPhi/R>      X28 = Nz <UzPhi/R>      *      *
!   X13 = B V1e            X21 = B V1i             X29 = B V1z             *      *
!   X30 = Nb                                                               *      *
!   X31 = Nb * UbrV                                                        0      *
!   X32 = <B Ubpara>                                                       *      *
!   X33 = Nb * <R UbPhi>                                                   *      *
!   X34 = Nb <UbPhi/R>                                                     *      *
!   X35 = B V1b                                                            *      *
!   X36 = Slow neutral                                                     N      D
!   X37 = Thermal neutral                                                  N      0
!   X38 = Halo Neutral                                                     N      0
!   X39 = Impurity neutral                                                 N      D
!   X40 = Nb ripple                                                        N      N
!
!   In case of MDFIXT /= 0, X10 = Te, X18 = Ti and X26 = Tz.
!   In case of MDBEAM == 0, X31 vanishes.
!
!   Nomenclature:
!     N : Neumann condition, 0 : Dirichlet condition with zero
!     D : Neumann or Dirichlet condition with finite value
!     * : No boundary condition imposed
  
program TASK_TX
  
  use tx_commons, only : SLID
  use tx_parameter_control, only : TXPARF, TXPARM_CHECK
  implicit none
  character(len=80) :: KPNAME

#ifndef nonGSAF
  call GSOPEN
#endif
  open(7,status='scratch',form='formatted')

  !     ***** Version ID *****
  !     SLID is used to identify data file.
  SLID = 'tx550.0'
  write(6,*) '######## TASK/TX V5.50.00 22/06/07 ########'

  call TXINIT
  KPNAME='txparm'
  call TXPARF(KPNAME)
  call TXPARM_CHECK

  call TXMENU

  close(7)
#ifndef nonGSAF
  call GSCLOS
#endif
end program TASK_TX
