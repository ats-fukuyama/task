! wminit.f90

MODULE wminit

  PRIVATE
  PUBLIC wm_init

CONTAINS

  SUBROUTINE wm_init

    USE wmcomm_parm
    IMPLICIT NONE
    INTEGER:: NA,i

!     *** MESH PARAMETERS ***

!     NRMAX  : Number of radial mesh points
!     NTHMAX : Number of poloidal mesh points (power of 2)
!     NHHMAX : Number of helically coupled toroidal modes (power of 2)
!                 =1 : axisymmetric calculation
!                 >1 : helical calculation
!     NPHMAX : Number of toroidal modes (power of 2)
!                 =1 : single toroidal mode calculation
!                 >1 : multi toroidal mode calculation (-NPHMAX/2+1..NPHMAX/2)
!                      NPHMAX >= NHHMAX*NHC
!     factor_nth : ratio of nthmax_f and nthmax [INTEGER]    
!     factor_nhh : ratio of nhhmax_f and nhhmax [INTEGER]    
!     factor_nph : ratio of nphmax_f and nphmax [INTEGER]    

      NRMAX   = 50
      NTHMAX  = 1
      NHHMAX  = 1
      NPHMAX  = 1
      factor_nth = 2 
      factor_nhh = 2
      factor_nph = 2

!     NSUMAX: Number of plasma surface plot points
!     NSWMAX: Number of wall surface plot points 
!     B0_FACT : B factor for equilibrum data

      NSUMAX  = 64
      NSWMAX  = 64
      B0_FACT = 1.D0
!     NTH0  : Central value of poloidal mode number
!     NPH0  : Central value of toroidal mode number
!     NHC   : Number of helical coils

      NTH0   = 0
      NPH0   = 8
      NHC    = 10

!     *** WAVE PARAMETER ***

!     RF    : Wave frequency                            (MHz)
!     RFI   : Wave growth rate                          (MHz)
!     RD    : Antenna minor radius                      (m)
!     NTH0  : Central value of poloidal mode number
!     NPH0  : Central value of toroidal mode number
!     NHC   : Number of helical coils
!     PRFIN : Input Power (0 for given antenna current) (W)

      RF     = 50.D0
      RFI    = 0.D0
      RD     = 1.1D0
      PRFIN  = 0.D0

!     NTH0  : Central value of poloidal mode number
!     NPH0  : Central value of toroidal mode number
!     NHC   : Number of helical coils

      NTH0   = 0
      NPH0   = 8
      NHC    = 10

!     *** ANTENNA PARAMETERS for multi mode analysis***

!        NAMAX : Number of vertical antennae
!        AJ    : Antenna current density                       (A/m)
!        AEWGT : Waveguide electric field (poloidal)           (V/m)
!        AEWGZ : Waveguide electric field (toroidal)           (V/m)
!        APH   : Antenna phase                              (degree)
!        THJ1  : Start poloidal angle of antenna            (degree)
!        THJ2  : End poloidal angle of antenna              (degree)
!        PHJ1  : Start toroidal angle of antenna            (degree)
!        PHJ2  : End toroidal angle of antenna              (degree)
!        ANTANG: Antenna angle: 0 for vertical antenna or perp WG
!        BETAJ : Antenna current profile parameter

      NAMAX  = 1
      DO NA=1,NAM
         AJ(NA)    = 1.D0
         AEWGT(NA) = 0.D0
         AEWGZ(NA) = 0.D0
         APH(NA)   = 0.D0
         THJ1(NA)  =-45.D0
         THJ2(NA)  = 45.D0
         PHJ1(NA)  = 0.D0
         PHJ2(NA)  = 0.D0
         ANTANG(NA)= 0.D0
         BETAJ(NA) = 0.D0
      ENDDO

!     **** ALPHA PARTICLE PARAMETERS ****

!        PNA   : Alpha density at center               (1.0E20/Mm*3)
!        PNAL  : Density scale length                            (m)
!        PTA   : Effective temperature                         (keV)

      PNA  = 0.02D0
      PNAL = 0.5D0
      PTA  = 3.5D3

!     **** ZEFF PARAMETERS ****

!     ZEFF  : Effective Z (sum n Z^2 / sum n Z)

      ZEFF  = 2.D0

!     *** CONTROL PARAMETERS ***

!        NPRINT: Control print output
!                   0: No print out
!                   1: Minimum print out (without input data)
!                   2: Minimum print out (with input data)
!                   3: Standard print out
!                   4: More print out
!                   5: and absorbed power detail
!        NGRAPH: Control graphic output
!                   0: No graphic out
!                   1: Standard graphic out (2D: Coutour)
!                   2: Standard graphic out (2D: Paint)
!                   3: Standard graphic out (2D: Birds eye)
!        MODELJ: Control antenna current model
!                   0: Loop antenna
!                   1: Waveguide
!                   2: Poloidal current
!                   3: Toroidal current
!                  2X: Vacuum eigen mode, poloidal current
!                  3X: Vacuum eigen mode, toroidal current
!        MODELA: Control alpha particle contribution
!                   0: No alpha effect
!                   1: Precession of alpha particles
!                   2: Precession of electrons
!                   3: Precession of both alpha particles and electrons
!                   4: Calculate alpha particle density using slowing down
!        MODELM: Control matrix solver
!                   0: Default
!        MDLWMK: toroidal effect of minimum parallel wave number
!                   0: OFF            
!                   1: ON
!        MDLWMX: model id
!                   0: wm
!                   1: wm_seki
!                   2: wmx
    

      NPRINT = 2
      NGRAPH = 1
      MODELJ = 0
      MODELA = 0
      MODELM = 0

      MDLWMK = 0
      MDLWMX = 2

!     *** EIGEN VALUE PARAMETERS ***

!        FRMIN : Minimum real part of frequency in amplitude scan
!        FRMAX : Maximum real part of frequency in amplitude scan
!        FIMIN : Minimum imag part of frequency in amplitude scan
!        FIMAX : Maximum imag part of frequency in amplitude scan
!        FI0   : Imag part of frequency in 1D amplitude scan

!        NGFMAX: Number of real freq mesh in 1D amplitude scan
!        NGXMAX: Number of real freq mesh in 2D amplitude scan
!        NGYMAX: Number of imag freq mesh in 2D amplitude scan

!        SCMIN : Minimum value in parameter scan
!        SCMAX : Maximum value in parameter scan
!        NSCMAX: Number of mesh in parameter scan
!        LISTEG: Listing in parameter scan

!        FRINI : Initial real part of frequency in Newton method
!        FIINI : Initial imag part of frequency in Newton method

!        DLTNW : Step size in evaluating derivatives in Newton method
!        EPSNW : Convergence criterion in Newton method
!        LMAXNW: Maximum iteration count in Newton method
!        LISTNW: Listing in Newton method
!        MODENW: Type of Newton method

!        NCONT : Number of contour lines
!        ILN1  : Line type of lower contours
!        IBL1  : Line boldness of lower contours
!        ICL1  : Line color of lower contours
!        ILN2  : Line type of higher contours
!        IBL2  : Line boldness of higher contours
!        ICL2  : Line color of higher contours

      FRMIN = 0.1D0
      FRMAX = 1.D0
      FIMIN =-0.1D0
      FIMAX = 0.1D0
      FI0   = 0.D0

      FRINI = RF
      FIINI = RFI

      NGFMAX= 11
      NGXMAX= 11
      NGYMAX= 11

      SCMIN = 0.1D0
      SCMAX = 1.D0
      NSCMAX= 11

      LISTEG= 1

      DLTNW = 1.D-6
      EPSNW = 1.D-6
      LMAXNW= 10
      LISTNW= 1
      MODENW= 0

      ILN1=0
      IBL1=2
      ICL1=7
      ILN2=0
      IBL2=2
      ICL2=7

      NCONT=30

!     *** ALFVEN FREQUENCY PARAMETERS ***

!        WAEMIN : Minimum frequency in Alfven frequency scan
!        WAEMAX : Maximum frequency in Alfven frequency scan

      WAEMIN = 0.001D0
      WAEMAX = 0.200D0

!    *** graphics ***

!        nthmax_g : number of poloidal mesh for graphics
!                   for modelg=4 or 6, nthmax_g is overwritten by nthmax

      nthmax_g = 256

!     *** dedub_info ***
!                41: wmemfp,wmsolv: nr1 and nr2 check
!                51: wmsetm1: IN/OUT mesh info
!                61: matrix coefficiewnts  knam_dump: 61+nrank
!                69: solution vector
!                71: nph0 mode dump
!                81: mode dependence of p_abs
      DO i=1,idebug_max
         idebuga(i)=0
      END DO
      knam_dump='wm.dump'
      
    RETURN
  END SUBROUTINE wm_init
END MODULE wminit
