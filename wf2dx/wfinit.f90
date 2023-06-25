! wfinit.f90

MODULE WFINIT

! **** INITIALIZE AND UPDATE INPUT PARAMETERS ******

CONTAINS

!     ****** INITIALIZE PARAMETERS ******

  SUBROUTINE WF_INIT

    use wfcomm
    implicit none
    integer :: NA,NM,NS

!     *** CONTROL PARAMETERS ***
!
!        MODELI = Definition of imaginary unit
!              * 0 : exp(-i omega t)
!                1 : exp( j omega t)

!        MODELD = Dielectric Tensor Model
!              * 0 : Cold

!        MODELP = Density Profile
!                0 : Flat
!                1 : Step function with radius R/RA 
!              * 2 : Parabolic with radius R/RA

!   model_dielectric: model of dielectric
!                1 : cold plasma
!                2 : warm plasma (k_para ~ n_phi/R)
!                3 : kinetic plasma
    
    CALL pl_allocate_ns

    MODELI=0
    MODELG=2
    MODELD=0
    MODELP=2

    model_dielectric=1

!        KFNAME: File name of element data
!        KFNAMA: File name of antenna data
!        KFNAMF: File name of field data
!        KFNAMB: File name of buffer

    KFNAME = 'elm-data'
    KFNAMA = 'ant-data'
    KFNAMF = 'fld-data'
    KFNAMB = '/tmp/wfx-buffer'

!     *** INITIAL PARAMETERS FOR DIVIDE ***

!     NRM,NZM: PARAMETER FOR WFDIV

    NRM = 1001
    NZM = 1001

!     *** CONFIGURATION PARAMETERS ***

!        BB    : Magnetic field at center                        (T)
!        RR    : Plasma major radius                             (m)
!        RA    : Plasma minor radius                             (m)

    BB     = 0.072D0
    RR     = 0.22D0
    RA     = 0.16D0

!     *** CIRCULAR COIL/STRAIGHT ROD PARAMETERS (MODELB=1) ***
!
!        NCOILMAX : Number of coils
!        RCOIL : Radial position of coil current               (m)
!        ZCOIL : Axial position of coil current                (m)
!        BCOIL : Magnetic field on axis, center of coil        (T)

    NCOILMAX  = 0
    RCOIL(1)  = 0.35D0
    ZCOIL(1)  = 0.D0
    BCOIL(1)  = 0.001D0
    RCOIL(2)  = 0.35D0
    ZCOIL(2)  = 0.05D0
    BCOIL(2)  =-0.001D0
    RCOIL(3)  = 0.35D0
    ZCOIL(3)  =-0.05D0
    BCOIL(3)  =-0.001D0

!     *** RF PARAMETERS ***
!
!        RF    : Wave frequency                               (MHz)
!        NPH   : Toroidal Mode Number for modelg=1,2
!        RKZ   : Vertical wave number for modelg=0,12         (1/m)

    RF     = 5000.D0
    NPH    = 0
    RKZ    = 0.D0

!     *** ANTENNA PARAMETERS ***
!
!        NAMAX : Number of antennae
!        AJ    : Antenna current density                       (A/m)
!        APH   : Antenna phase                              (degree)
!        AWD   : Antenna width in (z, phi, Z) direction     (degree)
!        APOS  : Antenna position in (z, phi, Z) direction  (degree)

    NAMAX=0
    DO NA=1,NAM
       AJ(NA)   = 1.D0
       APH(NA)  = 0.D0
       AWD(NA)  = 0.D0
       APOS(NA) = 0.D0
    ENDDO

!     *** PLASMA PARAMETERS ***
!
!        NSMAX : Number of particle species
!        PA    : Mass number
!        PZ    : Charge number
!        PN    : Density at center                     (1.0E20/m**3)
!        PNS   : Density on plasma surface             (1.0E20/m**3)
!        PTPR  : Parallel temperature at center                (keV)
!        PTPP  : Perpendicular temperature at center           (keV)
!        PTS   : Temperature on surface                        (keV)
!        PZCL  : Ratio of collision frequency to wave frequency

    NSMAX= 2
    IF(NSMAX.GT.NSM) NSMAX=NSM
  
    PA(1)  = AME/AMP
    PZ(1)  =-1.0D0
    PN(1)  = 1.0D-3
    PNS(1) = 0.0D0
    PTPR(1)= 1.0D0
    PTPP(1)= 1.0D0
    PTS(1) = 0.1D0
    PZCL(1)= 3.0D-3
 
    IF(NSM.GE.2) THEN
       PA(2)  = 1.D0
       PZ(2)  = 1.0D0
       PN(2)  = 1.D-3
       PNS(2) = 0.0D0
       PTPR(2)= 1.D0
       PTPP(2)= 1.D0
       PTS(2) = 0.1D0
       PZCL(2)= 3.0D-3
    ENDIF
  
    DO NS=3,NSM
       PA(NSM)  = 1.0D0
       PZ(NSM)  = 1.0D0
       PN(NSM)  = 0.0D0
       PNS(NSM) = 0.0D0
       PTPR(NSM)= 1.0D0
       PTPP(NSM)= 1.0D0
       PTS(NSM) = 0.1D0
       PZCL(NSM)= 3.0D-3
    END DO

!   Edge damping region (absorption ns=NSMAX)
!      MDAMP=0 : no damping region
!           <>0: damping region in all direction except wg layer
!             1: no damp in the region r=rmin, zdamp_min<z<zdamp_max 
!             2: no damp in the region r=rmax, zdamp_min<z<zdamp_max 
!             3: no damp in the region z=zmin, rdamp_min<r<rdamp_max 
!             4: no damp in the region z=zmax, rdamp_min<r<rdamp_max 
!      WDAMP   : width of damping region [m]
!      FDAMP   : damping factor
!                       epsilon=FDAMP*(WDAMP-DX)/(DX+CI*PZCL(NSMAX))
!      PZCL(NSMAX): damping factor: PN(NSMAX)=0

    MDAMP=0
    WDAMP=0.D0
    FDAMP=0.3D0

    ! *** collision enhancement near a layer ***
    !  model_coll_enhance:  0: no enhancement, 1: layer in x, 2: layer in y
    !  factor_coll_enhance: enhancement factor (0.D0: no enhancement)
    !  xpos_coll_enhance:   center position of the layer in x (R)
    !  xwidth_coll_enhance: Gaussian width in x (R)
    !  ypos_coll_enhance:   center position of the layer in y (Z)
    !  ywidth_coll_enhance: Gaussian width in y (Z)
    !       rzcl=rzcl_original*(1.D0+factor*exp(-(x-xpos)**2/xwidth**2))

    model_coll_enhance=0
    factor_coll_enhance=0.D0
    xpos_coll_enhance  =0.D0
    xwidth_coll_enhance=0.01D0
    ypos_coll_enhance  =0.D0
    ywidth_coll_enhance=0.01D0

  !     *** MEDIUM PARAMETERS ***

    call wfmed_allocate

    NMMAX=0
    DO NM=1,NMM
       EPSDM(NM)=1.D0
       AMUDM(NM)=1.D0
       SIGDM(NM)=0.D0
    ENDDO
  
    NBMAX=0

!     *** OUTPUT PARAMETERS ***
!
!        NPRINT: Print output parameter
!                0 : No output
!              * 1 : Parameter and global field data
!                2 : Local field data
!                3 : Element data
!        NDRAWD: Drawing parameter for elemendt divider
!                0 : Boundary shape
!              * 1 : Element shape
!                2 : Element shape + Element number
!                3 : Element shape + Element number + Node number
!        NDRAWA: Drawing parameter for antenna generater
!                0 : Antenna primary data
!                1 : Antenna secondary data
!              * 2 : Antenna secondary data + Element shape
!        NDRAWE: Drawing parameter for electric field profile
!              * 0 : cylinder (r,phi,z)
!                1 : toroidal (r,theta,phi)
!        NDRAWV: Vector field output parameter
!              * 0 : No output
!                1 : 2D Vector field in poloidal cross section
!        NGRAPH: Drawing parameter for field data
!               -1 : binary file output
!                0 : text file output
!              * 1 : Contour plot
!                2 : Paint plot
!                3-6 : Bird's eye view from four directions

    NPRINT = 1
    NDRAWD = 1
    NDRAWA = 2
    NDRAWE = 0
    NDRAWV = 0
    NGRAPH = 1

!     *** variable_sort parameters ***

   sort_weight_x=1.D0
   sort_weight_y=1.D0

!     *** DIVIDER PARAMETERS ***

!        RB    : Boundary radius (m)
!        DELR  : Typical element size in r direction (m)
!        DELZ  : Typical element size in z direction (m)

    RB     = 0.2D0
    DELR   = 0.05D0
    DELZ   = 0.05D0
    BDRMIN = 0.01d0
    BDRMAX = 0.5d0
    BDZMIN =-0.5d0
    BDZMAX = 0.5d0

!     *** ANTENNA SHAPE PARAMETERS ***

!        PIN   : Input Power (W)
!        RD    : Antenna radius (m)
!        THETJ1: Start angle of arc antenna (degree)
!        THETJ2: End angle of arc antenna (degree)
!        NJMAX : Number of primary grid points of antenna

    PIN    = 1.D0
    RD     = 0.1855D0
    THETJ1 = 40.D0
    THETJ2 =-40.D0
    NJMAX  = 41

    MODELWG=1         ! Profile: 0: step function, 1: gaussian
                      !         12: read file
    R1WG=0.5D0
    Z1WG=0.0D0
    R2WG=0.5D0
    Z2WG=0.05D0
    PH1WG=0.D0
    PH2WG=180.D0
    AMPWG=0.D0        ! Amplitude of waveguide electric field
    ANGWG=0.D0        ! Polarization angle wrt z axis
    ELPWG=0.D0        ! Polarization ellipticity
    DPHWG=0.D0        ! Phase curvature for focusing

!     PPN0   : Neutral pressure (Pa)  1 Torr = 1 mmHg = 133.322 Pa
!     PTN0   : Initial neutral temperarure (eV)

    PPN0   = 3.0D0
    PTN0   = 0.03D0

!     *** GRAPHICS PARAMETER ***

    NGXMAX = 101
    NGYMAX = 101
    NGVMAX = 101  

    GFACTOR= 0.5

    nxzone_max=100
    nyzone_max=100

!     *** Numerical computation parameter ***

    tolerance = 1.D-8

!     *** DEBUG CONTROL PARAMETER ***
    ! --- idebug=3 : wg e-field output ---

    IDEBUG = 0

!     *** LOAD FILE NAME ***

    KNAMPF=' '
    KNAMWG=' '

!     *** INITIALIZATION PARAMETERS (DO NOT MODIFY) ***

    NNMAX  = 0
    nemax=0
    nsdmax=0
    mlen=0
    rdamp_min=0.D0
    rdamp_max=0.D0
    zdamp_min=0.D0
    zdamp_max=0.D0

    NFOPEN = 0
    RNDMIN = 0.D0
    RNDMAX = 0.D0

    NDFILE=25

    MODELWF=0         ! side field: 0: positive 1: alternative

    RETURN
END SUBROUTINE WF_INIT
END MODULE WFINIT
