! wfinit.f90

MODULE WFINIT

  ! **** INITIALIZE AND UPDATE INPUT PARAMETERS ******

  PRIVATE
  PUBLIC wf_init

CONTAINS

!     ****** INITIALIZE PARAMETERS ******

  SUBROUTINE WF_INIT

    USE wfcomm
    USE plcomm
    implicit none
    integer :: nant,ns,nmed,id

    ! === Configuration parameters for wfdiv ===

    model_config = 1 ! coordinates type
    !                   1 : rectangular (x,y)
    !                   2 : torus (r,z)
    !                   3 : cylinder (z,r)
    
    model_shape = 1  ! boundary shape type
    !                   1 : rectabgular
    !                   2 : circle
    !                   3 : arc
    !                  11 : toroidal equilibrium (eq)
    !                  12 : rectangular mesh file

    ! --- model shape = 1 (rectangular) ---
    
    xdiv_min=0.D0
    xdiv_max=1.D0
    ydiv_min=0.D0
    ydiv_max=1.D0

    delx=0.01D0  ! typical division length
    dely=0.01D0  ! typical division length

    ! --- model shape = 2 (circular) --- DEFINED in plcomm

    ! RR=3.D0      ! major radius
    ! RB=1.2D0     ! wall radius
    ! RKAP=1.D0    ! elongation
    ! RDLT=0.D0    ! triangularity

    ! --- model shape = 3 (arc) ---
    
    rdiv_min=1.5D0
    rdiv_max=4.5D0
    thdiv_min=0.D0
    thdiv_max=90.D0

    ! === magnetic field parameter ===

    ! --- typcal magnetic field ---  DEFINED in pl/plcomm

    ! BB = 1.D0
    ! Q0 = 1.D0
    ! QA = 3.D0
    ! RIP   = 3.D0
    ! PROFJ = 2.D0

    ! --- coil parameters for cylindrical configuration --- DEFINED in plcomm

    !        ncoil_max : Number of coils
    !        rcoil(ncoilm) : r position of coil current [m]
    !        zcoil(ncoilm) : z position of coil current [m]
    !        bcoil(ncoim)  : maximum magnetic field on axis [T]

    ! ncoil_max = 3
    ! rcoil(1)  = 0.35D0
    ! zcoil(1)  = 0.D0
    ! bcoil(1)  = 0.001D0
    ! rcoil(2)  = 0.35D0
    ! zcoil(2)  = 0.05D0
    ! bcoil(2)  =-0.001D0
    ! rcoil(3)  = 0.35D0
    ! zcoil(3)  =-0.05D0
    ! bcoil(3)  =-0.001D0

    ! === RF PARAMETERS ===

    !        RF    : Wave frequency                               (MHz)
    !        NPH   : Toroidal Mode Number for modelg=1,2
    !        RKZ   : Vertical wave number for modelg=0,12         (1/m)

    RF     = 5000.D0
    RKZ    = 0.D0
    NPH    = 0

    ! --- ANTENNA PARAMETERS ---
    !
    !        nant_max : Number of antennae
    !        AJ(nantm)    : Antenna current density                       (A/m)
    !        APH(nantm)   : Antenna phase                              (degree)
    !        AWD(nantm)   : Antenna width in (z, phi, Z) direction     (degree)
    !        APOS(nantm)  : Antenna position in (z, phi, Z) direction  (degree)
    !        PIN(AM)   : Input Power (W)
    !        RD(nantm)    : Antenna radius (m)
    !        THETJ1(nantm): Start angle of arc antenna (degree)
    !        THETJ2(nantm): End angle of arc antenna (degree)
    !        nant_point_max(nantm) : Number of primary grid points of antenna

    nant_max=0
    DO nant=1,nantm
       AJ(nant)            = 1.D0
       APH(nant)           = 0.D0
       AWD(nant)           = 0.D0
       APOS(nant)          = 0.D0
       PIN(nant)           = 1.D0
       RD(nant)            = 0.1855D0
       THETJ1(nant)        = 40.D0
       THETJ2(nant)        =-40.D0
    ENDDO

    ! === waveguide antenna parameter ===

    model_wg=1   ! antenna model
    !              0   : line, step
    !              1   : line, gaussian
    !              6   : arc, step
    !              7   : arc, gaussian
    !             12:    read file
    
    xwg_min=0.5D0        ! waveguide range in xy
    xwg_max=0.5D0
    ywg_min=0.0D0
    ywg_max=0.05D0
    thwg_min=-30.D0      ! waveguide range in poloidal angle
    thwg_max= 30.D0
    
    phase_wg_min=0.D0    ! wave phase
    phase_wg_cen=90.D0
    phase_wg_max=180.D0

    amp_wg=0.D0      ! Amplitude of waveguide electric field
    angle_wg=0.D0    ! Polarization angle wrt z axis (degree)
    ellip_wg=0.D0    ! Polarization ellipticity

    !   Edge damping region (absorption ns=NSMAX)
    !      model_damp = 0 : no damping region
    !                damping region in all direction except waveguide mouth
    !             1: no damp on the wall x=xdiv_min, ydamp_min<y<ydamp_max 
    !             2: no damp on the wall x=xdiv_max, ydamp_min<y<ydamp_max 
    !             3: no damp on the wall y=ydiv_min, xdamp_min<x<xdamp_max 
    !             4: no damp on the wall y=ydiv_max, xdamp_min<x<rdamp_max
    !             5: no damp on the wall RS=RB, thdamp_min<theta<thdamp_max
    !      width_mdap  : width of damping region [m]
    !      factor_damp : damping factor (0 < DX=|x-wall| < width)
    !                       epsilon=factor*(width-DX)/(DX+CI*PZCL(NSMAX))

    model_damp=0
    xdamp_min=0.D0
    xdamp_max=0.D0
    ydamp_min=0.D0
    ydamp_max=0.D0
    thdamp_min=0.D0
    thdamp_max=0.D0

    width_damp=0.D0
    factor_damp=0.3D0

    ! === PLASMA PARAMETERS ===

    !     NSMAX : Number of particle species
    !                  (when model_damp.NE.0, ns=nsmax for wall damp)
    !     PA    : Mass number
    !     PZ    : Charge number
    !     PN    : Density at center                     (1.0E20/m**3)
    !     PNS   : Density on plasma surface             (1.0E20/m**3)
    !     PTPR  : Parallel temperature at center                (keV)
    !     PTPP  : Perpendicular temperature at center           (keV)
    !     PTS   : Temperature on surface                        (keV)
    !     PZCL  : Ratio of collision frequency to wave frequency

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
       PA(NS)  = 1.0D0
       PZ(NS)  = 1.0D0
       PN(NS)  = 0.0D0
       PNS(NS) = 0.0D0
       PTPR(NS)= 1.0D0
       PTPP(NS)= 1.0D0
       PTS(NS) = 0.1D0
       PZCL(NS)= 3.0D-3
    END DO

    !  PPN0   : Neutral pressure (Pa)  1 Torr = 1 mmHg = 133.322 Pa
    !  PTN0   : Initial neutral temperarure (eV)

    PPN0   = 3.0D0
    PTN0   = 0.03D0

    ! *** collision enhancement near a layer ***
    !  model_coll_enhance:  0: no enhancement, 1: layer in x, 2: layer in y
    !  factor_coll_enhance: enhancement factor (0.D0: no enhancement)
    !  xpos_coll_enhance:   center position of the layer in x
    !  xwidth_coll_enhance: Gaussian width in x
    !  ypos_coll_enhance:   center position of the layer in y
    !  ywidth_coll_enhance: Gaussian width in y
    !       rzcl=rzcl_original*(1.D0+factor*exp(-(x-xpos)**2/xwidth**2))

    model_coll_enhance=0
    factor_coll_enhance=0.D0
    xpos_coll_enhance  =0.D0
    xwidth_coll_enhance=0.01D0
    ypos_coll_enhance  =0.D0
    ywidth_coll_enhance=0.01D0

    ! *** interpolation level ***
    !   model_interpolation: 0: element center value
    !                        1: local linear interpolation

    model_interpolation=1

    ! *** MEDIUM PARAMETERS ***
    !   nmed_max : maximum number of medium type
    !   epsilon_nmed : electric permittivity
    !   amu_nmed     : magnetic permeability
    !   sig_nmed     : electric conductivity

    nmed_max=0
    DO nmed=1,nmedm
       model_nmed(nmedm)=0
       epsilon_nmed(nmed)=1.D0
       amu_nmed(nmed)=1.D0
       sigma_nmed(nmed)=0.D0
       xmin_nmed(nmed)=0.D0
       xmax_nmed(nmed)=0.D0
       ymin_nmed(nmed)=0.D0
       ymax_nmed(nmed)=0.D0
       rmin_nmed(nmed)=0.D0
       rmax_nmed(nmed)=0.D0
       thmin_nmed(nmed)=0.D0
       thmax_nmed(nmed)=0.D0
    ENDDO
  
    ! === OUTPUT PARAMETERS ===

    !   NPRINT: Print output parameter
    !            0 : No output
    !          * 1 : Parameter and global field data
    !            2 : Local field data
    !            3 : Element data
    !   NDRAWD: Drawing parameter for elemendt divider
    !            0 : Boundary shape
    !          * 1 : Element shape
    !            2 : Element shape + Element number
    !            3 : Element shape + Element number + Node number
    !   NDRAWA: Drawing parameter for antenna generater
    !            0 : Antenna primary data
    !            1 : Antenna secondary data
    !          * 2 : Antenna secondary data + Element shape
    !   NDRAWE: Drawing parameter for electric field profile
    !          * 0 : cylinder (r,phi,z)
    !            1 : toroidal (r,theta,phi)
    !   NDRAWV: Vector field output parameter
    !          * 0 : No output
    !            1 : 2D Vector field in poloidal cross section
    !   NGRAPH: Drawing parameter for field data
    !           -1 : binary file output
    !            0 : text file output
    !          * 1 : Contour plot
    !            2 : Paint plot
    !            3-6 : Bird's eye view from four directions

    nprint = 1
    ndrawd = 1
    ndrawa = 2
    ndrawe = 0
    ndrawv = 0
    ngraph = 1

    ! === variable_sort weight parameters ===

    sort_weight_x=1.D0
    sort_weight_y=1.D0


    ! === GRAPHICS PARAMETER ===

    ngxmax = 31
    ngymax = 31
    ngvmax = 31  

    gaspect= 1.0

    nxzone_max=100
    nyzone_max=100

    ! === Numerical computation parameter ===

    tolerance = 1.D-8

    ! === seg field parity ===

    model_wf = 0

    ! === CONTROL PARAMETERS ===

    !   MODELI = Definition of imaginary unit
    !      * 0 : exp(-i omega t)
    !        1 : exp( j omega t)


    MODELI=0

    !   KFNAME: File name of element data
    !   KFNAMA: File name of antenna data
    !   KFNAMF: File name of field data
    !   KFNAMB: File name of buffer

    KFNAME = 'elm-data'
    KFNAMA = 'ant-data'
    KFNAMF = 'fld-data'
    KFNAMB = '/tmp/wfx-buffer'

!     *** LOAD FILE NAME ***

    KNAMWG=' '

    ! === DEBUG CONTROL PARAMETER ===
    
    ! --- idebuga( 1) : wfdiv
    ! --- idebuga( 2) : wfant
    ! --- idebuga( 3) : wg e-field output ---
    ! --- idebuga( 4) : wfindex

    DO id=1,idebuga_max
       idebuga(id)=0
    END DO

    RETURN
  END SUBROUTINE WF_INIT
END MODULE WFINIT
