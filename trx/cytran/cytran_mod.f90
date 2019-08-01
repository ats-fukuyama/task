MODULE CYTRAN_MOD
!-------------------------------------------------------------------------------
!CYTRAN-CYclotron radiation TRANsport
!
!CYTRAN_MOD is an F90 module that calculates cyclotran radiation transport in
!  toroidal plasmas
!
!References:
!
!  S.Tamor, SAIC report SAI-023-81-110LJ/LAPS-71 (1981)
!  S.Tamor, SAIC report SAI-023-81-189LJ/LAPS-72 (1981)
!  D.C.Baxter, S.Tamor, SAIC report SAI-023-81-215LJ/LAPS-74 (1981)
!  D.C.Baxter, S.Tamor, Nucl Technol/Fusion 3 (1983) 181
!  W.A.Houlberg, F90 free format 8/2004
!
!Contains PUBLIC routine:
!
!  CYTRAN              -cyclotron radiation transport
!
!Comments:
!
!  CYTRAN is an approximation to the cyclotron radiation transport process that
!    requires benchmarking against more complete transport analysis to improve
!    confidence in its results.
!
!  The magnetic geometry is simplified by assuming a circular cross section in a
!    plane of constant toroidal angle with an average magnetic field that is
!    only a function of minor radius.
!
!  Tamor benchmarked CYTRAN against the more complete cyclotron transport
!    analysis of the SNECTR code.
!
!  More recent benchmarking of CYTRAN against another cyclotron radiation code
!    has been carried out by:
!    F.Albajar, M.Bornatici, F.Engelmann, submitted to Nucl Fusion (2001)
!
!  Because the transport is expressed in terms of plasma surface area and
!    volumes, it can be applied to non-circular axisymmetric plasmas while
!    preserving the energy balance; the same is true if applied to
!    non-axisymmetric plasmas, but these must be viewed with more caution
!    because there have been no benchmark calculations for plasmas in which
!    field curvature in the toroidal direction may be more important.
!
!  In many places, physical and empirical constants are merged - this was done
!    in the original coding - and have not been unraveled in the effort to
!    convert to a more modern and cleaner code.
!
!  The modernization of the code structure into an F90 module takes advantage of
!    some of the more attractive features of F90:
!    -use of KIND for precision declarations
!    -optional arguments for I/O
!    -generic names for all intrinsic functions
!    -compilation using either free or fixed form
!    -no common blocks or other deprecated Fortran features
!    -dynamic and automatic alocation of variables
!    -array syntax for vector operations
!-------------------------------------------------------------------------------
!USE SPEC_KIND_MOD
USE bpsd_kinds
IMPLICIT NONE

!-------------------------------------------------------------------------------
! Private procedures
!-------------------------------------------------------------------------------
PRIVATE :: &
  CYTRAN_OPACITY         !cyclotron opacity coefficients
                         !  called from CYTRAN

!-------------------------------------------------------------------------------
! Private data
!-------------------------------------------------------------------------------
!Constants
REAL(KIND=rkind), PRIVATE, PARAMETER :: &
  z_pi=3.141592654       !pi [-]

!-------------------------------------------------------------------------------
! Public procedures
!-------------------------------------------------------------------------------
CONTAINS

SUBROUTINE CYTRAN(ree,reo,roo,roe,nr_r,bmod_r,den_r,te_r,area_rm,dvol_r, &
                  psync_r, &
                  K_CYT_RES)
!-------------------------------------------------------------------------------
!CYTRAN estimates cyclotron radiation energy transport in toroidal plasmas
!
!References:
!  S.Tamor, SAIC report SAI-023-81-110LJ/LAPS-71 (1981)
!  S.Tamor, SAIC report SAI-023-81-189LJ/LAPS-72 (1981)
!  D.C.Baxter, S.Tamor, SAIC report SAI-023-81-215LJ/LAPS-74 (1981)
!  D.C.Baxter, S.Tamor, Nucl Technol/Fusion 3 (1983) 181
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  nr_r                   !no. of radial plasma cells [-]

REAL(KIND=rkind), INTENT(IN) :: &
  ree,                 & !refl of X mode from incident X mode [0-1]
  reo,                 & !refl of X mode from incident O mode [0-1]
  roe,                 & !refl of O mode from incident X mode [0-1]
  roo                    !refl of O mode from incident O mode [0-1]

REAL(KIND=rkind), INTENT(IN) :: &
  bmod_r(:),           & !<|B|> [T]
  den_r(:),            & !electron density in cell i [/m**3]
  dvol_r(:),           & !cell volume [m**3]
  te_r(:),             & !electron temperature in cell i [keV]
  area_rm(:)             !surface area at inner boundary of cell [m**2]
                         !plasma outer surface area is at nr_r+1


!Declaration of optional input variables
INTEGER, INTENT(IN), OPTIONAL :: &
  K_CYT_RES              !frequency resolution level [-]
                         !=1 delta(omega)=default of 100 frequncy intervals
                         !>1 100*K_CYT_RES frequency intervals
                         !10 recommended to keep graininess < a couple %


!Declaration of output variables
REAL(KIND=rkind), INTENT(OUT) :: &
  psync_r(:)             !net power source/loss in cell i [keV/m**3/s]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,j,k_res,nze,nzo

REAL(KIND=rkind) :: &
  areae,areao,blackb,delw,dwght,fwght,emass,omass,extray,ordray, &
  qee,qeo,qoe,qoo,srce,srco,taucrt,taue,tauo,testi,tse,tso,wfreq, &
  wght,wmax,wmin,xe,xo

REAL(KIND=rkind) :: &
  alphae(1:nr_r),alphao(1:nr_r),delta(1:nr_r),temp(1:nr_r)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Frequency resolution
k_res=1

IF(PRESENT(K_CYT_RES)) THEN

  IF(K_CYT_RES > 1) k_res=K_CYT_RES

ENDIF

!Minimum and maximum frequencies and interval
wmin=5.0e11*MINVAL(bmod_r(1:nr_r))
wmax=2.4e12*MAXVAL(bmod_r(1:nr_r))
delw=0.08*(wmax-wmin)/k_res

!Other
dwght=1
fwght=3
taucrt=0.6
testi=0
extray=1.0e36

psync_r(1:nr_r)=0
temp(1:nr_r)=te_r(1:nr_r)*dvol_r(1:nr_r)

DO i=1,nr_r-1 !Over radial zones

  delta(i)=2*dvol_r(i)/(area_rm(i)+area_rm(i+1))

ENDDO !Over radial zones

delta(nr_r)=dvol_r(nr_r)/area_rm(nr_r)

DO j=1,100*k_res !Over frequency

  IF(extray < 0.03*testi) EXIT

  wfreq=wmin+delw*REAL(j-1,rkind)
  dwght=-dwght
  wght=fwght+dwght
  IF(j == 1) wght=1
  wght=wght*delw/3

  !Blackbody radiation factor (/m**2/s)
  blackb=4.48e-20*wfreq*wfreq

  !Get opacity factors
  CALL CYTRAN_OPACITY(nr_r,te_r,bmod_r,wfreq,alphae,alphao)
  alphae(1:nr_r)=alphae(1:nr_r)*6.03e-17*den_r(1:nr_r)/bmod_r(1:nr_r)
  alphao(1:nr_r)=alphao(1:nr_r)*6.03e-17*den_r(1:nr_r)/bmod_r(1:nr_r)

  taue=0
  tauo=0
  nze=nr_r
  nzo=nr_r

  DO i=nr_r,1,-1 !Over radial zones

    taue=taue+alphae(i)*delta(i)
    IF(taue < taucrt) nze=i-1
    tauo=tauo+alphao(i)*delta(i)
    IF(tauo < taucrt) nzo=i-1

  ENDDO !Over radial zones 

  !Extraordinary mode
  IF(nze > 0) THEN

    areae=area_rm(nze+1)
    srce=areae*te_r(nze)
    emass=SUM(den_r(1:nr_r)*dvol_r(1:nr_r))

  ELSE

    areae=0
    srce=0
    emass=0

  ENDIF

  !Ordinary mode
  IF(nzo > 0) THEN

    areao=area_rm(nzo+1)
    srco=areao*te_r(nzo)
    omass=SUM(den_r(1:nr_r)*dvol_r(1:nr_r))

  ELSE

    areao=0
    srco=0
    omass=0

  ENDIF

!Extraordinary mode
  IF(nr_r > nze) THEN

    xe=4*SUM(dvol_r(nze+1:nr_r)*alphae(nze+1:nr_r))
    srce=srce+4*SUM(temp(nze+1:nr_r)*alphae(nze+1:nr_r))

  ELSE

    xe=0

  ENDIF

!Ordinary mode
  IF(nr_r > nzo) THEN

    xo=4*SUM(dvol_r(nzo+1:nr_r)*alphao(nzo+1:nr_r))
    srco=srco+4*SUM(temp(nzo+1:nr_r)*alphao(nzo+1:nr_r))

  ELSE

    xo=0

  ENDIF

  qee=areao+area_rm(nr_r+1)*(1.0-roo)+xo
  qoo=areae+area_rm(nr_r+1)*(1.0-ree)+xe
  qeo=reo*area_rm(nr_r+1)
  qoe=roe*area_rm(nr_r+1)
  extray=(qee*srce+qeo*srco)*blackb/(qee*qoo-qeo*qoe)
  ordray=(qoe*srce+qoo*srco)*blackb/(qee*qoo-qeo*qoe)
  testi=MAX(testi,extray)

  DO i=1,nr_r !Over radial zones

    !Extraordinary mode
    IF(i <= nze) THEN

      tse=z_pi*areae*(blackb*te_r(nze)-extray)*(den_r(i)/emass)*wght

    ELSE

      tse=4*z_pi*alphae(i)*(blackb*te_r(i)-extray)*wght

    ENDIF

    psync_r(i)=psync_r(i)-tse

    !Ordinary mode
    IF(i <= nzo) THEN

      tso=z_pi*areao*(blackb*te_r(nzo)-ordray)*(den_r(i)/omass)*wght

    ELSE

      tso=4*z_pi*alphao(i)*(blackb*te_r(i)-ordray)*wght

    ENDIF

    psync_r(i)=psync_r(i)-tso

  ENDDO !Over radial zones

ENDDO !Over frequency

END SUBROUTINE CYTRAN

SUBROUTINE CYTRAN_OPACITY(nr_r,te_r,bmod_r,wfreq,alphae,alphao)
!-------------------------------------------------------------------------------
!CYTRAN_OPACITY determines the opacity coefficients for cyclotron radiation
!
!References:
!  S.Tamor, SAIC report SAI-023-81-110LJ/LAPS-71 (1981)
!  S.Tamor, SAIC report SAI-023-81-189LJ/LAPS-72 (1981)
!  D.C.Baxter, S.Tamor, SAIC report SAI-023-81-215LJ/LAPS-74 (1981)
!  D.C.Baxter, S.Tamor, Nucl Technol/Fusion 3 (1983) 181
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  nr_r                   !no. of radial cells [-]

REAL(KIND=rkind), INTENT(IN) :: &
  wfreq                  !wave angular frequency [radians/s]

REAL(KIND=rkind), INTENT(IN) :: &
  te_r(:),             & !electron temperature [keV]
  bmod_r(:)              !<|B|> [T]

!Declaration of output variables
REAL(KIND=rkind), INTENT(OUT) :: &
  alphae(:),           & !extra-ordinary mode opacity [-]
  alphao(:)              !ordinary mode opacity [-]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i

REAL(KIND=rkind) :: &
  arg,tden,wloc

DO i=1,nr_r !Over zones

  !wloc=wfreq/electron cyclotron frequency
  wloc=wfreq/(1.7588e+11*bmod_r(i))
  tden=1.0_dp/(4.0+te_r(i)+25.0/te_r(i))
  arg=MAX(0.0_dp,0.045_dp+(wloc-2.0_dp)*tden)
  alphae(i)=10.0**(1.45-7.8*SQRT(arg))/(wloc**2)
  arg=MAX(0.0_dp,0.180_dp+(wloc-1.0_dp)*tden)
  alphao(i)=10.0**(2.45-8.58*SQRT(arg))/(wloc**2)

ENDDO !Over zones

END SUBROUTINE CYTRAN_OPACITY

END MODULE CYTRAN_MOD
      
