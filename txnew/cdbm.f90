!     ********************************************

!           CDBM Transport model (2009/03/06)
!              Modified by M. Honda (2016/05/12)

!     ********************************************

MODULE cdbm_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC:: cdbm

  integer,parameter :: rkind=selected_real_kind(12,100)
  integer,parameter :: ikind=selected_int_kind(8)
  integer,parameter :: dp=rkind

  real(rkind),parameter :: ZERO = 0.0_dp
  real(rkind),parameter :: HALF = 0.5_dp
  real(rkind),parameter :: ONE  = 1.0_dp
  real(rkind),parameter :: TWO  = 2.0_dp

  real(rkind),parameter :: PI   = 3.14159265358979323846_dp
  real(rkind),parameter :: TWOPI= PI+PI

  ! Physical constans, based on CODATA 2006
  real(rkind),parameter :: AEE  = 1.602176487E-19_dp ! elementary charge
  real(rkind),parameter :: AME  = 9.10938215E-31_dp  ! electron mass
  real(rkind),parameter :: AMP  = 1.672621637E-27_dp ! proton mass
  real(rkind),parameter :: VC   = 2.99792458E8_dp  ! speed of light
  real(rkind),parameter :: RMU0 = 4.E-7_dp*PI      ! permeability
  real(rkind),parameter :: EPS0 = ONE/(VC*VC*RMU0) ! permittivity

CONTAINS

  SUBROUTINE cdbm(bb,rr,rs,rkap,qp,shear,pne,rhoni,dpdr,dvexbdr, &
       &             calf,ckap,cexb,model,chi_cdbm,fsz,curvz,fez,omgexb)

    real(rkind),intent(in):: bb      ! Magnetic field strength [T]
    real(rkind),intent(in):: rr      ! Major radius [m]
    real(rkind),intent(in):: rs      ! Minor radius [m]
    real(rkind),intent(in):: rkap    ! elongation
    real(rkind),intent(in):: qp      ! Safety factor
    real(rkind),intent(in):: shear   ! Magnetic shear (r/q)(dq/dr)
    real(rkind),intent(in):: pne     ! Electron density [m^{-3}]
    real(rkind),intent(in):: dpdr    ! Pressure gradient [Pa/m]
    real(rkind),intent(in):: rhoni   ! Ion mass density [kg/m^3]
    !                                  (sum of ion-mass times ion-density)
    real(rkind),intent(in):: dvexbdr ! ExB drift velocity gradient [1/s]
    real(rkind),intent(in):: calf    ! Factor in s-alpha effects [1.0]
    real(rkind),intent(in):: ckap    ! Factor in magnetic curvature effects [1.0]
    real(rkind),intent(in):: cexb    ! Factor in ExB drift effects [1.0]
    integer(ikind),intent(in):: model! Model ID
    !                                    0: CDBM original
    !                                    1: CDBM05 including elongation
    !                                    2: CDBM original with ExB shear (Lorentzian)
    !                                    3: CDBM05 with ExB shear (Lorentzian)
    !                                    4: CDBM original with strong ExB shear (Exponent)
    !                                    5: CDBM05 with strong ExB shear (Exponent)
    !                                    6: CDBM original with strong ExB shear (Modified Lorentzian)
    !                                    7: CDBM05 with strong ExB shear (Modified Lorentzian)

    real(rkind),intent(out):: chi_cdbm! Thermal diffusion coefficient
    !                                   for both electrons and ions
    real(rkind),intent(out),optional :: fsz   ! Fitting function for output
    real(rkind),intent(out),optional :: curvz ! Magnetic curvature effect for output
    real(rkind),intent(out),optional :: fez   ! ExB shear reduction for output
    real(rkind),intent(in), optional :: omgexb ! ExB shearing rate for model = 6 or 7

    integer(ikind),parameter :: nexp = 2    ! Fixed numerical factor 
    real(rkind),parameter :: ckcdbm = 12.d0 ! Fixed numerical factor 
    real(rkind):: va,wpe2,delta2,alpha,curv,wexb,shearl,fs,fk,fe,gamcdbm,taua

    if(model.lt.0.or.model.gt.7) then
       write(6,*) 'XX cdbm: model: out of range'
       stop
    endif

    !     Alfven wave velocity
    va=SQRT(bb**2/(RMU0*rhoni))

    !     Square of Plasma frequency
    wpe2=pne*AEE*AEE/(AME*EPS0)

    !     Square of collisionless skin-denpth
    delta2=VC**2/wpe2

    !     Normalized pressure gradient (Shafranov shift factor)
    alpha=-2.D0*RMU0*qp**2*rr/bb**2*dpdr

    !     magnetic curvature
    curv=-(rs/rr)*(1.D0-1.D0/(qp*qp))

    !     rotational shear
    shearl=sqrt(shear**2+0.1D0**2)   !
    wexb = -qp*rr/(shearl*va)*dvexbdr
!    write(6,'(5ES15.7)') qp,rr,1.d0/(shearl*va),dvexbdr,wexb

    SELECT CASE(MOD(model,2))
    CASE(0)
       fk=1.d0
    CASE(1)
       fk=(2.D0*SQRT(rkap)/(1.D0+rkap**2))**1.5D0
    END SELECT

    fs=trcofs(shear,calf*alpha,ckap*curv)

    SELECT CASE(MOD((model/2),4))
    CASE(0)
       fe=1.D0
    CASE(1)
       fe=1.D0/(1.D0+cexb*wexb**2)
    CASE(2)
       fe=fexb(wexb,shear,alpha)
    CASE(3)
       taua=qp*rr/va
       gamcdbm = fs*sqrt(abs(alpha))/taua
       if(gamcdbm /= 0.d0) then
          fe=1.d0/(1.d0+(cexb*omgexb/gamcdbm)**nexp)
       else
          fe=1.d0
       end if
    END SELECT

    chi_cdbm=ckcdbm*fs*fk*fe*SQRT(ABS(alpha))**3*delta2*va/(qp*rr)

    IF(PRESENT(fsz))   fsz=fs
    IF(PRESENT(curvz)) curvz=curv
    IF(PRESENT(fez))   fez=fe
    
    RETURN
  END SUBROUTINE cdbm

! *** Form factor in CDBM model ***

  REAL(rkind) FUNCTION trcofs(shear,alpha,curv)

    IMPLICIT NONE
    real(rkind),intent(in):: shear ! Magnetic shear
    real(rkind),intent(in):: alpha ! Normalized pressure gradient
    real(rkind),intent(in):: curv  ! Magnetic curvature
    real(rkind):: fs1,fs2,sa

    IF(alpha.GE.0.D0) THEN
       sa=shear-alpha
       IF(sa.GE.0.D0) THEN
          fs1=(1.D0+9.0D0*SQRT(2.D0)*sa**2.5D0) &
               & /(SQRT(2.D0)*(1.D0+((-2.D0+(3.D0+2.0D0*sa)*sa)*sa)))
       ELSE
          fs1=1.D0/SQRT(2.D0*(1.D0-2.D0*sa)*(1.D0+(-2.D0+3.D0*sa)*sa))
       ENDIF
       IF(curv.GT.0.D0) THEN
          fs2=SQRT(curv)**3/(shear*shear)
       ELSE
          fs2=0.D0
       ENDIF
    ELSE
       sa=alpha-shear
       IF(sa.GE.0.D0) THEN
          fs1=(1.D0+9.0D0*SQRT(2.D0)*sa**2.5D0) &
               & /(SQRT(2.D0)*(1.D0+((-2.D0+(3.D0+2.0D0*sa)*sa)*sa)))
       ELSE
          fs1=1.D0/SQRT(2.D0*(1.D0-2.D0*sa)*(1.D0+(-2.D0+3.D0*sa)*sa))
       ENDIF
       IF(curv.LT.0.D0) THEN
          fs2=SQRT(-curv)**3/(shear*shear)
       ELSE
          fs2=0.D0
       ENDIF
    ENDIF
    trcofs=MAX(fs1,fs2)
    RETURN
  END FUNCTION trcofs

! *** ExB shearing effect for CDBM model ***
!   [M. Honda (2007) dissertation, Chapter 4]  

  REAL(rkind) FUNCTION fexb(wexb,shear,alpha)

    IMPLICIT NONE
    REAL(rkind),intent(in):: wexb  ! omega ExB
    REAL(rkind),intent(in):: shear ! Magnetic shear
    REAL(rkind),intent(in):: alpha ! Normalized pressure gradient
    REAL(8):: alpha1,alpha2,beta,gamma

    IF(ABS(alpha).LT.1.D-3) THEN
       alpha1=1.D-3
    ELSE
       alpha1=ABS(alpha)
    ENDIF
    beta=0.5D0*alpha1**(-0.602D0) &
         &   *(12.7D0+(-16.0D0+14.7D0*shear)*shear) &
         &   /( 1.0D0+( 0.27D0+1.60D0*shear)*shear)

    alpha2=-10.D0/3.D0*alpha+16.D0/3.D0
    IF(shear.LT.0.D0) THEN
       gamma = 1.D0/(1.1D0*SQRT(1.D0-(1.D0+(2.D0+3.D0*shear)*shear)*shear))+0.75D0
    ELSE
       gamma = (1.D0-0.5D0*shear) &
            & /(1.1D0+(-2.D0+(alpha2+4.D0*shear)*shear)*shear)+0.75D0
    ENDIF
    fexb=EXP(-beta*abs(wexb)**gamma)
    RETURN
  END FUNCTION FEXB
END MODULE cdbm_mod
