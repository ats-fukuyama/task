!     ********************************************

!           CDBM Transport model (2009/03/06)
!              Modified by M. Honda (2009/09/14)

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
  real(rkind),parameter :: VC   = 2.99792458E8_dp    ! speed of light
  real(rkind),parameter :: RMU0 = 4.E-7_dp*PI        ! permeability
  real(rkind),parameter :: EPS0 = ONE/(VC*VC*RMU0)   ! permittivity

CONTAINS

  SUBROUTINE cdbm(bb,rr,rs,rkap,qp,shear,pne,rhoni,dpdr,dvexbdr, &
       &      calf,ckap,cexb,model,chi_cdbm,                     &
       &      fsz,curvz,fez,mdl_tuned,c1_tuned,c2_tuned,k_tuned)

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
    !                               0: CDBM original
    !                               1: CDBM05 including elongation
    !                               2: CDBM original with weak ExB shear
    !                               3: CDBM05 with weak ExB shear
    !                               4: CDBM original with strong ExB shear
    !                               5: CDBM05 with strong ExB shear

    real(rkind),intent(out):: chi_cdbm! Thermal diffusion coefficient
    !                                   for both electrons and ions

    real(rkind),intent(out),optional :: fsz   ! Fitting function for output
    real(rkind),intent(out),optional :: curvz ! Magnetic curvature effect for output
    real(rkind),intent(out),optional :: fez   ! ExB shear reduction for output

    ! +++ implemented by Ikari 12/02/16 ++++++++++++++++++++++++++++++++++
    ! for Tuned CDBM
    integer(ikind),intent(in),optional :: mdl_tuned   ! Tuned model ID
    !                                      0: conventional CDBM/CDBM05
    !                                      1: Tuned CDBM by M.Yagi
    real(rkind),intent(in),optional :: c1_tuned, c2_tuned, k_tuned
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    real(rkind),parameter :: ckcdbm = 12.d0 ! Fixed numerical factor

    real(rkind):: va,wpe2,delta2,alpha,curv,wexb,shearl,fs,fk,fe,gamma_tuned

    if(model.lt.0.or.model.gt.5) then
       write(6,*) 'XX cdbm: model: out of range'
       stop
    endif

    if(present(mdl_tuned))then
       if(mdl_tuned.lt.0 .or. mdl_tuned.gt.1)then
          write(6,*) 'XX cdbm: mdl_tuned: out of range'
          stop
       endif
       if(.not. &
          (present(c1_tuned).and.present(c2_tuned).and.present(k_tuned)))then
          write(6,*) 'XX cdbm: factors for tuned CDBM are not found.'
          stop
       endif
    endif

    !     Alfven wave velocity
    va     = SQRT(bb**2/(RMU0*rhoni))

    !     Square of Plasma frequency
    wpe2   = pne*AEE*AEE/(AME*EPS0)

    !     Square of collisionless skin-denpth
    delta2 = VC**2/wpe2

    !     Normalized pressure gradient (Shafranov shift factor)
    alpha  = -2.D0*RMU0*qp**2*rr/bb**2*dpdr

    !     magnetic curvature
    curv   = -(rs/rr)*(1.D0-1.D0/(qp*qp))

    !     rotational shear
    shearl = sqrt(shear**2+0.1D0**2)
    wexb   = -qp*rr/(shearl*va)*dvexbdr

    SELECT CASE(MOD(model,2))
    CASE(0) ! CDBM original
       fk = 1.d0
    CASE(1) ! CDBM05
       fk = (2.D0*SQRT(rkap)/(1.D0+rkap**2))**1.5D0
    END SELECT

    SELECT CASE(MOD((model/2),3))
    CASE(0) ! without ExB shear ( model: 0, 1 )
       fe = 1.D0
    CASE(1) ! with weak ExB shear ( model: 2, 3 )
       fe = 1.D0/(1.D0+cexb*wexb**2)
    CASE(2) ! with strong ExB shear ( model: 4, 5 )
       shearl = sqrt(shear**2+0.1D0**2)
       fe     = cexb*fexb(wexb,shear,alpha)
    END SELECT


    IF(PRESENT(mdl_tuned).AND.mdl_tuned==1)THEN ! Tuned CDBM
       fs          = trcofs_tuned(shear,calf*alpha,ckap*curv)
       gamma_tuned = va/(qp*RR)*ABS(alpha)**0.5d0*fs

       chi_cdbm = ckcdbm*fs*fk*fe*SQRT(ABS(alpha))**3*delta2*va/(qp*rr) &
                  *c1_tuned/(1.d0+(c2_tuned*wexb/gamma_tuned)**k_tuned)
    ELSE ! original CDBM / CDBM05
       fs       = trcofs(shear,calf*alpha,ckap*curv)
       chi_cdbm = ckcdbm*fs*fk*fe*SQRT(ABS(alpha))**3*delta2*va/(qp*rr)
    END IF

    IF(PRESENT(fsz))   fsz   = fs
    IF(PRESENT(curvz)) curvz = curv
    IF(PRESENT(fez))   fez   = fe
    
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
               & /(SQRT(2.D0)*(1.D0-2.D0*sa+3.D0*sa*sa+2.0D0*sa*sa*sa))
       ELSE
          fs1=1.D0/SQRT(2.D0*(1.D0-2.D0*sa)*(1.D0-2.D0*sa+3.D0*sa*sa))
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
               & /(SQRT(2.D0)*(1.D0-2.D0*sa+3.D0*sa*sa+2.0D0*sa*sa*sa))
       ELSE
          fs1=1.D0/SQRT(2.D0*(1.D0-2.D0*sa)*(1.D0-2.D0*sa+3.D0*sa*sa))
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


  ! for Tuned CDBM model by M.Yagi (implemented by Ikari 12/02/16)
  REAL(rkind) FUNCTION trcofs_tuned(shear,alpha,curv)

    IMPLICIT NONE
    real(rkind),intent(in):: shear ! Magnetic shear
    real(rkind),intent(in):: alpha ! Normalized pressure gradient
    real(rkind),intent(in):: curv  ! Magnetic curvature
    real(rkind):: fs1,fs2,sa


    IF(alpha.GE.0.D0) THEN
       sa=shear-alpha
       IF(sa.GE.0.D0) THEN
          fs1=SQRT(3.d0)*(1.D0+2.D0*sa**2.5D0) &
               & /(1.D0-sa+1.35d0*sa*sa+0.7698d0*sa*sa*sa)
       ELSE
          fs1=SQRT(3.d0/((1.D0-sa)*(1.D0-sa+3.d0*sa*sa)))
       ENDIF
       IF(curv.GT.0.D0) THEN
          fs2=SQRT(curv)**3/(shear*shear)
       ELSE
          fs2=0.D0
       ENDIF
    ELSE
       sa=alpha-shear
       IF(sa.GE.0.D0) THEN
          fs1=SQRT(3.d0)*(1.D0+2.D0*sa**2.5D0) &
               & /(1.D0-sa+1.35d0*sa*sa+0.7698d0*sa*sa*sa)
       ELSE
          fs1=SQRT(3.d0/((1.D0-sa)*(1.D0-sa+3.d0*sa*sa)))
       ENDIF
       IF(curv.LT.0.D0) THEN
          fs2=SQRT(-curv)**3/(shear*shear)
       ELSE
          fs2=0.D0
       ENDIF
    ENDIF

    trcofs_tuned=MAX(fs1,fs2)
    RETURN
  END FUNCTION trcofs_tuned

! *** ExB shearing effect for CDBM model ***

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
         &   *(13.018D0-22.28915D0*shear+17.018D0*shear**2) &
         &   /(1.D0-0.277584D0*shear+1.42913D0*shear**2)

    alpha2=-10.D0/3.D0*alpha+16.D0/3.D0
    IF(shear.LT.0.D0) THEN
       gamma = 1.D0/(1.1D0*SQRT(1.D0-shear-2.D0*shear**2-3.D0*shear**3)) &
            & +0.75D0
    ELSE
       gamma = (1.D0-0.5D0*shear) &
            & /(1.1D0-2.D0*shear+alpha2*shear**2+4.D0*shear**3)+0.75D0
    ENDIF
    fexb=EXP(-beta*abs(wexb)**gamma)
    RETURN
  END FUNCTION FEXB
END MODULE cdbm_mod
