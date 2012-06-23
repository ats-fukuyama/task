MODULE trmbgb

  PRIVATE
  PUBLIC tr_mbgb

CONTAINS

! *************************************************************************
!
!  The interface between TASK/TR(trcoeftb) and mixed Bohm/gyro-Bohm model
!
! *************************************************************************

  SUBROUTINE tr_mbgb

    USE mixed_Bohm_gyro_Bohm, ONLY: mixed_model

    USE trcomm, ONLY: &
         ikind,rkind,ns_nsa,idnsa,pa,nsamax,nrmax,neqmax,RR,ra,BB, &
         abb1rho,rmnrho,rmjrho,nrmax,mdltr_tb,amp,rkev,rn,qp,      &
         dtr_tb,vtr_tb,cdtrn,cdtru,cdtrt !  ,nrd1,nrd2
    USE trcalv, ONLY: &
         rn_i,rn_ecl,rt_e,rt_i,rt_ecl,mshear,wexbp

    IMPLICIT NONE

    ! --- internal variables ---
    REAL(rkind) :: amm,zte_p8,zte_edge,pa_ave,sum_pan
    REAL(rkind) :: vti,gamma0,EXBfactor,SHRfactor,factor
    REAL(rkind),DIMENSION(1):: &
            rmajor,rminor,btor,tekev,tikev,q,aimass,charge,wexbs, &
            grdte,grdne,shear, &
            chi_i_mix,themix,thdmix,thigb,thegb,thibohm,thebohm

    REAL(rkind),DIMENSION(0:nrmax) :: &
         mbgb_chie,mbgb_chii
    REAL(rkind),DIMENSION(1:nrmax) :: &
         mbgb_chiem,mbgb_chiim

!    REAL(rkind),DIMENSION(0:nrmax) :: D_hyd,chie_b,chii_b,chie_gb,chii_gb
    REAL(rkind) :: FCTR
    INTEGER(ikind):: nr8,npoints,lflowshear

    INTEGER(ikind):: nr,nsa,ns,ierr

    dtr_tb(1:neqmax,1:neqmax,0:nrmax) = 0.D0
    vtr_tb(1:neqmax,1:neqmax,0:nrmax) = 0.D0
    mbgb_chie = 0.d0
    mbgb_chii = 0.d0

    amm      = amp          ! proton mass
    npoints  = 1            ! num. of values of 'jz' in all of arrays [int]
                            !  * 'jz' corresponds to 'nr' in TASK.

    nr8      = NINT(nrmax*0.8d0)
    zte_p8   = rt_e(nr8)    ! T_e(0.8a)
    zte_edge = rt_e(nrmax)  ! T_e(a)

    DO nr = 1, nrmax
       ! minor radius (half-width) of zone boundary [m]
       rminor(1)  = rmnrho(nr)
       ! major radius to geometric center of zone boundary [m] <- approx.
       rmajor(1)  = rmjrho(nr)
!       btor(1)    = BB     ! toroidal magnetic field [T] <- approx.
       btor(1) = abb1rho(nr) ! toroidal magnetic field [T]

       sum_pan = 0.d0
       DO nsa = 1, nsamax
          IF(idnsa(nsa)==1)THEN
             ns = ns_nsa(nsa)
             sum_pan = sum_pan+pa(ns)*rn(nsa,nr)
          END IF
       END DO
       pa_ave    = sum_pan / rn_i(nr) ! average ion mass [AMU]
       aimass(1) = pa_ave

       charge(1) = 1.0               ! charge number of main thermal ions
                                     ! ( 1.0 for hydrogenic ions )

       tekev(1)  = rt_e(nr)         ! T_e [keV]
       tikev(1)  = rt_i(nr)         ! T_i [keV] :effective ion temperature
       q(1)      = qp(nr)          ! safety-factor (half-mesh)

       ! *** There is a fault in 'mixed model' subroutine ***
       !      fault   : -RR ( d T_e / d r ) / T_e
       !                -RR ( d n_e / d r ) / n_e
       !      correct : -Ra ( d T_e / d r ) / T_e
       !                -Ra ( d n_e / d r ) / n_e

       grdte(1)  = - ra * rt_ecl(nr) ! -Ra ( d T_e / d r ) / T_e
       grdne(1)  = - ra * rn_ecl(nr) ! -Ra ( d n_e / d r ) / n_e
       shear(1)  = mshear(nr)        !  r ( d q   / d r ) / q
       
       wexbs(1)  = wexbp(nr)         ! ExB shearing rate [rad/s]
!      wexbs : ExB Rotation shear
!       "Effects of {ExB} velocity shear and magnetic shear
!           on turbulence and transport in magnetic confinement devices"
!       [Phys. of Plasmas, 4, 1499 (1997)]


       SELECT CASE(mdltr_tb)
       CASE(140)
          lflowshear = 0
          SHRfactor  = 1.d0
          EXBfactor  = 1.d0
       CASE(141)
          lflowshear = 1
          SHRfactor  = 1.d0
          EXBfactor  = 1.d0
       CASE(142)
          lflowshear = 0
          vti        = SQRT(2.d0*rt_i(nr)*rkev / (pa_ave*amm))
          gamma0     = vti / (qp(nr)*RR)
          SHRfactor  = 1.d0
          EXBfactor  = 1.d0 / (1.d0+(wexbp(nr)/gamma0)**2)
       CASE(143)
          lflowshear = 0
          vti        = SQRT(2.d0*rt_i(nr)*rkev/(pa_ave*amm))
          gamma0     = vti / (qp(nr)*RR)
          SHRfactor  = 1.d0 / MAX(1.d0,(mshear(nr)-0.5d0)**2)
          EXBfactor  = 1.d0 / (1.d0+(wexbp(nr)/gamma0)**2)
       CASE DEFAULT
          lflowshear = 0
          SHRfactor  = 1.d0
          EXBfactor  = 1.d0
       END SELECT

       call mixed_model ( &
            rminor,  rmajor,  tekev,   tikev,   q, &
            btor,    aimass,  charge,  wexbs, &
            grdte,   grdne,   shear, &
            zte_p8, zte_edge, npoints, &
            chi_i_mix,  themix,   thdmix, &
            thigb,   thegb,    thibohm, thebohm, &
            ierr, lflowshear)

       IF(ierr /= 0)THEN
          WRITE(*,*) 'error in subroutine "mixed_model"(mbgb): ierr =',ierr
          STOP
       END IF

! ***  lflowshear = 0 : for no magnetic and flow shear stabilization
! ***             = 1 : to use magnetic and flow shear stabilization
! ***  [ T.J.Tala et al. Plasma Phys. Controlled Fusion 44 (2002) A495]

!  << Output >>
!
! The following effective diffusivities are given in MKS units m^2/sec
!
!   chi_i_mix(jz) = total ion thermal diffusivity from the MIXED model
!   themix(jz)    = total electron thermal diffusivity from the MIXED model
!   thdmix(jz)    = total hydrogenic ion particle diffusivity 
!                                                      from the MIXED model
!
! The following contributions to the effective diffusivities are
!  for diagnostic purposes:
!
!   thigb(jz) = gyro-Bohm contribution to the ion thermal diffusivity
!   thegb(jz) = gyro-Bohm contribution to the electron thermal diffusivity
!
!   thibohm(jz) = Bohm contribution to the ion thermal diffusivity
!   thebohm(jz) = Bohm contribution to the electron thermal diffusivity
!
!   ierr    = returning with value .ne. 0 indicates error                  

       ! *** for diagnostic ***
!!$       ! Hydrogenic ion particle diffusivity [m^2/s]
!!$       D_hyd(nr)   = thdmix(1)  *factor 
!!$       ! Bohm contribution to electron thermal diffusivity [m^2/s]
!!$       chie_b(nr)  = thebohm(1) *factor
!!$       ! Bohm contribution to ion thermal diffusivity [m^2/s]
!!$       chii_b(nr)  = thibohm(1) *factor
!!$       ! gyro-Bohm contribution to electron thermal diffusivity [m^2/s]
!!$       chie_gb(nr) = thegb(1)   *factor
!!$       ! gyro-Bohm contribution to ion thermal diffusivity [m^2/s]
!!$       chii_gb(nr) = thigb(1)   *factor

       mbgb_chie(nr) = themix(1)
       mbgb_chii(nr) = chi_i_mix(1)
    END DO

    mbgb_chie(0) = &
         FCTR(rmnrho(1),rmnrho(2),mbgb_chie(1),mbgb_chie(2))
    mbgb_chii(0) = &
         FCTR(rmnrho(1),rmnrho(2),mbgb_chii(1),mbgb_chii(2))
    ! on grid -> on half grid
    mbgb_chiem(1:nrmax) = 0.5d0*(mbgb_chie(0:nrmax-1)+mbgb_chie(1:nrmax))
    mbgb_chiim(1:nrmax) = 0.5d0*(mbgb_chii(0:nrmax-1)+mbgb_chii(1:nrmax))

    factor   = EXBfactor * SHRfactor

    DO nr = 1, nrmax
       DO nsa = 1, nsamax
          IF(idnsa(nsa) == -1)THEN ! electron
             dtr_tb(1+3*nsa-2,1+3*nsa-2,nr) = cdtrn *mbgb_chiem(nr)*factor
             dtr_tb(1+3*nsa-1,1+3*nsa-1,nr) = cdtru *mbgb_chiem(nr)*factor
             dtr_tb(1+3*nsa  ,1+3*nsa  ,nr) = cdtrt *mbgb_chiem(nr)*factor
          ELSE IF(idnsa(nsa) /= 0)THEN ! ion
             dtr_tb(1+3*nsa-2,1+3*nsa-2,nr) = cdtrn *mbgb_chiim(nr)*factor
             dtr_tb(1+3*nsa-1,1+3*nsa-1,nr) = cdtru *mbgb_chiim(nr)*factor
             dtr_tb(1+3*nsa  ,1+3*nsa  ,nr) = cdtrt *mbgb_chiim(nr)*factor
          END IF
       END DO
    END DO

!    nrd1(1:nrmax) = chie_b(1:nrmax)
!    nrd2(1:nrmax) = chie_gb(1:nrmax)

    RETURN
  END SUBROUTINE tr_mbgb

END MODULE trmbgb
