MODULE trmbgb

  PRIVATE
  PUBLIC tr_mbgb

CONTAINS

  SUBROUTINE tr_mbgb

    USE mixed_Bohm_gyro_Bohm, ONLY: mixed_model

    USE trcomm, ONLY: &
         ikind,rkind,ns_nsa,idnsa,pa,nsamax,nrmax,     &
         RR,ra,BB,rm,nrmax,mdltr_tb,amp,rkev,rn,dtr_tb, &
         cdtrn,cdtru,cdtrt
    USE trcalv, ONLY: &
         rn_im,rn_ecl,rt_e,rt_em,rt_im,rt_ecl,qp_m,mshear,wexbp

    IMPLICIT NONE

    ! --- internal variables ---
    REAL(rkind) :: amm,zte_p8,zte_edge,pa_ave,sum_pan
    REAL(rkind) :: vti,gamma0,EXBfactor,SHRfactor,factor
    REAL(rkind),DIMENSION(1):: &
            rmajor,rminor,btor,tekev,tikev,q,aimass,charge,wexbs, &
            grdte,grdne,shear, &
            chi_i_mix,themix,thdmix,thigb,thegb,thibohm,thebohm
!    REAL(rkind),DIMENSION(0:nrmax) :: D_hyd,chie_b,chii_b,chie_gb,chii_gb
    INTEGER(ikind):: nr8,npoints,lflowshear

    INTEGER(ikind):: nr,nsa,ns,ierr

    amm      = amp          ! proton mass
    npoints  = 1            ! num. of values of 'jz' in all of arrays [int]
                            !  * 'jz' corresponds to 'nr' in TASK.

    nr8      = NINT(nrmax*0.8d0)
    zte_p8   = rt_e(nr8)    ! T_e(0.8a)
    zte_edge = rt_e(nrmax)  ! T_e(a)

    DO nr = 1, nrmax
       ! minor radius (half-width) of zone boundary [m]
       rminor(1)  = ra*rm(nr)
       ! major radius to geometric center of zone boundary [m] <- approx.
       rmajor(1)   = RR
       btor(1)     = BB          ! toroidal magnetic field [T] <- approx.

       sum_pan = 0.d0
       DO nsa = 1, nsamax
          IF(idnsa(nsa)==1)THEN
             ns = ns_nsa(nsa)
             sum_pan = sum_pan+pa(ns)*0.5d0*(rn(nsa,nr-1)+rn(nsa,nr))
          END IF
       END DO
       pa_ave    = sum_pan / rn_im(nr) ! average ion mass [AMU]
       aimass(1) = pa_ave

       charge(1) = 1.0               ! charge number of main thermal ions
                                     ! ( 1.0 for hydrogenic ions )

       tekev(1)  = rt_em(nr)         ! T_e [keV]
       tikev(1)  = rt_im(nr)         ! T_i [keV] :effective ion temperature
       q(1)      = qp_m(nr)          ! safety-factor (half-mesh)

       grdte(1)  = - RR / rt_ecl(nr) ! -R ( d T_e / d r ) / T_e
       grdne(1)  = - RR / rn_ecl(nr) ! -R ( d n_e / d r ) / n_e
       shear(1)  = mshear(nr)        !  r ( d q   / d r ) / q
       
       wexbs(1)  = wexbp(nr)         ! ExB shearing rate [rad/s]


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
          vti        = SQRT(2.d0*rt_im(nr)*rkev / (pa_ave*amm))
          gamma0     = vti / (qp_m(nr)*RR)
          SHRfactor  = 1.d0
          EXBfactor  = 1.d0 / (1.d0+(wexbp(nr)/gamma0)**2)
       CASE(143)
          lflowshear = 0
          vti        = SQRT(2.d0*rt_im(nr)*rkev/(pa_ave*amm))
          gamma0     = vti / (qp_m(nr)*RR)
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
!   thdmix(jz)    = total hydrogenic ion diffusivity from the MIXED model
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

       factor   = EXBfactor * SHRfactor

       DO nsa = 1, nsamax
          IF(idnsa(nsa) == -1)THEN ! electron
             dtr_tb(3*nsa-2,3*nsa-2,nr) = cdtrn *themix(1)*factor
             dtr_tb(3*nsa-1,3*nsa-1,nr) = cdtru *themix(1)*factor
             dtr_tb(3*nsa  ,3*nsa  ,nr) = cdtrt *themix(1)*factor
          ELSE IF(idnsa(nsa) /= 0)THEN ! ion
             dtr_tb(3*nsa-2,3*nsa-2,nr) = cdtrn *chi_i_mix(1)*factor
             dtr_tb(3*nsa-1,3*nsa-1,nr) = cdtru *chi_i_mix(1)*factor
             dtr_tb(3*nsa  ,3*nsa  ,nr) = cdtrt *chi_i_mix(1)*factor
             
          END IF
       END DO

       ! *** for diagnostic ***
       ! Hydrogenic ion particle diffusivity [m^2/s]
!       D_hyd(nr)   = thdmix(1)  *factor 
       ! Bohm contribution to electron thermal diffusivity [m^2/s]
!       chie_b(nr)  = thebohm(1) *factor
       ! Bohm contribution to ion thermal diffusivity [m^2/s]
!       chii_b(nr)  = thibohm(1) *factor
       ! gyro-Bohm contribution to electron thermal diffusivity [m^2/s]
!       chie_gb(nr) = thegb(1)   *factor
       ! gyro-Bohm contribution to ion thermal diffusivity [m^2/s]
!       chii_gb(nr) = thigb(1)   *factor

    END DO

    RETURN
  END SUBROUTINE tr_mbgb

END MODULE trmbgb
