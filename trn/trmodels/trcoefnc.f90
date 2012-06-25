MODULE trcoefnc

  USE trcomm, ONLY: ikind, rkind

  PRIVATE
  PUBLIC tr_coefnc,tr_calc_clseta,coulog

CONTAINS

  SUBROUTINE tr_coefnc
    USE trcomm, ONLY: &
         rkev,nsa_neq,nva_neq,neqmax,nsamax,nrmax,mdltr_nc,  &
         rmnrho,htr,eta,dtr_nc,vtr_nc,eta_nc,jbs_nc,jex_nc,  &
         rt,rn !,nrd1,nrd2
    USE trcalv, ONLY: &
         chi_ncp,chi_nct,d_ncp,d_nct,gfls,qfls,fls_tot, &
         vebs,qebs,dia_gdnc,dia_gvnc,cjbs_p,cjbs_t
    USE trncls, ONLY: tr_nclass

    IMPLICIT NONE
    REAL(rkind) :: deriv3

    ! internal parameters
    INTEGER(ikind) :: nr,nsa,nva,neq,nk,ierr
    INTEGER(ikind),SAVE :: etacls_save = 0
    REAL(rkind),DIMENSION(0:nrmax) :: rt_s,d_nc,v_nc

    dtr_nc(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    vtr_nc(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    jbs_nc(0:nrmax) = 0.d0

    ! CALL tr_calnc
    ! trapped particle fraction
    

    SELECT CASE(mdltr_nc) ! results are in dtr_nc and vtr_nc
    CASE(0)
       ! no transport
       
       ! resistivity
       CALL tr_calc_clseta
    CASE(1)
       ! preparation of eta for NCLASS calculation
       IF(etacls_save == 0)THEN
          CALL tr_calc_clseta
       END IF
       etacls_save = 1

       CALL tr_nclass(ierr)

       ! diffusion and convection coefficients
       ! ** diagonal term only
       fls_tot(3,1:nsamax,0:nrmax) = 0.d0
       DO neq = 1, neqmax
          nsa = nsa_neq(neq)
          nva = nva_neq(neq)
         
          IF(nva == 1)THEN! particle
             DO nk = 1, 5
                fls_tot(nva,nsa,0:nrmax) = fls_tot(nva,nsa,0:nrmax) &
                                            + gfls(nk,nsa,0:nrmax)
             END DO
             ! on half grid
             dtr_nc(neq,neq,1:nrmax) = &
                  0.5d0*(dia_gdnc(nsa,0:nrmax-1)+dia_gdnc(nsa,1:nrmax))
             vtr_nc(neq,neq,1:nrmax) = &
                  0.5d0*(dia_gvnc(nsa,0:nrmax-1)+dia_gvnc(nsa,1:nrmax))

          ELSE IF(nva == 2)THEN! velocity
!             dtr_nc(neq,neq,0:nrmax) = dia_dnc(nsa,0:nrmax)
!             vtr_nc(neq,neq,0:nrmax) = dia_vnc(nsa,0:nrmax)

          ELSE IF(nva == 3)THEN! energy
             DO nk = 1, 5
                fls_tot(nva,nsa,0:nrmax) = fls_tot(nva,nsa,0:nrmax) &
                                            + qfls(nk,nsa,0:nrmax)
             END DO             

             rt_s(0:nrmax) = rt(nsa,0:nrmax)
             DO nr = 0, nrmax
                d_nc(nr) = chi_nct(nsa,nsa,nr)+chi_ncp(nsa,nsa,nr)
                
                v_nc(nr) &
                  = fls_tot(nva,nsa,nr)                                    &
                    /(rn(nsa,nr)*rt(nsa,nr)*rkev)                          &
                  + d_nc(nr)*deriv3(nr,rmnrho,rt_s,nrmax,0) &
                    /rt(nsa,nr)
             END DO
             dtr_nc(neq,neq,1:nrmax) = &
                  0.5d0*(d_nc(0:nrmax-1)+d_nc(1:nrmax))
             ! interim way of substitution V_Es = V_Ks + (3/2)V_s
             vtr_nc(neq,neq,1:nrmax) =                  &
                  0.5d0*(v_nc(0:nrmax-1)+v_nc(1:nrmax)) &
                 +1.5d0*vtr_nc(neq-2,neq-2,1:nrmax)
          END IF
       END DO

!       nrd1(0:nrmax) = dia_gdnc(1,0:nrmax)
!       nrd2(0:nrmax) = dia_gvnc(1,0:nrmax)

       ! resistivity
       eta(0:nrmax) = eta_nc(0:nrmax)

       ! bootstrap current
!       htr(1,0:nrmax) = jbs_nc(0:nrmax) + jex_nc(0:nrmax)   

       ! *** off diagonal term...       

    END SELECT
    
    RETURN
  END SUBROUTINE tr_coefnc

!!$  SUBROUTINE tr_calc_bootstrap
!!$
!!$    RETURN
!!$  END SUBROUTINE tr_calc_bootstrap
!!$
!!$  SUBROUTINE tr_calc_nceta
!!$
!!$    RETURN
!!$  END SUBROUTINE tr_calc_nceta


! from trcalc.f90 in TASK/TR previous version
! **********************************************************************
  SUBROUTINE tr_calc_clseta
    USE trcomm, ONLY: aee,ame,nrmax,rn,rt,eta
    USE trcalv, ONLY: z_eff

    IMPLICIT NONE
    INTEGER(ikind) :: nr
    REAL(rkind) :: taue

      DO nr = 0, nrmax
!        ****** CLASSICAL RESISTIVITY (Spitzer) from JAERI Report ******

         ! electron collision time with ions
         taue = FTAUE(rn(1,nr),rn(2,nr),rt(1,nr),z_eff(nr))

         eta(nr) = ame/(rn(1,nr)*1.d20*aee**2*taue) &
                      *(0.29d0+0.46d0/(1.08d0+z_eff(nr)))

!         eta(nr) = 1.d-7
      END DO

    RETURN
  END SUBROUTINE tr_calc_clseta

! ----------------------------------------------------------------------
!           COULOMB LOGARITHM
! ----------------------------------------------------------------------
      REAL(8) FUNCTION COULOG(NS1,NS2,ANEL,TL)

!     ANEL : electron density [10^20 /m^3]
!     TL   : electron or ion temperature [keV]
!            in case of ion-ion collision, TL becomes ion temp.

      IMPLICIT NONE
      INTEGER(4):: NS1,NS2
      REAL(8)   :: ANEL,TL

      ! *** ref: 'Tokamaks 3rd Edition' p727 Coulomb logarithm ***
      IF(NS1.EQ.1.AND.NS2.EQ.1) THEN
         COULOG=14.9D0-0.5D0*LOG(ANEL)+LOG(TL)
      ELSE
         IF(NS1.EQ.1.OR.NS2.EQ.1) THEN
            COULOG=15.2D0-0.5D0*LOG(ANEL)+LOG(TL)
         ELSE
            COULOG=17.3D0-0.5D0*LOG(ANEL)+1.5D0*LOG(TL)
         ENDIF
      ENDIF

      RETURN
      END FUNCTION COULOG

! ----------------------------------------------------------------------
!           COLLISION TIME
! ----------------------------------------------------------------------
!     between electrons and ions
      REAL(8) FUNCTION FTAUE(ANEL,ANIL,TEL,ZL)

!     ANEL : electron density [10^20 /m^3]
!     ANIL : ion density [10^20 /m^3]
!     TEL  : electron temperature [kev]
!     ZL   : ion charge number

      USE TRCOMM, ONLY : AEE, AME, EPS0, PI, PZ, RKEV
      IMPLICIT NONE
      REAL(8) :: ANEL, ANIL, TEL, ZL
      REAL(8) :: COEF

      COEF = 6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)/(AEE**4*1.D20)
      IF(ZL-PZ(2).LE.1.D-7) THEN
         FTAUE = COEF*(TEL*RKEV)**1.5D0/(ANIL*ZL**2*COULOG(1,2,ANEL,TEL))
      ELSE
!     If the plasma contains impurities, we need to consider the
!     effective charge number instead of ion charge number.
!     From the definition of Zeff=sum(n_iZ_i^2)/n_e,
!     n_iZ_i^2 is replaced by n_eZ_eff at the denominator of tau_e.
         FTAUE = COEF*(TEL*RKEV)**1.5D0/(ANEL*ZL*COULOG(1,2,ANEL,TEL))
      ENDIF

      RETURN
      END FUNCTION FTAUE

!     between ions and ions
      REAL(8) FUNCTION FTAUI(ANEL,ANIL,TIL,ZL,PAL)

!     ANEL : electron density [10^20 /m^3]
!     ANIL : ion density [10^20 /m^3]
!     TIL  : ion temperature [kev]
!     ZL   : ion charge number
!     PAL  : ion atomic number

      USE TRCOMM, ONLY : AEE, AMP, EPS0, PI, RKEV
      IMPLICIT NONE
      REAL(8):: ANEL, ANIL, PAL, TIL, ZL
      REAL(8):: COEF

      COEF = 12.D0*PI*SQRT(PI)*EPS0**2*SQRT(PAL*AMP)/(AEE**4*1.D20)
      FTAUI = COEF*(TIL*RKEV)**1.5D0/(ANIL*ZL**4*COULOG(2,2,ANEL,TIL))

      RETURN
      END FUNCTION FTAUI

END MODULE trcoefnc
