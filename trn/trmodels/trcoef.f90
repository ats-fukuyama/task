MODULE trcoef

  USE bpsd_kinds
  USE trcalv, ONLY: tr_calv_nr_alloc,tr_calc_variables

  PRIVATE
  PUBLIC tr_coef, Pereverzev_check

CONTAINS

! ***** calculate transport coefficients and souce terms *****

  SUBROUTINE tr_coef
    USE trcomm, ONLY: mdltr_nc,mdltr_tb,mdltr_prv,dtr_nc,vtr_nc
    USE trsimple, ONLY: tr_simple
    USE trglf23, ONLY: tr_glf23
    USE trcdbm, ONLY: tr_cdbm
    IMPLICIT NONE

    call tr_calv_nr_alloc
    call tr_calc_variables

    SELECT CASE(mdltr_nc) ! results are in dtr_nc and vtr_nc
    CASE(0)
       dtr_nc = 0.d0
       vtr_nc = 0.d0
    END SELECT

    SELECT CASE(mdltr_tb) ! results are in dtr_tb and vtr_tb
    CASE(0:9)
       CALL tr_simple
    CASE(60:69)
       CALL tr_glf23
    CASE(130:139)
       CALL tr_cdbm
    END SELECT
    
    IF(mdltr_prv /= 0) CALL Pereverzev_method

    RETURN
  END SUBROUTINE tr_coef


! ----------------------- Pereverzev method -----------------------------
!
!   Numerical stabilization method for stiff turbulent transport models*
!    (especially for GLF23 model)
!   Pereverzev check subroutine is for calculating additional numerical
!    factor in element equations of FEM.
!   
!   Switch variable
!         mdltr_prv = 0 : off
!                   > 1 : on
!
! * G.V.Pereverzev, G.Corrigan, Computer Physics Comm. 179 (2008) 579-585
!------------------------------------------------------------------------
  SUBROUTINE Pereverzev_method
    USE trcomm, ONLY: ikind,rkind,nrmax,nsamax,ph0,rg,rt_prev, &
         mdltr_prv,dprv1,dprv2,dtr,vtr,nsa_neq,dtr_tb,vtr_tb,    &
         dtr_prv,vtr_prv,add_prv,cdtrn,cdtru,cdtrt,rn,ru,rt
    IMPLICIT NONE
    REAL(rkind) :: dtr_new,vtr_old,lt,drt,rtave,dvdrp,dvdrm,dtr_elm,vtr_elm
    INTEGER(ikind) :: nr,nsa

    DO nr=1,nrmax

       dvdrp = rg(nr  )
       dvdrm = rg(nr-1)

       DO nsa=1,nsamax
          drt=rt_prev(nsa,nr)-rt_prev(nsa,nr-1)

          rtave=( (2.d0*dvdrm +      dvdrp)*rt_prev(nsa,nr-1)  &
                 +(     dvdrm + 2.d0*dvdrp)*rt_prev(nsa,nr  )) &
                /(3.d0*(dvdrm+dvdrp))

          IF(drt > 0.D0) THEN
             lt = 0.D0
          ELSE
             lt = drt/(rg(nr)-rg(nr-1))
          END IF
          SELECT CASE(mdltr_prv)
          CASE(1)
             dtr_new = dprv1
          CASE(2)
             dtr_new = dprv2*dtr_tb(3*nsa,3*nsa,nr)
          CASE(3)
             dtr_new = dprv2*dtr_tb(3*nsa,3*nsa,nr)+dprv1
          END SELECT

          dtr_prv(3*nsa-2,nr) = cdtrn*dtr_new
          dtr_prv(3*nsa-1,nr) = cdtru*dtr_new
          dtr_prv(3*nsa  ,nr) = cdtrt*dtr_new

          dtr_tb(3*nsa-2,3*nsa-2,nr) = dtr_tb(3*nsa-2,3*nsa-2,nr) &
                                      +dtr_prv(3*nsa-2,nr)
          dtr_tb(3*nsa-1,3*nsa-1,nr) = dtr_tb(3*nsa-1,3*nsa-1,nr) &
                                      +dtr_prv(3*nsa-1,nr) 
          dtr_tb(3*nsa  ,3*nsa  ,nr) = dtr_tb(3*nsa  ,3*nsa  ,nr) &
                                      +dtr_prv(3*nsa  ,nr)

          vtr_old = dtr_new * lt / rtave

          vtr_prv(3*nsa-2,nr) = cdtrn*vtr_old
          vtr_prv(3*nsa-1,nr) = cdtru*vtr_old
          vtr_prv(3*nsa  ,nr) = cdtrt*vtr_old

          vtr_tb(3*nsa-2,3*nsa-2,nr) = vtr_tb(3*nsa-2,3*nsa-2,nr) &
                                      +vtr_prv(3*nsa-2,nr) 
          vtr_tb(3*nsa-1,3*nsa-1,nr) = vtr_tb(3*nsa-1,3*nsa-1,nr) &
                                      +vtr_prv(3*nsa-1,nr) 
          vtr_tb(3*nsa  ,3*nsa  ,nr) = vtr_tb(3*nsa  ,3*nsa  ,nr) &
                                      +vtr_prv(3*nsa  ,nr) 

       END DO
    END DO
    RETURN
  END SUBROUTINE Pereverzev_method

  SUBROUTINE Pereverzev_check
    USE trcomm, ONLY: ikind,rkind,neqmax,nsamax,nrmax,rg,dtr_prv,vtr_prv,&
                      add_prv,rn,ru,rt,dvrho,ar1rho,ar2rho

    REAL(rkind),DIMENSION(3*neqmax,0:nrmax) :: dtr_elm,vtr_elm
    REAL(rkind) :: dvdrm,dvdrp,gm1p,gm1m,gm2p,gm2m
    INTEGER(ikind) :: nr,nsa
    
    DO nr = 1, nrmax
       dvdrp = dvrho(nr)
       dvdrm = dvrho(nr-1)

       gm1p = dvdrp*ar1rho(nr)
       gm1m = dvdrm*ar1rho(nr-1)

       gm2p = dvdrp*ar2rho(nr)
       gm2m = dvdrm*ar2rho(nr-1)

!       gm1p = rg(nr)
!       gm1m = rg(nr-1)
!       gm2p = rg(nr)
!       gm2m = rg(nr-1)       
       
       DO nsa = 1, nsamax
          dtr_elm(3*nsa-2,nr) =                            &
               0.5D0*dtr_prv(3*nsa-2,nr)/(rg(nr)-rg(nr-1)) &
               *(gm2m+gm2p)*(rn(nsa,nr)-rn(nsa,nr-1))
          dtr_elm(3*nsa-1,nr) =                            &
               0.5D0*dtr_prv(3*nsa-1,nr)/(rg(nr)-rg(nr-1)) &
               *(gm2m+gm2p)*(ru(nsa,nr)-ru(nsa,nr-1))
          dtr_elm(3*nsa  ,nr) =                            &
               0.5D0*dtr_prv(3*nsa  ,nr)/(rg(nr)-rg(nr-1)) &
               *(gm2m+gm2p)*(rt(nsa,nr)-rt(nsa,nr-1))
          
          vtr_elm(3*nsa-2,nr) =                            &
            vtr_prv(3*nsa-2,nr)/6.D0                       &
            *((2.D0*gm1m+gm1p)*rn(nsa,nr-1) + (gm1m+2.D0*gm1p)*rn(nsa,nr))
          vtr_elm(3*nsa-1,nr) =                            &
            vtr_prv(3*nsa-1,nr)/6.D0                       &
            *((2.D0*gm1m+gm1p)*ru(nsa,nr-1) + (gm1m+2.D0*gm1p)*ru(nsa,nr))
          vtr_elm(3*nsa  ,nr) =                            &
            vtr_prv(3*nsa  ,nr)/6.D0                       &
            *((2.D0*gm1m+gm1p)*rt(nsa,nr-1) + (gm1m+2.D0*gm1p)*rt(nsa,nr))

          ! numerically additional term in nodal equation
          add_prv(3*nsa-2,nr) = dtr_elm(3*nsa-2,nr) - vtr_elm(3*nsa-2,nr)
          add_prv(3*nsa-1,nr) = dtr_elm(3*nsa-1,nr) - vtr_elm(3*nsa-1,nr)
          add_prv(3*nsa  ,nr) = dtr_elm(3*nsa  ,nr) - vtr_elm(3*nsa  ,nr)
          
       END DO
    END DO

  END SUBROUTINE Pereverzev_check
END MODULE trcoef
