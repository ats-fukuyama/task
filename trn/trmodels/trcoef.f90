MODULE trcoef

  USE bpsd_kinds

  PRIVATE
  PUBLIC tr_coef

CONTAINS

! ***** calculate transport coefficients and souce terms *****

  SUBROUTINE tr_coef
    USE trcomm, ONLY: mdltr_nc,mdltr_tb,mdltr_prv,dtr_nc,vtr_nc
    USE trsimple, ONLY: tr_simple
    USE trglf23, ONLY: tr_glf23
    USE trcdbm, ONLY: tr_cdbm
    IMPLICIT NONE

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

! -----Pereverzev method -----

  SUBROUTINE Pereverzev_method
    USE trcomm, ONLY: ikind,rkind,nrmax,nsamax,ph0,rg,rt_prev, &
         mdltr_prv,dprv1,dprv2,dtr,vtr,nsa_neq,dtr_tb,vtr_tb, &
         dtr_prv,vtr_prv,cdtrn,cdtru,cdtrt
    IMPLICIT NONE
    REAL(rkind) :: dtr_new,vtr_new,lt,drt,rtave
    INTEGER(ikind) :: nr,nsa

    DO nr=1,nrmax
       DO nsa=1,nsamax
          drt=rt_prev(nsa,nr)-rt_prev(nsa,nr-1)
          rtave=0.5D0*(rt_prev(nsa,nr)+rt_prev(nsa,nr-1))
!          rtave=rt_prev(nsa,nr)
          IF(drt > 0.D0) THEN
             lt = 0.D0
          ELSE
             lt = - drt/(rg(nr)-rg(nr-1))
          END IF
          SELECT CASE(mdltr_prv)
          CASE(1)
             dtr_new = dprv1
          CASE(2)
             dtr_new = dprv2*dtr_tb(3*nsa,3*nsa,nr)
          CASE(3)
             dtr_new = dprv2*dtr_tb(3*nsa,3*nsa,nr)+dprv1
          END SELECT
          vtr_new = dtr_new * lt / rtave

          dtr_prv(3*nsa-2,nr) = cdtrn*dtr_new
          dtr_prv(3*nsa-1,nr) = cdtru*dtr_new
          dtr_prv(3*nsa  ,nr) = cdtrt*dtr_new

          dtr_tb(3*nsa-2,3*nsa-2,nr) = dtr_tb(3*nsa-2,3*nsa-2,nr) &
                                      +dtr_prv(3*nsa-2,nr)
          dtr_tb(3*nsa-1,3*nsa-1,nr) = dtr_tb(3*nsa-1,3*nsa-1,nr) &
                                      +dtr_prv(3*nsa-1,nr) 
          dtr_tb(3*nsa  ,3*nsa  ,nr) = dtr_tb(3*nsa  ,3*nsa  ,nr) &
                                      +dtr_prv(3*nsa  ,nr)

          vtr_prv(3*nsa-2,nr) = cdtrn*vtr_new
          vtr_prv(3*nsa-1,nr) = cdtru*vtr_new
          vtr_prv(3*nsa  ,nr) = cdtrt*vtr_new

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
END MODULE trcoef
