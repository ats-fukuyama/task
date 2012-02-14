!  ************************************************************
!
!   The interface between TASK/TR(trcalv) and CDBM model(cdbm)
!
!  ************************************************************

MODULE trcdbm

  USE bpsd_kinds

  PRIVATE
  PUBLIC:: tr_cdbm

CONTAINS

  SUBROUTINE tr_cdbm

    USE cdbm_mod, ONLY: cdbm
    USE trcomm, ONLY: &
         bb,rr,rg,rm,rkap,rkev,qp,rn,rt,pa,amp,aee,cdtrn,cdtru,cdtrt, &
         dtr_tb,vtr_tb,nrmax,nsamax,ns_nsa,idnsa,mdltr_tb

    USE trcalv, ONLY: &
         rn_e,    &! the density of electron
         rp_totd, &! the total pressure
         qp_m,    &! safety factor (half-mesh)
         mshear,  &! magnetic shear (half-mesh)
         dvexbpdr  ! the gradient of ExBp drift velocity

    IMPLICIT NONE
    INTEGER(ikind):: nr,ns,nsa,model
    REAL(rkind):: calf,ckap,cexb, pnel,rhoni,chi_cdbm

    dtr_tb(1:3*nsamax,1:3*nsamax,0:nrmax) = 0.D0
    vtr_tb(1:3*nsamax,1:3*nsamax,0:nrmax) = 0.D0

! model : Model ID
!     0 : CDBM original
!     1 : CDBM05 including elongation
!     2 : CDBM original with weak ExB shear
!     3 : CDBM05 with weak ExB shear
!     4 : CDBM original with strong ExB shear
!     5 : CDBM05 with strong ExB shear       
    model = mdltr_tb - 130

    ! Factor in s-slpha effects [1.0]
    calf=1.D0
    ! Factor in magnetic curvature effects [1.0]
    ckap=1.D0
    ! Factor in ExB drift effects [1.0]
    cexb=1.D0

    DO nr = 1, nrmax
       ! electron density [m^-3]
       pnel = rn_e(nr)*1.d20
    
       rhoni = 0.D0
       DO nsa = 1, nsamax
          IF(idnsa(nsa)==1) THEN
             ns = ns_nsa(nsa)
             ! ion mass density [kg/m^3]
             ! (sum of ion-mass times ion-density)
             rhoni = rhoni + pa(ns)*amp*rn(nsa,nr)*1.D20
          END IF
       END DO

       CALL cdbm(BB,RR,rm(nr),rkap,qp_m(nr),mshear(nr),pnel,rhoni, &
            rp_totd(nr),dvexbpdr(nr),calf,ckap,cexb,model,chi_cdbm)

       DO nsa = 1, nsamax
          IF(idnsa(nsa) /= 0) THEN
             dtr_tb(3*nsa-2,3*nsa-2,nr) = cdtrn*chi_cdbm
             dtr_tb(3*nsa-1,3*nsa-1,nr) = cdtru*chi_cdbm
             dtr_tb(3*nsa  ,3*nsa  ,nr) = cdtrt*chi_cdbm
          ENDIF
       END DO
    END DO
    RETURN

  END SUBROUTINE tr_cdbm
END MODULE trcdbm
