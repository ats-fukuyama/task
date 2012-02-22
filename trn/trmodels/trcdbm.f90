MODULE trcdbm

  PRIVATE
  PUBLIC:: tr_cdbm

CONTAINS

!  ************************************************************
!
!   The interface between TASK/TR(trcalv) and CDBM model(cdbm)
!
!  ************************************************************

  SUBROUTINE tr_cdbm

    USE cdbm_mod, ONLY: cdbm
    USE trcomm, ONLY: &
         ikind,rkind,bb,rr,rmnrho,rkev,qp,rn,pa,amp,aee, &
         rkprho,cdtrn,cdtru,cdtrt, &
         dtr_tb,vtr_tb,nrmax,nsamax,ns_nsa,idnsa,mdltr_tb

    USE trcalv, ONLY: &
         rn_e,    &! the density of electron
         rp_totd, &! the total pressure
         qp_m,    &! safety factor (half-mesh)
         mshear,  &! magnetic shear (half-mesh)
         dvexbpdr  ! the gradient of ExBp drift velocity

    USE trlib, ONLY: mesh_convert_mtog,data_interpolate_gtom

    IMPLICIT NONE
    INTEGER(ikind):: nr,ns,nsa,model,model2
    REAL(rkind):: calf,ckap,cexb,qp_eng,pnel,rhoni
    ! on grid for numerical stability
    REAL(rkind),DIMENSION(0:nrmax) :: &
         mshearg,rp_totdg,dvexbpdrg,  &
         chi_cdbm,chim_cdbm,chi_tcdbm,chim_tcdbm

    dtr_tb(1:3*nsamax,1:3*nsamax,0:nrmax) = 0.D0
    vtr_tb(1:3*nsamax,1:3*nsamax,0:nrmax) = 0.D0

    ! model : Model ID
    !     0 : CDBM original
    !     1 : CDBM05 including elongation
    !     2 : CDBM original with weak ExB shear
    !     3 : CDBM05 with weak ExB shear
    !     4 : CDBM original with strong ExB shear
    !     5 : CDBM05 with strong ExB shear       
    model  = mdltr_tb - 130

    ! Tuned CDBM switch by M.Yagi ( 0: off, 1: on )
    model2 = 1

    ! Factor in s-slpha effects [1.0]
    calf=1.D0
    ! Factor in magnetic curvature effects [1.0]
    ckap=1.D0
    ! Factor in ExB drift effects [1.0]
    cexb=1.D0

    qp_eng = qp(nrmax) ! engineering safety factor for Tuned CDBM

    call mesh_convert_mtog(mshear(1:nrmax),mshearg(0:nrmax),nrmax)
    call mesh_convert_mtog(rp_totd(1:nrmax),rp_totdg(0:nrmax),nrmax)
    call mesh_convert_mtog(dvexbpdr(1:nrmax),dvexbpdrg(0:nrmax),nrmax)

    DO nr = 1, nrmax
       pnel = rn_e(nr)*1.d20  ! electron density [m^-3]
    
       rhoni = 0.D0
       DO nsa = 1, nsamax
          IF(idnsa(nsa)==1) THEN
             ns = ns_nsa(nsa)
             ! ion mass density [kg/m^3]
             ! (sum of ion-mass times ion-density)
             rhoni= rhoni + pa(ns)*amp*rn(nsa,nr)*1.D20
          END IF
       END DO

       CALL cdbm(BB,RR,rmnrho(nr),rkprho(nr),qp(nr),qp_eng,mshearg(nr), &
            pnel,rhoni,rp_totdg(nr),dvexbpdrg(nr), &
            calf,ckap,cexb,model,model2,chi_cdbm(nr),chi_tcdbm(nr))

    END DO

    call data_interpolate_gtom(chi_cdbm(0:nrmax),chim_cdbm(1:nrmax),nrmax)
    call data_interpolate_gtom(chi_tcdbm(0:nrmax),chim_tcdbm(1:nrmax),nrmax)

    DO nr = 1, nrmax
       IF(model2 == 0)THEN
          DO nsa = 1, nsamax
             IF(idnsa(nsa) /= 0) THEN
                dtr_tb(3*nsa-2,3*nsa-2,nr) = cdtrn*chim_cdbm(nr)
                dtr_tb(3*nsa-1,3*nsa-1,nr) = cdtru*chim_cdbm(nr)
                dtr_tb(3*nsa  ,3*nsa  ,nr) = cdtrt*chim_cdbm(nr)
!                write(*,*) chi_cdbm
             ENDIF
          END DO
       ELSE IF(model2 == 1)THEN
          DO nsa = 1, nsamax
             IF(idnsa(nsa) /= 0) THEN
                dtr_tb(3*nsa-2,3*nsa-2,nr) = cdtrn*chim_tcdbm(nr)
                dtr_tb(3*nsa-1,3*nsa-1,nr) = cdtru*chim_tcdbm(nr)
                dtr_tb(3*nsa  ,3*nsa  ,nr) = cdtrt*chim_tcdbm(nr)
             ENDIF
          END DO
       END IF

    END DO
    RETURN

  END SUBROUTINE tr_cdbm
END MODULE trcdbm
