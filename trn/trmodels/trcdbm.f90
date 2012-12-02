MODULE trcdbm

  PRIVATE
  PUBLIC tr_cdbm

CONTAINS

!  ************************************************************
!
!   The interface between TASK/TR(trcoeftb) and CDBM model(cdbm)
!
!  ************************************************************

  SUBROUTINE tr_cdbm

    USE cdbm_mod, ONLY: cdbm
    USE trcomm, ONLY: &
         ikind,rkind,rkev,amp,aee,pa,nrmax,nsamax,neqmax,ns_nsa,idnsa, &
         bb,rr,ra,abb1rho,rmnrho,rmjrho,rhog,rhom,rkprho,              &
         cdtrn,cdtru,cdtrt,qp,rn,dtr_tb,vtr_tb,mdltr_tb,               &
         rn_e,     &! the density of electron
         rn_i,     &! the effective density of hydrogenic ions
         rp_totd,  &! the total pressure gradient
         ai_ave,   &! mean atomic mass of tehrmal ions [AMU]
         mshear,   &! magnetic shear
         dvexbpdr   ! the gradient of ExBp drift velocity

    IMPLICIT NONE
    INTEGER(ikind):: nr,ns,nsa,model,model_t,ierr
    REAL(rkind):: calf,ckap,cexb,tc1,tc2,tk,pnel,rhoni
    REAL(rkind),DIMENSION(0:nrmax) :: chim_cdbm

    REAL(rkind) :: abb1rhom,rmjrhom,rmnrhom,rkprhom,qpm

    dtr_tb(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    vtr_tb(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    chim_cdbm(0:nrmax) = 0.d0

    ! model : Model ID
    !     0 : CDBM original
    !     1 : CDBM05 including elongation
    !     2 : CDBM original with weak ExB shear
    !     3 : CDBM05 with weak ExB shear
    !     4 : CDBM original with strong ExB shear
    !     5 : CDBM05 with strong ExB shear       
    model  = mdltr_tb - 130
    ! Tuned CDBM switch by M.Yagi ( 0: off, 1: on )
    model_t = 0

    ! Factor in s-slpha effects [1.0]
    calf=1.D0
    ! Factor in magnetic curvature effects [1.0]
    ckap=1.D0
    ! Factor in ExB drift effects [1.0]
    cexb=1.D0

    ! C1,C2,K for Tuned CDBM by M.Yagi
    tc1 = 0.5d0
    tc2 = 15.d0
    tk  = 3.d0

    DO nr = 1, nrmax
       pnel = 0.5d0*(rn_e(nr)+rn_e(nr-1))*1.d20  ! electron density [m^-3]

       ! ion mass density [kg/m^3]
       ! (sum of ion-mass times ion-density)
       rhoni = ai_ave(nr) * amp * 0.5d0*(rn_i(nr)+rn_i(nr-1))*1.d20

       abb1rhom = 0.5d0*(abb1rho(nr-1) + abb1rho(nr))
       rmjrhom  = 0.5d0*(rmjrho(nr-1)  +  rmjrho(nr))
       rmnrhom  = 0.5d0*(rmnrho(nr-1)  +  rmnrho(nr))
       rkprhom  = 0.5d0*(rkprho(nr-1)  +  rkprho(nr))
       qpm      = 0.5d0*(qp(nr-1)      +      qp(nr))

       ! engineering safety factor for Tuned CDBM
       ! *** for the time being ***
       IF(model_t == 1) qp(nr) = qp(nrmax) 

       CALL cdbm(abb1rhom,rmjrhom,rmnrhom,rkprhom,                   &
            qpm,mshear(nr),pnel,rhoni,rp_totd(nr),dvexbpdr(nr),      &
            calf,ckap,cexb,model,chim_cdbm(nr),                      &
            mdl_tuned=model_t,c1_tuned=tc1,c2_tuned=tc2,k_tuned=tk )

    END DO

    DO nr = 1, nrmax
          DO nsa = 1, nsamax
             IF(idnsa(nsa) /= 0) THEN
                dtr_tb(1+3*nsa-2,1+3*nsa-2,nr) = cdtrn*chim_cdbm(nr)
                dtr_tb(1+3*nsa-1,1+3*nsa-1,nr) = cdtru*chim_cdbm(nr)
                dtr_tb(1+3*nsa  ,1+3*nsa  ,nr) = cdtrt*chim_cdbm(nr)
             ENDIF
          END DO
    END DO

    RETURN

  END SUBROUTINE tr_cdbm
END MODULE trcdbm
