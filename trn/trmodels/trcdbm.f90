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
         ikind,rkind,bb,rr,abb1rho,rmnrho,rmjrho,rkev,pa,amp,aee, &
         rkprho,cdtrn,cdtru,cdtrt,qp,rn,                          &
         dtr_tb,vtr_tb,nrmax,nsamax,neqmax,ns_nsa,idnsa,mdltr_tb

    USE trcalv, ONLY: &
         rn_e,     &! the density of electron
         rp_totd , &! the total pressure
         mshear,   &! magnetic shear
         dvexbdr    ! the gradient of ExBp drift velocity

    USE trlib, ONLY: mesh_convert_mtog,data_interpolate_gtom

    IMPLICIT NONE
    INTEGER(ikind):: nr,ns,nsa,model,model_t
    REAL(rkind):: calf,ckap,cexb,tc1,tc2,tk,pnel,rhoni
    REAL(rkind),DIMENSION(0:nrmax) :: &
         chi_cdbm,chim_cdbm

    chi_cdbm(0:nrmax) = 0.d0
    dtr_tb(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    vtr_tb(1:neqmax,1:neqmax,0:nrmax) = 0.d0

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

       ! engineering safety factor for Tuned CDBM
       ! *** for the time being ***
       IF(model_t == 1) qp(nr) = qp(nrmax) 

       CALL cdbm(abb1rho(nr),rmjrho(nr),rmnrho(nr),rkprho(nr),       &
            qp(nr),mshear(nr),pnel,rhoni,rp_totd(nr),dvexbdr(nr),   &
            calf,ckap,cexb,model,chi_cdbm(nr),                       &
            mdl_tuned=model_t,c1_tuned=tc1,c2_tuned=tc2,k_tuned=tk )

    END DO

    ! on grid -> on half grid
    chim_cdbm(1:nrmax) = 0.5d0*(chi_cdbm(0:nrmax-1)+chi_cdbm(1:nrmax))

    DO nr = 1, nrmax
          DO nsa = 1, nsamax
             IF(idnsa(nsa) /= 0) THEN
                dtr_tb(1+3*nsa-2,1+3*nsa-2,nr) = cdtrn*chim_cdbm(nr)
                dtr_tb(1+3*nsa-1,1+3*nsa-1,nr) = cdtru*chim_cdbm(nr)
                dtr_tb(1+3*nsa  ,1+3*nsa  ,nr) = cdtrt*chim_cdbm(nr)
!                write(*,*) chi_cdbm
             ENDIF
          END DO
    END DO

    RETURN

  END SUBROUTINE tr_cdbm
END MODULE trcdbm
