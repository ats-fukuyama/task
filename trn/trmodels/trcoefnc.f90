MODULE trcoefnc

  USE trcomm, ONLY: ikind, rkind

  PRIVATE
  PUBLIC tr_coefnc

CONTAINS

  SUBROUTINE tr_coefnc
    USE trcomm, ONLY: &
         rkev,nsa_neq,nva_neq,nsab_nsa,neqmax,nsamax,nrmax,mdltr_nc,   &
         rhog,rhom,rmnrho,htr,eta,dtr_nc,vtr_nc,eta_nc,                &
         jbs_nc,jex_nc,rt,rn,vpol,vpar,vprp,cdtrn,cdtru,cdtrt!,nrd1,nrd2
    USE trncls, ONLY: tr_nclass, tr_ncls_allocate,                     &
         chi_ncp,chi_nct,d_ncp,d_nct,gfls,qfls,fls_tot,                &
         vebs,qebs,dia_gdnc,dia_gvnc,dia_qdnc,dia_qvnc,cjbs_p,cjbs_t,  &
         eta_ncls,jbs_ncls,jex_ncls,vpol_ncls,vpar_ncls,vprp_ncls

    USE libgrf, ONLY:grd1d

    IMPLICIT NONE
    REAL(rkind) :: deriv4,FCTR

    ! internal parameters
    INTEGER(ikind) :: nr,nsa,nsa1,nva,neq,ierr
    INTEGER(ikind) :: nclass_save = 0
    REAL(rkind)    :: cm,cp
    REAL(rkind),DIMENSION(0:nrmax) :: drhog,logetam,logetag

    dtr_nc(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    vtr_nc(1:neqmax,1:neqmax,0:nrmax) = 0.d0


    SELECT CASE(mdltr_nc) ! results are in dtr_nc and vtr_nc
    CASE(0)
       ! no transport
       
       nclass_save = 0 ! reset for NCLASS (input of eta)
    CASE(1)
       CALL tr_ncls_allocate

       IF(nclass_save==0)THEN
          eta_ncls(0:nrmax) = eta(0:nrmax)
       END IF

       CALL tr_nclass(ierr)

       ! diffusion and convection coefficients
       DO neq = 1, neqmax
          nsa  = nsa_neq(neq)
          nva  = nva_neq(neq)
          IF(nsa == 0) CYCLE
          IF(nsab_nsa(nsa) == 0) CYCLE ! only for bulk species

          IF(nva == 1)THEN! particle

             vtr_nc(neq,neq,1:nrmax) =                            &
                  0.5d0*(vebs(nsa,0:nrmax-1)+vebs(nsa,1:nrmax))   &
                + 0.5d0*(gfls(5,nsa,0:nrmax-1)/rn(nsa,0:nrmax-1)  &
                        +gfls(5,nsa,1:nrmax  )/rn(nsa,1:nrmax  ))
             DO nsa1 = 1, nsamax
                dtr_nc(neq,3*nsa1-1,1:nrmax) = cdtrn *                     &
                 0.5d0*(d_ncp(nsa,nsa1,0:nrmax-1)+d_ncp(nsa,nsa1,1:nrmax)  &
                       +d_nct(nsa,nsa1,0:nrmax-1)+d_nct(nsa,nsa1,1:nrmax))
             END DO

          ELSE IF(nva == 2)THEN! velocity
             CONTINUE

          ELSE IF(nva == 3)THEN! energy

             vtr_nc(neq,neq,1:nrmax) = &
                  0.5d0*(qebs(nsa,0:nrmax-1)+qebs(nsa,1:nrmax))  &
                + 0.5d0*(qfls(5,nsa,0:nrmax-1)/rn(nsa,0:nrmax-1) &
                        +qfls(5,nsa,1:nrmax  )/rn(nsa,1:nrmax  ))
             DO nsa1 = 1, nsamax
                dtr_nc(neq,3*nsa1+1,1:nrmax) = cdtrt *                       &
                0.5d0*(chi_ncp(nsa,nsa1,0:nrmax-1)+chi_ncp(nsa,nsa1,1:nrmax) &
                      +chi_nct(nsa,nsa1,0:nrmax-1)+chi_nct(nsa,nsa1,1:nrmax))
             END DO

          END IF
       END DO

       eta_nc(0:nrmax) = eta_ncls(0:nrmax)
       jbs_nc(0:nrmax) = jbs_ncls(0:nrmax)
       jex_nc(0:nrmax) = jex_ncls(0:nrmax)
       vpol(0:nrmax) = vpol_ncls(0:nrmax)
       vpar(0:nrmax) = vpar_ncls(0:nrmax)
       vprp(0:nrmax) = vprp_ncls(0:nrmax)

       nclass_save = 1
    END SELECT

    ! substitution
    IF(mdltr_nc /= 0)THEN
       ! resistivity
       eta(0:nrmax) = eta_nc(0:nrmax)
    END IF

!!$    write(6,*) eta
!!$    CALL PAGES
!!$    CALL GRD1D(1,rhog,eta,nrmax+1,nrmax+1,1,'eta',2)
!!$    CALL PAGEE

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

END MODULE trcoefnc
