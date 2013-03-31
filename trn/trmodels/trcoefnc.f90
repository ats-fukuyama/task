MODULE trcoefnc

  USE trcomm, ONLY: ikind, rkind

  PRIVATE
  PUBLIC tr_coefnc

CONTAINS

  SUBROUTINE tr_coefnc
    USE trcomm, ONLY: &
         rkev,nsa_neq,nva_neq,nsab_nsa,neqmax,nsamax,nrmax,mdltr_nc,   &
         rhog,rhom,rmnrho,htr,eta,dtr_nc,vtr_nc,eta_nc,                &
         jbs_nc,jex_nc,jbs,jex,rt,rn,vpol,vpar,vprp,cdtrn,cdtru,cdtrt
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
       eta_nc(1:nrmax) = 0.5d0*(eta(0:nrmax-1)+eta(1:nrmax))
       nclass_save = 0 ! reset for NCLASS (input of eta)
    CASE(1)
       CALL tr_ncls_allocate

       IF(nclass_save==0)THEN
          eta_ncls(1:nrmax) = 0.5d0*(eta(0:nrmax-1)+eta(1:nrmax))
       END IF

       CALL tr_nclass(ierr)

       ! diffusion and convection coefficients
       DO neq = 1, neqmax
          nsa  = nsa_neq(neq)
          nva  = nva_neq(neq)
          IF(nsa == 0) CYCLE
          IF(nsab_nsa(nsa) == 0) CYCLE ! only for bulk species

          IF(nva == 1)THEN! particle

             vtr_nc(neq,neq,1:nrmax) = vebs(nsa,1:nrmax) &
                                     + gfls(5,nsa,1:nrmax)/rn(nsa,1:nrmax)
             DO nsa1 = 1, nsamax
                dtr_nc(neq,3*nsa1-1,1:nrmax) = &
                     cdtrn*(d_ncp(nsa,nsa1,1:nrmax) + d_nct(nsa,nsa1,1:nrmax))
             END DO

          ELSE IF(nva == 2)THEN! velocity
             CONTINUE

          ELSE IF(nva == 3)THEN! energy

             vtr_nc(neq,neq,1:nrmax) = qebs(nsa,1:nrmax) &
                                     + qfls(5,nsa,1:nrmax)/rn(nsa,1:nrmax)
             DO nsa1 = 1, nsamax
                dtr_nc(neq,3*nsa1+1,1:nrmax) = &
                cdtrt*(chi_ncp(nsa,nsa1,1:nrmax) + chi_nct(nsa,nsa1,1:nrmax))
             END DO

          END IF
       END DO

       ! store for caluculation of magnetic diffusion coefficient (on HALF grid)
       eta_nc(1:nrmax) = eta_ncls(1:nrmax)
       logetam(1:nrmax) = LOG10(eta_ncls(1:nrmax)) ! for graphic

       jbs_nc(1:nrmax) = jbs_ncls(1:nrmax)
       jex_nc(1:nrmax) = jex_ncls(1:nrmax)

       ! on grid
       ! *** extrapolate center value ***                            
       logetag(0) = FCTR(rhom(1),rhom(2),logetam(1),logetam(2))
       eta(0) = FCTR(rhom(1),rhom(2),eta_ncls(1),eta_ncls(2))
       jbs(0) = FCTR(rhom(1),rhom(2),jbs_ncls(1),jbs_ncls(2))
       jex(0) = FCTR(rhom(1),rhom(2),jex_ncls(1),jex_ncls(2))
       vpol(0) = FCTR(rhom(1),rhom(2),vpol_ncls(1),vpol_ncls(2))
       vpar(0) = FCTR(rhom(1),rhom(2),vpar_ncls(1),vpar_ncls(2))
       vprp(0) = FCTR(rhom(1),rhom(2),vprp_ncls(1),vprp_ncls(2))
       
       ! *** interpolate values on grids linearly ***                
       drhog(1:nrmax) = rhog(1:nrmax) - rhog(0:nrmax-1)
       DO nr = 1, nrmax-1
          cm = drhog(nr+1) / (drhog(nr)+drhog(nr+1))
          cp = drhog(nr  ) / (drhog(nr)+drhog(nr+1))
          logetag(nr) = cm*logetam(nr) + cp*logetam(nr+1)
          eta(nr) = cm*eta_ncls(nr) + cp*eta_ncls(nr+1)
          jbs(nr) = cm*jbs_ncls(nr) + cp*jbs_ncls(nr+1)
          jex(nr) = cm*jex_ncls(nr) + cp*jex_ncls(nr+1)
          vpol(nr) = cm*vpol_ncls(nr) + cp*vpol_ncls(nr+1)
          vpar(nr) = cm*vpar_ncls(nr) + cp*vpar_ncls(nr+1)
          vprp(nr) = cm*vprp_ncls(nr) + cp*vprp_ncls(nr+1)
       END DO
       
       ! *** extrapolate edge value ***
       logetag(nrmax) = 2.d0*logetam(nrmax) - logetam(nrmax-1)
       eta(nrmax) = 2.d0*eta_ncls(nrmax) - eta_ncls(nrmax-1)
       jbs(nrmax) = 2.d0*jbs_ncls(nrmax) - jbs_ncls(nrmax-1)
       jex(nrmax) = 2.d0*jex_ncls(nrmax) - jex_ncls(nrmax-1)
       vpol(nrmax) = 2.d0*vpol_ncls(nrmax) - vpol_ncls(nrmax-1)
       vpar(nrmax) = 2.d0*vpar_ncls(nrmax) - vpar_ncls(nrmax-1)
       vprp(nrmax) = 2.d0*vprp_ncls(nrmax) - vprp_ncls(nrmax-1)

       nclass_save = 1
    END SELECT

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
