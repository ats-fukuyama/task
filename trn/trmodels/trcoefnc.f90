MODULE trcoefnc

  USE trcomm, ONLY: ikind, rkind

  PRIVATE
  PUBLIC tr_coefnc

CONTAINS

  SUBROUTINE tr_coefnc
    USE trcomm, ONLY: &
         rkev,nsa_neq,nva_neq,nsab_nsa,neqmax,nsamax,nrmax,mdltr_nc,       &
         rhog,rhom,rmnrho,htr,eta,dtr_nc,vtr_nc,eta_nc,jbs_nc,jex_nc,rt,rn,&
         vpol,vpar,vprp!,nrd1,nrd2
    USE trncls, ONLY: tr_nclass, tr_ncls_allocate,       &
         chi_ncp,chi_nct,d_ncp,d_nct,gfls,qfls,fls_tot,  &
         vebs,qebs,dia_gdnc,dia_gvnc,cjbs_p,cjbs_t,      &
         eta_ncls,jbs_ncls,jex_ncls,vpol_ncls,vpar_ncls,vprp_ncls

    IMPLICIT NONE
    REAL(rkind) :: deriv4,FCTR

    ! internal parameters
    INTEGER(ikind) :: nr,nsa,nsab,nva,neq,nk,ierr
    INTEGER(ikind) :: nclass_save = 0
    REAL(rkind)    :: cm,cp
    REAL(rkind),DIMENSION(0:nrmax) :: rt_s,d_nc,v_nc,drhog

    dtr_nc(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    vtr_nc(1:neqmax,1:neqmax,0:nrmax) = 0.d0

    ! CALL tr_calnc
    

    SELECT CASE(mdltr_nc) ! results are in dtr_nc and vtr_nc
    CASE(0)
       ! no transport
       
       nclass_save = 0 ! reset for NCLASS (input of eta)
    CASE(1)
       CALL tr_ncls_allocate

       IF(nclass_save==0)THEN
          eta_ncls(1:nrmax) = 0.5d0*(eta(0:nrmax-1)+eta(1:nrmax))
       END IF
       CALL tr_nclass(ierr)

       ! diffusion and convection coefficients
       ! ** For now, we consider diagonal terms only
       fls_tot(3,1:nsamax,0:nrmax) = 0.d0
       DO neq = 1, neqmax
          nsa  = nsa_neq(neq)
          nva  = nva_neq(neq)
          IF(nsa == 0)THEN
             CYCLE
          ELSE
             nsab = nsab_nsa(nsa)
          END IF
          IF(nsab == 0) CYCLE
         
          IF(nva == 1)THEN! particle
             DO nk = 1, 5
                fls_tot(nva,nsab,1:nrmax) = fls_tot(nva,nsab,1:nrmax) &
                                          + gfls(nk,nsab,1:nrmax)
             END DO
             ! on half grid
             dtr_nc(neq,neq,1:nrmax) = dia_gdnc(nsab,1:nrmax)
             vtr_nc(neq,neq,1:nrmax) = dia_gvnc(nsab,1:nrmax)

          ELSE IF(nva == 2)THEN! velocity
!             dtr_nc(neq,neq,0:nrmax) = dia_dnc(nsab,0:nrmax)
!             vtr_nc(neq,neq,0:nrmax) = dia_vnc(nsab,0:nrmax)

          ELSE IF(nva == 3)THEN! energy
             DO nk = 1, 5
                fls_tot(nva,nsab,1:nrmax) = fls_tot(nva,nsab,1:nrmax) &
                                          + qfls(nk,nsab,1:nrmax)
             END DO

             DO nr = 1, nrmax
!                drhog(nr) = rmnrho(nr) - rmnrho(nr-1)
                drhog(nr) = rhog(nr) - rhog(nr-1)
                d_nc(nr)  = chi_nct(nsab,nsab,nr) + chi_ncp(nsab,nsab,nr)

                v_nc(nr) &
                     = fls_tot(nva,nsab,nr)                          &
                        /(0.5d0*(rn(nsa,nr)+rn(nsa,nr-1))            &
                         *0.5d0*(rt(nsa,nr)+rt(nsa,nr-1))*rkev)      &
                      + d_nc(nr)*(rt(nsa,nr)-rt(nsa,nr-1))/drhog(nr) &
                        /(0.5d0*(rt(nsa,nr)+rt(nsa,nr-1)))
!                nrd1(nr) = v_nc(nr)
!                v_nc(nr) = qebs(nsab,nr)
!                nrd2(nr) = v_nc(nr)
             END DO
             dtr_nc(neq,neq,1:nrmax) = d_nc(1:nrmax)
             ! interim way of substitution V_Es = V_Ks + (3/2)V_s
             vtr_nc(neq,neq,1:nrmax) = v_nc(1:nrmax)
          END IF
       END DO

!       nrd1(0:nrmax) = dia_gdnc(1,0:nrmax)
!       nrd2(0:nrmax) = dia_gvnc(1,0:nrmax)

       ! *** extrapolate center value ***                            
       eta_nc(0) = FCTR(rhom(1),rhom(2),eta_ncls(1),eta_ncls(2))
       jbs_nc(0) = FCTR(rhom(1),rhom(2),jbs_ncls(1),jbs_ncls(2))
       jex_nc(0) = FCTR(rhom(1),rhom(2),jex_ncls(1),jex_ncls(2))
!       eta_nc(0) = eta_ncls(1)                                    
!       jbs_nc(0) = jbs_ncls(1)                                    
!       jex_nc(0) = jex_ncls(1)                                    
       vpol(0) = FCTR(rhom(1),rhom(2),vpol_ncls(1),vpol_ncls(2))
       vpar(0) = FCTR(rhom(1),rhom(2),vpar_ncls(1),vpar_ncls(2))
       vprp(0) = FCTR(rhom(1),rhom(2),vprp_ncls(1),vprp_ncls(2))
       
       ! *** interpolate values on grids linearly ***                
       DO nr = 1, nrmax-1
          cm = drhog(nr+1) / (drhog(nr)+drhog(nr+1))
          cp = drhog(nr  ) / (drhog(nr)+drhog(nr+1))
          eta_nc(nr) = cm*eta_ncls(nr) + cp*eta_ncls(nr+1)
          jbs_nc(nr) = cm*jbs_ncls(nr) + cp*jbs_ncls(nr+1)
          jex_nc(nr) = cm*jex_ncls(nr) + cp*jex_ncls(nr+1)
          vpol(nr) = cm*vpol_ncls(nr) + cp*vpol_ncls(nr+1)
          vpar(nr) = cm*vpar_ncls(nr) + cp*vpar_ncls(nr+1)
          vprp(nr) = cm*vprp_ncls(nr) + cp*vprp_ncls(nr+1)
       END DO
       
       ! *** extrapolate edge value ***                              
       eta_nc(nrmax) = 2.d0*eta_ncls(nrmax) - eta_ncls(nrmax-1)
       jbs_nc(nrmax) = 2.d0*jbs_ncls(nrmax) - jbs_ncls(nrmax-1)
       jex_nc(nrmax) = 2.d0*jex_ncls(nrmax) - jex_ncls(nrmax-1)
       vpol(nrmax) = 2.d0*vpol_ncls(nrmax) - vpol_ncls(nrmax-1)
       vpar(nrmax) = 2.d0*vpar_ncls(nrmax) - vpar_ncls(nrmax-1)
       vprp(nrmax) = 2.d0*vprp_ncls(nrmax) - vprp_ncls(nrmax-1)

       ! *** off diagonal term...       

       nclass_save = 1
    END SELECT

    ! resistivity
    eta(0:nrmax) = eta_nc(0:nrmax)

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
