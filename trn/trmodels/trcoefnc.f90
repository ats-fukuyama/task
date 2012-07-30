MODULE trcoefnc

  USE trcomm, ONLY: ikind, rkind

  PRIVATE
  PUBLIC tr_coefnc

CONTAINS

  SUBROUTINE tr_coefnc
    USE trcomm, ONLY: &
         rkev,nsa_neq,nva_neq,neqmax,nsamax,nrmax,mdltr_nc,       &
         rhog,rmnrho,htr,eta,dtr_nc,vtr_nc,eta_nc,jbs_nc,jex_nc,rt,rn,   &
         chi_ncp,chi_nct,d_ncp,d_nct,gfls,qfls,fls_tot,           &
         vebs,qebs,dia_gdnc,dia_gvnc,cjbs_p,cjbs_t
        !,nrd1,nrd2
    USE trncls, ONLY: tr_nclass

    IMPLICIT NONE
    REAL(rkind) :: deriv4

    ! internal parameters
    INTEGER(ikind) :: nr,nsa,nva,neq,nk,ierr
    REAL(rkind),DIMENSION(0:nrmax) :: rt_s,d_nc,v_nc,drhog

    dtr_nc(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    vtr_nc(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    jbs_nc(0:nrmax) = 0.d0

    ! CALL tr_calnc
    ! trapped particle fraction
    

    SELECT CASE(mdltr_nc) ! results are in dtr_nc and vtr_nc
    CASE(0)
       ! no transport
       
    CASE(1)
       CALL tr_nclass(ierr)

       ! diffusion and convection coefficients
       ! ** For now, we consider diagonal terms only
       fls_tot(3,1:nsamax,0:nrmax) = 0.d0
       DO neq = 1, neqmax
          nsa = nsa_neq(neq)
          nva = nva_neq(neq)
         
          IF(nva == 1)THEN! particle
             DO nk = 1, 5
                fls_tot(nva,nsa,1:nrmax) = fls_tot(nva,nsa,1:nrmax) &
                                            + gfls(nk,nsa,1:nrmax)
             END DO
             ! on half grid
             dtr_nc(neq,neq,1:nrmax) = dia_gdnc(nsa,1:nrmax)
             vtr_nc(neq,neq,1:nrmax) = dia_gvnc(nsa,1:nrmax)

          ELSE IF(nva == 2)THEN! velocity
!             dtr_nc(neq,neq,0:nrmax) = dia_dnc(nsa,0:nrmax)
!             vtr_nc(neq,neq,0:nrmax) = dia_vnc(nsa,0:nrmax)

          ELSE IF(nva == 3)THEN! energy
             DO nk = 1, 5
                fls_tot(nva,nsa,1:nrmax) = fls_tot(nva,nsa,1:nrmax) &
                                            + qfls(nk,nsa,1:nrmax)
             END DO             

             DO nr = 1, nrmax
!                drhog(nr) = rmnrho(nr) - rmnrho(nr-1)
                drhog(nr) = rhog(nr) - rhog(nr-1)
                d_nc(nr)  = chi_nct(nsa,nsa,nr) + chi_ncp(nsa,nsa,nr)
                
                v_nc(nr) &
                  = fls_tot(nva,nsa,nr)                          &
                    /(0.5d0*(rn(nsa,nr)+rn(nsa,nr-1))            &
                     *0.5d0*(rt(nsa,nr)+rt(nsa,nr-1))*rkev)      &
                  + d_nc(nr)*(rt(nsa,nr)-rt(nsa,nr-1))/drhog(nr) &
                    /(0.5d0*(rt(nsa,nr)+rt(nsa,nr-1)))
             END DO
             dtr_nc(neq,neq,1:nrmax) = d_nc(1:nrmax)
             ! interim way of substitution V_Es = V_Ks + (3/2)V_s
             vtr_nc(neq,neq,1:nrmax) = v_nc(1:nrmax) &
                                      +2.5d0*vtr_nc(neq-2,neq-2,1:nrmax)
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

END MODULE trcoefnc
