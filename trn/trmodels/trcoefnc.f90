MODULE trcoefnc

  USE trcomm, ONLY: ikind, rkind

  PRIVATE
  PUBLIC tr_coefnc

CONTAINS

  SUBROUTINE tr_coefnc
    USE trcomm, ONLY: &
         rkev,nsa_neq,nva_neq,neqmax,nsamax,nrmax,mdltr_nc,  &
         rmnrho,htr,eta,dtr_nc,vtr_nc,eta_nc,jbs_nc,jex_nc,  &
         rt,rn ,nrd1,nrd2
    USE trcalv, ONLY: &
         chi_ncp,chi_nct,d_ncp,d_nct,gfls,qfls,fls_tot, &
         vebs,qebs,dia_gdnc,dia_gvnc,cjbs_p,cjbs_t
    USE trncls, ONLY: tr_nclass

    IMPLICIT NONE
    REAL(rkind) :: deriv3

    ! internal parameters
    INTEGER(ikind) :: nr,nsa,nva,neq,nk,ierr
    REAL(rkind),DIMENSION(0:nrmax) :: rt_s,d_nc,v_nc

    dtr_nc(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    vtr_nc(1:neqmax,1:neqmax,0:nrmax) = 0.d0

    ! CALL tr_calnc
    ! trapped particle fraction
    

    SELECT CASE(mdltr_nc) ! results are in dtr_nc and vtr_nc
    CASE(0)
       ! no transport
    CASE(1)
       CALL tr_nclass(ierr)

       ! diffusion and convection coefficients
       ! ** diagonal term only
       fls_tot(3,1:nsamax,0:nrmax) = 0.d0
       DO neq = 1, neqmax
          nsa = nsa_neq(neq)
          nva = nva_neq(neq)
         
          IF(nva == 1)THEN! particle
             DO nk = 1, 5
                fls_tot(nva,nsa,0:nrmax) = fls_tot(nva,nsa,0:nrmax) &
                                            + gfls(nk,nsa,0:nrmax)
             END DO
             ! on half grid
             dtr_nc(neq,neq,1:nrmax) = &
                  0.5d0*(dia_gdnc(nsa,0:nrmax-1)+dia_gdnc(nsa,1:nrmax))
             vtr_nc(neq,neq,1:nrmax) = &
                  0.5d0*(dia_gvnc(nsa,0:nrmax-1)+dia_gvnc(nsa,1:nrmax))

          ELSE IF(nva == 2)THEN! velocity
!             dtr_nc(neq,neq,0:nrmax) = dia_dnc(nsa,0:nrmax)
!             vtr_nc(neq,neq,0:nrmax) = dia_vnc(nsa,0:nrmax)

          ELSE IF(nva == 3)THEN! energy
             DO nk = 1, 5
                fls_tot(nva,nsa,0:nrmax) = fls_tot(nva,nsa,0:nrmax) &
                                            + qfls(nk,nsa,0:nrmax)
             END DO             

             rt_s(0:nrmax) = rt(nsa,0:nrmax)
             DO nr = 0, nrmax
                d_nc(nr) = chi_nct(nsa,nsa,nr)+chi_ncp(nsa,nsa,nr)
                
                v_nc(nr) &
                  = fls_tot(nva,nsa,nr)                                    &
                    /(rn(nsa,nr)*rt(nsa,nr)*rkev)                          &
                  + d_nc(nr)*deriv3(nr,rmnrho,rt_s,nrmax,0) &
                    /rt(nsa,nr)
             END DO
             dtr_nc(neq,neq,1:nrmax) = &
                  0.5d0*(d_nc(0:nrmax-1)+d_nc(1:nrmax))
             ! interim way of substitution V_Es = V_Ks + (3/2)V_s
             vtr_nc(neq,neq,1:nrmax) =                  &
                  0.5d0*(v_nc(0:nrmax-1)+v_nc(1:nrmax)) &
                 +1.5d0*vtr_nc(neq-2,neq-2,1:nrmax)
          END IF
       END DO

       nrd1(0:nrmax) = dia_gdnc(1,0:nrmax)
       nrd2(0:nrmax) = dia_gvnc(1,0:nrmax)

       ! resistivity
       eta(0:nrmax) = eta_nc(0:nrmax)

       ! *** off diagonal term...       

    END SELECT

    ! bootstrap current
!    htr(1,0:nrmax) = jbs_nc(0:nrmax) + jex_nc(0:nrmax)   

    
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
