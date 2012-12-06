MODULE trcoeftb

  USE trcomm, ONLY: ikind, rkind

  PRIVATE
  PUBLIC tr_coeftb, Pereverzev_check

CONTAINS

! ***** calculate trubenlent transport coefficients *****

  SUBROUTINE tr_coeftb
    USE trcomm, ONLY: mdltr_tb,mdltr_prv,neqmax,nrmax, &
                      dtr_tb,vtr_tb,dtr_prv,vtr_prv
    USE trsimple, ONLY: tr_simple
    USE trglf23, ONLY: tr_glf23
    USE trcdbm, ONLY: tr_cdbm
    USE trmbgb, ONLY: tr_mbgb
    USE trmmm95,ONLY: tr_mmm95
    USE trmmm7_1,ONLY: tr_mmm7_1
    IMPLICIT NONE

    dtr_tb(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    vtr_tb(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    dtr_prv(1:neqmax,0:nrmax) = 0.d0
    vtr_prv(1:neqmax,0:nrmax) = 0.d0

    SELECT CASE(mdltr_tb) ! results are in dtr_tb and vtr_tb
    CASE(0:9)
       CALL tr_simple
    CASE(60:69)
       CALL tr_glf23
    CASE(130:139)
       CALL tr_cdbm
    CASE(140:149)
       CALL tr_mbgb
    CASE(150:159)
       CALL tr_mmm95
    CASE(160:169)
       CALL tr_mmm7_1
    END SELECT

    ! only for energy transport for now
    IF(mdltr_prv /= 0) CALL Pereverzev_method

    RETURN
  END SUBROUTINE tr_coeftb

!*************************************************************************

! ------------------------- Pereverzev method ----------------------------
!
!   Numerical stabilization method for stiff turbulent transport models*
!    (especially for GLF23 model and Weiland model etc...)
!   This method takes into account only energy trasport.
!
!   'Pereverzev check' subroutine is for calculating additional numerical
!    factor in element equations of FEM.
!   
!   Switch variable
!         mdltr_prv = 0 : off
!                   = 1 : D_add = dprv1
!                   = 2 : D_add = dprv2*dtr_tb(nr)
!                   = 3 : D_add = dprv2*dtr_tb(nr) + dprv1
!
! * G.V.Pereverzev, G.Corrigan, Computer Physics Comm. 179 (2008) 579-585
!
! ------------------------------------------------------------------------
  SUBROUTINE Pereverzev_method
    USE trcomm, ONLY: ikind,rkind,nrmax,nsamax,idnsa,ph0,rg,rhog,      &
         mdltr_prv,dprv1,dprv2,rhog_prv,dtr,vtr,nsa_neq,dtr_tb,vtr_tb, &
         dtr_prv,vtr_prv,cdtrn,cdtru,cdtrt,rn,ru,rt,                   &
         rn_prev,ru_prev,rt_prev,ar1rho,ar2rho,dvrho
    IMPLICIT NONE
    REAL(rkind) :: &
         dtr_new,vtr_old,lt,drnrt,rnrt_ave,dvdrp,dvdrm,dtr_elm,vtr_elm,&
         gm1p,gm1m,gm2p,gm2m,gm20
    INTEGER(ikind) :: nr,nsa

    DO nr = 1, nrmax

       IF(mdltr_prv > 10)THEN
          IF(rhog(nr) < rhog_prv) CYCLE
       END IF

       dvdrp = dvrho(nr  )
       dvdrm = dvrho(nr-1)

       gm1p = dvdrp*ar1rho(nr  )
       gm1m = dvdrm*ar1rho(nr-1)

       gm2p = dvdrp*ar2rho(nr  )
       gm2m = dvdrm*ar2rho(nr-1)
       gm20 = gm2p + gm2m

       DO nsa = 1, nsamax
          ! excluding neutral and fast ions
          IF(idnsa(nsa) == 0 .OR. idnsa(nsa) == 2) CYCLE

          drnrt = rn_prev(nsa,nr  )*rt_prev(nsa,nr  ) &
                 -rn_prev(nsa,nr-1)*rt_prev(nsa,nr-1)

          rnrt_ave=((2.d0*gm1m +      gm1p)                &
                     *rn_prev(nsa,nr-1)*rt_prev(nsa,nr-1)  &
                   +(     gm1m + 2.d0*gm1p)                &
                     *rn_prev(nsa,nr  )*rt_prev(nsa,nr ))  &
                   /(3.d0*gm20)

!          write(*,*) '*** drnrt ***', drnrt
!          write(*,*) '*** rnrt_ave ***', rnrt_ave

          IF(drnrt > 0.D0) THEN
             lt = 0.D0
          ELSE
             lt = drnrt/(rhog(nr)-rhog(nr-1))
          END IF
          SELECT CASE(mdltr_prv)
          CASE(1,11)
             dtr_new = dprv1
          CASE(2,12)
             dtr_new = dprv2*dtr_tb(3*nsa,3*nsa,nr)
          CASE(3,13)
             dtr_new = dprv2*dtr_tb(3*nsa,3*nsa,nr)+dprv1
          CASE(4,14)
             dtr_new = dprv1 + dprv2*(0.5d0*(rhog(nr)+rhog(nr-1)))**3
!!$          CASE(4)
!!$             IF(nr==1)THEN
!!$                dtr_new = dprv2 &
!!$                         *0.5d0**(dtr_tb(3*nsa,3*nsa,nr  )   &
!!$                                 +dtr_tb(3*nsa,3*nsa,nr+1))
!!$             ELSE IF(nr==nrmax)THEN
!!$                dtr_new = dprv2 &
!!$                         *0.5d0**(dtr_tb(3*nsa,3*nsa,nr-1)   &
!!$                                 +dtr_tb(3*nsa,3*nsa,nr  ))
!!$             ELSE
!!$                dtr_new = dprv2 &                  
!!$                         *0.33d0*(dtr_tb(3*nsa,3*nsa,nr-1)   &
!!$                                 +dtr_tb(3*nsa,3*nsa,nr  )   &
!!$                                 +dtr_tb(3*nsa,3*nsa,nr+1))
!!$             END IF
          END SELECT

!          dtr_prv(1+3*nsa-2,nr) = cdtrn*dtr_new
!          dtr_prv(1+3*nsa-1,nr) = cdtru*dtr_new
          dtr_prv(1+3*nsa  ,nr) = cdtrt*dtr_new

          dtr_tb(1+3*nsa-2,1+3*nsa-2,nr) = dtr_tb(1+3*nsa-2,1+3*nsa-2,nr) &
                                         +dtr_prv(1+3*nsa-2,nr)
          dtr_tb(1+3*nsa-1,1+3*nsa-1,nr) = dtr_tb(1+3*nsa-1,1+3*nsa-1,nr) &
                                         +dtr_prv(1+3*nsa-1,nr) 
          dtr_tb(1+3*nsa  ,1+3*nsa  ,nr) = dtr_tb(1+3*nsa  ,1+3*nsa  ,nr) &
                                         +dtr_prv(1+3*nsa  ,nr)

          vtr_old = dtr_new * lt / rnrt_ave
!          write(*,*) '*** nr, vtr_old ***', nr, vtr_old

!          vtr_prv(1+3*nsa-2,nr) = cdtrn*vtr_old
!          vtr_prv(1+3*nsa-1,nr) = cdtru*vtr_old
          vtr_prv(1+3*nsa  ,nr) = cdtrt*vtr_old

          vtr_tb(1+3*nsa-2,1+3*nsa-2,nr) = vtr_tb(1+3*nsa-2,1+3*nsa-2,nr) &
                                         +vtr_prv(1+3*nsa-2,nr) 
          vtr_tb(1+3*nsa-1,1+3*nsa-1,nr) = vtr_tb(1+3*nsa-1,1+3*nsa-1,nr) &
                                         +vtr_prv(1+3*nsa-1,nr) 
          vtr_tb(1+3*nsa  ,1+3*nsa  ,nr) = vtr_tb(1+3*nsa  ,1+3*nsa  ,nr) &
                                         +vtr_prv(1+3*nsa  ,nr) 

       END DO
    END DO
    RETURN
  END SUBROUTINE Pereverzev_method

  SUBROUTINE Pereverzev_check(add_prv)
    USE trcomm, ONLY: ikind,rkind,neqmax,nsamax,nrmax,idnsa,rg, &
         mdltr_prv,dtr_prv,vtr_prv,dtr,rn,ru,rt,dvrho,ar1rho,ar2rho,rhog

    REAL(rkind),DIMENSION(neqmax,0:nrmax) :: dtr_elm,vtr_elm,dtr_all
    REAL(rkind),DIMENSION(neqmax,0:nrmax),INTENT(OUT) :: add_prv
    REAL(rkind) :: term_n,term_u,term_t
    REAL(rkind) :: dvdrm,dvdrp,gm1p,gm1m,gm2p,gm2m,gm20
    INTEGER(ikind) :: nr,nsa

    add_prv(1:neqmax,0:nrmax) = 0.d0    
    IF(mdltr_prv == 0) RETURN
    
    DO nr = 1, nrmax
       dvdrp = dvrho(nr)
       dvdrm = dvrho(nr-1)

       gm1p = dvdrp*ar1rho(nr)
       gm1m = dvdrm*ar1rho(nr-1)

       gm2p = dvdrp*ar2rho(nr)
       gm2m = dvdrm*ar2rho(nr-1)
       gm20 = 0.5d0*(gm2p+gm2m)
       
       DO nsa = 1, nsamax
          ! excluding neutral and fast ions
          IF(idnsa(nsa) == 0 .OR. idnsa(nsa) == 2) CYCLE

          term_n = gm20*(rn(nsa,nr)-rn(nsa,nr-1))/(rhog(nr)-rhog(nr-1))
          term_u = gm20*(ru(nsa,nr)-ru(nsa,nr-1))/(rhog(nr)-rhog(nr-1))
          term_t = gm20*(rn(nsa,nr)*rt(nsa,nr)-rn(nsa,nr-1)*rt(nsa,nr-1)) &
                                                 /(rhog(nr)-rhog(nr-1))

          dtr_elm(1+3*nsa-2,nr) = dtr_prv(1+3*nsa-2,nr) * term_n
          dtr_elm(1+3*nsa-1,nr) = dtr_prv(1+3*nsa-1,nr) * term_u
          dtr_elm(1+3*nsa  ,nr) = dtr_prv(1+3*nsa  ,nr) * term_t

          ! diagonal element only
          dtr_all(1+3*nsa-2,nr) = dtr(1+3*nsa-2,1+3*nsa-2,nr) * term_n
          dtr_all(1+3*nsa-1,nr) = dtr(1+3*nsa-1,1+3*nsa-1,nr) * term_u
          dtr_all(1+3*nsa  ,nr) = dtr(1+3*nsa  ,1+3*nsa  ,nr) * term_t

          vtr_elm(1+3*nsa-2,nr) =                            &
          vtr_prv(1+3*nsa-2,nr)/6.D0                         &
            *((2.D0*gm1m+gm1p)*rn(nsa,nr-1) + (gm1m+2.D0*gm1p)*rn(nsa,nr))
          vtr_elm(1+3*nsa-1,nr) =                            &
          vtr_prv(1+3*nsa-1,nr)/6.D0                         &
            *((2.D0*gm1m+gm1p)*ru(nsa,nr-1) + (gm1m+2.D0*gm1p)*ru(nsa,nr))
          vtr_elm(1+3*nsa  ,nr) =                                 &
          vtr_prv(1+3*nsa  ,nr)/6.D0                              &
            *((2.D0*gm1m+     gm1p)*rn(nsa,nr-1)*rt(nsa,nr-1)     &
             +(     gm1m+2.D0*gm1p)*rn(nsa,nr  )*rt(nsa,nr  ))


          IF(dtr_all(1+3*nsa-2,nr)==0.d0) dtr_all(1+3*nsa-2,nr) = 0.01d0
          IF(dtr_all(1+3*nsa-1,nr)==0.d0) dtr_all(1+3*nsa-1,nr) = 0.01d0
          IF(dtr_all(1+3*nsa  ,nr)==0.d0) dtr_all(1+3*nsa  ,nr) = 0.01d0

          ! numerically additional term in nodal equation
          add_prv(1+3*nsa-2,nr) =                            &
               (dtr_elm(1+3*nsa-2,nr)-vtr_elm(1+3*nsa-2,nr)) &
               / dtr_all(1+3*nsa-2,nr)
          add_prv(1+3*nsa-1,nr) =                            &
               (dtr_elm(1+3*nsa-1,nr)-vtr_elm(1+3*nsa-1,nr)) &
               / dtr_all(1+3*nsa-1,nr)
          add_prv(1+3*nsa  ,nr) =                            &
               (dtr_elm(1+3*nsa  ,nr)-vtr_elm(1+3*nsa  ,nr)) &
               / dtr_all(1+3*nsa  ,nr)          
       END DO
    END DO

  END SUBROUTINE Pereverzev_check
END MODULE trcoeftb
