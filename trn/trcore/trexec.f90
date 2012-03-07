MODULE trexec
! This module solves the matrix equation

  USE trcomm,ONLY: rkind, ikind

  PUBLIC tr_exec
  PRIVATE

CONTAINS

  SUBROUTINE tr_exec

    USE trcomm,ONLY: &
         neqmax,neqrmax,nvmax,nvrmax,nrmax,dt,rhog, &
         dtr,vtr,ctr,str,htr,elmtx,limtx,rsimtx,rjimtx,r1imtx,r2imtx,r3imtx, &
         rimtx,lhmtx,rhv,xv,xv_new,xv_prev,neqr_neq,id_neq,id_neqnr
    IMPLICIT NONE
    REAL(rkind) :: dh0,dh1,dh2,dh3,dvdrp,dvdrm,ldrp,ldrm,rdrp,rdrm,  &
         gm1p,gm1m,gm2p,gm2m,cjexm,cjexp,coef1,coef2
    INTEGER(ikind) :: nr,neq,neq1,neqr,neqr1,nvrm,nvrp,ierr,nvm,nvp

    lhmtx(1:4*neqrmax-1,1:nvrmax) = 0.D0
    rhv(1:nvrmax)                 = 0.D0

    DO nr = 1, nrmax

       elmtx(1:2*neqrmax,1:2*neqrmax)=0.D0

       dh0 = (rhog(nr)-rhog(nr-1))/12.D0
       dh1 = dt/6.D0
       dh2 = dt*dh0
       dh3 = dt/(rhog(nr)-rhog(nr-1))

       DO neq = 1, neqmax

          call tr_neq_mtrc(nr,neq,dvdrp,dvdrm,ldrp,ldrm,rdrp,rdrm, &
               gm1p,gm1m,gm2p,gm2m,cjexm,cjexp,coef1,coef2)
          
          limtx(1,1,neq) = coef1*dh0*(3.D0*ldrm +      ldrp)
          limtx(2,1,neq) = coef1*dh0*(     ldrm +      ldrp)
          limtx(1,2,neq) = coef1*dh0*(     ldrm +      ldrp)
          limtx(2,2,neq) = coef1*dh0*(     ldrm + 3.D0*ldrp) 

          ! only for magnetic diffusion equation
          rjimtx(1,1,neq) = - dh1*(-4.D0*cjexm +      cjexp)
          rjimtx(2,1,neq) = - dh1*(-2.D0*cjexm -      cjexp)
          rjimtx(1,2,neq) = - dh1*(      cjexm + 2.D0*cjexp)
          rjimtx(2,2,neq) = - dh1*(-     cjexm + 4.D0*cjexp)

          rsimtx(1,1,neq) = dh2*(3.D0*dvdrm +      dvdrp)
          rsimtx(2,1,neq) = dh2*(     dvdrm +      dvdrp)
          rsimtx(1,2,neq) = dh2*(     dvdrm +      dvdrp)
          rsimtx(2,2,neq) = dh2*(     dvdrm + 3.D0*dvdrp) 

          DO neq1 = 1, neqmax
        
             ! at the center of element
             r1imtx(1,1,neq,neq1)=   dh3*dtr(neq,neq1,nr)   &
                                      *gm2m
             r1imtx(2,1,neq,neq1)= - dh3*dtr(neq,neq1,nr)   &
                                      *gm2m
             r1imtx(1,2,neq,neq1)= - dh3*dtr(neq,neq1,nr)   &
                                      *gm2p
             r1imtx(2,2,neq,neq1)=   dh3*dtr(neq,neq1,nr)   &
                                      *gm2p
           
             ! at the center of element
             r2imtx(1,1,neq,neq1)= - coef2*dh1*vtr(neq,neq1,nr) &
                                              *(2.d0*gm1m +      gm1p)
             r2imtx(2,1,neq,neq1)=   coef2*dh1*vtr(neq,neq1,nr) &
                                              *(2.d0*gm1m +      gm1p)
             r2imtx(1,2,neq,neq1)= - coef2*dh1*vtr(neq,neq1,nr) &
                                              *(     gm1m + 2.d0*gm1p)
             r2imtx(2,2,neq,neq1)=   coef2*dh1*vtr(neq,neq1,nr) &
                                              *(     gm1m + 2.d0*gm1p)
           
             r3imtx(1,1,neq,neq1)= dh2*(3.D0*ctr(neq,neq1,nr-1)*dvdrm &
                                      +      ctr(neq,neq1,nr  )*dvdrp)
             r3imtx(2,1,neq,neq1)= dh2*(     ctr(neq,neq1,nr-1)*dvdrm &
                                      +      ctr(neq,neq1,nr  )*dvdrp)
             r3imtx(1,2,neq,neq1)= dh2*(     ctr(neq,neq1,nr-1)*dvdrm &
                                      +      ctr(neq,neq1,nr  )*dvdrp)
             r3imtx(2,2,neq,neq1)= dh2*(     ctr(neq,neq1,nr-1)*dvdrm &
                                      + 3.D0*ctr(neq,neq1,nr  )*dvdrp)

             rimtx(1,1,neq,neq1)                           &
                =r1imtx(1,1,neq,neq1)-r2imtx(1,1,neq,neq1) &
                                     -r3imtx(1,1,neq,neq1)
             rimtx(2,1,neq,neq1)                           &
                =r1imtx(2,1,neq,neq1)-r2imtx(2,1,neq,neq1) &
                                     -r3imtx(2,1,neq,neq1)
             rimtx(1,2,neq,neq1)                           &
                =r1imtx(1,2,neq,neq1)-r2imtx(1,2,neq,neq1) &
                                     -r3imtx(1,2,neq,neq1)
             rimtx(2,2,neq,neq1)                           &
                =r1imtx(2,2,neq,neq1)-r2imtx(2,2,neq,neq1) &
                                     -r3imtx(2,2,neq,neq1)

          END DO ! End of NEQ1 loop

          CALL tr_bndry_set(neq,nr,nrmax)

          !--- calculate ELMTX with reduction ---

          IF(id_neq(neq) /= 0) THEN
             neqr=neqr_neq(neq)
             DO neq1=1,neqmax
                IF(id_neq(neq1) /= 0) THEN
                   neqr1=neqr_neq(neq1)
                   elmtx(neqr        ,neqr1        ) = rimtx(1,1,neq,neq1)
                   elmtx(neqr+neqrmax,neqr1        ) = rimtx(2,1,neq,neq1)
                   elmtx(neqr        ,neqr1+neqrmax) = rimtx(1,2,neq,neq1)
                   elmtx(neqr+neqrmax,neqr1+neqrmax) = rimtx(2,2,neq,neq1)
                END IF
             END DO
             elmtx(neqr        ,neqr        ) &
           = elmtx(neqr        ,neqr        ) + limtx(1,1,neq)

             elmtx(neqr+neqrmax,neqr        ) &
           = elmtx(neqr+neqrmax,neqr        ) + limtx(2,1,neq)

             elmtx(neqr        ,neqr+neqrmax) &
           = elmtx(neqr        ,neqr+neqrmax) + limtx(1,2,neq)

             elmtx(neqr+neqrmax,neqr+neqrmax) &
           = elmtx(neqr+neqrmax,neqr+neqrmax) + limtx(2,2,neq)
          END IF

       END DO ! End of NEQ loop

       !=== Build the LHS band matrix         

       DO neqr = 1, 2*neqrmax
          DO neqr1 = 1, 2*neqrmax
             lhmtx(2*neqrmax+neqr-neqr1,(nr-1)*neqrmax+neqr1) &
           = lhmtx(2*neqrmax+neqr-neqr1,(nr-1)*neqrmax+neqr1) &
                  + elmtx(neqr1,neqr)
          END DO
       END DO

       !=== Build the RHS vector

       nvm  = (nr-1)*neqmax
       nvp  =  nr   *neqmax
       nvrm = (nr-1)*neqrmax
       nvrp =  nr   *neqrmax
       DO neq = 1, neqmax
          IF(id_neq(neq) /= 0) THEN
             neqr = neqr_neq(neq)

             rhv(nvrm+neqr) = rhv(nvrm+neqr)          &
                  + limtx(1,1,neq) * xv_prev(nvm+neq) &
                  + limtx(1,2,neq) * xv_prev(nvp+neq) &
                  +rjimtx(1,1,neq) * htr(neq,nr-1)    &
                  +rjimtx(1,2,neq) * htr(neq,nr)      &
                  +rsimtx(1,1,neq) * str(neq,nr-1)    &
                  +rsimtx(1,2,neq) * str(neq,nr)

             rhv(nvrp+neqr) = rhv(nvrp+neqr)             &
                  + limtx(2,1,neq) * xv_prev(nvm+neq)    &
                  + limtx(2,2,neq) * xv_prev(nvp+neq)    &
                  +rjimtx(2,1,neq) * htr(neq,nr-1)       &
                  +rjimtx(2,2,neq) * htr(neq,nr)         &
                  +rsimtx(2,1,neq) * str(neq,nr-1)       &
                  +rsimtx(2,2,neq) * str(neq,nr)
          END IF
       END DO

    END DO ! END of NR loop

!    DO neq = 1, nvrmax
!       write(*,*) (lhmtx(nr,neq),nr=1,3)
!    END DO
!    DO neq = 1, nvrmax
!       write(*,*) rhv(neq)
!    END DO

    CALL BANDRD(lhmtx,rhv,nvrmax,4*neqrmax-1,4*neqrmax-1,ierr)
    IF(ierr /= 0) THEN
       WRITE(6,*) 'XX trexec: BANDRD error =',ierr
       STOP
    ENDIF

    DO neq=1,neqmax
       IF(id_neq(neq) /= 0) THEN
          neqr=neqr_neq(neq)
          DO nr=0,nrmax
             xv_new(nr*neqmax+neq) = rhv(nr*neqrmax+neqr)
          END DO
       ELSE
          DO nr=0,nrmax
             xv_new(nr*neqmax+neq) = xv(nr*neqmax+neq)
          END DO
       END IF
    END DO

    RETURN
  END SUBROUTINE tr_exec

! ***************************************************************************

  SUBROUTINE tr_neq_mtrc(nr,neq,dvdrp,dvdrm,ldrp,ldrm,rdrp,rdrm, &
                         gm1p,gm1m,gm2p,gm2m,cjexp,cjexm,coef1,coef2)
! ---------------------------------------------------------------------------
!   This subroutine decides coefficients of element equations of FEM.
! ---------------------------------------------------------------------------
    USE trcomm, ONLY: ar1rho,ar2rho,arrho,dvrho,abrho,ttrho,abb1rho,eta

    IMPLICIT NONE
    REAL(rkind),INTENT(OUT) :: &
         dvdrp,dvdrm,ldrp,ldrm,rdrp,rdrm,       &
         gm1p,gm1m,gm2p,gm2m,cjexp,cjexm,coef1,coef2

    INTEGER(ikind),INTENT(IN) :: nr,neq

    rdrm  = 0.d0
    rdrp  = 0.d0
    
    dvdrp = dvrho(nr)
    dvdrm = dvrho(nr-1)

    gm1p = dvdrp*ar1rho(nr) 
    gm1m = dvdrm*ar1rho(nr-1)
    
    IF(neq == 1)THEN ! magnetic field: d psi/d rho
       ldrp  = 1.d0   
       ldrm  = 1.d0
       coef1 = 1.d0
       coef2 = 1.d0
       gm2p  = dvdrp * abrho(nr  ) / ttrho(nr  ) ! C1 i+1
       gm2m  = dvdrm * abrho(nr-1) / ttrho(nr-1) ! C1 i

       eta = 0.d0
       cjexp = eta(nr  )*abb1rho(nr  )*arrho(nr  )/ttrho(nr  )
       cjexm = eta(nr-1)*abb1rho(nr-1)*arrho(nr-1)/ttrho(nr-1)
    ELSE IF(MOD((neq-1),3) == 2)THEN ! density
       ldrp  = dvrho(nr)
       ldrm  = dvrho(nr-1)
       coef1 = 1.d0
       coef2 = 1.d0
       gm2p  = 0.5d0*(dvdrp*ar2rho(nr) + dvdrm*ar2rho(nr-1))
       gm2m  = gm2p
    ELSE IF(MOD((neq-1),3) == 1)THEN ! toroidal velocity
       ldrp  = dvrho(nr)
       ldrm  = dvrho(nr-1)
       coef1 = 1.d0
       coef2 = 1.d0
       gm2p  = 0.5d0*(dvdrp*ar2rho(nr) + dvdrm*ar2rho(nr-1))
       gm2m  = gm2p
    ELSE IF(MOD((neq-1),3) == 0)THEN ! temperature
!       ldrp  = 3.d0/2.d0 * dvrho(nr  )**(5.d0/3.d0)
!       ldrm  = 3.d0/2.d0 * dvrho(nr-1)**(5.d0/3.d0)
       ldrp  = dvrho(nr)
       ldrm  = dvrho(nr-1)
!       rdrp  = dvrho(nr  )**(-2.d0/3.d0)
!       rdrm  = dvrho(nr-1)**(-2.d0/3.d0)
!       coef1 = 3.d0/2.d0 ! the coef. of LHS
!       coef2 = 3.d0/2.d0 ! the coef. of convection velocity
       coef1 = 1.d0
       coef2 = 1.d0
       gm2p  = 0.5d0*(dvdrp*ar2rho(nr) + dvdrm*ar2rho(nr-1))
       gm2m  = gm2p
    END IF

  END SUBROUTINE tr_neq_mtrc


  SUBROUTINE tr_bndry_set(neq,nr,nrmax)
! ---------------------------------------------------------------------------
!      This subroutine modifies element matrix 
!                             on the basis of the boundary conditions.
! ---------------------------------------------------------------------------
    USE trcomm, ONLY: id_neqnr,neqmax, limtx,rjimtx,rsimtx,rimtx
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN) :: neq,nr,nrmax
    INTEGER(ikind) :: neq1

    !--- fixed value ---
    IF(id_neqnr(neq,nr-1)==2)THEN
       IF(nr==1) limtx(1,1,neq) = 1.D0
       IF(nr/=1) limtx(1,1,neq) = 0.5D0
       limtx(1,2,neq) = 0.D0
       
       rjimtx(1,1,neq) = 0.D0
       rjimtx(1,2,neq) = 0.D0           
       
       rsimtx(1,1,neq) = 0.D0
       rsimtx(1,2,neq) = 0.D0
       
       DO neq1=1,neqmax
          rimtx(1,1,neq,neq1) = 0.D0
          rimtx(1,2,neq,neq1) = 0.D0
       END DO
       RETURN
    END IF    

    IF(id_neqnr(neq,nr)==2)THEN
       limtx(2,1,neq) = 0.D0
       IF(nr/=nrmax) limtx(2,2,neq) = 0.5D0
       IF(nr==nrmax) limtx(2,2,neq) = 1.D0
       
       rjimtx(2,1,neq) = 0.D0
       rjimtx(2,2,neq) = 0.D0
       
       rsimtx(2,1,neq) = 0.D0
       rsimtx(2,2,neq) = 0.D0
       
       DO neq1=1,neqmax
          rimtx(2,1,neq,neq1) = 0.D0
          rimtx(2,2,neq,neq1) = 0.D0
       END DO
    END IF

    RETURN
  END SUBROUTINE tr_bndry_set

END MODULE trexec
