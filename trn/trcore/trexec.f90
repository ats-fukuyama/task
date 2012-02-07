MODULE trexec

! This module solves the matrix equation

  PUBLIC tr_exec
  PRIVATE

CONTAINS

  SUBROUTINE tr_exec

    USE trcomm,ONLY: &
         rkind,ikind,pi,rr,rkap,neqmax,neqrmax,nvmax,nvrmax,nrmax,dt,rhog, &
         dtr,vtr,ctr,str,htr,elmtx,limtx,rsimtx,rjimtx,r1imtx,r2imtx,r3imtx, &
         rimtx,lhmtx,rhv,xv,xv_new,xv_prev,neqr_neq,id_neq,id_neqnr, &
         dvrho,ar1rho,ar2rho
    IMPLICIT NONE
    REAL(rkind) :: dh0,dh1,dh2,dh3,dvdrp,dvdrm,dvdr0,  &
                   gm1p,gm1m,gm10,gm2p,gm2m,gm20
    INTEGER(ikind) :: nr,neq,neq1,neqr,neqr1,nvrm,nvrp,ierr,nvm,nvp

    lhmtx(1:4*neqrmax-1,1:nvrmax) = 0.D0
    rhv(1:nvrmax)                 = 0.D0

    DO nr = 1, nrmax

       dvdrp = dvrho(nr)
       dvdrm = dvrho(nr-1)
       dvdr0 = 0.5D0*(dvdrm+dvdrp)

       gm1p = dvdrp*ar1rho(nr) 
       gm1m = dvdrm*ar1rho(nr-1)
       gm10 = 0.5d0*(gm1p+gm1m)

       gm2p = dvdrp*ar2rho(nr)
       gm2m = dvdrm*ar2rho(nr-1)
       gm20 = 0.5d0*(gm2p+gm2m)

       dh0 = (rhog(nr)-rhog(nr-1))/12.D0
       dh1 = dt/6.D0
       dh2 = dt*dh0
       dh3 = dt*gm20/(rhog(nr)-rhog(nr-1))

       elmtx(1:2*neqrmax,1:2*neqrmax)=0.D0

       DO neq = 1, neqmax

          limtx(1,1,neq) = dh0*(3.D0*dvdrm +      dvdrp)
          limtx(2,1,neq) = dh0*(     dvdrm +      dvdrp)
          limtx(1,2,neq) = dh0*(     dvdrm +      dvdrp)
          limtx(2,2,neq) = dh0*(     dvdrm + 3.D0*dvdrp) 

          rjimtx(1,1,neq) = dh1*(-4.D0*dvdrm +      dvdrp)
          rjimtx(2,1,neq) = dh1*(-2.D0*dvdrm -      dvdrp)
          rjimtx(1,2,neq) = dh1*(      dvdrm + 2.D0*dvdrp)
          rjimtx(2,2,neq) = dh1*(-     dvdrm + 4.D0*dvdrp)

          rsimtx(1,1,neq) = dh2*(3.D0*dvdrm +      dvdrp)
          rsimtx(2,1,neq) = dh2*(     dvdrm +      dvdrp)
          rsimtx(1,2,neq) = dh2*(     dvdrm +      dvdrp)
          rsimtx(2,2,neq) = dh2*(     dvdrm + 3.D0*dvdrp) 

          DO neq1 = 1, neqmax
        
             ! at the center of element
             r1imtx(1,1,neq,neq1)=   dh3*dtr(neq,neq1,nr)
             r1imtx(2,1,neq,neq1)= - dh3*dtr(neq,neq1,nr)
             r1imtx(1,2,neq,neq1)= - dh3*dtr(neq,neq1,nr)
             r1imtx(2,2,neq,neq1)=   dh3*dtr(neq,neq1,nr)
           
             ! at the center of element
             r2imtx(1,1,neq,neq1)= - dh1*vtr(neq,neq1,nr) &
                                    *(2.d0*gm1m +      gm1p)
             r2imtx(2,1,neq,neq1)=   dh1*vtr(neq,neq1,nr) &
                                    *(2.d0*gm1m +      gm1p)
             r2imtx(1,2,neq,neq1)= - dh1*vtr(neq,neq1,nr) &
                                    *(     gm1m + 2.d0*gm1p)
             r2imtx(2,2,neq,neq1)=   dh1*vtr(neq,neq1,nr) &
                                    *(     gm1m + 2.d0*gm1p)
           
             r3imtx(1,1,neq,neq1)= dh2*(3.D0*ctr(neq,neq1,nr-1)*dvdrm &
                                      +      ctr(neq,neq1,nr  )*dvdrp)
             r3imtx(2,1,neq,neq1)= dh2*(     ctr(neq,neq1,nr-1)*dvdrm &
                                      +      ctr(neq,neq1,nr  )*dvdrp)
             r3imtx(1,2,neq,neq1)= dh2*(     ctr(neq,neq1,nr-1)*dvdrm &
                                      +      ctr(neq,neq1,nr  )*dvdrp)
             r3imtx(2,2,neq,neq1)= dh2*(     ctr(neq,neq1,nr-1)*dvdrm &
                                      + 3.D0*ctr(neq,neq1,nr  )*dvdrp)

             rimtx(1,1,neq,neq1) &
                =r1imtx(1,1,neq,neq1)-r2imtx(1,1,neq,neq1)-r3imtx(1,1,neq,neq1)
             rimtx(2,1,neq,neq1) &
                =r1imtx(2,1,neq,neq1)-r2imtx(2,1,neq,neq1)-r3imtx(2,1,neq,neq1)
             rimtx(1,2,neq,neq1) &
                =r1imtx(1,2,neq,neq1)-r2imtx(1,2,neq,neq1)-r3imtx(1,2,neq,neq1)
             rimtx(2,2,neq,neq1) &
                =r1imtx(2,2,neq,neq1)-r2imtx(2,2,neq,neq1)-r3imtx(2,2,neq,neq1)
          END DO
        
          !--- fixed value ---

          IF(id_neqnr(neq,nr-1) == 2)THEN
             limtx(1,1,neq) = 1.D0
             limtx(1,2,neq) = 0.D0

             rjimtx(1,1,neq) = 0.D0
             rjimtx(1,2,neq) = 0.D0           

             rsimtx(1,1,neq) = 0.D0
             rsimtx(1,2,neq) = 0.D0

             DO neq1=1,neqMAX
                rimtx(1,1,neq,neq1) = 0.D0
                rimtx(1,2,neq,neq1) = 0.D0
             END DO
          END IF

          IF(id_neqnr(neq,nr) == 2)THEN
             limtx(2,1,neq) = 0.D0
             limtx(2,2,neq) = 1.D0

             rjimtx(2,1,neq) = 0.D0
             rjimtx(2,2,neq) = 0.D0

             rsimtx(2,1,neq) = 0.D0
             rsimtx(2,2,neq) = 0.D0

             DO neq1=1,neqMAX
                rimtx(2,1,neq,neq1) = 0.D0
                rimtx(2,2,neq,neq1) = 0.D0
             END DO
          END IF

          !--- calculate ELMTX with reduction ---

          IF(id_neq(neq) == 1) THEN
             neqr=neqr_neq(neq)
             DO neq1=1,neqmax
                IF(id_neq(neq1) == 1) THEN
                   neqr1=neqr_neq(neq1)
                   elmtx(neqr        ,neqr1        ) = rimtx(1,1,neq,neq1)
                   elmtx(neqr+neqrmax,neqr1        ) = rimtx(2,1,neq,neq1)
                   elmtx(neqr        ,neqr1+neqrmax) = rimtx(1,2,neq,neq1)
                   elmtx(neqr+neqrmax,neqr1+neqrmax) = rimtx(2,2,neq,neq1)
                END IF
             END DO
             elmtx(neqr        ,neqr        ) &
           = elmtx(neqr        ,neqr        )&
                  + limtx(1,1,neq)
             elmtx(neqr+neqrmax,neqr        ) &
           = elmtx(neqr+neqrmax,neqr        )&
                  + limtx(2,1,neq)
             elmtx(neqr        ,neqr+neqrmax) &
           = elmtx(neqr        ,neqr+neqrmax)&
                  + limtx(1,2,neq)
             elmtx(neqr+neqrmax,neqr+neqrmax) &
           = elmtx(neqr+neqrmax,neqr+neqrmax)&
                  + limtx(2,2,neq)
          END IF
       END DO

       !=== BUILD THE LEFT HAND BAND MATRIX         

       DO neqr = 1, 2*neqrmax
          DO neqr1 = 1, 2*neqrmax
             lhmtx(2*neqrmax+neqr-neqr1,(nr-1)*neqrmax+neqr1) &
           = lhmtx(2*neqrmax+neqr-neqr1,(nr-1)*neqrmax+neqr1) &
                  + elmtx(neqr1,neqr)
          END DO
       END DO

       !=== BUILD THE RHS VECTOR

       NVM=(nr-1)*neqmax
       NVP= nr   *neqmax
       NVRM=(nr-1)*neqrmax
       NVRP= nr   *neqrmax
       DO neq = 1, neqmax
          IF(id_neq(neq) == 1) THEN
             neqr=neqr_neq(neq)
             RHV(NVRM+neqr) = RHV(NVRM+neqr)          &
                  + limtx(1,1,neq) * xv_prev(NVM+neq) &
                  + limtx(1,2,neq) * xv_prev(NVP+neq) &
                  +rjimtx(1,1,neq) * htr(neq,nr-1)    &
                  +rjimtx(1,2,neq) * htr(neq,nr)      &
                  +rsimtx(1,1,neq) * str(neq,nr-1)    &
                  +rsimtx(1,2,neq) * str(neq,nr)

             RHV(NVRP+neqr) = RHV(NVRP+neqr)             &
                  + limtx(2,1,neq) * xv_prev(NVM+neq)    &
                  + limtx(2,2,neq) * xv_prev(NVP+neq)    &
                  +rjimtx(2,1,neq) * htr(neq,nr-1)       &
                  +rjimtx(2,2,neq) * htr(neq,nr)         &
                  +rsimtx(2,1,neq) * str(neq,nr-1)       &
                  +rsimtx(2,2,neq) * str(neq,nr)
          END IF
       END DO

    END DO ! nr LOOP END

    CALL BANDRD(lhmtx,rhv,nvrmax,4*neqrmax-1,4*neqrmax-1,ierr)
    IF(ierr /= 0) THEN
       WRITE(6,*) 'XX trexec: BANDRD error =',ierr
       STOP
    ENDIF

    DO neq=1,neqmax
       IF(id_neq(neq) == 1) THEN
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
END MODULE trexec
