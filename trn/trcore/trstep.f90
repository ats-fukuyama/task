MODULE trstep
! -------------------------------------------------------------------------
!  This module advances one time step.
!   - Full implicit method in time with iteration
!   - Finite element method in radial direction
! -------------------------------------------------------------------------
  USE trcomm, ONLY: ikind,rkind
  IMPLICIT NONE

  PUBLIC tr_step
  PRIVATE

CONTAINS

  SUBROUTINE tr_step(ierr)

    USE trcomm, ONLY: nrmax,neqmax,nvmax,xv,xv_prev,xv_new,       &
         rn,ru,rt,rp,rp_tot,dpdrho,lmaxtr,epsltr,nsa_neq,nva_neq, &
         error_it,nitmax,mdltr_prv,dtr_tb  ! ,nrd1
    USE trcalc2, ONLY: tr_calc2
    USE trexec, ONLY: tr_exec

    IMPLICIT NONE
    INTEGER(ikind),INTENT(OUT):: ierr
    INTEGER(ikind):: nr,neq,nv,nit
    REAL(rkind):: dif(neqmax),ave(neqmax),difmax

    ierr = 0

    ! xv <-- dpdrho,rn,ru,rt
    CALL tr_get_xv(xv,dpdrho,rn,ru,rt)

    xv_prev(1:nvmax)   = xv(1:nvmax)
    error_it(1:lmaxtr) = epsltr
       
    non_linear: DO nit = 1, lmaxtr

       CALL tr_calc2

       CALL tr_exec

       !--- CONVERGENCE CHECK ---
       DO neq = 1, neqmax
          dif(neq) = 0.d0
          ave(neq) = 0.d0
          DO nr = 0, nrmax
             nv = nr*neqmax + neq
!             write(6,*) nr, xv_new(nv), xv(nv)
             dif(neq) = dif(neq) + (xv_new(nv)-xv(nv))**2
             ave(neq) = ave(neq) + xv_new(nv)**2
          END DO
          ave(neq) = SQRT(ave(neq))
          dif(neq) = SQRT(dif(neq))
          IF(ave(neq) > 0.D0) THEN
             dif(neq) = dif(neq)/ave(neq)
          ENDIF
       END DO

       difmax = MAXVAL(dif)
!       write(*,*) difmax
       error_it(nit) = difmax
       xv(1:nvmax) = xv_new(1:nvmax)

       IF(difmax < epsltr) THEN
          nitmax=MAX(nit,nitmax)
!          nitmax = nit
          EXIT
       END IF

       ! dpdrho,rn,ru,rt <-- xv
       CALL tr_set_xv(xv,dpdrho,rn,ru,rt,rp,rp_tot)

       CALL tr_check_negative(ierr)
       IF(ierr /= 0) EXIT

       CALL tr_calc_bpqpj

    END DO non_linear

    IF(nit==lmaxtr+1) nitmax = lmaxtr + 1

    CALL tr_set_xv(xv,dpdrho,rn,ru,rt,rp,rp_tot)
!    nrd1(0:nrmax) = rt(1,0:nrmax)*rn(1,0:nrmax)


    RETURN
  END SUBROUTINE tr_step

! *************************************************************************

  SUBROUTINE tr_get_xv(xv,dpdrho,rn,ru,rt)
    USE trcomm, ONLY:nsamax,nrmax,neqmax,nvmax,nsa_neq,nva_neq
    
    REAL(rkind),DIMENSION(nsamax,0:nrmax),INTENT(IN) :: rn,ru,rt
    REAL(rkind),DIMENSION(0:nrmax),INTENT(IN) :: dpdrho
    REAL(rkind),DIMENSION(1:nvmax),INTENT(OUT) :: xv

    INTEGER(ikind) :: nr,neq,nv,nsa,nva

    DO nr = 0, nrmax
       DO neq = 1, neqmax
          nv = nr*neqmax + neq
          IF(nsa_neq(neq) == 0) THEN
             xv(nv)=dpdrho(nr)
          ELSE
             nsa=nsa_neq(neq)
             nva=nva_neq(neq)
             SELECT CASE(nva)
             CASE(1)
                xv(nv)=rn(nsa,nr)
             CASE(2)
                xv(nv)=ru(nsa,nr)
             CASE(3)
                xv(nv)=rn(nsa,nr)*rt(nsa,nr)
!                xv(nv)=rt(nsa,nr)
             END SELECT
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE tr_get_xv


  SUBROUTINE tr_set_xv(xv,dpdrho,rn,ru,rt,rp,rp_tot)
    USE trcomm,ONLY:rkev,nsamax,nrmax,neqmax,nvmax,nsa_neq,nva_neq,id_neq
    
    REAL(rkind),DIMENSION(nsamax,0:nrmax),INTENT(OUT) :: rn,ru,rt,rp
    REAL(rkind),DIMENSION(0:nrmax),INTENT(OUT) :: dpdrho,rp_tot
    REAL(rkind),DIMENSION(1:nvmax),INTENT(IN) :: xv

    INTEGER(ikind) :: nr,neq,nv,nsa,nva

    DO nr = 0, nrmax
       DO neq = 1, neqmax
          nv = nr*neqmax + neq
          IF(nsa_neq(neq) == 0) THEN
             dpdrho(nr) = xv(nv)
          ELSE
             nsa = nsa_neq(neq)
             nva = nva_neq(neq)
             SELECT CASE(nva)
             CASE(1)
                rn(nsa,nr) = xv(nv)
             CASE(2)
                ru(nsa,nr) = xv(nv)
             CASE(3)
                IF(id_neq(neq) /= 0)THEN
                   rt(nsa,nr) = xv(nv)/rn(nsa,nr)
                END IF
             END SELECT
          END IF
       END DO
    END DO

    rp_tot(0:nrmax) = 0.d0
    DO nr = 0, nrmax
       DO nsa = 1, nsamax
          ! the pressure of each species
          rp(nsa,nr) = rn(nsa,nr)*1.d20 * rt(nsa,nr)*rkev
          rp_tot(nr) = rp_tot(nr) + rp(nsa,nr)
       END DO
    END DO

    RETURN
  END SUBROUTINE tr_set_xv


  SUBROUTINE tr_check_negative(ierr)
! ------------------------------------------------------------------------
! check temperature and density profiles whether they are negative or not
! ------------------------------------------------------------------------
    USE trcomm,ONLY: nrmax,nsamax,t,dt,rn,rt,eta
    IMPLICIT NONE

    INTEGER(ikind),INTENT(OUT) :: ierr

    INTEGER(ikind)    :: nsa,nr
    CHARACTER(LEN=32) :: fmt_1, fmt_2, fmt_3

    ierr = 0
    fmt_1 = '(1X,A,I8,A)'
    fmt_2 = '(1X,A10,I2,A5,I3,A3,ES12.4)'
    fmt_3 = '(1X,A11,I3,A3,ES12.4)'

    DO nsa = 1, nsamax
       DO nr = 0, nrmax
          IF(rn(nsa,nr) < 0.d0)THEN
             WRITE(6,fmt_1) &
       'XX tr_negative_check: negative density value at step (',INT(t/dt),' )'
             WRITE(6,fmt_2) 'XX n (nsa=',nsa,', nr=',nr,') =', rt(nsa,nr)
             ierr = 1
          END IF
          IF(rt(nsa,nr) < 0.d0)THEN
             WRITE(6,fmt_1) &
   'XX tr_negative_check: negative temperature value at step (',INT(t/dt),' )'
             WRITE(6,fmt_2) 'XX T (nsa=',nsa,', nr=',nr,') =', rt(nsa,nr)
             ierr = 1
          END IF
          
          IF(ierr /= 0) RETURN
       END DO
    END DO

!!$    DO nr = 0, nrmax
!!$       IF(eta(nr) < 0.d0)THEN
!!$          WRITE(6,fmt_1) &
!!$       'XX tr_negative_check: negative eta value at step (',INT(t/dt),' )'
!!$          WRITE(6,fmt_3) 'XX eta (nr=',nr,') =',eta(nr)
!!$          ierr = 1
!!$       END IF
!!$    END DO

    RETURN
  END SUBROUTINE tr_check_negative

! **************************************************************************

  SUBROUTINE tr_calc_bpqpj
! --------------------------------------------------------------------------
!  calculate following conversions:  d psi/d rho --> bp,qp,jtor,jtot
! --------------------------------------------------------------------------
    USE trcomm, ONLY: pi,rmu0,nrmax,RR,ar1rho,ttrho,rmjrho,arrho,dvrho,  &
         abb1rho,abrho,abvrho,rhog,dpdrho,rdpvrho,qp,q0,qa,bp,rip,       &
         jtot,joh,jtor,jbs_nc,jex_nc,jcd_nb,jcd_ec,jcd_lh,jcd_ic,        &
         id_neq,mdlijq
!         ,nrd3,nrd4

    IMPLICIT NONE
    INTEGER(ikind) :: nr, mdlid
    REAL(rkind) :: FCTR    ! function defined in TASK/lib
!    REAL(rkind) :: FCTR4pt ! function defined in TASK/lib
    REAL(rkind) :: deriv3  ! function defined in TASK/lib
!    REAL(rkind) :: ipl,dr
    REAL(rkind),DIMENSION(0:nrmax) :: factor1,factor2

    factor1(0:nrmax) = dvrho(0:nrmax)*abrho(0:nrmax)*dpdrho(0:nrmax)
    factor2(0:nrmax) = factor1(0:nrmax)/ttrho(0:nrmax)

    ! dpdrho --> bp
    bp(0:nrmax) = ar1rho(0:nrmax)*dpdrho(0:nrmax)/rmjrho(0:nrmax)

    IF(id_neq(1) == 0)THEN

       mdlid = MOD(mdlijq,2)
       SELECT CASE(mdlid)
       CASE(0) ! in the case (qp --> dpdrho)
          DO nr = 1, nrmax
             ! dpdrho --> jtot(j_para)
             jtot(nr) = ttrho(nr)**2/(rmu0*abb1rho(nr)*dvrho(nr)) &
                        * deriv3(nr,rhog,factor2,nrmax,0)
          END DO
!          jtot(0) = FCTR4pt(rhog(1),rhog(2),rhog(3),jtot(1),jtot(2),jtot(3))
          jtot(0) = FCTR(rhog(1),rhog(2),jtot(1),jtot(2))
          
       CASE(1) ! in the case (jtot --> dpdrho)
          ! dpdrho --> qp
          qp(1:nrmax) = ttrho(1:nrmax)*arrho(1:nrmax)*dvrho(1:nrmax)    &
                        /(4.d0*pi**2 * dpdrho(1:nrmax))
!          qp(0)       = FCTR4pt(rhog(1),rhog(2),rhog(3),qp(1),qp(2),qp(3))
          qp(0)       = FCTR(rhog(1),rhog(2),qp(1),qp(2))
!          qp(0) = qp(1)
          q0 = qp(0)
          qa = qp(nrmax)

       END SELECT

    ELSE IF(id_neq(1) == 2)THEN
       ! dpdrho --> qp
       qp(1:nrmax) = ttrho(1:nrmax)*arrho(1:nrmax)*dvrho(1:nrmax)    &
                     /(4.d0*pi**2 * dpdrho(1:nrmax))
!       qp(0)       = FCTR4pt(rhog(1),rhog(2),rhog(3),qp(1),qp(2),qp(3))
       qp(0)       = FCTR(rhog(1),rhog(2),qp(1),qp(2))
!       qp(0) = qp(1)
       q0 = qp(0)
       qa = qp(nrmax)

       DO nr = 1, nrmax
          ! dpdrho --> jtot(j_para)
          jtot(nr) = ttrho(nr)**2/(rmu0*abb1rho(nr)*dvrho(nr)) &
                     * deriv3(nr,rhog,factor2,nrmax,0)
       END DO
!       jtot(0) = FCTR4pt(rhog(1),rhog(2),rhog(3),jtot(1),jtot(2),jtot(3))
       jtot(0) = FCTR(rhog(1),rhog(2),jtot(1),jtot(2))

    END IF

    ! only for graphic output -----------------------------------------
    rdpvrho(0:nrmax) = ttrho(0:nrmax)*arrho(0:nrmax)          &
                      /(4.d0*pi**2*qp(0:nrmax)) ! d psi/d V 
    DO nr = 1, nrmax
          ! dpdrho --> jtor(j_toroidal)    
          jtor(nr) = rmjrho(nr)/(rmu0*dvrho(nr))               &
                     * deriv3(nr,rhog,factor1,nrmax,0)

    END DO
!    jtor(0) = FCTR4pt(rhog(1),rhog(2),rhog(3),jtor(1),jtor(2),jtor(3))
    jtor(0) = FCTR(rhog(1),rhog(2),jtor(1),jtor(2))
    ! -----------------------------------------------------------------


    joh(0:nrmax) = jtot(0:nrmax)                &
         - (jbs_nc(0:nrmax) + jex_nc(0:nrmax) + &
            jcd_nb(0:nrmax) + jcd_ec(0:nrmax) + &
            jcd_ic(0:nrmax) + jcd_lh(0:nrmax))

    ! ***** inverse conversion for confirmation *****
!!$    ipl = 0.d0
!!$    DO nr = 1, nrmax
!!$       dr = rhog(nr) - rhog(nr-1)
!!$       ipl = ipl + 0.5d0*(dvrho(nr)+dvrho(nr-1))  &
!!$                  *0.5d0*(jtor(nr)+jtor(nr-1))*dr
!!$    END DO
!!$    ipl = ipl/(2.d0*pi*RR)
!!$    write(*,*) rip,ipl
    ! ***********************************************

!    nrd3(0:nrmax) = dpdrho(0:nrmax)
!    nrd4(0:nrmax) = jtot(0:nrmax)

    RETURN
  END SUBROUTINE tr_calc_bpqpj

END MODULE trstep
