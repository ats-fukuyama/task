MODULE trstep

  USE trcomm, ONLY: ikind,rkind
! This module advances one time step.
! Full implicit method in time with iteration
! Finite element method in radial direction

  PUBLIC tr_step
  PRIVATE

CONTAINS

  SUBROUTINE tr_step(ierr)

    USE trcomm, ONLY: nrmax,neqmax,xv,xv_prev,xv_new, &
         rn,ru,rt,dpdrho,lmaxtr,epsltr,nsa_neq,nva_neq,nvmax,error_it,&
         nitmax,mdltr_prv  ! ,nrd1
    USE trcalc, ONLY: tr_calc
    USE trexec, ONLY: tr_exec

    IMPLICIT NONE
    INTEGER(ikind),INTENT(OUT):: ierr
    INTEGER(ikind):: nr,neq,nv,nit
    REAL(rkind):: dif(neqmax),ave(neqmax),difmax

    ! xv <-- dpdrho,rn,ru,rt
    CALL tr_get_xv(xv,dpdrho,rn,ru,rt)

    xv_prev(1:nvmax)   = xv(1:nvmax)
    error_it(1:lmaxtr) = epsltr
       
    DO nit = 1, lmaxtr

!       write(*,*) 'in nonlinear iteration'

       CALL tr_calc

       CALL tr_exec

       !--- CONVERGENCE CHECK ---
       DO neq = 1, neqmax
          dif(neq) = 0.d0
          ave(neq) = 0.d0
          DO nr = 0, nrmax
             nv = nr*neqmax + neq
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
       error_it(nit) = difmax
!       write(*,*) difmax

       xv(1:nvmax) = xv_new(1:nvmax)

       ! dpdrho,rn,ru,rt <-- xv
       CALL tr_set_xv(xv,dpdrho,rn,ru,rt)

       IF(difmax < epsltr) THEN
          nitmax=MAX(nit,nitmax)
          GO TO 700
       END IF
    END DO

    nitmax=lmaxtr+1

700 CONTINUE

    CALL tr_set_xv(xv,dpdrho,rn,ru,rt)
!    nrd1(0:nrmax) = rt(1,0:nrmax)*rn(1,0:nrmax)
       
!   --- error check here ---
    ierr=0

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

  SUBROUTINE tr_set_xv(xv,dpdrho,rn,ru,rt)
    USE trcomm, ONLY:nsamax,nrmax,neqmax,nvmax,nsa_neq,nva_neq
    
    REAL(rkind),DIMENSION(nsamax,0:nrmax),INTENT(OUT) :: rn,ru,rt
    REAL(rkind),DIMENSION(0:nrmax),INTENT(OUT) :: dpdrho
    REAL(rkind),DIMENSION(1:nvmax),INTENT(IN) :: xv

    INTEGER(ikind) :: nr,neq,nv,nsa,nva

    DO nr=0,nrmax
       DO neq=1,neqmax
          nv=nr*neqmax+neq
          IF(nsa_neq(neq) == 0) THEN
             dpdrho(nr)=xv(nv)
          ELSE
             nsa=nsa_neq(neq)
             nva=nva_neq(neq)
             SELECT CASE(nva)
             CASE(1)
                rn(nsa,nr)=xv(nv)
             CASE(2)
                ru(nsa,nr)=xv(nv)
             CASE(3)
                rt(nsa,nr)=xv(nv)/rn(nsa,nr)
!                rt(nsa,nr)=xv(nv)
             END SELECT
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE tr_set_xv


END MODULE trstep
