MODULE trstep

! This module advances one time step.
! Full implicit method in time with iteration
! Finite element method in radial direction

  PUBLIC tr_step
  PRIVATE

CONTAINS

  SUBROUTINE tr_step(ierr)

    USE trcomm, ONLY: rkind,ikind,nrmax,neqmax,xv,xv_prev,xv_new,qp,rn,ru,rt, &
         lmaxtr,epsltr,nsa_neq,nva_neq,nvmax,error_it,nitmax,mdltr_prv
    USE trcoef, ONLY: Pereverzev_check
    USE trcalc, ONLY: tr_calc
    USE trexec, ONLY: tr_exec

    IMPLICIT NONE
    INTEGER(ikind),INTENT(OUT):: ierr
    INTEGER(ikind):: nr,neq,nv,nsa,nva,nit
    REAL(rkind):: dif(neqmax),ave(neqmax),difmax

    DO nr=0,nrmax
       DO neq=1,neqmax
          nv=nr*neqmax+neq
          IF(nsa_neq(neq) == 0) THEN
             xv(nv)=qp(nr)
          ELSE
             nsa=nsa_neq(neq)
             nva=nva_neq(neq)
             SELECT CASE(nva)
             CASE(1)
                xv(nv)=rn(nsa,nr)
             CASE(2)
                xv(nv)=ru(nsa,nr)
             CASE(3)
                xv(nv)=rt(nsa,nr)
             END SELECT
          END IF
       END DO
    END DO

    xv_prev(1:nvmax)=xv(1:nvmax)
    error_it(1:lmaxtr) = epsltr
       
    DO NIT = 1, lmaxtr

       CALL tr_calc

       CALL tr_exec

       !--- CONVERGENCE CHECK ---
       DO neq=1,neqmax
          dif(neq)=0.D0
          ave(neq)=0.d0
          DO nr=0,nrmax
             nv=nr*neqmax+neq
             dif(neq)=dif(neq)+(xv_new(nv)-xv(nv))**2
             ave(neq)=ave(neq)+xv_new(nv)**2
          END DO
          ave(neq)=SQRT(ave(neq))
          dif(neq)=SQRT(dif(neq))
          IF(ave(neq) > 0.D0) THEN
             dif(neq)=dif(neq)/ave(neq)
          ENDIF
       END DO
       difmax=MAXVAL(dif)
       error_it(nit) = difmax

       xv(1:nvmax)=xv_new(1:nvmax)

       DO nr=0,nrmax
          DO neq=1,neqmax
             nv=nr*neqmax+neq
             IF(nsa_neq(neq) == 0) THEN
                qp(nr)=xv(nv)
             ELSE
                nsa=nsa_neq(neq)
                nva=nva_neq(neq)
                SELECT CASE(nva)
                CASE(1)
                   rn(nsa,nr)=xv(nv)
                CASE(2)
                   ru(nsa,nr)=xv(nv)
                CASE(3)
                   rt(nsa,nr)=xv(nv)
                END SELECT
             END IF
          END DO
       END DO

       IF(difmax < epsltr) THEN
          nitmax=MAX(nit,nitmax)
          GO TO 700
       END IF
    END DO

    nitmax=lmaxtr+1

700 CONTINUE

    DO nr=0,nrmax
       DO neq=1,neqmax
          nv=nr*neqmax+neq
          IF(nsa_neq(neq) == 0) THEN
             qp(nr)=xv(nv)
          ELSE
             nsa=nsa_neq(neq)
             nva=nva_neq(neq)
             SELECT CASE(nva)
             CASE(1)
                rn(nsa,nr)=xv(nv)
             CASE(2)
                ru(nsa,nr)=xv(nv)
             CASE(3)
                rt(nsa,nr)=xv(nv)
             END SELECT
          END IF
       END DO
    END DO

    IF(mdltr_prv /= 0)THEN
       CALL Pereverzev_check
!       WRITE(*,*)
    END IF
       
!   --- error check here ---
    ierr=0

    RETURN
  END SUBROUTINE tr_step
END MODULE trstep
