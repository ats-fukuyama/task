MODULE trsimple

  USE bpsd_kinds

  PRIVATE
  PUBLIC tr_simple

CONTAINS

! ***** calculate simple transport coefficients *****

  SUBROUTINE tr_simple

    USE trcomm, ONLY: nrmax,nsamax,neqmax,mdltr_tb,rg,rm,rhom,rt, &
         ltcr,dtr0,dtr1,dtr_tb,vtr_tb,cdtrn,cdtru,cdtrt
    IMPLICIT NONE
    INTEGER(ikind) :: NR, NEQ, nsa
    REAL(rkind) :: LT
    REAL(rkind),DIMENSION(nsamax,nrmax):: DTR_DIAG

    dtr_tb(1:neqmax,1:neqmax,1:nrmax)=0.D0
    vtr_tb(1:neqmax,1:neqmax,1:nrmax)=0.D0

    dtr_diag(1:nsamax,1:nrmax) = 0.d0

    SELECT CASE(mdltr_tb)
    ! --- No transport
    CASE(0)
       RETURN
    ! --- Flat profile ---
    CASE(1)
       DO nr = 1, nrmax
          DO nsa = 1, nsamax
             dtr_diag(nsa,nr) = dtr0
          END DO
       END DO

    ! ---   ---
    CASE(2)
       DO nr = 1, nrmax
          DO nsa = 1, nsamax
             dtr_diag(nsa,nr) = dtr0 + dtr1*rhom(nr)**2
          END DO
       END DO

    ! --- Typical stiff model ---
    CASE(3)
       DO nr = 1, nrmax
          DO nsa = 1, nsamax
             IF(rt(nsa,nr)-rt(nsa,nr-1) > 0) THEN
                lt = 0.D0
             ELSE
                lt = - (rt(nsa,nr)-rt(nsa,nr-1))/(rg(nr)-rg(nr-1))
             END IF
             IF(lt > ltcr)THEN
                dtr_diag(nsa,nr) = dtr0 + dtr1 - dtr1*ltcr/lt
             else
                dtr_diag(nsa,nr) = dtr0
             END IF
          END DO
       END DO
    END SELECT
         
    DO nr = 1, nrmax
       DO nsa = 1, nsamax
          dtr_tb(1+3*nsa-2,1+3*nsa-2,nr) = cdtrn*dtr_diag(nsa,nr)
          dtr_tb(1+3*nsa-1,1+3*nsa-1,nr) = cdtru*dtr_diag(nsa,nr)
          dtr_tb(1+3*nsa  ,1+3*nsa  ,nr) = cdtrt*dtr_diag(nsa,nr)
       END DO
    END DO
    RETURN
  END SUBROUTINE tr_simple

END MODULE trsimple
