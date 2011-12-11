MODULE trsimple

  USE bpsd_kinds

  PRIVATE
  PUBLIC tr_simple

CONTAINS

! ***** calculate simple transport coefficients *****

  SUBROUTINE tr_simple

    USE trcomm, ONLY: nrmax,nsamax,mdltr_tb,ltcr,d0,d1,rg,rt,dtr_tb,vtr_tb, &
         lt_save
    IMPLICIT NONE
    INTEGER(ikind) :: NR, NEQ, nsa
    REAL(rkind) :: LT
    REAL(rkind),DIMENSION(nsamax,nrmax):: DTR_DIAG

    dtr_tb(1:3*nsamax,1:3*nsamax,1:nrmax)=0.D0
    vtr_tb(1:3*nsamax,1:3*nsamax,1:nrmax)=0.D0

    SELECT CASE(mdltr_tb)
       ! --- No transport
    CASE(0)

       ! --- FLAT PROFILE ---
    CASE(1)
       DO nr = 1, nrmax
          DO nsa = 1, nsamax
             dtr_diag(nsa,nr) = d0
             if(nsa /= 0) lt_save(nsa,nr)=0.d0
          END DO
       END DO
      
       ! --- Typical stiff model (Pereverzev) ---
    CASE(2)
       DO nr = 1, nrmax
          DO nsa = 1, nsamax
             IF(rt(nsa,nr)-rt(nsa,nr-1) > 0) THEN
                lt = 0.D0
             ELSE
                lt = - (rt(nsa,nr)-rt(nsa,nr-1))/(rg(nr)-rg(nr-1))
             END IF
             lt_save(nr,nsa)=lt
             IF(lt > ltcr)THEN
                dtr_diag(nsa,nr) = d0 + d1 - d1*ltcr/lt
             else
                dtr_diag(nsa,nr) = d0
             END IF
          END DO
       END DO
    END SELECT
         
    DO nr = 1, nrmax
       DO nsa = 1, nsamax
          dtr_tb(3*nsa-2,3*nsa-2,nr) = dtr_diag(nsa,nr)
          dtr_tb(3*nsa-1,3*nsa-1,nr) = dtr_diag(nsa,nr)
          dtr_tb(3*nsa  ,3*nsa  ,nr) = dtr_diag(nsa,nr)
       END DO
    END DO
    RETURN
  END SUBROUTINE tr_simple

END MODULE trsimple
