!    $Id$

MODULE wiexec
  PRIVATE
  PUBLIC wi_exec
 
CONTAINS

  SUBROUTINE wi_exec(iprint,ratea,ierr)

    USE wicomm,ONLY: qkind,ikind,modelg,modelp
    USE wicold,ONLY: wi_cold
    USE wiwarm,ONLY: wi_warm
    USE wihot,ONLY: wi_hot
    
    IMPLICIT NONE
    INTEGER(ikind),INTENT(IN):: iprint
    REAL(qkind),INTENT(OUT):: ratea
    INTEGER,INTENT(OUT):: ierr

    SELECT CASE(modelg)
    CASE(0)
       SELECT CASE(modelp)
       CASE(0)
          CALL wi_cold(iprint,ratea,ierr)
       CASE(1)
          CALL wi_warm(iprint,ratea,ierr)
       CASE(2)
          CALL wi_hot(iprint,ratea,ierr)
       END SELECT
    END SELECT

    RETURN
  END SUBROUTINE wi_exec
END MODULE wiexec
