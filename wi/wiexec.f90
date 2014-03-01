!    $Id$

MODULE wiexec
  PRIVATE
  PUBLIC wi_exec
 
CONTAINS

  SUBROUTINE wi_exec(iprint,ratea,ierr)

    USE wicomm,ONLY: rkind,ikind,modelg,modelp
    USE wiunmag,ONLY: wi_unmag
    USE wifluid,ONLY: wi_fluid
    
    IMPLICIT NONE
    INTEGER(ikind),INTENT(IN):: iprint
    REAL(rkind),INTENT(OUT):: ratea
    INTEGER,INTENT(OUT):: ierr

    SELECT CASE(modelg)
    CASE(0)
       SELECT CASE(modelp)
       CASE(0)
          CALL wi_fluid(iprint,ratea,ierr)
       CASE(1)
          CALL wi_unmag(iprint,ratea,ierr)
       END SELECT
    END SELECT

    RETURN
  END SUBROUTINE wi_exec
END MODULE wiexec
