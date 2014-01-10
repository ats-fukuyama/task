!    $Id$

MODULE wiexec
  PRIVATE
  PUBLIC wi_exec
 
CONTAINS

  SUBROUTINE wi_exec(ierr)

    USE wicomm,ONLY: modewi
    USE wiunmag,ONLY: wi_unmag
    
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    SELECT CASE(modewi)
    CASE(0,1)
       CALL wi_unmag(ierr)
    END SELECT

    RETURN
  END SUBROUTINE wi_exec
END MODULE wiexec
