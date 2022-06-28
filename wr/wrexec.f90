! wrexec.f90

MODULE wrexec

  PRIVATE
  PUBLIC wr_exec

CONTAINS

  SUBROUTINE wr_exec(ierr)

    USE wrcomm
    USE wrexecr,wrexecb
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    IF(mode_beam.EQ.0) THEN
       CALL wr_exec_rays(ierr)
    ELSE
       CALL wr_exec_beams(ierr)
    END IF
  END SUBROUTINE wr_exec
END MODULE wrexec
