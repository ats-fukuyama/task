! wrsetup.f90

MODULE wrsetup

  PRIVATE
  PUBLIC wr_setup

CONTAINS

  SUBROUTINE wr_setup(ierr)

    USE wrcomm
    USE wrsetupr,wrsetupb
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    IF(mode_beam.EQ.0) THEN
       CALL wr_setup_rays(ierr)
    ELSE
       CALL wr_setup_beams(ierr)
    END IF
  END SUBROUTINE wr_setup
END MODULE wrsetup

