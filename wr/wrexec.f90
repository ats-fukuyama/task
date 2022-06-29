! wrexec.f90

MODULE wrexec

  PRIVATE
  PUBLIC wr_exec

CONTAINS

  SUBROUTINE wr_exec(ierr)

    USE wrcomm
    USE wrexecr
    USE wrexecb
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    IF(Rmax_wr.EQ.0.D0) Rmax_wr=RR+1.2*RA
    IF(Rmin_wr.EQ.0.D0) Rmin_wr=RR-1.2*RA
    IF(Rmin_wr.LT.0.D0) Rmin_wr=0.D0
    IF(Zmax_wr.EQ.0.D0) Zmax_wr=RR+1.2*rkap*RA
    IF(Zmin_wr.EQ.0.D0) Zmin_wr=RR-1.2*rkap*RA

    IF(mode_beam.EQ.0) THEN
       CALL wr_exec_rays(ierr)
    ELSE
       CALL wr_exec_beams(ierr)
    END IF
  END SUBROUTINE wr_exec
END MODULE wrexec
