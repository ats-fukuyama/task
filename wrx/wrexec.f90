! wrexec.f90

MODULE wrexec

  PRIVATE
  PUBLIC wr_exec

CONTAINS

  SUBROUTINE wr_exec(nstat,ierr)

    USE wrcomm
    USE wrexecr
    USE wrexecb
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: nstat,ierr
    INTEGER:: NS,NSS
    REAL(rkind):: ZA1,ZA2

    IF(Rmax_wr.EQ.0.D0) Rmax_wr=RR+bdr_threshold*RA
    IF(Rmin_wr.EQ.0.D0) Rmin_wr=RR-bdr_threshold*RA
    IF(Rmin_wr.LT.0.D0) Rmin_wr=0.D0
    IF(Zmax_wr.EQ.0.D0) Zmax_wr= bdr_threshold*rkap*RA
    IF(Zmin_wr.EQ.0.D0) Zmin_wr=-bdr_threshold*rkap*RA

!     --- eliminate disp factor for same z/a species ---

    DO NS=1,NSMAX
       NSDP(NS)=1
       DO NSS=1,NS-1
          ZA1=PZ(NS)/PA(NS)
          ZA2=PZ(NSS)/PA(NSS)
          IF(ABS(ZA1-ZA2).LE.1.D-8) NSDP(NS)=0
       ENDDO
    ENDDO

    IF(mode_beam.EQ.0) THEN
       CALL wr_exec_rays(ierr)
       nstat=1
    ELSE
       CALL wr_exec_beams(ierr)
       nstat=2
    END IF
  END SUBROUTINE wr_exec
END MODULE wrexec
